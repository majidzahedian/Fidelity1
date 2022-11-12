#abort in principle is weaker than Ctrl+c --> replace 'abort' by 'Ctrl+c'

#the stage is commanded via the venus programming language. It has an input memory of 256 signs with a
#FIFO structure which only can be surpassed by the Ctrl+C command. 
#So called non-Blocking commands however can be executed during a positioning is executed.

from __future__ import print_function, absolute_import, division
__metaclass__ = type

import visa
import numpy as np

class Micos:
    def __init__(self, logger=None):
        rm = visa.ResourceManager(r"C:\Windows\System32\visa64.dll")
        self.device_xy = rm.open_resource('COM1', write_termination=' ', read_termination='\r\n', baud_rate=57600, query_delay=0.05, timeout=5000)
        self.device_zphi = rm.open_resource('COM4', write_termination=' ', read_termination='\r\n', baud_rate=57600, query_delay=0.05, timeout=5000)
        self.logger = logger

    __MOVEMENT_TIMEOUT__ = 60.

    axes = ['x', 'y', 'z', 'phi']

    positioning_accuracy = {'x': 6, 'y': 6, 'z': 6, 'phi': 6}  #number of significant decimal places for the magnet positioning accuracy. 3 (i.e. accuracy of 1micron) is a conservative estimate

    ignore_range = False

    @property
    def range(self):
        return self.range_r

    @property
    def range_r(self):
        if self.ignore_range:
            range = dict(zip(['x', 'y', 'z'], [{'min': -1e4, 'max': 1e4}] * 3))
        else:
            range = {}
            range['x'] = {'min': 138., 'max': 148.}  # xmin should be 26.5
            range['y'] = {'min': 20., 'max': 104.45}  # physical limit is 159.2984
            range['z'] = {'min': 0.0, 'max': 32.32}
        return range

    @property
    def range_dc(self):
        if self.ignore_range:
            range = dict(zip(['x', 'y', 'z', 'phi'], [{'min': -1e4, 'max': 1e4}] * 4))
        else:
            range = {}
            range['x'] = {'min': 0.0, 'max': 100.}
            range['y'] = {'min': 0.0, 'max': 1000.0}
            range['z'] = {'min': 0.0, 'max': 32.32}
            range['phi'] = {'min': 44., 'max': 46.}
        return range

    block = False

    @property
    def pos_dc(self):
        out = {}
        pxy = self.device_xy.query('pos').split()
        out['x'] = pxy[0]
        out['y'] = pxy[1]
        pzphi = self.device_zphi.query('pos').split()
        out['z'] = pzphi[0]
        out['phi'] = pzphi[1]
        for axis in self.axes:
            out[axis] = round(float(out[axis]), self.positioning_accuracy[axis])
        if not self.is_in_range(out, self.range_dc):
            raise Exception("The current magnet position {} is outside the chosen range {}".format(out, self.range_dc))
        self.cut_to_range(out, self.range_dc)
        return self.round2pa(out)

    @pos_dc.setter
    def pos_dc(self, value):
        target_pos = self.complement_target_pos(value, self.pos_dc)
        if not self.is_in_range(target_pos, self.range_dc):
            raise Exception("The chosen magnet position {} is outside the chosen range {}".format(target_pos, self.range))
        if self.block:
            raise Exception('Magnet movement was blocked by the user.')
        else:
            delta = {}
            pos_dc = self.pos_dc
            for axis in ['x', 'y', 'z', 'phi']:
                delta[axis] = target_pos[axis] - pos_dc[axis]
            self.device_xy.write('%f %f move'%(target_pos['x'], target_pos['y']))
            self.device_zphi.write('%f %f move'%(target_pos['z'], target_pos['phi']))
            info = "Moving magnet to position {}, delta: {} (rc: {})".format(self.round2pa(target_pos), self.round2pa(delta), self.round2pa(self.dc2r(target_pos)))
            if self.logger is not None:
                self.logger.info(info)
            else:
                print(info)

    @property
    def pos_r(self):
        return self.round2pa(self.dc2r(self.pos_dc))

    @pos_r.setter
    def pos_r(self, val):
        val = self.complement_target_pos(val, self.pos_r)
        if not self.is_in_range(val, self.range_r):
            raise Exception("The chosen magnet position {} is outside the chosen range {}".format(val, self.range))
        self.pos_dc = self.r2dc(val)

    @property
    def pos(self):
        return self.pos_r

    @pos.setter
    def pos(self, val):
        self.pos_r = val

    def dc2r(self, pos_dc):
        posxyz_arr = np.dot(self.t_matrix, np.array([pos_dc['x'], pos_dc['y'], pos_dc['z']]))
        return dict(zip(['x', 'y', 'z'], posxyz_arr))

    def r2dc(self, pos_r):
        posxyz_arr_t = np.dot(self.t_matrix.T, np.array([pos_r[i] for i in ['x', 'y', 'z']]))
        return dict(zip(['x', 'y', 'z', 'phi'], list(posxyz_arr_t) + [self.pos_dc['phi']]))

    def complement_target_pos(self, target_pos, current_pos):
        return dict([[axis, float(target_pos.get(axis, current_pos[axis]))] for axis in current_pos.keys()])

    def is_in_range(self, pos, range):
        in_range = True
        for axis in pos:
            tol = 10 ** -self.positioning_accuracy[axis]
            d_min = range[axis]['min'] - pos[axis]
            d_max = pos[axis] - range[axis]['max']
            if d_min > tol or d_max > tol:
                print('Axis {} is out of range ({}, {}).'.format(axis, d_min, d_max))
                in_range = False
        return in_range

    def cut_to_range(self, pos, range):
        for axis in pos:
            if pos[axis] < range[axis]['min']:
                pos[axis] = range[axis]['min']
                print("Position {} of axis {} is out of the range {} but within the positioning tolerance of {} digits. It was set to {}".format(pos[axis], axis, range[axis], self.positioning_accuracy[axis], range[axis]['min']))
            elif pos[axis] > range[axis]['max']:
                pos[axis] = range[axis]['max']
                print("Position {} of axis {} is out of the range {} but within the positioning tolerance of {} digits. It was set to {}".format(pos[axis], axis, range[axis], self.positioning_accuracy[axis], range[axis]['max']))
        return pos

    @property
    def t_matrix(self):
        por = np.deg2rad(self.pos_dc['phi'])
        return np.array([[np.cos(por), np.sin(por), 0],
                         [-np.sin(por), np.cos(por), 0],
                         [0, 0, 1]])

    def round2pa(self, pos):
        for axis in pos.keys():
            pos[axis] = round(float(pos[axis]), self.positioning_accuracy[axis])
        return pos

    def make_step(self, axis, stepsize):
        """
        Moves the Magnet a distance of stepsize [mm or degree] along the chosen axis 'x','y','z','phi'
        """
        self.pos = self.complement_target_pos({axis: self.pos[axis] + stepsize}, self.pos)

    def make_step_dc(self, axis, stepsize):
        """
        Moves the Magnet a distance of stepsize [mm or degree] along the chosen axis 'x','y','z','phi'
        """
        self.pos_dc = self.complement_target_pos({axis: self.pos_dc[axis] + stepsize}, self.pos_dc)

    def calibrate_all(self):
        self.calibrate_x()
        self.calibrate_y()
        self.calibrate_z()
        self.calibrate_phi()

    def calibrate_x(self):
        self.device_xy.write('1 1 setaxis')
        self.device_xy.write('4 2 setaxis')
        self.device_xy.write('cal')
        self.limits = self.range_dc

    def calibrate_y(self):
        self.device_xy.write('4 1 setaxis')
        self.device_xy.write('1 2 setaxis')
        self.device_xy.write('cal')
        self.limits = self.range_dc


    def calibrate_z(self):
        self.device_zphi.write('1 1 setaxis')
        self.device_zphi.write('4 2 setaxis')
        self.device_zphi.write('cal')
        self.limits = self.range_dc

    def calibrate_phi(self):
        self.device_zphi.write('4 1 setaxis')
        self.device_zphi.write('1 2 setaxis')
        self.device_zphi.write('cal')
        self.limits = self.range_dc

    @property
    def limits(self):
        x = self.device_xy.query("getlimit").split()  # two line answer
        y = self.device_xy.query("").split()
        z = self.device_zphi.query("getlimit").split()  # two line answer
        phi = self.device_zphi.query("").split() #two line answer
        range = {}
        range['x'] = {'min': float(x[0]), 'max': float(x[1])}  # xmin should be 26.5
        range['y'] = {'min': float(y[0]), 'max': float(y[1])}
        range['z'] = {'min': float(z[0]), 'max': float(z[1])}
        range['phi'] = {'min': float(phi[0]), 'max': float(phi[1])}
        return range

    @limits.setter
    def limits(self, val):
        self.device_xy.write(
            "{} {} {} {} setlimit".format(val['x']['min'], val['y']['min'], val['x']['max'], val['y']['max']))
        self.device_zphi.write(
            "{} {} {} {} setlimit".format(val['z']['min'], val['phi']['min'], val['z']['max'], val['phi']['max']))

    def abort(self):
        self.device_xy.write('abort')
        self.device_zphi.write('abort')


if __name__ == '__main__':
    m = Micos()
