# from __future__ import print_function, absolute_import, division #careful, probably unicode_literals make parts disfunctional
__metaclass__ = type
import visa
import time

class Amp:

    def __init__(self):
        rm = visa.ResourceManager()
        self.device = rm.open_resource('GPIB0::1', timeout=5000, read_termination='\n') #GPBI
        # self.device = rm.open_resource('USB0::0x0547::0x1B58::333074::INSTR', send_end=False, query_delay=0.1, read_termination='\n') #USB

    __ERROR_DICT__ = {'00': 'No Fault', '01': 'AC Interlock', '02': 'Interlock','03': 'PS3','04': 'RS 485 Error', '05': 'PS4', '06': 'ALC', '07': 'Under Current', '08': 'Over Current',
                      '09': 'No Resp. to External Fault', '0A': '50V lost on module', '0B': 'Thermal A3', '0C': 'Thermal A4', '0D': 'Thermal A5', '0E': 'Thermal A6'}

    @property
    def errors(self):
        return self.__ERROR_DICT__[self.device.query('FSTA?')[-2:]]

    def check_errors(self):
        e = self.errors
        if e != 'No Fault':
            raise Exception('Error {} occured on the device.'.format(e))

    def query(self, cmd):
        e = self.errors
        if e != 'No Fault':
            raise Exception('Error {} occured on the device PRIOR TO querying'.format(e))
        q = self.device.query(cmd).replace("{}= ".format(cmd.replace('?', '')), '')
        e = self.errors
        if e != 'No Fault':
            print(q)
            raise Exception('Error {} occured on the device AFTER querying.'.format(e))
        else:
            return q

    def write(self, command):
        self.device.write(command)
        if not command in 'POWER':
            self.device.read()
        self.check_errors()

    @property
    def state(self):
        def b(si):
            return format(int(bin(int(si))[2:]), '04d')

        s = self.query('STATE?')
        out = {}
        xb = b(s[0])
        out['pulse_status'] = {0: 'OFF', 1: 'PULSE'}[int(xb[3])]
        out['blank_selection'] = {0: 'OFF', 1: 'INHIBITED'}[int(xb[2])]
        out['blank'] = {0: 'DISABLED', 1: 'ENABLED'}[int(xb[1])]
        out['remote_control'] = {0: 'DISABLED', 1: 'ENABLED'}[int(xb[0])]
        yb = b(s[1])
        out['power_status'] = {0: 'OFF', 1: 'ON'}[int(yb[3])]
        out['standby_status'] = {0: 'OFF', 1: 'STANDBY'}[int(yb[2])]
        out['rf_on'] = {0: 'OFF', 1: 'ON'}[int(yb[1])]
        out['fault_status'] = {0: 'OFF', 1: 'FAULT EXISTS'}[int(yb[0])]
        zb = b(s[2])
        out['keylock_inhibit'] = {0: 'OFF', 1: 'INHIBITED'}[int(zb[3])]
        ab = b(s[3])
        out['cw_mode'] = {0: 'DISABLED', 1: 'ENABLED'}[int(ab[3])]
        return out

    def check_set(self, prop, val):
        a = getattr(self, prop)
        if a != val:
            print("a: {}, val:{}".format(a,val))
            raise Exception('{} could not be set to value {}'.format(prop, val))

    @property
    def rf_gain(self):
        rg = int(self.query('RFG?'))
        return rg

    @rf_gain.setter
    def rf_gain(self, val):
        if type(val) is not int:
            raise Exception('Error: val = {}'.format(val))
        self.write('LEVEL:GAIN{}'.format(val))
        self.check_set('rf_gain', val)

    @property
    def rf_on(self):
        return {'ON': True, 'OFF': False}[self.state['rf_on']]

    @rf_on.setter
    def rf_on(self, val):
        self.write('RF:{}'.format({True: 'ON', False: 'OFF'}[val]))
        import numpy as np
        for i in np.arange(0, 20, 0.3): # the display updates very fast, but the query for rf_on takes some seconds to updated
            time.sleep(i)
            if self.rf_on == val:
                break
        self.check_set('rf_on', val)

    @property
    def power_on(self):
        return {'ON': True, 'OFF': False}[self.state['power_status']]

    @power_on.setter
    def power_on(self, val):
        self.write('POWER:{}'.format({True: 'ON', False: 'OFF'}[val]))
        self.check_set('power_on', val)

    @property
    def forward_power(self):
        return float(self.query('FPOW?'))

    @property
    def reverse_power(self):
        return float(self.query('RPOW?'))

    def start(self):
        self.power_on = True
        self.rf_on = True

    def stop(self):
        self.power_on = False

if __name__ == '__main__':
    a = Amp()