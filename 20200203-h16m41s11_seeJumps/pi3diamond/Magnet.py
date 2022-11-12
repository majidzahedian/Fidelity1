from __future__ import print_function, absolute_import, division


from utility import GetSetItemsMixin
from pi3diamond import pi3d

import misc

from traits.api import *
from traitsui.api import *


import datetime

class PositionLabelCloseHandler(Handler):
    """
    This handler makes sure, that no multiple labels for Magnet positions can be assigned.
    """

    def close(self, info, is_ok):
        magnet = pi3d.magnet
        for pos in magnet.positions_list:
            if pos.label == info.object.label:
                print("This label already exists. Chose another one.")
                return
            else:
                pi3d.logger.info("Magnet position x: {},y: {}, z: {} saved as {}".format(info.object.x, info.object.y, info.object.z, info.object.label))
        return Handler.close(self, info, is_ok)

class PositionsRow(HasTraits):
    def __init__(self, label=None):
        super(PositionsRow, self).__init__()
        now = datetime.datetime.now()
        if label is None:
            self.edit_traits()
        else:
            self.label=label
        self.date = now.strftime("%Y-%m-%d %H:%M")
        self.x = pi3d.magnet.x
        self.y = pi3d.magnet.y
        self.z = pi3d.magnet.z

    may_exist = Bool(True)
    label = Str(label='Label', auto_set=False, enter_set=False)
    date = Str(label='Date')
    x = Float(label='x')
    y = Float(label='y')
    z = Float(label='z')
    view_label = View(
        Item('label'),
        title='Enter a label for this magnet position',
        kind='modal',
        resizable=True,
        buttons=['OK', 'Cancel'],
        handler=PositionLabelCloseHandler
    )
    model_selection = 0


positions_editor = TableEditor(
    columns=[ObjectColumn(name='label', width=0.30, editable=False),
             ObjectColumn(name='date', width=0.4, editable=False),
             ObjectColumn(name='x', width=0.20, editable=False),
             ObjectColumn(name='y', width=0.20, editable=False),
             ObjectColumn(name='z', width=0.20, editable=False),
    ],
    deletable=True,
    auto_size=False,
    orientation='vertical',
    selection_mode='row',
    show_toolbar=True,
    row_factory=PositionsRow,
    selected='selected_row'
)

class Magnet(SingletonHasTraits, GetSetItemsMixin):
    """Magnet positioning software.
    Enables precise positioning of the magnetic field with the micos stage around x,y,z.
    Provides automatized alignment of the magnetic field with respect to the NV axis depending on the fluorescence intensity (rough) or on T1 of the nitrogen nuclear spin (fine).
    """

    def __init__(self):
        super(Magnet, self).__init__()
        self.do_refocus_on_step = False
        self._write_pos_to_gui(pi3d.magnet_stage_micos.pos)


    # Set Magnet range view
    x = Range(pi3d.magnet_stage_micos.range['x']['min'], pi3d.magnet_stage_micos.range['x']['max'], 0.0, label='x [mm]', mode='slider')
    y = Range(pi3d.magnet_stage_micos.range['y']['min'], pi3d.magnet_stage_micos.range['y']['max'], 0.0, label='y [mm]', mode='slider')
    z = Range(pi3d.magnet_stage_micos.range['z']['min'], pi3d.magnet_stage_micos.range['z']['max'], 0.0, label='z [mm]', mode='slider')


    do_refocus_on_step = misc.ret_property_typecheck('do_refocus_on_step', bool)

    # General
    block_magnet_bool = Bool(value=False, label='Block all magnet movement')
    abort_button = Button('Abort')
    move_button = Button('Move')
    load_position_button = Button('Load selected Position')
    get_position_button = Button('Get magnet position')

    positions_list = List(PositionsRow)
    selected_row = Any()

    # Calibration
    calibrate_button = Button('Calibrate magnet position')

    # Steps
    stepsize_x = Float(0.)
    stepsize_y = Float(0.)
    stepsize_z = Float(0.)
    make_step_x_plus_button = Button('x+')
    make_step_x_minus_button = Button(label='x-')
    make_step_y_plus_button = Button('y+')
    make_step_y_minus_button = Button('y-')
    make_step_z_plus_button = Button('z+')
    make_step_z_minus_button = Button('z-')

    def move(self, target_pos=None):
        if target_pos is None:
            pi3d.magnet_stage_micos.pos = {'x': self.x, 'y': self.y, 'z': self.z}
        else:
            pi3d.magnet_stage_micos.pos = pi3d.magnet_stage_micos.complement_target_pos(target_pos, pi3d.magnet_stage_micos.pos)
        self.save_magnet_pos_single_values()


    def save_magnet_pos_single_values(self):
        pi3d.save_values_to_file([pi3d.magnet_stage_micos.pos_dc[axis] for axis in ['x', 'y', 'z', 'phi']], 'magnet_pos')

    def _move_button_fired(self):
        self.move()

    def _get_position_button_fired(self):
        self._write_pos_to_gui(pi3d.magnet_stage_micos.pos)

    def _load_position_button_fired(self):
        self._write_pos_to_gui({'x': self.selected_row.x, 'y': self.selected_row.y, 'z': self.selected_row.z})

    def _calibrate_button_changed(self):
        pi3d.magnet_stage_micos.calibrate_all()

    @on_trait_change('make_step_x_plus_button, make_step_x_minus_button, \
                      make_step_y_plus_button, make_step_y_minus_button, \
                      make_step_z_plus_button, make_step_z_minus_button')
    def test(self, name, new):
        for axis in pi3d.magnet_stage_micos.axes:
            for sign in ['plus', 'minus']:
                if axis in name and sign in name:
                    stepsize = getattr(self, "stepsize_{}".format(axis))
                    if sign == 'minus':
                        stepsize *= -1
                    pi3d.magnet_stage_micos.make_step(axis, stepsize)
        if self.do_refocus_on_step:
            pi3d.confocal.state='refocus_crosshair'



    def _write_pos_to_gui(self, pos):
        pos = pi3d.magnet_stage_micos.complement_target_pos(pos, pi3d.magnet_stage_micos.pos)
        for axis in pos.keys():
            setattr(self, axis, pos[axis])

    def _block_magnet_bool_changed(self):
        if self.block_magnet_bool:
            pi3d.magnet_stage_micos.abort()
            pi3d.magnet_stage_micos.block = True
        else:
            pi3d.magnet_stage_micos.block = False

    def _abort_button_fired(self):
        self.block_magnet_bool = True
        pi3d.magnet_stage_micos.abort()  #just in case there is a bug in _block_magnet_bool_changed

    def __getstate__(self):
        state = SingletonHasTraits.__getstate__(self)
        for key in []:
            if state.has_key(key):
                del state[key]
        return state

    traits_view = View(
        VGroup(
            Tabbed(
                VGroup(
                    HGroup(
                        Item('x', width=650),
                        spring,
                        Item('make_step_x_minus_button', label='step x', show_label=True, width=-30),
                        Item('stepsize_x', show_label=False, width=-30),
                        Item('make_step_x_plus_button', label='x+', show_label=False, width=-30)
                    ),
                    HGroup(
                        Item('y', width=650),
                        spring,
                        Item('make_step_y_minus_button', label='step y', show_label=True, width=-30),
                        Item('stepsize_y', show_label=False, width=-30),
                        Item('make_step_y_plus_button', label='y+', show_label=False, width=-30)
                    ),
                    HGroup(
                        Item('z', width=650),
                        spring,
                        Item('make_step_z_minus_button', label='step z', show_label=True, width=-30),
                        Item('stepsize_z', show_label=False, width=-30),
                        Item('make_step_z_plus_button', label='z+', show_label=False, width=-30)
                    ),
                    Item('move_button', show_label=False),
                    Item('get_position_button'),
                    # Item('positions_list', show_label=False, editor=positions_editor),
                    Item('load_position_button', label='Load selected Position', show_label=False),
                    label='Position'
                ),
            ),
            HGroup(
                Item('calibrate_button', show_label=False),
                Item('block_magnet_bool'),
                Item('abort_button', width=600)
            ),
        ),
        title='Magnet',
        kind='live',
        width=900,
        height=700,
        resizable=True,
        buttons=['OK', 'Undo'],
        # handler=CloseHandler
    )