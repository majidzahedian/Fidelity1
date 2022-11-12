from __future__ import print_function, absolute_import, division
__metaclass__ = type

import numpy as np
np.set_printoptions(linewidth=500, suppress=True)
import os
os.chdir(r'D:/Python/pi3diamond')
import sys

from traits.etsconfig.api import ETSConfig
ETSConfig.toolkit = 'wx' #possible values: 'wx', 'qt4'

import PyQt5.QtWidgets

app = PyQt5.QtWidgets.QApplication(sys.argv)
if 'pi3d' not in globals():
    from pi3diamond_custom import pi3d

from qutip_enhanced import *

if __name__ == '__main__':
    pi3d.gated_counter.gui.show()
    pi3d.gui.show()
    pi3d.logger.info('Starting pi3diamond...')

    module_names = ['mcas_dict', 'tt', 'confocal', 'odmr', 'gated_counter']#, 'magnet', 'pulsed']

    pi3d.show_gui(module_names)
    import UserScripts.helpers.snippets_awg as sna
    pi3d.oxxius_laser._set_laser_diode_current_fast(sna.__CURRENT_POL_RED__)
    pi3d.add_to_queue('standard_awg_sequences', folder='D:/Python/pi3diamond/UserScripts/helpers')
    pi3d.add_to_queue('refocus_confocal_odmr', folder='D:/Python/pi3diamond/UserScripts')
    # pi3d.user_script_folder = r'D:\Python\pi3diamond\UserScripts\QFT'
    # pi3d.selected_user_script = 'qutrit_qft'
    # # pi3d.user_script_folder = r'D:\Python\pi3diamond\UserScripts\parameters'
    # # pi3d.selected_user_script = 'electron_rabi_pulsed'
    # pi3d.add_to_queue()