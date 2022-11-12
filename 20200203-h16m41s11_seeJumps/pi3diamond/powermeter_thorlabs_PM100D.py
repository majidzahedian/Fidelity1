from __future__ import print_function, absolute_import, division

from pi3diamond import pi3d

from traits.api import *
from traitsui.api import *

nidaq = pi3d.nidaq


class Powermeter(HasTraits):
    
    volt_min = Float(0.,label = 'Minimum voltage [V]') # given in Volts
    volt_max = Float(2.,label = 'Maximum voltage [V]') 
    power_min = Float(0.,label = 'Minimum power [mW]') #given in mW
    power_max = Float(16.,label = 'Maximum Power [mW]')
    
    power = Float(0., label = 'Power [mW]')    
    get_power_button = Button()
    
    traits_view = View(
        VGroup(
            HGroup(
                Item('volt_min'),
                Item('volt_max'),
            ),
            HGroup(
                Item('power_min'),
                Item('power_max'),
            ),
        ),
        Item('get_power_button'),
        Item('power'),
        title='Powermeter',
        kind = 'live',
        width=400, 
        height=400,
        resizable=True, 
        buttons = ['OK','Undo'],
    )
    
    def _get_power_button_fired(self):
        self.power = self.get_power()
    
    def get_power(self):
        """returns measured power in mW"""
        voltage = nidaq.read_powermeter_ai()
        power = voltage[0][0]/(self.volt_max - self.volt_min) * (self.power_max - self.power_min)
        return power
