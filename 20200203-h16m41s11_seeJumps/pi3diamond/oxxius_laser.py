from __future__ import print_function, absolute_import, division

import serial
import visa

class Laser:

    def __init__(self):
        self.laser = serial.Serial('COM3', baudrate=19200, timeout=5.) #serial needs port (i.e. COM3), not the device name, which is ok for visa (see magnet_stage_micos)

    @property
    def analog_control_mode(self):
        """
        Laser mode

        0 Power (APC)
        1 Current (ACC)
        """
        return int(self.ask('ANA'))

    @analog_control_mode.setter
    def analog_control_mode(self, value):
        self.write('ANA', str(value))

    @property
    def analog_modulation_state(self):
        """
        External analog modulation input state

        1 external
        0 internal
        """
        return int(self.ask('AM'))

    @analog_modulation_state.setter
    def analog_modulation_state(self, value):
        self.write('AM', str(value))

    @property
    def base_plate_temp(self):
        """
        Retrieves the measured temperature of the laser base plate
        """
        return float(self.ask('BT'))

    @property
    def laser_diode_current(self):
        """
        Measured current of the laser diode
        Ranges from 0 to self.maximum_laser_current
        """
        return float(self.ask('C'))

    @laser_diode_current.setter
    def laser_diode_current(self, value):
        self.write('C', str(value))

    def _set_laser_diode_current_fast(self, value):
        self.write('CM', str(value))

    laser_diode_current_fast = property(fset=_set_laser_diode_current_fast, doc="Sets the laser diode current without saving the value in memory. \n"
                                                                                "Use this command instead of 'C', i.e. the 'laser_diode_current' attribute \n"
                                                                                "if you want to modify the current in real-time.\n"
                                                                                "Maximum frequency : 50Hz \n"
                                                                                "(RS-232) or 1kHz (USB)")

    @property
    def cdrh_state(self):
        """
        CDRH state

        0 Delay Off
        1 Delay On
        """
        return int(self.ask('CDRH'))

    @cdrh_state.setter
    def cdrh_state(self, value):
        self.write('CDRH', str(value))

    @property
    def modulation_state(self):
        """
        Modulation state (SMB connector)

        1 CW
        0 Modulated
        """
        return int(self.ask('CW'))

    @modulation_state.setter
    def modulation_state(self, value):
        self.write('CW', str(value))

    @property
    def diode_temp_set_point(self):
        """
        Retrieves the set temperature of the laser diode
        """

        return float(self.ask('DST'))

    @property
    def measured_diode_temp(self):
        """
        Retrieves the measured temperature of the laser diode
        """
        return float(self.ask('DT'))

    @property
    def measured_module_temp(self):
        """
        Retrieves the measured temperature of the module
        """
        return float(self.ask('ET'))

    @property
    def fault_number(self):
        """
        Retrieves current fault
        0 : No alarm
        1 : Diode current
        2 : Laser power
        3 : Power supply
        4 : Diode temperature
        5 : Base temperature
        6 : Warning end of Life
        """
        return int(self.ask('F'))

    @property
    def cumulated_time_of_operation(self):
        """
        Cumulated time of operation
        Tracking interval : 1 minute.
        """
        return float(self.ask('HH'))

    @property
    def serial_number_and_wavelength(self):
        """
        Retrieves the head serial
        number and nominal
        wavelength, comma delimited.
        """
        return self.ask('HID')

    @property
    def interlock_state(self):
        """
        Interlock state
        0 Open
        1 Closed
        """
        return int(self.ask('INT'))

    @property
    def input_voltage(self):
        """
        Allowed values from 5 to 6.5
        """
        return float(self.ask('IV'))

    @property
    def laser_emission_activation(self):
        """
        Laser emission activation state
        """
        return int(self.ask('L'))

    @laser_emission_activation.setter
    def laser_emission_activation(self, value):
        self.write('L', str(value))

    @property
    def maximum_laser_current(self):
        """
        Maximum laser current
        Range 0 to 1000
        """
        return float(self.ask('MAXLC'))

    @maximum_laser_current.setter
    def maximum_laser_current(self, value):
        self.write('MAXLC', str(value))

    @property
    def maximum_laser_power(self):
        """
        Maximum laser power
        Range from 0 to 500
        """
        return float(self.ask('MAXLP'))

    @maximum_laser_power.setter
    def maximum_laser_power(self, value):
        self.write('MAXLP', str(value))

    @property
    def laser_output_power(self):
        """
        Laser output power
        """
        return float(self.ask('P'))


    @laser_output_power.setter
    def laser_output_power(self, value):
        self.write('P', str(value))

    def _set_laser_output_power_fast(self, value):
        self.write('PM', str(value))

    laser_diode_current_fast = property(fset=_set_laser_output_power_fast, doc="Sets the laser power without saving the value in memory. \n"
                                                                                "Use this command instead of 'P', i.e. the 'laser_output_power' attribute \n"
                                                                                "if you want to modify the current in real-time.\n"
                                                                                "Maximum frequency : 50Hz \n"
                                                                                "(RS-232) or 1kHz (USB)")

    @property
    def processor_temperature(self):
        """
        Retrieves the measured
        temperature of the
        microcontroller
        """
        return float(self.ask('PST'))

    @property
    def current_set_point(self):
        """
        Retrieves the laser current
        setting
        """
        return float(self.ask('SC'))

    @property
    def power_set_point(self):
        """
        Retrieves the laser power
        setting
        """
        return float(self.ask('SP'))

    @property
    def operation_status(self):
        """
        Retrieves the system operating
        status
        1 : Warm Up
        2 : Standby
        3 : Laser ON
        4 : Error
        5 : Alarm
        6 : Sleep
        7 : Searching SLM point
        """
        return int(self.ask('STA'))

    def alarm_reset(self):
        print(self.write('RST', '1'))

    @property
    def software_version(self):
        """
        Laser firmware version
        """
        return self.ask('SV')

    @property
    def tec_state(self):
        """
        TEC status

        0 TEC disabled
        1 TEC enabled
        If TEC = 0, ?STA = Sleep
        """
        return int(self.ask('T'))

    @tec_state.setter
    def tec_state(self, value):
        self.write('T', str(value))

    def _set_digital_modulation_state(self, value):
        self.write('TTL', str(value))

    digital_modulation_state = property(fset=_set_digital_modulation_state, doc="TTL 1 External\nTTL 0 Internal. This is readonly, as 'CW', i.e.the 'modulation_state' property is equivalent (opposite).")

    def ask(self, command):
        self.laser.write(str("?{}\n".format(command)))
        return self.laser.readline().replace('\r\n', '')


    def write(self, command, value):
        self.laser.write(str("{} {}\n".format(command, value)))
        return self.laser.readline().replace(command, '').replace('=', '').replace('\r\n', '')
