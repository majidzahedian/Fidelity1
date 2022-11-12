# coding=utf-8
from __future__ import print_function, absolute_import, division  # unicode_literals leads to an undetermined bug

__metaclass__ = type

import ctypes

# lib = ctypes.CDLL('./libfoo.so')


import numpy
import time
import numpy as np
import PyDAQmx

class NIDAQ():
    """NIDAQ measurement card.
    2 counters are used for standard photon counting (time trace).
    The other 2 counters are used on demand for functions like scanning, odmr etc."""
    _MaxCounts = 1e7
    _RWTimeout = 60

    def __init__(self, photon_source, function_counter_in, function_counter_out,
                 scanner_ao_channels, scanner_xrange, scanner_yrange, scanner_zrange, odmr_trig_channel=None, scanner_ai_channels=None):
        self._photon_source = photon_source  # channel where apd is conntected
        self._odmr_trig_channel = odmr_trig_channel  # channel for mw trigger during odmr
        self._function_counter_in = function_counter_in  # counter for seperate function / photon counting (eg scanning, odmr)
        self._function_counter_out = function_counter_out  # counter for pulse generation for seperate function
        self._scanner_xrange = scanner_xrange  # x range of scanner as (xmin, xmax)
        self._scanner_yrange = scanner_yrange
        self._scanner_zrange = scanner_zrange
        self._scanner_xvolt_range = (10, -10)
        self._scanner_yvolt_range = (-10, 10)
        self._scanner_zvolt_range = (-10, 10)
        self.calibrated = False
        self._scanner_ao_channels = scanner_ao_channels  # scanner ao channels for x,y,z control
        self._scanner_ai_channels = scanner_ai_channels  # scanner ai channels for x,y,z readout

        self.function_state = 'idle'  # set to 'scan', 'odmr', during scanning, odmr. before starting function, check if state is idle

        self.scanner_x = 0.0  # current x, y, z values of scanner
        self.scanner_y = 0.0
        self.scanner_z = 0.0

        self.scanner_ao_task = PyDAQmx.Task()
        self.scanner_ao_task.CreateAOVoltageChan(self._scanner_ao_channels, '',  # use sanncer ao_channels, name = ''
                                                 -10.,  # min voltage
                                                 10.,  # max voltage
                                                 PyDAQmx.DAQmx_Val_Volts, '')
        self.scanner_ao_task.SetSampTimingType(PyDAQmx.DAQmx_Val_OnDemand)

        # init ai channels / task to read scanner position
        if self._scanner_ai_channels != None:
            self.scanner_ai_task = PyDAQmx.Task()
            self.scanner_ai_task.CreateAIVoltageChan(
                self._scanner_ai_channels,  # use scanner ao_channels
                '',  # use scanner ao_channels name = ''
                PyDAQmx.DAQmx_Val_RSE,  # measuring against ground? DAQmx_Val_Cfg_Default ?
                -10.,  # min voltage
                10.,  # max voltage
                PyDAQmx.DAQmx_Val_Volts,  # units is Volt
                '')  # use no costum scaling
            self.scanner_ai_task.SetSampTimingType(PyDAQmx.DAQmx_Val_OnDemand)  # set task timing to on demand, i.e. when demanded by software

    # -------------- Scanner ----------------------------------
    def scanner_pos_to_volt(self, Pos):
        x = self._scanner_xrange
        y = self._scanner_yrange
        z = self._scanner_zrange
        return numpy.vstack(
            (
                -1. * (20.0 / (x[1] - x[0]) * (Pos[0] - x[0]) - 10),
                20.0 / (y[1] - y[0]) * (Pos[1] - y[0]) - 10,
                -1. * (20.0 / (z[1] - z[0]) * (Pos[2] - z[0]) - 10)
            )
        )

    def scanner_volt_to_pos(self, Volt):
        xpoint1 = (self._scanner_xvolt_range[0], self._scanner_xrange[0])  # (u1,x1)
        xpoint2 = (self._scanner_xvolt_range[1], self._scanner_xrange[1])  # (u2,x2)

        ypoint1 = (self._scanner_yvolt_range[0], self._scanner_yrange[0])  # (u1,z1)
        ypoint2 = (self._scanner_yvolt_range[1], self._scanner_yrange[1])  # (u2,z2)

        zpoint1 = (self._scanner_zvolt_range[0], self._scanner_zrange[0])  # (u1,z1)
        zpoint2 = (self._scanner_zvolt_range[1], self._scanner_zrange[1])  # (u2,z2)

        return numpy.vstack((xpoint1[1] + ((Volt[0] - xpoint1[0]) * (xpoint2[1] - xpoint1[1]) / (xpoint2[0] - xpoint1[0])),
                             ypoint1[1] + ((Volt[1] - ypoint1[0]) * (ypoint2[1] - ypoint1[1]) / (ypoint2[0] - ypoint1[0])),
                             zpoint1[1] + ((Volt[2] - zpoint1[0]) * (zpoint2[1] - zpoint1[1]) / (zpoint2[0] - zpoint1[0]))
                             ))

    def calibrate_scanner_remote(self):
        voltrange = []
        self.scanner_set_pos(self._scanner_xrange[0], self._scanner_yrange[0], self._scanner_zrange[0])
        time.sleep(0.5)

        voltrange.append(self.read_scanner_ai())

        self.scanner_set_pos(self._scanner_xrange[1], self._scanner_yrange[1], self._scanner_zrange[1])

        time.sleep(0.5)

        voltrange.append(self.read_scanner_ai())

        self._scanner_xvolt_range = []
        self._scanner_yvolt_range = []
        self._scanner_zvolt_range = []

        self._scanner_xvolt_range.extend([voltrange[0][0], voltrange[1][0]])
        self._scanner_yvolt_range.extend([voltrange[0][1], voltrange[1][1]])
        self._scanner_zvolt_range.extend([voltrange[0][2], voltrange[1][2]])
        self.calibrated = True

    def write_scanner_ao(self, data, length=1, start=False):
        self.scanner_ao_task.WriteAnalogF64(
            length,
            start,
            self._RWTimeout,
            PyDAQmx.DAQmx_Val_GroupByChannel,
            data,
            None,
            None
        )

    def read_scanner_ai(self):
        line_points = 1
        data = numpy.empty((3, line_points), dtype=numpy.float64)
        self.scanner_ai_task.ReadAnalogF64(

            line_points,
            self._RWTimeout,
            PyDAQmx.DAQmx_Val_GroupByChannel,
            data,
            3 * line_points,
            None,
            None
        )
        return data

    def get_scanner_pos(self):
        if not self.calibrated:
            self.calibrate_scanner_remote()
        return self.scanner_volt_to_pos(self.read_scanner_ai())

    def scanner_setx(self, x):
        """Move stage to x, y, z"""
        if self.function_state != 'scan':
            self.write_scanner_ao(
                self.scanner_pos_to_volt((x, self.scanner_y, self.scanner_z)),
                start=True
            )
        self.scanner_x = x

    def scanner_sety(self, y):
        """Move stage to x, y, z"""
        if self.function_state != 'scan':
            self.write_scanner_ao(
                self.scanner_pos_to_volt((self.scanner_x, y, self.scanner_z)),
                start=True
            )
        self.scanner_y = y

    def scanner_setz(self, z):
        """Move stage to x, y, z """
        if self.function_state != 'scan':
            self.write_scanner_ao(
                self.scanner_pos_to_volt((self.scanner_x, self.scanner_y, z)),
                start=True
            )
        self.scanner_z = z

    def scanner_set_pos(self, x, y, z):
        """Move stage to x, y, z"""
        if self.function_state != 'scan':
            self.write_scanner_ao(self.scanner_pos_to_volt((x, y, z)), start=True)
        self.scanner_x, self.scanner_y, self.scanner_z = x, y, z

    def init_scan(self, settling_time, count_time):
        """set up tasks for scanning."""
        self.function_state = 'scan'
        self.scanner_pulse_task = PyDAQmx.Task()  # ctypes.c_uint64()  #create handle for task, this task will generate pulse signal for photon counting
        self.scanner_count_task = PyDAQmx.Task()  # ctypes.c_uint64() #this task will count photons with binning defined by pulse task
        self.scan_freq = 1. / (count_time + settling_time)  # counting freqeuency
        duty_cycle = self.scan_freq * count_time  # width of pulse (relative to count_interval) during which photons are counted
        self.scan_duty_cycle = duty_cycle
        # now create pulse signal defined by f and duty_cycle
        self.scanner_pulse_task.CreateCOPulseChanFreq(
            self._function_counter_out,  # use this counter
            '',                       #name
            PyDAQmx.DAQmx_Val_Hz,
            PyDAQmx.DAQmx_Val_Low,
            0,  # units, idle state, inital delay
            self.scan_freq,   #pulse frequency
            duty_cycle
        )  # duty cycle of pulses
        # set up pulse width measurement in photon ticks, i.e. the width of each pulse generated by pulse_out_task is measured in photon ticks.
        self.scanner_count_task.CreateCIPulseWidthChan(self._function_counter_in,  # use this counter
                                                       '',  # name
                                                       0,  # expected minimum value
                                                       self._MaxCounts * duty_cycle / self.scan_freq,  # expected maximum value
                                                       PyDAQmx.DAQmx_Val_Ticks,  # units of width measurement, here photon ticks
                                                       PyDAQmx.DAQmx_Val_Rising, '')  # start pulse width measurement on rising edge
        # set the pulses to counter self._trace_counter_in
        self.scanner_count_task.SetCIPulseWidthTerm(self._function_counter_in, self._function_counter_out + 'InternalOutput')
        # set the timebase for width measurement as self._photon_source
        self.scanner_count_task.SetCICtrTimebaseSrc(self._function_counter_in, self._photon_source)
        self.line_points = -1  # reset variable, scan length is checked and adjusted for each line scan
        return 0

    def stop_scan(self):
        "Clear tasks, so that counters are not in use any more."
        self.function_state = 'idle'
        self.scanner_ao_task.SetSampTimingType(PyDAQmx.DAQmx_Val_OnDemand)

    def set_scan_length(self, N):
        """Set length for line scan, i.e. number of clock pulses, number of count samples etc."""
        self.scanner_ao_task.SetSampTimingType(PyDAQmx.DAQmx_Val_SampClk)  # set task timing to use a sampling clock
        if N < numpy.inf:
            self.scanner_ao_task.CfgSampClkTiming(
                self._function_counter_out + 'InternalOutput',  # use these pulses as clock
                self.scan_freq,  # maximum expected clock frequency
                PyDAQmx.DAQmx_Val_Falling,
                PyDAQmx.DAQmx_Val_FiniteSamps,  # genarate sample on falling edge, generate finite number of samples
                N)  # samples to generate
            # set timing for scanner pulse and count task.
            # Because clock pulse starts high, but position is set on falling edge, first count will be ignored
            self.scanner_pulse_task.CfgImplicitTiming(
                PyDAQmx.DAQmx_Val_FiniteSamps,
                N + 1
            )
            self.scanner_count_task.CfgImplicitTiming(
                PyDAQmx.DAQmx_Val_FiniteSamps,
                N + 1
            )
        # read samples from beginning of acquisition, do not overwrite
        self.scanner_count_task.SetReadRelativeTo(PyDAQmx.DAQmx_Val_CurrReadPos)
        self.scanner_count_task.SetReadOffset(1)  # do not read first sample
        self.scanner_count_task.SetReadOverWrite(PyDAQmx.DAQmx_Val_DoNotOverwriteUnreadSamps)
        self._scan_count_timeout = 2. * (N + 1) / self.scan_freq
        self.line_points = N

    def config_readout_sampling(self, N, trend_is_positive):
        self.scanner_ai_task.SetSampTimingType(PyDAQmx.DAQmx_Val_SampClk)  # set task timing to use a sampling clock

        if trend_is_positive:
            active_edge = PyDAQmx.DAQmx_Val_Falling
        else:
            active_edge = PyDAQmx.DAQmx_Val_Rising

            self.scanner_ai_task.CfgSampClkTiming(  # set up sample clock for task timing
                self._function_counter_out + 'InternalOutput',  # use these pulses as clock
                self.scan_freq,  # maximum expected clock frequency
                active_edge,
                PyDAQmx.DAQmx_Val_FiniteSamps,  # genarate sample on falling edge, generate finite number of samples
                N)  # samples to read

    def scan_line(self, Line):
        """Perform a line scan."""
        # check if length setup is correct, if not, adjust.
        line_points = Line.shape[1]
        if self.line_points != line_points:
            self.set_scan_length(line_points)

        # set up scanner ao channel
        self.write_scanner_ao(self.scanner_pos_to_volt(Line), length=self.line_points)
        # start tasks
        # CHK( dll.DAQmxStartTask(self.scanner_ao_task) )
        self.scanner_ao_task.StartTask()
        self.scanner_count_task.StartTask()
        self.scanner_pulse_task.StartTask()

        self.scanner_count_task.WaitUntilTaskDone(self._scan_count_timeout)

        self._scan_data = numpy.empty((self.line_points,), dtype=numpy.uint32)  # count data will be written here
        import time
        t0=time.time()
        self.scanner_count_task.ReadCounterU32(
            self.line_points,  # read number of "line_points" samples
            self._RWTimeout,
            self._scan_data,  # write into this array
            self.line_points,  # length of array to write into
            None, # ctypes.byref(n_read_samples)
            None
        )
        time.sleep(0.1)
        self.scanner_count_task.StopTask()
        self.scanner_pulse_task.StopTask()
        self.scanner_ao_task.StopTask()

        return self._scan_data * self.scan_freq / self.scan_duty_cycle  # normate to counts per second

    def scan_read_line(self, Line):
        """Perform a line scan and read."""
        # check if length setup is correct, if not, adjust.
        line_points = Line.shape[1]
        if Line[0][0] <= Line[0][-1]:
            trend_is_positive = True
        else:
            trend_is_positive = False

        if self.line_points != line_points:
            self.set_scan_length(line_points)

        self.config_readout_sampling(line_points, trend_is_positive)

        # set up scanner ao channel
        self.write_scanner_ao(self.scanner_pos_to_volt(Line), length=self.line_points)
        # start tasks
        # CHK( dll.DAQmxStartTask(self.scanner_ao_task) )
        self.scanner_ao_task.StartTask()
        self.scanner_ai_task.StartTask()
        self.scanner_count_task.StartTask()
        self.scanner_pulse_task.StartTask()

        self.scanner_count_task.WaitUntilTaskDone(self._scan_count_timeout)

        self._scan_data = numpy.empty((self.line_points,), dtype=numpy.uint32)  # count data will be written here
        self._read_pos_data = numpy.empty((3, self.line_points), dtype=numpy.float64)  # read position data will be written here
        self.scanner_count_task.ReadCounterU32(
            self.line_points,  # read number of "line_points" samples
            self._RWTimeout,
            self._scan_data,  # write into this array
            self.line_points,  # length of array to write into
            None, # ctypes.byref(n_read_samples),
            None
        )
        # readout voltages positions here
        self.scanner_ai_task.ReadAnalogF64(
            self.line_points,
            self._RWTimeout,
            PyDAQmx.DAQmx_Val_GroupByChannel,
            self._read_pos_data,
            3 * self.line_points,
            None, #ctypes.byref(n_read_samples),
            None
        )
        time.sleep(0.05)
        self.scanner_count_task.StopTask()
        self.scanner_pulse_task.StopTask()
        self.scanner_ao_task.StopTask()
        self.scanner_ai_task.StopTask()

        return self.scanner_volt_to_pos(self._read_pos_data), self._scan_data * self.scan_freq / self.scan_duty_cycle  # normalize to counts per second, return also voltages/positions

    # --------------- ODMR -------------------------------------
    def init_odmr(self, settling_time, count_time):
        """Start task for odmr measurement.
        Use pulse width measurment with two counters: Counter 1 generates pulses for gate of counter two,
        which measures width of pulses in terms of photon ticks."""

        self.function_state = 'odmr'
        self.odmr_pulse_task = PyDAQmx.Task()
        self.odmr_count_task = PyDAQmx.Task()
        # create pulses (clock) for odmr.
        self.odmr_freq = 1. / (count_time + settling_time)  # counting freqeuency
        duty_cycle = self.odmr_freq * count_time  # width of pulse (relative to count_interval) during which photons are counted
        self.odmr_duty_cycle = duty_cycle
        # now create pulse signal defined by f and duty_cycle
        self.odmr_pulse_task.CreateCOPulseChanFreq(
            self._function_counter_out,  # use this counter
            '',  # name
            PyDAQmx.DAQmx_Val_Hz,
            PyDAQmx.DAQmx_Val_High,
            0,  # units, idle state, inital delay
            self.odmr_freq,  # pulse frequency
            1 - duty_cycle
        )  # duty cycle of pulses
        # set up pulse width measurement in photon ticks, i.e. the width of each pulse generated by pulse_out_task is measured in photon ticks.
        self.odmr_count_task.CreateCIPulseWidthChan(
            self._function_counter_in,  # use this counter
            '',  # name
            0,  # expected minimum value
            self._MaxCounts * duty_cycle / self.odmr_freq,  # expected maximum value
            PyDAQmx.DAQmx_Val_Ticks,  # units of width measurement, here photon ticks
            PyDAQmx.DAQmx_Val_Falling, ''
        )  # start pulse width measurement on falling edge
        self.odmr_count_task.SetCIPulseWidthTerm(self._function_counter_in, self._function_counter_out + 'InternalOutput')
        # set the timebase for width measurement as self._photon_source
        self.odmr_count_task.SetCICtrTimebaseSrc(self._function_counter_in, self._photon_source)
        PyDAQmx.DAQmxConnectTerms(self._function_counter_out + 'InternalOutput', self._odmr_trig_channel, PyDAQmx.DAQmx_Val_DoNotInvertPolarity)
        self.odmr_points = -1  # reset variable, odmr length is checked and adjusted for each count operation
        return 0

    def stop_odmr(self):
        self.function_state = 'idle'
        PyDAQmx.DAQmxDisconnectTerms(self._function_counter_out + 'InternalOutput', self._odmr_trig_channel)

    def set_odmr_length(self, N):
        """Set length for line scan, i.e. number of clock pulses, number of count samples etc."""
        # set timing for odmr pulse and count task.
        # Contrary to scanner, here first frequency points should be set without trigger, i.e. first count sample can be used. After last count sample pulse goes low, i.e. frequency generator will get trigger and return to list pos 1
        self.odmr_pulse_task.CfgImplicitTiming(PyDAQmx.DAQmx_Val_FiniteSamps, N)
        self.odmr_count_task.CfgImplicitTiming(PyDAQmx.DAQmx_Val_FiniteSamps, N)
        # read samples from beginning of acquisition, do not overwrite
        self.odmr_count_task.SetReadRelativeTo(PyDAQmx.DAQmx_Val_CurrReadPos)
        self.odmr_count_task.SetReadOffset(0)
        self.odmr_count_task.SetReadOverWrite(PyDAQmx.DAQmx_Val_DoNotOverwriteUnreadSamps)
        self._odmr_count_timeout = 2. * (N) / self.odmr_freq
        self.odmr_points = N

    def count_odmr(self, points):
        """Count one run of odmr."""
        # check if length setup is correct, if not, adjust.
        if self.odmr_points != points:
            self.set_odmr_length(points)

        # start tasks
        self.odmr_count_task.StartTask()
        self.odmr_pulse_task.StartTask()

        self.odmr_count_task.WaitUntilTaskDone(self._odmr_count_timeout)

        self._odmr_data = numpy.empty((self.odmr_points,), dtype=numpy.uint32)  # count data will be written here
        self.odmr_count_task.ReadCounterU32(
            self.odmr_points,  # read number of "line_points" samples
            self._RWTimeout,
            self._odmr_data,  # write into this array
            self.odmr_points,  # length of array to write into
            None,
            None
        )  # number of samples which were read

        time.sleep(0.05)
        self.odmr_count_task.StopTask()
        self.odmr_pulse_task.StopTask()

        return self._odmr_data * self.odmr_freq / self.odmr_duty_cycle  # normate to counts per second
