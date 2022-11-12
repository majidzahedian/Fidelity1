# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import, division

import numpy
import numpy as np
import os
import pickle

from traits.api import *
from traitsui.api import *
from enable.api import ComponentEditor
from chaco.api import ArrayPlotData, Plot
from chaco.tools.api import PanTool, ZoomTool
import time
import traceback
import sys

try:
    from pi3diamond import pi3d
except:
    pi3d = None
    print('pi3diamond could not be imported')
import itertools
import operator
from qutip_enhanced import *
import pym8190a
from utility import GetSetItemsMixin

class Pulsed(HasTraits, GetSetItemsMixin):
    """This program is used to perform standard pulsed measurements.
    The used pulse sequence has to be provided by PulsePattern.
    From this sequence, it gets the starting index ('laser_on_indices') and length of the laser pulses.
    If use_as_tau is set for the sequence, it also gets the tau-values (x-axis).
    Alternatively, tau are set manually via 'tau_start' and 'tau_step'.
    It also sorts the 'laser_on_indices' based on tau, resp. depending on 'first_is_last' etc.
    It sets up the FastComtecCard (time-resolved photon counting), and starts the sequence (with DTG).
    The data is pulled from the FastComtecCard and saved in 'count_data' as [[data_pulse_tau_1], [data_pulse_tau_2], ...], i.e. it only pulls data from the pulse starting index, and length depending on pulse length and 'max_offset'.
    """

    sequence_name_list = List([''], transient=True)
    sequence_name = Str(label='Sequence')  # name of current sequence
    update_sequence_name_list_button = Button(label='Get Sequences')

    # device settings
    # mw_freq = Float(2870, desc='Frequency [MHz]', label='frequency [MHz]')
    # mw_power = Range(low=-100., high=25., value=-40, desc='Power [dBm]', label='power [dBm]', mode='text',
    #                  auto_set=False, enter_set=True)
    self.count_binnings = [1.,2.,3.,4.,5.,6.,7.,8.,9.,10.]
    self.bin_width = self.count_binnings[0]
    counter_safety = Int(1500,
                         desc='time [ns] subtracted from sequence length for photon counting. This is needed for some cards, as they will miss triggers otherwise.')

    #pulse analysis
    tau_auto = Bool(True, label='Tau defined in sequence?',
                    desc='If tau (i.e. mw-time for rabi) is defined in the sequence (this can be done by the PulsePattern block generator), this program can set up x-axis and pulse association automatically.')
    tau_start = Float(0, label='Tau start [micro-s]')  # If tau is not defined in sequence, the x-axis is made by the following parameters.
    tau_step = Float(0.1, label='Tau step [micro-s]')  # ...
    ignore_first = Bool(False, label='Ignore first laser?')
    ignore_last = Bool(False, label='Ignore last laser?')
    first_is_last = Bool(False, label='Make first laser last')
    alternating_sequence = Bool(False, label='Alternating sequence?',
                                desc='Splits results in two subsets of results, by alternating picking the pulses. If tau is set manually, uses the given parameters for both results individually')
    substract_results = Bool(False, label='Subtract results?', desc='Substracts results of an alternating sequence.')
    count_data = Array  # The data returned by the count card for each laser pulse.
    count_sum = Array  # Sum of all laser pulses.
    #pulses_enum = Property(trait=List, depends_on='count_data')     #
    pulses_list = List(['Sum'], desc='List 1...n, n=number of pulses')
    pulse_to_show = Str('Sum', label='Show pulse')
    analyze_time_window = Int(300, desc='time window for pulse analysis [ns]', label='analyze time window [ns]')
    auto_analyze_time = Button(label='Auto',
                               desc='This will search for the best analyze_time_window (best signal/noise)')
    normalize_time_window = Int(1000, label='Normalize time window[ns]',
                                desc='The time window at the end of the pulse which is used for normalization.')
    normalize_safety = Int(5, label='Normalize safety [ns]', desc='Last x [ns] of pulse are not used for normalization')
    max_offset = Int(600, label='Max offset [ns]',
                     desc='flourecence pulse is expected to start somewhere in start_offset + max_offset')
    start_offset = Int(0, label="Start offset [ns]",
                       desc='expected offset [ns] for the search for the fluorescence pulse, i.e. pulse is expected to start somewhere in start_offset + max_offset')
    do_normalization = Bool(True, label='Normalize?')
    average_normalize = Bool(True, label='Do averaged normalization?')
    calc_weight_func = Button(label='Calc weight function')
    weight_func = Array

    #refocus
    initial_refocus = Bool(False, label='Refocus/ODMR at start?')
    refocus_interval = Int(1, desc='Refocus before each n-th ODMR (0=None)', label='Refocus interval [n-th ODMR]')
    odmr_interval = Int(90, desc='Minimum time interval between ODMR (for sweep-like measurements: nuclear rabi)',
                        label='ODMR interval [s]')

    #run
    state = Enum('idle', 'run', transient=True)
    planned_sweeps = Range(low=1., high=1e10, value=1e9, desc='number of sweeps', label='Planned sweeps', mode='text',
                           auto_set=False, enter_set=True)
    sweeps = Float(value=0, desc='Elapsed Sweeps', label='Elapsed Sweeps', mode='text')
    get_data = Button(label='Pull and analyze')
    #analyze_button = Button(label='Analyze')
    file_path = Str('D:\data\Pulsed')
    file_name = Str('enter filename')
    file_notes = Str('', label='Notes')
    save_button = Button(label='Save')

    #fit
    fit_function = Enum('None', 'Rabi decay', 'CosineMultiDet', label='Function')
    do_fit = Button(label='Fit')
    centerfreq = Float(0, desc='Center frequency [MHz]', label='f0 [MHz]')
    rabi_contrast = Float(0, label='Contrast')
    rabi_period = Float(0, label='Period [ns]')
    rabi_decay = Float(0, label='Decay [micro-s]')
    rabi_offset = Float(0, label='Offset (tau=t+offset) [ns]')

    #plots
    results = Array()
    tau = Array()  #x-axis for results
    plot_results_data = Instance(ArrayPlotData, transient=True)
    plot_results = Instance(Plot, transient=True)
    colormap = ['blue', 'green', 'orange']
    plot_pulse_data = Instance(ArrayPlotData, transient=True)
    plot_pulse = Instance(Plot, transient=True)
    plot_snr_data = Instance(ArrayPlotData, transient=True)
    plot_snr = Instance(Plot, transient=True)
    plot_matrix_data = Instance(ArrayPlotData, transient=True)
    plot_matrix = Instance(Plot, transient=True)

    def __init__(self):
        HasTraits.__init__(self)
        self.count_binnings = [1.,2.,3.,4.,5.,6.,7.,8.,9.,10.]
        self.bin_width = self.count_binnings[0]

    def run(self, abort):
        try:
            self.init_run()
            pi3d.microwave.CW(pi3d.odmr.f_mhz * 1e6, pi3d.odmr.power)
            pi3d.mcas_dict[self.mcas.name] = self.mcas
            self.start_pulsed()
            time.sleep(1)
            pi3d.timetagger.pulsed.start()
            time.sleep(0.2)
            self.mcas.initialize()
            self.mcas.start_awgs()

            start_time = time.time()
            while self.sweeps < self.planned_sweeps:  # while True:
                if abort.is_set():
                    break
                self.do_refocusodmr(abort)
                if abort.is_set(): break
                time.sleep(2.)
                self.elapsed_time = time.time() - start_time
                self.sweeps = int(pi3d.timetagger.pulsed.getCounts()) # int(self.counter.getCounts())
            self._get_data_changed()
        except:
            abort.set()
            exc_type, exc_value, exc_tb = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_tb)
        finally:
            try:
                pi3d.timetagger.pulsed.stop()
                # self.counter.stop()
            except:
                pass
            pi3d.mcas_dict.stop_awgs()
            self.state = 'idle'

    def do_refocusodmr(self, abort):
        delta_t = time.time() - self.last_odmr
        if self.odmr_interval != 0 and (delta_t >= self.odmr_interval):
            pi3d.timetagger.pulsed.stop()
            if self.refocus_interval != 0 and self.odmr_count % self.refocus_interval == 0:
                pi3d.confocal.run_refocus()
            time.sleep(0.2)
            pi3d.odmr.run(abort)
            self.odmr_count += 1
            self.last_odmr = time.time()
            time.sleep(0.2)
            pi3d.timetagger.pulsed.start()
            time.sleep(0.2)
            self.mcas.initialize()
            self.mcas.start_awgs()

    def start_pulsed(self):
        n_bins=self.laser_length_bins + int(round(self.max_offset*1./self.bin_width))
        if self.alternating_sequence:
            n_histos=len(self.tau) + len(self.tau2)
        else:
            n_histos=len(self.tau)

        self.counter = pi3d.timetagger.init_counter(
            'pulsed',
            bw=int(self.bin_width * 1000.),
            n_bins=int(n_bins),
            n_histos=int(n_histos)
        )

    def init_run(self):
        self.state = 'run'
        self.sweeps = 0
        self.odmr_count = 0
        self.last_odmr = -np.inf
        self.get_parameters_from_sequence()

    def get_parameters_from_sequence(self):

        def consecutive_list_list(l):
            """
            example:

            l = [3,4, 5, 8, 10, 11]
            out = [[3,4], [5], [8], [10,11]]

            """
            return [map(operator.itemgetter(1), g) for k, g in itertools.groupby(enumerate(l), lambda (i, x): i-x)]

        asd = self.mcas.sequences['2g']
        wsl = asd[2].wavestep_list
        self.sequence_length = asd[1].repeated_length_smpl
        self.sequence_length_ns = int(asd[1].repeated_length_mus * 1e3)

        wavestep_length_smpl_list = np.array([step.length_smpl for step in wsl])
        wavestep_length_mus_list = np.array([step.length_mus for step in wsl])

        laser_on_wavestep_indices = np.array(consecutive_list_list([i for i, step in enumerate(wsl) if step.smpl_marker])) #during which wavsteps is the laser on
        self.laser_length_ns_list = pym8190a.elements.round_length_mus_full_sample([1e3*sum(wavestep_length_mus_list[i]) for i in laser_on_wavestep_indices]) #list of duration of all laserpulses. all should be equal
        self.laser_on_indices_raw = np.cumsum([[0.0] + list(wavestep_length_smpl_list)])[laser_on_wavestep_indices[:,0]] # after which samples does laser start
        tau_wavestep_indices = [i for i, step in enumerate(wsl) if '_pulsed_tau_' in step.name] # Which wavesteps are used as tau.
        tau_list_length_mus = np.array([sum(wavestep_length_mus_list[i]) for i in consecutive_list_list(tau_wavestep_indices)]) #How long (mus) are these wavesteps combined per incident
        self.tau_raw = pym8190a.elements.round_length_mus_full_sample(tau_list_length_mus)
        if self.sequence_length_ns - self.laser_on_indices_raw[-1]  / (pym8190a.elements.__SAMPLE_FREQUENCY__*1e-3) + \
                self.laser_length_ns_list[-1] < self.max_offset + self.start_offset + self.counter_safety:
            raise RuntimeError(
                'duration after last laser must be at least start_offset (%s ns) + least max_offset (%s ns) + counter_safety (%s ns) = %s ns' % (
                    self.start_offset, self.max_offset, self.counter_safety,
                    self.start_offset + self.max_offset + self.counter_safety))
        if self.check_laser_length() is False:
            raise RuntimeError('Laser pulses with non-equal length in sequence.')
        self._get_tau()

    def _get_tau(self):
        tau = self.tau_raw[:]
        self.laser_on_indices = self.laser_on_indices_raw[:]
        if not self.tau_auto:  #if tau are not defined in sequence, create x-axis here.
            if self.first_is_last:  #if first laser is (logically) last point in sequence.
                temp = self.laser_on_indices[0]
                self.laser_on_indices.append(temp)  # put first laser_on_index at end of sequence
                self.laser_on_indices = self.laser_on_indices[1:]  # remove first laser_on_index
                #temp = self.laser_length_ns[0]
                #self.laser_length_ns.append(temp)
                #self.laser_length_ns = self.laser_length_ns[1:]
            tau = numpy.arange(len(self.laser_on_indices))  # tau = [0,1,2,3..], multiply tau_step, add tau_start later
            if self.alternating_sequence:  #i alternating sequence, tau = [0,0,1,1,2,2,...], divide into two lists later
                for i in range(len(tau)):
                    if i % 2 == 1:
                        tau[i] = tau[i - 1]
            tau *= self.tau_step
            tau += self.tau_start
        if self.ignore_last:  # ignore last laser pulse
            self.laser_on_indices = self.laser_on_indices[:-1]
            #self.laser_length_ns = self.laser_length_ns[:-1]
            tau = tau[:-1]
        if self.ignore_first:  # ignore first laser pulse
            self.laser_on_indices = self.laser_on_indices[1:]
            #self.laser_length_ns = self.laser_length_ns[1:]
            tau = tau[1:]
        if self.tau_auto and not self.alternating_sequence:  # (Col(D)/0.00272-1)/0.0039083sort x-axis and corresponding pulses
            tau_sorted, self.laser_on_indices = pi3d.misc.sort_lists(tau, self.laser_on_indices)
            #tau_sorted, self.laser_length_ns = pi3d.misc.sort_lists(tau, self.laser_length_ns)
            tau = tau_sorted
        if self.alternating_sequence:  #divide into two lists
            tau1 = tau[::2]
            lon1 = self.laser_on_indices[::2]
            #llns1 = self.laser_length_ns[::2]
            tau1_sorted, lon1_sorted = pi3d.misc.sort_lists(tau1, lon1)
            #tau1_sorted, llns1_sorted = pi3d.misc.sort_lists(tau1, llns1)
            tau2 = tau[1::2]
            lon2 = self.laser_on_indices[1::2]
            #llns2 = self.laser_length_ns[1::2]
            tau2_sorted, lon2_sorted = pi3d.misc.sort_lists(tau2, lon2)
            #tau2_sorted, llns2_sorted = pi3d.misc.sort_lists(tau2, llns2)
            self.laser_on_indices = []
            #self.laser_length_ns = []
            for i in range(len(tau1_sorted)):
                self.laser_on_indices.append(lon1_sorted[i])
                self.laser_on_indices.append(lon2_sorted[i])
                #self.laser_length_ns.append(llns1_sorted[i])
                #self.laser_length_ns.append(llns2_sorted[i])
            self.tau = tau1_sorted
            self.tau2 = tau2_sorted
            self.results2 = numpy.zeros(len(self.tau2))
        else:
            self.tau = tau
        if self.alternating_sequence and self.substract_results:
            for i in range(len(self.tau)):
                if self.tau[i] != self.tau2[i]:
                    print('Warning: Subtracting results of alternating sequence with different tau!!!')
                    break
        self.results = numpy.zeros(len(self.tau))

    def _get_data_changed(self):

        self._get_tau()
        self._init_plot_results()

        self.laser_length_bins = int(round(self.laser_length_ns * 1. / self.bin_width))
        # length = self.laser_length_bins + int(round(self.max_offset * 1. / self.bin_width))
        bin_index = []
        for i, index in enumerate(self.laser_on_indices):
            bin_index.append(int(round((index  / (pym8190a.elements.__SAMPLE_FREQUENCY__*1e-3) + self.start_offset) * 1. / self.bin_width)))
        new_data = pi3d.timetagger.pulsed.getData() #self.counter.getData()
        self.count_data = numpy.array(new_data)
        self.pulses_list = list(numpy.arange(len(self.count_data)))
        self.count_sum = numpy.sum(self.count_data, 0)  #sum of all pulses
        overlap = []
        for i in range(len(self.count_sum) - self.laser_length_bins):
            overlap.append(numpy.sum(self.count_sum[i:i + self.laser_length_bins]))
        self.offset = numpy.array(overlap).argmax()
        average = max(overlap) * 1. / self.laser_length_bins
        average_array = numpy.zeros(self.offset)
        average_array = numpy.append(average_array, numpy.ones(self.laser_length_bins) * average)
        average_array = numpy.append(average_array, numpy.zeros(len(self.count_sum) - len(average_array)))
        self.plot_pulse_data = self._plot_pulse_data_default()
        self.plot_pulse_data.set_data('overlap', average_array)
        self.plot_pulse = self._plot_pulse_default()
        self.plot_pulse.request_redraw()
        self.analyze()

    def _analyze_button_changed(self):
        self.analyze()

    def analyze(self):
        offset = self.offset
        self.results_all = []
        signal = []
        normalization = []
        for i, pulse in enumerate(self.count_data):
            signal.append(1. * sum(pulse[offset:int(round(offset + self.analyze_time_window * 1. / self.bin_width))]))
            normalization.append(self.analyze_time_window * 1. / self.normalize_time_window *
                                 sum(pulse[offset + self.laser_length_bins - int(
                                     round((self.normalize_time_window - self.normalize_safety) * 1. / self.bin_width)):
                                 offset + self.laser_length_bins - int(
                                     round(self.normalize_safety * 1. / self.bin_width))]))
            if self.do_normalization and not self.average_normalize:
                self.results_all.append(signal[-1] * 1. / normalization[-1])
        self.normalization = normalization
        if self.do_normalization and self.average_normalize:
            self.normalization_average = sum(self.normalization) * 1. / len(self.normalization)
            self.results_all = numpy.array(signal) / self.normalization_average
        if not self.do_normalization:
            self.results_all = numpy.array(signal)
        if self.alternating_sequence:
            results = numpy.array(self.results_all[::2])
            results2 = numpy.array(self.results_all[1::2])
            if self.substract_results:
                self.results = results2 - results
            else:
                self.results = results
                self.results2 = results2
        else:
            self.results = numpy.array(self.results_all)

        self.plot_results_data.set_data('y', self.results)
        if self.alternating_sequence:
            self.plot_results_data.set_data('y2', self.results2)
        self.plot_results.request_redraw()

    def _auto_analyze_time_changed(self):
        if len(self.count_data) == 0:
            print('no data')
            return
        offset = self.offset
        pulse_bright = self.count_data[self.results_all.argmax()]
        pulse_dark = self.count_data[self.results_all.argmin()]
        snr = []
        signal_bright = numpy.cumsum(pulse_bright[offset:offset + self.laser_length_bins], dtype=int)
        signal_dark = numpy.cumsum(pulse_dark[offset:offset + self.laser_length_bins], dtype=int)
        for i in range(len(signal_bright)):
            snr.append(abs(signal_bright[i] - signal_dark[i]) * 1. / ((signal_bright[i] + signal_dark[i]) ** 0.5))
        self.snr = numpy.array(snr)
        self.analyze_time_window = int(self.snr.argmax() * self.bin_width)
        self.plot_snr_data.set_data('x', numpy.arange(min(self.snr.argmax() * 3, self.laser_length_bins)))
        self.plot_snr_data.set_data('snr', self.snr[:min(self.snr.argmax() * 3, self.laser_length_bins)])
        self.plot_snr.request_redraw()

    def _calc_weight_func_changed(self):
        if len(self.count_data) == 0:
            print('no data')
            return
        offset = self.offset

        def smoothing(data, binning=2):
            new_data = []
            for i in range(len(data)):
                new_data.append(sum(
                    data[max(0, i - 3):min(i + 3, len(data))] * 1. / len(data[max(0, i - 3):min(i + 3, len(data))])))
            return numpy.array(new_data)

        #calculate with signal from brightest and darkest pulse
        pulse_bright = smoothing(
            numpy.array(self.count_data[self.results_all.argmax()][offset:offset + self.laser_length_bins], dtype=int))
        pulse_dark = smoothing(
            numpy.array(self.count_data[self.results_all.argmin()][offset:offset + self.laser_length_bins], dtype=int))
        weight_func = []
        for i in range(len(pulse_bright)):
            #formula from lukins repetitive readout paper (Science 2009)
            weight_func.append(abs(pulse_bright[i] - pulse_dark[i]) * 1. / (pulse_bright[i] + pulse_dark[i]))
        self.weight_func = numpy.array(weight_func)
        self.snr_weight = numpy.cumsum((pulse_bright - pulse_dark) ** 2 / (pulse_bright + pulse_dark), dtype=int) ** 0.5
        self._auto_analyze_time_changed()
        self.plot_snr_data.set_data('snr_weight', self.snr_weight[:min(self.snr.argmax() * 3, self.laser_length_bins)])
        self.plot_snr.plot(('x', 'snr_weight'), style='line', color='green', name='snr_weight')
        self.plot_snr.request_redraw()

    def check_laser_length(self):
        laser_length = self.laser_length_ns_list[0]
        for length in self.laser_length_ns_list:
            if length != laser_length:
                return False
        self.laser_length_ns = laser_length
        self.laser_length_bins = int(round(self.laser_length_ns * 1. / self.bin_width))

    def _do_fit_changed(self):
        if self.fit_function == "None":
            print('select fit function')
            return
        if self.fit_function != "None":
            if self.fit_function == 'Rabi decay':
                x = self.tau
                y = self.results
                try:
                    m = lmfit_models.CosineModel()
                    p = m.guess(data=y, x=x)
                    r = m.fit(data=y, params=p, x=x)
                except:
                    print('Fitting failed.')
                    return None
                # if p[0] < 0:
                #     p[0] = -p[0]
                #     p[2] = ( ( p[2] / p[1] + 0.5 ) % 1 ) * p[1]
                #     p = pi3d.Fit.Fit(x, y, pi3d.Fit.CosinusNoOffset, p)
                # p = (p[0], p[1], p[2], y_offset)
                # p = pi3d.Fit.Fit(x, self.results, pi3d.Fit.Cosinus, p)
                # while (p[2] > 0.5 * p[1]):
                #     p[2] -= p[1]
                # p = pi3d.Fit.Fit(x, self.results, pi3d.Fit.Cosinus, p)
                # p = list(p)
                # p.append(10 * max(x))
                # p = pi3d.Fit.Fit(x, self.results, pi3d.Fit.Cosinus_dec, p)
                # self.rabi_contrast = 200 * p[0] / p[3]
                # self.rabi_period = p[1] * 1e3
                # self.rabi_decay = p[4]
                # self.rabi_offset = p[2] * 1e3
                # function = pi3d.Fit.Cosinus_dec(*p)(x)
                # self.plot_results_data.set_data('fit', function)
            if self.fit_function == 'CosineMultiDet':
                x = self.tau
                y = self.results
                try:
                    min_chisqr = np.inf
                    for method in ['leastsq', 'powell', 'cg', 'slsqp']:
                        m = lmfit_models.CosineMultiDetModel([pi3d.tt.get_f('14n_hf'), pi3d.tt.get_f('13c414_hf'), pi3d.tt.get_f('13c90_hf')])
                        p = m.guess(data=y, x=x)
                        p['f0'].vary=getattr(self, 'vary_f0', False)
                        new_r = m.fit(data=y, params=p, x=x, method=method)
                        if new_r.chisqr < min_chisqr:
                            min_chisqr = new_r.chisqr
                            r = new_r
                except:
                    print('Fitting failed.')
                    return None
            self.rabi_contrast = 200 * r.params['amplitude'].value / r.params['c'].value
            self.rabi_period = r.params['T'].value * 1e3
            self.rabi_decay = r.params['t2'].value
            self.rabi_offset = r.params['x0'].value * 1e3
            x_fit_plot = np.linspace(x[0], x[-1], 1000)
            y_fit_plot = r.eval(params=r.params, x=x_fit_plot)
            self.plot_results_data.set_data('x_fit_plot', x_fit_plot)
            self.plot_results_data.set_data('y_fit_plot', y_fit_plot)
            self.plot_results.plot(('x_fit_plot', 'y_fit_plot'), style='line', color='red')
            self.plot_results.request_redraw()

    def _get_pulses_enum(self):
        pulses = ['Sum']
        for i in range(len(self.count_data)):
            pulses.append(i)
        return pulses

    def _pulse_to_show_changed(self):
        if self.pulse_to_show == 'Sum':
            pulse = self.count_sum
            laser_on = numpy.zeros(self.offset)
            laser_on = numpy.append(laser_on, numpy.ones(self.laser_length_bins) * max(pulse))
            laser_on = numpy.append(laser_on, numpy.zeros(len(self.count_sum) - len(laser_on)))
        else:
            pulse_index = int(self.pulse_to_show)
            pulse = self.count_data[pulse_index]
            laser_on = numpy.zeros(self.offset)
            laser_on = numpy.append(laser_on, numpy.ones(self.laser_length_bins) * max(pulse))
            laser_on = numpy.append(laser_on, numpy.zeros(len(pulse) - len(laser_on)))
        self.plot_pulse_data = self._plot_pulse_data_default()
        self.plot_pulse_data.set_data('overlap', laser_on)
        self.plot_pulse_data.set_data('y', pulse)
        self.plot_pulse = self._plot_pulse_default()
        self.plot_pulse.request_redraw()

    def _update_sequence_name_list_button_fired(self):
        self.update_sequence_name_list()

    def update_sequence_name_list(self):
        self.sequence_name_list = list(
            set(pi3d.awg2g.ch[1].sequence_name_list) & set(pi3d.awg2g.ch[2].sequence_name_list))

    def _tau_default(self):
        return numpy.array([0, 1])

    def _results_default(self):
        return numpy.zeros(self.tau.shape)

    def _plot_results_data_default(self):
        return ArrayPlotData(x=self.tau, y=self.results)

    def _plot_results_default(self):
        plot = Plot(self.plot_results_data, padding_left=60, padding_top=10, padding_right=10, padding_bottom=30)
        plot.plot(('x', 'y'), style='line', color='blue', name='result')
        plot.index_axis.title = 'tau [micro-s]'
        plot.value_axis.title = 'Intensity [a.u.]'
        plot.tools.append(PanTool(plot))
        plot.overlays.append(ZoomTool(plot))
        return plot

    def _init_plot_results(self):
        plots = []
        for plot in self.plot_results.plots:
            plots.append(plot)
        for plot in plots:
            self.plot_results.delplot(plot)
            self.plot_results.datasources.clear()
        data = self.plot_results_data.list_data()
        for d in data:
            self.plot_results_data.del_data(d)
        self.plot_results_data.set_data('x', self.tau)
        self.plot_results_data.set_data('y', self.results)
        if self.alternating_sequence and not self.substract_results:
            self.plot_results_data.set_data('x2', self.tau2)
            self.plot_results_data.set_data('y2', self.results2)
            self.plot_results.plot(('x2', 'y2'), style='line', color='green', name='result2')
        self.plot_results.plot(('x', 'y'), style='line', color='blue', name='result')
        self.plot_results.request_redraw()

    def _count_sum_default(self):
        return numpy.zeros(4000)

    def _plot_pulse_data_default(self):
        return ArrayPlotData(x=numpy.arange(len(self.count_sum)), y=self.count_sum, overlap=self.count_sum)

    def _plot_pulse_default(self):
        plot = Plot(self.plot_pulse_data, padding_left=60, padding_top=10, padding_right=10, padding_bottom=30)
        plot.plot(('x', 'y'), style='line', color='blue', name='pulse')
        plot.plot(('x', 'overlap'), style='line', color='red', name='overlap')
        plot.index_axis.title = 'bin'
        plot.value_axis.title = 'number of photons'
        plot.tools.append(PanTool(plot))
        plot.overlays.append(ZoomTool(plot))
        return plot

    def _init_plot_pulse(self):
        plots = []
        for plot in self.plot_pulse:
            plots.append(plot)
        for plot in plots:
            self.plot_pulse.delplot(plot)
            self.plot_pulse.datasources.clear()
        data = self.plot_pulse_data.list_data()
        for d in data:
            self.plot_pulse_data.del_data(d)
        self.plot_pulse_data.set_data('x', numpy.arange(len(self.count_sum)))
        self.plot_pulse_data.set_data('y', self.count_sum)
        self.plot_pulse_data.set_data('overlap', self.count_sum)
        self.plot_pulse.plot(('x', 'y'), style='line', color='blue', name='pulse')
        self.plot_pulse.plot(('x', 'overlap'), style='line', color='red', name='overlap')
        self.plot_pulse.request_redraw()

    def _plot_snr_data_default(self):
        x_temp = numpy.array(range(1000))
        return ArrayPlotData(x=x_temp, snr=numpy.zeros(x_temp.shape))

    def _plot_snr_default(self):
        plot = Plot(self.plot_snr_data, padding_left=60, padding_top=10, padding_right=10, padding_bottom=30)
        plot.plot(('x', 'snr'), style='line', color='blue', name='snr')
        plot.index_axis.title = 'bin'
        plot.value_axis.title = 'SNR'
        plot.tools.append(PanTool(plot))
        plot.overlays.append(ZoomTool(plot))
        return plot

    def _save_button_changed(self):
        #os.chdir(self.file_path)
        if self.file_name in os.listdir(self.file_path):
            print('File already exists! Data NOT saved!')
            print('Choose other filename!')
            return
        os.mkdir(self.file_path + '/' + self.file_name)
        filename = self.file_path + '/' + self.file_name + '/' + self.file_name
        file_notes = open(filename + '_notes.txt', 'w')
        file_notes.write(self.file_notes)
        file_notes.close()
        #try:
        #    fil = open(filename+'.pyd','w')
        #    pickle.dump(self.__getstate__(), fil)
        #    fil.close()
        #except Exception:
        #    print 'failed to dump pulsed'
        file = open(filename + '.pyd', 'w')
        data = {}
        state = SingletonHasTraits.__getstate__(self)
        for key in state:
            if key == 'tau': data['tau [mu-s]'] = self.tau
            if key == 'results': data['result'] = self.results
            if self.alternating_sequence:
                if key == 'tau2': data['tau2 [mu-s]'] = self.tau2
                if key == 'results2': data['result2'] = self.results2
            if key == 'sequence_name': data['sequence_name'] = self.sequence_name
            if key == 'mw_freq': data['mw_freq [MHz]'] = pi3d.odmr.f_mhz
            if key == 'mw_power': data['mw_power [dB]'] = pi3d.odmr.power
            if key == 'bin_width': data['bin_width of counter [ns]'] = self.bin_width
            if key == 'count_data': data['raw count_data'] = self.count_data
            if key == 'sweeps': data['number of sweeps'] = self.sweeps
        pickle.dump(data, file)
        file.close()

        file = open(filename + '_result.dat', 'w')
        if self.alternating_sequence:
            file.write('tau [micro-s] \t result \t tau2 [micro-s] \t result2\n')
            for i, x in enumerate(self.tau):
                file.write('%s\t%s\t%s\t%s\n' % (self.tau[i], self.results[i], self.tau2[i], self.results2[i]))
            file.close()
        else:
            file.write('tau [micro-s] \t result \n')
            for i, x in enumerate(self.tau):
                file.write('%s \t %s \n' % (x, self.results[i]))
            file.close()
        #file = open(path + 'pulses.dat', 'w')
        self.save_figure(self.plot_results, filename + '_plot.png')
        print('saved pulsed ' + self.file_name)

    # def _state_changed(self):
    #     self.stop()
    #     if self.state == 'run' or self.state == 'continue':
    #         self.start()
    #
    # def __getstate__(self):
    #     """Returns current state of a selection of traits.
    #     Overwritten HasTraits.
    #     """
    #     state = SingletonHasTraits.__getstate__(self)
    #     for key in ['is_tracker_client', 'pause_order', 'stop_request', 'thread', 'AWG', 'pp', 'sequence',
    #                 'sequence_length', 'sequence_length_ns']:
    #         if key in state:
    #             del state[key]
    #     return state

    def save_line_plot(self, filename=None):
        if filename is None:
            filename = pi3d.get_filename() + '_odmr_line_plot.png'
        else:
            filename = pi3d.get_filename(filename)
        self.save_figure(self.plot_results, filename)

    traits_view = View(
        Tabbed(
            VGroup(
                HGroup(
                    VGroup(
                        # Item('mw_freq', width=-70),
                        # Item('mw_power', width=-70),
                        Item('bin_width', width=-70),
                    ),

                    VGroup(
                        Item('update_sequence_name_list_button', show_label=False),
                        Item('sequence_name', editor=EnumEditor(name='sequence_name_list'), width=-100),
                    ),
                    VGroup(
                        Item('ignore_first'),
                        Item('ignore_last'),
                        Item('alternating_sequence'),
                        Item('substract_results', visible_when="alternating_sequence"),
                    ),
                    VGroup(
                        Item('tau_auto'),
                        Item('tau_start', visible_when="not tau_auto"),
                        Item('tau_step', visible_when="not tau_auto"),
                        Item('first_is_last', visible_when="not tau_auto"),
                    ),
                ),
                HGroup(
                    Item('refocus_interval', width=-40),
                    Item('odmr_interval', width=-40),
                ),
                HGroup(
                    Item('state', style='custom', show_label=False, ),
                    Item('planned_sweeps', editor=TextEditor(evaluate=float, format_func=lambda x: '%.3e' % x),
                         width=40),
                    # Item('expected_duration', style='readonly', width=40,
                    #      editor=TextEditor(evaluate=float, format_func=lambda x: '%s' % (round(x, 0)))),
                    Item('sweeps', style='readonly', width=40,
                         editor=TextEditor(evaluate=float, format_func=lambda x: '%.2e' % x)),
                ),
                HGroup(
                    Item('get_data', show_label=False),
                    Item('file_name'),
                    Item('file_notes'),
                    Item('save_button', show_label=False),
                ),
                HGroup(
                    Item('do_fit', show_label=False),
                    Item('fit_function'),
                    Item('centerfreq', width=35, visible_when='fit_function=="CosineMultiDet"',
                         editor=TextEditor(evaluate=float, format_func=lambda x: '%s' % (round(x, 1)))),
                    Item('rabi_period', style='readonly', width=35, visible_when='fit_function in ["Rabi decay", "CosineMultiDet"]',
                         editor=TextEditor(evaluate=float, format_func=lambda x: '%s' % (round(x, 2)))),
                    Item('rabi_contrast', style='readonly', width=35, visible_when='fit_function in ["Rabi decay", "CosineMultiDet"]',
                         editor=TextEditor(evaluate=float, format_func=lambda x: '%s' % (round(x, 1)))),
                    Item('rabi_decay', style='readonly', width=35, visible_when='fit_function in ["Rabi decay", "CosineMultiDet"]',
                         editor=TextEditor(evaluate=float, format_func=lambda x: '%s' % (round(x, 2)))),
                    Item('rabi_offset', style='readonly', width=35, visible_when='fit_function in ["Rabi decay", "CosineMultiDet"]',
                         editor=TextEditor(evaluate=float, format_func=lambda x: '%s' % (round(x, 1)))),
                ),
                Item('plot_results', editor=ComponentEditor(), show_label=False, width=400, height=300, resizable=True),
                label='Sequence'
            ),
        ),
        VGroup(
            HGroup(
                VGroup(
                    Item('analyze_time_window'),
                    Item('auto_analyze_time'),
                    Item('average_normalize'),
                    Item('do_normalization'),
                    Item('normalize_time_window'),
                    Item('normalize_safety'),
                    Item('start_offset'),
                    Item('max_offset'),
                ),
                Item('plot_snr', editor=ComponentEditor(), show_label=False, width=200, height=250, resizable=True),
            ),
            HGroup(
                VGroup(
                    Item(label='Show pulse:'),
                    Item('pulse_to_show', show_label=False, editor=EnumEditor(name='pulses_list')),
                ),
                Item('plot_pulse', editor=ComponentEditor(), show_label=False, width=400, height=300, resizable=True),
            ),
            label='Analysis'
        ),
        title='Pulsed',
        width=800,
        height=600,
        resizable=True,
        buttons=['OK'],
        kind='live',
        id='PulsedView'
    )
