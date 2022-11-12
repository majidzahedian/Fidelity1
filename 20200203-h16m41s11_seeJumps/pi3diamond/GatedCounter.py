from __future__ import print_function, absolute_import, division
__metaclass__ = type

import misc
import traceback
import time
import sys
import os
import threading
import multi_channel_awg_seq as MCAS
import Analysis
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
#import pyqtgraph as pg
from matplotlib.backends.backend_qt5 import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import PyQt5.QtCore
import numpy as np
import logging
import zmq

# from qutip_enhanced import *
import qutip_enhanced.qtgui.gui_helpers
from numbers import Number
from pi3diamond import pi3d

class GatedCounter(qutip_enhanced.qtgui.gui_helpers.WithQt):

    def __init__(self, title='Gated Counter', parent=None, gui=True, **kwargs):
        super(GatedCounter, self).__init__(title=title, parent=parent, gui=gui, QtGuiClass=GatedCounterQt)
        self.readout_interval = 1e6 # DO NOT SET BELOW 0.2, plotting runs in main loop, and takes around 0.15, which then makes the console annoying to work with.
        self.readout_duration = 10e6
        self.trace = Analysis.Trace()
        if title is not None:
            self.update_window_title(title)

    readout_interval = misc.ret_property_typecheck('readout_interval', Number)

    progress = misc.ret_property_typecheck('progress', int)

    @property
    def points(self):
        if hasattr(self, 'n_values'):
            return self.n_values/sum([step[3] for step in self.trace.analyze_sequence])
        return self._points

    @points.setter
    def points(self, val):
        self._points = misc.check_type(val, 'points', int)

    @property
    def n_values(self):
        return self._n_values

    @n_values.setter
    def n_values(self, val):
        self._n_values = val


    def set_n_values(self, mcas, analyze_sequence=None):
        analyze_sequence = self.trace.analyze_sequence if analyze_sequence is None else analyze_sequence
        self.n_values = int(self.readout_duration / mcas.length_mus * sum([step[3] for step in analyze_sequence]))

    def read_trace(self):
        self.gated_counter_data = pi3d.timetagger.gated_counter_countbetweenmarkers.getData().astype(np.int16)
        self.set_progress()
        self.trace_rep = Analysis.TraceRep(trace=self.gated_counter_data[:self.progress], analyze_sequence=self.trace.analyze_sequence, number_of_simultaneous_measurements=self.trace.number_of_simultaneous_measurements)
        self.trace.trace = np.array(self.trace_rep.df.groupby(['run', 'sm', 'step', 'memory']).agg({'n': np.sum}).reset_index().n)
        self.update_plot_data()

    def count(self, abort, ch_dict, turn_off_awgs=True):
        if hasattr(self, '_gui'):
            self.gui.clear_plot(len(self.trace.analyze_sequence))
        try:
            self.set_counter()
            MCAS.start_awgs(pi3d.awgs, ch_dict=ch_dict)
            self.progress = 0
            while True:
                if abort.is_set():
                    break
                ready = pi3d.timetagger.gated_counter_countbetweenmarkers.ready()
                self.read_trace()
                self.update_plot()
                if ready:
                    break
                else:
                    time.sleep(self.readout_interval/1e6)
        except Exception as e:
            abort.set()
            exc_type, exc_value, exc_tb = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_tb)
        finally:
            if turn_off_awgs:
                pi3d.mcas_dict.stop_awgs()
            pi3d.timetagger.gated_counter_countbetweenmarkers.stop()
            self.state = 'idle'

    def set_counter(self):

        def f():
            nlp_per_point = sum([step[3] for step in self.trace.analyze_sequence])
            pi3d.timetagger.init_counter(
                'gated_counter_countbetweenmarkers',
                n_values=self.n_values - self.n_values % nlp_per_point if hasattr(self, '_n_values') else nlp_per_point * self.points,
            )
        for i in range(2):
            t = threading.Thread(target=f)
            t.start()
            t.join(2)
            if t.is_alive():
                logging.getLogger().info('Setting up counter failed!')
                logging.getLogger().info('Trying to restart timetagger..')
                pi3d.restart_timetagger()
            else:
                break
        else:
            raise Exception('Error: timeout.')

    def set_progress(self):
        if pi3d.timetagger.gated_counter_countbetweenmarkers.ready():
            self.progress = int(len(self.gated_counter_data))
        else:
            self.progress = int(len(self.gated_counter_data) - np.argmax(self.gated_counter_data[::-1] != 0) - 1)

    def init_plot(self):
        if hasattr(self, '_gui'):
            self.gui.init_plot()

    def update_plot_data(self):
        self.effective_subtrace_list = []
        self.hist_list = []
        df = self.trace.df_extended()
        for idx, step in enumerate(self.trace.analyze_sequence):
            try:
                self.effective_subtrace_list.append(df.loc[df.step == idx, range(step[5])].astype(int).values)
                self.hist_list.append([])
                estl = self.effective_subtrace_list[-1].T
                if estl.shape[0] == 2:
                    self.hist_list[-1].append(self.trace.hist(estl[0, :] - estl[1, :]))
                else:
                    for t in estl:
                        self.hist_list[-1].append(self.trace.hist(t))
            except:
                pass

    def update_plot(self):
        if hasattr(self, '_gui'):
            self.update_plot_data()
            self.gui.update_plot(self.effective_subtrace_list, self.hist_list)

class GatedCounterQt(qutip_enhanced.qtgui.gui_helpers.QtGuiClass):
    def __init__(self, parent=None, no_qt=None):
        super(GatedCounterQt, self).__init__(parent=parent, no_qt=no_qt, ui_filepath=os.path.join(pi3d.app_dir, 'qtgui/gated_counter.ui'))

    clear_plot_signal = PyQt5.QtCore.pyqtSignal(int)
    update_plot_signal = PyQt5.QtCore.pyqtSignal(list, list)

    def clear_plot(self, number_of_subtraces):
        self.clear_plot_signal.emit(number_of_subtraces)

    def clear_plot_signal_emitted(self, number_of_subtraces):
        try:
            self.fig.clear()
            self.axes = self.fig.subplots(number_of_subtraces, 2, squeeze=False)
            self.canvas.draw()
        except:
            exc_type, exc_value, exc_tb = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_tb)

    def clear_signal_emitted(self):
        pass

    def init_gui(self):
        super(GatedCounterQt, self).init_gui()
        for name in [
            'clear_plot',
            'update_plot',
        ]:
            getattr(getattr(self, "{}_signal".format(name)), 'connect')(getattr(self, "{}_signal_emitted".format(name)))

        # Figure
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)

        self.ax = self.fig.add_subplot(111)
        self.canvas.draw()
        self.plot_layout.addWidget(self.canvas, 1, 1, 20, 20)
        self.toolbar_layout.addWidget(self.toolbar, 21, 1, 1, 20)

    def update_plot(self, effective_subtrace_list, hist_list):
        self.update_plot_signal.emit(effective_subtrace_list, hist_list)

    def update_plot_signal_emitted(self, effective_subtrace_list, hist_list):
        try:
            t0 = time.time()
            for nst, st in enumerate(effective_subtrace_list):
                self.axes[nst, 0].cla()
                self.axes[nst, 0].plot(st)
                self.axes[nst, 1].cla()
                if len(hist_list[nst]) > 0:
                    for h in hist_list[nst]:
                        self.axes[nst, 1].plot(h['bin_edges'], h['hist'])
                else:
                    self.axes[nst, 1].text(
                        0.5,
                        0.5,
                        'HISTOGRAM FAILED',
                        horizontalalignment='center',
                        verticalalignment='center',
                        transform=self.axes[nst, 1].transAxes,
                        color='r',
                        fontsize=30
                    )
                self.fig.tight_layout()
                self.fig.canvas.draw_idle()
        except:
            exc_type, exc_value, exc_tb = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_tb)


class GatedCounterQtGraph(qutip_enhanced.qtgui.gui_helpers.QtGuiClass):
    """
    The same as GatedCounterQt but with Pyqtgraph'''
    """
    def __init__(self, parent=None, no_qt=None):
        super(GatedCounterQt, self).__init__(parent=parent, no_qt=no_qt, ui_filepath=os.path.join(pi3d.app_dir, 'qtgui/gated_counter.ui'))

    clear_plot_signal = PyQt5.QtCore.pyqtSignal(int)
    update_plot_signal = PyQt5.QtCore.pyqtSignal(list, list)

    def clear_plot(self, number_of_subtraces):
        self.clear_plot_signal.emit(number_of_subtraces)

    def clear_plot_signal_emitted(self, number_of_subtraces):
        try:
            self.fig.clear()
            self.axes = self.fig.subplots(number_of_subtraces, 2, squeeze=False)
            self.canvas.draw()
        except:
            exc_type, exc_value, exc_tb = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_tb)

    def clear_signal_emitted(self):
        pass

    def init_gui(self):
        super(GatedCounterQt, self).init_gui()
        for name in [
            'clear_plot',
            'update_plot',
        ]:
            getattr(getattr(self, "{}_signal".format(name)), 'connect')(getattr(self, "{}_signal_emitted".format(name)))

        # Figure




        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)




        self.ax = self.fig.add_subplot(111)
        self.canvas.draw()



        self.plot_layout.addWidget(self.canvas, 1, 1, 20, 20)
        self.toolbar_layout.addWidget(self.toolbar, 21, 1, 1, 20)

    def update_plot(self, effective_subtrace_list, hist_list):
        self.update_plot_signal.emit(effective_subtrace_list, hist_list)

    def update_plot_signal_emitted(self, effective_subtrace_list, hist_list):
        try:
            t0 = time.time()
            for nst, st in enumerate(effective_subtrace_list):
                self.axes[nst, 0].cla()
                self.axes[nst, 0].plot(st)
                self.axes[nst, 1].cla()
                if len(hist_list[nst]) > 0:
                    for h in hist_list[nst]:
                        self.axes[nst, 1].plot(h['bin_edges'], h['hist'])
                else:
                    self.axes[nst, 1].text(
                        0.5,
                        0.5,
                        'HISTOGRAM FAILED',
                        horizontalalignment='center',
                        verticalalignment='center',
                        transform=self.axes[nst, 1].transAxes,
                        color='r',
                        fontsize=30
                    )
                self.fig.tight_layout()
                self.fig.canvas.draw_idle()
        except:
            exc_type, exc_value, exc_tb = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_tb)