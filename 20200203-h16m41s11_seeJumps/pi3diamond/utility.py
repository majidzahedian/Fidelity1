from __future__ import print_function, absolute_import, division

import Queue
import cPickle
import threading
import time

from chaco.api import PlotGraphicsContext
from chaco.tools.simple_zoom import SimpleZoom

try:
    from pi3diamond import pi3d
except:
    pi3d = None
    CloseHandler = None
    print('pi3diamond could not be imported')

class Singleton:
    def __new__(cls, *a, **k):
        if not hasattr(cls, '_inst'):
            cls._inst = super(Singleton, cls).__new__(cls)
        return cls._inst

class GetSetItemsMixin:
    """Provides save, load, save figure methods. Useful with HasTraits models.
    Data is stored in a dictionary with keys that are strings and identical to
    class attribute names. To save, pass a list of strings that denote attribute names.
    Load methods accept a filename. The dictionary is read from file and attributes
    on the class are set (if necessary created) according to the dictionary content. 
    """

    _file_mode_map = {'asc':'', 'bin':'b'}
    _pickle_mode_map = {'asc':0, 'bin':1}
    
    def set_items(self, d):
        for key, val in d.items():
            try:
                setattr(self, key, val)
            except:
                pi3d.logger.warning("failed to set item '"+key+"'")

    def get_items(self, keys):
        d = {}
        for key in keys:
            d[key] = getattr(self, key)
        return d
    
    def copy_items(self, keys):
        d = {}
        for key in keys:
            item = getattr(self, key)
            if hasattr(item,'copy'):
                d[key] = item.copy()
            else:
                d[key] = item
        return d
    
    def dump_items(self, keys, filename, mode='asc'):
        """Copy class attributes into a dictionary and save it to disc using cPickle.
        'keys' is a list of strings denoting the attributes to be saved.
        """
        f = open(filename, 'w'+self._file_mode_map[mode])
        cPickle.dump(self.get_items(keys), f, self._pickle_mode_map[mode])
        f.close()

    def load_items(self, filename, mode='asc'):
        f = open(filename, 'r'+self._file_mode_map[mode])
        self.set_items(cPickle.load(f))
        f.close()

    def save_figure(self, figure, filename):
        gc = PlotGraphicsContext(figure.outer_bounds, dpi=72)
        gc.render_component(figure)
        gc.save(filename)


class AspectZoomTool(SimpleZoom):

    def _do_zoom(self):
        """ Does the zoom operation.
        """
        # Sets the bounds on the component using _cur_stack_index
        low, high = self._current_state()
        orig_low, orig_high = self._history[0]
    
        if self._history_index == 0:
            if self.tool_mode == "range":
                mapper = self._get_mapper()
                mapper.range.low_setting = self._orig_low_setting
                mapper.range.high_setting = self._orig_high_setting
            else:
                x_range = self.component.x_mapper.range
                y_range = self.component.y_mapper.range
                x_range.low_setting, y_range.low_setting = \
                    self._orig_low_setting
                x_range.high_setting, y_range.high_setting = \
                    self._orig_high_setting

                # resetting the ranges will allow 'auto' to pick the values
                x_range.reset()
                y_range.reset()
               
        else:   
            if self.tool_mode == "range":
                mapper = self._get_mapper()
                if self._zoom_limit_reached(orig_low, orig_high, low, high, mapper):
                    self._pop_state()
                    return
                mapper.range.low = low
                mapper.range.high = high
            else:
                for ndx in (0, 1):
                    mapper = (self.component.x_mapper, self.component.y_mapper)[ndx]
                    if self._zoom_limit_reached(orig_low[ndx], orig_high[ndx],
                                                low[ndx], high[ndx], mapper):
                        # pop _current_state off the stack and leave the actual
                        # bounds unmodified.
                        self._pop_state()
                        return
                x_range = self.component.x_mapper.range
                y_range = self.component.y_mapper.range
                x_range.low, y_range.low = low
                x_range.high, y_range.high = high

        plot = self.component.container
        plot.aspect_ratio = (x_range.high - x_range.low) / (y_range.high - y_range.low)
                
        self.component.request_redraw()
        return

def get_timestamp_string():
    return time.strftime('%y-%m-%d_%Hh%Mm%S', time.localtime())


class TrackerManagerSingleton( Singleton ):
    def compare_order(self, a, b):
        return cmp(a.pause_order, b.pause_order)

    def get_clients_sorted(self):
        clients = []
        for thread in threading.enumerate():
            if hasattr(thread, 'is_tracker_client') and thread.is_tracker_client and thread.isAlive():
                clients.append(thread)
        return sorted(clients, self.compare_order)

    #def request_pause_old(self):
    #    for thread in self.get_clients_sorted():
    #        thread.toggle_pause.set()
    #        thread.pausing.wait()

    def request_pause(self):
        for thread in self.get_clients_sorted():
            thread.command_queue.put('pause_request')
            thread.command_queue.join()

    def release_pause(self):
        clients = self.get_clients_sorted()
        clients.reverse() #reverse() returns None. from python import bad surprise at runtimes!
        for thread in clients:
            thread.command_queue.put('resume')
            thread.command_queue.join()

class ThreadedMeasurementMixin():
    """Provides thread management, interaction with the global tracker_manager and logging.
    * Set ignore_tracker_manager to True to disable offering measurement pauses to the tracker.
    * Override run() to implement your measurement.
    ** Don't forget to offer_pause() once in a while.
    ** Also include the line "if self.stop_request.is_set(): break" whenever you can offer aborting your thread worker loop.
    * Use pi3d.logger.info(), pi3d.logger.error() etc. if you have to say something.
    """
    default_join_timeout = 20
    default_track_timeout = 600
    tracker_manager = TrackerManagerSingleton()

    def __init__(self, is_tracker_client = True, pause_order = 5):
        self.is_tracker_client = True
        self.pause_order = 5
        self.stop_request = threading.Event()
        self.stop_request.clear()
        self.thread = None

    def start(self, stop_request=None):
        """Start the worker thread with target=run()"""
        ## try to clean up
        self.stop()
        if stop_request is None:
            self.stop_request.clear()
        else:
            self.stop_request = stop_request
        pi3d.logger.debug('Active threads: %i' % threading.activeCount())

        ## allocate new thread
        for i in range(1,10):
            try:
                self.thread = threading.Thread(target=self.run,
                                        name=self.__class__.__name__ + get_timestamp_string())
                break
            except Exception as e:
                ## we've had problems with a "can't start new thread" exception. therefore we retry here
                pi3d.logger.exception(str(e))
                pi3d.logger.debug('Active threads: %i' % threading.activeCount())
                time.sleep(5)
        ## set flags for tracker
        # TODO: maybe subclass thread and add these parameters, here we're storing data twice. but not now...
        self.thread.is_tracker_client = self.is_tracker_client
        self.thread.pause_order = self.pause_order
        self.thread.command_queue = Queue.Queue()

        ## offer a pause, and then start
        self.offer_pause()
        self.thread.start()

    def stop(self):
        """Stop the worker thread"""
        pi3d.logger.debug('setting stop request...')
        self.stop_request.set() ##if run() is not executed in extra thread, stop() still should abort run() by setting this flag
        if self.thread is None:
            pi3d.logger.debug('no thread to stop, returning, ...')
            return
        elif not self.thread.isAlive():
            pi3d.logger.debug('thread no longer alive, returning, ...')
            return
        elif self.thread is threading.current_thread():
            pi3d.logger.debug('stop request from current thread, returning...')
            return
        pi3d.logger.debug('waiting for thread to finish ...')
        self.thread.join(self.default_join_timeout)
        pi3d.logger.debug('Active threads: %i' % threading.activeCount())

    def offer_pause(self):
        """If the tracker_manager wants to do some tracking, pause until it is done"""
        if self.thread is None or not self.thread.is_tracker_client: return
        q = self.thread.command_queue
        if not q.empty():
            cmd = q.get_nowait()
            if cmd == 'pause_request':
                pi3d.logger.debug('preparing pause...')
                self.prepare_pause()
                q.task_done()
                cmd = q.get(timeout=self.default_track_timeout)
                if cmd == 'resume':
                    pi3d.logger.debug('releasing pause...')
                    self.prepare_resume()
                    q.task_done()
                else:
                    raise RuntimeError("Expected 'resume' but got:" + cmd)
            else:
                raise RuntimeError('Expected pause_request but got:' + cmd)

    def prepare_pause(self):
        """Prepare optimal conditions for the tracker (i.e. switch off MW pulses etc.)
        """
        pass

    def prepare_resume(self):
        """Resume normal conditions after tracking pause
        """
        pass

    def run(self):
        """NEVER CALL THIS METHOD DIRECTLY. Use start() and stop() instead. Override it to implement your measurement.
        Don't forget to offer_pause() once in a while.
        Also, include the line "if self.stop_request.is_set(): break" in your main loop if you can accept thread stop.
        """
        while(True):
            pi3d.logger.debug("Yeah, still taking data like the LHC!")
            if self.stop_request.is_set():
                pi3d.logger.debug("Stopping")
                break
            pi3d.logger.debug("Offering pause")
            time.sleep(1)
            self.offer_pause()
            time.sleep(1)
            pi3d.logger.debug("Back from pause")
            time.sleep(2)

class PickleMixin:
    class Dump:
        pass

    def __init__(self):
        self.filename_base = self.__class__.__name__ + '_' + get_timestamp_string()
        self.attributes_to_pickle = []

    def pickle(self):
        d = PickleMixin.Dump()
        for a in self.attributes_to_pickle:
            setattr(d, a, getattr(self, a))
        file = open(self.filename_base + '.pkl', 'wb')
        cPickle.dump(d, file, -1)

    def unpickle(self, filename_base=None):
        if filename_base is None: filename_base = self.filename
        file = open(filename_base + '.pkl', 'rb')
        d = cPickle.load(file)
        for a in self.attributes_to_pickle:
            setattr(self, a, getattr(d, a))


class History(object):
    """History of length 'length'."""
    def __init__(self, length):
        self.length = length
        self.items = [ ]
        self.i = 0

    def get(self):
        return self.items[self.i]

    def back(self):
        if self.i != 0:
            self.i = self.i - 1
        return self.items[self.i]

    def forward(self):
        if self.i != len(self.items) - 1:
            self.i = self.i + 1
        return self.items[self.i]

    def put(self, item):
        while self.i < len(self.items) - 1:
            self.items.pop()
        if self.i == self.length - 1:
            self.items.pop(0)
        self.items.append(item)
        self.i = len(self.items) - 1

    def setlength(self, length):
        while len(self.items) > length:
            self.items.pop(0)
            self.i = self.i - 1
        self.length = length