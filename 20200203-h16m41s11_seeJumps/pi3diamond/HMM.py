"""
Program to analyze fluorescence data from single shot readout time traces
Using Hidden Markov Model (HMM) "hidden states" can be extracted from a timetrace.
The class HMM initializes and maintaines a HMM for given data (fluorescence trace).

usually the function "hmm = run_prog( data_list, time_increment, **kwargs )" should be used. It returns a fully trained HMM containing all necessary information.

programm code is based on the book "numerical recipes"
"""
from __future__ import print_function, absolute_import, division
__metaclass__ = type

import pickle, Fitting
import numpy as np
import matplotlib.pyplot as plt
import os



# ================================================
class HMM:
    """
    Hidden Markov model
    contains parameters decsribing fluorescence data with a hidden markov model
    These parameters can be refined by using forward-backward algorithm
    the fluorescence data provided as a list of integers should exhibit N more or less destinct Poissonians in its histogram (N fluorescence levels belonging to different states)
    
    data   - list of photon numbers (int)
    N      - number of states (default=2)
    M      - number of symbols/occuring photon numbers (default - calculated by program)
    A      - transition matrix (N x N, list A_ij) (default - calculated by program)
    B      - symbol probability matrix (N x M, list B_ij) (default - calculated by program)
    K      - list of all used symbols (default - calculated by program)
    P      - list of probability of state occupation for each time step (default - calculated by forward-backward algorithm)
    Llhood - log of likelihood for given Probabilities P
    
    dwell_hist_x - data for dwell time histogram
    dwell_hist_y
    trajectory   - state trajectory
    fl_n         - avg fluorescence levels for each state
    viterbi      - avg fluorescence at each time according to trajectory and fl_n
    lifetimes    - avg lifetimes for each state
    
    !!! up to now working for:
    1d data
    N = 2
    """

    def __init__(self, data, A=None, B=None, N=2, dt=0.001, tau=0.05, low=None, high=None):
        """
        initialize all necessary variables
        """

        self.N = N
        self.dt = dt  # time increment in seconds
        self.tau = tau
        self.P = []
        self.Llhood = None
        self.data = data[:]
        self.A = A
        self.B = B
        self.B_init = B

        self.histx = []  # histogram of measurement data (x)
        self.histn = []  # histogram of measurement data (n)
        self.fl_n = []  # fluorescence levels for states
        self.hist_fit = []  # fit to histogram
        self.hist_fit_areas = []  # fit areas for fits to histogram (B_init)
        self.trajectory = []  # most likely state trajectory (state numbers)
        self.viterbi = []  # most likely fluorescence trajectory (avg fluorescence)
        self.dwell_hist_x = []  # dwell time histogram (x)
        self.dwell_hist_y = []  # dwell time histogram (y)
        self.lifetimes = []  # liefetimes for each state

        # create symbol list
        self.make_symbol_list()
        self.M = len(self.K)

        # make symbol dictionary
        self.symbol_map = []
        self.create_symb_dict()

        # make data list according to symbol map
        self.data_ix = []
        self.symbol2index()

        # check validity of A or create A
        self.check_A()

        # check validity of B or create B
        self.check_B(low, high)


    def check_A(self):
        """
        check validity of A
        creates A if it is None
        """
        if self.A == None:
            self.create_A_2()
        elif len(np.shape(self.A)) != 2:
            print('A needs to be a 2 dimensional array')
        elif np.shape(self.A) != (self.N, self.N):
            print('A needs to be a square array N x N')
        else:
            for a in self.A:
                sum = 0.
                for aa in a:
                    sum += aa
                if abs(sum - 1) > 0.001:
                    print('Transition matrix A is not normalized')
        self.A = self.A[:]

    def check_B(self, low, high):
        """
        check validity of B
        creates B if it is None
        """
        if self.B == None:
            if low == None or high == None:
                self.create_B_3()
            else:
                self.create_B_4(low, high)
        elif len(np.shape(self.B)) != 2:
            print('B needs to be a 2-dimensional array')
        elif np.shape(self.B) != (self.N, self.M):
            print('number of rows x columns in B needs to be equal to N x M')
        else:
            for b in self.B:
                sum = 0.
                for bb in b:
                    sum += bb
                if abs(sum - 1) > 0.001:
                    print('symbol probability matrix B is not normalized')
        self.B = self.B[:]

    def make_symbol_list(self):
        """
        create a list of symbols occuring in data list
        """
        self.K = list(set(self.data))
        # K = []
        # data_array = np.array(self.data)
        # K_range = np.arange(data_array.min(), data_array.max() + 1, 1)
        # for k in K_range:
        #     if k in data_array:
        #         K.append(k)
        # self.K = K[:]

    def create_A_2(self):
        """
        create a matrix using timesteps dt and lifetime tau
        """
        a = np.exp(-self.dt / self.tau)
        b = 1 - a
        A = [[a, b], [b, a]]
        self.A = A[:]

    def create_B_3(self):
        """
        create B matrix with two a priori gaussian distributions
        """
        histn, histx = hist(self.data)
        self.histx = histx
        self.histn = histn
        Area = histn.sum()
        Center = (histx.max() + histx.min()) / 2
        center1 = Fitting.Parameter(Center * 0.85)
        center2 = Fitting.Parameter(Center * 1.15)
        area1 = Fitting.Parameter(Area * 0.6)
        area2 = Fitting.Parameter(Area * 0.4)
        f = Fitting.Two_Pois(center1, area1, center2, area2)
        fit_res, dum = Fitting.fit(f, [center1, area1, center2, area2], histn, histx)
        dist1 = Fitting.Poisson(center1, Fitting.Parameter(1))
        dist2 = Fitting.Poisson(center2, Fitting.Parameter(1))
        B = [dist1(np.array(self.K)), dist2(np.array(self.K))]
        self.B = B[:]
        self.B_init = B[:]
        self.hist_fit_areas = [area1(), area2()]
        self.hist_fit = fit_res

    def create_B_4(self, low, high):
        """
        create B matrix with two a priori gaussian distributions
        """
        histn, histx = hist(self.data)
        self.histx = histx
        self.histn = histn
        Area = histn.sum()
        # Center       = (histx.max() + histx.min())/2
        center1 = Fitting.Parameter(low)
        center2 = Fitting.Parameter(high)
        area1 = Fitting.Parameter(Area * 0.9)
        area2 = Fitting.Parameter(Area * 0.1)
        f = Fitting.Two_Pois(center1, area1, center2, area2)

        def f3(area1, area2):
            f1 = Fitting.Poisson(center1, area1)
            f2 = Fitting.Poisson(center2, area2)

            def f(n):
                return ( f1(n) + f2(n) )

            return f

        fit_res, dum = Fitting.fit(f, [center1, area1, center2, area2], histn, histx)
        # fit_res = histn
        dist1 = Fitting.Poisson(center1, Fitting.Parameter(1))
        dist2 = Fitting.Poisson(center2, Fitting.Parameter(1))
        B = [dist1(np.array(self.K)), dist2(np.array(self.K))]
        self.B = B[:]
        self.B_init = B[:]
        self.hist_fit_areas = [area1(), area2()]
        self.hist_fit = fit_res

    def create_symb_dict(self):
        """
        create dictionary for the observable symbols
        it returns the index of each symbol ijn the list of symbols K
        """
        symbol_map = dict(zip(self.K, range(len(self.K))))
        self.symbol_map = symbol_map.copy()

    def symbol2index(self):
        """
        Converts an obeservation symbol sequence into a sequence
        of indices for accessing distribution matrices.
        """
        data_ix = [0. for k in range(len(self.data))]
        for ix, o in enumerate(self.data):
            data_ix[ix] = self.symbol_map[o]
        self.data_ix = data_ix[:]

    def plot_viterbi(self, plot_range=[0, 1000]):
        if len(self.viterbi) == 0:
            pass
        else:
            if len(self.viterbi) < plot_range[1]:
                plot_range = [0, len(self.viterbi)]
            plt.figure()
            plt.plot(list(np.arange(plot_range[0] * self.dt, plot_range[1] * self.dt, self.dt)),
                     self.data[plot_range[0]:plot_range[1]], c=(0.7, 0.7, 0.7), linewidth=5, label='fluor. data')
            plt.plot(list(np.arange(plot_range[0] * self.dt, plot_range[1] * self.dt, self.dt)),
                     self.viterbi[plot_range[0]:plot_range[1]], c=(0, 0, 0.7), linewidth=2, label='viterbi path')
            plt.legend()
            plt.xlabel('time [s]')
            plt.ylabel('photons per bin')

    def plot_hist(self, plot_B=0):
        """
        plot histogram of measurement data plus fit functions B_init
        """
        if self.histx == None:
            pass
        else:
            plt.figure()
            plt.plot(self.histx, self.histn, 'k-', label='histogram of data')
            plt.plot(self.histx, self.hist_fit, c=(0.6, 0.6, 0.6), label='fit of histogram')
            plt.xlabel('# of photons per %.1f ms' % (self.dt * 1000))
            plt.ylabel('# of events')
            for k, (bb, area) in enumerate(zip(self.B_init, self.hist_fit_areas)):
                y = []
                for b in bb:
                    y.append(b * area)
                plt.plot(self.K, y, c=((0.6 + k * 0.3), 0, 0), label='Poissonfit state %i' % k)
            plt.legend()
            if plot_B == True:
                for k, (bb, area) in enumerate(zip(self.B, self.hist_fit_areas)):
                    y = []
                    for b in bb:
                        y.append(b * area)
                    plt.plot(self.K, y, c=(0, 0, (0.6 + k * 0.3)), label='emission fct. state %i' % k)
            plt.legend()


            # def check_A_matrix( self ):
            # """
            # check validity of A
            # creates A if it is None
            # """
            # if self.A == None:
            # self.create_A_1()
            # elif self.A.ndim != 2:
            # print 'A needs to be a 2-dimensional array'
            # elif self.A.shape[0] != self.A.shape[1]:
            # print 'A needs to be a square array'
            # else:
            # for s in self.A.sum(1):
            # if abs(s-1)>0.001:
            # print 'Transition matrix A is not normalized'
            # self.A = self.A.copy()

            # def check_B_matrix( self ):
            # """
            # check validity of B
            # creates B if it is None
            # """
            # if self.B == None:
            # self.create_B_2()
            # elif self.B.ndim != 2:
            # print 'B needs to be a 2-dimensional array'
            # elif self.B.shape != (self.N,self.M):
            # print 'number of rows x columns in B needs to be equal to N x M'
            # else:
            # for s in self.B.sum(1):
            # if abs(s-1)>0.001:
            # print 'symbol probability matrix B is not normalized'
            # self.B = self.B.copy()

            # def create_A_1( self ):
            # """
            # create a matrix using timesteps dt and lifetime tau
            # """
            # a = np.exp(-self.dt/self.tau)
            # b = 1-a
            # A = sc.matrix( [ [ a, b ], [ b, a ] ] )
            # self.A = A.copy()

            # def create_B_1( self ):
            # """
            # create B matrix with two a priori gaussian distributions
            # """
            # K = np.array(self.K)
            # width = K.max() - K.min()
            # N1 = (K.max()+K.min())/2.*1.1
            # N2 = (K.max()+K.min())/2.*0.9
            # dist_y1 = 1/(np.sqrt(N1*2.*np.pi))*np.exp(-1./2*(K-N1)**2/N1)
            # dist_y2 = 1/(np.sqrt(N2*2.*np.pi))*np.exp(-1./2*(K-N2)**2/N2)
            # dist_y1 = dist_y1/dist_y1.sum()
            # dist_y2 = dist_y2/dist_y2.sum()
            # B = sc.matrix( [ dist_y2, dist_y1 ] )
            # self.B = B.copy()

            # def create_B_2( self ):
            # """
            # create B matrix with two a priori gaussian distributions
            # """
            # histn, histx = hist( self.data )
            # Area    = histn.sum()
            # Center  = histx.max() - histx.min()
            # center1 = Fitting.Parameter( Center*0.8 )
            # center2 = Fitting.Parameter( Center*1.2 )
            # area1   = Fitting.Parameter( Area*0.6 )
            # area2   = Fitting.Parameter( Area*0.4 )
            # f       = Fitting.Two_Pois( center1, center2, area1, area2 )
            # Fitting.fit(f, [ center1, center2, area1, area2 ], histn, histx )
            # dist1   = Fitting.Normal( center1, Fitting.Parameter(1), Fitting.Parameter(np.sqrt(center1())), Fitting.Parameter(0) )
            # dist2   = Fitting.Normal( center2, Fitting.Parameter(1), Fitting.Parameter(np.sqrt(center2())), Fitting.Parameter(0) )
            # B       = sc.matrix( [ dist1(self.K), dist2(self.K) ] )
            # self.B = B.copy()


# ================================================
# ================================================
# ================================================


def run_prog(data, dt, **kwargs):
    """
    supply data, time increment and optional parameters and get everything done
    return is a fully trained hmm including all necessary data
    """
    if 'A' in kwargs:
        A = kwargs['A']
    else:
        A = None
    if 'B' in kwargs:
        B = kwargs['B']
    else:
        B = None
    if 'N' in kwargs:
        N = kwargs['N']
    else:
        N = 2
    # if 'dt'      in kwargs: dt      = kwargs['dt']
    # else: dt      = 0.001
    if 'tau' in kwargs:
        tau = kwargs['tau']
    else:
        tau = 0.05
    if 'f_plot' in kwargs:
        f_plot = kwargs['f_plot']
    else:
        f_plot = 0
    if 'f_print' in kwargs:
        f_print = kwargs['f_print']
    else:
        f_print = 0
    if 'dL' in kwargs:
        dL = kwargs['dL']
    else:
        dL = 0.001
    if 'steps' in kwargs:
        steps = kwargs['steps']
        dL = None
    else:
        steps = 10
    if 'low' in kwargs:
        low = kwargs['low']
    else:
        low = None
    if 'high' in kwargs:
        high = kwargs['high']
    else:
        high = None

    hmm = HMM(data, A, B, N, dt, tau, low, high)
    if f_plot == True: hmm.plot_hist(0)
    iterations(hmm, f_plot, f_print, dL=dL, steps=steps)
    if f_plot == True: hmm.plot_hist(1)
    most_likely_state(hmm, 0)
    dwell_times(hmm, f_plot)
    results(hmm, f_plot, f_print)

    return hmm


# ================================================
# ================================================
# ================================================


def fwd_bwd_3(hmm):
    """
    forward backward algorithm to calculate the likelyhood of each step
    """
    data_ix = hmm.data_ix[:]
    N = hmm.N
    M = hmm.M
    L = len(data_ix)
    alpha = [[0. for l in range(N)] for k in range(L)]  # [ [0.]*N ]*L
    beta = [[0. for l in range(N)] for k in range(L)]  # [ [0.]*N ]*L
    arnrm = [0 for k in range(L)]  # [ 0 ]*L
    brnrm = [0 for k in range(L)]  # [ 0 ]*L
    P = [[0. for l in range(N)] for k in range(L)]  # [ [0.]*N ]*L
    A = hmm.A[:]
    B = hmm.B[:]
    A_new = A[:]
    B_new = [[0. for l in range(M)] for k in range(N)]  # [ [0.]*M ]*N
    BIG = 1e10
    BIGI = 1 / BIG

    # forward
    for i in xrange(N):
        alpha[0][i] = B[i][data_ix[0]]
    for t in xrange(1, L, 1):
        asum = 0.
        for j in xrange(N):
            sum = 0.
            for i in xrange(N):
                sum += ( alpha[t - 1][i] * A[i][j] * B[j][data_ix[t]] )
            alpha[t][j] = sum
            asum += sum
        arnrm[t] = arnrm[t - 1]
        if asum < BIGI:
            arnrm[t] += 1
            for j in xrange(N):
                alpha[t][j] *= BIG

    # backward
    for i in xrange(N):
        beta[L - 1][i] = 1.
    for t in xrange(L - 2, -1, -1):
        bsum = 0.
        for i in xrange(N):
            sum = 0.
            for j in xrange(N):
                sum += ( A[i][j] * B[j][data_ix[t + 1]] * beta[t + 1][j] )
            beta[t][i] = sum
            bsum += sum
        brnrm[t] = brnrm[t + 1]
        if bsum < BIGI:
            brnrm[t] += 1
            for j in xrange(N):
                beta[t][j] *= BIG

    # likelihood
    lhood = 0.
    for i in xrange(N):
        lhood += ( alpha[0][i] * beta[0][i] )
    lrnrm = arnrm[0] + brnrm[0]
    while lhood < BIGI and lhood != 0:
        lhood *= BIG
        lrnrm += 1
    Llhood = np.log10(lhood) - lrnrm * np.log10(BIG)
    # print 'log10(likelihood) =', Llhood

    # probability of states
    for t in xrange(L):
        sum = 0.
        for i in xrange(N):
            P[t][i] = p = ( alpha[t][i] * beta[t][i] )
            sum += p
        for i in xrange(N):
            P[t][i] /= sum

    # baumwelch
    powtab = [0. for k in range(30)]  # [ 0. ]*30
    pt_offset = 9
    for i in xrange(len(powtab)):
        powtab[i] = pow(BIGI, i - pt_offset)
    for i in xrange(N):
        denom = 0.
        for k in xrange(M):
            B_new[i][k] = 0.
        for t in xrange(L - 1):
            term = ( alpha[t][i] * beta[t][i] / lhood ) * powtab[arnrm[t] + brnrm[t] - lrnrm + pt_offset]
            denom += term
            B_new[i][data_ix[t]] += term
        for j in xrange(N):
            num = 0.
            for t in xrange(L - 1):
                num += (alpha[t][i] * B[j][data_ix[t + 1]] * beta[t + 1][j] * powtab[
                    arnrm[t] + brnrm[t + 1] - lrnrm + pt_offset] / lhood )
            A_new[i][j] *= ( num / denom )
        for k in xrange(M):
            B_new[i][k] /= denom

    hmm.A = A_new[:]
    hmm.B = B_new[:]
    hmm.P = P[:]
    hmm.Llhood = Llhood


def iterations(hmm, f_plot=0, f_print=0, **kwargs):
    """
    do forward-backward steps including baumwelch
    number of steps is "steps" or, if "dL" is supplied, runs until improvement of Llhood, is less than dL
    """
    if 'dL' in kwargs:
        dL = kwargs['dL']
    else:
        dL = None
    if 'steps' in kwargs:
        steps = kwargs['steps']
    else:
        steps = 10

    Llhood_list = []
    x = []
    if f_plot == True:
        plt.figure()
        plt.xlabel('# of iterations')
        plt.ylabel('log( likelihood )')
    if dL == None:
        for k in xrange(steps):
            fwd_bwd_3(hmm)
            Llhood_list.append(hmm.Llhood)
            if f_print == True:
                print('LogLikelihood:')
                print('\t%f' % hmm.Llhood)
            if f_plot == True:
                x.append(k + 1)
                plt.plot(x[-1], Llhood_list[-1], 'bo')
    else:
        fwd_bwd_3(hmm)
        if f_print == True:
            print('LogLikelihood:')
            print('\t%f' % hmm.Llhood)
        Llhood_list.append(hmm.Llhood)
        x.append(1)
        if f_plot == True: plt.plot(x[-1], Llhood_list[-1], 'bo')

        fwd_bwd_3(hmm)
        if f_print == True: print('\t%f' % hmm.Llhood)
        Llhood_list.append(hmm.Llhood)
        x.append(2)
        if f_plot == True: plt.plot(x[-1], Llhood_list[-1], 'bo')

        k = 1
        DL = 1.
        D0 = Llhood_list[1] - Llhood_list[0]
        while DL > dL and k < (steps * 5):
            k += 1
            fwd_bwd_3(hmm)
            if f_print == True: print('\t%f' % hmm.Llhood)
            Llhood_list.append(hmm.Llhood)
            DL = abs(( Llhood_list[k] - Llhood_list[k - 1] ) / D0)
            x.append(k + 1)
            if f_plot == True: plt.plot(x[-1], Llhood_list[-1], 'bo')


def most_likely_state(hmm, f_plot=0):
    """
    calculates the most likely state trajectory
    """
    P = hmm.P[:]
    trajectory = []
    for p in P:
        state = 0
        value = p[0]
        for v in p[1:]:
            if v > value:
                value = v
                state += 1
        trajectory.append(state)
    if f_plot == True:
        plt.figure()
        plt.xlabel('time [bin]')
        plt.ylabel('state')
        plt.plot(trajectory[:1000])
    hmm.trajectory = trajectory[:]


def dwell_times(hmm, f_plot=0):
    """
    extract the occurences of dwell times for each state from the state trajectory
    """
    dwell_hist_y = []
    dwell_hist_x = []
    # prepare dwell_hist_y
    for n in range(hmm.N):
        dwell_hist_y.append([])
    # extract dwell times from trajectory
    state0 = hmm.trajectory[0]
    n = 1
    for state in hmm.trajectory[1:]:
        if state == state0:
            n += 1
        else:
            length = len(dwell_hist_y[state0])
            while length < n:
                dwell_hist_y[state0].append(0)
                length = len(dwell_hist_y[state0])
            dwell_hist_y[state0][n - 1] += 1
            state0 = state
            n = 1
    # generate dwell_hist_x values
    for n, k in enumerate(dwell_hist_y):
        x = []
        for i, l in enumerate(k):
            x.append((i + 1) * hmm.dt * 1000)
        dwell_hist_x.append(x)
    # plotting
    if f_plot == True:
        plt.figure()
        plt.xlabel('lifetime [ms]')
        plt.ylabel('# of events')
        for k, (x, y) in enumerate(zip(dwell_hist_x, dwell_hist_y)):
            plt.semilogy(x, y, label='state %i' % k)
        plt.legend()
    hmm.dwell_hist_y = dwell_hist_y[:]
    hmm.dwell_hist_x = dwell_hist_x[:]


def results(hmm, f_plot=0, f_print=0):
    """
    calculate results from fully trained hmm
    """
    # fluorescence levels for each state
    fl_n = []
    for k, b in enumerate(hmm.B):
        maxi = 0
        max = 0.
        for i, bb in enumerate(b):
            if bb > max:
                max = bb
                maxi = i
        maxx = hmm.K[maxi]
        if maxx == 0:
            maxx = 0.1
        center = Fitting.Parameter(maxx)
        area = Fitting.Parameter(1.)
        # sigma  = Fitting.Parameter( np.sqrt( maxx ) )
        # offset = Fitting.Parameter( 0. )
        #f      = Fitting.Normal( center, area, sigma, offset)
        f = Fitting.Poisson(center, area)
        Fitting.fit(f, [center, area], b, hmm.K)
        fl_n.append(center())
    hmm.fl_n = fl_n[:]
    # viterbi
    viterbi = hmm.trajectory[:]
    for k, state in enumerate(viterbi):
        viterbi[k] = fl_n[state]
    hmm.viterbi = viterbi[:]
    # lifetimes
    lifetimes = []
    for k, a in enumerate(hmm.A):
        lifetimes.append(-hmm.dt / np.log(a[k]))
    hmm.lifetimes = lifetimes[:]
    # print results
    if f_print == True:
        print('results:')
        for k, (fl, tau) in enumerate(zip(fl_n, lifetimes)):
            print('\tstate %i:' % k)
            print('\t\t fluor. level %.1f [photons/dt], lifetime %.2f [dt]' % (fl / hmm.dt, tau))
    # plot results
    if f_plot == True:
        hmm.plot_viterbi([0, 10000])


def load_pyd(fname):
    """
    load data from *.pyd files
    """
    tracefile = open(fname)
    timetrace = pickle.load(tracefile)
    tracefile.close()
    return timetrace


def load_asc(fname):
    """
    load data from ascii files
    """
    # read data from file
    file_in = open(fname, 'r')
    lines = file_in.readlines()
    file_in.close()
    # look for header "data"
    start = 0
    for i, line in enumerate(lines[0:10]):
        if 'Data' in line:
            start = i + 1
        i += 1
    # print i, start, line
    # create data list
    data = []
    for line in lines[start:]:
        data.append(int(line))
    return data


def load_asc_2(fname):
    """
    load data from ascii files out of T1 matrix measurement
    5 header lines
    """
    file_in = open(fname, 'r')
    lines = file_in.readlines()
    file_in.close()
    x = float(lines.pop(0))
    y = float(lines.pop(0))
    z = float(lines.pop(0))
    phi = float(lines.pop(0))
    dummy = lines.pop(0)
    data = []
    for line in lines:
        data.append(int(line))
    return data, x, y, z, phi


def hist(y):
    """
    compute histogramm of data with an entry for every integer within the valid symbol range
    """
    y_array = np.array(y)
    x_min = y_array.min()
    x_max = y_array.max()
    bins = range(x_min, x_max + 2, 1)
    ( histn, histx ) = np.histogram(y_array, bins)
    return histn, histx[:-1]


def binning(data_in, N):
    """
    binning of data
    """
    data_out = []
    a = 0
    for i, data in enumerate(data_in):
        a += data
        if i % N == 0:
            data_out.append(a)
            a = 0
    return data_out


def analyze_trace(trace, binlength, tau_d, tau_b, f_print=True, f_plot=False, binning=1):
    if binning != 1:
        data_new = []
        i = 0
        while i < len(trace) - binning:
            data_new.append(sum(trace[i:i + binning]))
            i += binning
        trace = data_new
        binlength = binlength * binning
    A = [[np.exp(-binlength / tau_d), 1 - np.exp(-binlength / tau_d)],
         [1 - np.exp(-binlength / tau_b), np.exp(-binlength / tau_b)]]
    result = run_prog(trace, binlength, **{'A': A, 'f_print': f_print, 'f_plot': f_plot})
    return result


def analyze(filename, dt, tau_d=0.3, tau_b=0.003, low=None, high=None, f_print=True, f_plot=False, binning=1):
    data = load_asc('D:\\Data\\Timetrace\\' + filename + '_Trace.dat')
    if binning != 1:
        data_new = []
        i = 0
        while i < len(data) - binning:
            data_new.append(sum(data[i:i + binning]))
            i += binning
        data = data_new
        dt = dt * binning
    A = [[np.exp(-dt / tau_d), 1 - np.exp(-dt / tau_d)], [1 - np.exp(-dt / tau_b), np.exp(-dt / tau_b)]]
    result = run_prog(data, dt, **{'A': A, 'low': low, 'high': high, 'f_print': f_print, 'f_plot': f_plot})
    return filename, result


def analyze2(filename, dt, tau_d=0.3, tau_b=0.003, low=None, high=None, f_print=True, f_plot=False, binning=1):
    data = load_asc('D:\\Data\\Timetrace\\' + filename + '_Trace.dat')
    if binning != 1:
        data_new = []
        i = 0
        while i < len(data) - binning:
            data_new.append(sum(data[i:i + binning]))
            i += binning
        data = data_new
        dt = dt * binning
    A = [[np.exp(-dt / tau_d), 1 - np.exp(-dt / tau_d)], [1 - np.exp(-dt / tau_b), np.exp(-dt / tau_b)]]

    length = int(len(data) / 10.)
    traces = []
    results = []
    t1_bright = []
    t1_dark = []
    print('now')
    for i in range(10):
        traces.append(data[i * length:(i + 1) * length])
    for trace in traces:
        results.append(run_prog(trace, dt, **{'A': A, 'low': low, 'high': high, 'f_print': f_print, 'f_plot': f_plot}))
        dark = 0
        bright = 1
        if results[-1].fl_n[0] > results[-1].fl_n[1]:
            dark = 1
            bright = 0
        t1_bright.append(results[-1].lifetimes[bright])
        t1_dark.append(results[-1].lifetimes[dark])
    result = run_prog(data, dt, **{'A': A, 'low': low, 'high': high, 'f_print': f_print, 'f_plot': True})
    error_dark = (1. / (len(t1_dark) * (len(t1_dark) - 1)) * sum((np.array(t1_dark) - np.mean(t1_dark)) ** 2)) ** 0.5
    print('error_dark: ' + str(error_dark * 1000))
    error_bright = (1. / (len(t1_bright) * (len(t1_bright) - 1)) * sum(
        (np.array(t1_bright) - np.mean(t1_bright)) ** 2)) ** 0.5
    print('error_bright: ' + str(error_bright * 1000))
    print('t1_bright:')
    for t in t1_bright:
        print(t * 1000)
    print('t1_dark:')
    for t in t1_dark:
        print(t * 1000)
    filename = 'D:\\Data\\NVswitch\\' + filename
    export_viterbi(result, filename)
    export_dwelltimes(result, filename)
    export_distributions(result, filename)
    file = open(filename + '_results.txt', 'w')
    file.write('results:\n')
    error = [error_dark, error_bright]
    for k, (fl, tau) in enumerate(zip(result.fl_n, result.lifetimes)):
        file.write('\tstate %i:\n' % k)
        file.write('\t\t fluor. level %.1f kcounts, lifetime %.2f ms +- %.3f ms\n' % (
            fl / result.dt / 1000, tau * 1000, error[k] * 1000))
    file.close()
    return filename, result


def export_viterbi(hmm, filename_base):
    """
    export time, data and viterbi vector as three columns into tab separated file
    """
    file = open(filename_base + '_viterbi.dat', 'w')
    i = 0
    for x1, x2 in zip(hmm.data, hmm.viterbi):
        file.write('%g\t%i\t%g\n' % (i * hmm.dt, x1, x2))
        i += 1
    file.close()


def export_dwelltimes(hmm, filename_base):
    """
    export dwelltime histograms
    """
    file = open(filename_base + '_dwelltimes.dat', 'w')
    for xx, yy in zip(hmm.dwell_hist_x, hmm.dwell_hist_y):
        for x, y in zip(xx, yy):
            file.write('%g\t%i\n' % (x, y))
        file.write('\n')
    file.close()


def export_results(hmm, filename_base):
    """
    export results
    """
    file = open(filename_base + '_results.txt', 'w')
    file.write('results:\n')
    for k, (fl, tau) in enumerate(zip(hmm.fl_n, hmm.lifetimes)):
        file.write('\tstate %i:\n' % k)
        file.write('\t\t fluor. level %.1f kcounts, lifetime %.2f ms\n' % (fl / hmm.dt / 1000, tau * 1000))
    file.close()


def export_distributions(hmm, filename_base):
    """
    export distributions of states (emission functions)
    """
    file = open(filename_base + '_distributions.dat', 'w')
    B = np.array(hmm.B)
    C = B.transpose()
    for k, b in zip(hmm.K, C):
        file.write('%g' % k)
        for bb in b:
            file.write('\t%g' % bb)
        file.write('\n')
    file.close()


def save_all(result, folder=None, filename='result'):
    """filenames extensions are dat and generated automatically"""
    folder = os.getcwd() if folder is None else folder
    filename = os.path.join(folder, filename)
    export_viterbi(result, filename)
    export_dwelltimes(result, filename)
    export_results(result, filename)
    export_distributions(result, filename)