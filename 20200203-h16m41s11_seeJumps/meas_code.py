# coding=utf-8
from pi3diamond import pi3d
import numpy as np
import os
import UserScripts.helpers.sequence_creation_helpers as sch; reload(sch)
import UserScripts.helpers.shared as shared; reload(shared)
import multi_channel_awg_seq as MCAS; reload(MCAS)
import UserScripts.helpers.snippets_awg as sna; reload(sna)
import UserScripts.helpers.shared as ush;reload(ush)
from qutip_enhanced import *

from collections import OrderedDict

seq_name = os.path.basename(__file__).split('.')[0]
nuclear = sch.create_nuclear(__file__)
with open(os.path.abspath(__file__).split('.')[0] + ".py", 'r') as f:
    meas_code = f.read()

__TAU_HALF__ = 2*192/12e3

ael = 1.0

def ret_ret_mcas(pdc):
    def ret_mcas(current_iterator_df):
        mcas = MCAS.MultiChSeq(seq_name=seq_name, ch_dict={'2g': [1, 2], '128m': [1]})
        for idx, _I_ in current_iterator_df.iterrows():
            n = _I_['ssr_reps']
            if _I_['state_result'] == 'n+':
                sna.ssr(mcas, frequencies=[pi3d.tt.mfl({'14n': [+1, 0, -1], '13c414': [+.5]}),
                                           pi3d.tt.mfl({'14n': [+1, 0, -1], '13c414': [-.5]})], nuc='13c414', robust=True, repetitions=int(n), mixer_deg=-90, step_idx=0)

            elif _I_['state_result'] == 'nn+':
                sna.ssr(mcas, frequencies=[pi3d.tt.mfl({'14n': [+1, 0, -1], '13c414': [+.5, -.5], '13c90': [+.5]}),
                                           pi3d.tt.mfl({'14n': [+1, 0, -1], '13c414': [+.5, -.5], '13c90': [-.5]})],
                                            nuc='13c90', robust=True, repetitions=int(n), mixer_deg=-90, step_idx=0)

            elif _I_['state_result'] == '+':
                freq1 = pi3d.tt.mfl({'14N': [+1]}, ms_trans=_I_['ms'])
                freq2 = pi3d.tt.mfl({'14N': [0]}, ms_trans=_I_['ms'])
                #freq3 = pi3d.tt.mfl({'14N': [-1]}, ms_trans=_I_['ms'])
                sna.ssr(mcas, frequencies=[freq1, freq2], nuc='14N+1', robust=True, repetitions=int(n), mixer_deg=-90, step_idx=0)

            elif _I_['state_result'] in ["".join(i) for i in itertools.product(['+', '0', '-'], ['+', '-'], ['+', '-'])]:
                sna.ssr_single_state(mcas, state=_I_['state_result'], step_idx=0)

            pi3d.gated_counter.set_n_values(mcas)
        return mcas
    return ret_mcas
# __LASER_DUR_DICT__ = {'14n+1': .175, '14n-1': .175, '14n': .175, '14n0': .2, '13c414': .2, '13c90': .21, 'single_state': .2, 'charge_state': 2000.0}
def settings(pdc={}):
    ana_seq=[
        ['result', '<', 0, 0, 10, 2],
        #['result', '<', 0, 0, 10, 2]
        #['result', '<', 'auto', 123123, 1, 1],
    ]
    sch.settings(
        nuclear=nuclear,
        ret_mcas=ret_ret_mcas(pdc),
        analyze_sequence=ana_seq,
        pdc=pdc,
        meas_code=meas_code
    )
    nuclear.x_axis_title = 'tau_half [mus]'
    nuclear.analyze_type = 'consecutive'


    pi3d.gated_counter.trace.analyze_type = 'consecutive'
    pi3d.gated_counter.trace.consecutive_valid_result_numbers = [0]
    pi3d.gated_counter.trace.average_results = True

    nuclear.parameters = OrderedDict(
        (
            ('sweeps', range(100)),
            ('rabi_period', [0.1]),
            ('state_result', ['n+','+','nn+']),
            #('state_result', ["".join(i) for i in itertools.product(['+','0','-'], ['+', '-'], ['+', '-'])]+['+','n+','nn+']),
            # ('state_result', ["".join(i) for i in itertools.product(['+', '0'], ['+', '-'], ['+','-'])]),
            # ('state_init', ["".join(i) for i in itertools.product(['+'], ['+','-'], ['+','-'])]),
            #('state_init', ['+++']),
            #('state_result', ['+++', '0++']),
            ('ms', [-1]),
            ('ddt', ['hahn']),# 'fid','hahn', 'xy4', 'xy16', 'kdd4', 'kdd16']),
            ('n_rep_dd', [1]),
            ('ssr_reps',range(400,2000,200))
        )
    )
    nuclear.number_of_simultaneous_measurements = 1#len(nuclear.parameters['phase_pi2_2'])

def run_fun(abort, **kwargs):
    pi3d.readout_duration = 150e6
    nuclear.debug_mode = False
    settings()
    nuclear.run(abort)
