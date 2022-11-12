from lmfit import Model
import numpy as np
import qutip_enhanced
reload(qutip_enhanced)
def twocosinefunc(x,a1,a2,f1,f2,T2,phi1,phi2):
    return np.exp(-(x/T2)**2)*(a1*np.cos(x*f1+phi1) + a2*np.cos(x*f2+phi2))

def fourcosine(x,a1,a2,a3,a4,f1,f2,f3,f4,T2,phi1,phi2,phi3,phi4):
    return np.exp(-(x / T2) ** 2) * (a1 * np.cos(x * f1 + phi1) +
                                     a2 * np.cos(x * f2 + phi2)+
                                     a3 * np.cos(x * f3 + phi3) +
                                     a4 * np.cos(x * f4 + phi4))
def run(pi3d):

    gmodel = Model(twocosinefunc)
    #gmodel = Model(fourcosine)
    def guesfun(*args, **kwargs):
        return gmodel.make_params(a1=0.3, a2=0.3, f1=1, f2=1.5, phi1=0, phi2=0, T2=15.0)
        #return gmodel.make_params(a1=0.1, a2=0.1, f1=1, f2=1.4, phi1=0, phi2=0, T2=12.0,
         #                         a3=0.1, a4=0.1, f3=2, f4=2.1, phi3=np.pi, phi4=np.pi)
    gmodel.guess = guesfun

    pi3d.NuclearOPselectron_t2_20181212h01m20s06.pld.custom_model = gmodel


def run_echo_fit(pi3d):

    gmodel = qutip_enhanced.lmfit_custom_models.ExpPowerDecayModel()

    pi3d.NuclearOPselectron_t2_20181212h16m57s52.pld.custom_model = gmodel