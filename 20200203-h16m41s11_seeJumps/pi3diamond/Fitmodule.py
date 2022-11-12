from __future__ import print_function, absolute_import, division

import numpy as np
from lmfit import Model, minimize, Parameters, report_fit

def Fit(x, y, f, d, **kwargs):
    params = Parameters()
    mod = Model(f)
    for key, val in d.items():
        v = True
        v = False if key == 'f0' else v
        params.add(key, value=val, vary=v)

    result = []
    redchi = []
    for method in ['nelder', 'leastsq', 'lbfgsb', 'powell']:
        result.append(mod.fit(y, params, x=x, method=method))
        redchi.append(result[-1].redchi)
    result = result[np.argmin(redchi)]
    return result.best_fit, result.best_values
