# fitting stuff
# from www.scipy.org
# ==============================
#
# I.   define parameters
#		height = Parameter(10)
#		sigma  = Parameter(2)
#		mu     = Parameter(0)
# II.  define your function:
#		e.g.
# 		def f(x): return height() * exp(-((x-mu())/sigma())**2)
#		or
#		f = Gauss(center,height,width,offset)
# III. do the fit
#		data and x have to be an numpy.array
#		fitfunction(x), x = fit(f, [mu, sigma, height], data, x )
#
# ...	use adjusted parameters, fitfunction(x) and/or x
# 
# ==============================
from __future__ import print_function, absolute_import, division
from scipy import optimize
import numpy as np
from scipy import stats


class Parameter:
    '''Class for handling a Parameter'''

    def __init__(self, value):
        '''Parameter initialization'''
        self.value = value

    def set(self, value):
        '''set/change Parameter'''
        self.value = value

    def __call__(self):
        '''recall Parameter'''
        return self.value


def fit(function, parameters, data, x=None):
    '''performing a least square fit using function with parameters (list of Parameter instances). Data is y. x is optional.'''

    def f(params):
        '''return difference between function and data'''
        i = 0
        for p in parameters:
            p.set(params[i])
            i += 1
        return data - function(x)

    if x is None:
        x = np.arange(data.shape[0])

    p = [param() for param in parameters]

    optimize.leastsq(f, p)

    return function(x), x


def fit2d(function, parameters, data, x=None, y=None):
    '''performing a least square fit using function with parameters (list of Parameter instances). Data is y. x is optional.'''

    def f(params):
        '''return difference between function and data'''
        i = 0
        for p in parameters:
            p.set(params[i])
            i += 1
        result = data - function(x, y)
        return result.flatten()

    if x is None and y is None: x, y = np.meshgrid(np.arange(data.shape[0]), np.arange(data.shape[1]))
    if x is None: x, dummy = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[1]))
    if y is None: dummy, y = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[1]))

    p = [param() for param in parameters]

    optimize.leastsq(f, p)

    return function(x, y), x, y


# =================================================================

def Gauss(center, height, FWHM, offset):
    '''
	returns a 1d gauss function f(x) with given parameters (instances of class Parameter)
	center - center of gauss
	height - height of gauss at center
	FWHM   - full width at half maximum
	offset - y offset of gauss
	'''

    def f(x):
        return offset() + height() * np.exp(-((x - center()) * 2 * np.log(2) / FWHM()) ** 2)

    return f


def Normal(center, area, sigma, offset):
    '''
	returns a 1d gauss function f(x) with given parameters (instances of class Parameter)
	center - center of gauss
	area   - whole area under gauss
	sigma  - standart deviation of gaussian distribution
	offset - y offset of gauss
	'''

    def f(x):
        return offset() + area() / sigma() / np.sqrt(2 * np.pi) * np.exp(-1 / 2 * ((x - center()) / sigma()) ** 2)

    return f


def Two_Gauss(center1, center2, area1, area2):
    '''
	returns a 1d function f(x) that's made of two gaussians/poissonians, parameters (instances of class Parameter)
	center1 - center of gauss 1
	area1   - whole area under gauss 1
	center2 - center of gauss 2
	area2   - whole area under gauss 2
	'''

    def f(x):
        return area1() / np.sqrt(center1()) / np.sqrt(2 * np.pi) * np.exp(-1 / 2 * ((x - center1()) / np.sqrt(center1())) ** 2) + area2() / np.sqrt(center2()) / np.sqrt(2 * np.pi) * np.exp(-1 / 2 * ((x - center2()) / np.sqrt(center2())) ** 2)

    return f


def Poisson(avg, area):
    def f(n): return area() * stats.distributions.poisson.pmf(n, avg())

    return f


def Two_Pois(avg1, area1, avg2, area2):
    f1 = Poisson(avg1, area1)
    f2 = Poisson(avg2, area2)

    def f(n):
        return (f1(n) + f2(n))

    return f


def Gauss_2D_1(x0, y0, height, FWHM, offset):
    '''
	resturns a 2D Gauss, symmetric
	x0
	y0
	height
	FWHM   - full width half maximum for all directions
	offset
	'''
    param0 = Parameter(0)
    param1 = Parameter(1)
    Gauss_x = Gauss(x0, param1, FWHM, param0)
    Gauss_y = Gauss(y0, param1, FWHM, param0)

    def f(x, y):
        return height() * (Gauss_x(x) * Gauss_y(y)) + offset()

    return f


def Gauss_2D_2(x0, y0, height, wx, wy, offset):
    '''
	resturns a 2D Gauss, with 2 axis (x,y)
	x0
	y0
	height
	wx     - full width half maximum in x direction
	wy     - full width half maximum in y direction
	offset
	'''
    param0 = Parameter(0)
    param1 = Parameter(1)
    Gauss_x = Gauss(x0, param1, wx, param0)
    Gauss_y = Gauss(y0, param1, wy, param0)

    def f(x, y):
        return height() * (Gauss_x(x) * Gauss_y(y)) + offset()

    return f


def Gauss_2D_3(x0, y0, height, w1, w2, theta, offset):
    '''
	resturns a 2D Gauss, with 2 axis, rotated
	x0
	y0
	height
	w1     - full width half maximum in 1 direction
	w2     - full width half maximum in 2 direction
	theta  - tilt of ellipsoidal gauss
	offset
	'''
    param0 = Parameter(0)
    param1 = Parameter(1)
    Gauss_1 = Gauss(x0, param1, w1, param0)
    Gauss_2 = Gauss(y0, param1, w2, param0)

    def f(X, Y):
        x = X.copy()
        y = Y.copy()
        x -= x0()
        y -= y0()
        x, y = x * np.cos(theta()) - y * np.sin(theta()), x * np.sin(theta()) + y * np.cos(theta())
        x += x0()
        y += y0()
        return height() * (Gauss_1(x) * Gauss_2(y)) + offset()

    return f
