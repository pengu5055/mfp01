import matplotlib.pyplot as plt
from mpmath import *
import numpy as np


def f(x, n):
    c = mpf(0)  # Counter
    for k in range(0, n+1):
        step = rf(1/3, k) * mpf("3")**mpf(k) * mpf(x)**mpf(3*k)/fac(3*k)
        c += step
    return c


def f_arr(array, n):
    return np.array([f(mpf(x), n) for x in array])


def arrayize(array, function, *args,):
    return np.array([function(x, *args) for x in array])



x_lin = np.linspace(0, 1, 20)
print(f_arr(x_lin, 10))
print(arrayize(x_lin, f, 10))
