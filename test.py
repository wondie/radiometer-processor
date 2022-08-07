'''Import Python math library'''
import math
import matplotlib.pyplot as plt
'''try and except ImportError handler will be printing message if you haven't install required python library'''
try:
    import scipy.optimize as opt
    import scipy
    from scipy.optimize import minimize
except ImportError:
    print ("You must install Python scipy first before running this program")

try:
    import numpy as np
except ImportError:
    print ("You must install Python numpy first before running this program")

'''function f(x) = sin(x)'''
def f(x):
    print (x)
    #return -x**(2)
    return np.sin(x)

'''
check here for fmin_l_bgfs_b docs : https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fmin_l_bfgs_b.html#scipy.optimize.fmin_l_bfgs_b
and what params need to used.

or you can take a look at here : https://docs.scipy.org/doc/scipy/reference/optimize.html for other method and algorithm
'''
initial_guess = 1.0

minbound = -9
maxbound = 9

max_bounds_area = (minbound,maxbound)
max_x = opt.fmin_l_bfgs_b(lambda x: -f(x), initial_guess, bounds=[(max_bounds_area)],approx_grad=True)

# I want to plot the result. I try this but I get TypeError: only length-1 arrays can be converted to Python scalars
t = np.arange(minbound, maxbound, initial_guess)
s = f(t)
plt.plot(t, s)
plt.show()