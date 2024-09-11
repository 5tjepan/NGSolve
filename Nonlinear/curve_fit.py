import numpy as np
from scipy.optimize import curve_fit


B=np.array([1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.74,1.77,1.8,1.83])
H=np.array([3.33148997, 3.38979934, 3.46729718, 3.56473246, 3.69037833,
       3.8537581 , 4.06645936, 4.34626991, 4.67282883, 5.07517382,
       5.59842196, 6.2146081 , 6.90775528])

#H=np.array([27.98,29.66,32.05,35.33,40.06,47.17,58.35,77.19,107,160,270,500,1000])

#B=np.array([1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.74,1.77])
#H=np.array([27.98,29.66,32.05,35.33,40.06,47.17,58.35,77.19,107,160,270])


""" def fun(B,a,b,c):
    return a*B/(B*b+c) """

def fun(B,a,b,c,d):
    #return a + b*B + c*B**2 + d*B**3
    e=2.718281828459
    return e**(a + b*B + c*B**2 + d*B**3)

bounds = ([-100, -2,-5], [100, 2,5])
initial_guess = [-14, 1.865,1]
params, covariance = curve_fit(fun, B, H)
#params, covariance = curve_fit(fun, B, H, p0=None, bounds=None)

a, b, c, d =params

Hfit=fun(B,a,b,c,d)

import matplotlib.pyplot as plt
plt.scatter(B, H)
plt.plot(B, Hfit)
plt.xlabel("B")
plt.ylabel("H")
plt.show()


a=-144.54207430769281
b= 304.6915669443657
c=-209.76353623263986
d=48.34162716175372

#e = np.e
e=2.718281828459
e**(a + b*Babs + c*Babs**2 + d*Babs**3)