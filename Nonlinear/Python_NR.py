import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
#from scipy.misc import derivative

xos=np.linspace(0,4,20)
y=xos**2 - 1
f=interp1d(xos, y, kind='cubic')

dydx = np.gradient(y, xos)
dfdx = interp1d(xos,dydx, kind='cubic')
#dydx=derivative(f, xos, dx=1e-5)

#plt.plot(xos,y)
plt.xlabel('xos')
plt.ylabel('funkcija')
plt.title('graf $y=x^2-1$')
plt.grid(True)

#plt.subplot(2,1,1)
plt.plot(xos,y,linestyle='', marker='o')
plt.plot(xos,f(xos))

#plt.subplot(2,1,2)
plt.plot(xos,dydx, linestyle='', marker='*', markersize=10)
plt.plot(xos,dfdx(xos))

#plt.show()


x0=4*np.random.rand()
print(f'(x0, dfdx) je ({x0},{dfdx(x0)})')
xi=x0
tol=1e-5
for i in range(1,10):
    print(f"u {i}-toj iteraciji xi={xi}")
    xj=-f(xi)/(dfdx(xi)) + xi
    if ((abs(xi-xj)/xi) < tol):
        break
    xi=xj

print(f'nultocka funkcije f je x={xj}')