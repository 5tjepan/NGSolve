import numpy as np
import matplotlib.pyplot as plt
from extremaseek import FindDominantExtrema
from everett import Everett_atan,Everett_exp, Everett_exp2
from inverseEverett import invEverett

# i n p u t funkcija:
#--------------
def fun(t):
    y= -1000* np.cos(5*t) * np.cos(3*t) /(t+1)
    return y


Everett = Everett_atan
#---------------
def preisach_output(u, domindex, u_max):
    domex=u[domindex]

    if np.abs(domex[-1]) == u_max: 
        output= 0.5* np.sign(domex[-1]) * Everett(u_max,-u_max)
        return output

    output= np.sign(domex[0]) * 0.5*Everett(np.abs(domex[0]),-np.abs(domex[0]))
    #print('domex=',domex[0])
    #print('----->', output)
    if len(domex)>1: 
        flag=domex[1]>domex[0] #domex[0] je globalni ekstrem i zanima nas je li globalni min ili globalni max
        if len(domex)%2==0: domex=np.append(domex,domex[-1]) #appendanjem domex[-1] na domex ne mijenja se output jer je Everett(a,a)=0
        if flag:
            for i in range(1,len(domex),2):
                output+= Everett(domex[i],domex[i-1]) - Everett(domex[i], domex[i+1])
                #print('output2',output)
        else:
            for i in range(1,len(domex),2):
                output+= -Everett(domex[i-1],domex[i]) + Everett(domex[i+1], domex[i])
                #print('output3=',output)
    return output    
#-----------------


#==================
#  P O Z I V
#==================

if __name__ == "__main__":

    t= np.linspace(0,2.0,400)
    H=fun(t)
    #H=39*np.sin(6.2832*t)
    #H= -300*np.cos(6.2832*t) /(0.5*t+1)

    #print(H)

    Hmax=400
    H=np.minimum(H,Hmax)
    H=np.maximum(H,-Hmax)

    #---------------

    Everett = Everett_atan

    dominantni=FindDominantExtrema(H,Hmax)
    print('B =',preisach_output(H,dominantni,Hmax))

    B=[]
    for i in range(1,len(H)+1):
        domIndeksi=FindDominantExtrema(H[:i],Hmax)
        #print('domIndeksi',domIndeksi)
        B.append(preisach_output(H[:i],domIndeksi,Hmax))

    #print('B=',B)
    Bscaled=[B[i]*100 for i in range(len(B))]
    plt.plot(t,H, label=r'$H(t)$')
    plt.plot(t[dominantni], H[dominantni], marker='o', linestyle='')
    plt.plot(t,Bscaled, label=r'$B(t)$')
    plt.xlabel(r'$t$', fontsize=15)
    plt.ylabel(r'$f(t)$', fontsize=15)
    plt.legend()
    plt.grid()
    plt.show()

    plt.plot(H,B)
    plt.grid()
    plt.axhline(0, color='black', linewidth=0.8, linestyle='--')  # Oznaka x-osi
    plt.axvline(0, color='black', linewidth=0.8, linestyle='--')  # Oznaka y-osi
    plt.xlabel(r'$H$', fontsize=15)
    plt.ylabel(r'$B$', fontsize=15)
    plt.show()





#:::::::::::::::::::::
""" t= np.linspace(0,1.0,100)
#I=fun(t)
#I=1.99*np.sin(6.2832*t)
#I= -1.5*np.cos(6.2832*t) /(0.2*t+1)

I=np.array(B)
Imax=1.99
I=np.minimum(I,Imax)
I=np.maximum(I,-Imax)

Everett = invEverett

#M is magnetization
#I is mag polarization

M=[]
for i in range(1,len(I)+1):
    domIndeksi=FindDominantExtrema(I[:i],Imax)
    #print('domIndeksi',domIndeksi)
    M.append(preisach_output(I[:i],domIndeksi,Imax))

Mscaled=[M[i]/100 for i in range(len(M))]
plt.plot(t,I)

dominantni=FindDominantExtrema(I,Imax)

plt.plot(t[dominantni], I[dominantni], marker='o', linestyle='')
plt.plot(t,Mscaled)
plt.show()

plt.plot(I,M)
plt.grid()
plt.axhline(0, color='black', linewidth=0.8, linestyle='--')  # Oznaka x-osi
plt.axvline(0, color='black', linewidth=0.8, linestyle='--')  # Oznaka y-osi
plt.xlabel('I')
plt.ylabel('M')
plt.show() """




