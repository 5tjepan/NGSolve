import numpy as np
import matplotlib.pyplot as plt
from extremaseek import FindDominantExtrema
from everett import Everett_atan,Everett_exp, Everett_exp2
from inverseEverett import invEverett

# i n p u t funkcija:
#--------------
def fun(t):
    #y= 100*t 
    y=1.8*(100*np.sin(6.28*t) - 20*np.sin(3*6.28*(t-0)) - 35*np.sin(5*6.28*(t-0)) ) #+20 #*(1+np.cos(0.05*6.28*t))
    #y= 100*np.sin(4*t) - 15*np.sin(3*4*t-2) + 35*np.sin(5*4*t-5) - 0
    #y= 100*(1*np.sin(8.5*t)+3*np.sin(t)-2) #ovaj signal pokazuje presjecanje uzlaznih grana malih histereznih krivulja...potvrdjuje hystory-dependance
    #y= 100*(np.sin(18*t)+4*np.sin(t)-2)
    #y= -100000* np.cos(4*t) * np.cos(90*t) *np.sin(t)/(t+4) *(t-0.25)**2
    #y= -1000* np.cos(5*t) * np.cos(3*t) /(t+1)
    return y


Everett = Everett_exp #Everett_atan
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

    t= np.linspace(0,2,1500)
    H=fun(t)
    #H=39*np.sin(6.2832*t)
    #H= -300*np.cos(6.2832*t) /(0.5*t+1)

    #print(H)

    Hmax=500
    H=np.minimum(H,Hmax)
    H=np.maximum(H,-Hmax)

    #---------------

    Everett = Everett_exp #Everett_atan

    dominantni=FindDominantExtrema(H,Hmax)
    print('B =',preisach_output(H,dominantni,Hmax))

    B=[]
    for i in range(1,len(H)+1):
        domIndeksi=FindDominantExtrema(H[:i],Hmax)
        #print('domIndeksi',domIndeksi)
        B.append(preisach_output(H[:i],domIndeksi,Hmax))

    #print('B=',B)
    """Bscaled=[B[i]*100 for i in range(len(B))]
    plt.plot(t,H, linewidth=3, label=r'$H(t)$')
    #plt.plot(t[dominantni], H[dominantni], marker='o', linestyle='') #UGASIO SAM MARKERE DOMINANTNIH EKSTREMA
    plt.plot(t,Bscaled, linewidth=3, label=r'$B(t)$')
    plt.xlabel(r'$t$', fontsize=20)
    plt.ylabel(r'$f(t)$', fontsize=20)
    plt.legend()
    plt.grid()
    plt.show() """
    

    # A. Postavljanje fonta za standardni tekst
    plt.rcParams["font.family"] = "serif"
    # Navođenje Times New Romana kao prvog izbora unutar serifne obitelji
    plt.rcParams["font.serif"] = ["Times New Roman", "Times"]

    # B. Postavljanje fonta za matematičke izraze (ako koristite LaTeX stil r'$...$')
    # Ovo je ključno za ujednačavanje fonta u svim labelama
    plt.rcParams["mathtext.fontset"] = "custom"
    plt.rcParams["mathtext.rm"] = "Times New Roman" 
    # Opcionalno, za talike i podebljane simbole:
    plt.rcParams["mathtext.it"] = "Times New Roman:italic"
    plt.rcParams["mathtext.bf"] = "Times New Roman:bold"



    fig, ax1 = plt.subplots()

    # --- Lijeva y-osa (H) ---
    ax1.plot(t, H, color='tab:blue', linewidth=2.5, linestyle = '--', label=r'$H(t)$')
    ax1.set_xlabel(r'$t$, (s)', fontsize=22)
    ax1.set_ylabel(r'$H$, (A/m)', fontsize=22)
    ax1.tick_params(axis='both', labelsize=20)

    # --- Desna y-osa (B) ---
    ax2 = ax1.twinx()  # kreira drugu y-osu koja dijeli istu x-os
    ax2.plot(t, B, color='tab:red', linewidth=2.5, linestyle = '-', label=r'$B(t)$')
    ax2.set_ylabel(r'$B$, (T)', fontsize=22)
    ax2.tick_params(axis='both', labelsize=20)
    ax2.set_ylim(-2.32, 2.32)

    # --- Ostalo ---
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper right', fontsize=18, handlelength=1.0, handletextpad=0.5, labelspacing=0.1)

    fig.tight_layout()
    ax1.grid(True)
    plt.show()

    plt.plot(H,B, linewidth=2.5)
    plt.grid()
    plt.axhline(0, color='black', linewidth=0.8, linestyle='--')  # Oznaka x-osi
    plt.axvline(0, color='black', linewidth=0.8, linestyle='--')  # Oznaka y-osi
    plt.xlabel(r'$H$, (A/m)', fontsize=22)
    plt.ylabel(r'$B$, (T)', fontsize=22)

    plt.tick_params(
    axis='both',       # Primijeni na obje osi (x i y)
    which='major',     # Primijeni na glavne tikove
    labelsize=20)       # Postavite željenu veličinu fonta
    
    plt.ylim(-2.32, 2.32)
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




