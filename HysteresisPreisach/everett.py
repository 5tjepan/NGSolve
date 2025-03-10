import numpy as np
import matplotlib.pyplot as plt


def fun(x,a,b,c,d,e):
    #y=np.arctan(a*x)
    y= c/(b+ np.exp(a*x + e)) + d
    return y

def fun2(x,a,b,c,d,e):
    #y=np.arctan(a*x)
    y= c/(0+ 4*np.sinh(e)*np.cosh(a*x)) + d
    return y


def plot_BHloop():
    x=np.linspace(-200,200,10)
    #y=np.sin(4*x)
    #y1=fun(x,0.03,4.5,9,-1,0)

    y1=fun(x,-0.03,1,   2,-1,-1.5)
    y2=fun(x,-0.03,1,   2,-1, 1.5)
    
    y3=fun(x,-0.03,1.2, 2,-1,-1.5)
    y4=fun(x, 0.03,1.2,-2, 1,-1.5)
    
    y5=fun(x, 0.03,1.4,-2, 1,-1.5)
    #print('y3',y3)
    #print('y5',y5)
    #print('x',x)
    
    plt.plot(x,y1,color='blue')
    plt.plot(x,y2,color='red')
    plt.plot(x,y3,color='blue',linestyle='--')
    plt.plot(x,y4,color='red', linestyle='--')
    plt.plot(x,y5,color='green', linestyle='--')
    plt.grid()

    plt.axhline(0, color='black', linewidth=0.8, linestyle='--')  # Oznaka x-osi
    plt.axvline(0, color='black', linewidth=0.8, linestyle='--')  # Oznaka y-osi

    plt.xlim(-280, 280)
    plt.ylim(-1.6, 1.6)

    plt.show()


def mu1(a,b):
    y=(0.00171038 *np.exp(-0.03*a - 0.03*b))/(1 + 4.25856 *np.exp(-0.03*a) + \
                                              0.22313* np.exp(-0.03* b) )**3
    return y

def mu(a,b):
    y= - (0.0036 *np.exp(0.03*a - 0.03*b))/(1 + np.exp(0.03*a) + np.exp(-0.03*b))**3 + \
         (0.0036 *np.exp(-0.03*a + 0.03*b))/(1 + np.exp(-0.03*a) + np.exp(0.03*b))**3
    return y


def plot_contour_2dplane(b,a, F):
    
    # Konturirani graf
    plt.contour(b,a, F, levels=20, cmap='viridis')
    plt.colorbar()  # Bočna traka za boje

    plt.title('Contour Plot of Preisach function P=mu(a, b)')
    plt.xlabel('b')
    plt.ylabel('a')

    plt.axhline(0, color='black', linestyle='--', linewidth=1)  
    plt.axvline(0, color='black', linestyle='--', linewidth=1)

    plt.show()

def Everett1(a,b):
    y= 2/(1 + np.exp(-0.03*a + 1.5)) - \
          2/(1 + 2*np.sinh(1.5)*np.exp(-0.03*a) + np.exp(-0.03*b - 1.5))
    return y

def Everett2(a,b):
    y= -2/(1+ np.exp(0.03*a - 1.5)) + 1 - \
          2/(1 + 2*np.sinh(1.5)*np.exp(-0.03*a) + np.exp(-0.03*b - 1.5)) + 1
    return y

def Everett_exp(a,b):
    c=0.25
    y= 1/(c + np.exp(-0.03*a) + np.exp(+0.03*b) ) - \
          1/(c + np.exp(+0.03*a) + np.exp(-0.03*b) )
    return y

def Everett_exp2(a,b):
    y= 1/(0.25 + np.exp(-0.03*a) + np.exp(+0.03*b) + np.exp(-0.04*a) + np.exp(+0.04*b)) - \
          1/(0.25 + np.exp(+0.03*a) + np.exp(-0.03*b) +np.exp(+0.04*a) + np.exp(-0.04*b))
    return y

def Everett_atan(x,y):
    a=0.0196483
    b=2.95329554
    c=0.02211744
    d=1.04359946
    
    valid= x>=y

    alpha=valid*x
    beta=valid*y

    #E= (np.arctan(a*x) - np.arctan(a*y))**b + (np.arctan(c*x)**3 - np.arctan(c*y)**3)**d
    E= 0.1*((np.arctan(a*alpha) - np.arctan(a*beta))**b + (np.arctan(c*alpha)**3 - np.arctan(c*beta)**3)**d)
    return E


if __name__=='__main__':

    x=np.linspace(-400,400,2000)
    y=np.linspace(-400,400,2000)

    a,b = np.meshgrid(x,y)

    P=mu1(a,b)
    E1=Everett1(a,b)
    Et=Everett_exp(a,b)
    Ea=Everett_atan(a,b)

    #plt.plot(x, Everett_atan(x,np.full_like(y,1000)))
    #plt.plot(x, Everett_atan(x,np.full_like(y,500)))
    #plt.plot(x, Everett_atan(x,np.full_like(y,100)))

    #plt.plot(x,2/(1 + np.exp(0.03*x) + np.exp(-0.03*x)))
    plt.show()


    #plot_BHloop()
    #plot_contour_2dplane(b,a,Ea)
    #plot_contour_2dplane(b,a,E1)
    plot_contour_2dplane(b,a,Et)
    #plot_contour_2dplane(b,a,P)

    print('everett=', Everett_atan(100,-100))


