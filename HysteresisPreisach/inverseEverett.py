import numpy as np
import matplotlib.pyplot as plt


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


def invEverett2(a,b):
    y= (-2/0.03)*np.log( -0.25 + (1/(a+2.0)) + (1/(-b+2.0)) ) - (-2/0.03)*np.log( -0.25 + (1/(-a+2.0)) + (1/(b+2.0)) )
    return y

def invEverett(a,b):
    c=0.25
    d=1
    y= (-2/0.03)*np.log( c+ (d/(a+2.0)) + (d/(-b+2.0)) ) - (-2/0.03)*np.log( c + (d/(-a+2.0)) + (d/(b+2.0)) )
    return y

if __name__=='__main__':

    x=np.linspace(-1.99,1.99,100)
    y=np.linspace(-1.99,1.99,100)
    a,b = np.meshgrid(x,y)
    invEt=invEverett(a,b)


    def fun(x,a,b):
        y = np.log(( -b + (1/(x+2.0)) )**(-1/a) )
        return y

    x=np.linspace(-1.999,1.999,100)
    y1=fun(x,0.03,0.25)
    #plt.plot(x,y1,color='blue')
    #plt.show()
    #print(y1)

    plot_contour_2dplane(b,a,invEt)


""" def invEverett(a,b,c,d):
    y= (-2/0.03)*np.log( c+ (d/(a+2.0)) + (d/(-b+2.0)) ) - (-2/0.03)*np.log( c + (d/(-a+2.0)) + (d/(b+2.0)) )
    return y
 """