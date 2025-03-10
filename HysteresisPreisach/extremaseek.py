import numpy as np
import matplotlib.pyplot as plt
#import time

def fun(x):
    y= -1010* np.cos(5*x) * np.cos(3*x) /(x+1)
    return y

#-----------------------

def localmaxima(u):
    if len(u)<3: 
        lex=np.arange(len(u))
    else: #maksimum je ako je >= od prethodnog i iskljucivo veci od iduceg
        lex= np.nonzero((u[1:-1]>=u[0:-2]) & (u[1:-1]>u[2:]))[0] + 1 #np.nonzero() vraca tuple of arrays pa zato [0]
    
        if (u[0]>u[1]) : lex = np.concatenate(([0], lex))
        n=len(u)
        if (u[-1]>=u[-2]): lex = np.concatenate((lex, [n-1]))
    return lex


#------------------------

#pronalazi niz dominantnih ekstrema
def dominant_extrema(u, ekstremi):
    #pretpostavljamo da je len(ekstremi)>2, tj imamo barem tri lokalna ekstrema
    domEkstremi=[]  #inicijalizacija liste za dominantne ekstreme
    domEkstremi.append(ekstremi[-1])  #posljednji lokalni ekstrem je uvijek i dominantni ekstrem
    #if (u[ekstremi[-1]] == u_max): return domEkstremi

    globExtrInd=np.nonzero(np.abs(u[ekstremi])==np.max(np.abs(u[ekstremi])))[0][-1]
    ekstremi=ekstremi[globExtrInd:]
    
    if (len(ekstremi)==1): return domEkstremi

    flag = u[ekstremi[-2]] < u[domEkstremi[-1]]  #flag=True if domEkstrema[-1] je lokalni maksimum i obratno
    for i in range(len(ekstremi)-2,-1,-1):
        if (flag and (u[ekstremi[i]] > u[domEkstremi[-1]])):  #if True, pronadji dominantni minimum
            #domEkstremi.append(ekstremi[i] + np.argmin(u[ekstremi[i]:domEkstremi[-1]])) #daje krivi index kod konst funkcije, al tocan iznos ekstrema
            domindex=ekstremi[i] + np.nonzero(u[ekstremi[i]:domEkstremi[-1]]==np.min(u[ekstremi[i]:domEkstremi[-1]]))[0][-1]
            domEkstremi.append(domindex)
            flag=False
        elif (not flag and (u[ekstremi[i]] < u[domEkstremi[-1]])):  #if True, pronadji dom maksimum
            #domEkstremi.append(ekstremi[i] + np.argmax(u[ekstremi[i]:domEkstremi[-1]]))
            domindex=ekstremi[i] + np.nonzero(u[ekstremi[i]:domEkstremi[-1]]==np.max(u[ekstremi[i]:domEkstremi[-1]]))[0][-1]
            domEkstremi.append(domindex)
            flag=True
    #...ova petlja ne moze appendat najstariji dominantni maksimum/minimum M0 pa moramo jos njega pronaci i provjeriti je li on globalni ekstrem
    #Cim smo u ovoj funkciji znamo da je len(ekstremi)>2, pa zbog prethodne petlje znamo da postoji \
    #...dominantni maksimum/minimum M0 stariji od u[domEkstremi[-1]] al ne znamo je li on ujedno i dominantni ekstrem jer svojstvo najstarijeg
    #dominantnog ekstrema je to da je on ujedno i globalni ekstrem!
    #Dakle, da bi M0 appendali u domEkstremi, mora vrijediti da je abs(M0)>abs(u[domEkstremi])
    #ono sto sigurno znamo je da postoji ekstrem E0 takav da je abs(E0) > abs(domEkstremi[-2]), ako postoji abs(domEkstremi[-2])
    #tj. ako je len(domEkstremi)>1
    #domEkstremi[-1] je trenutno vremenski najstariji dominantni ekstrem
    #ako je domEkstremi[-1] maksimum, moguce je da prebrise najstariji dominantni ekstrem (minimum) i obratno.
    #...to je moguce u slucaju da je np.abs(u[domEkstremi[-1]]) > od apsolutne vrijendosti najstarijeg dominantnog ekstrema

    domEkstremi.append(ekstremi[0])

    domEkstremi.reverse()
    return domEkstremi


#---------------------------

#sortiranje minimuma i maksimuma u listu ekstrema::
def localextrema(u):
    ekstremi=[]
    if len(u)<3: 
        ekstremi = localmaxima(u) #jer u ima samo dva clana
    #    print('ekstremi=',ekstremi)
        return ekstremi
    
    max_indeksi = localmaxima(u)
    min_indeksi = localmaxima(-u)
    #print('max_indeksi=', max_indeksi)
    #print('min_indeksi=', min_indeksi)
    if (min_indeksi[0]<=max_indeksi[0]):
        for i in range(len(max_indeksi)+len(min_indeksi)):
            if i%2==0:
                ekstremi.append(min_indeksi[i//2])
            else: 
                ekstremi.append(max_indeksi[i//2])
    else:
        for i in range(len(max_indeksi)+len(min_indeksi)):
            if i%2==0:
                ekstremi.append(max_indeksi[i//2])
            else: 
                ekstremi.append(min_indeksi[i//2])
    #print('ekstremi=',ekstremi)
    return ekstremi




""" def DominantExtrema(u, max_indeksi, min_indeksi):
    domMax_ind=[]
    domMin_ind=[]
    if max_indeksi[-1]>min_indeksi[-1]:
        Nmax=len(max_indeksi)
        domMax_ind.append(max_indeksi[-1]) #u(t) je sigurni clan
        domMin_ind.append(max_indeksi[-1]) #privremeni clan
        j=-1
        for i in range(Nmax-2,-1,-1):
            if (u[max_indeksi[i]] > u[max_indeksi[j]]):
                pMinInd=max_indeksi[i] + \
                              np.argmin( u[max_indeksi[i] : max_indeksi[j]] )
                if (u[pMinInd] < u[domMin_ind[-1]]):
                    domMin_ind.append(pMinInd)
                    domMax_ind.append(max_indeksi[i])
                    j=i
                else:
                    domMax_ind.pop()
                    domMax_ind.append(max_indeksi[i])    
        domMax_ind.reverse()
        domMin_ind.reverse()
        domMin_ind.pop()
    
    else:
        Nmin=len(min_indeksi)
        domMin_ind.append(min_indeksi[-1])
        domMax_ind.append(min_indeksi[-1])
        j=-1
        for i in range(Nmin-2,-1,-1):
            if (u[min_indeksi[i]] < u[min_indeksi[j]]):
                pMaxInd=min_indeksi[i] + \
                              np.argmin( u[min_indeksi[i] : min_indeksi[j]] )
                if (u[pMaxInd] > u[domMax_ind[-1]]):
                    domMax_ind.append(pMaxInd)
                    domMin_ind.append(min_indeksi[i])
                    j=i
                else:
                    domMin_ind.pop()
                    domMin_ind.append(min_indeksi[i])
                    
        domMax_ind.reverse()
        domMin_ind.reverse()
        domMax_ind.pop()
    
    return domMax_ind, domMin_ind """



#----------------------------

def FindDominantExtrema(u,u_max):

    ekstremi=localextrema(u)
    
    if abs(u[ekstremi[-1]])==abs(u_max):
        domindex=[ekstremi[-1]] #lista s jednim clanom
    elif len(ekstremi)==1:
        domindex = ekstremi
    else:
        domindex = dominant_extrema(u,ekstremi)
    #print('domindex=', domindex)

    return domindex



#===============
#  P O Z I V   
#===============

if __name__ == "__main__":

    x= np.linspace(0,9,251)
    H=fun(x)

    Hmax=430
    H=np.minimum(H,Hmax)
    H=np.maximum(H,-Hmax)
    #H0=np.minimum(H,Hmax)
    #H0=np.maximum(H0,-Hmax)

    
    dominantni = FindDominantExtrema(H,Hmax)
    #dominantni = FindDominantExtrema(H0,Hmax)

    plt.plot(x,H,label='graf u(t)')
    plt.plot(x[dominantni], H[dominantni], marker='o', linestyle='',label='dominantni ekstremi')
    #plt.plot(x[dominantni], H0[dominantni], marker='o', linestyle='')

    x_ticks = np.linspace(min(x), max(x), 10) 
    y_ticks = np.linspace(-420, 420, 11)  
    plt.xticks(x_ticks)
    plt.yticks(y_ticks)

    plt.xlabel('Vrijeme t')  # Dodajte opis osi x
    plt.ylabel('Ulazni signal u')  # Dodajte opis osi y
    plt.legend()
    plt.grid()
    plt.show()



""" 
#pronalazi niz dominantnih ekstrema
def dominant_extrema(u, ekstremi,u_max):
    #pretpostavljamo da je len(ekstremi)>2, tj imamo barem tri lokalna ekstrema
    domEkstremi=[]  #inicijalizacija liste za dominantne ekstreme
    domEkstremi.append(ekstremi[-1])  #posljednji lokalni ekstrem je uvijek i dominantni ekstrem
    #if (u[ekstremi[-1]] == u_max): return domEkstremi

    globExtrInd=np.nonzero(np.abs(u[ekstremi])==np.max(np.abs(u[ekstremi])))[0][-1]
   
    flag = u[ekstremi[-2]] < u[domEkstremi[-1]]  #flag=True if domEkstrema[-1] je lokalni maksimum i obratno
    for i in range(len(ekstremi)-2,-1,-1):
        if (flag and (u[ekstremi[i]] > u[domEkstremi[-1]])):
            #domEkstremi.append(ekstremi[i] + np.argmin(u[ekstremi[i]:domEkstremi[-1]])) #daje krivi index kod konst funkcije, al tocan iznos ekstrema
            domindex=ekstremi[i] + np.nonzero(u[ekstremi[i]:domEkstremi[-1]]==np.min(u[ekstremi[i]:domEkstremi[-1]]))[0][-1]
            domEkstremi.append(domindex)
            flag=False
        elif (not flag and (u[ekstremi[i]] < u[domEkstremi[-1]])):
            #domEkstremi.append(ekstremi[i] + np.argmax(u[ekstremi[i]:domEkstremi[-1]]))
            domindex=ekstremi[i] + np.nonzero(u[ekstremi[i]:domEkstremi[-1]]==np.max(u[ekstremi[i]:domEkstremi[-1]]))[0][-1]
            domEkstremi.append(domindex)
            flag=True
    #...ova petlja ne moze appendat najstariji dominantni maksimum/minimum M0 pa moramo jos njega pronaci i provjeriti je li on globalni ekstrem
    #Cim smo u ovoj funkciji znamo da je len(ekstremi)>2, pa zbog prethodne petlje znamo da postoji \
    #...dominantni maksimum/minimum M0 stariji od u[domEkstremi[-1]] al ne znamo je li on ujedno i dominantni ekstrem jer svojstvo najstarijeg
    #dominantnog ekstrema je to da je on ujedno i globalni ekstrem!
    #Dakle, da bi M0 appendali u domEkstremi, mora vrijediti da je abs(M0)>abs(u[domEkstremi])
    #ono sto sigurno znamo je da postoji ekstrem E0 takav da je abs(E0) > abs(domEkstremi[-2]), ako postoji abs(domEkstremi[-2])
    #tj. ako je len(domEkstremi)>1
    #domEkstremi[-1] je trenutno vremenski najstariji dominantni ekstrem
    #ako je domEkstremi[-1] maksimum, moguce je da prebrise najstariji dominantni ekstrem (minimum) i obratno.
    #...to je moguce u slucaju da je np.abs(u[domEkstremi[-1]]) > od apsolutne vrijendosti najstarijeg dominantnog ekstrema


    olderEks=u[ekstremi[:ekstremi.index(domEkstremi[-1])]] #ovo su svi preostali lokalni ekstremi stariji od u[domEkstremi[-1]]

    if (np.abs(u[domEkstremi[-1]])>=np.max(np.abs(olderEks))):  #u ovom slucaju svi stariji ekstremi su prebrisani pa izlazimo iz funkcije
        domEkstremi.reverse()  
        return domEkstremi

    #if (np.abs(u[domEkstremi[-1]])==u_max):  #u ovom slucaju svi stariji ekstremi su prebrisani pa izlazimo iz funkcije
    #    domEkstremi.reverse()  
    #    return domEkstremi

    if (flag and len(domEkstremi)==1):
        domEkstremi.append( np.nonzero(olderEks == np.min(olderEks))[0][-1] )
    elif (not flag and len(domEkstremi)==1):     
        olderEks=u(ekstremi[:ekstremi.index(domEkstremi[-1])])
        domEkstremi.append( np.nonzero(olderEks == np.max(olderEks))[0][-1] )
    
    print('domEkstremi',domEkstremi)
    if(flag and any(u[:domEkstremi[-1]]<u[domEkstremi[-2]])): 
        domEkstremi.append( np.nonzero(u[:domEkstremi[-1]] == np.min(u[:domEkstremi[-1]]))[0][-1] )
    elif(not flag and any(u[:domEkstremi[-1]]>u[domEkstremi[-2]])): 
        domEkstremi.append( np.nonzero(u[:domEkstremi[-1]] == np.max(u[:domEkstremi[-1]]))[0][-1] )
    domEkstremi.reverse()
    return domEkstremi


#--------------------------- """