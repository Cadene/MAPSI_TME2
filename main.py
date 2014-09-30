# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import pickle as pkl
import random as r
import math as m
from mpl_toolkits.mplot3d import Axes3D
plt.close('all')    



def bernoulli(p):
    x = r.random()
    if x<=p :
        return 1
    else :
        return 0
        
#print bernoulli(0.5)

def binomiale (n,p) :
    cpt = 0
    for i in range(n):
        if bernoulli(p) > 0:
           cpt+=1
    return cpt
    
    
#print binomiale(10,0.9)


def Galton(n, billes):
    tab=np.zeros(billes, int)
    tab_possibilite = np.zeros(n, int)
    intervalle = 0;
    for i in range(billes):
        tab[i]=binomiale(n,0.5)
        tab_possibilite[tab[i]]+=1
    for i in range (n): 
        if tab_possibilite[i] > 0:
            intervalle += 1
    #print tab_possibilite
    #print intervalle
    plt.hist(tab, intervalle)
    plt.show()
        
#Galton(20, 1000)


def normale(k, sigma ):
    if k % 2 == 0 : 
        raise ValueError('le nombre k doit etre impair')
    x=np.linspace(-2*sigma, 2*sigma, k)
    tab=np.zeros(len(x))
    
    for i in range (len(x)):
        tab[i]= (1/m.sqrt(2*m.pi*sigma))*m.exp((-1/2)* pow((x[i]/sigma), 2))
    
    #plt.plot(x,tab)
    #plt.show()    
    return tab
    
#normale(51,12)


def proba_affine(k, slope) :
    if k % 2 == 0:
        raise ValueError ( 'le nombre k doit etre impair' )
    if abs ( slope  ) > 2.0 / ( k * k ):
        raise ValueError ( 'la pente est trop raide : pente max = ' + 
        str ( 2.0 / ( k * k ) ) )
    tab=np.zeros(k)
    for i in range (k) :
        tab[i]= (1.0/k) + (i - (k-1)/2)*slope
    #plt.plot(tab)
    #plt.show()
    return tab
    
#proba_affine(7,0.03 )

def Pxy(A,B) :
    tab=np.zeros((len(A), len(B)), float)
    for i in range(len(A)) :
        for j in range(len(B)) :
            tab[i,j]=A[i]*B[j]
    return tab
        
PxyA = np.array ( [0.2, 0.7, 0.1] )
PxyB = np.array ( [0.4, 0.24, 0.2] )
#print Pxy(PA,PB)        


def dessine ( P_jointe ):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x = np.linspace ( -3, 3, P_jointe.shape[0] )
    y = np.linspace ( -3, 3, P_jointe.shape[1] )
    X, Y = np.meshgrid(x, y)
    ax.plot_surface(X, Y, P_jointe, rstride=1, cstride=1 )
    ax.set_xlabel('A')
    ax.set_ylabel('B')
    ax.set_zlabel('P(A) * P(B)')
    plt.show ()
    

P_jointe = Pxy( normale(51,2), proba_affine(51,0.0003) )
#dessine(p)


 
def project1Var ( P, index ):
    """
    Supprime une variable de proba jointe
    Param P : une distribution de proba jointe sous forme d'un array
    a 1 dimension (toutes les variables aleatoires sont supposees binaires)
    Param index : represente la variable aleatoire a marginaliser
    """
    length = 2**( index + 1 )
    reste = 2**index
    vect = np.zeros ( P.size / 2 )
    for i in range ( P.size ):
        j = m.floor ( i / length ) * length / 2 + ( i % reste )
        vect[j] += P[i]
    return vect

proj1_P = np.array([0.05, 0.1, 0.15, 0.2, 0.02, 0.18, 0.13, 0.17])
#print project1Var( proj1_P, 1)   

    
def project ( P, ind_to_remove ):
    """
    Calcul de la projection d'une distribution de probas

    Param P une distribution de proba sous forme d'un array à 1 dimension
    Param ind_to_remove un array d'index representant les variables à
    supprimer. 0 = 1ère variable, 1 = 2ème var, etc.
    """
    v = P
    ind_to_remove.sort ()
    for i in range ( ind_to_remove.size - 1, -1, -1 ):
        v = project1Var ( v, ind_to_remove[i] )
    return v

#print project( proj1_P, np.array([1,2]))  

def expanse1Var ( P, index ):
    """
    duplique une distribution de proba |X| fois, où X est une des variables
    aléatoires de la probabilité jointe P. Les variables étant supposées
    binaires, |X| = 2. La duplication se fait à l'index de la variable passé
    en argument.
    Par exemple, si P = [0,1,2,3] et index = 0, expanse1Var renverra
    [0,0,1,1,2,2,3,3]. Si index = 1, expanse1Var renverra [0,1,0,1,2,3,2,3].

    Param P : une distribution de proba sous forme d'un array à 1 dimension
    Param index : représente la variable à dupliquer (0 = 1ère variable,
       1 = 2ème variable, etc).
    """
    length = 2**(index+1)
    reste = 2**index
    vect = np.zeros ( P.size * 2 )
    for i in range ( vect.size ):
        j = m.floor ( i / length ) * length / 2 + ( i % reste )
        vect[i] = P[j]
    return vect

#print expanse1Var (np.array([0.2, 0.3, 0.15, 0.35]), 1)    
   

def expanse ( P, ind_to_add ):
    """
    Expansion d'une probabilité projetée

    Param P une distribution de proba sous forme d'un array à 1 dimension
    Param ind_to_add un array d'index representant les variables permettant
    de dupliquer la proba P. 0 = 1ère variable, 1 = 2ème var, etc.
    """
    v = P
    ind_to_add.sort ()
    for ind in ind_to_add.size:
        v = expanse1Var ( v, ind )
    return v
    
def nb_vars ( P ):
    i = P.size
    nb = 0
    while i > 1:
        i /= 2
        nb += 1
    return nb
    
    
    
# def proba_conditionnelle(P) :
#     nb_vars 
    
    
    
    
#p=Pxy(normale(51,2),proba_affine(51,0.0003 ))
