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

#print project1Var( np.array([0.05, 0.1, 0.15, 0.2, 0.02, 0.18, 0.13, 0.17]), 1)   

    
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

print project( np.array([0.05, 0.1, 0.15, 0.2, 0.02, 0.18, 0.13, 0.17]), np.array([1,2]))  

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
    for ind in ind_to_add:
        v = expanse1Var ( v, ind )
    return v
    
def nb_vars ( P ):
    i = P.size
    nb = 0
    while i > 1:
        i /= 2
        nb += 1
    return nb
    
    
    
def proba_conditionnelle ( P ) :
    """
    P(A|B) = P(A,B) / P(B)
    P(A|B,C) = P(A,B,C) / P(B,C)
    Note : La somme n'est pas égale à 1
    """
    n = nb_vars ( P ) - 1
    P_XnlXi = project1Var ( P, n )
    P_XnlXi_double = expanse1Var ( P_XnlXi, n )
    #res = P
    # L'EGALITE EST UN PUTAIN DE POITEUR EN PYTHON, QUAND TU MODIFIAIS RES, TU MODIFIAIS P
    # GG SALE MERDE, TU M'AS FAIT PERDRE 2H DE MA VIE
    res=np.zeros(len(P_XnlXi_double))
    for i in range ( len(res) ) :
        if P_XnlXi_double[i] != 0 :
            res[i] = P[i] / P_XnlXi_double[i]
        else :
            res[i] = float(0)
    return res

#print proba_conditionnelle ( np.array([0.05, 0.1, 0.15, 0.2, 0.02, 0.18, 0.13, 0.17]) )
    

def is_indep ( P, index, epsilon = m.exp(-6) ) :
    if nb_vars ( P ) - 1 == index :
        raise ValueError ( 'L\indice i d\'une variable Xi doit être différent de n')
    P_cond = proba_conditionnelle ( P )
    P_no_Xi = project1Var ( P, index )
    P_no_Xi_double = expanse1Var ( P_no_Xi, index )
    P_cond_no_Xi = proba_conditionnelle ( P_no_Xi_double )
    for ind in range ( len(P)) :
        if ( abs ( P_cond_no_Xi[ind] - P_cond[ind] ) > epsilon ) :
            return False
    return True

def test_is_indep ():
    P_asia = np.loadtxt ( 'asia_2014.txt' )
    for i in range ( 7 ) :    
        print str(i) + " : " + str(is_indep ( P_asia , i ))
       
#test_is_indep ()
#print is_indep ( np.array([0.25, 0.25, 0.25, 0.25]), 0)


def find_indep ( P, epsilon = m.exp(-6) ) :
    """
    E: P: np.array proba jointe taille n
    E: epsilon : écart
    S: n: nb variables dans la proba jointe passée en argument
    S: proba_cond: proba cond réduisant la taille à n/2
    S: indep_i: array composé des i des Xi qui constituent la proba cond (!isIndep)
    """
    n = nb_vars ( P ) - 1
    
    isIndep = np.zeros ( n )
    nbTrue = 0
    for i in range ( n ) : 
        isIndep[i] = is_indep ( P, i, epsilon )
        if isIndep[i] :
            nbTrue += 1
     
    indep_i = np.zeros(nbTrue)
    ind_indep = 0
    for i in range(n) :
        if isIndep[i] :
            indep_i[ind_indep] = i
            ind_indep += 1
    proba_cond = project(P, indep_i)
    
    return n, proba_cond, indep_i.astype(int)
    
def test_find_indep ():
    P_asia = np.loadtxt ( 'asia_2014.txt' )   
    print find_indep ( P_asia )

#test_find_indep ()


def find_all_indep ( P, epsilon = m.exp(-6) ) :
 
    nb_v = nb_vars ( P ) #nb vars initial used
    n = nb_v - 1
 
    list_tab_i = [[]] #[[n],[n,n-1],...[n,...1]]
    for j in range ( n ) :
        tab_i = np.zeros ( j+1, dtype=int ) #because begin at j=0
        for i in range ( 1, j+2 ) : #because (n-1)+2=n+1 excluded
            tab_i[i-1] = nb_v - i
        list_tab_i.append(tab_i)
    mat_i = np.array(list_tab_i)
 
    string = 'mat_i: ' + str(mat_i) + '\n'
 
    list_cond = []
    conso_totale = 0 #nb vars (total) used for this new formula
    i = n
    for tab_i in mat_i :
        
        if i == n :
            cond_tmp = np.copy(proba_conditionnelle (P))
            nb_v_distrib = nb_v
            
        else :
            proj_tmp = project ( P, tab_i )
            nb_v_distrib = i+1
            cond_tmp = np.copy(proba_conditionnelle ( proj_tmp ))
            
        indep_nb, indep_tab, indep_var = find_indep ( cond_tmp, epsilon ) 
        nb_v_compact = nb_vars ( indep_tab ) 
        conso_totale += nb_v_compact
        
        cond = expanse ( indep_tab, tab_i )
        
        list_cond.append ( cond ) 

        string += 'i: ' + str ( i ) + '\n'
        string += ' :nb_v_distrib: ' + str ( nb_v_distrib ) + '\n'
        string += ' :nb_v_compact: ' + str ( nb_v_compact ) + '\n'
        string += ' :cond:  ' + str(cond) + '\n'     
        string += ' :tab_i: ' + str(tab_i) + '\n'
        i -= 1
        
    P_result = np.ones ( P.size )
    for cond in list_cond :
        P_result *= cond
 
    string += 'P_result: ' + str(P_result) + '\n'
    string += 'P       : ' + str(P) + '\n'
    string += 'sum P_res: ' + str(sum(P_result)) + '\n'
    string += 'sum P    : ' + str(sum(P)) + '\n'
    string += 'conso P: ' + str(nb_v) + '\n'
    string += 'conso compact: ' + str(conso_totale)

    print string


find_all_indep ( np.array([0.05, 0.1, 0.15, 0.2, 0.02, 0.18, 0.13, 0.17]) )
    
