#!/opt/anaconda/bin/python

import numpy as np
#import scipy.stats as st
import matplotlib.pyplot as plt
#from sklearn.neighbors import KernelDensity
from  scipy.stats import gaussian_kde as kde

Xa=np.genfromtxt("Xa.dat")
Xb=np.genfromtxt("Xb.dat")

for P in [0,1,2,3,4]:

    ka=kde(Xa[P,:])
    kb=kde(Xb[P,:])

    mn=min(np.min(Xa[P,:]),np.min(Xb[P,:]))
    mx=max(np.max(Xa[P,:]),np.max(Xb[P,:]))

    rng=mx-mn
    d=rng*0.3
    
    positions=np.linspace(mn-d,mx+d,200)

    #plt.plot(positions,ka(positions),label="posterior")
    plt.fill_between(positions,ka(positions),label="posterior",alpha=0.5)
    #plt.plot(positions,kb(positions),label="prior")
    plt.fill_between(positions,kb(positions),label="prior",alpha=0.5)
    plt.plot(Xa[P,:],np.zeros(len(Xa[P,:])),"ko",label="post. samples")
    plt.title("KDE for parameter "+str(P+1))
    plt.xlabel("parameter value")
    plt.ylabel("likliehood")
    plt.legend()
    plt.savefig("kde_"+str(P+1)+".png")
    plt.clf()
