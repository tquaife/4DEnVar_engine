from copy import copy
import numpy as np
from matplotlib import pyplot as plt
from genLinear import *

class linearModelEnsemble(linearModelEnsemble):
    
    def R_inv(self):
        R_inv=np.eye(len(self.obs_y))*1./self.obs_uncert
        return R_inv
    
    def B_inv(self):
        B_inv=np.eye(len(self.uncert_prior))
        for i in range(len(self.uncert_prior)):
            B_inv[i,i]=1./self.uncert_prior[i]
        return B_inv


    def get_J_4DVar(self, params):
        """get the value of the 4DVar cost function
        at position x
        """
        R_inv=self.R_inv()
        B_inv=self.B_inv()
        
        #background term
        t1=np.matmul(B_inv, np.asarray(params)-np.asarray(self.coefs_prior))
        t1=np.matmul(np.asarray(params)-np.asarray(self.coefs_prior), t1)

        #observation term
        hx=linearModel(params)
        t2=np.matmul(R_inv, np.asarray(hx.eval(self.obs_x))-np.asarray(self.obs_y))
        t2=np.matmul(np.asarray(hx.eval(self.obs_x))-np.asarray(self.obs_y), t2)

        return((t1+t2)/2.)
    
    def plot_JSurface(self, dims=(0,1),range1=(0,2,0.1),range2=(0,2,0.1)):
        
        start,stop,step=range1
        range1=np.arange(start,stop,step)
        start,stop,step=range2
        range2=np.arange(start,stop,step)
        
        surf=np.zeros((len(range1),len(range2)))
        
        #print(surf)
        
        for (n,i) in enumerate(range1):
            for (m,j) in enumerate(range2):
                params=copy(self.truth.coefs)
                params[dims[0]]=i
                params[dims[1]]=j
                surf[n,m]=self.get_J_4DVar(params)
                #print(i,j,self.get_J_4DVar(params))
        
        plt.imshow(surf)
        plt.show()
    
if __name__=="__main__":
    import subprocess

    truth=[1.,1.,0.]
    #coefs_truth, coefs_prior, uncert_prior, nens, nobs, obs_uncert
    l=linearModelEnsemble(truth,[1.,0.5,0.3],[0.02,0.02,0.02],20,5,0.01)

    print(l.get_J_4DVar([0,0,0]))
    print(l.R_inv())
    print(l.B_inv())
    
    l.plot_JSurface(range1=(0,2,0.01),range2=(0,2,0.01))
    
    
