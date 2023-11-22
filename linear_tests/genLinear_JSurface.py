from copy import copy
import numpy as np
from matplotlib import pyplot as plt
from genLinear import *

import sys

class linearModelEnsemble_JSurf(linearModelEnsemble):

    def __init__(self, coefs_truth, coefs_prior, uncert_prior, nens, nobs, obs_uncert, rand_obs_y=True ):
        super().__init__( coefs_truth, coefs_prior, uncert_prior, nens, nobs, obs_uncert, rand_obs_y )
        self.R_inv=self.get_R_inv()
        self.B_inv=self.get_B_inv()

    def surface_looper(self, range1, range2):
        """an iterator to loop across the error surface"""
        start,stop,step=range1
        range1=np.arange(start,stop,step)
        start,stop,step=range2
        range2=np.arange(start,stop,step)
        for (n,i) in enumerate(range1):
            for (m,j) in enumerate(range2):                
                yield (n,m,i,j)       
    
    def get_blank_surface(self,range1, range2):
        start,stop,step=range1
        range1=np.arange(start,stop,step)
        start,stop,step=range2
        range2=np.arange(start,stop,step)
        surf=np.zeros((len(range1),len(range2)))
        return surf
    
    def get_R_inv(self):
        R_inv=np.eye(len(self.obs_y))*1./np.power(self.obs_uncert,2)
        return R_inv
    
    def get_B_inv(self):
        B_inv=np.eye(len(self.uncert_prior))
        for i in range(len(self.uncert_prior)):
            B_inv[i,i]=1./np.power(self.uncert_prior[i],2)
        return B_inv

    def get_J_4DVar(self, params):
        """get the value of the 4DVar cost function
        at position given by params"""
        #background term:
        t1=np.matmul(self.B_inv, np.asarray(params)-np.asarray(self.coefs_prior))
        t1=np.matmul(np.asarray(params)-np.asarray(self.coefs_prior), t1)
        #observation term:
        hx=linearModel(params)
        t2=np.matmul(self.R_inv, np.asarray(hx.eval(self.obs_x))-np.asarray(self.obs_y))
        t2=np.matmul(np.asarray(hx.eval(self.obs_x))-np.asarray(self.obs_y), t2)
        
        #return(t1/2.)
        return((t1+t2)/2.)
    
    def get_JSurface_4DVar(self, dims=(0,1),range1=(0,2,0.1),range2=(0,2,0.1)):
        surf=self.get_blank_surface(range1,range2)
        for (n,m,i,j) in self.surface_looper(range1,range2):
            params=copy(self.truth.coefs)
            params[dims[0]]=i
            params[dims[1]]=j
            surf[n,m]=self.get_J_4DVar(params)
        return surf

    def get_JSurface_4DEnVar(self, data, dims=(0,1),range1=(0,2,0.1),range2=(0,2,0.1)):
        surf=self.get_blank_surface(range1,range2)
        for (n,m,i,j) in self.surface_looper(range1,range2):
            params=copy(self.truth.coefs)
            params[dims[0]]=i
            params[dims[1]]=j
            surf[n,m]=float(data.pop(0).split()[-1])                   
        return surf

    def write_JSurface_xeval_file(self, dims=(0,1),range1=(0,2,0.1),range2=(0,2,0.1)):
        with open("0x_eval.dat","w") as f:
            for (n,m,i,j) in self.surface_looper(range1,range2):
                params=copy(self.truth.coefs)
                params[dims[0]]=i
                params[dims[1]]=j
                for k in range(len(params)):
                    f.write("%f "%params[k])
                f.write("\n")

def plot_JSurface_compare(linear, data, analysis, dims=(0,1),range1=(0,2,0.1),range2=(0,2,0.1), filename=None):
    
    surf_4DEnVar=linear.get_JSurface_4DEnVar(data, dims,range1,range2)
    surf_4DVar=linear.get_JSurface_4DVar(dims,range1,range2)

    jmin = np.min([np.min(surf_4DEnVar), np.min(surf_4DVar)])
    jmax = np.max([np.max(surf_4DEnVar), np.max(surf_4DVar)])    

    fig, ax = plt.subplots(2, 2)

    #contour plots    
    cmap="inferno_r"
    #ax[0,0].imshow(surf_4DVar, cmap=cmap, vmin=jmin, vmax=jmax, extent=[range1[0],range1[1],range2[0],range2[1]])
    cs=ax[0,0].contourf(np.log(surf_4DVar), cmap=cmap, extent=[range1[0],range1[1],range2[0],range2[1]])
    ax[0,0].axis('square')
    ax[0,0].set_title("4DVar") 
    #ax[0,1].imshow(surf_4DEnVar, cmap=cmap, vmin=jmin, vmax=jmax, extent=[range1[0],range1[1],range2[0],range2[1]])
    ax[0,1].contourf(np.log(surf_4DEnVar), cmap=cmap, levels=cs.levels, extent=[range1[0],range1[1],range2[0],range2[1]])
    ax[0,1].axis('square')
    ax[0,1].set_title("4DEnVar") 

    #scatter plot
    ax[1,0].plot(np.ndarray.flatten(surf_4DEnVar),np.ndarray.flatten(surf_4DVar),"k.", markersize=0.1)
    ax[1,0].plot([jmin,jmax],[jmin,jmax],"k--", linewidth=0.1)
    ax[1,0].set_xlabel("J - 4DEnVar")
    ax[1,0].set_ylabel("J - 4DVar")
    ax[1,0].axis('square')
    
    #model fit summary:
    a=linear.range[0]
    b=linear.range[1]
    x=np.arange(a,b,(b-a)/100.)
    for n in range(linear.nens):
        ax[1,1].plot(x,linear.ensemble[n].eval(x),'k',alpha=0.2)
    ax[1,1].plot(x,linear.truth.eval(x),'r')
    ax[1,1].plot(linear.obs_x,linear.obs_y,'.r')
    l=linearModel(analysis)
    ax[1,1].plot(x,l.eval(x),'g')            
    ax[1,1].set_xlim(linear.range)
    ax[1,1].set_ylim((linear.truth.eval(a),linear.truth.eval(b)))
    ax[1,1].set_aspect(abs( (a-b)/(linear.truth.eval(a)-linear.truth.eval(b))))

    plt.tight_layout()

    if filename is None:
        plt.show()
    else:
        plt.savefig(filename)
    plt.clf()
    
if __name__=="__main__":
    import subprocess

    #ranges for error surface calculation
    range1=(0.,2,0.01)
    range2=(0.,2,0.01)

    truth=[0.75,1.25,2.]
    #coefs_truth, coefs_prior, uncert_prior, nens, nobs, obs_uncert
    l=linearModelEnsemble_JSurf(truth,[1.,1.,2.0],[0.5,0.5,0.5],100,10,0.01)
    
    l.write_files()
    l.write_JSurface_xeval_file(range1=range1,range2=range2)

    #run the 4DEnVar surface generator via a subprocess    
    data=subprocess.run(["../4DEnVar_surf","0xb.dat","0hx.dat","0y.dat","0R.dat","0hxbar.dat","0x_eval.dat"],capture_output=True)
    surface=data.stdout.decode("utf-8").rstrip().split("\n")
            
    #run the 4DEnVar optimiser via a subprocess
    out=subprocess.run(["../4DEnVar","0xb.dat","0hx.dat","0y.dat","0R.dat","0hxbar.dat"],capture_output=True)
    out=out.stdout.decode("utf-8").rstrip().split("\n")

    #read the results of the analysis
    analysis=[]    
    for i in range(len(truth)):
        analysis.append(float(out[i]))

    plot_JSurface_compare(l, surface, analysis, range1=range1,range2=range2, filename='linear_fit_summary.png')
            
    

    
    
