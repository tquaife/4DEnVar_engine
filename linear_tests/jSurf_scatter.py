from genLinear_JSurface import *
import numpy as np
from matplotlib import pyplot as plt

if __name__=="__main__":
    import subprocess

    truth=[1.,1.,1.]

    #ranges for error surface calculation
    dims=[0,1]
    range1=(0.,2,0.02)
    range2=(0.,2,0.02)

    #coefs_truth, coefs_prior, uncert_prior, nens, nobs, obs_uncert
    l=linearModelEnsemble_JSurf(truth,truth,[0.05,0.05,0.05],100,10,0.01,rand_obs_y=True)
    
    l.write_files()
    l.write_JSurface_xeval_file(range1=range1,range2=range2)
    
    data=subprocess.run(["../4DEnVar_surf","0xb.dat","0hx.dat","0y.dat","0R.dat","0hxbar.dat","0x_eval.dat"],capture_output=True)
    data=data.stdout.decode("utf-8").rstrip().split("\n")
        

    surf_4DEnVar=l.get_JSurface_4DEnVar(data, dims,range1,range2)
    surf_4DVar=l.get_JSurface_4DVar(dims,range1,range2)

    jmin = np.min([np.min(surf_4DEnVar), np.min(surf_4DVar)])
    jmax = np.max([np.max(surf_4DEnVar), np.max(surf_4DVar)])    

    print(jmin,jmax)

    plt.plot(np.ndarray.flatten(surf_4DEnVar),np.ndarray.flatten(surf_4DVar),"y.")
    plt.plot([jmin,jmax],[jmin,jmax],"k--")
    plt.xlabel("J - 4DEnVar")
    plt.ylabel("J - 4DVar")
    plt.show()

