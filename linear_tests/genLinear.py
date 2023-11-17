import numpy as np
from matplotlib import pyplot as plt


class linearModel:

    def __init__(self, coefs):
        self.coefs=coefs

    def eval(self, x):
        return np.polyval(self.coefs,x)

class linearModelEnsemble:

    def __init__(self, coefs_truth, coefs_prior, uncert_prior, nens, nobs, obs_uncert, rand_obs_y ):

        self.truth=linearModel(coefs_truth) 
        self.prior=linearModel(coefs_prior) 

        self.coefs_prior=coefs_prior 
        self.uncert_prior=uncert_prior 
        self.obs_uncert=obs_uncert
        
        self.nobs=nobs
        self.nens=nens
              
        self.range=(0,1)
        self.rand_obs_x=False
        self.rand_obs_y=rand_obs_y

        self.gen_prior_ensemble()
        self.gen_obs()
    
    def gen_prior_ensemble(self):
        self.ensemble=[]
        for n in range(self.nens):
            coefs=np.random.randn(len(self.coefs_prior))
            coefs=coefs*self.uncert_prior
            coefs=self.coefs_prior+coefs
            self.ensemble.append(linearModel(coefs))                
    
    def gen_obs(self):
        if self.rand_obs_x:
            self.obs_x=np.random.rand(self.nobs)
            self.obs_x=self.obs_x*(self.range[1]-self.range[0])+self.range[0]
        else:
            a=self.range[0]
            b=self.range[1]
            step=(b-a)/self.nobs
            self.obs_x=np.arange(a+step/2.,b,step)
        if self.rand_obs_y:
            self.obs_y=self.truth.eval(self.obs_x)+np.random.randn(self.nobs)*self.obs_uncert
        else:
            self.obs_y=self.truth.eval(self.obs_x)
            
    def plot(self,filename=None,analysis=None):
        a=self.range[0]
        b=self.range[1]
        x=np.arange(a,b,(b-a)/100.)
    
        for n in range(self.nens):
            plt.plot(x,self.ensemble[n].eval(x),'k',alpha=0.2)
    
        plt.plot(x,self.truth.eval(x),'r')
        plt.plot(self.obs_x,self.obs_y,'.r')
    
        if analysis is not None:
            l=linearModel(analysis)
            plt.plot(x,l.eval(x),'g')            
    
        plt.xlim(self.range)
        plt.ylim((self.truth.eval(a),self.truth.eval(b)))
        if filename is None:
            plt.show()
        else:    
            plt.savefig(filename)
        plt.clf()
               
    def write_files(self):
    
        #observations:
        with open("0y.dat","w") as f:
            for n in range(self.nobs):
                f.write("%f\n"%self.obs_y[n])

        #R matrix:
        with open("0R.dat","w") as f:
            for m in range(self.nobs):
                for n in range(self.nobs):
                    if n==m:
                        f.write("%f "%np.power(self.obs_uncert,2))
                    else:
                        f.write("0.0 ")
                f.write("\n")

        #parameter ensemble:
        with open("0xb.dat","w") as f:
            for m in range(len(self.coefs_prior)):
                for n in range(self.nens):
                    f.write("%f "%self.ensemble[n].coefs[m])
                f.write("\n")

        #predicted observations from ensemble:
        with open("0hx.dat","w") as f:
            for m in range(self.nobs):
                for n in range(self.nens):
                    f.write("%f "%self.ensemble[n].eval(self.obs_x[m]))
                f.write("\n")
             
        #predicted observations from expected value of 
        #the prior distribution:
        with open("0hxbar.dat","w") as f:
            for m in range(self.nobs):
                f.write("%f\n"%self.prior.eval(self.obs_x[m]))
             


if __name__=="__main__":

    import subprocess

    #coefs_truth, coefs_prior, uncert_prior, nens, nobs, obs_uncert
    truth=[2.,1.1,0.]
    l=linearModelEnsemble(truth,[1.,0.5,0.3],[0.2,0.2,0.2],20,10,0.01)
    l.write_files()
        
    #run the 4DEnVar via a subprocess
    out=subprocess.run(["../4DEnVar","0xb.dat","0hx.dat","0y.dat","0R.dat","0hxbar.dat"],capture_output=True)
    out=out.stdout.decode("utf-8").rstrip().split("\n")

    #read the results of the analysis
    analysis=[]    
    for i in range(len(truth)):
        analysis.append(float(out[i]))
            
    print(analysis)
    #l.plot(filename="linear_example.png",analysis=analysis)
    l.plot(analysis=analysis)


