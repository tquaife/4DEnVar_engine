import numpy as np
import matplotlib.pyplot as plt
import subprocess
import sys

class simpleEcosystem:

    def __init__(self,s1=100,p1=0.5,p2=0.001):
        self.s1_ini=s1
        self.s1=s1
        self.p1=p1
        self.p2=p2
        self.ts=0
                
    def integrate_ts(self,f):
        self.s1+=f*self.p1-self.s1*self.p2
        self.ts+=1
        
    def forcing(self):
        f=np.sin(self.ts/365.*np.pi*2)+1.0
        return f

    def run_model(self, nts):
        out_s1=[]
        out_f=[]
        self.s1=self.s1_ini
        self.ts=0
        for n in range(nts):
            f=self.forcing()
            self.integrate_ts(f)
            out_s1.append(self.s1)
            out_f.append(f)
        return (np.asarray(out_f),np.asarray(out_s1))


class simpleEcosystem_newF(simpleEcosystem):
    """change forcing for the model
    """
    def forcing(self):
        phase_shift=0.5
        f=np.sin((self.ts/365+phase_shift)*np.pi*2)+1.0
        if f<1.0:
            f=0
        else:
            f*=3.
        return f


class simpleEcosystemEnsemble:

    def __init__(self, n_ts ):

        self.n_ts=n_ts
        self.read_4DEnVar_inputs( )
        
        self.coefs_truth=(100,0.4,0.001)
        self.coefs_prior=(250,0.5,0.001)
        
        #not used, but for prosperity:
        self.uncert_prior=(50,0.05,0.0002)

        self.truth=simpleEcosystem_newF(s1=self.coefs_truth[0],p1=self.coefs_truth[1],p2=self.coefs_truth[2]) 
        self.prior=simpleEcosystem(s1=self.coefs_prior[0],p1=self.coefs_prior[1],p2=self.coefs_prior[2])


    def read_4DEnVar_inputs(self):
        
        #prior ensemble
        self.xb=np.genfromtxt("0xb.dat")       
        self.nens=np.shape(self.xb)[1]
        self.ensemble=[]
        for n in range(self.nens):
            self.ensemble.append(simpleEcosystem(s1=self.xb[0,n],p1=self.xb[1,n],p2=self.xb[2,n]))                

        self.hxb_bar=np.genfromtxt("0hxbar.dat")
        self.obs_x=np.genfromtxt("0obs_x.dat")
        self.obs_y=np.genfromtxt("0y.dat")


    def plot(self,analysis=None):

        x=np.arange(self.n_ts)                

        for n in range(self.nens):
            (_,state)=self.ensemble[n].run_model(self.n_ts)        
            if n==0:
                plt.plot(x,state,"y-",alpha=0.2,label="HX'b")
            else:    
                plt.plot(x,state,"y-",alpha=0.2)

        (_,prior_state)=self.prior.run_model(self.n_ts)        
        plt.plot(x,prior_state,"y-",label="h(xb)",linewidth=2.0)
 
        (_,true_state)=self.truth.run_model(self.n_ts)        
        plt.plot(x,true_state,"g-",label="truth")
        plt.plot(self.obs_x,self.obs_y,"go",label="y",linewidth=2.0)


        if analysis is not None:
            posterior=simpleEcosystem(s1=analysis[0],p1=analysis[1],p2=analysis[2])
            (_,posterior_state)=posterior.run_model(self.n_ts)        
            plt.plot(x,posterior_state,"r-",label="post",linewidth=2)            

        plt.xlabel("time step")
        plt.ylabel("model state")
        plt.legend()

        plt.plot()
        plt.show()

def do_4denvar():
    #run the 4DEnVar via a subprocess
    out=subprocess.run(["../../4DEnVar","0xb.dat","0HXb.dat","0y.dat","0R.dat","0hxbar.dat"],capture_output=True)
    out=out.stdout.decode("utf-8").rstrip().split("\n")

    #read the results of the analysis
    analysis=[]    
    for i in range(3):
        analysis.append(float(out[i]))

    return analysis



if __name__=="__main__":
  
    nyears=2
    n_ts=int(365*nyears)
    s=simpleEcosystemEnsemble(n_ts)
    #s.plot()

    analysis=do_4denvar()
    print(s.coefs_prior)
    print(s.coefs_truth)
    print(analysis)
    s.plot(analysis=analysis)


