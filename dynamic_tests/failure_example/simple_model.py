import numpy as np
import matplotlib.pyplot as plt
import subprocess


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

    def __init__(self, n_ts, coefs_truth, coefs_prior, uncert_prior, nens, nobs, obs_uncert ):

        self.n_ts=n_ts
        #self.truth=simpleEcosystem(s1=coefs_truth[0],p1=coefs_truth[1],p2=coefs_truth[2]) 
        self.truth=simpleEcosystem_newF(s1=coefs_truth[0],p1=coefs_truth[1],p2=coefs_truth[2]) 
        self.prior=simpleEcosystem(s1=coefs_prior[0],p1=coefs_prior[1],p2=coefs_prior[2])

        self.coefs_prior=coefs_prior 
        self.uncert_prior=uncert_prior 
        self.obs_uncert=obs_uncert
        
        self.nobs=nobs
        self.nens=nens

        self.gen_prior_ensemble()
        self.gen_obs()

    def gen_obs(self):    
        self.obs_x=np.random.rand(self.nobs)
        #sort not needed, but makes looking at
        #the outputs a bit more intuitive:
        self.obs_x.sort()
        self.obs_x=np.asarray(self.obs_x*self.n_ts,dtype=int)
        (_,state)=self.truth.run_model(self.n_ts)        
        self.obs_y=state[self.obs_x]+np.random.randn(self.nobs)*self.obs_uncert
        
    def gen_prior_ensemble(self):
        self.ensemble=[]
        for n in range(self.nens):
            coefs=np.random.randn(len(self.coefs_prior))
            coefs=coefs*self.uncert_prior
            coefs=self.coefs_prior+coefs
            self.ensemble.append(simpleEcosystem(s1=coefs[0],p1=coefs[1],p2=coefs[2]))                


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
            for n in range(self.nens):
                f.write("%f "%self.ensemble[n].s1)
            f.write("\n")
            for n in range(self.nens):
                f.write("%f "%self.ensemble[n].p1)
            f.write("\n")
            for n in range(self.nens):
                f.write("%f "%self.ensemble[n].p2)
            f.write("\n")

        #predicted observations from ensemble:
        with open("0HXb.dat","w") as f:
            ens=[]
            for n in range(self.nens):
                (_,state)=self.ensemble[n].run_model(self.n_ts)
                ens.append(state)   
            for mx in self.obs_x:
                for n in range(self.nens):
                    f.write("%f "%ens[n][mx])
                f.write("\n")
             
        #predicted observations from expected value of 
        #the prior distribution:
        with open("0hxbar.dat","w") as f:
            (_,state)=self.prior.run_model(self.n_ts)
            for mx in self.obs_x: 
                f.write("%f\n"%state[mx])
             


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
    for i in range(len(coefs_truth)):
        analysis.append(float(out[i]))

    return analysis


def simple_model_run():
    """simple demo function for testing
    """
    nts=365*10
    s=simpleEcosystem()
    (forcing,state)=s.run_model(nts)
    x=np.arange(nts)
    plt.plot(x,state)
    plt.xlabel("time step")
    plt.ylabel("model state")
    plt.show()

if __name__=="__main__":

    nyears=2
    n_ts=int(365*nyears)
    coefs_truth=(100,0.4,0.001)
    coefs_prior=(250,0.5,0.001)
    uncert_prior=(50,0.05,0.0002)
    #uncert_prior=(5,0.005,0.00002)
    nobs=10
    nens=20
    obs_uncert=1 #stddev
    s=simpleEcosystemEnsemble(n_ts,coefs_truth,coefs_prior,uncert_prior,nens,nobs,obs_uncert)
    #turned off writing to preserve the 
    #specific case previously generated
    #s.write_files()

    analysis=do_4denvar()
    print(coefs_prior)
    print(coefs_truth)
    print(analysis)
    s.plot(analysis=analysis)


