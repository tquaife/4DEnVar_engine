import numpy as np
from matplotlib import pyplot as plt
from  scipy.stats import gaussian_kde as kde
import scipy.stats as stats


class zeroOrderModelEnsemble:

    def __init__(self, truth, prior, uncert_prior, nens, nobs, obs_uncert ):

        self.truth=truth 
        self.prior=prior 
        self.uncert_prior=uncert_prior 
        self.obs_uncert=obs_uncert
        
        self.perturb_obs=True
        
        self.nobs=nobs
        self.nens=nens
              
        self.gen_prior_ensemble()
        self.gen_obs()
    
    def gen_prior_ensemble(self):
        self.ensemble=np.zeros([self.nens,len(self.truth)])
        for n in range(self.nens):
            self.ensemble[n,:]=self.prior+np.sqrt(self.uncert_prior)*np.random.randn(len(self.prior))
        
    def gen_obs(self):
        
        #note - self.nobs is an array containing
        # the number of observation per dimension
        self.obs=[]
        for n in range(len(self.nobs)):
            tmp=np.ones(self.nobs[n])*self.truth[n]
            if self.perturb_obs:
                self.obs.append(tmp+np.sqrt(self.obs_uncert[n])*np.random.randn(len(tmp)))
            else:
                self.obs.append(tmp)
        
    def plot(self,filename=None,analysis=None,full_posterior=None):

        for n in range(len(self.truth)):
            k_prior=kde(self.ensemble[:,n],bw_method='silverman')
    
            #work out plot positions for KDE
            mn=np.min(self.ensemble[:,n])
            mx=np.max(self.ensemble[:,n])            
            rng=mx-mn
            d=rng*0.7            
            positions=np.linspace(mn-d,mx+d,200)
   
            plt.fill_between(positions,k_prior(positions),label="prior",alpha=0.5)

            #add analytical distribution:
            plt.plot(positions, stats.norm.pdf(positions, prior[n], np.sqrt(prior_uncert[n])))
            
            if full_posterior is not None:
                k_post=kde(full_posterior[:,n],bw_method='silverman')
                plt.fill_between(positions,k_post(positions),label="posterior",alpha=0.5)

                #calculate analytical posterior
                p=self.uncert_prior[n]
                r=self.obs_uncert[n]
                x=self.prior[n]
                for m in range(self.nobs[n]):
                    y=self.obs[n][m]
                    k=p/(p+r)
                    x=x+k*(y-x)
                    p=(1-k)*p
                #plot analytical posterior
                plt.plot(positions, stats.norm.pdf(positions, x, np.sqrt(p)))


            #add prior mean and truth
            ymn=plt.gca().get_ylim()[0]
            ymx=plt.gca().get_ylim()[1]
            plt.plot([self.prior[n],self.prior[n]],[0,ymx],"k--",label="prior")
            plt.plot([self.truth[n],self.truth[n]],[0,ymx],"k-.",label="truth")

            if analysis is not None:
                plt.plot([analysis[n],analysis[n]],[0,ymx],"r:",label="analysis")


            #add obs
            plt.plot(self.obs[n],np.zeros(len(self.obs[n])),"ko",label="obs")

            plt.gca().set_ylim([0,ymx])            
            plt.legend()
            plt.show()
            plt.clf()
                             
    def write_files(self):
    
        #observations:
        with open("0y.dat","w") as f:
            for n in range(len(self.nobs)):
                for m in range(self.nobs[n]):
                    f.write("%f\n"%self.obs[n][m])

        #R matrix:
        #(got to be a better way of doing this!
        #problem is that there can be different
        #numbers of obs per dimension)
        with open("0R.dat","w") as f:
            i=0
            j=0
            for n in range(len(self.nobs)):
                for m in range(self.nobs[n]):
                    for nn in range(len(self.nobs)):
                        for mm in range(self.nobs[nn]):
                            if i==j:
                                f.write("%f "%self.obs_uncert[nn])
                            else:
                                f.write("0.0 ")
                            j+=1
                    f.write("\n")
                    j=0
                    i+=1


        #background ensemble:
        with open("0xb.dat","w") as f:
            for m in range(len(self.prior)):
                for n in range(self.nens):
                    f.write("%f "%self.ensemble[n,m])
                f.write("\n")


        #predicted observations from ensemble:
        with open("0hx.dat","w") as f:
            for n in range(len(self.nobs)):
                for nn in range(self.nobs[n]):
                    for m in range(self.nens):                        
                        f.write("%f "%self.ensemble[m,n])
                    f.write("\n")
             
        #predicted observations from expected value of 
        #the prior distribution:
        with open("0hxbar.dat","w") as f:
            for n in range(len(self.nobs)):
                for nn in range(self.nobs[n]):
                    f.write("%f\n"%np.mean(self.ensemble[:,n]))
                    #or:
                    #f.write("%f\n"%self.prior[n])
             

if __name__=="__main__":

    import subprocess

    #all uncertainties are in variances
    truth=[5.0,5.0]
    prior=[1.0,1.0]
    prior_uncert=[4.0,4.0]
    nens=200
    nobs=[1,4]
    obs_uncert=[4.0,4.0]

    z=zeroOrderModelEnsemble(truth,prior,prior_uncert,nens,nobs,obs_uncert)
    z.perturb_obs=False
    z.gen_obs()
    z.write_files()

    #run the 4DEnVar via a subprocess
    out=subprocess.run(["../4DEnVar","0xb.dat","0hx.dat","0y.dat","0R.dat","0hxbar.dat"],capture_output=True)
    out=out.stdout.decode("utf-8").rstrip().split("\n")

    #parse output
    analysis=[]
    for i in range(len(truth)):
        analysis.append(float(out[i]))

    posterior=np.zeros([nens,len(truth)])
    for i in range(len(truth)):
        line=out[i+len(truth)+1].split()
        for j in range(nens):
            posterior[j,i]=float(line[j])

        print(np.std(z.ensemble[:,i]),np.std(posterior[:,i]))    

    z.plot(analysis=analysis,full_posterior=posterior)



