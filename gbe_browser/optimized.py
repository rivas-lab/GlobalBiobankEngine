
from __future__ import print_function
from __future__ import division
from random import shuffle

# Written by Manuel A. Rivas
# Updated 04.21.2017

import sys,re
import os
import numpy
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

def Polygenic(key, betas, ses, labels, chains = 8, iter = 200, warmup = 100, cores = 8):
    nstudies = len(labels)
    labels = numpy.array(labels)
    betas = numpy.array(betas)
    ses = numpy.array(ses)
    fn = 'static/images/Polygenic/Polygenic_' + str(key) + '.svg'
    ro.r('''
  covarr<-function(nstudies, fn, betas, ses, labels, chains, iter, warmup, cores){       


       library(rstan)
       require(ggplot2)
       require(corpcor)
       require(RColorBrewer)
       library(ggthemes)

       # Create data for Stan
      stan.data <- list(
           N = nrow(betas),
           M = ncol(betas),
           B = betas,
           SE = ses,
           K = 2
          )

 
sm <- stan_model(file = "model2_mixture.stan")

fit2 = optimizing(sm, data = stan.data, hessian = TRUE, as_vector=FALSE)
#fit2 = vb(sm, data = stan.data)
#print(fit2)     
  
     #  fit2 <- stan(
   #        file = "model2_mixture.stan",  # Stan program
   #        data = stan.data,    # named list of data
   #        chains = 8,             # number of Markov chains
   #        warmup = 100,          # number of warmup iterations per chain
   #        iter = iter,            # total number of iterations per chain
   #        cores = 8,              # number of cores (using 2 just for the vignette)
   #        refresh = 50          # show progress every 'refresh' iterations
   #       )

#print(fit2, pars=c("L_Omega", "L_Theta", "tau", "pi","Omegacor","Thetacor"), probs=c(0.025, 0.5, 0.975), digits_summary=5)
#print(fit2, pars=c("L_Omega", "L_Theta", "tau", "pi","Omegacor","Thetacor"), digits_summary=5)
#print(fit2, pars = c("L_omega","L_Theta","tau","pi","Omegacor","Thetacor"))
print(fit2$par$pi)
print(fit2$par$Omegacor)
print(fit2$par$Thetacor)
print(labels)
#mylabels = c("tau","pi", rev(labels),"Thetacor")
#mylabels = c(rep(1,nstudies),rep(1,2),rep(1,nstudies*nstudies),rep(1,nstudies*nstudies))

#d = cbind(rep("tau_",nstudies),rev(labels))
#colnames(d) = c("id","val")
#d = data.frame(d)
#mylabels[1:nstudies] = with(d,paste0(id,val))
#mylabels[nstudies+1] = "pi_m1"
#mylabels[nstudies+2] = "pi_null"

#labelsrev = rev(labels)

#init = nstudies + 2
#for(i in 1:nstudies){
#   for(j in 1:nstudies){
#      init = init + 1
#      mylabels[init] = paste(c("cor_e_", labelsrev[j], "_", labelsrev[i]), sep = "", collapse="")
#}

#}

#for(i in 1:nstudies){
#   for(j in 1:nstudies){
#      init = init + 1
#      mylabels[init] = paste(c("cor_g_", labelsrev[j], "_", labelsrev[i]), sep = "", collapse="")
#}

#}



#print(mylabels)
#p.1 <- plot(fit2, pars=c("Omegacor", "Thetacor", "pi", "tau"))  
#p.1 <- p.1  +  scale_y_continuous(labels = mylabels, breaks = seq(1,length(mylabels))) +  theme_wsj() + scale_colour_wsj("colors6", "")
#ggsave(fn, plot = p.1, width = 8, height = 8, device = "svg") 
#p.2 <- traceplot(fit2, pars = c( "tau","pi","Omegacor","Thetacor"), inc_warmup = TRUE, nrow = 5) +  theme_wsj() + scale_colour_wsj("colors6", "")
#ggsave(paste(c(fn,2,".svg"), sep = "", collapse=""), plot = p.2, width = 8, height = 8, device = "svg")

## extract the simulated draws from the posterior and note the number for nsims
#theta = extract(fit2)

#print(names(theta))
#nsims = length(theta$Sigmas)
#print(nsims)
#print(dim(theta$Sigmas))


#save.image(paste(fn, ".RData", sep=""))

}

       ''')
    covarf = ro.globalenv['covarr']
    covarf = ro.r['covarr']
    res = covarf(nstudies, fn, betas, ses, labels, chains, iter, warmup, cores)




## START NEW 

def PolygenicCoding(key, betas, ses, labels, chains = 8, iter = 200, warmup = 100, cores = 8):
    nstudies = len(labels)
    labels = numpy.array(labels)
    betas = numpy.array(betas)
    ses = numpy.array(ses)
    fn = 'static/images/PolygenicCoding/PolygenicCoding_' + str(key) + '.svg'
    ro.r('''
  covarr<-function(nstudies, fn, betas, ses, labels, chains, iter, warmup, cores){       

       library(rstan)
       require(ggplot2)
       require(RColorBrewer)
       library(ggthemes)
    require(svglite)
       # Create data for Stan
      stan.data <- list(
           N = nrow(betas),
           M = ncol(betas),
           B = betas,
           SE = ses,
           K = 2
          )

 
sm <- stan_model(file = "model2_mixture.stan")

#fit2 = optimizing(sm, data = stan.data, hessian = TRUE, as_vector=FALSE)
#fit2 = vb(sm, data = stan.data)
#print(fit2)     
  
      fit2 <- stan(
           file = "model2_mixture.stan",  # Stan program
           data = stan.data,    # named list of data
           chains = 4,             # number of Markov chains
           warmup = 100,          # number of warmup iterations per chain
           iter = iter,            # total number of iterations per chain
           cores = 8,              # number of cores (using 2 just for the vignette)
           refresh = 50          # show progress every 'refresh' iterations
          )



print(fit2, pars=c("L_Omega", "L_Theta", "tau", "pi","Omegacor","Thetacor"), probs=c(0.025, 0.5, 0.975), digits_summary=5)

#print(fit2$par$pi)
#print(fit2$par$Omegacor)
#print(fit2$par$Thetacor)
#print(labels)
mylabels = c("tau","pi", rev(labels),"Thetacor")
mylabels = c(rep(1,nstudies),rep(1,2),rep(1,nstudies*nstudies),rep(1,nstudies*nstudies))

d = cbind(rep("tau_",nstudies),rev(labels))
colnames(d) = c("id","val")
d = data.frame(d)
mylabels[1:nstudies] = with(d,paste0(id,val))
mylabels[nstudies+1] = "pi_m1"
mylabels[nstudies+2] = "pi_null"

labelsrev = rev(labels)

init = nstudies + 2
for(i in 1:nstudies){
   for(j in 1:nstudies){
      init = init + 1
      mylabels[init] = paste(c("cor_e_", labelsrev[j], "_", labelsrev[i]), sep = "", collapse="")
}

}

for(i in 1:nstudies){
   for(j in 1:nstudies){
      init = init + 1
      mylabels[init] = paste(c("cor_g_", labelsrev[j], "_", labelsrev[i]), sep = "", collapse="")
}

}



print(mylabels)
p.1 <- plot(fit2, pars=c("Omegacor", "Thetacor", "pi", "tau"))  
p.1 <- p.1  +  scale_y_continuous(labels = mylabels, breaks = seq(1,length(mylabels))) +  theme_wsj() + scale_colour_wsj("colors6","")
#p.1 <- p.1  +  scale_y_continuous(labels = mylabels, breaks = seq(1,length(mylabels))) +  theme_hc() + scale_colour_hc("colors6", "")
#    svg(file = fn, width = 8, height = 8)
#    p.1
#    dev.off()

p.2 <- traceplot(fit2, pars = c( "tau","pi","Omegacor","Thetacor"), inc_warmup = TRUE, nrow = 5) +  theme_wsj() + scale_colour_wsj("colors6", "")
#    svg(file = paste(c(fn,2,".svg"), sep = "", collapse=""), width = 8, height = 8)
#    p.2
#    dev.off()
ggsave(fn, plot = p.1, width = 8, height = 8, device = "svg") 
ggsave(paste(c(fn,2,".svg"), sep = "", collapse=""), plot = p.2, width = 8, height = 8, device = "svg")

## extract the simulated draws from the posterior and note the number for nsims
#theta = extract(fit2)

#print(names(theta))
#nsims = length(theta$Sigmas)
#print(nsims)
#print(dim(theta$Sigmas))


#save.image(paste(fn, ".RData", sep=""))

}

       ''')
    covarf = ro.globalenv['covarr']
    covarf = ro.r['covarr']
    res = covarf(nstudies, fn, betas, ses, labels, chains, iter, warmup, cores)

