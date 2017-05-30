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

def meta(key, betas, ses, labels, chains = 4, iter = 5000, warmup = 1000, cores = 2):
    print(key,betas,ses,labels)
    nstudies = len(labels)
    labels = numpy.array(labels)
    betas = numpy.array(betas)
    ses = numpy.array(ses)
    print(betas,ses)
    fn = 'static/images/' + str(key) + '.svg'
    ro.r('''
  metar<-function(nstudies, fn, betas, ses, labels, chains, iter, warmup, cores){       
library(svglite)
       library(rstan)
       require(GGally)
       require(network)
       require(sna)
       require(ggplot2)
       require(corpcor)
       require(RColorBrewer)

 meta_data <- list(
  J = nstudies,
  betahat = betas,
  sigma = ses
)

print(meta_data)
print(nstudies)
print(betas)
print(ses)
library(rstan)
fit1 <- stan(
  file = "stan/meta.stan",  # Stan program
  data = meta_data,    # named list of data
  chains = chains,             # number of Markov chains
  warmup = warmup,          # number of warmup iterations per chain
  iter = iter,            # total number of iterations per chain
  cores = cores,              # number of cores (using 2 just for the vignette)
  refresh = 1000          # show progress every refresh iterations
  ) 

my_labels = c("tau","mu", rev(labels))
print(my_labels)
print(names(fit1))
print(fit1, pars=c("beta", "mu", "tau", "lp__"), probs=c(.1,.5,.9))
p.1 <- plot(fit1, pars=c("beta", "mu", "tau"))  
p.1 <- p.1  + scale_y_continuous(labels = my_labels, breaks = seq(1,length(my_labels))) + theme(axis.text.y = element_text(face="bold", 
                           size=8))
ggsave(fn, plot = p.1, width = 6.5, height = 8, device = "svg") 
p.2 <- traceplot(fit1, pars = c("mu", "tau"), inc_warmup = TRUE, nrow = 2)
ggsave(paste(c(fn,2,".svg"), sep = "", collapse=""), plot = p.2, width = 6.5, height = 8, device = "svg")
}

       ''')
    metaf = ro.globalenv['metar']
    metaf = ro.r['metar']
    res = metaf(nstudies, fn, betas, ses, labels, chains, iter, warmup, cores)

