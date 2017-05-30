from __future__ import print_function
from __future__ import division
from random import shuffle

# Written by Manuel A. Rivas
# Updated 11.28.2016

from optparse import OptionParser
from collections import Counter
import array
import itertools
import math
import sys,re
import os
import logging
from scipy.stats import binom as binomial
from scipy.stats import wishart
import numpy as np
import numpy.matlib
import time
from scipy.stats import invgamma
import sklearn
import sklearn.covariance
# Set up basic logging
logger = logging.getLogger('Log')
from scipy import stats
from scipy.stats import multivariate_normal
import random

def is_pos_def(x):
    i = 0
    x = np.matrix(x)
    if np.all(np.linalg.eigvals(x) > 0):
        return True
    else:
       return False

# return BIC -2*log(p(Data | theta that maximizes C, Mc)) + vc log(n) : vc is the number of parameters (K+J)*(C-1), K is the number of phenotypes, J is the number of genes, C is the number of clusters
def mr(betas,ses,vymat,annotvec,genevec,protvec,chroffvec,clusters,fout,Rphen,Rpheninv,phenidarr, Rphenuse=True,  niter=200,burn=100,thinning=1,verbose=True, outpath='/Users/mrivas/'):
    print("Running MCMC algorithm...")
    epsilon = .0000000000000001
    storephensvar = []
    S = vymat
    xi0 = 1 # hyperparameter to control spread of proposals for annotation
    xialpha0 = 1
    betas = numpy.matrix(betas)
    ses = numpy.matrix(ses)
    S = numpy.matrix(S)
    Sinv = numpy.linalg.inv(S)
    # Let k be the number of clusters, where cluster 1 is the null model cluster
    C = clusters
    maxloglkiter = np.zeros((niter+2,1))
    # Let k be the number of phenotypes
    k = betas.shape[1]
    # Let m be the number of variants
    m = betas.shape[0]
    # Initialize
    #Sigma0 for alternative clusters
    if Rphenuse:
        if is_pos_def(Rphen):
            Theta0 = Rphen
            Theta0inv = Rpheninv
            Omega0 = 0.2*0.2*Theta0
            Omega0inv = numpy.linalg.inv(Omega0)
        else:
            Theta0  = sklearn.covariance.shrunk_covariance(Rphen)
            Theta0inv = numpy.linalg.inv(Theta0)
            Omega0 = 0.2*0.2*Theta0
            Omega0inv = numpy.linalg.inv(Omega0)
    else:
        Theta0  = numpy.eye(Rphen.shape[0])
        Theta0inv = numpy.linalg.inv(Theta0)
        Omega0 = 0.2*0.2*Theta0
        Omega0inv = numpy.linalg.inv(Omega0)
    #scale matrix
    geneset = set(genevec)
    genemap = list(geneset)
    annotset = set(annotvec)
    annotlen = len(annotset)
    annotmap = list(annotset)
    scales = numpy.zeros((niter+2,annotlen))
    # store the mean trait value across the clusters for individuals that are members
    bc = numpy.zeros((niter+2,C,k))
    # store the probabilities (proportions) of cluster memberships
    pc = numpy.zeros((niter+2,1,C))
    # for each iteration keep record of the variant membership
    deltam = numpy.zeros((niter+2,m))
    ###### Why are these stored separately?
    # non-normalized probabilities for each individual variant
    uc = numpy.zeros((niter+2,m,C))
    # normalized probabilities for each individual variant
    ws = numpy.zeros((niter+2,m,C))
    # for each iteration keep record of the variant membership
    tm = numpy.zeros((niter+2,m))
    sigmainvdict = {}
    sigmadict = {}
    thetadict = {}
    thetainvdict = {}
    for clusteriter in range(2,C+1):
        sigmadict[0,clusteriter] = S
        sigmainvdict[0,clusteriter] = Sinv
        thetadict[0,clusteriter] = Omega0
        thetainvdict[0,clusteriter] = Omega0inv
    pc[0,0,:] = np.random.dirichlet([1]*C)
    bc[0,0,:] = np.array([0]*k)
    for clusteridx in range(0,C):
        bc[0,clusteridx,:] = bc[0,0,:]
    # initialize variant membership across clusters
    deltam[0,:] = np.random.randint(1,C+1,m)
    # Iterations MCMC samplers
    for iter in range(1,niter+1):
        if iter % 100 == 0:
            print(iter)
        ## a) Update \pi_0 : Proposal centered around the current value, 
        cnts = [1]*C
        for clusteridx in range(1,C+1):
            for varidx in range(0,m):
                if deltam[iter-1,varidx] == clusteridx:
                    cnts[clusteridx-1] += 1
        pcproposal = np.random.dirichlet(pc[iter-1,0,:] + cnts)
        pc[iter,0,:] = pcproposal
        # c) Update delta_jm
        xk = numpy.arange(1,C+1)
        for varidx in range(0,m):
            probmjc = [0]*C
            lprobmjcu = [0]*C
            uc = [0]*C
            varannot = annotvec[varidx]
            annotidx = [i for i in range(0,annotlen) if annotmap[i] == varannot][0]
            genevar = genevec[varidx]
            geneid = [i for i in range(0,len(genemap)) if genemap[i] == genevar][0]
            atmp = np.array(ses[varidx,:])[0]
            dtmp = numpy.matlib.eye(len(atmp))
            np.fill_diagonal(dtmp,atmp)
            Vjm = dtmp*S*dtmp + np.matlib.eye(S.shape[0])*.000001
            # Gives covariance matrix of variant  effect on sets of phenotypes (after fixed effect meta-analysis has been applied across all studies available)
            for cidx in range(1,C+1):
                if cidx == 1:
                    llk2 = multivariate_normal.logpdf(betas[varidx,:],bc[0,0,:],Vjm)  + np.log(pc[iter,0,cidx-1])
                else:
                    llk2 = multivariate_normal.logpdf(betas[varidx,:],bc[0,0,:],Vjm + thetadict[iter-1,cidx]) + np.log(pc[iter,0,cidx-1])
                if deltam[iter-1,varidx] == cidx:
                    maxloglkiter[iter-1,0] += llk2
                lprobmjcu[cidx-1] += llk2
            #normalize uc - set to wc
            maxloglk = numpy.max(lprobmjcu)
            for cidx in range(0,C):
                uc[cidx] = numpy.exp(lprobmjcu[cidx] - maxloglk)
            for cidx in range(0,C):
                probmjc[cidx] = uc[cidx]/numpy.sum(uc)
            if numpy.isnan(probmjc[0]):
                wstmp = numpy.random.dirichlet(numpy.repeat(numpy.array([1]),C,axis = 0))
                custm = stats.rv_discrete(name='custm',values=(xk,wstmp))
            else:
                custm = stats.rv_discrete(name='custm',values=(xk,probmjc))
            deltam[iter,varidx] = custm.rvs(size=1)[0]
       # print(deltam[iter,:])
        # d) Update Sigma_c using a Gibbs update from a Gaussian distribution
        for cidx in range(2,C+1):
            cnt = 0
            betastmp = numpy.matlib.zeros((m,k))
            for varidx in range(0,m):
                if deltam[iter,varidx] == cidx:
                    cnt += 1
                    betastmp[cnt-1,:] = betas[varidx,:]
            SS = betastmp[0:cnt,:]
            SS = SS.T*SS
            # Use prior as 50 variants
            SSsample = wishart.rvs(50 + cnt, SS + Omega0inv)
            thetainvdict[iter,cidx] = SSsample
            thetadict[iter,cidx] = numpy.linalg.inv(SSsample)
    ## Write output for input files
    mcout = open(outpath + str(fout) + '.mcmc.posteriors','w+')
    varprobdict = {}
    for varidx in range(0,m):
        mcout.write(chroffvec[varidx] + '\t' + annotvec[varidx] + '\t' + protvec[varidx] + '\t' + genevec[varidx] + '\t' + str(genevec[varidx] + ':' + annotvec[varidx] + ':' +  protvec[varidx]))
        for cidx in range(1,C+1):
            probclustervar = numpy.where(deltam[burn+1:niter+1,varidx] == cidx)[0].shape[0]/(niter - burn)
            varprobdict[chroffvec[varidx],cidx] = probclustervar
            mcout.write('\t'  + str(probclustervar))
        mcout.write('\n')
    mcout.close()
    maxllkiter = np.max(maxloglkiter[burn+1:niter:thinning,0])
    BIC = -2*maxllkiter + (k)*(C-1)*np.log(m)
    AIC = -2*maxllkiter + (k)*(C-1)*2
    returndict = {}
    returndict['bic'] = BIC
    returndict['aic'] = AIC
    returndict['thetainv'] = thetainvdict
    returndict['iter'] = iter
    returndict['c'] = C
    return returndict
#    return [BIC,AIC,thetainvdict,iter,C]
    if verbose:
        probout = fout + '.mcmc.probs'
        numpy.savetxt(outpath + probout, deltam, fmt='%1.3f')
        bcout = open(outpath + str(fout) +  '.mcmc.bc','w+')
        bcout.write('cluster')
        for i in range(0,len(phenidarr)):
            print(("\t%s\t%s\t%s") % (phenidarr[i] + 'm50',phenidarr[i] + 'l95', phenidarr[i] + 'u95'), end = '', file = bcout)
        bcout.write('\n')
        for cidx in range(0,C):
            mean = numpy.mean(bc[burn+1:niter+1:thinning,cidx,:],axis = 0)
            l95ci = numpy.percentile(bc[burn+1:niter+1:thinning,cidx,:],2.5, axis = 0)
            u95ci = numpy.percentile(bc[burn+1:niter+1:thinning,cidx,:],97.5, axis = 0)
            bcout.write(str(cidx))
            for phenidx in range(0,mean.shape[0]):
                print(("\t%2.2f\t%2.2f\t%2.2f") % (mean[phenidx], l95ci[phenidx], u95ci[phenidx]), end = '', file = bcout)
            bcout.write('\n')
        bcout.close()
        scaleout = open(outpath + str(fout) + '.mcmc.scale','w+')
        for annotidx in range(0,annotlen):
            mean = numpy.mean(np.sqrt(scales[burn+1:niter+1:thinning,annotidx]),axis = 0)
            l95ci = numpy.percentile(np.sqrt(scales[burn+1:niter+1:thinning,annotidx]),2.5, axis = 0)
            u95ci = numpy.percentile(np.sqrt(scales[burn+1:niter+1:thinning,annotidx]),97.5, axis = 0)
            print(("%s\t%s\t%2.2f\t%2.2f\t%2.2f") % (str(annotidx),annotmap[annotidx],mean,l95ci,u95ci), file = scaleout)
        scaleout.close()
        tmpbc = open(outpath + str(fout)  + '.theta.bc', 'w+')
        for jidx in range(0,k):
            for kidx in range(0,k):
                print(Theta0[jidx,kidx], file = tmpbc,end = ' ')
            print('\n',end='',file=tmpbc)
        tmpbc.close()
        pc[0,0,:]
        print('geneset', np.mean(pcj[burn+1:niter+1:thinning,:],axis=0))
        # initialize pcj (proportions for each gene j)
        genesdict = {}
        for geneidx in range(0,genenum):
            genesdict[genemap[geneidx]] = genemap[geneidx]
            genedatm50[genemap[geneidx]] = np.mean(pcj[burn+1:niter+1:thinning,geneidx,:],axis=0)
            genedatl95[genemap[geneidx]] = np.percentile(pcj[burn+1:niter+1:thinning,geneidx,:], 2.5, axis=0)
            genedatu95[genemap[geneidx]] = np.percentile(pcj[burn+1:niter+1:thinning,geneidx,:], 97.5, axis=0)
    alphaout = open(outpath + str(fout) + '.mcmc.alpha','w+')
    mean = numpy.mean(alpha[burn+1:niter+1:thinning,0],axis = 0)
    l95ci = numpy.percentile(alpha[burn+1:niter+1:thinning,0],2.5, axis = 0)
    u95ci = numpy.percentile(alpha[burn+1:niter+1:thinning,0],97.5, axis = 0)
    print(mean)
    print(("%2.2f\t%2.2f\t%2.2f") % (mean,l95ci,u95ci), file = alphaout)
    alphaout.close()
