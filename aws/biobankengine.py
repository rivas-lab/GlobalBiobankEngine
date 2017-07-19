from __future__ import print_function
import logging
import os
import sys
import math
import zipfile

import boto3
import numpy
import sklearn
from scipy.stats import invgamma, multivariate_normal, rv_discrete

logger = logging.getLogger('log')

s3_client = boto3.client('s3')

SAMPLE_DATA_PREFIX_PATH = '/tmp/t-mrp-lof_variant_missense_variant_RH160_HC153_IL33_IL10'


def open_files(base_name):
    '''
    Load files required for processing.
    '''
    file_postfixes = ('.betas.npy', '.se.npy', '.pvalues.npy', '.annotations.npy', '.protein_annotations.npy',
        '.variant_ids.npy', '.icd.npy', '.gene_return.npy', '.rsids.npy', '.alts.npy', '.allele_frequencies.npy')
    file_names = (base_name + postfix for postfix in file_postfixes)
    return (numpy.load(file_name) for file_name in file_names)

def is_pos_def(x):
    i = 0
    x = numpy.matrix(x)
    if numpy.all(numpy.linalg.eigvals(x) > 0):
        return True
    else:
       return False

def extract_zip(file, to_dir):
    zf = zipfile.ZipFile('/tmp/'+file, 'r')
    zf.extractall(to_dir)
    zf.close()

def test(event, context):
    print('Hello World')

def run_mrp_lambda(event, context):
    """
    Process a file upload.
    """

    # Get the uploaded file's information
    bucket = event['Records'][0]['s3']['bucket']['name'] # Will be `my-bucket`
    key = event['Records'][0]['s3']['object']['key'] # Will be the file path of whatever file was uploaded.

    # Get the bytes from S3
    s3_client.download_file(bucket, key, '/tmp/' + key) # Download this file to writable tmp space.
    extract_zip(key, '/tmp')
    try:
        os.mkdir('/tmp/MRP_out')
    except:
        pass
    print(run_mrp())

def run_mrp(lof=True, missense=True, genes=['RH160', 'HC153', 'IL33', 'IL10'], fdr=5, phenidarr = ['ICD1462','ICD1463']):
    fdr = int(fdr)/100
    annotations = []

    # Add in the selected annotations to the category list
    if lof:
        annotations.append('lof_variant')
    if missense:
        annotations.append('missense_variant')
    # If the input file has genes
    if genes != None:
        # Find variants with the given annotation
        key = 't-mrp-' + '_'.join(annotations) + '_' + '_'.join(phenidarr) + '_' + '_'.join(genes)
        # Generate relevant files
        betas, se, pvalues, annotations, protein_annotations, variant_ids, icd, gene_return, rsids, alts, allele_frequencies = open_files(SAMPLE_DATA_PREFIX_PATH)
        # Reshape betas from genome query
        C = numpy.matlib.eye(betas.shape[1])
        annotvec = [str(annotations[i].strip('"').strip('[').strip(']').strip("'")) for i in range(0,len(annotations))]

        # Run MRP with 2 clusters, output to the MRP_out subdirectory in the gbe_browser directory
        bicarr = []
        cmax = 4
        fail = 0
        for tmpc in range(1,cmax+1):
            [BIC, AIC, genedat] = mrpmm(betas,se, C, annotvec, gene_return, rsids, variant_ids,tmpc,key,C, numpy.linalg.inv(C),icd, fdr=fdr, niter=51,burn=10,thinning=1,verbose=True, outpath = '/tmp/MRP_out/')
            if tmpc > 2 and BIC > bicarr[len(bicarr)-1]:
                fail += 1
            if fail >= 2:
                break
            bicarr.append(BIC)

            print(tmpc,BIC,AIC)
#        if tmpc == 5:
#            with PyCallGraph(output=GraphvizOutput()):
#                mrpmm(betas,se, C, annotvec, gene_return, rsids, variant_ids,5,key,C, numpy.linalg.inv(C),icd, fdr=fdr, niter=201,burn=50,thinning=1,verbose=True, outpath = '/tmp/MRP_out/')
        cminarr = bicarr[1:]
        clustminidx = cminarr.index(min(cminarr))
        clustminval = clustminidx + 2
        clustmaxidx = cminarr.index(max(cminarr))
        clustmaxval = clustmaxidx + 2
        if bicarr[0] > 0 and bicarr[1] > 0:
            lbf2 = 0.5*(bicarr[0] - cminarr[0] )
            lbf = 0.5*( bicarr[0] - cminarr[clustminidx])
        else:
            lbf2 = 0.5*(bicarr[0] - cminarr[0])
            lbf = 0.5*(bicarr[0] - cminarr[clustminidx])
        print("bicstuff",bicarr[0],cminarr[clustminidx])
        lbfout = open("/tmp/MRP_out/" + key + ".lbf",'w')
        lbfout.write(str(lbf))
        lbfout.close()
        clustvalue = 1
        if lbf > 1 and (lbf - lbf2) > 1 and lbf2 > 1:
#        if lbf > 1:
            [BIC, AIC, genedat] = mrpmm(betas,se, C, annotvec, gene_return, rsids, variant_ids,clustminval,key,C, numpy.linalg.inv(C), icd, fdr=fdr, niter=201,burn=100,thinning=1,verbose=True, outpath = '/tmp/MRP_out/', protectivescan=True)
            print("log bayes factor: ",lbf)
            clustvalue = clustminval
        elif lbf2 > 1:
            clustminval = 2
            [BIC, AIC, genedat] = mrpmm(betas,se, C, annotvec, gene_return, rsids, variant_ids,clustminval,key,C, numpy.linalg.inv(C), icd, fdr=fdr, niter=201,burn=100,thinning=1,verbose=True, outpath = '/tmp/MRP_out/', protectivescan=True)
            print("log bayes factor: ",lbf)
            clustvalue = clustminval
        else:
            [BIC, AIC, genedat] = mrpmm(betas,se, C, annotvec, gene_return, rsids, variant_ids,1,key,C, numpy.linalg.inv(C), icd, fdr=fdr, niter=51,burn=10,thinning=1,verbose=True, outpath = '/tmp/MRP_out/', protectivescan=True)
        print(bicarr)
        print(icd)
        print("genedat",genedat)
    return([key, lbf, clustvalue, genedat])

def mrpmm(betas,ses,vymat,annotvec,genevec,protvec,chroffvec,clusters,fout,Rphen,Rpheninv,phenidarr, Rphenuse=True, fdr=.05, niter=1000,burn=100,thinning=1,verbose=True, protectivescan = False, outpath='/Users/mrivas/', maxlor = 0.693):
    print("Running MCMC algorithm...")
    print(sys.flags.optimize)
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
    maxloglkiter = numpy.zeros((niter+2,1))
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
        else:
            Theta0  = sklearn.covariance.shrunk_covariance(Rphen)
            Theta0inv = numpy.linalg.inv(Theta0)
    else:
        Theta0  = numpy.eye(Rphen.shape[0])
        Theta0inv = numpy.linalg.inv(Theta0)


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
    # store the probabilities (proportions) of cluster memberships for each gene
    genenum = len(set(genevec))
    pcj = numpy.zeros((niter+2,genenum,C))
    # for each iteration keep record of the variant membership
    deltam = numpy.zeros((niter+2,m))

    ###### Why are these stored separately?
    # non-normalized probabilities for each individual variant
    uc = numpy.zeros((niter+2,m,C))
    # normalized probabilities for each individual variant
    ws = numpy.zeros((niter+2,m,C))


    # for each iteration keep record of the variant membership
    tm = numpy.zeros((niter+2,m))
    #sharing parameter
    alpha = numpy.zeros((niter+2,1))
    ks = numpy.arange(1,C+1)
    sigmainvdict = {}
    sigmadict = {}
    thetadict = {}
    thetainvdict = {}
    for clusteriter in range(2,C+1):
        sigmadict[0,clusteriter] = S
        sigmainvdict[0,clusteriter] = Sinv
        thetadict[0,clusteriter] = Theta0
        thetainvdict[0,clusteriter] = Theta0inv
    # For Metropolois Hastings sub-step : keep track of acceptance rate
    acceptmh1 = 0
    rejectmh1 = 0
    acceptmh1_postburnin = 0
    rejectmh1_postburnin = 0
    acceptmh3 = 0
    rejectmh3 = 0
    acceptmh3_postburnin = 0
    rejectmh3_postburnin = 0
    acceptmh2 = [0]*annotlen
    rejectmh2 = [0]*annotlen
    acceptmh2_postburnin = [0]*annotlen
    rejectmh2_postburnin = [0]*annotlen
    # initialize \alpha : sharing of clusters across genes
    alpha[0,:] = invgamma.rvs(1,0,1,size = 1)
    # initialize pc (proportions across all variants)
    pc[0,0,:] = numpy.random.dirichlet([1]*C)
    # initialize pcj (proportions for each gene j)
    for geneidx in range(0,genenum):
        pcj[0,geneidx,:] = numpy.random.dirichlet(alpha[0,0]*pc[0,0,:])
    bc[0,0,:] = numpy.array([0]*k)
    for clusteridx in range(1,C):
        bc[0,clusteridx,:] = numpy.random.multivariate_normal(numpy.array([0]*k).T,Theta0)
    for scaleidx in range(0,annotlen):
        scales[0,scaleidx] = numpy.power(0.2,2)
    # initialize variant membership across clusters
    deltam[0,:] = numpy.random.randint(0,C,m)
    # protective candidate alleles
    protind = numpy.zeros((niter+2,m))
    # Iterations MCMC samplers
    for iter in range(1,niter+1):
        gamma = 1
        if iter % 500 == 0:
            print(iter)
        ## a) Update \pi_0 : Proposal centered around the current value, Set gamma to 1 , how to set gamma?
        ## mhstep1
        pcproposal = numpy.random.dirichlet(alpha[iter-1,0]*pc[iter-1,0,:])
      #  lnormDprop = math.lgamma(numpy.sum([alpha[iter-1,0]*i for i in pcproposal])) - numpy.sum([math.lgamma(max(alpha[iter-1,0]*i,epsilon)) for i in pcproposal])
        lnormDprop = math.lgamma(numpy.sum([gamma*i for i in pcproposal])) - numpy.sum([math.lgamma(max(gamma*i,epsilon)) for i in pcproposal])
        # second part of density
    #    densitypropb = numpy.sum([(alpha[iter-1,0]*pcproposal[i] - 1)*numpy.log(pc[iter-1,0,i]) for i in range(0,C)])
        densitypropb = numpy.sum([(gamma*pcproposal[i] - 1)*numpy.log(pc[iter-1,0,i]) for i in range(0,C)])
        lpdirprop = lnormDprop + densitypropb
        #go through each gene
        lpdirpropgene = 0
        lnormDprop = math.lgamma(numpy.sum([alpha[iter-1,0]*i for i in pcproposal])) - numpy.sum([math.lgamma(max(alpha[iter-1,0]*i,epsilon)) for i in pcproposal])
        for geneidx in range(0,genenum):
            # second part of density
            densitypropb = numpy.sum([(alpha[iter-1,0]*pcproposal[i] - 1)*numpy.log(pcj[iter-1,geneidx,i]) for i in range(0,C)])
            lpdirpropgene += densitypropb + lnormDprop
        lpdirnum = lpdirprop + lpdirpropgene
        # denominator, iteration - 1 pc
     #   lnormD = math.lgamma(numpy.sum([alpha[iter-1,0]*i for i in pc[iter-1,0,:]])) - numpy.sum([math.lgamma(max(alpha[iter-1,0]*i,epsilon)) for i in pc[iter-1,0,:]])
        lnormD = math.lgamma(numpy.sum([gamma*i for i in pc[iter-1,0,:]])) - numpy.sum([math.lgamma(max(gamma*i,epsilon)) for i in pc[iter-1,0,:]])
        # second part of density
        densityb = numpy.sum([(gamma*pc[iter-1,0,i] - 1)*numpy.log(pcproposal[i]) for i in range(0,C)])
        lpdir = lnormD + densityb
        #go through each gene
        lpdirgene = 0
        lnormD = math.lgamma(numpy.sum([alpha[iter-1,0]*i for i in pc[iter-1,0,:]])) - numpy.sum([math.lgamma(max(alpha[iter-1,0]*i,epsilon)) for i in pc[iter-1,0,:]])
        for geneidx in range(0,genenum):
            # second part of density
            densityb = numpy.sum([(alpha[iter-1,0]*pc[iter-1,0,i] - 1)*numpy.log(pcj[iter-1,geneidx,i]) for i in range(0,C)])
            lpdirgene += densityb + lnormD
        lpdirdenom = lpdir + lpdirgene
        lpdir = lpdirnum - lpdirdenom
        ## Metropolis-Hastings step
        if numpy.log(numpy.random.uniform(0,1,size = 1)[0]) < min(0, lpdir):
            acceptmh1 += 1
            pc[iter,0,:] = pcproposal
            if iter > burn:
                acceptmh1_postburnin += 1
        else:
            rejectmh1 += 1
            pc[iter,0,:] = pc[iter-1,0,:]
            if iter > burn:
                rejectmh1_postburnin += 1
        # b) For each gene j = 1, ..., J update \pi_j
        for geneidx in range(0,genenum):
            paramvecshared = alpha[iter-1,0]*pc[iter,0,:]
            for geneiter in range(0,len(genevec)):
                if genevec[geneiter] == genemap[geneidx]:
                    paramvecshared[int(deltam[iter-1,geneiter])] += 1
            pcj[iter,geneidx,:] = numpy.random.dirichlet(paramvecshared)
        # c) Update delta_jm
        xk = numpy.arange(0,C)
        for varidx in range(0,m):
            probmjc = [0]*C
            lprobmjcu = [0]*C
            uc = [0]*C
            varannot = annotvec[varidx]
            annotidx = [i for i in range(0,annotlen) if annotmap[i] == varannot][0]
            genevar = genevec[varidx]
            geneid = [i for i in range(0,len(genemap)) if genemap[i] == genevar][0]
            atmp = numpy.array(ses[varidx,:])[0]
            dtmp = numpy.matlib.eye(len(atmp))
            numpy.fill_diagonal(dtmp,atmp)
            Vjm = dtmp*S*dtmp + numpy.matlib.eye(S.shape[0])*.000001
            # Gives covariance matrix of variant  effect on sets of phenotypes (after fixed effect meta-analysis has been applied across all studies available)
            for cidx in range(0,C):
                llk2 = multivariate_normal.logpdf(betas[varidx,:],numpy.sqrt(scales[iter-1,annotidx])*bc[iter-1,cidx,:],Vjm) + numpy.log(pcj[iter,geneid,cidx])
                if deltam[iter-1,varidx] == cidx:
                    maxloglkiter[iter-1,0] += llk2
                lprobmjcu[cidx] += llk2
                        #normalize uc - set to wc
            maxloglk = numpy.max(lprobmjcu)
            for cidx in range(0,C):
                uc[cidx] = numpy.exp(lprobmjcu[cidx] - maxloglk)
            for cidx in range(0,C):
                probmjc[cidx] = uc[cidx]/numpy.sum(uc)
            if numpy.isnan(probmjc[0]):
                wstmp = numpy.random.dirichlet(numpy.repeat(numpy.array([1]),C,axis = 0))
                custm = rv_discrete(name='custm',values=(xk,wstmp))
            else:
                custm = rv_discrete(name='custm',values=(xk,probmjc))
            deltam[iter,varidx] = custm.rvs(size=1)[0]
            if protectivescan:
                protbool = 0
                protadverse = 0
                for tmptidx in range(0, k):
                    if numpy.sqrt(scales[iter-1,annotidx])*bc[iter-1,int(deltam[iter,varidx]),tmptidx] >= maxlor:
                        protadverse = 1
                        if numpy.sqrt(scales[iter-1,annotidx])*bc[iter-1,int(deltam[iter,varidx]),tmptidx] < -.1:
                            protbool = 1
                if protbool == 1 and protadverse == 0:
                    protind[iter,varidx] = 1
        # d) Update b_c using a Gibbs update from a Gaussian distribution
        for cidx in range(1,C):
            cnt = 0
            mucurrenttmp1 = 0
            varcurrenttmp1 = 0
            mucurrenttmp2 = 0*betas[0,:]
            mucurrenttmp2 = mucurrenttmp2.T
            for varidx in range(0,m):
                if deltam[iter,varidx] == cidx:
                    cnt += 1
                    if cnt == 1:
                        varannot = annotvec[varidx]
                        annotidx = [i for i in range(0,annotlen) if annotmap[i] == varannot][0]
                        atmp = numpy.array(ses[varidx,:])[0]
                        dtmp = numpy.matlib.eye(len(atmp))
                        numpy.fill_diagonal(dtmp,atmp)
                        Vjmtmp = dtmp*S*dtmp + numpy.matlib.eye(S.shape[0])*.000001
                        Vjminvtmp = numpy.linalg.inv(Vjmtmp)
                        mucurrenttmp1 = scales[iter-1,annotidx]*Vjminvtmp
                        mucurrenttmp2 = numpy.sqrt(scales[iter-1,annotidx])*Vjminvtmp*betas[varidx,:].T
                        varcurrenttmp1 = scales[iter-1,annotidx]*Vjminvtmp
                    else:
                        varannot = annotvec[varidx]
                        annotidx = [i for i in range(0,annotlen) if annotmap[i] == varannot][0]
                        atmp = numpy.array(ses[varidx,:])[0]
                        dtmp = numpy.matlib.eye(len(atmp))
                        numpy.fill_diagonal(dtmp,atmp)
                        Vjmtmp = dtmp*S*dtmp + numpy.matlib.eye(S.shape[0])*.000001
                        Vjminvtmp = numpy.linalg.inv(Vjmtmp)
                        mucurrenttmp1 += scales[iter-1,annotidx]*Vjminvtmp
                        mucurrenttmp2 += numpy.sqrt(scales[iter-1,annotidx])*Vjminvtmp*betas[varidx,:].T
                        varcurrenttmp1 += scales[iter-1,annotidx]*Vjminvtmp
            mucurrenttmp1 += Theta0inv
            varcurrenttmp1 += Theta0inv
            meanparam = numpy.ravel(numpy.linalg.inv(mucurrenttmp1)*mucurrenttmp2)
            varparam = numpy.linalg.inv(varcurrenttmp1)
            bc[iter,cidx,:] = numpy.random.multivariate_normal(meanparam,varparam)
        # e) Update scale sigma^2 annot.
        for annotidx in range(0,annotlen):
            scaleprop = abs(numpy.random.normal(numpy.sqrt(scales[iter-1,annotidx]),xi0,size = 1)[0])
            annotdata = annotmap[annotidx]
            probnum1 = invgamma.logpdf(numpy.power(scaleprop,2),1,scale=1)
            probdenom1 = invgamma.logpdf(scales[iter-1,annotidx],1,scale=1)
            lnum2 = 0
            ldenom2 = 0
            for varidx in range(0,m):
                if annotvec[varidx] == annotdata:
                    atmp = numpy.array(ses[varidx,:])[0]
                    dtmp = numpy.matlib.eye(len(atmp))
                    numpy.fill_diagonal(dtmp,atmp)
                    Vjm = dtmp*S*dtmp + numpy.matlib.eye(S.shape[0])*.000001
                    cidx = int(deltam[iter,varidx])
                   # print(cidx,iter,varidx)
                    lnum2 += multivariate_normal.logpdf(betas[varidx,:],scaleprop*bc[iter,cidx,:],Vjm)
                    ldenom2 += multivariate_normal.logpdf(betas[varidx,:],numpy.sqrt(scales[iter-1,annotidx])*bc[iter,cidx,:],Vjm)
            ## Metropolis-Hastings step
            if numpy.log(numpy.random.uniform(0,1,size = 1)[0]) < min(0, (lnum2 + probnum1) - (probdenom1 + ldenom2)):
                acceptmh2[annotidx] += 1
                scales[iter,annotidx] = numpy.power(scaleprop,2)
                if iter > burn:
                    acceptmh2_postburnin[annotidx] += 1
            else:
                rejectmh2[annotidx] += 1
                scales[iter,annotidx] = scales[iter-1,annotidx]
                if iter > burn:
                    rejectmh2_postburnin[annotidx] += 1
        # f) alpha
        alphaprop = abs(numpy.random.normal(alpha[iter-1,0],xialpha0,size = 1)[0])
        alphanum = -2*numpy.log(alphaprop) - 1/alphaprop
        alphadenom =  -2*numpy.log(alpha[iter-1,0]) - 1/alpha[iter-1,0]
        alphanum2 = 0
        alphadenom2 = 0
        lnormDprop = 0
        lpdirpropgene = 0
        lnormDliter = 0
        lpdirlgene = 0
        densitypropa = 0
        densitya = 0
        lnormDprop = math.lgamma(numpy.sum([alphaprop*i for i in pc[iter,0,:]])) - numpy.sum([math.lgamma(max(alphaprop*i,epsilon)) for i in pc[iter,0,:]])
        lnormDliter = math.lgamma(numpy.sum([alpha[iter-1,0]*i for i in pc[iter,0,:]])) - numpy.sum([math.lgamma(max(alpha[iter-1,0]*i,epsilon)) for i in pc[iter,0,:]])
        for geneidx in range(0,genenum):
            densitypropa = numpy.sum([(alphaprop*pc[iter,0,:] - 1)*numpy.log(pcj[iter-1,geneidx,i]) for i in range(0,C)])
            lpdirpropgene += densitypropa + lnormDprop
            densitya = numpy.sum([(alpha[iter-1,0]*pc[iter,0,:] - 1)*numpy.log(pcj[iter-1,geneidx,i]) for i in range(0,C)])
            lpdirlgene += densitya + lnormDliter
        ladirnum = alphanum + lpdirpropgene
        ladirdenom = alphadenom + lpdirlgene
        ladir = ladirnum - ladirdenom
        ## Metropolis-Hastings step
        if numpy.log(numpy.random.uniform(0,1,size = 1)[0]) < min(0, ladir):
            acceptmh3 += 1
            alpha[iter,:] = alphaprop
            if iter > burn:
                acceptmh3_postburnin += 1
        else:
            rejectmh3 += 1
            alpha[iter,:] = alpha[iter-1,:]
            if iter > burn:
                rejectmh3_postburnin += 1
    ## Write output for input files
    mcout = open(outpath + str(fout) + '.mcmc.posteriors','w+')
    varprobdict = {}
    for varidx in range(0,m):
        mcout.write(chroffvec[varidx] + '\t' + annotvec[varidx] + '\t' + protvec[varidx] + '\t' + genevec[varidx] + '\t' + str(genevec[varidx] + ':' + annotvec[varidx] + ':' +  protvec[varidx]))
        for cidx in range(0,C):
            probclustervar = numpy.where(deltam[burn+1:niter+1,varidx] == cidx)[0].shape[0]/(niter - burn)
            varprobdict[chroffvec[varidx],cidx + 1] = probclustervar
            mcout.write('\t'  + str(probclustervar))
        mcout.write('\n')
    mcout.close()
    ## Write output for protective scan
    if protectivescan:
        protout = open(outpath + str(fout) + '.mcmc.protective','w+')
        for varidx in range(0,m):
            protout.write(chroffvec[varidx] + '\t' + annotvec[varidx] + '\t' + protvec[varidx] + '\t' + genevec[varidx] + '\t' + str(genevec[varidx] + ':' + annotvec[varidx] + ':' +  protvec[varidx]))
            protdattmp = numpy.where(protind[burn+1:niter+1,varidx] == 1)[0].shape[0]/(niter - burn)
            protout.write('\t'  + str(protdattmp))
            protout.write('\n')
        protout.close()
    fdrout = open(outpath + str(fout) + '.fdr','w+')
    print(str(fdr),file = fdrout)
    varprobnull = []
    varfdrid = []
    for varidx in range(0,m):
        varfdrid.append(chroffvec[varidx])
        varprobnull.append(varprobdict[chroffvec[varidx],1])
    idxsort = sorted(range(len(varprobnull)), key=lambda k: varprobnull[k])
    varprobnullsort = [varprobnull[i] for i in idxsort]
    varfdridsort = [varfdrid[i] for i in idxsort]
    numfdrtmp = 0
    counter = 0
    varlfdr = []
    for i in range(0,len(varprobnullsort)):
        counter += 1
        numfdrtmp += varprobnullsort[i]
        fdrtmp = numfdrtmp/counter
        if fdrtmp <= fdr:
            print(varfdridsort[i], file = fdrout)
    fdrout.close()
    rejectionrate = rejectmh1_postburnin/(acceptmh1_postburnin + rejectmh1_postburnin)
    logger.info(("Your acceptance rate is %2.2f")  % ( rejectmh1_postburnin/(acceptmh1_postburnin + rejectmh1_postburnin)))
    genedatm50 = {}
    genedatl95 = {}
    genedatu95 = {}
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
            mean = numpy.mean(numpy.sqrt(scales[burn+1:niter+1:thinning,annotidx]),axis = 0)
            l95ci = numpy.percentile(numpy.sqrt(scales[burn+1:niter+1:thinning,annotidx]),2.5, axis = 0)
            u95ci = numpy.percentile(numpy.sqrt(scales[burn+1:niter+1:thinning,annotidx]),97.5, axis = 0)
            print(("%s\t%s\t%2.2f\t%2.2f\t%2.2f") % (str(annotidx),annotmap[annotidx],mean,l95ci,u95ci), file = scaleout)
        scaleout.close()
        tmpbc = open(outpath + str(fout)  + '.theta.bc', 'w+')
        for jidx in range(0,k):
            for kidx in range(0,k):
                print(Theta0[jidx,kidx], file = tmpbc,end = ' ')
            print('\n',end='',file=tmpbc)
        tmpbc.close()
        genesdict = {}
        for geneidx in range(0,genenum):
            genesdict[genemap[geneidx]] = genemap[geneidx]
            genedatm50[genemap[geneidx]] = numpy.percentile(pcj[burn+1:niter+1:thinning,geneidx,:],50, axis=0)
            genedatl95[genemap[geneidx]] = numpy.percentile(pcj[burn+1:niter+1:thinning,geneidx,:], 2.5, axis=0)
            genedatu95[genemap[geneidx]] = numpy.percentile(pcj[burn+1:niter+1:thinning,geneidx,:], 97.5, axis=0)
    alphaout = open(outpath + str(fout) + '.mcmc.alpha','w+')
    mean = numpy.mean(alpha[burn+1:niter+1:thinning,0],axis = 0)
    l95ci = numpy.percentile(alpha[burn+1:niter+1:thinning,0],2.5, axis = 0)
    u95ci = numpy.percentile(alpha[burn+1:niter+1:thinning,0],97.5, axis = 0)
    print(("%2.2f\t%2.2f\t%2.2f") % (mean,l95ci,u95ci), file = alphaout)
    alphaout.close()
    maxllkiter = numpy.max(maxloglkiter[burn+1:niter:thinning,0])
    BIC = -2*maxllkiter + (k+ genenum)*(C-1)*numpy.log(m)
    AIC = -2*maxllkiter + (k+ genenum)*(C-1)*2
    geneout = open(outpath + str(fout) + '.mcmc.gene.posteriors','w+')
    for genekey in genesdict.keys():
        print(genekey, file = geneout, end = '')
        for i in range(0,len(genedatm50[genekey])):
            print(("\t%2.2f") % (genedatm50[genekey][i]), file = geneout, end = '')
        for i in range(0,len(genedatl95[genekey])):
            print(("\t%2.2f") % (genedatl95[genekey][i]), file = geneout, end = '')
        for i in range(0,len(genedatu95[genekey])):
            print(("\t%2.2f") % (genedatu95[genekey][i]), file = geneout, end = '')
        geneout.write("\n")
    geneout.close()
    return [BIC,AIC,genedatm50]
