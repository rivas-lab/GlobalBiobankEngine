data {
    int<lower=0> N; // number of samples
    int<lower=1> M; // dimension
    matrix[N, M] B; // data observed
    matrix[N, M] SE; // standard errors
    int<lower=1> K; // number of mixture components
}
transformed data{
    vector[M] zeros; 
    matrix[M,M] SE_mat[N];
    
    zeros = rep_vector(0, M);
    
    // fill in the matrix of standard errors
    for (n in 1:N) {
        SE_mat[n] = diag_matrix(to_vector(SE[n]));
    }
}

parameters {
    simplex[K] pi; // mixing proportions
    cholesky_factor_corr[M] L_Omega;
    vector<lower=0>[M] tau; 
}

transformed parameters{

    matrix[M, M] Sigma;
    matrix[M, M] Sigmas[K];
    Sigma = diag_pre_multiply(tau, L_Omega)*diag_pre_multiply(tau, L_Omega)';
    
    Sigmas[1] = diag_matrix(rep_vector(0,M));
    Sigmas[2] = Sigma;
}

model {
    vector[K] ps; // contributions of each
    tau ~ cauchy(0, 2.5);
    L_Omega ~ lkj_corr_cholesky(2);
    pi ~ dirichlet(rep_vector(1, K));

    for (n in 1:N){

       // two components
       for (k in 1:K){
           ps[k] = log(pi[k]) + multi_normal_lpdf(B[n] | zeros, SE_mat[n] + Sigmas[k]);
       }
       target += log_sum_exp(ps);
     }
 
}
