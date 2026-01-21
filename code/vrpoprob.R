
# variable response propensity + oprobit, extension of Peress (2010)
# reference: https://www.jstor.org/stable/27920175


# notation:
# M: number of outcome categories
# R: number of response categories
# J: number of vars in x
# K: number of vars in z
# y: ordered outcome, scalar \in {1...M}
# r: response status, scalar \in {1...R}
# x[1,2...J]: explanatory variable for y (row vector)
# z[1,2...K]: explanatory variable for r (row vector)
# alpha[1...J]: coefficients for x
# beta[1...K]: coefficients for z
# lambda[1...(M-1)]: thresholds for outcome oprobit, includes fixed lambda[M-1] = 0
# theta[1...R]: thresholds for response oprobit, includes fixed theta[R] = 0
# rho \in (-1,1): correlation of shocks
# xi[1...(J+K+(M-2)+(R-1)+1)]: vector of transformed/packed parameters
# Xpop[1...NXpop,1...J]: matrix of x values in population
# WXpop[1...NXpop]: populations weights for Xpop
# Zpop[1...NZpop,1...J]: matrix of z values in population
# WZpop[1...NZpop]: populations weights for Zpop
# Nmiss: number of nonresponses
# ydata[1...Nobs]: vector of observed y
# rdata[1...Nobs]: vector of observed r
# xdata[1...Nobs,1...J]: matrix of observed x
# zdata[1...Nobs,1...K]: matrix of observed z


# dependencies ------------------------------------------------------------

library(numDeriv)
library(mnorm)
library(maxLik)


# main functions ----------------------------------------------------------

vrpoprob_estim = function(ydata, rdata, xdata, zdata, Nmiss, WXpop, Xpop, WZpop, Zpop) {
  # main estimation routine
  # returns list with elements: alpha, beta, lambda, theta, rho, alpha_se, beta_se
  # also pphat, pphat_se are estimates of population proportions for outcomes
  
  # get dimensions
  J = ncol(xdata)
  K = ncol(zdata)
  M = max(ydata)
  R = max(rdata)
  
  # initial guess
  xi0 = rep(0, J+K+R+M-2)
  
  # call maxLik
  objfun = function(xi) vrpoprob_loglik_vec(xi, ydata, rdata, xdata, zdata, Nmiss, WZpop, Zpop)
  sol = maxLik(objfun, start=xi0, method="BFGS", control=list(printLevel=1))
  converged = (returnCode(sol)==0)  # this for BFGS, other methods may have different codes
  if (!converged) {warning("no convergence")} 
  
  # get pphat and its s.e.
  pphat = vrpoprob_xi_to_pphat(coef(sol), WXpop, Xpop, J, K, M, R)
  pphat_se = vrpoprob_delta_se(function(u) vrpoprob_xi_to_pphat(u, WXpop, Xpop, J, K, M, R), 
                               coef(sol), vcov(sol))
  
  # get pphat_nonresp its s.e.
  pphat_nonresp = vrpoprob_xi_to_pphat_resp_nonresp(coef(sol), WXpop, Xpop, Zpop, J, K, M, R)$pphat_nonresp
  pphat_nonresp_se = vrpoprob_delta_se(function(u) vrpoprob_xi_to_pphat_resp_nonresp(u, WXpop, Xpop, Zpop, J, K, M, R)$pphat_nonresp,
                                      coef(sol), vcov(sol))
  # get pphat_resp its s.e.
  pphat_resp = vrpoprob_xi_to_pphat_resp_nonresp(coef(sol), WXpop, Xpop, Zpop, J, K, M, R)$pphat_resp
  pphat_resp_se = vrpoprob_delta_se(function(u) vrpoprob_xi_to_pphat_resp_nonresp(u, WXpop, Xpop, Zpop, J, K, M, R)$pphat_resp,
                                       coef(sol), vcov(sol))
  
  # output
  out = list()
  out[["alpha"]] = vrpoprob_unpack(coef(sol), J, K, M, R)[["alpha"]]
  out[["beta"]] = vrpoprob_unpack(coef(sol), J, K, M, R)[["beta"]]
  out[["lambda"]] = vrpoprob_unpack(coef(sol), J, K, M, R)[["lambda"]]
  out[["theta"]] = vrpoprob_unpack(coef(sol), J, K, M, R)[["theta"]]
  out[["rho"]] = vrpoprob_unpack(coef(sol), J, K, M, R)[["rho"]]
  out[["alpha_se"]] = vrpoprob_unpack(stdEr(sol), J, K, M, R)[["alpha"]]
  out[["beta_se"]] = vrpoprob_unpack(stdEr(sol), J, K, M, R)[["beta"]]
  out[["pphat"]] = pphat
  out[["pphat_se"]] = pphat_se
  out[["pphat_nonresp"]] = pphat_nonresp
  out[["pphat_nonresp_se"]] = pphat_nonresp_se
  out[["pphat_resp"]] = pphat_resp
  out[["pphat_resp_se"]] = pphat_resp_se
  out[["converged"]] = converged
  out
}


# internal functions ------------------------------------------------------

vrpoprob_delta_se = function(f, x, V) {
	# quick delta method for standard errors
	Df = jacobian(f, x)  # from numDeriv
	Vf = Df %*% V %*% t(Df)
	sqrt(diag(Vf))
}


vrpoprob_pack = function(alpha, beta, lambda, theta, rho) {
  # pack parameters
  xi = c(alpha, beta, log(diff(lambda)), log(diff(theta)), atanh(rho/0.99))
}


vrpoprob_unpack = function(xi, J, K, M, R) {
  # unpack parameters
  alpha = xi[1:J]
  beta = xi[(J+1):(J+K)]
  lambda = rep(0, M-1)
  if (M>2) {
    for (i in (M-2):1) {
      lambda[i] = lambda[i+1] - exp(xi[J+K+i])
    }
  }
  theta = rep(0, R)
  for (i in (R-1):1) {
    theta[i] = theta[i+1] - exp(xi[J+K+(M-2)+i])
  }
  rho = 0.99*tanh(xi[J+K+M+R-2])
  list(alpha=alpha, beta=beta, lambda=lambda, theta=theta, rho=rho)
}


vrpoprob_loglik_vec = function(xi, ydata, rdata, xdata, zdata, Nmiss, WZpop, Zpop) {
	# compute log-likelihood, vectorized
	
	# get dimensions
	J = ncol(xdata)
	K = ncol(zdata)
	M = max(ydata) #?
	R = max(rdata) #?
	Nobs = length(ydata)
	
	# unpack parameters
	parms = vrpoprob_unpack(xi, J, K, M, R)
	alpha = parms$alpha
	beta = parms$beta
	lambda = parms$lambda
	theta = parms$theta
	rho = parms$rho
	
	# predictions for latent vars
	ystarhat = xdata %*% alpha
	rstarhat = zdata %*% beta
	
	# compute bounds A_epsilon < epsilon < B_epsilon
	LAMBDA = c(-Inf, lambda, Inf)
	A_epsilon = LAMBDA[(ydata-1)+1] - ystarhat
	B_epsilon = LAMBDA[(ydata)  +1] - ystarhat
	
	# compute bounds A_eta < eta < B_eta
	THETA = c(-Inf, theta, Inf)
	A_eta = THETA[(rdata-1)+1] - rstarhat
	B_eta = THETA[(rdata)  +1] - rstarhat
	
	# loglikelihood
	Sigma = matrix(1.0, nrow=2, ncol=2)
	Sigma[1,2] = rho
	Sigma[2,1] = rho
	ll = pmnorm(lower=cbind(A_epsilon,A_eta), upper=cbind(B_epsilon, B_eta), mean=c(0.0,0.0), sigma=Sigma, log=TRUE)
	
	# evaluate unconditional probability of nonresponse
	p_nonresp = sum(WZpop * pnorm(-(Zpop %*% beta), lower.tail=FALSE))
	
	# output
	sum(ll$prob) + Nmiss*log(p_nonresp)
}


vrpoprob_xi_to_pphat = function(xi, WXpop, Xpop, J, K, M, R) {
  # compute population proportions for outcome from parameters
  parms = vrpoprob_unpack(xi, J, K, M, R)
  alpha = parms$alpha
  lambda = parms$lambda
  ystarhat = as.matrix(Xpop) %*% alpha
  P = matrix(0, nrow(Xpop), M)
  P[,1] = pnorm(lambda[1] - ystarhat)
  if (M>2) {
    for (i in 2:(M-1)) {
      P[,i] = pnorm(lambda[i] - ystarhat) - pnorm(lambda[i-1] - ystarhat)
    }
  }
  P[,M] = pnorm(lambda[M-1] - ystarhat, lower.tail=FALSE)
  pphat = colSums(P * WXpop)
  pphat
}


vrpoprob_xi_to_pphat_resp_nonresp = function(xi, Wpop, Xpop, Zpop, J, K, M, R){
  # compute population proportions for nonresponders from parameters
  # Note: in general, this requires joint distribution of (X,Z) over population,
  #       whereas previous code required only marginal distributions which could
  #       be defined over different partitions of population. For our application,
  #       here we assume just a single partition that captures both X and Z.
  parms = vrpoprob_unpack(xi, J, K, M, R)
  alpha = parms$alpha
  lambda = parms$lambda
  ystarhat = as.matrix(Xpop) %*% alpha
  
  # we also need response variables
  beta = parms$beta
  rstarhat = as.matrix(Zpop) %*% beta
  
  rho = parms$rho
  
  Sigma = matrix(1.0, nrow=2, ncol=2)
  Sigma[1,2] = rho
  Sigma[2,1] = rho
  
  P_nonresp = matrix(0, nrow(Xpop), M)
  P_resp    = matrix(0, nrow(Xpop), M)
  
  for (j in 1:nrow(Xpop)){
    
    # category 1
    p1 = pmnorm(lower=c(-Inf, -Inf),
                 upper=c(lambda[1] - ystarhat[j], 0 - rstarhat[j]),
                 mean=c(0,0), sigma=Sigma)$prob
    
    p0 = pmnorm(lower=c(-Inf, 0 - rstarhat[j]),
                 upper=c(lambda[1] - ystarhat[j], Inf),
                 mean=c(0,0), sigma=Sigma)$prob
    
    if (!is.finite(p0) | p0 <= 0.0) p0 = 1e-16
    if (!is.finite(p1) | p1 <= 0.0) p1 = 1e-16
    
    P_nonresp[j,1] = p0
    P_resp[j,1]    = p1
    
    if (M>2) {
      for (i in 2:(M-1)) {
        
        p1 = pmnorm(lower=c(lambda[i-1] - ystarhat[j], -Inf),
                     upper=c(lambda[i] - ystarhat[j], 0 - rstarhat[j]),
                     mean=c(0,0), sigma=Sigma)$prob
        
        p0 = pmnorm(lower=c(lambda[i-1] - ystarhat[j], 0 - rstarhat[j]),
                     upper=c(lambda[i] - ystarhat[j], Inf),
                     mean=c(0,0), sigma=Sigma)$prob
        
        if (!is.finite(p0) | p0 <= 0.0) p0 = 1e-16
        if (!is.finite(p1) | p1 <= 0.0) p1 = 1e-16
        
        P_nonresp[j,i] = p0
        P_resp[j,i]    = p1
      }
    }
    
    # category M
    p1 = pmnorm(lower=c(lambda[M-1] - ystarhat[j], -Inf),
                 upper=c(Inf, 0 - rstarhat[j]),
                 mean=c(0,0), sigma=Sigma)$prob
    
    p0 = pmnorm(lower=c(lambda[M-1] - ystarhat[j], 0 - rstarhat[j]),
                 upper=c(Inf, Inf),
                 mean=c(0,0), sigma=Sigma)$prob
    
    if (!is.finite(p0) | p0 <= 0.0) p0 = 1e-16
    if (!is.finite(p1) | p1 <= 0.0) p1 = 1e-16
    
    P_nonresp[j,M] = p0
    P_resp[j,M]    = p1
  }
  
  Pr_resp    = pnorm(-rstarhat)
  Pr_nonresp = 1-pnorm(-rstarhat)
  
  P_nonresp = P_nonresp / matrix(Pr_nonresp,
                                 nrow = nrow(P_nonresp),
                                 ncol = ncol(P_nonresp))
  
  P_resp = P_resp / matrix(Pr_resp,
                           nrow = nrow(P_resp),
                           ncol = ncol(P_resp))
  
  pphat_nonresp = colSums(P_nonresp * Wpop)
  pphat_resp    = colSums(P_resp    * Wpop)
  
  list(pphat_nonresp = pphat_nonresp,
       pphat_resp    = pphat_resp)
}

