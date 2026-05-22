
# variable response propensity + oprobit, extension of Peress (2010)
# reference: https://www.jstor.org/stable/27920175


# dependencies ------------------------------------------------------------

library(numDeriv)
library(mnorm)
library(maxLik)


# main functions ----------------------------------------------------------

vrpoprob_estim = function(ydata, rdata, xdata, zdata, Nmiss, Wpop, Xpop, Zpop) {
  # Estimate model parameters and model-implied proprtions
  # given data, nonresponse number and population distribution of covariates.
  #
  # INPUTS
  # * ydata: vector of observed y 
  # * rdata: vector of observed r (same length as ydata)
  # * xdata: matrix of observed x (nrow=nrow(ydata), ncol is nr. of X covariates)
  # * zdata: matrix of observed z (nrow=nrow(ydata), ncol is nr. of Z covariates)
  # * Nmiss: number of nonresponses, integer
  # * Wpop: vector of population weights for x, z
  # * Xpop: matrix of population values for x (nrow=length(Wpop), ncol=ncol(xdata))
  # * Zpop: matrix of population values for z (nrow=length(Wpop), ncol=ncol(zdata))
  #
  # OUTPUT: list with elements:
  # * alpha, beta, alpha_se, beta_se: covariate coefficients and their s.e., vectors
  # * rho, rho_se: correlation parameter and its s.e., scalars
  # * lambda, theta: threshold estimates, vectors
  # * pphat, pphat_se: implied population outcome probabilities and s.e., vectors
  # * pphat_nonresp, pphat_nonresp_se: outcome probabilities and s.e. for nonresponders, vectors
  # * pphat_resp, pphat_resp_se: outcome probabilities and s.e. for responders, vectors
  # * converged: TRUE if optimization converged, else FALSE, logical
  
  # get dimensions
  J = ncol(xdata)
  K = ncol(zdata)
  M = max(ydata)
  R = max(rdata)
  
  # multiple starts in the rho direction to avoid spurious local optima
  # three starting values: rho in {-0.5, 0, 0.5}; all other params start at 0
  rho_idx  <- J + K + M + R - 2
  objfun   <- function(xi) vrpoprob_loglik(xi, ydata, rdata, xdata, zdata, Nmiss, Wpop, Zpop)
  xi0_base <- rep(0, J+K+R+M-2)
  rho_starts <- c(-0.5, 0, 0.5)

  sols <- lapply(rho_starts, function(r0) {
    xi0 <- xi0_base
    xi0[rho_idx] <- r0
    maxLik(objfun, start=xi0, method="BFGS", control=list(printLevel=0))
  })

  sol       <- sols[[which.max(vapply(sols, logLik, numeric(1)))]]
  converged <- (returnCode(sol) == 0)
  if (!converged) warning("no convergence")

  # covariance matrix: exact inverse of -H; NA SEs if H is non-finite or singular
  # Recompute Hessian via finite differences on the summed log-likelihood;
  # maxLik's BFGS-approximated sol$hessian can be rank-deficient and yield zero SEs.
  H <- tryCatch(
    numDeriv::hessian(function(xi) sum(objfun(xi)), coef(sol)),
    error = function(e) sol$hessian
  )
  H <- (H + t(H)) / 2
  V <- if (any(!is.finite(H))) {
    NULL
  } else {
    tryCatch(solve(-H), error = function(e) NULL)
  }
  se_ok <- !is.null(V)

  # point estimates (always available)
  parms  <- vrpoprob_unpack(coef(sol), J, K, M, R)
  pphat  <- vrpoprob_xi_to_pphat(coef(sol), Wpop, Xpop, J, K, M, R)
  rnonresp <- vrpoprob_xi_to_pphat_resp_nonresp(coef(sol), Wpop, Xpop, Zpop, J, K, M, R)

  # standard errors (NA when V unavailable)
  if (se_ok) {
    xi_se            <- sqrt(pmax(diag(V), 0))
    # delta method: rho = 0.99*tanh(xi), d(rho)/d(xi) = 0.99*(1-(rho/0.99)^2)
    rho_se           <- xi_se[rho_idx] * 0.99 * (1 - (parms$rho / 0.99)^2)
    pphat_se         <- vrpoprob_delta_se(function(u) as.numeric(vrpoprob_xi_to_pphat(u, Wpop, Xpop, J, K, M, R)), coef(sol), V)
    pphat_nonresp_se <- vrpoprob_delta_se(function(u) as.numeric(vrpoprob_xi_to_pphat_resp_nonresp(u, Wpop, Xpop, Zpop, J, K, M, R)$pphat_nonresp), coef(sol), V)
    pphat_resp_se    <- vrpoprob_delta_se(function(u) as.numeric(vrpoprob_xi_to_pphat_resp_nonresp(u, Wpop, Xpop, Zpop, J, K, M, R)$pphat_resp),    coef(sol), V)
  } else {
    xi_se            <- rep(NA_real_, length(coef(sol)))
    rho_se           <- NA_real_
    pphat_se         <- rep(NA_real_, M)
    pphat_nonresp_se <- rep(NA_real_, M)
    pphat_resp_se    <- rep(NA_real_, M)
  }

  # output
  out = list()
  out[["alpha"]]           <- parms$alpha
  out[["beta"]]            <- parms$beta
  out[["lambda"]]          <- parms$lambda
  out[["theta"]]           <- parms$theta
  out[["rho"]]             <- parms$rho
  out[["rho_se"]]          <- rho_se
  out[["alpha_se"]]        <- vrpoprob_unpack(xi_se, J, K, M, R)[["alpha"]]
  out[["beta_se"]]         <- vrpoprob_unpack(xi_se, J, K, M, R)[["beta"]]
  out[["pphat"]]           <- pphat
  out[["pphat_se"]]        <- pphat_se
  out[["pphat_nonresp"]]   <- rnonresp$pphat_nonresp
  out[["pphat_nonresp_se"]]<- pphat_nonresp_se
  out[["pphat_resp"]]      <- rnonresp$pphat_resp
  out[["pphat_resp_se"]]   <- pphat_resp_se
  out[["converged"]]       <- converged
  out
}


# internal functions ------------------------------------------------------


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


vrpoprob_delta_se = function(f, x, V) {
	# quick delta method for standard errors
	# method="simple" uses p+1 evaluations vs ~144 for Richardson default
	Df = jacobian(f, x, method="simple")
	Vf = Df %*% V %*% t(Df)
	sqrt(pmax(diag(Vf), 0))
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


vrpoprob_loglik = function(xi, ydata, rdata, xdata, zdata, Nmiss, WZpop, Zpop) {
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
  # compute population proportions for nonresponders from parameters.
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



