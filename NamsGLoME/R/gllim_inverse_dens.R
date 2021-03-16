#' @export
gllim_inverse_dens = function(y,theta,x,verb=0){
# %%%%%%%%%%%%%%%%% Inverse Mapping from Gllim Parameters %%%%%%%%%%%%%%%%%%%
# %%%% Author: Antoine Deleforge (July 2012) - antoine.deleforge@inria.fr %%%
# %% Converted to R: Emeline Perthame (2016) - perthame.emeline@gmail.com %%%
# %% Extented: TrungTin Nguyen (14-03-2021) - tinnguyen0495@gmail.com     %%%

# % Description:
# 1. Map N observations y using the inverse conditional
# expectation E[x|y;theta] of the gllim model with parameters theta.
# 2. Calculate Gaussian mixture parameters of the inverse
# conditional density p(x|y;theta) (i.e theta^star) in space Y using gllim
#  parameters theta.
# 3. Evaluate conditional log-likelihood with theta^star
# given sample (x_i,y_i)_{i = 1,...,N}.
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%% Input %%%%
# % - y (DxN)               % Input observations to map (covariate variable)
# % - x (LxN)               % Output observations (response variable)
# % for conditional log-likelihood
# % - theta  (struct)       % Gllim model parameters
# %   - theta.c (LxK)       % Gaussian means of X's prior
# %   - theta.Gamma (LxLxK) % Gaussian covariances of X's prior
# %   - theta.pi (1xK)      % Gaussian weights of X's prior
      # %   - theta.A (DxLxK)     % Affine transformation matrices
      # %   - theta.b (DxK)       % Affine transformation vectors
# %   - theta.Sigma (DxDxK) % Error covariances
# % - verb {0,1,2}          % Verbosity (def 1)
# %%%% Output %%%%
# % - x_exp (LxN)           % Posterior mean estimates E[xn|yn;theta]
# % - alpha (NxK)           % Weights of the posterior GMMs
# % thetaStar               % The estimated parameter for gllim_inverse
# % CNLL (1x1)               % Conditional negative-log-likelihood with theta^star
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	D = nrow(y) ; N = ncol(y) ;
	L = nrow(theta$c) ; K = ncol(theta$c) ;
	cs <- matrix(0,D,K); Gammas <- array(NaN,dim=c(D,D,K));
	As <- array(NaN,dim=c(L,D,K)); bs <- matrix(0,L,K); Sigmas <- array(NaN,dim=c(L,L,K));
	pXyz <- matrix(0,N,K)

# % ======================Inverse density parameters=========================
	if(verb>=1) print('Compute K projections to X space and weights')
	proj=array(NaN,dim=c(L,N,K));
  logalpha=matrix(0,N,K);

	for (k in 1:K){
    if(verb>=2) print(paste("k=",k,sep=""));
    if(verb>=2) print('AbcG ');
    if (L==1) {Ak=theta$A[,,k,drop=FALSE];} else {Ak=matrix(theta$A[,,k],ncol=L,nrow=D);} # % DxL
    bk=theta$b[,k]; # % Dx1
    Sigmak=theta$Sigma[,,k]; # %DxD  ## OK
    if (L==1) {ck=theta$c[,k,drop=FALSE];} else {ck=theta$c[,k];} # % Lx1
    Gammak=theta$Gamma[,,k]; # % LxL ## OK

    # the forward parameter vector theta^{star}:
    if(verb>=2) print('cks ');
    cks=Ak%*%ck+bk; cs[,k] <-cks; ## OK

    if(verb>=2) print('Gks ');
    Gammaks=Sigmak+Ak%*%tcrossprod(Gammak,Ak); Gammas[,,k] <- Gammaks ## OK

    if(verb>=2) print('iSks ');
    iSk = solve(Sigmak) #diag(1/diag(Sigmak))
    invSigmaks2=diag(L)+Gammak%*%crossprod(Ak, iSk)%*%Ak #Sigmaks=(invSigmaks2)^-1%*%Gammak
    Sigmaks <- solve(invSigmaks2,Gammak); Sigmas[,,k] <- Sigmaks

    if(verb>=2) print('Aks ');
    Aks=solve(invSigmaks2,Gammak)%*%crossprod(Ak, iSk); As[,,k] <- Aks;

    if(verb>=2) print('bks ');
    bks= solve(invSigmaks2) %*% (ck-Gammak%*%crossprod(Ak, iSk%*%bk)); bs[,k] <- bks; # Lx1
    # Note that bks = Sigmaks%*%Gammak^(-1) %*% %*% (ck-Gammak%*%crossprod(Ak, iSk%*%bk))
    # and Sigmaks%*%Gammak^(-1) = (invSigmaks2)^(-1)=solve(invSigmaks2) (definition of invSigmaks2).

    if(verb>=2) print('projections ');
    proj[,,k]=sweep(Aks%*%y,1,bks,"+"); # LxN

    if(verb>=2) print('logalpha ');
      logalpha[,k]=log(theta$pi[k])+loggausspdf(y,cks,Gammaks);
    if(verb>=2) print('');

    # Calculate probability p(x|y,Z=k;theta)
    pXyz[,k] <- exp(loggausspdf(x,matrix(proj[,,k],nrow = L),matrix(Sigmas[,,k],nrow = L))) # NxK
    #pXyz[,k] <- exp(loggausspdf(x,sweep(Aks%*%y,1,bks,"+"),Sigmas[,,k])) # NxK
	}

# % ===Conditional expectation E[x|y;theta] of the gllim model with parameters theta===========
	den=logsumexp(logalpha,1) # Nx1
	logalpha= sweep(logalpha,1,den,"-") # N xK
	alpha=exp(logalpha) # NxK

	x_exp = lapply(1:L,function(l) proj[l,,]*alpha) # list: Lx(NxK)
	x_exp = lapply(x_exp,rowSums) # Lx(Nx1)
	x_exp = do.call(rbind, x_exp) # LxN

# % ======================Conditional log-likelihood with theta^star===========================
	# Vector of pdf values over all sample \hat{s}_{\hat{m}}(Y_i|X_i), i = 1,..., n*ny.
	CLL_vec <- rowSums(pXyz*alpha)

	# Summation of conditional negative log-likelihood value over all samples.
	CNLL <- -sum(log(rowSums(pXyz*alpha)))

return(list(CNLL = CNLL, CLL_vec = CLL_vec, x_exp=x_exp,alpha=alpha,
            cs=cs,Gammas=Gammas,pis = theta$pi,As = As,bs = bs,Sigmas = Sigmas,proj = proj))
}
