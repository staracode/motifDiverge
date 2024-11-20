
#############################################
#- 1. Same length of the two series of trials
#############################################


dcbern.sl = function(k,n,p10,p01,err=1E-20,g.approx=TRUE,details=FALSE){
########################################################################
	
	RES = vector(length=4,mode="list")
	names(RES) = c("value","error","nterms","g.approx")

	#- special easy case: no hits
	if((p10 ==0) & (p01==0)) {
		if(k==0){
			RES$value    = 1
			RES$error    = 0
			RES$nterms   = 0
			RES$g.approx = 0
			return(RES)
		} 
		RES$value    = 0
                RES$error    = 0
                RES$nterms   = 0
                RES$g.approx = 0
                return(RES)
	}
 
	if(!g.approx) {
		res          = rcpp_prob_sym(k,n,p10,p01,err)
		RES$value    = res[1]
		RES$error    = res[2]
		RES$nterms   = res[3]
		RES.g.approx = FALSE
	}
	else {
		if(cbern_sym_approx(n,p10,p01)){
			pars         = cbern_sym_approx(n,p10,p01,pars=TRUE)
			RES$value    = pnorm(k+.5,mean=pars[1],sd=pars[2]) - pnorm(k-.5,mean=pars[1],sd=pars[2])
			RES$error    = NA
			RES$nterms   = NA
			RES$g.approx = TRUE
		} else {
			res          = rcpp_prob_sym(k,n,p10,p01,err)
                	RES$value    = res[1]
                	RES$error    = res[2]
                	RES$nterms   = res[3]
                	RES$g.approx = FALSE
		}

	}
        if(details) return(RES)
	return(RES$value)
}


#################################################
#- 2. Differnt length of the two series of trials
#################################################


dcbern.dl = function(n1,n2,p10,p01,p11,g.approx=TRUE){
######################################################

	#- no k, because we always get whole range

	#- special easy case: no hits
	if((p10 ==0) & (p01==0)){
		RES           = cbind(0,1)
		colnames(RES) = c("k","p.k")
		return(RES)
	}
        
	#- EQUAL LENGTH:
	#===============
	if(n1==n2){
		#- TODO: shorten domain
		kmin = -n1
		kmax = +n1
		#- May be slow:
		res  = sapply(kmin:kmax,function(k) dcbern.sl(k,n1,p10,p01,g.approx))
		return(cbind(kmin:kmax,res))
	}

	#- DIFFERENT LENGTH:
	#===================

	nmin = min(n1,n2) #- smaller number of trials
        nmax = max(n1,n2) #- larger number of trials
        ndif = abs(n1-n2) #- binom number of trials
        p    = ifelse(n1>n2,p10+p11,p01+p11) #- binom success

	#- 1. Approximation?
	#-------------------
	if(g.approx){
		#- TODO: shorten domain
		if(cbern_asym_approx(n1,n2,p10,p01,p11)){
			pars = cbern_asym_approx(n1,n2,p10,p01,p11,pars=TRUE)
			rnge = (-nmin):nmax
			res  = pnorm(rnge+0.5,mean=pars[1],sd=pars[2]) - pnorm(rnge-0.5,mean=pars[1],sd=pars[2])
			fac  = ifelse(n1>n2,1,-1)
			rnge = fac*rnge
			ord  = order(rnge)
			return(cbind(rnge[ord],res[ord]))
		}			
	}	

	#- 2. Convolution!
	#-----------------

	rnge.ps = (-nmin):nmin 	#- symmetric part
	rnge.bn = 0:ndif	#- binomial part
	
	#- Can we shorten the ranges?
	#............................

	PMIN = 1E-20 #- min prob

	rnge.bn.max = rcpp_bin_quantile_rec(ndif,p,0,ndif,PMIN, 1)	
	rnge.bn.min = rcpp_bin_quantile_rec(ndif,p,0,ndif,PMIN,-1)	
	rnge.bn     = rnge.bn.min:rnge.bn.max

	#- do p-sym nested intervals here by hand
	#........................................
	
	#- find mode
	psfu <- function(x) dcbern.sl(x,nmin,p10,p01)
	curr.x = round(nmin*(p10-p01))
	curr.y = psfu(curr.x)
	while(1){
		if(psfu(curr.x+1) > psfu(curr.x)) {
			curr.x = curr.x+1
			next
		}
		if(psfu(curr.x-1) > psfu(curr.x)) {
                        curr.x = curr.x-1
                        next
                }
		break
	}	
	ps.mode = curr.x

	#- find maximum support
	cmin = ps.mode 
	cmax = nmin
	while( cmin + 1< cmax){
		cmid = ceiling(cmin + (cmax-cmin)/2)
		if(dcbern.sl(cmid,nmin,p10,p01,g.approx=FALSE) >= PMIN){
			cmin = cmid
		} else {
			cmax = cmid
		}
	}
	rnge.ps.max = cmax
	
	#- find minimum support
	cmin = -nmin
        cmax = ps.mode  
        while( cmin + 1< cmax){
                cmid = floor(cmin + (cmax-cmin)/2)
                if(dcbern.sl(cmid,nmin,p10,p01,g.approx=FALSE) >= PMIN){
                        cmax= cmid
                } else {
                        cmin = cmid
                }
        }
	rnge.ps.min = cmin
	rnge.ps     = rnge.ps.min:rnge.ps.max

	#- finally do the work
	#.....................

	vec.ps = sapply(rnge.ps,function(x) dcbern.sl(x,nmin,p10,p01))
	vec.bn = dbinom(rnge.bn,size=ndif,prob=p)

	if(n1>n2){
		res    = rcpp_convolve(vec.ps,vec.bn) #- look up this line 
		k.res  = (rnge.ps.min+rnge.bn.min):(rnge.ps.max+rnge.bn.max)	

	}
	if(n2>n1){
		res    = rcpp_convolve(vec.ps,rev(vec.bn)) #- look up this line
		k.res  = (rnge.ps.min-rnge.bn.max):(rnge.ps.max-rnge.bn.min)

	}
	res = res/sum(res) #- usually small error ~ 1E-12 because (1) p-sym too small and (2) cut terms
	
	r.res = cbind(k.res,res)
	colnames(r.res) = c("k","p.k")
	return(r.res)

}

################################
#- 3. Summarize to one interface
################################


dcbern = function(k, n1, n2, p10, p01, p11,g.approx=TRUE){
##########################################################

	#- FIXME: k needs to be kmin:kmax type for n1 != n2

	if( 	length(n1)  > 1 ||
		length(n2)  > 1 ||
		length(p10) > 1 ||
		length(p01) > 1 ||
		length(p11) > 1  )   stop("Illegal arguments. Only k can be a vector.\n" )
	 
	isnt.anint = function(x, tol = .Machine$double.eps^0.5) !(abs(x - round(x)) < tol)
	
	if( 	any(isnt.anint(k)) ||
		isnt.anint(n1)     ||	
		isnt.anint(n2)      ) stop("k, n1 and n2 need to be interger-valued.\n"   )	
	
	if( p10<0 || p01<0 || p11<0) stop("Negative probabilities in arguments\n"        )
	if( p10>1 || p01>1 || p11>1) stop("Probabilities > 1\n in arguments.\n"          )
	if( abs(1-p10-p01-p11) < 0 ) stop("Probabilities in arguments sum to > 1.\n"     )

	#- special easy case: no hits
	if((p10 ==0) & (p01==0) & (length(k) ==1)){
		if(k==0) return(1)
		return(0)
	}
 
	if(n1==n2){
		if(length(k)==1){
		if(k > n1)  return(NA)
		if(k < -n1) return(NA)
			return(dcbern.sl(k,n1,p10,p01))
		} else {
			return(sapply(k,function(x)dcbern.sl(x,n1,p10,p01,g.approx=g.approx)))
		}
	} else 	{
		tmp = dcbern.dl(n1,n2,p10,p01,p11,g.approx=g.approx)
		#- pick out the ks queried
		min.k = tmp[1,1]
		max.k = tmp[nrow(tmp),1]
		afu <- function(x){ 
			if(x>max.k) return(0)
			if(x<min.k) return(0)
			return(tmp[x-min.k+1,2])
		}
		res = sapply(k,afu)	
		return(res)
	}
	
}

pcbern.single = function(k, n1, n2, p10, p01, p11, lower.tail=TRUE,g.approx=TRUE){
##################################################################################

        nmin = min(n1,n2)
        nmax = max(n1,n2)
	
	#- easy specieal case
	if((p10==0) & (p01==0)){
		if(k==0) return(1)
		return(0)
	}	

        if(n1>=n2){
                kmin = -nmin
                kmax =  nmax
        }
        else{
                kmin = -nmax
                kmax = nmin
        }

	#- can't get that low
	if( k < kmin) return(NA)
	if( k > kmax) return(NA)

	#- handle gaussian approximation here 
	if(g.approx){
		if(n1 == n2){
			if(cbern_sym_approx(n1,p10,p01)){
				pars = cbern_sym_approx(n1,p10,p01,pars=TRUE)
				if(lower.tail) return(pnorm(k+0.5,mean=pars[1],sd=pars[2]))
				return(pnorm(k-0.5,mean=pars[1],sd=pars[2],lower.tail=FALSE))
			}
		}
		if(n1 != n2){
			if(cbern_asym_approx(n1,n2,p10,p01,p11)){
				pars = cbern_asym_approx(n1,n2,p10,p01,p11,pars=TRUE)
				if(lower.tail) return(pnorm(k+0.5,mean=pars[1],sd=pars[2]))
				return(pnorm(k-0.5,mean=pars[1],sd=pars[2],lower.tail=FALSE))
			}
		}
	}

	#- no approximation - use densities and sum
	if(lower.tail) return(sum(dcbern(kmin:k,n1,n2,p10,p01,p11,g.approx=g.approx)))	
	return(sum(dcbern(k:kmax,n1,n2,p10,p01,p11,g.approx=g.approx)))

}

pcbern = function(k, n1, n2, p10, p01, p11,lower.tail=TRUE,g.approx=TRUE){
##########################################################################

	if( 	length(n1)  > 1 ||
		length(n2)  > 1 ||
		length(p10) > 1 ||
		length(p01) > 1 ||
		length(p11) > 1  )   stop("Illegal arguments. Only k can be a vector.\n" )
	 
	isnt.anint = function(x, tol = .Machine$double.eps^0.5) !(abs(x - round(x)) < tol)
	
	if( 	any(isnt.anint(k)) ||
		isnt.anint(n1)     ||	
		isnt.anint(n2)      ) stop("k, n1 and n2 need to be interger-valued.\n"   )	
	
	if( p10<0 || p01<0 || p11<0) stop("Negative probabilities in arguments\n"        )
	if( p10>1 || p01>1 || p11>1) stop("Probabilities > 1\n in arguments.\n"          )
	if( abs(1-p10-p01-p11) < 0 ) stop("Probabilities in arguments sum to > 1.\n"     )

	kmin = ifelse(n1>=n2,-n2,-n1)
	kmax = ifelse(n1>=n2,+n1,+n2)
	if(min(k)<kmin) stop("too small k \n")
	if(max(k)>kmax) stop("too large k \n")

	if(length(k)==1) return(pcbern.single(k,n1,n2,p10,p01,p11,lower.tail))

	if(n1==n2){
		res = sapply(k,function(x) pcbern.single(x,n1,n2,p10,p01,p11,lower.tail,g.approx=g.approx))
	} else {
		
		#- gaussian approximation
		if(g.approx){
			if(cbern_asym_approx(n1,n2,p10,p01,p11)){
				pars = cbern_asym_approx(n1,n2,p10,p01,p11,pars=TRUE)
				if(lower.tail) {
					pvf <- function(k) pnorm(k+0.5,mean=pars[1],sd=pars[2])
					res = sapply(k,pvf)
					return(res)
				}
				pvf <- function(k) pnorm(k-0.5,mean=pars[1],sd=pars[2],lower.tail=FALSE)
				res = sapply(k,pvf)
				return(res)
			}
		}

		#- use densities
		tmp = dcbern.dl(n1,n2,p10,p01,p11,g.approx=g.approx)
		kmin = tmp[1,1]
		kmax = tmp[nrow(tmp),1]
		fu <- function(k){
			if(lower.tail){
				if(k<kmin) return(0)
				return(sum(tmp[1:(k-kmin+1),2]))
			} else {
				if(k>kmax) return(0)
				return(sum(tmp[(k:kmax)-kmin+1,2]))
			}
		}
		res = sapply(k,fu)
		return(res)
	}
	
}
	

rcbern = function(n,n1,n2,p10,p01,p11){
#######################################

	nmin = min(n1,n2)
	nmax = max(n1,n2)	

	#- draw from multinomial 
	tmp1  = rmultinom(n,nmin,prob=c(p10,p01,p11,1-p10-p01-p11))
	tmp1  = apply(tmp1,2,function(x) x[1] - x[2])

	#- draw from binomial if necessary
	tmp2 = 0
	if(nmax>nmin){
		if(n1>n2) p = p10 + p11
		if(n2>n1) p = p01 + p11
		tmp2  = rbinom(n,nmax-nmin,p)

		if(n1>n2) return(tmp1+tmp2)
		return(tmp1-tmp2)
	} else {
		return(tmp1)
	}
}



