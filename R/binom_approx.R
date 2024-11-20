
binom_approx <- function(n,p,pars=FALSE){
#########################################

	res = FALSE
	if(n*p > 10) 					      res = TRUE
	if(abs(1/sqrt(n)*(sqrt((1-p)/p)-sqrt(p/(1-p)))) > .3) res = FALSE
	if(3*sqrt(n*p*(1-p)) > n)		 	      res = FALSE
	if(3*sqrt(n*p*(1-p)) < 0) 			      res = FALSE
	
	#- parameter results
	if(pars) {
		if(res) return(c(n*p,sqrt(n*p*(1-p))))
		return(c(NA,NA))
	}

	#- boolean results
	return(res)
}

cbern_sym_approx <- function(n,p10,p01,pars=FALSE){
###################################################

	res       = FALSE
	if( binom_approx(n,p10)) res = TRUE
	if(!binom_approx(n,p01)) res = FALSE
	
	#- parameter results
	if(pars){
		if(res){
			mn = n*(p10-p01)
			vr = n*p10*(1-p10)+n*p01*(1-p01)-2*n*p10*p01
			return(c(mn,sqrt(vr)))
		}
		return(c(NA,NA))
	}

	#- boolean result
	return(res)

}

cbern_asym_approx <- function(n1,n2,p10,p01,p11,pars=FALSE){
############################################################

	n.min = min(n1,n2)
	n.dif = abs(n1-n2)
	p     = ifelse(n1>n2,p10+p11,p01+p11)
	res   = FALSE
	if(cbern_sym_approx(n.min,p10,p01)) res = TRUE
	if(!binom_approx(n.dif,p)) res = FALSE
	
	#-parameter results
	if(pars){
		if(res){
			prs.1 =	cbern_sym_approx(n.min,p10,p01,pars=TRUE) 
			prs.2 = binom_approx(n.dif,p,pars=TRUE)
			return(c(prs.1[1]+prs.2[1],sqrt(prs.1[2]^2+prs.2[2]^2)))
		}
		return(c(NA,NA))
	}

	#- boolean result
	return(res)
}





