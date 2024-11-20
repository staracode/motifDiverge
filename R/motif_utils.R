
scoreDist <- function(pssm,pspm=rep(1/4,4),gran=0.01){
######################################################

	pssm = round(pssm/gran)
	pssm = pssm*gran

	if(is.null(ncol(pspm))) pspm = matrix(rep(pspm,ncol(pssm)),nrow=4,ncol=ncol(pssm)) 

	tmin = min(pssm[,1])
	tmax = max(pssm[,1])
	tvec = rep(0,(tmax-tmin)/gran+1)
	tvec[pssm[,1]/gran+1-tmin/gran] = pspm[,1]
	for(i in 2:ncol(pssm)){
		nmin = min(pssm[,i])
		nmax = max(pssm[,i])
		nvec = rep(0,(nmax-nmin)/gran+1)
		nvec[pssm[,i]/gran+1-nmin/gran] = pspm[,i]
                tvec = rcpp_convolve(tvec,nvec)
		tmin = tmin + nmin
		tmax = tmax + nmax
        }

	tmp = tvec
	tmp = zapsmall(tmp)
        tmp = tmp/sum(tmp)
        prb = tmp[tmp>0]

	vls = (which(tmp>0) -1) * gran + tmin
	return(cbind(vls,prb))
}


scoreCutFromErr <- function(err,pssm,pspm,bg=rep(1/4,4),type="type1",gran=0.01){
################################################################################

	if(type=="type1"){
		bgdist = scoreDist(pssm,pspm=bg,gran=gran)
		ind    = min(which(1-cumsum(bgdist[,2]) <= err))
		return(bgdist[ind,1])
	} else if(type=="type2") {
		fgdist = scoreDist(pssm,pspm=pspm,gran=gran)
                ind    = max(which(cumsum(fgdist[,2])  <= err))
		return(fgdist[ind,1])
	} else {
		bgdist = scoreDist(pssm,pspm=bg  , gran=gran)
		fgdist = scoreDist(pssm,pspm=pspm, gran=gran)
		cdf.bg = stepfun(bgdist[,1],c(0,cumsum(bgdist[,2])))
		cdf.fg = stepfun(fgdist[,1],c(0,cumsum(fgdist[,2])))
		minval = min(bgdist[1,1],fgdist[1,1])	
		maxval = max(max(bgdist[,1]),max(fgdist[,1]))
		t1err  = 1-cdf.bg(seq(minval,maxval,by=gran))
		t2err  = cdf.fg(seq(minval,maxval,by=gran))
		ind    = which.min(abs(t1err-t2err))
		return(seq(minval,maxval,by=gran)[ind])
	}
}


errFromScoreCut <- function(cut,pssm,pspm,bg=rep(1/4,4),type="type1",gran=0.01){
################################################################################

	if(type=="type1"){
		bgdist = scoreDist(pssm,pspm=bg,gran=gran)
		return(sum(bgdist[bgdist[,1]>=cut,2]))
	} else if(type=="type2"){
		fgdist = scoreDist(pssm,pspm=pspm,gran=gran)
		return(sum(fgdist[fgdist[,1]<=cut,2]))
	} else {
		stop("Type argument makes no sense.\n")
	}
}



pscmToPspm <- function(pscm, bg=NULL, digits=4, reg=TRUE){
#=========================================================

	#- get background from motif count matrix
	if(is.null(bg)) bg = rowSums(pscm)/sum(pscm) 

	#- add a (as in: one) pseudocount distributed accroding to the background
	if(reg) pscm = pscm + bg

	#- snap to digits and divide (pscm can still contain zeros, also when reg)
	#---------------------------

	DIGITS = digits
        pscm.t = t(t(pscm) / rowSums(t(pscm)))
        min = 1/10^DIGITS
        pscm.t[pscm.t <= min] = min + 1* 1/10^(DIGITS+1)
        pscm.t = round(pscm.t,DIGITS)

        pscm.t = pscm.t * 10^DIGITS

        while(sum(colSums(pscm.t) - .4 > 10^DIGITS) > 0){
                ind = which(colSums(pscm.t) -.4 > 10^DIGITS)
                for(i in ind){
                        pscm.t[which.max(pscm.t[,i]),i] =  pscm.t[which.max(pscm.t[,i]),i] - 1
                }
        }

        while(sum(colSums(pscm.t) +.4 < 10^DIGITS) > 0){
                ind = which(colSums(pscm.t) + .4 < 10^DIGITS)
                for(i in ind){
                        pscm.t[which.min(pscm.t[,i]),i] =  pscm.t[which.min(pscm.t[,i]),i] + 1
                }
        }
        pscm.t = round(pscm.t/10^DIGITS,DIGITS)

        return(pscm.t)
}


pspmToPssm <- function(pspm, bg=NULL, digits=4, return.pspm=FALSE){
#==================================================================

	#- regularize pspm if necessary (should be avoided)
	if(any(pspm< 10^(-digits))){

		#- this is a mess but (hopefully) a good guess
		pscm = round( pspm * 1/min(pspm[pspm>10^(-digits)]) )
		pspm = pscmToPspm(pscm,bg=bg,digits=digits)
	}

	if(is.null(bg)) bg = rowSums(pspm)/sum(pspm)

	pssm = log2(pspm) - log2(bg)

	if(return.pspm == TRUE) return(list(pssm=pssm,pspm=pspm))
	return(pssm)
}


pscmToPssm <- function(pscm, bg=NULL, digits=4, reg=TRUE){
#=========================================================

	pspm = pscmToPspm(pscm,bg=bg,digits=digits,reg=reg)
	pssm = pspmToPssm(pspm,bg=bg)

	return(pssm)
}


makeRandomPspm = function(len=8){
#================================

	#- colums sampled from unit simplex 
	pspm = matrix(-log(runif(len*4)),nrow=4)
	pspm = t(t(pspm)/colSums(pspm))
	return(pspm)

}



snapVec <- function(vec, vecSum=1, DIGITS=4, pC=FALSE){
#======================================================

	if(any(vec<0)) stop("snapVec does not cope with negative entries\n")

	vec                 = vec/sum(vec)*vecSum
        vec                 = round(vec, DIGITS)
        if(pC) vec[vec ==0] = 1/10^DIGITS

        vec            = vec * 10^DIGITS
        while(sum(vec) - .4 > (10^DIGITS)*vecSum){
                vec[which.max(vec)] =  vec[which.max(vec)] - 1
        }
        while(sum(vec) + .4 < (10^DIGITS)*vecSum){
                vec[which.min(vec)] =  vec[which.min(vec)] + 1
        }
        vec = round(vec/10^DIGITS, DIGITS)

        return(vec)
}

