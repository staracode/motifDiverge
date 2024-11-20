
cbernEstimateModelPars <- function(seqs.x,seqs.y,pssm,pspm,cut.fw,cut.rc,indep=TRUE,useCounts=TRUE,modList=NULL,NSAMP=100000,aligned=NULL){
#================================================================================================================

	#- checks go here

	if(useCounts==FALSE){
		if(is.null(modList)) stop("Tree models are needed for model-based estimates.")
		res = estimate.pars.Model(	seqs.x=seqs.x	,
						seqs.y=seqs.y	,
						pssm=pssm	,
						pspm=pspm	,
						cut.fw=cut.fw	,
						cut.rc=cut.rc	,
						indep=indep	,
						mods=modList	,
						NSAMP=NSAMP	,
						aligned=aligned )
		return(res)
	} else {
		res = estimate.pars.noModel(	seqs.x=seqs.x	,
						seqs.y=seqs.y	,
						pssm=pssm	,
						cut.fw=cut.fw	,
						cut.rc=cut.rc	,
						indep=indep	,
						aligned=aligned )
		return(res)
	}
}

estimate.pars.noModel <- function(seqs.x,seqs.y,pssm,cut.fw,cut.rc,indep=TRUE,aligned=NULL){
#===========================================================================================

        l = dim(pssm)[2]

        count.fu    = function(DNAString) {
        #----------------------------------
                a = start(matchPWM(pssm,DNAString,min.score=cut.fw))
                b = start(matchPWM(reverseComplement(pssm),DNAString,min.score=cut.rc))
                return(sort(unique(c(a,b))))
        }

        if(indep){
        #=========

        starts.x = lapply((seqs.x),count.fu) #- FIXME: list conversion slow
        starts.y = lapply((seqs.y),count.fu) #- FIXME: list conversion slow
        hits.x = unlist(lapply(starts.x,length))
        hits.y = unlist(lapply(starts.y,length))

        p.hat = (hits.x + hits.y) / (width(seqs.x)+width(seqs.y)-2*l+2)

        n11.x = unlist(lapply(starts.x,function(x) sum(diff(x)==1)))
        n11.y = unlist(lapply(starts.y,function(x) sum(diff(x)==1)))

        p.hat.1.1 = (n11.x+n11.y)/(hits.y+hits.x)
        p.hat.1.1[hits.y+hits.x==0] = 0 #- ROUGH for short sequences (little hits)

        rhofu <- function(p.hat,lam.hat, kx,ky){
		
		if(lam.hat == p.hat) return(0)
                if(p.hat == 0) return(0)

                afu <- function(p,lam,n){
                        2*p*(1-p)*(lam-p)/(1-lam)*(( (n-1) - (lam-p)/(1-lam)*(1-((lam-p)/(1-p))^n) ))
                }
                #rho = 1/2/min(kx,ky)*( afu(p.hat,lam.hat,kx) + afu(p.hat,lam.hat,ky) )
                rho = (-1) * 1 / 2 / min(kx,ky) / p.hat / (1-p.hat) * ( afu(p.hat,lam.hat,kx) + afu(p.hat,lam.hat,ky) )
                return(rho)
        }

        afu<- function(index) rhofu(p.hat[index],p.hat.1.1[index],width(seqs.x)[index]-l+1, width(seqs.y)[index]-l+1)
        rhos = sapply(1:length(seqs.x),afu)

        #p11 = p.hat^2 - pmin(rhos,0) ; p11 = pmin(p.hat,p11)
        #p10 = p.hat - p11
        #p01 = p10

	rmin = -p.hat/(1-p.hat)
	rmax = 1
	rhos[rhos < rmin]  = rmin[rhos < rmin]
	rhos[rhos > rmax]  = rmax 
	p01 = p.hat*(1-p.hat)*(1-rhos)
	p01 = pmin(p01,p.hat) #- not larger than p.hat
	p01 = pmax(0,p01)     #- not smaller than zero
	p10 = p01
	p11 = p.hat - p10

        res = cbind(hits.x-hits.y,width(seqs.x)-l+1,width(seqs.y)-l+1,p10,p01,p11,p.hat,hits.x,hits.y)
        colnames(res) = c("nxy","kx","ky","p10","p01","p11","p.hat","nx","ny")
        return(res)

        } else {
        #=======
	
	if(is.null(aligned)) stop("need aligend sequences for useCounts=TRUE and indep=FALSE")
	seqs.x = aligned$seqs.x #- aligned version
	seqs.y = aligned$seqs.y #- aligned version

        starts.x = lapply((seqs.x),count.fu) #- FIXME: list conversion slow
        starts.y = lapply((seqs.y),count.fu) #- FIXME: list conversion slow
        hits.x = unlist(lapply(starts.x,length))
        hits.y = unlist(lapply(starts.y,length))

        p.hat = (hits.x + hits.y) / (width(seqs.x)+width(seqs.y)-2*l+2)

        p.hat.1.1 = sapply(1:length(seqs.x),function(index) sum(starts.x[index] %in% starts.y[index])/hits.x[index])
        p.hat.1.1[hits.x==0] =0
        p.hat.1.1[p.hat.1.1 < p.hat] = (p.hat)[p.hat.1.1 < p.hat]

        p11 = p.hat * p.hat.1.1
        p10 = p.hat - p11
        p01 = p10

        res = cbind(hits.x-hits.y,width(seqs.x)-l+1,width(seqs.y)-l+1,p10,p01,p11,p.hat,hits.x,hits.y)
        colnames(res) = c("nxy","kx","ky","p10","p01","p11","p.hat","nx","ny")
        return(res)

        }

}

estimate.pars.Model <- function(seqs.x,seqs.y,pssm,pspm,cut.fw,cut.rc,indep=FALSE,mods,hitProbs=NULL,NSAMP=100000,aligned=NULL){
#==================================================================================================================

        l = dim(pssm)[2]

	#- prepare sequences
	#===================
	if(is.null(aligned)){
		#- no alignments; need to do them
		message("aligning sequences; keeping stand")
		nseqs      = length(seqs.x)
		
		afu <- function(ind){
			message("seqence ",ind, " of ", nseqs)
			tmp = pairwiseAlignment(seqs.x[[ind]],seqs.y[[ind]])
			return(list(sx = pattern(tmp),sy=subject(tmp)))
		}
	
		tmp = lapply(1:nseqs,afu)
		tmp.x = lapply(tmp,function(x) x$sx)
		tmp.y = lapply(tmp,function(x) x$sy)
		seqs.x.ali = unlist(DNAStringSetList(tmp.x))
		seqs.y.ali = unlist(DNAStringSetList(tmp.y))
	} else {
		#- take passed-on alignments
		seqs.x.ali = aligned$seqs.x
		seqs.y.ali = aligned$seqs.y
	}

        #- 1. get estimate for p.hat
        #----------------------------

        count.fu <- function(seq) {
                a = start(matchPWM(pssm,seq,min.score=cut.fw))
                b = start(matchPWM(reverseComplement(pssm),seq,min.score=cut.rc))
                return(sort(unique(c(a,b))))
        }

        starts.x = lapply(as.list(seqs.x),count.fu) #- FIXME: list conversion slow
        starts.y = lapply(as.list(seqs.y),count.fu) #- FIXME: list conversion slow

        hits.x = unlist(lapply(starts.x,length))
        hits.y = unlist(lapply(starts.y,length))

        p.hat = (hits.x + hits.y) / (width(seqs.x)+width(seqs.y)-2*l+2)

        #- 2. get estimate for zeta in the  HMM
        #--------------------------------------

	bg        = alphabetFrequency(seqs.x) + alphabetFrequency(seqs.y)
	bg        = colSums(bg[,c("A","C","G","T")]) ; bg = bg/sum(bg)
	bg        = snapVec(.5*(bg + reverseComplement(t(bg)))) ; bg = as.vector(bg)
	ePspm     = get.extendedPspm(pspm, bg)
        allMotifs = get.allMotifs(pspm,ePspm)

        if(is.null(hitProbs)){
        hitProbs = apply(allMotifs,1,function(x){
                dist.fw = scoreDist(pssm,ePspm[,x])
                dist.rc = scoreDist(reverseComplement(pssm),ePspm[,x])
                sum(dist.fw[dist.fw[,1]>=cut.fw,2]) + sum(dist.rc[dist.rc[,1]>=cut.rc,2])
                })
        }

        get.p = function(zeta,allMotifs,hitProbs){
        #-----------------------------------------
        
		eq.freq.bg =  (1-zeta)/((1-zeta)+l*zeta)
        	eq.freq.mt =  zeta/2/((1-zeta)+l*zeta)
        	eq.freqs   =  c(eq.freq.bg,rep(eq.freq.mt,2*l))
        	trans.mat  =  get.hmm(pspm,zeta/2,zeta/2,transmat.only=TRUE)
        	trans.mat  =  get.hmm(pspm,zeta/2,zeta/2,transmat.only=TRUE)
        	p1         = sapply(allMotifs[,1],function(x) eq.freqs[x])

        	i.beg = allMotifs[,1:(l-1)]
        	i.end = allMotifs[,2:l]
        	lp2   = sapply(1:nrow(i.beg) ,function(x)  sum(log(trans.mat[cbind(i.beg[x,],i.end[x,])])) )

        return(sum(exp(log(p1)+lp2+log(hitProbs))))

        }

 	#- what if zeta is zero?:q
	p.min   = get.p(0,allMotifs,hitProbs)
	min.ind = p.hat <= p.min 	
        zetas = sapply(p.hat,function(x){
		message("x",appendLF=FALSE)
                ofu = function(zeta){
                        abs(get.p(zeta,allMotifs,hitProbs) - x)
                }
                zeta.hat = optimize(ofu,interval=c(1E-5,.25))$minimum
	        return(zeta.hat)
        })
	zetas[min.ind] = 1E-5 #- otherwise rphast unstable FIXME: easier sampling possible

	#- 3. get the parameter estimates by simulation
	#----------------------------------------------

        afu <- function(index){

		#- FIXME: afu should be modularized

                if(index%%10==0) message(".",appendLF=FALSE)
                zeta       = zetas[index]
                ##eq.freq.bg = (1-zeta)/((1-zeta)+l*zeta)
                ##eq.freq.mt =  zeta/2/((1-zeta)+l*zeta)
                ##eq.freqs   = c(eq.freq.bg,rep(eq.freq.mt,2*l))
                moiHmm     = get.hmm(pspm,zeta/2,zeta/2)

                #- get the time estimate
                minlen = min(width(seqs.x[index]),width(seqs.y[index]))
                sx     = substr(seqs.x.ali[index],start=1,stop=minlen)
                sy     = substr(seqs.y.ali[index],start=1,stop=minlen)
                msa    = msa(c(as.character(sx),as.character(sy)),names=c("sx","sy"))


                get.lik <- function(time){
                        time = exp(time)
                        tmp.mods      = mods
                        for(i in seq_along(tmp.mods)){
                                tmp.mods[[i]]$tree = rescale.tree(      tmp.mods[[i]]$tree,
                                                                        scale=1/branchlength.tree(tmp.mods[[i]]$tree)*time)
                        }
                        tmpres = score.hmm(msa=msa,mod=tmp.mods,hmm=moiHmm,quiet=TRUE)
                        return(-tmpres$likelihood)
                }
                time = exp(optimize(get.lik,interval=c(-5,5))$minimum)

                #- NSAMP sample
                seqs = get.seqs (mods,moiHmm,time,(l+1)*NSAMP,toBiostring=TRUE)
                a  = matchPWM(pwm=pssm,subject= seqs[[1]],min.score=cut.fw)
                b  = matchPWM(pwm=reverseComplement(pssm),subject= seqs[[1]],min.score=cut.rc)
                starts.1 = sort(unique(c(start(a),start(b))))
                a  = matchPWM(pwm=pssm,subject= seqs[[2]],min.score=cut.fw)
                b  = matchPWM(pwm=reverseComplement(pssm),subject= seqs[[2]],min.score=cut.rc)
                starts.2 = sort(unique(c(start(a),start(b))))

                p.1.1.ind = sum(diff(starts.1) == 1)/length(starts.1)
                p.1.1.cor = sum(starts.1 %in% starts.2)/length(starts.1)

                res = (c(p.1.1.ind,p.1.1.cor,time))
                return(res)
        }

        pars = sapply(1:length(zetas), afu)
        pars = t(pars)

        p11.dep  = p.hat * pars[,2]
        p10.dep  = p.hat - p11.dep
        p01.dep  = p10.dep
	times    = pars[,3] 

	if(indep==FALSE){
        	res = cbind(hits.x-hits.y,width(seqs.x)-l+1,width(seqs.y)-l+1,p10.dep,p01.dep,p11.dep,hits.x,hits.y,p.hat,zetas,times)
        	colnames(res) = c("nxy","kx","ky","p10","p01","p11","nx","ny","p.hat","zeta.hat","that")
		return(res)
	} else if(indep==TRUE){

        	rhofu <- function(p.hat,lam.hat, kx,ky){

                	afu <- function(p,lam,n){
                    	    2*p*(1-p)*(lam-p)/(1-lam)*(( (n-1)*-(lam-p)/(1-lam)*(1-((lam-p)/(1-p))^n) ))
                	}
                	rho = 1/2/min(kx,ky)*( afu(p.hat,lam.hat,kx) + afu(p.hat,lam.hat,ky) )
                	return(rho)
       		}

        	afu<- function(index) rhofu(p.hat[index],pars[index,1],width(seqs.x)[index]-l+1, width(seqs.y)[index]-l+1)
        	rhos = sapply(1:length(seqs.x),afu)

        	p11 = p.hat^2 - pmin(rhos,0) ; p11 = pmin(p.hat,p11)
        	p10 = p.hat - p11
        	p01 = p10

        	res = cbind(hits.x-hits.y,width(seqs.x)-l+1,width(seqs.y)-l+1,p10,p01,p11,hits.x,hits.y,p.hat,zetas,times)
        	colnames(res) = c("nxy","kx","ky","p10","p01","p11","nx","ny","p.hat","zeta.hat","that")
        	return(res)
	} else {

		return(NA)
	}

}

get.extendedPspm = function(pspm,bg=NULL){
#=========================================

        if(is.null(bg)) bg = rowSums(pspm)/sum(pspm)
        e.PSPM = cbind(bg,pspm,reverseComplement(pspm))
}

get.allMotifs = function(pspm,ePspm=NULL,bg=NULL){
#=================================================

        if(is.null(ePspm)) ePspm = get.extendedPspm(bg,pspm)

        l = ncol(pspm)
        i.bg = 1
        i.fw = 2
        i.rc = l+2

        #- pure
        A1 = rep(i.bg,l)
        A2 = (i.fw:(i.fw+l-1))
        A3 = (i.rc:(i.rc+l-1))

        #- 1-mix
        A4 = t(sapply(0:(l-2),function(x) c(rep(i.bg,(l-1)-x),i.fw:(i.fw+x))))
        A5 = t(sapply(0:(l-2),function(x) c(rev((i.fw+l-1):(i.fw+l-1-x)),rep(i.bg,(l-1)-x))))
        A6 = t(sapply(0:(l-2),function(x) c(rep(i.bg,(l-1)-x),i.rc:(i.rc+x))))
	A7 = t(sapply(0:(l-2),function(x) c(rev((i.rc+l-1):(i.rc+l-1-x)),rep(i.bg,(l-1)-x))))

        #- 2-mix
        n.bg = 1:(l-2)
        n.rc = 1:(l-2)
        f.a <- function(n.bg,n.rc) c((i.fw+n.bg+n.rc):(i.fw+l-1),rep(i.bg,n.bg),i.rc:(i.rc+(n.rc-1)))
        f.b <- function(n.bg,n.rc) c((i.rc+n.bg+n.rc):(i.rc+l-1),rep(i.bg,n.bg),i.fw:(i.fw+(n.rc-1)))

        afu <- function(n.bg) {
                n.rc = (l-1-n.bg):1
                a = t(sapply(n.rc,function(x) f.a(n.bg,x)))
                b = t(sapply(n.rc,function(x) f.b(n.bg,x)))
                rbind(a,b)
        }

        A8 = sapply(n.bg,afu)
        A8 = matrix(unlist(lapply(A8,t)),ncol=l,byrow=T)

        A = rbind(A1,A2,A3,A4,A5,A6,A7,A8)
        rownames(A) = paste("M",1:nrow(A),sep="")
        return(A)
}

get.neutralMod <- function(filename, download=FALSE){
#====================================================

        #- get the model
	if(download){
        	tfile = tempfile()
        	download.file(filename, tfile)
        	neutral.mod = read.tm(tfile)
        	unlink(tfile)
	} else {
		neutral.mod = read.tm(filename)
	}
	
        #- only two sequences and combined length of 1 expected subs
        names            = leafnames.tree(neutral.mod$tree)[1:2]
        neutral.mod$tree = prune.tree(neutral.mod$tree,seqs=names,all.but=TRUE)
        neutral.mod$tree = rescale.tree(neutral.mod$tree,scale=1/branchlength.tree(neutral.mod$tree))
        neutral.mod$tree = rename.tree(neutral.mod$tree,old.names=names,new.names=c("sx","sy"))

        return(neutral.mod)
}

get.modList <- function(neutral.mod,pspm,names=NULL){
#====================================================

        #- make the models for each motif position
	l = ncol(pspm)
        motif.mod.fw = vector(length=l,mode="list")
        for(i in seq_along(motif.mod.fw)){
                motif.mod.fw[[i]] = mod.backgd.tm(neutral.mod,new.backgd = pspm[,i])
        }

        motif.mod.rc = vector(length=l,mode="list")
        for(i in seq_along(motif.mod.rc)){
                motif.mod.rc[[i]] = mod.backgd.tm(neutral.mod,new.backgd = pspm[4:1,l-i+1])
        }

        #- one list for background and motif
        mods                                  = vector(length=2*l+1,mode="list")
        mods[[1]]                             = neutral.mod
        mods[2:(ncol(pspm)+1)]                = motif.mod.fw
        mods[(ncol(pspm)+2):(2*ncol(pspm)+1)] = motif.mod.rc

	if(is.null(names)) names = c("neutral", paste("site.",1:l,".fw",sep=""), paste("site.",1:l,".rc",sep=""))
        names(mods)                           = names

        return(mods)
}

get.hmm <- function(pspm,zeta.fw,zeta.rc,transmat.only=FALSE,...){
#=================================================================

        l = ncol(pspm)
        trans.mat = matrix(0,ncol=2*l+1,nrow=2*l+1)
        colnames(trans.mat) = c("neutral", paste("site.",1:l,".fw",sep=""), paste("site.",1:l,".rc",sep=""))
        rownames(trans.mat) = colnames(trans.mat)

        trans.mat["neutral","neutral"] = 1-zeta.fw - zeta.rc
        trans.mat["neutral","site.1.fw"] = zeta.fw
        trans.mat["neutral","site.1.rc"] = zeta.rc
        trans.mat[cbind(2:(ncol(pspm)),3:(ncol(pspm)+1))] = 1
        trans.mat[(ncol(pspm)+1),"neutral"] = 1 - zeta.fw - zeta.rc
        trans.mat[(ncol(pspm)+1), 2] = zeta.fw
        trans.mat[(ncol(pspm)+1),(ncol(pspm)+2)] = zeta.rc
        trans.mat[cbind((ncol(pspm)+2):(2*ncol(pspm)), (ncol(pspm)+3):(2*ncol(pspm)+1))] = 1
        trans.mat[2*ncol(pspm)+1,"neutral"] = 1 - zeta.fw - zeta.rc
        trans.mat[2*ncol(pspm)+1,2] = zeta.fw
        trans.mat[2*ncol(pspm)+1,ncol(pspm)+2] = zeta.rc

        if(transmat.only) return(trans.mat)

        #- rphast don't do this right
	zeta 	   = zeta.fw + zeta.rc
        eq.freq.bg = (1-zeta)/((1-zeta)+l*zeta)
        eq.freq.mt =  zeta/2/((1-zeta)+l*zeta)
        eq.freqs   = c(eq.freq.bg,rep(eq.freq.mt,2*l))

        pspm.hmm = hmm(trans.mat,eq.freq=eq.freqs,begin.freq=eq.freqs,...)

        return(pspm.hmm)
}

get.seqs <- function(mods,hmm,t,len,toBiostring=FALSE){
#======================================================

        for(i in seq_along(mods)){
                #- makes copy
                mods[[i]]$tree = rescale.tree(mods[[i]]$tree,scale=1/branchlength.tree(mods[[i]]$tree)*t)
        }
        seqs = simulate.msa(mods, len, hmm=hmm, get.features=TRUE)
        if(!toBiostring) return(seqs)
        res = DNAStringSet(seqs[[1]][[1]])
        names(res) = c("sx","sy")
        return(res)

}

