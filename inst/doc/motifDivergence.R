## ----include=FALSE,echo=FALSE--------------------------------------------
  require(knitr)
  opts_chunk$set(comment=NA, echo=TRUE, tidy=FALSE, eval=TRUE, prompt=TRUE, warning=FALSE)

## ----message=FALSE, warning=FALSE----------------------------------------
require(motifDiverge)
require(Biostrings)
require(MotifDb)
enh.hg.file = system.file( "extdata", "enh_human.fa", 
			   package="motifDiverge"   )
enh.mm.file = system.file( "extdata", "enh_mouse.fa", 
			   package="motifDiverge"   )  
enh.hg = readDNAStringSet(enh.hg.file)
enh.mm = readDNAStringSet(enh.mm.file)
enh.hg
enh.mm

## ------------------------------------------------------------------------
providerId = "MA0063.1" #- Jaspar NKX-2.5
index      = grep(providerId,values(MotifDb)$providerId)
pspm       = MotifDb[index][[1]]
pssm       = pspmToPssm(pspm)
pspm       = pspmToPssm(pspm,return.pspm=TRUE)$pspm

## ------------------------------------------------------------------------
bg  = colSums(alphabetFrequency(c(enh.hg,enh.mm))[,1:4])
bg  = bg/sum(bg)
cut = scoreCutFromErr(err=.01, pssm=pssm, pspm=pspm, bg=bg, type="type1")

## ------------------------------------------------------------------------
pars.nomodel = cbernEstimateModelPars(	seqs.x = enh.mm,
					seqs.y = enh.hg,
					pssm   = pssm  ,
					pspm   = pspm  ,
					cut.fw = cut   ,
					cut.rc = cut   )  

## ----echo=FALSE----------------------------------------------------------
evo.mod.file = system.file( "extdata", "primates.mod", package="motifDiverge" )

## ------------------------------------------------------------------------
require(rphast)
neutral.mod = get.neutralMod(evo.mod.file)       #- point ucscURL to model
neutral.mod = mod.backgd.tm(neutral.mod,bg) #- adjust background
pspm.mods   = get.modList(neutral.mod,pspm) #- models for pspm columns
pars.model  = cbernEstimateModelPars(  seqs.x    = enh.mm   ,
                                       seqs.y    = enh.hg   ,
                                       pssm      = pssm     ,
                                       pspm      = pspm     ,
                                       cut.fw    = cut      ,
                                       cut.rc    = cut      ,
				       indep     = FALSE    ,
				       useCounts = FALSE    ,
				       modList   = pspm.mods)

## ------------------------------------------------------------------------
pvals.enr = apply(pars.model[,1:6],1,function(x) pcbern(x[1],x[2],x[3],
							x[4],x[5],x[6],
							lower.tail=F))
pvals.dep = apply(pars.model[,1:6],1,function(x) pcbern(x[1],x[2],x[3],
                                                        x[4],x[5],x[6],
                                                        lower.tail=T))

