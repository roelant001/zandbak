### function to compare RPMs
### IN: 
###  -RPM A as logical matrix: voxels as rows, cols as pts
###  -RPM B as logical matrix: voxels as rows, cols as pts
### OUT:
###  -sA.nii: summed locations of rpm A
###  -sB.nii: summed locations of rpm B
###  -sAB.nii: summed locations of rpm A and B combined
###  -diff.nii: difference between proportion of tumors in deficit and non-deficit maps
###  -p.nii: randomized p-value map (corrects for auto-correlation)
###  -qnull.nii: q-values based on null hypothesis distribution of randomizations 
###  -qBH.nii: q-values based on Benjamini Hochberg method 
### additional standard parameters for testing:
 
debug <- FALSE 
if (debug) {
      min.pt.info          <- 1
      nr.randomizations    <- 200
      results.dir          <- '/Users/philip/Downloads/'
      #results.dir          <- '/home/ubuntu/results/'
      #mni.file             <- '/Users/philip/Documents/projects/probabilisticmap/atlases/MNI152/sMNI152_T1_1mm.nii.gz'
      mni.file             <- '/Users/philip/Documents/projects/probabilisticmap/atlases/MNI152/MNI152_T1_1mm_brain.nii.gz'
      #mni.file             <- '/home/ubuntu/data/MNI152_T1_1mm_brain.nii.gz'
      false.discovery.rate <- 0.40
      my.result.txt        <- 'test_'
      my.dimensions        <- dim( example.2mm )
      mylogfile              <- "log.txt"
      #min.pts=min.pt.info; R=nr.randomizations; rdir=results.dir; mfile=mni.file; fdr=false.discovery.rate; mytxt=my.result.txt; mydim=my.dimensions
}

### version history

### v7-6:
###  - randomization p values based on relative risk row (equal to cols)
###  - p based on global null distribution (code for local retained)
###  - q0 using global distribution for Vt (code for local retained)
### v7-5:
###  - randomization p values based on relative risk row
###  - p based on local null distribution
###  - storey's q
### v7-4:
###  - randomization p values based on relative risk row + col
###  - added continuity correction (+1/2) to avoid infinity
### v7-3:
###  - Barnard extact p values instead of Fisher exact p values
### v7-2:
###  - no randomization p value, use simple Fisher extact p values
###  - no q0, use Storey qvalue
### v7-1:
###  - p values calculated by local (per voxel) randomizations instead of global
### v7-0:
###  - replaced diffcol and diffrow by percol and perrow
###  - remove tdiffcol and tdiffrow
###  - write nii after resampling to MNI space
### v5-8:
###  - fixed mistake in plotting of nA and nB, tB and ntB
###  - combined diffcol and diffrow in one run
###  - q according to Storey and null distribution
### v5-6:
###  - test statistic is difference in rows (comparing population w/o and with deficit)
###  - p-value based on global diff null distribution, not on null distribution per voxel
###  - q-value according to Storey
### v5-5:
###  problem: q-values of lowest p-values become 1.
###  - q-values based on p-values off null distribution instead of diff-values: added R.p.diff
### v5-4:
###  - to compare qnull's: qnullg1, qnullg2, qnulll1, qnulll2 
###  - removed Fp, logOR
### v5-3:
###  - solved mistake in  qnull
###  - show variations in qnull of diff (global null, voxel null, 2-sided, 1-sided)
###  - re-inserted qBH as v4-3 to v5-2
### v5-2:
###  - re-inserted Fp, logOR as test statistic from v4-2
### v5-1:
###  - solved mistake in randomization calculation which introduced a spatial artefact
### v5-0:
###  - only diff1 and qnull and qBH
###  - parallel randomization
### v4-3:
###  - calculate q variants
### v4-2:
###  - all test statistics (Fisher's p, logOR, diff1, diff2, chisquare)
### v4-1:
###  - added logging function of output
### v4-0:
###  - for use local
###  - changed movie part to plot reduced dimensions to MNI@1x1x1mm
### v3-1: 
###  - for use with ubuntu
###  - only fisher test for p
###  - only null distr  for q 

#install.packages('oro.nifti')
#install.packages('data.table')
#install.packages('pbapply')
#biocLite('qvalue')
 library(oro.nifti)
 library(data.table)
 library(pbapply)
 library(qvalue)
 library(parallel)
 library(Exact)
 
compare.RPM.bin <- function(
 A,
 B, 
 min.pts = min.pt.info, 
 R       = nr.randomizations, 
 rdir    = results.dir, 
 mfile   = mni.file, 
 fdr     = false.discovery.rate, #
 mytxt   = my.result.txt, 
 mydim   = my.dimensions,
 mylog   = mylogfile               # opened connection to file in apend mode
 ) {

# message outputs to console; cat to logfile
message("====================")
message(paste('processing ',mytxt,' ...',sep=""))
message(paste('FDR : ',fdr,sep=""))
cat('==================== \n',file=mylog)
cat(date(),'\n',file=mylog)
cat('processing ...       : ',mytxt,'\n',file=mylog)
cat('min pts in a voxel   : ',min.pts,'\n',file=mylog)
cat('nr of randomizations : ',R,'\n',file=mylog)
cat('results directory    : ',rdir,'\n',file=mylog)
cat('false discovery rate : ',fdr,'\n',file=mylog)
cat('dimensions           : ',mydim,'\n',file=mylog)

########################################################################
### RANDOMIZATION TEST OF DIFFERENCE ON TUMOR DISTRIBUTION VOXELWISE ###
########################################################################

# data array AB with information (TRUE/FALSE), pts in cols, vox in rows
 # nrows is all mni voxels
 AB <- data.table( matrix(as.logical(A),ncol=ncol(A)), matrix(as.logical(B),ncol=ncol(B)) )
 AB <- data.table(1:nrow(AB), AB) 
 mycolnames <- c("mni.id",colnames(A),colnames(B))
 setnames(AB,1:ncol(AB),mycolnames)

# get subset of voxels in DT with any tumor information as population of voxels for statistics, alternative is to use mnibrain voxels, but less power
 # nrows is nr of voxels containing patient info
 # alternatively use 'brain' instead of ss
 nrs <- AB[,rowSums(.SD),.SDcols=2:ncol(AB)]
 #hist(nrs)
 cat('nr of voxels with min pts nr \n', file=mylog)
 (tab.nrs <- table(nrs))
 cat(paste(tab.nrs,'\n',sep=""), file=mylog)
 ss <- nrs >= min.pts # info from more than min.pts patients 
 DT <- data.table(1:sum(ss),AB[ss,])
 setnames(DT,1,"DT.id")
 message(paste(nrow(DT),' of ',nrow(AB),' voxels contain information of at least ',min.pts,' patient(s)',sep=''))
 message(paste('A contains ',ncol(A),' patients',sep=''))
 message(paste('B contains ',ncol(B),' patients',sep=''))
 cat(nrow(DT),' of ',nrow(AB),' voxels contain information of at least ',min.pts,' patient(s) \n', file=mylog)
 cat('A contains ',ncol(A),' patients \n', file=mylog)
 cat('B contains ',ncol(B),' patients \n', file=mylog)

# original observations using fast matrix operations with keyed subset of DT
 # note t() to get pts in rows and unique voxels in cols for fast matrix calculations
 # DT is sorted by pat data columns, order is similar to uq.id
 # T/F * 1 is easy trick to get numeric  
 # logical matrix of voxels in columns with pts in rows
 # order of matrix is as 3:ncol(DT) which we need later on
 last.pt.col <- ncol(DT)
 setkeyv(DT,colnames(DT)[3:last.pt.col])
 DT <- data.table(DT,1:nrow(DT))
 setnames(DT,ncol(DT),"M.id")
 M <- as.matrix( t(subset(DT,,key(DT))))
 
 # original group labeling of pts
 grp.org <- c(rep(TRUE,ncol(A)),rep(FALSE,ncol(B)))

 message('started calculating DT ...')

 # calculate Ax, Ay, Bx, By at each voxel
 DT[ ,Ax := as.vector( ( grp.org*1) %*% ( M*1) ) ] # nr of pts with tumor in rpm A
 DT[ ,Bx := as.vector( (!grp.org*1) %*% ( M*1) ) ] # nr of pts with tumor in rpm B
 DT[ ,Ay := as.vector( ( grp.org*1) %*% (!M*1) ) ] # nr of pts w/o  tumor in rpm A
 DT[ ,By := as.vector( (!grp.org*1) %*% (!M*1) ) ] # nr of pts w/o  tumor in rpm B

# An <- ncol(A)
# Bn <- ncol(B)
 
 #  atlas            
 #        +----+----+                                  +----+----+
 #   A    | Ax | Ay | An : pts in A (e.g. no-deficits) |  a |  b |
 #        +----+----+                                  +----+----+
 #   B    | Bx | By | Bn : pts in B (e.g. deficits)    |  c |  d |
 #        +----+----+                                  +----+----+
 #          t    !t   n pts
 #
 # A = no deficit  ;  B = deficit
 # t = with tumor  ; !t = w/o tumor
 # Ax = a ; Ay = b ; Bx = c ; By = d

 # get relative risk for tumor rate among deficits divided by tumor rate among non-deficits
 get.RRrow <- function(a,b,c,d,CC=TRUE) {
 if (CC) { # continuity correction
	 	a <- a + 0.5
	  b <- b + 0.5
	  c <- c + 0.5
	  d <- d + 0.5	}  
  res <- (c / (c + d)) / (a / (a + b))
  return(res)
 }

 # select unique combinations of Ax and Bx for efficient calculation
 setkeyv(DT,c("Ax","Ay","Bx","By"))
 DT.uq <- unique(DT)
 DT.uq <- data.table(1:nrow(DT.uq), DT.uq)
 setnames(DT.uq,1,"uq.id")
 DT.uq[ ,RRr  := get.RRrow(Ax,Ay,Bx,By), by=key(DT.uq)]
 # lookup RR for DT in DT.uq
 DT[ ,RRr  := DT[DT.uq,list(DT.id,RRr)][,RRr] ]

 message(paste('started calculating R.DT with ',R,' randomizations...',sep=''))
 message(paste('using ',R*nrow(DT)/10^6,' million entries',sep=""))

# create R randomizations of group label, first and matrix, then as list
 grp.rnd <- t(replicate(R, sample(grp.org)))

# function to calculate null distribution of randomization (in rows) per voxel (in cols)
## RR rows
 randomizeR.row <- function(myrow,M) {
  # store one randomization in data.table temporary data table tDT to calculate efficiently 
  # tDT in same order as M which is based on DT 3:last.pt.col
  tDT <- data.table(1:ncol(M),
   as.vector(( myrow*1) %*% ( M*1)), # R.Ax
   as.vector(( myrow*1) %*% (!M*1)), # R.Ay
   as.vector((!myrow*1) %*% ( M*1)), # R.Bx
   as.vector((!myrow*1) %*% (!M*1))) # R.By
  setnames(tDT,1:5,c("M.id","R.Ax","R.Ay","R.Bx","R.By"))
  # order according to R.Ax and R.Bx for unique and efficient calculation
  setkeyv(tDT,c("R.Ax","R.Ay","R.Bx","R.By"))
  # identify unique combinations of Ax Bx in tDT.uq to calculate these only once
  tDT.uq <- unique(tDT)
  tDT.uq[ ,RRr    := get.RRrow(R.Ax,R.Ay,R.Bx,R.By), by=key(tDT.uq)]
  # lookup calculations in tDT
  # result lookup first in order of R.Ax and R.Bx, original order of M and DT 3:last.pt.col is stored as tDT.id
  result <- tDT[tDT.uq,list(M.id,RRr)]
  # result second in order tDT.id as of M and DT 3:last.pt.row
  setkeyv(result,c("M.id"))
  # write results to storage matrix by row
  return( result[,RRr] )	
 }
# R.perrow with vox in rows and R in cols in order of M i.e. DT 3:last.pt.row
 startpb(max=R)
 R.RRr <- pbapply(grp.rnd,MARGIN=1,FUN=randomizeR.row,M=M)

# message('started verification null distribution') # test once
#
#  # verify null distribution p-values
#  # p's based on local null (per voxel)
#  # input v: vector (i.e. a row of R.RRr)
#  startpb(min=1,max=nrow(R.RRr))
#  pl.R.RRr <- pbapply(R.RRr, 1, function(v) {sapply(v, function(x) { sum( v >= x ) / length(v) } )} )
#  dim(pl.R.RRr)
#	pl.R.RRr <- t(pl.R.RRr)
#	hist(pl.R.RRr)
#  # verify null distribution p-values
#  # p's based on global null (all voxels)
#  # input e: element (i.e. of R.RRr) and m: matrix
#  pg.R.RRr <- matrix(NA,ncol=ncol(R.RRr),nrow=nrow(R.RRr))
#	startpb(min=1,max=prod(dim(R.RRr)))
#  pg.R.RRr[] <- pbsapply(R.RRr, function(e,m=R.RRr) {sum( m >= e )} )
#  pg.R.RRr <- pg.R.RRr / prod(dim(R.RRr))
#  hist(pg.R.RRr)

 message('started calculating randomized p-values ...')

 # calculate randomization p-value of RRr using global null distribution R.RRr and put in DT
 # testing extremer values using global null
 setkeyv(DT,"RRr")
 DT.uq <- unique(DT)
 denom <- prod(dim(R.RRr))
 DT.uq[ , p.RRr := sum( R.RRr >= RRr ) / denom, by=key(DT.uq)]
 DT[    , p.RRr := DT[DT.uq,list(DT.id,p.RRr)][,p.RRr] ]

# # calculate randomization p-value of RRr using local (per voxel) null distribution R.RRr and put in DT setkeyv(DT,c("M.id"))
# # testing extremer using local null (per voxel)
# setkeyv(DT,"M.id")
# DT[ ,p.RRr := ( sum( R.RRr[M.id,] >= RRr ) / R ), by=key(DT)]   

 message('started calculating randomized q-values ...')

# ## calculate q-value of perrow and percol using Storey and put in DT
# setkeyv(DT,"DT.id")
# #qobj <- qvalue(DT[,p.RRr],pi0.method="bootstrap",robust=TRUE)
# qobj <- qvalue(DT[,p.RRr])
# message(paste("q storey p.RRr pi0 : ",round(qobj$pi0,3),sep=''))
# cat("q storey p.RRr pi0 : ",qobj$pi0,'\n', file=mylog)
# DT[ ,qs.RRr := qobj$qvalues] 

 ## qnull RRr
 ## calculate q-value of perroq using null distribution R.RRr and put in DT
 ## Vt and Rt as table 1 in e.g. Chen, Hum Brain Mapp 2009
  ## Rt: nr of declared significant elements for each randomized p-value in DT
  setkeyv(DT,c("p.RRr")) # because p.RRr's are not unique for each c("Ax","Bx")
  DT.uq <- unique(DT) # rewrite DT.uq
  # calculate weight 
  DT.uq[ ,w  := as.vector( table( DT[ ,p.RRr] ))]
  # calculate number of significant voxels in DT.uq and lookup in DT
  DT.uq[ ,Rt := DT.uq[,sum( DT.uq[,w] * (DT.uq[,p.RRr] <= p.RRr) ), by=key(DT.uq)][,V1] ]
  DT[    ,Rt := DT[DT.uq,list(DT.id,Rt)][,Rt] ]
  ## Vt: estimated nr of false significant elements in DT for each randomized RRr-value by frequency of RRr-values (R.RRr) in whole null distribution 
  # ffs: fraction of false significants in whole null distribution = ( sum( R.perrow >= perrow ) / prod(dim(R.perrow)) )
  # calculate Vt based on p.diff global
  setkeyv(DT,c("RRr")) # because perrow's are not unique for each c("Ax","Bx")
  DT.uq <- unique(DT) # rewrite DT.uq
  # determine Vt by global null distribution
  DT.uq[ ,Vt := ( sum( R.RRr >= RRr ) / prod(dim(R.RRr)) ) * nrow(DT), by=key(DT.uq)]  
  DT[    ,Vt := DT[DT.uq,list(DT.id,Vt)][,Vt] ]
# # determine Vt based on each voxel's local null distribution 
#  setkeyv(DT,c("M.id"))
#  DT[    ,Vt := ( sum( R.RRr[M.id,] >= RRr ) / R ) * nrow(DT), by=key(DT) ]

  setkeyv(DT,c("DT.id"))
  DT[ ,q.RRr := min(Vt/Rt, 1), by=DT.id]

 message(paste('writing result file ',rdir,mytxt,'result.dat',sep=''))

# store results for mni.id's in Res array using lookup for DT in DT.uq
 Res <- data.table(1:nrow(AB))
 setnames(Res,1,"mni.id")
 # lookup logOR in DT from DT.uq using uq.id
 Res[ ,          Ax:=as.double(NA)] # col 2
 Res[ ,          Ay:=as.double(NA)] # col 3
 Res[ ,          Bx:=as.double(NA)] # col 4
 Res[ ,          By:=as.double(NA)] # col 5
 Res[ ,        row1:=as.double(NA)] # col 6
 Res[ ,        row2:=as.double(NA)] # col 7 
 Res[ ,         RRr:=as.double(NA)] # col 8
 Res[ ,       p.RRr:=as.double(NA)] # col 9
 Res[ ,       q.RRr:=as.double(NA)] # col 10
 Res[ ,  chk.mni.id:=as.double(NA)] # col 11
 
 # rearrange DT as original nii volume order
 setkeyv(DT,"mni.id")
 DT[,row1 := Ax/(Ax+Ay),by=key(DT)] #row1: normalized proportion of lesions given no-deficit [reference]
 DT[,row2 := Bx/(Bx+By),by=key(DT)] #row2: normalized proportion of lesions given deficit 
 #DT[,col1 := By/(Ay+By),by=key(DT)] #col1: normalized proportion of deficits given no-lesion [reference]
 #DT[,col2 := Bx/(Ax+Bx),by=key(DT)] #col2: normalized proportion of deficits given lesion
 set(Res,which(ss),2:11,DT[,list(
  Ax,Ay,Bx,By,
  row1,row2,RRr,p.RRr,q.RRr,mni.id)])
 Res <- cbind(Res,nrs)
 write.table(Res,paste(rdir,mytxt,'result.dat',sep=''),sep='\t')

 tot.vox <- nrow(Res)
 inf.vox <- sum(!is.na(Res[,p.RRr]))
 sig1.vox  <- sum(Res[,q.RRr]<0.40, na.rm=TRUE)
 sig2.vox  <- sum(Res[,q.RRr]<0.20, na.rm=TRUE)
 sig3.vox  <- sum(Res[,q.RRr]<0.05, na.rm=TRUE)
 sig4.vox  <- sum(Res[,p.RRr]<0.20, na.rm=TRUE)
 sig5.vox  <- sum(Res[,p.RRr]<0.10, na.rm=TRUE)
 sig6.vox  <- sum(Res[,p.RRr]<0.05, na.rm=TRUE)
 message(paste("total nr of voxels       : ",tot.vox,sep=''))
 message(paste("informative nr of voxels : ",inf.vox," ( ",round( inf.vox/tot.vox*100, 2) ," % of total ) ", sep=''))
 message(paste("q<0.40 nr of voxels   : ",sig1.vox," ( ",round( sig1.vox/inf.vox*100, 2) ," % of informative ) ", sep=''))
 message(paste("q<0.20 nr of voxels   : ",sig2.vox," ( ",round( sig2.vox/inf.vox*100, 2) ," % of informative ) ", sep=''))
 message(paste("q<0.05 nr of voxels   : ",sig3.vox," ( ",round( sig3.vox/inf.vox*100, 2) ," % of informative ) ", sep=''))
 message(paste("p<0.20 nr of voxels   : ",sig4.vox," ( ",round( sig4.vox/inf.vox*100, 2) ," % of informative ) ", sep=''))
 message(paste("p<0.10 nr of voxels   : ",sig5.vox," ( ",round( sig5.vox/inf.vox*100, 2) ," % of informative ) ", sep=''))
 message(paste("p<0.05 nr of voxels   : ",sig6.vox," ( ",round( sig6.vox/inf.vox*100, 2) ," % of informative ) ", sep=''))
 cat("total nr of voxels       : ",tot.vox,'\n', file=mylog)
 cat("informative nr of voxels : ",inf.vox," ( ",round( inf.vox/tot.vox*100, 2) ," % of total ) \n", file=mylog)
 cat("q < 0.40 nr of voxels : ",sig1.vox," ( ",round( sig1.vox/inf.vox*100, 2) ," % of informative ) \n", file=mylog)
 cat("q < 0.20 nr of voxels : ",sig2.vox," ( ",round( sig2.vox/inf.vox*100, 2) ," % of informative ) \n", file=mylog)
 cat("q < 0.05 nr of voxels : ",sig3.vox," ( ",round( sig3.vox/inf.vox*100, 2) ," % of informative ) \n", file=mylog)
 cat("p < 0.20 nr of voxels : ",sig4.vox," ( ",round( sig4.vox/inf.vox*100, 2) ," % of informative ) \n", file=mylog)
 cat("p < 0.10 nr of voxels : ",sig5.vox," ( ",round( sig5.vox/inf.vox*100, 2) ," % of informative ) \n", file=mylog)
 cat("p < 0.05 nr of voxels : ",sig6.vox," ( ",round( sig6.vox/inf.vox*100, 2) ," % of informative ) \n", file=mylog)

############################
### plotting starts here ###
############################

 message('plotting started ...')

 mni   <- readNIfTI(mfile,reorient=F)
 #cal.min(mni) <- min(mni)
 #cal.max(mni) <- max(mni)

####################
### MAKE NIFTI's ###
####################

#mydim <- dim( example.8mm )

# verification file for use in fslview
 mni.address.nii <- nifti(
   array(Res[,chk.mni.id],dim=mydim), 
   datatype = 16)
 writeNIfTI(mni.address.nii,paste(rdir,mytxt,'mni.address',sep=""))

# convert result to nifti file and save
 message(paste("max %   row1 : ",round(max(Res[,row1],na.rm=TRUE),3),sep=''))
 message(paste("max n   row1 : ",round(max(Res[,Ax],na.rm=TRUE),3),sep=''))
 message(paste("max sum row1 : ",round(max(Res[,Ax+Ay],na.rm=TRUE),3),sep=''))
 message(paste("max %   row2 : ",round(max(Res[,row2],na.rm=TRUE),3),sep=''))
 message(paste("max n   row2 : ",round(max(Res[,Bx],na.rm=TRUE),3),sep=''))
 message(paste("max sum row2 : ",round(max(Res[,Bx+By],na.rm=TRUE),3),sep=''))
 cat("max %   row1 : ",round(max(Res[,row1],na.rm=TRUE),3),'\n', file=mylog)
 cat("max n   row1 : ",round(max(Res[,Ax],na.rm=TRUE),3),'\n', file=mylog)
 cat("max sum row1 : ",round(max(Res[,Ax+Ay],na.rm=TRUE),3),'\n', file=mylog)
 cat("max %   row2 : ",round(max(Res[,row2],na.rm=TRUE),3),'\n', file=mylog)
 cat("max n   row2 : ",round(max(Res[,Bx],na.rm=TRUE),3),'\n', file=mylog)
 cat("max sum row2 : ",round(max(Res[,Bx+By],na.rm=TRUE),3),'\n', file=mylog)

 row1.nii <- nifti(
   array(Res[,row1]/max(Res[,row1],na.rm=TRUE),dim=mydim), # normalized
   datatype = 16)

 row2.nii <- nifti(
   array(Res[,row2]/max(Res[,row2],na.rm=TRUE),dim=mydim),  # normalized
   datatype = 16)

 RRr.nii <- nifti(
   array(Res[,RRr],dim=mydim),
   datatype = 16)

 p.RRr.nii <- nifti(
   array(Res[,p.RRr],dim=mydim),
   datatype = 16)

 q.RRr.nii <- nifti(
   array(Res[,q.RRr],dim=mydim), 
   datatype = 16)

#######################
### CONVERT NIFTI's ###
### TO MNI 1 MM     ###
#######################

message('resampling dims ...')

 resample <- function(nii) {
  # function to resample nii to specified dimensions using ants:ResampleImage
  # nii: nii object in R created by read.NIfTI() to change size 
  # returns nii object in R with desired dimensions
  # fix NAs that become 0's 
    cal.min.old  <- cal.min(nii)
    nii[is.na(nii)] <- -998 # NA
    nii[nii == 0]   <- -999 # 0
    cal.min(nii) <- min(nii)
    writeNIfTI(nii,"tmp")
    system(paste(	
	 "/Users/philip/data/projects/probabilisticmap/ANTS/antsbin/bin/ResampleImage 3 ",
     paste(wd,"tmp.nii.gz",sep=''), " ",
     "tmp_mni.nii.gz ",
     "182x218x182 1 1 ; ", # size:1 and interpolation:1 [nn] MNI 1x1x1mm: 182x218x182
     sep=""))
    out <- readNIfTI(paste("tmp_mni.nii.gz",sep=""),reorient=F)
    out[out == -998] <- NA # replace NAs
    out[out == -999] <- 0  # replace 0s  
    cal.min(out) <- cal.min.old
    return(out)
 }

# resample to MNI 1x1x1 mm 182x218x182 pix
 row1.nii       <- resample(row1.nii)
 row2.nii       <- resample(row2.nii)
 RRr.nii        <- resample(RRr.nii)
 p.RRr.nii      <- resample(p.RRr.nii)
 q.RRr.nii      <- resample(q.RRr.nii)

 cal.min(row1.nii) <- ifelse(is.infinite(min(row1.nii,na.rm=T)),0,min(row1.nii,na.rm=T))
 cal.max(row1.nii) <- ifelse(is.infinite(max(row1.nii,na.rm=T)),1,max(row1.nii,na.rm=T))
 writeNIfTI(row1.nii,paste(rdir,mytxt,'row1',sep=""))

 cal.min(row2.nii) <- ifelse(is.infinite(min(row2.nii,na.rm=T)),0,min(row2.nii,na.rm=T))
 cal.max(row2.nii) <- ifelse(is.infinite(max(row2.nii,na.rm=T)),1,max(row2.nii,na.rm=T))
 writeNIfTI(row2.nii,paste(rdir,mytxt,'row2',sep=""))

 cal.min(RRr.nii) <- ifelse(is.infinite(min(RRr.nii,na.rm=T)),0,min(RRr.nii,na.rm=T))
 cal.max(RRr.nii) <- ifelse(is.infinite(max(RRr.nii,na.rm=T)),1,max(RRr.nii,na.rm=T))
 writeNIfTI(RRr.nii,paste(rdir,mytxt,'RRr',sep=""))

 cal.min(p.RRr.nii) <- ifelse(is.infinite(min(p.RRr.nii,na.rm=T)),0,min(p.RRr.nii,na.rm=T))
 cal.max(p.RRr.nii) <- ifelse(is.infinite(max(p.RRr.nii,na.rm=T)),1,max(p.RRr.nii,na.rm=T))
 writeNIfTI(p.RRr.nii,paste(rdir,mytxt,'p.RRr',sep=""))

 cal.min(q.RRr.nii) <- ifelse(is.infinite(min(q.RRr.nii,na.rm=T)),0,min(q.RRr.nii,na.rm=T))
 cal.max(q.RRr.nii) <- ifelse(is.infinite(max(q.RRr.nii,na.rm=T)),1,max(q.RRr.nii,na.rm=T))
 writeNIfTI(q.RRr.nii,paste(rdir,mytxt,'q.RRr',sep=""))

######################
### COLOR PALETTES ###
######################

 # colorbrewer scheme purple - grey - green
 divcol.cold <- adjustcolor(colorRampPalette( c(
    rgb(208, 28,139,1,maxColorValue=255),
    rgb(208, 28,139,1,maxColorValue=255), # purple
    rgb(241,182,218,1,maxColorValue=255),
    rgb(247,247,247,1,maxColorValue=255), # grey
    rgb(184,225,134,1,maxColorValue=255),
    rgb( 77,172, 38,1,maxColorValue=255), # green
    rgb( 77,172, 38,1,maxColorValue=255)))(1000),alpha.f=0.7)
  # colorbrewer scheme purple - grey - green
  divcol.rr <- adjustcolor(colorRampPalette( c(
     rgb(216,179,101,1,maxColorValue=255), # brown       - 0
     rgb(254,254,254,1,maxColorValue=255), # grey        - 1
     rgb( 90,180,172,1,maxColorValue=255), # turqoise    - 2
     rgb(  1,102, 94,1,maxColorValue=255))  #            - 3
		 )(1000),alpha.f=0.7)
 # colorbrewer scheme blue - yellow - red
 divcol.warm <- adjustcolor(colorRampPalette( c(
    rgb( 44,123,182,1,maxColorValue=255),
    rgb(171,217,233,1,maxColorValue=255),
    rgb(255,255,191,1,maxColorValue=255),
    rgb(253,174, 97,1,maxColorValue=255),
    rgb(215, 25, 28,1,maxColorValue=255)))(1000),alpha.f=0.7)
 # colorbrewer scheme warm greens
 mycol <- adjustcolor(colorRampPalette(c(
    rgb(  0,104, 55,255,maxColorValue=255), # green  - 0.00
    rgb( 49,163, 84,255,maxColorValue=255), #        - 0.20
    rgb(120,198,121,255,maxColorValue=255), #        - 0.40
    rgb(173,221,142,255,maxColorValue=255), #        - 0.60
    rgb(217,240,163,255,maxColorValue=255), #        - 0.80
    rgb(255,255,204,255,maxColorValue=255)) # yellow - 1.00
   )(1000),alpha.f=0.7)
  mycol <- rev(mycol)
	mycol[1] <- adjustcolor("white",alpha.f=0.7)
 # colorbrewer scheme yellow - red - darkred - 0.
 statred <- colorRampPalette( c(
    rgb(255,255,178,255,maxColorValue=255), # yellow - 0.00
    rgb(255,255,178,255,maxColorValue=255), # yellow - 0.03
    rgb(254,204, 92,255,maxColorValue=255), # red    - 0.06
    rgb(253,141, 60,255,maxColorValue=255), # red    - 0.09
    rgb(240, 59, 32,255,maxColorValue=255), # dkred  - 0.12
    rgb(189,  0, 38,255,maxColorValue=255), # dkred  - 0.15
    rep(rgb(189,  0, 38,255,maxColorValue=255),29)) # dkred  - 0.18:0.99
   )(1000)
 # colorbrewer scheme yellow/red/darkred - blue - blue
 statq  <- colorRampPalette( c(
    rgb(255,255,178,255,maxColorValue=255), # yellow - 0.00
  	rgb(254,204, 92,255,maxColorValue=255), # red    - 0.10
    rgb(253,141, 60,255,maxColorValue=255), # red    - 0.20
    rgb(240, 59, 32,255,maxColorValue=255), # dkred  - 0.30
    rgb(189,  0, 38,255,maxColorValue=255), # dkred  - 0.40
    rgb( 44,123,182,255,maxColorValue=255), # blue   - 0.50
    rgb( 44,123,182,255,maxColorValue=255), # blue   - 0.60
    rgb( 44,123,182,255,maxColorValue=255), # blue   - 0.70
    rgb( 44,123,182,255,maxColorValue=255), # blue   - 0.80
    rgb( 44,123,182,255,maxColorValue=255), # blue   - 0.90
    rgb( 44,123,182,255,maxColorValue=255)) # blue   - 1.00
   )(1000)
 # colorbrewer scheme yellow/red - blue - black
 statq1  <- colorRampPalette( c(
    rgb(255,255,178,255,maxColorValue=255), # yellow - 0.00
    rgb(254,204, 92,255,maxColorValue=255), # red    - 0.10
    rgb(253,141, 60,255,maxColorValue=255), # red    - 0.20
    rgb(240, 59, 32,255,maxColorValue=255), # dkred  - 0.30
    rgb(189,  0, 38,255,maxColorValue=255)) # dkred  - 0.40
   )(400)
 statq2  <- colorRampPalette( c(
    rgb( 44,123,182,255,maxColorValue=255), # blue   - 0.50
    rgb( 44,123,182,255,maxColorValue=255), # blue   - 0.75
    rgb( 44,123,182,255,maxColorValue=255)) # blue   - 1.00
   )(600)
 statq <- c(statq1,statq2)
 transseq <- function(a,b){ # transparency vector
    # length 1000
    # a: opaque until eg 200
    # b: transparent from eg 800
    # b>a
    mytransp = 0.75
 	myvec <- vector(length=1000)
    myvec[1:a] <- mytransp
    myvec[b:1000] <- 0
 	for (i in 1:(b-a)) {
 		myvec[a+i] <- mytransp * ((b-a)-i)/(b-a)
 	}
    return(myvec)
 }
 statreda <- statred # add gradual transparency
 transvecp <- transseq(150,500)
 for (i in 1:1000) {
 	statreda[i] <- adjustcolor(statred[i],alpha.f=transvecp[i])
 }
 statqa   <- statq # add gradual transparency
 transvecq <- transseq(400,600)
 for (i in 1:1000) {
 	statqa[i] <- adjustcolor(statq[i],alpha.f=transvecq[i])
 }
 # verify colors
 #pie(rep(1,1000),col=statred,lty=0,labels=NULL)

##################
### MAKE PNG's ###
##################

# png(type = "cairo-png") to avoid color problem in ffmpeg

 slices <- seq(20,150,1)

 message('processing segmentation pngs ...')

### n segmentations
 #cal.min(row1.nii) <- 0.0001 # 0 transparant
 cal.min(row1.nii) <- 0
 cal.max(row1.nii) <- 1
 for (i in 1:length(slices)) {
	png(filename=paste(rdir,mytxt,"row1_",formatC(slices[i],digits=0,wid=3,format="d",flag="0"),".png",sep="" ),width=960,height=960,type = "cairo-png")
	 overlay(mni,row1.nii,z=slices[i],plot.type="single",col.y=mycol)
	dev.off()
 }
 cal.min(row2.nii) <- 0
 cal.max(row2.nii) <- 1
 for (i in 1:length(slices)) {
	png(filename=paste(rdir,mytxt,"row2_",formatC(slices[i],digits=0,wid=3,format="d",flag="0"),".png",sep="" ),width=960,height=960,type = "cairo-png")
	 overlay(mni,row2.nii,z=slices[i],plot.type="single",col.y=mycol)
	dev.off()
 }
 # legend for summary of tumors 
 png(filename=paste(rdir,mytxt,"legend.rowcol.png",sep=""), bg = "black", width = 128, height = 256,type = "cairo-png")
  lg      <- as.matrix(seq(0,1,0.001))
  par(mar=c(0,0,0,6)+0.1,
    oma=c(0,0,0,0)+0.1)
  image(t(lg),
    axes=F,	
    col=mycol)
  axis(4, at=seq(0,1),labels=c("0","max"),las=1,col.axis="white",cex=2)
 dev.off()

 message('processing RR-value pngs ...')

 ### RRr
  tRRr     <- 3 # threshold for RR
  max.RRr  <- max(RRr.nii[!is.na(RRr.nii) & !is.infinite(RRr.nii)]) 
  prc.tRRr <- sum(RRr.nii>tRRr,na.rm=TRUE)/sum(!is.na(RRr.nii))
  message(paste('max for RRr: ',round(max.RRr,2)," , ",round(prc.tRRr*100,1)," %",sep=""))
  tRRr.nii <- RRr.nii
	tRRr.nii[tRRr.nii>tRRr] <- tRRr
  cal.min(tRRr.nii) <- 0
  cal.max(tRRr.nii) <- tRRr
  for (i in 1:length(slices)) {
 	png(filename=paste(rdir,mytxt,"RRr_",formatC(slices[i],digits=0,wid=3,format="d",flag="0"),".png",sep="" ),width=960,height=960,type = "cairo-png")
 	 overlay(mni,tRRr.nii,z=slices[i],plot.type="single",col.y=divcol.rr)
 	dev.off()
  }
  # legend diff
  png(filename=paste(rdir,mytxt,"legend.RRr.png",sep=""), bg = "black", width = 128, height = 256,type = "cairo-png")
   lg      <- as.matrix(seq(0,max.RRr,length.out=1001))
   brks    <- seq(0,max.RRr,length.out=1001)
   par(mar=c(0,0,0,6)+0.3,
     oma=c(0,0,0,0)+0.3)
   image(t(lg),
     axes=F,	
     col=divcol.rr,
     breaks=brks)
   axis(4, at=seq(0,1,length.out=4),labels=c("0","1","2",">3"),las=1,col.axis="white",cex=2)
  dev.off()

  message('processing p-value pngs ...')

### p.RRr
 cal.min(p.RRr.nii) <- 0
 cal.max(p.RRr.nii) <- 1
 for (i in 1:length(slices)) {
	png(filename=paste(rdir,mytxt,"p.RRr_",formatC(slices[i],digits=0,wid=3,format="d",flag="0"),".png",sep=""),width=960,height=960,type = "cairo-png")
	 overlay(mni,p.RRr.nii,z=slices[i],plot.type="single",col.y=statreda)
	dev.off()
 }
 # legend p
  png(filename=paste(rdir,mytxt,"legend.p.png",sep=""), bg = "black", width = 128, height = 256,type = "cairo-png")
   mymax   <- 500
   lg      <- as.matrix(seq(1,mymax,1))
   brks    <- seq(0,mymax,1)
   par(mar=c(0,0,0,6)+0.3,
     oma=c(0,0,0,0)+0.3)
   image(t(lg[1:mymax,]),
     axes=F,	
     col=statreda[1:mymax],
     breaks=brks)
   axis(4, 
    at=seq(0,1,length.out=6),
    labels=formatC(seq(0,mymax/1000,0.1),digits=2,width=3,format="f"),
    las=1,col.axis="white",cex=2)
  dev.off()

 message('processing q-value pngs ...')

### q
 cal.min(q.RRr.nii) <- 0
 cal.max(q.RRr.nii) <- 1
 for (i in 1:length(slices)) {
	png(filename=paste(rdir,mytxt,"q.RRr_",formatC(slices[i],digits=0,wid=3,format="d",flag="0"),".png",sep=""),width=960,height=960,type = "cairo-png")
	 overlay(mni,q.RRr.nii,z=slices[i],plot.type="single",col.y=statqa)
	dev.off()
 }
 # legend q
 png(filename=paste(rdir,mytxt,"legend.q.png",sep=""), bg = "black", width = 128, height = 256,type = "cairo-png")
  mymax <- 600
  lg      <- as.matrix(seq(1,mymax,1))
  brks    <- seq(0,mymax,1)
  par(mar=c(0,0,0,6)+0.3,
    oma=c(0,0,0,0)+0.3)
  image(t(lg[1:mymax,]),
    axes=F,	
    col=statqa[1:mymax],
    breaks=brks)
  axis(4, 
   at=seq(0,1,length.out=7),
   labels=formatC(seq(0,mymax/1000,0.1),digits=2,width=3,format="f"),
   las=1,col.axis="white",cex=2)
 dev.off()

 message('done ...')
 cat(' \n \n', file=mylog)

}