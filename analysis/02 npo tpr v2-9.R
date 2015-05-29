 wd <- '/Users/philip/data/projects/probabilisticmap/analysis_cogn/'
 setwd(wd)
 debug=FALSE

 #install.packages('foreign')
  library(foreign)
  library(xlsx)
  library(oro.nifti)

# load 2 mm data for baseline pts
 load(paste(wd,"workspace_01_20150107.RData",sep=""))	

# load example 2 mm data file for dimensions
 example.2mm <- readNIfTI(paste('/Users/philip/data/projects/probabilisticmap/data/npo/2mm/2mm.N001.tum.nii.gz',sep=''),reorient=F)

## get NPO data ###

# read mri data
 fmri <- '/Users/philip/data/projects/cognition/data/NPO-final EHE 20150107.xls'
 dmri <- read.xlsx(fmri,sheetName="Blad1",
  colIndex=c(1,22,23,27,68,71,76,77,78,79,124),
	colClasses=c("character","integer","integer","integer","numeric","character","numeric","numeric","numeric","numeric","character"))
 colnames(dmri) <- c("N","include.b","include.d","zisnr","grade","side","motorpre","languagepre","motordet","languagedet","languagedominance")

# read npo data
 fnpo <- '/Users/philip/data/projects/cognition/data/NPO pre post domeinscores_3105.sav'
 dnpo <- read.spss(fnpo,to.data.frame=TRUE)
 #cbind(colnames(dnpo))
 dnpo <- dnpo[,
  c(2,     # zis
    20,21, # vem: verbal memory
    27,28, # vic: visual construction
    34,35, # vim: visual memory
    41,42, # isp: infospeed
    48,49, # exe: execution
    55,56, # att: attention
    62,63, # wom: working memory 
    71,72, # bnt: boston naming
    78,79  # mot: motor
    ,22,23,24,25,26 # to verify wom calculations by EHA with following in R
   )] 
 colnames(dnpo) <- c("zisnr",
  "vem1","vem2","vic1","vic2","vim1","vim2",
  "isp1","isp2","exe1","exe2","att1","att2",
  "wom1","wom2","bnt1","bnt2","mot1","mot2",
  "vem_delta","vem_deficitpre","vem_deficitpost","vem_declinepost","vem_change")
 dnpo$zisnr <- as.integer(as.character(dnpo$zisnr)) # make integer for merge
 #dnpo$zisnr 

 mri.b <- dmri[dmri$include.b==1 & !is.na(dmri$include.b),]
 mri.b <- mri.b[mri.b$grade>=1,]
 (bdat <- merge(dnpo,mri.b,by="zisnr"))
 dim(bdat) # 72 pts

 mri.d <- dmri[dmri$include.d==1 & !is.na(dmri$include.d),]
 mri.d <- mri.d[mri.d$grade>=2,]
 (ddat <- merge(dnpo,mri.d,by="zisnr"))
 dim(ddat) # 64 pts

# load compare function for binary maps
#TODO local=TRUE to run sourced parallel: http://stackoverflow.com/questions/16028671/why-does-sourcing-a-script-within-rs-parallel-functions-fail
#TODO alternative: start parallel from global environment with parRapply in sourced function: http://stackoverflow.com/questions/18035711/environment-and-scope-when-using-parallel-functions?rq=1

# parameters
 myR      = 4000 # randomizations
 myFDR    = 0.4  # false discovery rate
 myminpts = 3    # min nr of pts with info in voxel
 mydir    = wd   # dir to write to 
 myxyz    = dim( example.2mm )
 mymni    = '/Users/philip/data/projects/probabilisticmap/atlases/MNI152/MNI152_T1_1mm_brain.nii.gz'

 # start log file
  logname <- paste(wd,"log_",Sys.time(),".txt",sep="")
  logfile <- file(logname,open="a") # opens a connection to this file using [a]ppend mode

  # load source files
  SF1 <- paste(wd,'compare_RPM_bin_v7-7.R',sep='')
	SF2 <- paste(wd,'make_movies_v1-21.R',sep='')
  cat(paste(SF1,'\n'), file=logfile)
  cat(paste(SF2,'\n'), file=logfile)
  source(SF1)
  source(SF2)

################
### baseline ###
################

# preop npo deficit [0: no deficit, 1: deficit, NA] 
# values below -1.5 SD considered deficit
 th <- -1.5
 bdat$vem.pre <- ifelse(bdat$vem1 < th, 1, 0)
 bdat$vic.pre <- ifelse(bdat$vic1 < th, 1, 0)
 bdat$vim.pre <- ifelse(bdat$vim1 < th, 1, 0)
 bdat$isp.pre <- ifelse(bdat$isp1 < th, 1, 0)
 bdat$exe.pre <- ifelse(bdat$exe1 < th, 1, 0)
 bdat$att.pre <- ifelse(bdat$att1 < th, 1, 0)
 bdat$wom.pre <- ifelse(bdat$wom1 < th, 1, 0)
 bdat$bnt.pre <- ifelse(bdat$bnt1 < th, 1, 0)
 bdat$mot.pre <- ifelse(bdat$mot1 < th, 1, 0)
 bdat$nem.pre <- bdat$motorpre # neurological examination preoperative motor findings
 bdat$nel.pre <- bdat$languagepre # neurological examination preoperative language findings

 mytable_baseline <- rbind(
  table(bdat$vem.pre),
  table(bdat$vic.pre),
  table(bdat$vim.pre),
  table(bdat$isp.pre),
  table(bdat$exe.pre),
  table(bdat$att.pre),
  table(bdat$wom.pre),
  table(bdat$bnt.pre),
  table(bdat$mot.pre),
  table(bdat$nem.pre),
  table(bdat$nel.pre))
 rownames(mytable_baseline) <- c("vem","vic","vim","isp","exe","att","wom","bnt","mot","nem","nel")
 colnames(mytable_baseline) <- c("no deficit","deficit [<1.5 sd]")
 mytable_baseline

# verify scores by EHA with these calculations 
 #table(bdat$vem.pre,bdat$vem_deficitpre)
 #cbind(bdat$vem1,bdat$vem.pre)

# test
 if (debug) {
  A=Ntpr[,bdat$nel.pre==0 & !is.na(bdat$nel.pre)] # no deficit 
  B=Ntpr[,bdat$nel.pre==1 & !is.na(bdat$nel.pre)] # deficit
  min.pts=3;R=1000;rdir=mydir;fdr=myFDR;mydim=myxyz;mylog=logfile;mfile=mymni;mytxt='test_'
 compare.RPM.bin(
  Ntpr[,bdat$nel.pre==0 & !is.na(bdat$nel.pre)], # no deficit 
  Ntpr[,bdat$nel.pre==1 & !is.na(bdat$nel.pre)], # deficit
  min.pts=4,R=1000,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
  mytxt='test_')  
 make.movies('test_')
 }
 
 
###################################
### POSITIVE REFERENCE BASELINE ###
###################################

# nel with min.pts 2,3,4,5 all hav q<0.4 for language regions (Broca, Wernicke, Geschwindt, AF/SLF), most (absolute nr of voxels) at min.pts 3
# vem only demonstrates signal q<0.4 at min.pts 3,4, most at 4; but brain region with information in reduced more and more
# therefore min.pts 3 seems reasonable


# baseline verbal memory
  compare.RPM.bin(
   Ntpr[,bdat$vem.pre==0 & !is.na(bdat$vem.pre)], # no deficit 
   Ntpr[,bdat$vem.pre==1 & !is.na(bdat$vem.pre)], # deficit
   min.pts=5,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
   mytxt='baseline_vem_5_')  
  make.movies('baseline_vem_5_')
	compare.RPM.bin(
	 Ntpr[,bdat$vem.pre==0 & !is.na(bdat$vem.pre)], # no deficit 
	 Ntpr[,bdat$vem.pre==1 & !is.na(bdat$vem.pre)], # deficit
	 min.pts=4,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
	 mytxt='baseline_vem_4_')  
	make.movies('baseline_vem_4_')
	compare.RPM.bin(
	 Ntpr[,bdat$vem.pre==0 & !is.na(bdat$vem.pre)], # no deficit 
	 Ntpr[,bdat$vem.pre==1 & !is.na(bdat$vem.pre)], # deficit
	 min.pts=3,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
	 mytxt='baseline_vem_3_')  
	make.movies('baseline_vem_3_')
	compare.RPM.bin(
	 Ntpr[,bdat$vem.pre==0 & !is.na(bdat$vem.pre)], # no deficit 
	 Ntpr[,bdat$vem.pre==1 & !is.na(bdat$vem.pre)], # deficit
	 min.pts=2,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
	 mytxt='baseline_vem_2_')  
	make.movies('baseline_vem_2_')
	compare.RPM.bin(
	 Ntpr[,bdat$vem.pre==0 & !is.na(bdat$vem.pre)], # no deficit 
	 Ntpr[,bdat$vem.pre==1 & !is.na(bdat$vem.pre)], # deficit
	 min.pts=1,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
	 mytxt='baseline_vem_1_')  
	make.movies('baseline_vem_1_')

# baseline motor
  compare.RPM.bin(
   Ntpr[,bdat$mot.pre==0 & !is.na(bdat$mot.pre)], # no deficit 
   Ntpr[,bdat$mot.pre==1 & !is.na(bdat$mot.pre)], # deficit
   min.pts=5,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
   mytxt='baseline_mot_5_')  
  make.movies('baseline_mot_5_')
	compare.RPM.bin(
	 Ntpr[,bdat$mot.pre==0 & !is.na(bdat$mot.pre)], # no deficit 
	 Ntpr[,bdat$mot.pre==1 & !is.na(bdat$mot.pre)], # deficit
	 min.pts=4,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
	 mytxt='baseline_mot_4_')  
	make.movies('baseline_mot_4_')
	compare.RPM.bin(
	 Ntpr[,bdat$mot.pre==0 & !is.na(bdat$mot.pre)], # no deficit 
	 Ntpr[,bdat$mot.pre==1 & !is.na(bdat$mot.pre)], # deficit
	 min.pts=3,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
	 mytxt='baseline_mot_3_')  
	make.movies('baseline_mot_3_')
	compare.RPM.bin(
	 Ntpr[,bdat$mot.pre==0 & !is.na(bdat$mot.pre)], # no deficit 
	 Ntpr[,bdat$mot.pre==1 & !is.na(bdat$mot.pre)], # deficit
	 min.pts=2,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
	 mytxt='baseline_mot_2_')  
	make.movies('baseline_mot_2_')
	compare.RPM.bin(
	 Ntpr[,bdat$mot.pre==0 & !is.na(bdat$mot.pre)], # no deficit 
	 Ntpr[,bdat$mot.pre==1 & !is.na(bdat$mot.pre)], # deficit
	 min.pts=1,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
	 mytxt='baseline_mot_1_')  
	make.movies('baseline_mot_1_')

# baseline neurological examination motor
  compare.RPM.bin(
   Ntpr[,bdat$nem.pre==0 & !is.na(bdat$nem.pre)], # no deficit 
   Ntpr[,bdat$nem.pre==1 & !is.na(bdat$nem.pre)], # deficit
   min.pts=5,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
   mytxt='baseline_nem_5_')  
  make.movies('baseline_nem_5_')
	compare.RPM.bin(
	 Ntpr[,bdat$nem.pre==0 & !is.na(bdat$nem.pre)], # no deficit 
	 Ntpr[,bdat$nem.pre==1 & !is.na(bdat$nem.pre)], # deficit
	 min.pts=4,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
	 mytxt='baseline_nem_4_')  
	make.movies('baseline_nem_4_')
	compare.RPM.bin(
	 Ntpr[,bdat$nem.pre==0 & !is.na(bdat$nem.pre)], # no deficit 
	 Ntpr[,bdat$nem.pre==1 & !is.na(bdat$nem.pre)], # deficit
	 min.pts=3,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
	 mytxt='baseline_nem_3_')  
	make.movies('baseline_nem_3_')
	compare.RPM.bin(
	 Ntpr[,bdat$nem.pre==0 & !is.na(bdat$nem.pre)], # no deficit 
	 Ntpr[,bdat$nem.pre==1 & !is.na(bdat$nem.pre)], # deficit
	 min.pts=2,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
	 mytxt='baseline_nem_2_')  
	make.movies('baseline_nem_2_')
	compare.RPM.bin(
	 Ntpr[,bdat$nem.pre==0 & !is.na(bdat$nem.pre)], # no deficit 
	 Ntpr[,bdat$nem.pre==1 & !is.na(bdat$nem.pre)], # deficit
	 min.pts=1,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
	 mytxt='baseline_nem_1_')  
	make.movies('baseline_nem_1_')

# baseline neurological examination language
  compare.RPM.bin(
   Ntpr[,bdat$nel.pre==0 & !is.na(bdat$nel.pre)], # no deficit 
   Ntpr[,bdat$nel.pre==1 & !is.na(bdat$nel.pre)], # deficit
   min.pts=5,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
   mytxt='baseline_nel_5_')  
  make.movies('baseline_nel_5_')
	compare.RPM.bin(
	 Ntpr[,bdat$nel.pre==0 & !is.na(bdat$nel.pre)], # no deficit 
	 Ntpr[,bdat$nel.pre==1 & !is.na(bdat$nel.pre)], # deficit
	 min.pts=4,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
	 mytxt='baseline_nel_4_')  
	make.movies('baseline_nel_4_')
	compare.RPM.bin(
	 Ntpr[,bdat$nel.pre==0 & !is.na(bdat$nel.pre)], # no deficit 
	 Ntpr[,bdat$nel.pre==1 & !is.na(bdat$nel.pre)], # deficit
	 min.pts=3,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
	 mytxt='baseline_nel_3_')  
	make.movies('baseline_nel_3_')
	compare.RPM.bin(
	 Ntpr[,bdat$nel.pre==0 & !is.na(bdat$nel.pre)], # no deficit 
	 Ntpr[,bdat$nel.pre==1 & !is.na(bdat$nel.pre)], # deficit
	 min.pts=2,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
	 mytxt='baseline_nel_2_')  
	make.movies('baseline_nel_2_')
	compare.RPM.bin(
	 Ntpr[,bdat$nel.pre==0 & !is.na(bdat$nel.pre)], # no deficit 
	 Ntpr[,bdat$nel.pre==1 & !is.na(bdat$nel.pre)], # deficit
	 min.pts=1,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
	 mytxt='baseline_nel_1_')  
	make.movies('baseline_nel_1_')

# baseline bnt
  compare.RPM.bin(
   Ntpr[,bdat$bnt.pre==0 & !is.na(bdat$bnt.pre)], # no deficit 
   Ntpr[,bdat$bnt.pre==1 & !is.na(bdat$bnt.pre)], # deficit
   min.pts=5,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
   mytxt='baseline_bnt_5_')  
  make.movies('baseline_bnt_5_')
	compare.RPM.bin(
	 Ntpr[,bdat$bnt.pre==0 & !is.na(bdat$bnt.pre)], # no deficit 
	 Ntpr[,bdat$bnt.pre==1 & !is.na(bdat$bnt.pre)], # deficit
	 min.pts=4,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
	 mytxt='baseline_bnt_4_')  
	make.movies('baseline_bnt_4_')
	compare.RPM.bin(
	 Ntpr[,bdat$bnt.pre==0 & !is.na(bdat$bnt.pre)], # no deficit 
	 Ntpr[,bdat$bnt.pre==1 & !is.na(bdat$bnt.pre)], # deficit
	 min.pts=3,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
	 mytxt='baseline_bnt_3_')  
	make.movies('baseline_bnt_3_')
	compare.RPM.bin(
	 Ntpr[,bdat$bnt.pre==0 & !is.na(bdat$bnt.pre)], # no deficit 
	 Ntpr[,bdat$bnt.pre==1 & !is.na(bdat$bnt.pre)], # deficit
	 min.pts=2,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
	 mytxt='baseline_bnt_2_')  
	make.movies('baseline_bnt_2_')
	compare.RPM.bin(
	 Ntpr[,bdat$bnt.pre==0 & !is.na(bdat$bnt.pre)], # no deficit 
	 Ntpr[,bdat$bnt.pre==1 & !is.na(bdat$bnt.pre)], # deficit
	 min.pts=1,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
	 mytxt='baseline_bnt_1_')  
	make.movies('baseline_bnt_1_')

##################
#### BASELINE ####
##################

# baseline verbal memory
 compare.RPM.bin(
  Ntpr[,bdat$vem.pre==0 & !is.na(bdat$vem.pre)], # no deficit 
  Ntpr[,bdat$vem.pre==1 & !is.na(bdat$vem.pre)], # deficit
  min.pts=myminpts,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
  mytxt='baseline_vem_')  
 make.movies('baseline_vem_')

# baseline visual construction
 compare.RPM.bin(
  Ntpr[,bdat$vic.pre==0 & !is.na(bdat$vic.pre)], # no deficit 
  Ntpr[,bdat$vic.pre==1 & !is.na(bdat$vic.pre)], # deficit
  min.pts=myminpts,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
  mytxt='baseline_vic_')  
 make.movies('baseline_vic_')

# baseline visual memory
 compare.RPM.bin(
  Ntpr[,bdat$vim.pre==0 & !is.na(bdat$vim.pre)], # no deficit 
  Ntpr[,bdat$vim.pre==1 & !is.na(bdat$vim.pre)], # deficit
  min.pts=myminpts,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
  mytxt='baseline_vim_')  
 make.movies('baseline_vim_')

# baseline infospeed
 compare.RPM.bin(
  Ntpr[,bdat$isp.pre==0 & !is.na(bdat$isp.pre)], # no deficit 
  Ntpr[,bdat$isp.pre==1 & !is.na(bdat$isp.pre)], # deficit
  min.pts=myminpts,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
  mytxt='baseline_isp_')  
 make.movies('baseline_isp_')

# baseline execution
 compare.RPM.bin(
  Ntpr[,bdat$exe.pre==0 & !is.na(bdat$exe.pre)], # no deficit 
  Ntpr[,bdat$exe.pre==1 & !is.na(bdat$exe.pre)], # deficit
  min.pts=myminpts,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
  mytxt='baseline_exe_')  
 make.movies('baseline_exe_')

# baseline attention
 compare.RPM.bin(
  Ntpr[,bdat$att.pre==0 & !is.na(bdat$att.pre)], # no deficit 
  Ntpr[,bdat$att.pre==1 & !is.na(bdat$att.pre)], # deficit
  min.pts=myminpts,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
  mytxt='baseline_att_')  
 make.movies('baseline_att_')

# baseline working memory
 compare.RPM.bin(
  Ntpr[,bdat$wom.pre==0 & !is.na(bdat$wom.pre)], # no deficit 
  Ntpr[,bdat$wom.pre==1 & !is.na(bdat$wom.pre)], # deficit
  min.pts=myminpts,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
  mytxt='baseline_wom_')  
 make.movies('baseline_wom_')

 # baseline neurological examination language : POSITIVE REFERENCE
compare.RPM.bin(
 Ntpr[,bdat$nel.pre==0 & !is.na(bdat$nel.pre)], # no deficit 
 Ntpr[,bdat$nel.pre==1 & !is.na(bdat$nel.pre)], # deficit
 min.pts=myminpts,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
 mytxt='baseline_nel_')  
make.movies('baseline_nel_')
 
# too few pts for bnt and mot.pre
## baseline boston naming test, only in L language dominant pts
# compare.RPM.bin(
#  Ntpr[,bdat$bnt.pre==0 & !is.na(bdat$bnt.pre) & bdat$language == "L"], # no deficit 
#  Ntpr[,bdat$bnt.pre==1 & !is.na(bdat$bnt.pre) & bdat$language == "L"], # deficit
#  min.pts=myminpts,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
#  mytxt='baseline_bnt_')  
# make.movies('baseline_bnt_')

## baseline motor
# compare.RPM.bin(
#  Ntpr[,bdat$mot.pre==0 & !is.na(bdat$mot.pre)], # no deficit 
#  Ntpr[,bdat$mot.pre==1 & !is.na(bdat$mot.pre)], # deficit
#  min.pts=myminpts,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
#  mytxt='baseline_mot_')  
# make.movies('baseline_mot_')
 

#####################
### deterioration ###
#####################

# load 2 mm data for deterioration pts
 load(paste(wd,"workspace_01_20150121.RData",sep=""))	

# deterioration postop - preop npo deficit [0: no deficit, 1: deficit, NA] 
# values more than -1.5 SD considered deficit
 th <- -1.5
 ddat$vem.det <- ifelse(ddat$vem2 - ddat$vem1 < th, 1, 0)
 ddat$vic.det <- ifelse(ddat$vic2 - ddat$vic1 < th, 1, 0)
 ddat$vim.det <- ifelse(ddat$vim2 - ddat$vim1 < th, 1, 0)
 ddat$isp.det <- ifelse(ddat$isp2 - ddat$isp1 < th, 1, 0)
 ddat$exe.det <- ifelse(ddat$exe2 - ddat$exe1 < th, 1, 0)
 ddat$att.det <- ifelse(ddat$att2 - ddat$att1 < th, 1, 0)
 ddat$wom.det <- ifelse(ddat$wom2 - ddat$wom1 < th, 1, 0)
 ddat$bnt.det <- ifelse(ddat$bnt2 - ddat$bnt1 < th, 1, 0)
 ddat$mot.det <- ifelse(ddat$mot2 - ddat$mot1 < th, 1, 0)
 ddat$nem.det <- ddat$motordet # neurological examination preoperative motor findings
 ddat$nel.det <- ddat$languagedet # neurological examination preoperative language findings

#cbind(ddat$vem1,ddat$vem2,ddat$vem.det)

 mytable_deterioration <- rbind(
  table(ddat$vem.det),
  table(ddat$vic.det),
  table(ddat$vim.det),
  table(ddat$isp.det),
  table(ddat$exe.det),
  table(ddat$att.det),
  table(ddat$wom.det),
  table(ddat$bnt.det),
  table(ddat$mot.det),
  table(ddat$nem.det),
  table(ddat$nel.det))
		
 rownames(mytable_deterioration) <- c("vem","vic","vim","isp","exe","att","wom","bnt","mot","nem","nel")
 colnames(mytable_deterioration) <- c("no deterioration","deterioration [drop more than -1.5 sd]")
 mytable_deterioration

 # start log file
  logname <- paste(wd,"log_",Sys.time(),".txt",sep="")
  logfile <- file(logname,open="a") # opens a connection to this file using [a]ppend mode

  # load source files
  SF1 <- paste(wd,'compare_RPM_bin_v7-7.R',sep='')
	SF2 <- paste(wd,'make_movies_v1-21.R',sep='')
  cat(paste(SF1,'\n'), file=logfile)
  cat(paste(SF2,'\n'), file=logfile)
  source(SF1)
  source(SF2)

 ########################################
 ### POSITIVE REFERENCE DETERIORATION ###
 ########################################

 # bnt (7 deficits, 26 no-deficits) with min.pts 2,3,4,5 all have q<0.4 for language regions (Wernicke, insula, AF/SLF), most (absolute nr of voxels) at min.pts 4
 # nel (8 deficits, 56 no-deficits) with min.pts 5 with some q<0.4
 # mot (2 deficits, 25 no-deficits) no q<0.4
 # nem (9 deficits, 55 no-deficits) with min.pts 3,4,5 all have q<0.4 for scattered regions
 # vem (10 deficits, 53 no-deficits) no q<0.4
 # therefore positive reference: bnt and nem with min.pts 3 seems reasonable
 

 # deterioration verbal memory
   compare.RPM.bin(
    Ntpr[,ddat$vem.det==0 & !is.na(ddat$vem.det)], # no deficit 
    Ntpr[,ddat$vem.det==1 & !is.na(ddat$vem.det)], # deficit
    min.pts=5,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
    mytxt='deterioration_vem_5_')  
   make.movies('deterioration_vem_5_')
 	compare.RPM.bin(
 	 Ntpr[,ddat$vem.det==0 & !is.na(ddat$vem.det)], # no deficit 
 	 Ntpr[,ddat$vem.det==1 & !is.na(ddat$vem.det)], # deficit
 	 min.pts=4,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
 	 mytxt='deterioration_vem_4_')  
 	make.movies('deterioration_vem_4_')
 	compare.RPM.bin(
 	 Ntpr[,ddat$vem.det==0 & !is.na(ddat$vem.det)], # no deficit 
 	 Ntpr[,ddat$vem.det==1 & !is.na(ddat$vem.det)], # deficit
 	 min.pts=3,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
 	 mytxt='deterioration_vem_3_')  
 	make.movies('deterioration_vem_3_')
 	compare.RPM.bin(
 	 Ntpr[,ddat$vem.det==0 & !is.na(ddat$vem.det)], # no deficit 
 	 Ntpr[,ddat$vem.det==1 & !is.na(ddat$vem.det)], # deficit
 	 min.pts=2,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
 	 mytxt='deterioration_vem_2_')  
 	make.movies('deterioration_vem_2_')
 	compare.RPM.bin(
 	 Ntpr[,ddat$vem.det==0 & !is.na(ddat$vem.det)], # no deficit 
 	 Ntpr[,ddat$vem.det==1 & !is.na(ddat$vem.det)], # deficit
 	 min.pts=1,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
 	 mytxt='deterioration_vem_1_')  
 	make.movies('deterioration_vem_1_')

 # deterioration motor
   compare.RPM.bin(
    Ntpr[,ddat$mot.det==0 & !is.na(ddat$mot.det)], # no deficit 
    Ntpr[,ddat$mot.det==1 & !is.na(ddat$mot.det)], # deficit
    min.pts=5,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
    mytxt='deterioration_mot_5_')  
   make.movies('deterioration_mot_5_')
 	compare.RPM.bin(
 	 Ntpr[,ddat$mot.det==0 & !is.na(ddat$mot.det)], # no deficit 
 	 Ntpr[,ddat$mot.det==1 & !is.na(ddat$mot.det)], # deficit
 	 min.pts=4,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
 	 mytxt='deterioration_mot_4_')  
 	make.movies('deterioration_mot_4_')
 	compare.RPM.bin(
 	 Ntpr[,ddat$mot.det==0 & !is.na(ddat$mot.det)], # no deficit 
 	 Ntpr[,ddat$mot.det==1 & !is.na(ddat$mot.det)], # deficit
 	 min.pts=3,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
 	 mytxt='deterioration_mot_3_')  
 	make.movies('deterioration_mot_3_')
 	compare.RPM.bin(
 	 Ntpr[,ddat$mot.det==0 & !is.na(ddat$mot.det)], # no deficit 
 	 Ntpr[,ddat$mot.det==1 & !is.na(ddat$mot.det)], # deficit
 	 min.pts=2,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
 	 mytxt='deterioration_mot_2_')  
 	make.movies('deterioration_mot_2_')
 	compare.RPM.bin(
 	 Ntpr[,ddat$mot.det==0 & !is.na(ddat$mot.det)], # no deficit 
 	 Ntpr[,ddat$mot.det==1 & !is.na(ddat$mot.det)], # deficit
 	 min.pts=1,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
 	 mytxt='deterioration_mot_1_')  
 	make.movies('deterioration_mot_1_')

 # deterioration neurological examination motor
   compare.RPM.bin(
    Ntpr[,ddat$nem.det==0 & !is.na(ddat$nem.det)], # no deficit 
    Ntpr[,ddat$nem.det==1 & !is.na(ddat$nem.det)], # deficit
    min.pts=5,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
    mytxt='deterioration_nem_5_')  
   make.movies('deterioration_nem_5_')
 	compare.RPM.bin(
 	 Ntpr[,ddat$nem.det==0 & !is.na(ddat$nem.det)], # no deficit 
 	 Ntpr[,ddat$nem.det==1 & !is.na(ddat$nem.det)], # deficit
 	 min.pts=4,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
 	 mytxt='deterioration_nem_4_')  
 	make.movies('deterioration_nem_4_')
 	compare.RPM.bin(
 	 Ntpr[,ddat$nem.det==0 & !is.na(ddat$nem.det)], # no deficit 
 	 Ntpr[,ddat$nem.det==1 & !is.na(ddat$nem.det)], # deficit
 	 min.pts=3,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
 	 mytxt='deterioration_nem_3_')  
 	make.movies('deterioration_nem_3_')
 	compare.RPM.bin(
 	 Ntpr[,ddat$nem.det==0 & !is.na(ddat$nem.det)], # no deficit 
 	 Ntpr[,ddat$nem.det==1 & !is.na(ddat$nem.det)], # deficit
 	 min.pts=2,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
 	 mytxt='deterioration_nem_2_')  
 	make.movies('deterioration_nem_2_')
 	compare.RPM.bin(
 	 Ntpr[,ddat$nem.det==0 & !is.na(ddat$nem.det)], # no deficit 
 	 Ntpr[,ddat$nem.det==1 & !is.na(ddat$nem.det)], # deficit
 	 min.pts=1,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
 	 mytxt='deterioration_nem_1_')  
 	make.movies('deterioration_nem_1_')

 # deterioration neurological examination language
   compare.RPM.bin(
    Ntpr[,ddat$nel.det==0 & !is.na(ddat$nel.det)], # no deficit 
    Ntpr[,ddat$nel.det==1 & !is.na(ddat$nel.det)], # deficit
    min.pts=5,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
    mytxt='deterioration_nel_5_')  
   make.movies('deterioration_nel_5_')
 	compare.RPM.bin(
 	 Ntpr[,ddat$nel.det==0 & !is.na(ddat$nel.det)], # no deficit 
 	 Ntpr[,ddat$nel.det==1 & !is.na(ddat$nel.det)], # deficit
 	 min.pts=4,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
 	 mytxt='deterioration_nel_4_')  
 	make.movies('deterioration_nel_4_')
 	compare.RPM.bin(
 	 Ntpr[,ddat$nel.det==0 & !is.na(ddat$nel.det)], # no deficit 
 	 Ntpr[,ddat$nel.det==1 & !is.na(ddat$nel.det)], # deficit
 	 min.pts=3,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
 	 mytxt='deterioration_nel_3_')  
 	make.movies('deterioration_nel_3_')
 	compare.RPM.bin(
 	 Ntpr[,ddat$nel.det==0 & !is.na(ddat$nel.det)], # no deficit 
 	 Ntpr[,ddat$nel.det==1 & !is.na(ddat$nel.det)], # deficit
 	 min.pts=2,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
 	 mytxt='deterioration_nel_2_')  
 	make.movies('deterioration_nel_2_')
 	compare.RPM.bin(
 	 Ntpr[,ddat$nel.det==0 & !is.na(ddat$nel.det)], # no deficit 
 	 Ntpr[,ddat$nel.det==1 & !is.na(ddat$nel.det)], # deficit
 	 min.pts=1,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
 	 mytxt='deterioration_nel_1_')  
 	make.movies('deterioration_nel_1_')

 # deterioration bnt
   compare.RPM.bin(
    Ntpr[,ddat$bnt.det==0 & !is.na(ddat$bnt.det)], # no deficit 
    Ntpr[,ddat$bnt.det==1 & !is.na(ddat$bnt.det)], # deficit
    min.pts=5,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
    mytxt='deterioration_bnt_5_')  
   make.movies('deterioration_bnt_5_')
 	compare.RPM.bin(
 	 Ntpr[,ddat$bnt.det==0 & !is.na(ddat$bnt.det)], # no deficit 
 	 Ntpr[,ddat$bnt.det==1 & !is.na(ddat$bnt.det)], # deficit
 	 min.pts=4,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
 	 mytxt='deterioration_bnt_4_')  
 	make.movies('deterioration_bnt_4_')
 	compare.RPM.bin(
 	 Ntpr[,ddat$bnt.det==0 & !is.na(ddat$bnt.det)], # no deficit 
 	 Ntpr[,ddat$bnt.det==1 & !is.na(ddat$bnt.det)], # deficit
 	 min.pts=3,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
 	 mytxt='deterioration_bnt_3_')  
 	make.movies('deterioration_bnt_3_')
 	compare.RPM.bin(
 	 Ntpr[,ddat$bnt.det==0 & !is.na(ddat$bnt.det)], # no deficit 
 	 Ntpr[,ddat$bnt.det==1 & !is.na(ddat$bnt.det)], # deficit
 	 min.pts=2,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
 	 mytxt='deterioration_bnt_2_')  
 	make.movies('deterioration_bnt_2_')
 	compare.RPM.bin(
 	 Ntpr[,ddat$bnt.det==0 & !is.na(ddat$bnt.det)], # no deficit 
 	 Ntpr[,ddat$bnt.det==1 & !is.na(ddat$bnt.det)], # deficit
 	 min.pts=1,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
 	 mytxt='deterioration_bnt_1_')  
 	make.movies('deterioration_bnt_1_')


# deterioration verbal memory
 compare.RPM.bin(
  Ncav[,ddat$vem.det==0 & !is.na(ddat$vem.det)], # no deficit 
  Ncav[,ddat$vem.det==1 & !is.na(ddat$vem.det)], # deficit
  min.pts=myminpts,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
  mytxt='deterioration_vem_')  
 make.movies('deterioration_vem_')

# deterioration visual construction
 compare.RPM.bin(
  Ncav[,ddat$vic.det==0 & !is.na(ddat$vic.det)], # no deficit 
  Ncav[,ddat$vic.det==1 & !is.na(ddat$vic.det)], # deficit
  min.pts=myminpts,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
  mytxt='deterioration_vic_')  
 make.movies('deterioration_vic_')

# deterioration visual memory
 compare.RPM.bin(
  Ncav[,ddat$vim.det==0 & !is.na(ddat$vim.det)], # no deficit 
  Ncav[,ddat$vim.det==1 & !is.na(ddat$vim.det)], # deficit
  min.pts=myminpts,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
  mytxt='deterioration_vim_')  
 make.movies('deterioration_vim_')

# deterioration infospeed
 compare.RPM.bin(
  Ncav[,ddat$isp.det==0 & !is.na(ddat$isp.det)], # no deficit 
  Ncav[,ddat$isp.det==1 & !is.na(ddat$isp.det)], # deficit
  min.pts=myminpts,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
  mytxt='deterioration_isp_')  
 make.movies('deterioration_isp_')

# deterioration execution
 compare.RPM.bin(
  Ncav[,ddat$exe.det==0 & !is.na(ddat$exe.det)], # no deficit 
  Ncav[,ddat$exe.det==1 & !is.na(ddat$exe.det)], # deficit
  min.pts=myminpts,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
  mytxt='deterioration_exe_')  
 make.movies('deterioration_exe_')

# deterioration attention
 compare.RPM.bin(
  Ncav[,ddat$att.det==0 & !is.na(ddat$att.det)], # no deficit 
  Ncav[,ddat$att.det==1 & !is.na(ddat$att.det)], # deficit
  min.pts=myminpts,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
  mytxt='deterioration_att_')  
 make.movies('deterioration_att_')

# deterioration working memory
 compare.RPM.bin(
  Ncav[,ddat$wom.det==0 & !is.na(ddat$wom.det)], # no deficit 
  Ncav[,ddat$wom.det==1 & !is.na(ddat$wom.det)], # deficit
  min.pts=myminpts,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
  mytxt='deterioration_wom_')  
 make.movies('deterioration_wom_')

# deterioration boston naming test, only in L language dominant pts
 compare.RPM.bin(
  Ncav[,ddat$bnt.det==0 & !is.na(ddat$bnt.det) & ddat$language == "L"], # no deficit 
  Ncav[,ddat$bnt.det==1 & !is.na(ddat$bnt.det) & ddat$language == "L"], # deficit
  min.pts=myminpts,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
  mytxt='deterioration_bnt_')  
 make.movies('deterioration_bnt_')

# deterioration motor
 compare.RPM.bin(
  Ncav[,ddat$mot.det==0 & !is.na(ddat$mot.det)], # no deficit 
  Ncav[,ddat$mot.det==1 & !is.na(ddat$mot.det)], # deficit
  min.pts=myminpts,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
  mytxt='deterioration_mot_')  
 make.movies('deterioration_mot_')

 # deterioration neurological examination motor
  compare.RPM.bin(
   Ncav[,ddat$nem.det==0 & !is.na(ddat$nem.det)], # no deficit 
   Ncav[,ddat$nem.det==1 & !is.na(ddat$nem.det)], # deficit
   min.pts=myminpts,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
   mytxt='deterioration_nem_')  
  make.movies('deterioration_nem_')

	# deterioration neurological examination language
	 compare.RPM.bin(
	  Ncav[,ddat$nel.det==0 & !is.na(ddat$nel.det)], # no deficit 
	  Ncav[,ddat$nel.det==1 & !is.na(ddat$nel.det)], # deficit
	  min.pts=myminpts,R=myR,rdir=mydir,fdr=myFDR,mydim=myxyz,mylog=logfile,mfile=mymni,
	  mytxt='deterioration_nel_')  
	 make.movies('deterioration_nel_')

# end log file
 close(logfile) 
