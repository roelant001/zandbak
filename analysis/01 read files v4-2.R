wd <- '/Users/philip/data/projects/probabilisticmap/analysis_cogn/'
setwd(wd)

# v4-2 : separated baseline and deterioration 
# v4-1 : only pts with grade 2-3-4

# TODO
# - 1 mm dimensions?

################
### read nii ###
################

#install.packages('oro.nifti')
#install.packages('xlsx')
#install.packages('lubridate')
library(oro.nifti)
library(xlsx)
library(lubridate)

rdir <- '/Users/philip/data/projects/cognition/data/'
dir2 <- '/Users/philip/data/projects/probabilisticmap/data/npo/2mm/'
mni         <- readNIfTI('/Users/philip/data/projects/probabilisticmap/atlases/MNI152/MNI152_T1_1mm_brain.nii.gz',reorient=F)
example.2mm <- readNIfTI(paste(dir2,'2mm.N001.tum.nii.gz',sep=''),reorient=F)

# read imaging info for baseline data
#info_old <- read.xlsx(paste(rdir,'NPO-final EHE 20140525.xls',sep=''),1,colIndex=c(1,22,69,81,82,94))
info <- read.xlsx(paste(rdir,'NPO-final EHE 20150107.xls',sep=''),1,
 colIndex=c(1,22,23,68,86,87,99), # N nr; include baseline; include deterioration; grade; DOB; gender; OR date
 colClasses=c("character","numeric","numeric","numeric","Date","character","Date"))
info$age <- new_interval(ymd(info$DOB), ymd(info$OR.date)) / duration(num = 1, units = "years")
colnames(info)[1] <- "id"
dim(info)

# function to read nifti
readnii <- function(mydir,x,y='tum') {
 message(paste('N',x,'.',y,' processing ...', sep=""))
 readNIfTI(paste(mydir,'N',x,'.',y,'.nii.gz',sep=''),reorient=F)
}

### BASELINE

# select baseline pts
b <- info[info$include.PWH.baseline==1 & !is.na(info$include.PWH.baseline),]
# exclude pts N029 N047 N048 N066 for deterioration because radiological progression/ further treatment before follow-up NPO, those are included for baseline
# only pts with grade 1-2-3-4 (exclude 2 cavernomas)
b <- b[b$grade >= 1,]
table(b$grade)
dim(b) # 72 pts

### preop tumor based on postop segmentations (res+cav) : tpo ###
Ntpo <- matrix(nrow=prod(dim( example.2mm )),ncol=nrow(b))
for (i in 1:nrow(b)) {
 Ntpo[,i] <- readnii(paste(dir2,"2mm.",sep=''),b$id[i],"tum")
}
colnames(Ntpo) <- as.character(paste("N",b$id,sep=''))
message("Ntpo done")

### preop tumor based on preop segmentations (tmr) : tpr ###
Ntpr <- matrix(nrow=prod(dim( example.2mm )),ncol=nrow(b))
for (i in 1:nrow(b)) {
 Ntpr[,i] <- readnii(paste(dir2,"2mm.",sep=''),b$id[i],"tmr")
}
colnames(Ntpr) <- as.character(paste("N",b$id,sep=''))
message("Ntpr done")

### residual tumor based on postop segmentations (res) : res ###
Nres <- matrix(nrow=prod(dim( example.2mm )),ncol=nrow(b))
for (i in 1:nrow(b)) {
 Nres[,i] <- readnii(paste(dir2,"2mm.",sep=''),b$id[i],"res")
}
colnames(Nres) <- as.character(paste("N",b$id,sep=''))
message("Nres done")

### cavity of tumor resected based on postop tum and res : cav
Ncav <- Ntpo - Nres

mode(Ntpr) <- "logical"
mode(Ntpo) <- "logical"
mode(Nres) <- "logical"
mode(Ncav) <- "logical"

save.image("workspace_01_20150107.RData")

### DETERIORATION

# select deterioration pts
d <- info[info$include.PWI.deterioration==1 & !is.na(info$include.PWI.deterioration),]
# exclude pts N029 N047 N048 N066 for deterioration because radiological progression/ further treatment before follow-up NPO, those are included for baseline
# only pts with grade 1-2-3-4 (exclude 2 cavernomas)
d <- d[d$grade >= 2,]
table(d$grade)
dim(d) # 64 pts

### preop tumor based on postop segmentations (res+cav) : tpo ###
Ntpo <- matrix(nrow=prod(dim( example.2mm )),ncol=nrow(d))
for (i in 1:nrow(d)) {
 Ntpo[,i] <- readnii(paste(dir2,"2mm.",sep=''),d$id[i],"tum")
}
colnames(Ntpo) <- as.character(paste("N",d$id,sep=''))
message("Ntpo done")

### preop tumor based on preop segmentations (tmr) : tpr ###
Ntpr <- matrix(nrow=prod(dim( example.2mm )),ncol=nrow(d))
for (i in 1:nrow(d)) {
 Ntpr[,i] <- readnii(paste(dir2,"2mm.",sep=''),d$id[i],"tmr")
}
colnames(Ntpr) <- as.character(paste("N",d$id,sep=''))
message("Ntpr done")

### residual tumor based on postop segmentations (res) : res ###
Nres <- matrix(nrow=prod(dim( example.2mm )),ncol=nrow(d))
for (i in 1:nrow(d)) {
 Nres[,i] <- readnii(paste(dir2,"2mm.",sep=''),d$id[i],"res")
}
colnames(Nres) <- as.character(paste("N",d$id,sep=''))
message("Nres done")

### cavity of tumor resected based on postop tum and res : cav
Ncav <- Ntpo - Nres # alternative: Ntpr - Nres

mode(Ntpr) <- "logical"
mode(Ntpo) <- "logical"
mode(Nres) <- "logical"
mode(Ncav) <- "logical"

save.image("workspace_01_20150121.RData")


