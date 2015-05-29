#!/bin/bash

# function myprocess uses 2 arguments: <folder to process> <o:no_overwrite 1:overwrite>
# executed in Terminal, e.g.: 
# cd /Users/philip/Documents/projects/probabilisticmap/bash
# ./myprocess_v1-5.sh /Users/philip/Documents/projects/probabilisticmap/brainlab/20131011 1 | tee log_20131028.txt
# ./myprocess_v1-9.sh /Users/philip/Documents/projects/probabilisticmap/brainlab/20131109 1 | tee log_20131110.txt
# output to terminal als log.txt file

### history
# v1-9c
# - for npo_pre folder
# - A000.tmr.nii.gz copied and changed to N000.tmr.nii.gz
# v1-9b
# - changed for mprpre and tmr (as preop tumor)
# v1-9
# - inserted intermediate myA_template_trimmed atlas to correct cortex mismatch
# - NearestNeighbor instead of Linear for mask (tum,res) warps, so that threshold step is unnecessary
# - separate doWarpMPR function outputs mprpost.mya and mprpost.mni
# - changed SyN metric form Mattes to CC for precision
# - tumor dilated from 4 to 3
# - more preciion for affine 1.e-6 
# v1-8
# - tum.inv dilated 4
# - CC[,,4]
# - affine 4 stages
# - use myA_template instead of MNI

### initialize some parameters
 myDataDir=$1
 # if $2 is missing or not 1 do not overwrite, else if $2 is 1 overwrite
 myResDir='/Users/philip/Documents/projects/probabilisticmap/data/npo'

 myAtlas='/Users/philip/Documents/projects/probabilisticmap/atlases/buildtemplate/myA_template_trimmed.nii.gz'
 myMNI='/Users/philip/Documents/projects/probabilisticmap/atlases/MNI152/MNI152_T1_1mm.nii.gz'
 myAtlasBrain='/Users/philip/Documents/projects/probabilisticmap/data/results/EM_Map_mask_edit.nii.gz'
 antsDir='/Users/philip/Documents/projects/probabilisticmap/ANTS/antsbin/bin'

### functions

DICOMtoNII() {
# convert dicoms to nii using default settings
# with 1 argument : mpr, cav, res, tmr

   dcm2nii \
    -a n -c y -d y -e n -f n -g y -i n -m n -n y \
    -o $myOutDir -p y -r n -s n -v y -x n "$myFolder/$1"

   # rename last modified file
   outDirList=`ls $myOutDir`
   arr=($outDirList)
   latestfile=${arr[0]}
   for filename in $outDirList ; do     
    if test "$myOutDir/$filename" -nt "$myOutDir/$latestfile" ; then
       latestfile=$filename
    fi
   done    
   mv "$myOutDir/$latestfile" "$myOutDir/$1.org.nii.gz"
   echo "renamed $myOutDir/$latestfile into $myOutDir/$1.org.nii.gz"

}   

makeThresholds() {
# function makeThresholds
# with 1 argument : cav, res, tmr
   
   # create binary mask
   #c3d $myOutDir/$1.org.nii.gz -threshold 65535 inf 1 0 -o $myOutDir/$1.mpr.nii.gz
   # alternative using fslmaths
   # fslmaths $myOutDir/$1.org.nii.gz -thr 65535 -bin $myOutDir/$1.mpr.nii.gz 
   # $antsDir/ImageMath 3 $myOutDir/$1.mpr.nii.gz \
   # ReplaceVoxelValue $myOutDir/$1.org.nii.gz 0 65534 0 \
   # ReplaceVoxelValue $myOutDir/$1.org.nii.gz 65535 65535 1 
   $antsDir/ThresholdImage 3 \
    $myOutDir/$1.org.nii.gz \
    $myOutDir/$1.mprpre.nii.gz \
    0 65534 0 1

   echo "binarized $myOutDir/$1.mprpre.nii.gz"

}

oneMNI() {
# function create null mni

   # create empty volume with mni atlas dimensions
   #c3d $myAtlas -threshold -inf inf 0 0 -o $myOutDir/mni.nul.nii.gz 
   #fslmaths $myAtlas -thr 100000 $myOutDir/mni.nul.nii.gz 
   $antsDir/ImageMath 3 $myOutDir/mya.one.nii.gz \
    ReplaceVoxelValue $myAtlas 0 65535 1 

   echo "created empty $myOutDir/mya.one.nii.gz"

}

addNII() {
# function to add two niftis and create a binarized output and inverted binarized output
# 3 arguments : $1 + $2 = $3, eg cav + res = tum
   # alternative in fsl, including binarization
   #fslmaths $myOutDir/$1.mpr.nii.gz -add $myOutDir/$2.mpr.nii.gz \
   # -bin $myOutDir/$3.mpr.nii.gz
   # alternative in fsl, including inverted binarization
   #fslmaths $myOutDir/$1.mpr.nii.gz -add $myOutDir/$2.mpr.nii.gz \
   # -binv $myOutDir/$3.inv.nii.gz

#   # add $1 and $2
#   $antsDir/ImageMath 3 \
#    $myOutDir/$3.mprpost.nii.gz \
#    + $myOutDir/$1.mprpost.nii.gz $myOutDir/$2.mprpost.nii.gz
#   # binarize $3
#   $antsDir/ThresholdImage 3 \
#    $myOutDir/$3.mprpost.nii.gz \
#    $myOutDir/$3.mprpost.nii.gz \
#    1 65535 1 0
#
#   # dilate with 3 voxels for tum.inv
#   $antsDir/ImageMath 3 \
#    $myOutDir/$3.inv.nii.gz \
#    MD $myOutDir/$3.mprpost.nii.gz 3
#   # inverse binary so that 0=mask of dilated tumor and 1=the rest
#   $antsDir/ImageMath 3 \
#    $myOutDir/$3.inv.nii.gz \
#    - $myOutDir/$3.inv.nii.gz 1
#   $antsDir/ImageMath 3 \
#    $myOutDir/$3.inv.nii.gz \
#    abs $myOutDir/$3.inv.nii.gz

   # dilate with 3 voxels for tum.inv
   $antsDir/ImageMath 3 \
    $myOutDir/$1.inv.nii.gz \
    MD $myOutDir/$1.mprpre.nii.gz 3
   # inverse binary so that 0=mask of dilated tumor and 1=the rest
   $antsDir/ImageMath 3 \
    $myOutDir/$1.inv.nii.gz \
    - $myOutDir/$1.inv.nii.gz 1
   $antsDir/ImageMath 3 \
    $myOutDir/$1.inv.nii.gz \
    abs $myOutDir/$1.inv.nii.gz

   echo "created $myOutDir/$1.mprpre.nii.gz" 

}

doN4() {

     $antsDir/N4BiasFieldCorrection 3 \
      -i $myOutDir/mprpre.org.nii.gz \
      -b [200] -s 3 -c [50x50x30x20,1e-6] \
      -o [$myOutDir/mprpre.tmp.nii.gz, $myOutDir/mprpre.n4b.nii.gz]

}

doants() {

    f=$myAtlas
    m=$myOutDir/mprpre.org.nii.gz

    $antsDir/antsRegistration \
     -d 3 \
     -r [ $f, $m ,1] \
     -m mattes[  $f, $m , 1 , 32, regular, 0.1 ] \
      -t translation[ 0.1 ] \
      -c [1000x1000x1000,1.e-7,10]  \
      -s 4x2x1vox  \
      -f 6x4x2 -l 1 \
     -m mattes[  $f, $m , 1 , 32, regular, 0.1 ] \
      -t rigid[ 0.1 ] \
      -c [1000x1000x1000,1.e-7,10]  \
      -s 4x2x1vox  \
      -f 6x4x2 -l 1 \
     -m mattes[  $f, $m , 1 , 32, regular, 0.25 ] \
      -t affine[ 0.1 ] \
      -c [ 1000x1000x1000x50,1.e-6,10 ] \
      -s 6x2x1x1vox  \
      -f 8x4x2x1 -l 1 \
     -m cc[  $f, $m , 1 , 4] \
      -t SyN[ 0.2, 3, 0 ] \
      -c [ 100x100x100, 1.e-4, 5 ] \
      -s 6x2x1vox \
      -f 8x4x2 -l 1 -u 1 -z 1 \
     -m cc[  $f, $m , 1 , 4] \
      -t BSplineSyN[0.25,10x10x10,0x0x0,3] \
      -c [ 20, 1.e-4, 5 ] \
      -s 1vox \
      -f 2 -l 1 -u 1 -z 1 \
     -x [$myOutDir/mya.one.nii.gz,$myOutDir/tmr.inv.nii.gz] \
     -o [$myOutDir/mprpre.mya,$myOutDir/mprpre.mya_diff.nii.gz,$myOutDir/mprpre.mya_inv.nii.gz]

}

doWarpMPR() {

    $antsDir/antsApplyTransforms \
     -d 3 -i $myOutDir/mprpre.org.nii.gz -r $myAtlas -n linear \
     -t $myOutDir/mprpre.mya1Warp.nii.gz \
     -t $myOutDir/mprpre.mya0GenericAffine.mat \
     -o $myOutDir/mprpre.mya.nii.gz

    $antsDir/antsApplyTransforms \
     -d 3 -i $myOutDir/mprpre.org.nii.gz -r $myMNI -n linear \
     -t /Users/philip/Documents/projects/probabilisticmap/atlases/buildtemplate/mya.mni1Warp.nii.gz \
     -t /Users/philip/Documents/projects/probabilisticmap/atlases/buildtemplate/mya.mni0GenericAffine.mat \
     -t $myOutDir/mprpre.mya1Warp.nii.gz \
     -t $myOutDir/mprpre.mya0GenericAffine.mat \
     -o $myOutDir/mprpre.mni.nii.gz

}

doSimilarity() {

    $antsDir/MeasureImageSimilarity 3 1 \
     $f $myOutDir/mprpre.mni.nii.gz $myOutDir/similarity

}

doWarp() {

    $antsDir/antsApplyTransforms \
     -d 3 -i $myOutDir/$1.mprpre.nii.gz -r $myMNI -n NearestNeighbor \
     -t /Users/philip/Documents/projects/probabilisticmap/atlases/buildtemplate/mya.mni1Warp.nii.gz \
     -t /Users/philip/Documents/projects/probabilisticmap/atlases/buildtemplate/mya.mni0GenericAffine.mat \
     -t $myOutDir/mprpre.mya1Warp.nii.gz \
     -t $myOutDir/mprpre.mya0GenericAffine.mat \
     -o $myOutDir/$1.mni.nii.gz
   echo "warped $myOutDir/$1.mprpre.nii.gz to $myOutDir/$1.mni.nii.gz"

}

makeResult() {
# limit info to brain mask of atlas only by multiplying result with brain mask
# copy result file to result folder

   $antsDir/ImageMath 3 \
    $myOutDir/$PatID.$1.nii.gz \
    m $myOutDir/$1.mni.nii.gz $myAtlasBrain

   cp $myOutDir/$PatID.$1.nii.gz $myResDir/$PatID.$1.nii.gz
   echo "saved $myOutDir/$PatID.$1.nii.gz and copied to $myResDir/"

}

# function doIt
# as wrapper
function doIt { 

   DICOMtoNII mprpre
   DICOMtoNII tmr
   makeThresholds tmr
   oneMNI
   addNII tmr
   doants
   doWarpMPR
   doSimilarity
   doWarp tmr
   makeResult tmr
   
}

### MAIN EXECUTABLE STARTS HERE

### date & time stamp
date

### cycle through folders of myDataDir
dirList=`ls $myDataDir`
for mySampleDir in $dirList ; do
   
   PatID=$mySampleDir
   myFolder="$myDataDir/$mySampleDir"
   myOutDir="$myFolder/OUT"
   if test ! -d $myOutDir ; then
    if test ! -f $myOutDir ; then
      mkdir $myOutDir
    fi
   fi

   if test -d $myOutDir ; then

    echo "################"  
    echo "PROCESSING: $PatID"
    t0=$SECONDS

    # test if overwrite is required
    if [[ $2 != 1 && -f $myOutDir/$PatID.tmr.nii.gz ]] ; then
     echo "$myOutDir/ already processed"
     cp $myOutDir/*.tmr.nii.gz $myResDir/$PatID.tmr.nii.gz
    else
     doIt
    fi

    # print time
    t1=$SECONDS
    tsum=`echo $t1 - $t0  | bc -l`
    echo "time: $tsum secs"
    echo "################"
    echo " "

   fi

done

### copy log file to analysis dir
cp ./log.txt $1/
