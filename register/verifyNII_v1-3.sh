#!/bin/bash

# function to cycle through folders to verify tum res in mni space
# executed in Terminal, e.g.: 
# first make the file executable : chmod 755 ./verifyNII_v1-2.sh
# cd /Users/philip/Documents/projects/probabilisticmap/bash
# /Users/philip/Documents/projects/probabilisticmap/bash/verifyNII_v1-2.sh /Users/philip/Documents/projects/probabilisticmap/brainlab/20131031 A022 
# /Users/philip/Documents/projects/probabilisticmap/bash/verifyNII_v1-2.sh /Users/philip/Documents/projects/probabilisticmap/images_expert_atlas/ATLASvumcA001-063 A026 
# ./verifyNII_v1-2.sh /Users/philip/Documents/projects/probabilisticmap/images_expert_atlas/ATLASfranceB001-059 B001 
# ./verifyNII_v1-2.sh /Users/philip/Documents/projects/probabilisticmap/brainlab/20131031 A022 
# /Users/philip/Documents/projects/probabilisticmap/bash/verifyNII_v1-2.sh /Users/philip/Documents/projects/probabilisticmap/brainlab/20131109 A015 
# /Users/philip/Documents/projects/probabilisticmap/bash/verifyNII_v1-2.sh /Users/philip/Documents/projects/probabilisticmap/brainlab/20131110 A015 
# /Users/philip/Documents/projects/probabilisticmap/bash/verifyNII_v1-2.sh /Users/philip/Documents/projects/probabilisticmap/images_expert_atlas A001
# /Users/philip/Documents/projects/probabilisticmap/bash/verifyNII_v1-2.sh /Users/philip/Documents/projects/probabilisticmap/images_expert_atlas A051
# /Users/philip/Documents/projects/probabilisticmap/bash/verifyNII_v1-2.sh /Volumes/RPM/images_expert_atlas_post A093
# /Users/philip/data/projects/probabilisticmap/bash/verifyNII_v1-3.sh /Users/philip/data/presentations/PICTURE/ A052

myDataDir=$1
myFirst=$2
init=0
myAtlas='/Users/philip/data/projects/probabilisticmap/atlases/MNI152/MNI152_T1_1mm.nii.gz'
SnapApp='/Applications/ITK-SNAP.app/Contents/MacOS/InsightSNAP'
antsDir='/Users/philip/data/projects/probabilisticmap/ANTS/antsbin/bin'

dirList=`ls $myDataDir`
for mySampleDir in $dirList ; do

   PatID=$mySampleDir
   myOutDir=$myDataDir/$mySampleDir/OUT

   if [ "$PatID" = "$myFirst" ] ; then
    init=1
   fi
   
   if [[ $init = 1 && -d $myOutDir ]] ; then

    echo "$PatID"
   
    # correct cal-min/cal-max for visualization
    min=`fslstats $myOutDir/mprpost.org.nii.gz -r | awk '{print $1}'`
    max=`fslstats $myOutDir/mprpost.org.nii.gz -r | awk '{print $2}'`

    fslview -m ortho $myAtlas -l "Blue-Lightblue" -t 1.0 \
     $myOutDir/$PatID.tum.nii.gz -l "Green" -t 0.5 \
     $myOutDir/$PatID.res.nii.gz -l "Red"  -t 0.5 \
     $myOutDir/mprpost.mni.nii.gz -l "Copper" -t 0.8 -b $min,$max 
   
   
   fi
    
done     
	