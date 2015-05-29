#wd <- '~/data/projects/probabilisticmap/analysis_cogn/results'
#setwd(wd)
#mytxt = "baseline_vem_"
# the pngs should be in rgba color format as created by png(type='cairo-png') 

# brew install ffmpeg
# brew update && brew update ffmpeg
do.m4v <- function(x) {
 system(paste("ffmpeg ",
                       # function to make m4v from sequentially numbered pngs
 "-loglevel warning ", # Limit verbose to warnings/errors/panic
 "-r 15 ",             # Read input at framerate of 15 frames per second
 "-vsync 1 ",          # Do not skip frames to match the input rate
 "-f image2 ",         # Input codec – sets it to still images
 "-start_number 020 ", # seq nr to start with instead of _001.png
 "-i ",x,"_%03d.png ", # input png's sequentially numbered _001, _002, _003, etc
 "-pix_fmt yuv420p ",  # generalized color scheme
 "-c:v libx264 ",      # codec to use
 "-preset slower ",    # preset quality protocol
 "-crf 18 ",           # target quality (lower is better) usually 18-25
 "-y ",                # overwrite files without asking
 x,".m4v",             # name of output file
                       #
 sep=""))
 
}

make.movies <- function(mytxt) {

message('make single m4vs ...')

do.m4v(paste(mytxt,"row1",sep=""))
do.m4v(paste(mytxt,"row2",sep=""))
do.m4v(paste(mytxt,"RRr",sep=""))
do.m4v(paste(mytxt,"p.RRr",sep=""))
do.m4v(paste(mytxt,"q.RRr",sep=""))

message('make composite m4v for RRr ...')

 # composite 8 streams using perrow results 
 system(paste("ffmpeg ",
  "-loglevel warning ", # Limit verbose to warnings/errors/panic
  "-i ",mytxt,"row1.m4v ",
  "-i ",mytxt,"row2.m4v ",
  "-i ",mytxt,"RRr.m4v ",
  "-i ",mytxt,"p.RRr.m4v ",
  "-i ",mytxt,"q.RRr.m4v ",
  "-filter_complex \" ",
   "nullsrc=size=4800x960 [base]; ",
   "[0:v] setpts=PTS-STARTPTS, scale=960x960 [U1]; ",
   "[1:v] setpts=PTS-STARTPTS, scale=960x960 [U2]; ",
   "[2:v] setpts=PTS-STARTPTS, scale=960x960 [U3]; ",
   "[3:v] setpts=PTS-STARTPTS, scale=960x960 [U4]; ",
   "[4:v] setpts=PTS-STARTPTS, scale=960x960 [U5]; ",
   "[base][U1] overlay=shortest=1:x=0 [tmp1]; ",
   "[tmp1][U2] overlay=shortest=1:x=960 [tmp2]; ",
   "[tmp2][U3] overlay=shortest=1:x=1920 [tmp3]; ",
   "[tmp3][U4] overlay=shortest=1:x=2880 [tmp4]; ",
   "[tmp4][U5] overlay=shortest=1:x=3840  ",
  " \" ",
  "-c:v libx264 ",
  "-preset slower ", 
  "-crf 18 ",
   "-y ",
  mytxt,"result_RRr.m4v",
  sep="")) 

 system(paste("ffmpeg ",
  "-loglevel warning ", # Limit verbose to warnings/errors/panic
  "-i ",mytxt,"result_RRr.m4v ",
  "-vf scale=1024:768 ",
  "-preset slower ",
  "-crf 18 ", 
  "-y ",
  mytxt,"result_RRr_1024x768.m4v ",
  sep="")) 
   
### cleanup
 
message('make tar.gz  ...')
 
 # tar.gz all x.png's
 system(paste("tar ",
  "-z -c -f ", #-v # g[z]ip [c]reate [v]erbose [f]ile
  mytxt,"pngs.tar.gz ", 
  mytxt,"*.png ", 
  sep=""))
	
 ## brew install trash
 # send x.png's to trash
 system(paste("trash ",
  mytxt,"*.png ", 
  sep=""))
 
}