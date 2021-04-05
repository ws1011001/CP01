#!/bin/bash
## ---------------------------
## [script name] psmeta_visualization_TEP_dynamics.sh
##
## SCRIPT to combine TEP surface figures (that output by HCP workbench), and to make a movie for showing TEP dynamics.
##
## By Shuai Wang, [date] 2021-01-13
##
## ---------------------------
## Notes:
##   
##
## ---------------------------

## set environment (packages, functions, working path etc.)
# setup path
mdir='/media/wang/BON/Projects/CP01'
ddir="$mdir/SEEG_LectureVWFA/derivatives/mia_SEEG_LectureVWFA/sub-08/Figures"
wdir="$mdir/manuscript/visualization"
# visualization parameters
declare -a dmodes=("ERP" "HGA1" "HGA2")
scale='full'  # full or fixed
windows=$(seq 1 39)
fontsize=150
framerate=5
## ---------------------------

## combine figures
for dmode in ${dmodes[@]};do
  echo -e "Make frames for data modality : $dmode ......"
  Adir = "$ddir/sub-08_dynamic_networks_language_${dmode}_Ap"
  Vdir = "$ddir/sub-08_dynamic_networks_language_${dmode}_Vp"
  Cdir = "$ddir/sub-08_dynamic_networks_language_${dmode}_Combined"
  if [ ! -d $Cdir ];then mkdir -p $Cdir;fi
  for w in $windows;do
    echo -e "combine views for time-window $w."
    imgw=`printf "Network%02d.png" $w`          # figure name
    imgl=`printf "Labeled_Network%02d.png" $w`  # new figure name
    # add text to show condition
    convert -pointsize $fontsize -fill red -draw 'text 350,200 "Auditory"' $Adir/$imgw $Adir/$imgl
    convert -pointsize $fontsize -fill blue -draw 'text 350,200 "Visual"' $Vdir/$imgw $Vdir/$imgl
    # combine two conditions vertically
    convert $Adir/$imgl $Vdir/$imgl -append $Cdir/$imgw
  done
  echo -e "Data $dmode is ready for creating a movie."
  # create movies
  echo -e "Make movie for data : $dmode ......"
  ffmpeg -framerate $framerate -i $Cdir/Network%02d.png -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p $wdir/sub-08_dynamic_networks_language_${dmode}_movie.mp4
done
## ---------------------------

