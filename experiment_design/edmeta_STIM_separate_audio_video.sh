#!/bin/bash

sdir="AVI_files_updated"

cd $sdir
for w in *.avi;do
	echo "separate audio and video files for the word $w ..."
	ffmpeg -i $w -acodec copy A${w:3:-4}.wav
	ffmpeg -i $w -vcodec copy -an VF${w:3:-4}.wmv
done
