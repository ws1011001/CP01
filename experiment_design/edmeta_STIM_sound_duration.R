## SCRIPT to calculate the duratuon of video and audio files.
# By WS, 2019-12-10.

# clean up
rm(list=ls())

# set environments
library('tuneR')
library('seewave')
wdir <- '/data/agora/Chotiga_VOTmultimod/experiment_design/stimuli_selection_SEEG/'
#sdir <- file.path(wdir,'selected_AV_files/')
sdir <- file.path(wdir,'Stimuli_70dB/')
ndir <- file.path(wdir,'Stimuli_SEEG_70dB/')

# read the names (words) of files
#sf <- read.table(file=file.path(wdir,'selected_AV_files.txt'))
sf <- read.csv(file.path(wdir,'Stimuli_Auditory_Localizer.csv'))
nf <- dim(sf)[1]

# calculate duration
adur <- rep(0,nf)
#vdur <- rep(0,nf)
for (i in 1:nf){
  #wavf <- paste0(sf$stimulus[i],'.wav')    # file name in .wav format
  wavf <- sf$stimulus[i]
  wavo <- readWave(filename=file.path(sdir,wavf))  # wave object
  wavd <- duration(wavo)                           # duration in seconds
#  sprintf('The duration of audio file %s is %f seconds.\n',sf$V1[i],wavd)
  adur[i] <- wavd
#  wmvf <- paste0('/usr/bin/ffprobe -show_entries format=duration -v quiet -i ',file.path(sdir,paste0(sf$V1[i],'.avi')))
#  wmvo <- system(wmvf,intern=TRUE)
#  wmvd <- as.numeric(substr(wmvo[2],10,20))
#  sprintf('The duration of video file %s is %f seconds.\n',sf$V1[i],wmvd)
#  vdur[i] <- wmvd
}

# update sound files info
#sf$V2 <- adur
#sf$V3 <- vdur
#names(sf) <- c('Words','AudioDuraion','VideoDuration')
#write.csv(sf,file=file.path(wdir,'Audio-Video_Files_Duration2.csv'),row.names=FALSE)

sf$duration <- adur
write.csv(sf,file=file.path(wdir,'Stimuli_Auditory_Localizer_Duration.csv'),row.names=FALSE)

# create auditory catch trials
catch.aud <- read.csv(file=file.path(wdir,'selected_words_catch_trials.txt'))
catch.num <- length(catch.aud$pseudowords)
beep <- readWave(filename=file.path(sdir,'BEEP.wav'))
beep <- resamp(beep,g=48000,output = 'Wave')
for (i in 1:catch.num){
  wavf <- paste0(catch.aud$pseudowords[i],'.wav')
  wavo <- readWave(filename=file.path(sdir,wavf),from=0,to=0.3,units='seconds')
  wavn <- bind(wavo,beep)
  writeWave(wavn,filename=file.path(ndir,paste0(catch.aud$pseudowords[i],'_catch.wav')))
}