## SCRIPT to calculate the duratuon of video and audio files.
# By WS, 2019-12-10.

# clean up
rm(list=ls())

# set environments
library('tuneR')
library('seewave')
wdir <- '.'
sdir <- file.path(wdir,'audioOnly/')
ndir <- file.path(wdir,'audioCut/')

# read the names (words) of files
sf <- read.csv(file.path(wdir,'AudioOnsetOffset.csv'))
nf <- dim(sf)[1]

# calculate original duration and cut off sound waves
adur <- rep(0,nf)  # original duration
cdur <- rep(0,nf)  # cut duration
for (i in 1:nf){
  # load up audio file
  wavf <- paste0('A-',sf$stimulus[i],'.wav')       # file name in .wav format
  wavo <- readWave(filename=file.path(sdir,wavf))  # wave object
  # calculate original duration (in seconds)
  adur[i] <- duration(wavo)                           
  message(sprintf('The original duration of audio file %s is %f seconds.\n',sf$stimulus[i],adur[i]))
  # show acoustic onset and offset (in seconds)
  onset <- sf$onset[i]/1000
  offset <- sf$offset[i]/1000
  message(sprintf('The onset and offset of spoken word %s are %f and %f respectively. Cut off it in this range.\n',sf$stimulus[i],onset,offset))
  # cut off pre-onset and post-offset
  wav1 <- cutw(wavo,channel=1,from=onset,to=offset,output='Wave')
  wav2 <- cutw(wavo,channel=2,from=onset,to=offset,output='Wave') 
  wavn <- stereo(wav1,wav2)
  cdur[i] <- duration(wavn)
  writeWave(wavn,filename=file.path(ndir,paste0('Ac-',sf$stimulus[i],'.wav')))
  # output sound wave plot
  png(filename=file.path(ndir,paste0('Ac-',sf$stimulus[i],'.png')))
    cutw(wavo,from=onset,to=offset,plot=TRUE)
  dev.off()
}

# update sound files info
sf$duration_original <- adur
sf$duration_cut <- cdur
write.csv(sf,file=file.path(wdir,'AudioCut_Files_Duration.csv'),row.names=FALSE)

