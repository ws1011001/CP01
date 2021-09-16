# replay old word by new words
rm(list=ls())

words <- read.csv(file='AudioCut_Files_Duration_selected.csv',as.is=TRUE)

for (i in 1:5){
  blk <- read.csv(file=paste0('RSElist_Eprime_v4_block',i,'.csv'),as.is=TRUE)
  nstim <- dim(blk)[1]
  blk$VDOduration <- rep(NA,nstim)
  for (j in 1:nstim){
    if (blk$Condition[j]=='WV'){
      wold <- blk$StimulusV[j]
      wnew <- words[words$ortho==wold,]
      blk$StimulusV[j] <- wnew$word_new
    } else if (blk$Condition[j]=='WA'){
      wold <- strsplit(blk$StimulusA[j],'.wav',fixed=TRUE)[[1]]
      wnew <- words[words$word_old==wold,]
      blk$StimulusA[j] <- paste0('Ac-',wnew$stimulus,'.wav')
      blk$AUDduration[j] <- wnew$duration_cut*1000 
    } else if (blk$Condition[j]=='FV'){
      wold <- strsplit(blk$StimulusVDO[j],'.wmv',fixed=TRUE)[[1]]
      wnew <- words[words$word_old==wold,]
      blk$StimulusVDO[j] <- paste0('VF-',wnew$stimulus,'.wmv')
      blk$VDOduration[j] <- wnew$offset200
    } else if (blk$Condition[j]=='FA'){
      wold <- strsplit(blk$StimulusVDO[j],'.avi',fixed=TRUE)[[1]]
      wnew <- words[words$word_old==wold,]
      blk$StimulusVDO[j] <- paste0('AVF-',wnew$stimulus,'.avi')     
      blk$VDOduration[j] <- wnew$offset200     
    }
  }
write.csv(blk,file=paste0('RSElist_Eprime_v5_block',i,'.csv'),row.names=FALSE,quote=FALSE,na='')
}