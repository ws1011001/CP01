## SCRIPT to select words as stimuli for the SEEG frequency tagging experiment (CP01).
# By WS, 2019-12-31.

# clean up
rm(list=ls())

# setup environments
library('psych')                 # describeBy()
library('anticlust')             # anticlustering()
`%!in%`=Negate(`%in%`)

# words selection parameters
WordInfo <- c('ortho','phon','sfname','freqlivres','freqfilms2','nbsyll','nblettres','nbphons','old20','pld20','puphon')
VarMatched <- WordInfo[4:length(WordInfo)]
ft.sf.dur <- 0.625  # sound files for the FT task should be no longer than 625ms, that is ~1.6Hz at baseline
ft.sf.num <- 200    # select 200 sound files for the FT task
ft.vw.num <- 200    # select 200 visual words for the FT task
ft.sf.units <- 4    # 50 words per unit; one unit for the VVVA condition, the other three units for the AAAV condition
ft.trials <- 50     # 50 trials per block: the VVVA block and the AAAV block

# read the v4 set of words
load('stimuli_selection/Lexique_Words_v4.Rdata')

# group all words that have sound files
words.al.info <- read.csv('LabView_Scripts/Selected_Words_AduitoryLocalizer.csv')
names(words.al.info)[2] <- 'sfname'
words.mt.info <- read.csv('LabView_Scripts/Selected_Words_MainTask_byConditions.csv')
words.mt.sf <- read.csv('LabView_Scripts/Stimuli_MainTaskAuditory.csv')
words.mt.info <- merge(words.mt.info,words.mt.sf[,c('ortho','sfname')],by='ortho')
words.info <- rbind(words.al.info[,WordInfo],words.mt.info[,WordInfo])
words.sf <- read.csv('LabView_Scripts/Sound_Files/SoundLength_v2.csv')
names(words.sf)[2] <- 'sfname'
words.info <- merge(words.info,words.sf[,c('sfname','v2length')],by.x='sfname')
save(words.info,file='Words_wSoundFiles_Info.Rdata')

# select sound files for the frequency-tagging task
task.ft.sf.pool <- words.info[words.info$v2length<=ft.sf.dur,]
task.ft.sf.pool$v2lengthOrd <- abs(task.ft.sf.pool$v2length - median(task.ft.sf.pool$v2length))
task.ft.sf <- task.ft.sf.pool[order(task.ft.sf.pool$v2lengthOrd)[1:ft.sf.num],]

# complete the words selection for the FT task
UnitComparisonP <- rep(0,length(VarMatched))
while (any(UnitComparisonP < 0.3)){
  # split the selected sound files into $ft.sf.units matched units
  task.ft.sf.units <- anticlustering(features=task.ft.sf[,VarMatched],K=ft.sf.units,objective='variance')
  # compare target variables across all units with Kruskal-Wallis Rank Sum Test
  UnitComparison <- lapply(task.ft.sf[,VarMatched], function(x) kruskal.test(x,as.factor(task.ft.sf.units)))
  UnitComparisonP <- unlist(lapply(UnitComparison, function(x) x$p.value))
}
print(paste('#####Quick Check - Units##### The Comparison results are: P =',UnitComparisonP))
task.ft.sf$units <- task.ft.sf.units  # auditory units: 1,2,3,4
task.ft.vw.pool <- words.all[words.all$ortho %!in% task.ft.sf$ortho,]  # remove the selected sound files
FTComparisonP <- rep(0,length(VarMatched))
while (any(FTComparisonP < 0.3)){
  randvw <- sample(c(rep(1,ft.vw.num),rep(0,dim(task.ft.vw.pool)[1]-ft.vw.num)))
  task.ft.pre <- rbind(task.ft.vw.pool[randvw!=0,VarMatched],task.ft.sf[,VarMatched])
  # compare target variables between sound files and visual words
  FTComparison <- lapply(task.ft.pre, function(x) wilcox.test(x[1:ft.sf.num],x[(ft.sf.num+1):(ft.sf.num+ft.vw.num)]))
  FTComparisonP <- unlist(lapply(FTComparison, function(x) x$p.value))
}
print(paste('#####Quick Check - Units##### The Comparison results are: P =',FTComparisonP))
task.ft.vw <- task.ft.vw.pool[as.logical(randvw),]  # the selected visual words
UnitVWComparisonP <- rep(0,length(VarMatched))
while (any(UnitVWComparisonP < 0.3)){
  # split the selected sound files into $ft.sf.units matched units
  task.ft.vw.units <- anticlustering(features=task.ft.vw[,VarMatched],K=ft.sf.units,objective='variance')
  # compare target variables across all units with Kruskal-Wallis Rank Sum Test
  UnitVWComparison <- lapply(task.ft.vw[,VarMatched], function(x) kruskal.test(x,as.factor(task.ft.vw.units)))
  UnitVWComparisonP <- unlist(lapply(UnitVWComparison, function(x) x$p.value))
}
task.ft.vw$units <- task.ft.vw.units+4  # visual units: 5,6,7,8

# output the trials list and the words list for the FT task
task.ft.vw$sfname <- task.ft.vw$ortho
task.ft <- rbind(task.ft.sf[,c(WordInfo,'units')],task.ft.vw[,c(WordInfo,'units')])
task.ft$modality <- rep(c('Aud','Vis'),each=ft.sf.num)
task.ft.aaav <- data.frame(base1=task.ft$sfname[task.ft$units==1],
                           base2=task.ft$sfname[task.ft$units==2],
                           base3=task.ft$sfname[task.ft$units==3],
                           target=task.ft$ortho[task.ft$units==8])
task.ft.vvva <- data.frame(base1=task.ft$ortho[task.ft$units==5],
                           base2=task.ft$ortho[task.ft$units==6],
                           base3=task.ft$ortho[task.ft$units==7],
                           target=task.ft$sfname[task.ft$units==4])
task.ft$block <- rep('VVVA',ft.sf.num+ft.vw.num)
task.ft$block[task.ft$units %in% c(1,2,3,8)] <- 'AAAV'
write.csv(task.ft.aaav,file='Selected_Words_FrequencyTagging_AAAV_Trials.csv',row.names=FALSE)
write.csv(task.ft.vvva,file='Selected_Words_FrequencyTagging_VVVA_Trials.csv',row.names=FALSE)
write.csv(task.ft,file='Selected_Words_FrequencyTagging_byBlocks.csv',row.names=FALSE)

# compare target variables between VVVA block and AAAV block
BlockComparison <- lapply(task.ft[,VarMatched], function(x) wilcox.test(x[task.ft$block=='VVVA'],x[task.ft$block=='AAAV']))