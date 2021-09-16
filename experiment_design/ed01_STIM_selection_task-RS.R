## SCRIPT to select stimuli for the RS task of the SEEG study
# By WS, 2020-01-09

# clean up
rm(list=ls())

# setup environments
`%!in%`=Negate(`%in%`)
wdir <- '/data/agora/Chotiga_VOTmultimod/experiment_design/stimuli_selection_SEEG/'

# read various stimuli
# real words used for the RS trials
stim.words <- read.csv(file.path(wdir,'Selected_Words_SEEG_RepetitionSuppression_byConditions.csv'))
# pseudowords and scrambled sounds used as fillers (A&V)
stim.sounds <- read.csv(file.path(wdir,'Stimuli_Auditory_Localizer_Duration.csv'))
stim.pwords <- stim.sounds[stim.sounds$type=='pseudoword',]
stim.scrams <- stim.sounds[stim.sounds$type=='vocoded',]
stim.consos <- read.csv(file=file.path(wdir,'Stimuli_Visual_Localizer.csv'))
stim.consos <- stim.consos[stim.consos$type=='consonant','stimulus']
stim.consos <- stim.consos[sample(1:length(stim.consos))[1:50]]
# real word with video files
stim.videos <- read.csv(file=file.path(wdir,'Audio-Video_Files_Duration.csv'))
# other words that have sound files for the catch trials
load(file=file.path(wdir,'Words_wSoundFiles_Info.Rdata'))  # words.info

# create the stimuli list: ortho, sfname, condition, trigger, proc, duration
# select real words for the RS trials
stim.list.words <- stim.words[,c('ortho','sfname','condition')]
stim.list.words$sfname <- paste0(stim.words$sfname,'.wav')
stim.list.words$trigger <- rep('',dim(stim.words)[1])
stim.list.words$proc <- rep('',dim(stim.words)[1]) 
stim.list.words$duration <- rep('',dim(stim.words)[1]) 
# select baselines (pseudowords and scrambles) for the fillers
stim.pwords <- stim.pwords[order(stim.pwords$duration)[1:50],]
stim.list.bases <- data.frame(ortho=c(rep(as.character(stim.pwords$stimulus),3),as.character(stim.consos)),
                              sfname=c(paste0(rep(stim.pwords$stimulus,2),'.wav'),paste0(stim.pwords$stimulus,'_vocoded.wav'),rep('',50)),
                              condition=rep(c('PA','PV','SA','SV'),each=50),
                              trigger=rep(c(220,230,240,250),each=50),
                              proc=rep(c('TrialProcAud','TrialProcVis','TrialProcAud','TrialProcVis'),each=50),
                              duration=rep('',200))
# select VDO words (WA, WV, FA, FV) for the fillers
stim.list.vdeos <- data.frame(ortho=rep(stim.videos$Words,4),
                              sfname=c(paste0(stim.videos$Words,'.wav'),rep('',50),paste0(stim.videos$Words,'.avi'),paste0(stim.videos$Words,'.wmv')),
                              condition=rep(c('WA','WV','FA','FV'),each=50),
                              trigger=rep('',200),
                              proc=c(rep('TrialProcAud',50),rep('TrialProcVis',50),rep('TrialProcVDO',50),rep('TrialProcVDO',50)),
                              duration=c(rep('',100),rep(stim.videos$AudioDuraion*1000,2)))
# select catch words
stim.catch <- words.info[words.info$ortho %!in% stim.list.words$ortho,]
stim.catchA <- stim.catch[order(stim.catch$v2length)[1:40],]
stim.catchV <- sub("(.{2})(.*)","\\1###\\2",stim.catchA$ortho)                              
stim.list.catch <- data.frame(ortho=c(stim.catchA$ortho,stim.catchV),
                              sfname=c(paste0(stim.catchA$sfname,'_catch.wav'),rep('',40)),
                              condition=rep(c('catchA','catchV'),each=40),
                              trigger=rep(255,80),
                              proc=rep(c('TrialProcAud','TrialProcVis'),each=40),
                              duration=rep('',80))
# combine all stimuli into one stimuli list
stim.list <- rbind(stim.list.words,stim.list.bases,stim.list.vdeos,stim.list.catch)
stim.list$stimIdx <- seq(1:dim(stim.list)[1])
write.csv(stim.list,file=file.path(wdir,'Selected_Stimuli_SEEG_RS.csv'),row.names=FALSE)

# 

