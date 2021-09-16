## SCRIPT to generate the EPrime2 stimuli table for the RS task of the SEEG study.
# By WS, 2020-01-09

# clean up
rm(list=ls())

# setup environments
insert <- function(a,bs,pos){
  as <- split(a,cumsum(seq(a)%in%(pos+1)))
  idx <- order(c(seq_along(as),seq_along(bs)))
  unname(unlist(c(as,bs)[idx]))
}
wdir <- '/data/agora/Chotiga_VOTmultimod/experiment_design/stimuli_selection_SEEG/'
blks.num <- 5
trls.num <- 10  # number of trials per condition in each block
cond.rs.num <- 4  # number of conditions (4 RS conditions and 8 filler conditions)
cond.fillers.num <- 8
catch.num <- trls.num*(cond.rs.num*2+cond.fillers.num)*0.1  # 10% catch trials 

# read the stimuli lists
stimuli <- read.csv(file=file.path(wdir,'Selected_Stimuli_SEEG_RS.csv'))


# generate the Eprime2 table 
# Weight Nested Procedure Bloc StimulusV StimulusA StimulusVDO VDOduration RSE Condition Trigger CorrPesp Jitter
conditions.rs <- rep(c('AA','VV','AV','VA'),each=trls.num)                           # should be 40
conditions.fillers <- rep(c('WA','WV','FA','FV','PA','PV','SA','SV'),each=trls.num)  # should be 80
conditions.fillers.pos <- c(rep(c(1,3),each=10),rep(2,20))
conditions.catch <- rep(c('catchA','catchV'),each=catch.num/2)  
tmp.stimuli <- stimuli
tmp.stimuli.num <- 1
# trial-wise triggers
trigger.wa <- 11  # to 60
trigger.wv <- 61  # to 110
trigger.fa <- 111  # to 160
trigger.fv <- 161  # to 210
for (blk in 1:blks.num){
  ep2.tab <- data.frame(Weight=1,
                        Nested='',
                        Procedure='TrialProcXXX',
                        Bloc='',
                        StimulusV='word',
                        StimulusA='word.wav',
                        StimulusVDO='word.avi',
                        VDOduration=1000,
                        RSE='',
                        Condition='',
                        Trigger=0,
                        CorrPesp='SSS',
                        Jitter=500)
  # generate the conditions sequence
  tmp.rs <- conditions.rs[sample(1:length(conditions.rs))]  # randomized rs sequence
  tmp.fillers <- conditions.fillers[sample(1:length(conditions.fillers))]
  tmp.fillers.pos <- conditions.fillers.pos[sample(1:length(conditions.fillers.pos))]
  tmp.fillers.tag <- list()
  i <- 1
  for (pos in tmp.fillers.pos){
    tmp.fillers.tag[[i]] <- c(tmp.fillers[1:pos])
    tmp.fillers <- tmp.fillers[(pos+1):length(tmp.fillers)]
    i <- i+1
  }
  tmp.catch <- conditions.catch[sample(1:length(conditions.catch))]
  tmp.conditions <- insert(tmp.rs,tmp.fillers.tag,1:length(tmp.rs))
  tmp.conditions <- insert(tmp.conditions,tmp.catch,sample(1:length(tmp.conditions))[1:length(tmp.catch)])
  write.table(tmp.conditions,file=file.path(wdir,sprintf("Conditions_Sequence_RS_Block%02d.txt",blk)),row.names=FALSE,quote=FALSE,col.names=FALSE)
  # convert the conditions sequence into the eprime2 table
  for (j in 1:length(tmp.conditions)){
    tmp.cond.tag <- tmp.conditions[j]
    tmp.cond.stim <- tmp.stimuli[tmp.stimuli$condition==tmp.cond.tag,]
    tmp.stim <- tmp.cond.stim[sample(1:dim(tmp.cond.stim)[1])[1],]
    tmp.stimuli.num <- tmp.stimuli.num+1
#    tmp.stimuli <- tmp.stimuli[tmp.stimuli$stimIdx!=tmp.stim$stimIdx,]
    tmp.stimuli <- subset(tmp.stimuli,stimIdx!=tmp.stim$stimIdx)
    jitter <- sample(500:600)[1]
    if (tmp.cond.tag %in% c('PV','SV')){
      tmp.ep2 <- data.frame(Weight=1,
                            Nested='',
                            Procedure=tmp.stim$proc,
                            Bloc='',
                            StimulusV=tmp.stim$ortho,
                            StimulusA='',
                            StimulusVDO='',
                            VDOduration='',
                            RSE='',
                            Condition=tmp.cond.tag,
                            Trigger=tmp.stim$trigger,
                            CorrPesp='',
                            Jitter=jitter)     
    } else if (tmp.cond.tag %in% c('PA','SA')){
      tmp.ep2 <- data.frame(Weight=1,
                            Nested='',
                            Procedure=tmp.stim$proc,
                            Bloc='',
                            StimulusV='',
                            StimulusA=tmp.stim$sfname,
                            StimulusVDO='',
                            VDOduration='',
                            RSE='',
                            Condition=tmp.cond.tag,
                            Trigger=tmp.stim$trigger,
                            CorrPesp='',
                            Jitter=jitter)  
    } else if (tmp.cond.tag=='WA'){
      tmp.ep2 <- data.frame(Weight=1,
                            Nested='',
                            Procedure=tmp.stim$proc,
                            Bloc='',
                            StimulusV='',
                            StimulusA=tmp.stim$sfname,
                            StimulusVDO='',
                            VDOduration='',
                            RSE='',
                            Condition=tmp.cond.tag,
                            Trigger=trigger.wa,
                            CorrPesp='',
                            Jitter=jitter)      
      trigger.wa <- trigger.wa+1
    } else if (tmp.cond.tag=='WV'){
      tmp.ep2 <- data.frame(Weight=1,
                            Nested='',
                            Procedure=tmp.stim$proc,
                            Bloc='',
                            StimulusV=tmp.stim$ortho,
                            StimulusA='',
                            StimulusVDO='',
                            VDOduration='',
                            RSE='',
                            Condition=tmp.cond.tag,
                            Trigger=trigger.wv,
                            CorrPesp='',
                            Jitter=jitter)      
      trigger.wv <- trigger.wv+1
    } else if (tmp.cond.tag=='FA'){
      tmp.ep2 <- data.frame(Weight=1,
                            Nested='',
                            Procedure=tmp.stim$proc,
                            Bloc='',
                            StimulusV='',
                            StimulusA='',
                            StimulusVDO=tmp.stim$sfname,
                            VDOduration=tmp.stim$duration,
                            RSE='',
                            Condition=tmp.cond.tag,
                            Trigger=trigger.fa,
                            CorrPesp='',
                            Jitter=jitter)      
      trigger.fa <- trigger.fa+1
    } else if (tmp.cond.tag=='FV'){
      tmp.ep2 <- data.frame(Weight=1,
                            Nested='',
                            Procedure=tmp.stim$proc,
                            Bloc='',
                            StimulusV='',
                            StimulusA='',
                            StimulusVDO=tmp.stim$sfname,
                            VDOduration=tmp.stim$duration,
                            RSE='',
                            Condition=tmp.cond.tag,
                            Trigger=trigger.fv,
                            CorrPesp='',
                            Jitter=jitter)      
      trigger.fv <- trigger.fv+1
    } else if (tmp.cond.tag=='AA'){
      jitter.plus <- sample(500:600)[1]
      tmp.ep2 <- data.frame(Weight=c(1,1),
                            Nested=c('',''),
                            Procedure=rep('TrialProcAud',2),
                            Bloc=c('',''),
                            StimulusV=c('',''),
                            StimulusA=rep(tmp.stim$sfname,2),
                            StimulusVDO=c('',''),
                            VDOduration=c('',''),
                            RSE=c('prime','target'),
                            Condition=rep(tmp.cond.tag,2),
                            Trigger=c(1,2),
                            CorrPesp=c('',''),
                            Jitter=c(jitter,jitter.plus))      
    } else if (tmp.cond.tag=='AV'){
      jitter.plus <- sample(500:600)[1]
      tmp.ep2 <- data.frame(Weight=c(1,1),
                            Nested=c('',''),
                            Procedure=c('TrialProcAud','TrialProcVis'),
                            Bloc=c('',''),
                            StimulusV=c('',as.character(tmp.stim$ortho)),
                            StimulusA=c(as.character(tmp.stim$sfname),''),
                            StimulusVDO=c('',''),
                            VDOduration=c('',''),
                            RSE=c('prime','target'),
                            Condition=rep(tmp.cond.tag,2),
                            Trigger=c(3,4),
                            CorrPesp=c('',''),
                            Jitter=c(jitter,jitter.plus)) 
    } else if (tmp.cond.tag=='VA'){
      jitter.plus <- sample(500:600)[1]
      tmp.ep2 <- data.frame(Weight=c(1,1),
                            Nested=c('',''),
                            Procedure=c('TrialProcVis','TrialProcAud'),
                            Bloc=c('',''),
                            StimulusV=c(as.character(tmp.stim$ortho),''),
                            StimulusA=c('',as.character(tmp.stim$sfname)),
                            StimulusVDO=c('',''),
                            VDOduration=c('',''),
                            RSE=c('prime','target'),
                            Condition=rep(tmp.cond.tag,2),
                            Trigger=c(5,6),
                            CorrPesp=c('',''),
                            Jitter=c(jitter,jitter.plus)) 
    } else if (tmp.cond.tag=='VV'){
      jitter.plus <- sample(500:600)[1]
      tmp.ep2 <- data.frame(Weight=c(1,1),
                            Nested=c('',''),
                            Procedure=rep('TrialProcVis',2),
                            Bloc=c('',''),
                            StimulusV=rep(tmp.stim$ortho,2),
                            StimulusA=c('',''),
                            StimulusVDO=c('',''),
                            VDOduration=c('',''),
                            RSE=c('prime','target'),
                            Condition=rep(tmp.cond.tag,2),
                            Trigger=c(7,8),
                            CorrPesp=c('',''),
                            Jitter=c(jitter,jitter.plus))
    } else if (tmp.cond.tag=='catchA'){
      tmp.ep2 <- data.frame(Weight=1,
                            Nested='',
                            Procedure=tmp.stim$proc,
                            Bloc='',
                            StimulusV='',
                            StimulusA=tmp.stim$sfname,
                            StimulusVDO='',
                            VDOduration='',
                            RSE='',
                            Condition=tmp.cond.tag,
                            Trigger=tmp.stim$trigger,
                            CorrPesp='SSS',
                            Jitter=jitter)
    } else if (tmp.cond.tag=='catchV'){
      tmp.ep2 <- data.frame(Weight=1,
                            Nested='',
                            Procedure=tmp.stim$proc,
                            Bloc='',
                            StimulusV=tmp.stim$ortho,
                            StimulusA='',
                            StimulusVDO='',
                            VDOduration='',
                            RSE='',
                            Condition=tmp.cond.tag,
                            Trigger=tmp.stim$trigger,
                            CorrPesp='SSS',
                            Jitter=jitter)
    }
    ep2.tab <- rbind(ep2.tab,tmp.ep2)
  }# finish one block
  write.table(ep2.tab,file=file.path(wdir,sprintf('Eprime2_Block%2d.txt',blk)),sep='\t',row.names=FALSE,quote=FALSE)
}

# add duration info (AUDduration) for audio stimuli 
rm(list=ls())
duration.video <- read.csv(file='Audio-Video_Files_Duration.csv')
duration.pword <- read.csv(file='Stimuli_Auditory_Localizer_Duration.csv')
load('Words_wSoundFiles_Info.Rdata')
StimulusA.duration <- data.frame(stimulus=c(as.character(duration.video$Words),as.character(duration.pword$stimulus[duration.pword$type!='word']),as.character(words.info$sfname)),
                                 duration=c(duration.video$AudioDuraion,duration.pword$duration[duration.pword$type!='word'],words.info$v2length))  # seconds
for (k in 1:5){
  blk.table <- read.csv(file=sprintf('Eprime_RSE_Block%d.csv',k))
  blk.n <- dim(blk.table)[1]
  blk.table$AUDduration <- rep('',blk.n)
  for (q in 1:blk.n){
    tmp.wav <- as.character(blk.table$StimulusA[q])
    if (nchar(tmp.wav)>0){
      if (blk.table$Condition[q] %in% c('catchA','catchV')){
        tmp.wav.dur <- 0.604
      } else {
        tmp.stim <- strtrim(tmp.wav,nchar(tmp.wav)-4)  # remove '.wav'
        tmp.wav.dur <- StimulusA.duration$duration[StimulusA.duration$stimulus==tmp.stim]
      }
      blk.table$AUDduration[q] <- ceiling(tmp.wav.dur*1000)
    }
  }
  write.csv(blk.table,file=sprintf('Eprime_RSE_Block%d_wAUDduration.csv',k),row.names=FALSE,quote=FALSE)
}

# generate expected E-prime timings
rm(list=ls())
mdir <- '/media/wang/BON/Projects/CP01/experiment_design/'
blocks <- c(1,2,3,4,5)
dur_vis <- 550  # in ms
for (i in 1:length(blocks)){
  # read design
  fblk <- file.path(mdir,sprintf('task-RS_Eprime-list-v4_block-%02d.csv',i))
  dblk <- read.csv(fblk, stringsAsFactors = FALSE)
  ntrials <- dim(dblk)[1]
  # calculate expected timings in E-prime
  edat <- data.frame(Trial = double(ntrials), Condition = character(ntrials), FixationCross.OnsetTime = double(ntrials),
                     RSE = character(ntrials), StimAud.OnsetTime = double(ntrials), StimVDO.OnsetTime = double(ntrials),
                     StimVis.OnsetTime = double(ntrials), Trigger = double(ntrials), stringsAsFactors = FALSE) 
  onset <- 0  # in ms
  for (itrial in 1:ntrials){
    # stimulus onset
    if (!is.na(dblk$VDOduration[itrial])){
      edat$StimVDO.OnsetTime[itrial] <- onset
      onset <- onset + dblk$VDOduration[itrial]
    } else if (!is.na(dblk$AUDduration[itrial])){
      edat$StimAud.OnsetTime[itrial] <- onset
      onset <- onset + dblk$AUDduration[itrial]
    } else {
      edat$StimVis.OnsetTime[itrial] <- onset
      onset <- onset + dur_vis
    }
    # fixation onset
    edat$FixationCross.OnsetTime[itrial] <- onset
    # next stimulus onset
    onset <- onset + dblk$Jitter[itrial]
    # other info
    edat$Trial[itrial] <- itrial
    edat$Condition[itrial] <- dblk$Condition[itrial]
    edat$RSE[itrial] <- dblk$RSE[itrial]
    edat$Trigger[itrial] <- dblk$Trigger[itrial]
  }
  # output E-prime list
  fepl <- file.path(mdir, sprintf('sub-xx_ses-01_task-RS_run-01_block-%02d_edat.txt', i))
  write.table(edat, file=fepl, sep='\t', row.names=FALSE, quote=FALSE, fileEncoding = 'UTF-16')
}