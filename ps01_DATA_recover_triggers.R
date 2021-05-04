## ---------------------------
## [script name] ps01_DATA_recover_triggers.R
##
## SCRIPT to recover missing triggers based on E-prime logs.
##
## By Shuai Wang, [date]
##
## ---------------------------
## Notes: - This script shold be performed on certain subject whose triggers are incomplete.
##   
##
## ---------------------------

## clean up
rm(list=ls())
## ---------------------------

## set environment (packages, functions, working path etc.)
# setup working path
subj <- 'sub-test'
mdir <- '/media/wang/BON/Projects/CP01/SEEG_LectureVWFA/sourcedata/eprime_logs'
wdir <- file.path(mdir,subj)
frec <- sprintf('%s_ses-01_task-RS_run-01_events_recovered.csv',subj)
# setup parameters
blocks <- c(1)  # sub-04: c(2,3,4,5); sub-06: c(1,3,4,5)
delay_limit <- 0.8    # in seconds
## ---------------------------

## recover trigger for each block
# initialize data frames
delays_recorded <- data.frame(conditions=character(),delays=double())
triggers_all <- data.frame(conditions=character(),onsets=double())
for (i in 1:length(blocks)){
  ib <- blocks[i]
  # clean recorded triggers
  triggers <- read.csv(file.path(wdir,sprintf('%s_ses-01_task-RS_run-01_block-%02d_events.csv',subj,ib)),
                       header=FALSE,stringsAsFactors=FALSE)
  names(triggers) <- c('conditions','onsets','durations')
  triggers <- triggers[!triggers$conditions %in% c('press'),]  # remove "press"
  triggers_working <- triggers  # initialize working triggers
  # read edat events
  edat <- read.table(file=file.path(wdir,sprintf('%s_ses-01_task-RS_run-01_block-%02d_edat.txt',subj,ib)),
                     header=TRUE,fileEncoding='UTF-16',sep='\t',stringsAsFactors=FALSE) 
  # change conditions names
  edat$conditions <- paste0(edat$Condition,substring(edat$RSE,1,1))
  edat$conditions <- gsub('catch[AV]','catch',edat$conditions)
  # calculate absolute timing for edat
  edat[is.na(edat)] <- 0
  edat$onsets_recorded <- edat$StimAud.OnsetTime + edat$StimVDO.OnsetTime + edat$StimVis.OnsetTime
  edat$onsets_absolute <- (edat$onsets_recorded - edat$onsets_recorded[1])/1000  # start with 0, in seconds
  # align absolute timing 
  trigger_first <- triggers_working$conditions[1]
  edat_align_first <- which(edat$conditions==trigger_first)[1]  # align the first recorded trigger
  triggers_working$onsets_absolute <- (triggers_working$onsets - triggers_working$onsets[1]) + edat$onsets_absolute[edat_align_first]
  # recover missing triggers
  ntrials <- dim(edat)[1]
  edat$delays <- rep(0,ntrials)
  edat$onsets_trigger <- rep(0,ntrials)
  edat$recorded_trigger <- rep(0,ntrials) 
  edat$recovered_trigger <- rep(0,ntrials)
  for (j in 1:ntrials){
    jcond <- edat$conditions[j]
    jtime <- edat$onsets_absolute[j]  # E-prime timing
    jtrig <- which(triggers_working$conditions == jcond)[1]
    if (!is.na(jtrig)){
      jdelay <- triggers_working$onsets_absolute[jtrig] - jtime
      if (abs(jdelay) < delay_limit){
        # correct E-prime timings according to this recorded trigger
        edat$delays[j] <- jdelay*1000  # in ms
        #edat$onsets_absolute[j:ntrials] <- edat$onsets_absolute[j:ntrials] + jdelay
        edat$onsets_trigger[j] <- triggers_working$onsets[jtrig]
        edat$recorded_trigger[j] <- 1
        triggers_working <- tail(triggers_working,-jtrig)
      } else {
        # no recorded trigger
        edat$recovered_trigger[j] <- 1
      }
    } else {
      # no recorded trigger
      edat$recovered_trigger[j] <- 1
    }
  }
  for (j in 1:ntrials){
    if (edat$recovered_trigger[j] == 1){
      jpos <- which(edat$recovered_trigger==0)-j
      jref <- jpos[which.min(abs(jpos))]+j
      edat$onsets_trigger[j] <- edat$onsets_trigger[jref]-edat$onsets_absolute[jref]+edat$onsets_absolute[j]
    }
  }
  # report the sum and mean delays
  delay_max <- max(abs(edat$delays))
  delay_sum <- sum(abs(edat$delays))
  delay_avg <- delay_sum/sum(edat$onsets_trigger != 0)
  cat(sprintf('\nNr1 / Nr2 = %d / %d. The max, mean and summed delays of Block %g are : %f ms, %f ms, and %f ms.\n',
              dim(triggers)[1], sum(edat$recorded_trigger), ib, delay_max, delay_avg, delay_sum))
  delays_recorded <- rbind(delays_recorded,edat[edat$onsets_trigger != 0,c('conditions','delays')])
  # combine triggers of all runs
  triggers_all <- rbind(triggers_all,edat[,c('conditions','onsets_trigger')])
  # output the information of recovered timings
  write.table(edat,file=file.path(wdir,sprintf('%s_ses-01_task-RS_run-01_block-%02d_events_recovered-info.csv',subj,ib)),
              row.names=FALSE,quote=FALSE,sep=',')
}
# output the complete timings for BST
#triggers_all$durations <- rep(0,dim(triggers_all)[1])
write.table(triggers_all,file=file.path(wdir,frec),col.names=FALSE,row.names=FALSE,quote=FALSE,sep=',')
write.table(delays_recorded,file=file.path(wdir,sprintf('%s_ses-01_task-RS_run-01_events_recovered-delays.csv',subj,ib)),
            row.names=FALSE,quote=FALSE,sep=',')
## ---------------------------