## ---------------------------
## [script name] ps02_DATA_examine_trigger_delays.R
##
## SCRIPT to examine the timing differences between recorded triggers and E-prime logs, and to validate the temporal
##           precision of recovered triggers.
##
## By Shuai Wang, [date] 2021-03-01
##
## ---------------------------
## Notes:
##   
##
## ---------------------------

## clean up
rm(list=ls())
## ---------------------------

## set environment (packages, functions, working path etc.)
# load up packages
library(ggplot2)
library(cowplot)
# setup working path
mdir <- '/media/wang/BON/Projects/CP01/SEEG_LectureVWFA/sourcedata/eprime_logs'
# define parameters
subjects <- c('sub-04','sub-05','sub-06','sub-07')
settings <- c('bedside','lab','bedside','lab')
isRemove0s <- TRUE
## ---------------------------

## analyze timing delays for each subject
timings_all <- data.frame(subject=character(),conditions=character(),delays=double(),settings=character())
for (i in 1:length(subjects)){
  subj <- subjects[i]
  sdir <- file.path(mdir,subj)
  fdel <- file.path(sdir,sprintf('%s_ses-01_task-RS_run-01_events_recovered-delays.csv',subj))
  ffig <- file.path(mdir,sprintf('%s_ses-01_task-RS_run-01_events_recovered-delays.png',subj))
  timings <- read.csv(file=fdel)
  if (isRemove0s){
    timings <- timings[timings$delays > 1 | timings$delays < -1,]  # only remain timings with delay greater than 1 ms
  }
  plist <- lapply(unique(timings$conditions), 
               function(x) qplot(timings$delays[timings$conditions==x],geom="histogram",binwidth=1,xlab=x,ylim=c(0,10)))
  ps <- plot_grid(plotlist=plist)
  save_plot(filename=ffig,ps,base_height=5)
  # combine subjects
  timings$subject <- rep(subj,dim(timings)[1])
  timings$settings <- rep(settings[i],dim(timings)[1])
  timings_all <- rbind(timings_all,timings)
}
## ---------------------------

## analyze timing delays for all subjects
# by subjects
plist_subject <- lapply(unique(timings_all$subject), 
                           function(x) qplot(timings_all$delays[timings_all$subject==x],geom="histogram",binwidth=1,xlab=x))
ps_subject <- plot_grid(plotlist=plist_subject)
save_plot(filename=file.path(mdir,'all-subject_task-RS_run-01_events_recovered-delays.png'),ps_subject,base_height=5)
# by conditions
plist_conditions <- lapply(unique(timings_all$conditions), 
                           function(x) qplot(timings_all$delays[timings_all$conditions==x],geom="histogram",binwidth=1,xlab=x,ylim=c(0,30)))
ps_conditions <- plot_grid(plotlist=plist_conditions)
save_plot(filename=file.path(mdir,'all-conditions_task-RS_run-01_events_recovered-delays.png'),ps_conditions,base_height=5)
# by settings
plist_settings <- lapply(unique(timings_all$settings), 
                           function(x) qplot(timings_all$delays[timings_all$settings==x],geom="histogram",binwidth=1,xlab=x))
ps_settings <- plot_grid(plotlist=plist_settings)
save_plot(filename=file.path(mdir,'all-settings_task-RS_run-01_events_recovered-delays.png'),ps_settings,base_height=5)
## ---------------------------