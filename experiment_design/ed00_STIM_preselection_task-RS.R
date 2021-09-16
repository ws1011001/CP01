## SCRIPT to select words as stimuli for the SEEG repetition suppression experiment (CP01).
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
rs.sf.num <- 150
rs.sf.units <- 3
rs.vw.num <- 50

# read the v4 set of words
load('stimuli_selection/Lexique_Words_v4.Rdata')

# select sound files for the RS task
load('Words_wSoundFiles_Info.Rdata')
task.rs.sf <- words.info[order(words.info$v2length)[1:rs.sf.num],]
UnitComparisonP <- rep(0,length(VarMatched))
while (any(UnitComparisonP < 0.3)){
  # split the selected sound files into $ft.sf.units matched units
  task.rs.sf.units <- anticlustering(features=task.rs.sf[,VarMatched],K=rs.sf.units,objective='variance')
  # compare target variables across all units with Kruskal-Wallis Rank Sum Test
  UnitComparison <- lapply(task.rs.sf[,VarMatched], function(x) kruskal.test(x,as.factor(task.rs.sf.units)))
  UnitComparisonP <- unlist(lapply(UnitComparison, function(x) x$p.value))
}
print(paste('#####Quick Check - Units##### The Comparison results are: P =',UnitComparisonP))
task.rs.sf$units <- task.rs.sf.units

# complete the word selection for the RS task
task.rs.vw.pool <- words.all[words.all$ortho %!in% task.rs.sf$ortho,]  # remove the selected sound files
RSComparisonP <- rep(0,length(VarMatched))
while (any(RSComparisonP < 0.3)){
  randvw <- sample(c(rep(1,rs.vw.num),rep(0,dim(task.rs.vw.pool)[1]-rs.vw.num)))
  task.rs.pre <- rbind(task.rs.vw.pool[randvw!=0,VarMatched],task.rs.sf[,VarMatched])
  # compare target variables across all units with Kruskal-Wallis Rank Sum Test
  RSComparison <- lapply(task.rs.pre, function(x) kruskal.test(x,as.factor(c(rep(4,rs.vw.num),task.rs.sf$units))))
  RSComparisonP <- unlist(lapply(RSComparison, function(x) x$p.value))
}
print(paste('#####Quick Check - Units##### The Comparison results are: P =',RSComparisonP))
task.rs.vw <- task.rs.vw.pool[as.logical(randvw),]  # the selected visual words

# output the trials list and the words list for the FT task
task.rs.vw$units <- rep(4,rs.vw.num)
task.rs.vw$sfname <- task.rs.vw$ortho
task.rs <- rbind(task.rs.sf[,c(WordInfo,'units')],task.rs.vw[,c(WordInfo,'units')])
task.rs$condition <- c('AA','AV','VA','VV',task.rs$units)[match(task.rs$units,c(1:4,task.rs$units))]
write.csv(task.rs,file='Selected_Words_RepetitionSuppression_byConditions.csv',row.names=FALSE)