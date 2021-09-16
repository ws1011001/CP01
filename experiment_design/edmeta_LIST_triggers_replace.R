rm(list=ls())

for (i in 1:5){
  triggers <- read.table(sprintf('Triggers_RSE%d.txt',i))

  triggers.new <- triggers$V1
  
  # replace old triggers 
 
  triggers.new[triggers.new %in% c(11:60)] <- 41
  triggers.new[triggers.new %in% c(61:110)] <- 42
  triggers.new[triggers.new %in% c(111:160)] <- 43
  triggers.new[triggers.new %in% c(161:210)] <- 44
  triggers.new[triggers.new==255] <- 63
  triggers.new[triggers.new==220] <- 20
  triggers.new[triggers.new==230] <- 30
  triggers.new[triggers.new==240] <- 40
  triggers.new[triggers.new==250] <- 50
  triggers.new[triggers.new==1] <- 31
  triggers.new[triggers.new==2] <- 32
  triggers.new[triggers.new==3] <- 33
  triggers.new[triggers.new==4] <- 34
  triggers.new[triggers.new==5] <- 35
  triggers.new[triggers.new==6] <- 36
  triggers.new[triggers.new==7] <- 37
  triggers.new[triggers.new==8] <- 38
  write.table(triggers.new,file=sprintf('Triggers_v03_RSE%d.txt',i),row.names=FALSE,col.names=FALSE,quote=FALSE) 
}
