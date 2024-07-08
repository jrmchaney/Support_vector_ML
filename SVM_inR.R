#this script runs SVM on badaga FFR data
#trains on a subset of quiet trials, then tests on the remaining quiet trials
#written by Jacie McHaney 8/24/22s

rm(list = ls())

library(tidyverse)
library(e1071)
library(foreach)
library(doParallel)
library(tictoc)
library(rstatix)

####----set up cluster for parallel processing----
cl <- makeCluster(future::availableCores()) #must use future package for this to work properly on Pitt CRC
registerDoParallel(cl)

#load libraries inside the cluster
clusterEvalQ(cl, {
  library(tidyverse)
  library(e1071)
})

#####----start of script----
#paths
maindir <- '/bgfs/bchandrasekaran/jbm90/NSAA'
inpath <- paste(maindir,'/FFR_data/processed_FFR/spin/FFR/sub_svm',sep='') #location of dat files
outpath <- paste(maindir, '/Radialsvm_NSAA_0to70_trainQ_testQ',sep='') #where to store svm info (BE SPECIFIC WITH THIS FOLDER BECAUSE YOU WILL RUN MANY DIFFERENT ITERATIONS)
ifelse(!dir.exists(file.path(outpath)), #automatically make the above outpath if doesn't already exist
       dir.create(file.path(outpath)), FALSE) #will output TRUE to console if it didn't exist and was created successfully

#list of dat files
files <- dir(inpath,'.dat')

#read single file
#data <- read.table(paste(inpath,'/',files[1],sep  = ''), sep = '\t')


# parameters --------------------------------------------------------------

stims <- c('BA','DA','GA') #change this to the OG, HP, LP, ILL conditions
conds <- c('quiet','SSN') #change to laurel / variable distinction

spoint <- 1000/25000 #how many milliseconds per sampling point
prestim <- 20 #how many milliseconds prestim is there (check with Anoop)
ipstart <- (prestim + 0)/spoint #where to start svm analysis. currently set to sound onset
ipend <- (prestim + 70)/spoint #where to end analysis.  currently set to length of stimulus. Change 70 to match length of your stimulus in ms
iprange <- seq(ipstart,ipend) #range of time (in sampling points) to run SVM on
iprange.char <- as.character(foreach(i = 1:length(iprange)) %do% paste('sp',iprange[i],sep = '')) #list of column names to trim to 

allresp <- data.frame()
lenFiles <- seq(1,length(files)) #vector of file numbers
outacc <- data.frame()
today <- Sys.Date()

#parameters for averaging across trials. Single trial FFRs are often too noisy. Here, it is set to average across 25 trials for a new "single" trial
bin <- 25 #number of trials to average across
#create vector of bin numbers
binlab <- 1000/bin #here, change 1000 to match the total number of trials each subject has per listening stimulus
binvec1 <- seq(1,binlab)
binvec <-  rep(binvec1,bin)

###----send variables to cluster----
vars <- c('allresp','stims','conds','iprange','files','inpath','lenFiles','outpath','outacc','bin','binvec')
clusterExport(cl,vars)

##----define function to compile all info for each subject from dat files----
#Ja Young will need to modify the below function to call in the HUGE .mat file to get a final dataframe with the following layout:
#subject,group assignment,stimulus, trial number, then ffr waveworm
compiledat <- function(i){
  
  fname <- files[i]
  fnow <- read.table(paste(inpath,'/',fname,sep  = ''), sep = '\t')
  
  #get file info
  nameparts <- strsplit(fname, split = '_')
  stim <- as.vector(nameparts[[1]][2])
  sub <- as.vector(nameparts[[1]][1])
  sub <- substring(sub,4)
  cond <- as.vector(nameparts[[1]][3])
  cond <- substring(cond,1,nchar(cond)-4)
  
  #extract subject data
  subdat <- fnow %>% 
    mutate(sub = sub,
           stim = stim,
           cond = cond,
           trial = seq(1,nrow(fnow))) %>%
    select(sub,stim,cond,trial,everything())
  
  #add together
  allresp <- as.data.frame(rbind(allresp,subdat))
  
  #make sure allresp is spit out from function
  return(allresp)
  
} #end compiledat function


##----now actually compile subject FFR data from dat files----
tic()
#allresp <- foreach(i = 1:length(files)) %dopar% compiledat(i,files,inpath)
allresp <- parLapply(cl,lenFiles,compiledat)
toc() #83s for 78files with 8 cores

#convert from weird list to a nice data frame
allresp.now <- data.frame()
tic()
for(i in 1:length(allresp)){
  
  allresp.now <- rbind(allresp.now,allresp[[i]])
  
}
toc() #183s for 78 files with 8 cores


##----prepare for SVM----
#make list of new colnames
oldnames <- colnames(allresp.now)[5:5754] #change these numbers to match your dataset
newnames <- as.character(foreach(i = 1:length(oldnames)) %do% paste('sp',substring(oldnames[i],2),sep = ''))

#1. trim all responses to just time of interest and reshape
# trim to 1000 trials per sub, stim, condition (I chose 1000 trials for NSAA subs. Change to fit your needs)
# group by bin size and find the average FFR per bin
#first remove any subjects who don't even have 1000 trials
bad.subs <- allresp.now %>% 
  group_by(sub,cond,stim) %>% 
  count() %>% 
  ungroup() %>% 
  subset(n < 1000) %>% 
  select(sub) %>% 
  distinct() %>% 
  pull(sub)

#now average remaining subjects across 25 trials (i call them bins)
allresp.svm <- allresp.now %>% 
  rename_at(vars(oldnames),~newnames) %>% 
  mutate(cond = case_when(cond == 'ssn' ~ 'SSN',
                          cond == 'Quiet' ~ 'quiet',
                          cond == 'SSNeeg'~ 'SSN',
                          TRUE ~ cond)) %>% 
  select(sub,stim,cond,trial,iprange.char)%>% 
  subset(!sub %in% bad.subs) %>% 
  group_by(sub,cond,stim) %>% 
  slice_head(n=1000) %>% 
  mutate(bin = binvec) %>% 
  select(sub:trial,bin,everything()) %>% 
  ungroup() %>% 
  group_by(sub,cond,stim,bin) %>% 
  summarise(across(iprange.char,mean))

sublist <- unique(allresp.svm$sub)

####need to trim to bins
# #1. how many trials does each sub have?
# sub.trials <- allresp.svm  %>% 
#   subset(grepl('quiet',cond, ignore.case = TRUE)) %>% 
#   group_by(sub,stim) %>% 
#   count() %>% 
#   ungroup() %>% 
#   arrange(n)
# min <- min(sub.trials$n)

#subtrain.out <- list()


##----define a function to run SVM----
vars2 <- c('allresp.svm','sublist')
clusterExport(cl,vars2)
#clusterExport(cl,'allresp.svm')

runsvm <- function(sub.now){ #instead of sub.now, change to stim.now
  
  idx <- allresp.svm %>% 
    subset(sub == sub.now) %>% #change to stim.now
    subset(grepl('quiet',cond,ignore.case = TRUE)) #i think you can delete this line
  
  n <- nrow(idx)
  
  if(n == 0){
    next
  }
  
  #get number of trials, in 75% for training and 25% for testing
  tot_trials <- nrow(idx)
  train.trials <- sample(1:tot_trials, round(tot_trials*.75), replace = FALSE) #generate random trial numbers for testing (75%)
  
  #split to training and testing sets
  train_set <- idx %>% 
    filter(row_number() %in% train.trials)
  test_set <- idx %>% 
    filter(!row_number() %in% train.trials)
  
  #without stim column
  train_setx <- subset(train_set, select = -c(sub,stim,cond,bin)) #ffr info only. 1 row per trial. quiet only
  train_sety <- as.factor(train_set$stim) #just list of stimuli
  test_setx <- subset(test_set, select = -c(sub,stim,cond,bin)) #ssn only. 1 row per trial
  test_sety <- as.factor(test_set$stim) #list of stimuli for testing
  
  #classification mode. predict stimulus (y) based on FFRs defined in x
  #tic()
  model <- svm(x = train_setx,
               y = train_sety,
               kernel = 'radial') #can also try linear kernel
  #toc() #1 sub takes ~ 8min
  
  #test model with testing data
  y_pred <- predict(model,newdata = test_setx)
  
  #make confusion matrix
  cm <- table(test_sety,y_pred) #rows = what it actually is, cols = what model classified as
  
  #save outputs
  #write.table(model, paste(outpath,'/sub',sub.now,'_svmModel.txt',sep = ''),row.names = FALSE)
  write.csv(y_pred, paste(outpath,'/sub',sub.now,'_y_pred_bin',bin,'.csv',sep = ''),row.names = FALSE)
  write.csv(test_setx, paste(outpath,'/sub',sub.now,'_testset_x_bin',bin,'.csv',sep = ''),row.names = FALSE)
  write.csv(cm, paste(outpath,'/sub',sub.now,'_confusionmatrix_bin',bin,'.csv',sep = ''))
  write.csv(test_sety, paste(outpath,'/sub',sub.now,'_testset_y_bin',bin,'.csv',sep = ''),row.names = FALSE)
  
  #calculate decoding accuracy
  
  
  
  return(outacc)
  
}

##----now actually run svm----
tic()
outacc <- parLapply(cl,sublist,runsvm)
toc()

#convert from weird list to a nice data frame
outacc.now <- data.frame()
tic()
for(i in 1:length(outacc)){
  
  outacc.now <- rbind(outacc.now,outacc[[i]])
  
}
toc()

write.csv(outacc.now,paste(maindir,'/svm_outacc','/NSAA_radial_trainQ_testQ_0to70_bin',bin,'_',today,'.csv',sep = ''),row.names = FALSE)

#shutdown the cluster
stopCluster(cl)
