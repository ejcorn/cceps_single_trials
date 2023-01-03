# make a table of patient characteristics

rm(list=ls())

basedir <- '/Users/Eli/Dropbox/CNTProjects/CCEPs_Projects/cceps_single_trials/'
source(paste0(basedir,'cceps_files.R'))
source(paste0(basedir,'trial_variability/miscfxns/packages.R'))
source(paste0(basedir,'trial_variability/miscfxns/miscfxns.R'))

# load redcap - ids + age + sex + Presurgical Clinical Info all fields
#redcap.df <- read.csv(file='~/Documents/Research/CNT/CCEPs/CNTSurgicalRepositor-EliCCEPs_DATA_2022-06-10_1519.csv',stringsAsFactors = F)
redcap.df <- read.csv(file='~/Documents/Research/CNT/CCEPs/CNTSurgicalRepositor-EliCCEPs_DATA_LABELS_2022-06-10_1603.csv',stringsAsFactors = F,check.names = F)
redcap.df <- redcap.df[redcap.df$`HUP Number` != '',]
rownames(redcap.df) <- paste0('HUP',redcap.df$`HUP Number`) # link row names to HUP ID

hup.ids <- locations$subjects # get hup IDs of patients I studied

# get med names from redcap indicator column names
get.med.names <- function(redcap.df,hup.ids,med.col.start='Current AED medications'){
  medcols <- colnames(redcap.df)[grep(med.col.start,colnames(redcap.df))]
  med.names <- sapply(strsplit(medcols,'='), function(X) strsplit(X[2],'/|)'))
  med.names <- sapply(med.names, function(X) X[1])
  hup.id.meds <- sapply(hup.ids, function(hup.id) med.names[which(redcap.df[hup.id,medcols]=='Checked')])
  hup.id.meds <- sapply(hup.id.meds, function(meds) paste(meds,collapse = ', '))
  return(hup.id.meds)
}

hup.id.current.meds <- get.med.names(redcap.df,hup.ids,'Current AED medications')
hup.id.past.meds <- get.med.names(redcap.df,hup.ids,'AED Trials')

# load SOZ localization file from Erin on google drive
soz <- read.csv(paste0(locations$data_folder,'SOZ.csv'),row.names = 1,stringsAsFactors = F)

pt.tbl <- data.frame(`HUP ID` = hup.ids,
           `Age (y)` = redcap.df[hup.ids,'Age at implant'],
           Sex = redcap.df[hup.ids,'Sex assigned at birth'],
           `Seizure Onset Zone` = tools::toTitleCase(soz[hup.ids,'SOZ.localization']),
           `AEDs during CCEPs` = hup.id.current.meds,
           `Home AEDs` = hup.id.past.meds,
           check.names = F,stringsAsFactors = F)

# manually correct some things
pt.tbl$`Epilepsy Localization`[is.na(pt.tbl$`Epilepsy Localization`)] <- 'Unavailable'
pt.tbl$`Home AEDs` <- gsub(pattern = 'LEV',replacement = 'Levetiracetam',x=pt.tbl$`Home AEDs`)
pt.tbl$`Home AEDs` <- gsub(pattern = 'OXC',replacement = 'Oxcarbazepine',x=pt.tbl$`Home AEDs`)
pt.tbl$`Home AEDs` <- gsub(pattern = 'ZNS',replacement = 'Zonisamide',x=pt.tbl$`Home AEDs`)
pt.tbl$`Home AEDs` <- gsub(pattern = 'TPM',replacement = 'Topiramate',x=pt.tbl$`Home AEDs`)
pt.tbl$`Home AEDs` <- gsub(pattern = 'LTG',replacement = 'Lamotrigine',x=pt.tbl$`Home AEDs`)
     
# Load erin patient data file
erin <- read.csv(paste0(locations$data_folder,'ErinPatientData.csv'),row.names = 1,stringsAsFactors = F)

if(identical(rownames(erin),rownames(pt.tbl))){
  pt.tbl$`Home AEDs` <- tools::toTitleCase(erin$Home.AEDs)
  pt.tbl$`AEDs during CCEPs` <- tools::toTitleCase(erin$Stim.AEDs)
  pt.tbl$`Seizures with Stim` <- erin$Seizure.with.stim.
}

write.csv(file = paste0(locations$results_folder,'PatientCharacteristicsCCEPs.csv'),x=pt.tbl)
library(xtable)
print(xtable(pt.tbl,caption='Patient characteristics',label='table:table1'),include.rownames=FALSE)

