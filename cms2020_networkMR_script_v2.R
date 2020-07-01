########################
########################
# PREPARING R ENVIRONMENT
########################
########################

########################
# Install packages (first time only)
########################

#Run this block seperately first before running the rest of the script
#install.packages("devtools")
#install.packages("usethis")
#library(devtools)

#Install packages - uncomment and run the first time you run the script, or to update packages
# install_github("MRCIEU/TwoSampleMR")
# install_github('qingyuanzhao/mr.raps')
# install_github("WSpiller/RadialMR")
# install_github("explodecomputer/tryx")
# install_github("explodecomputer/simulateGP")

########################
# Load packages
########################

# rm(list=ls(all=TRUE))
# library(devtools)
# library(TwoSampleMR)
# library(ggplot2)
# library(knitr)
# library(mr.raps)
# library(MRInstruments)
# library(readr)
# library(dplyr)
# library(tidyr)
# library(tryx)

########################
# Initialise variables
########################

# #Build a list of tests already conducted list / variable for exposures/outcome 
# #data/harmonisation already obtained/conducted in previous runtimes of the script 
# #to minimise restart time (OPTIONAL)
# 
# if (exists("exposure_dat")){
#   if(!is.null(exposure_dat)){
#     if(nrow(exposure_dat)>0){
#       dataExistsFor_exposures <- unique(exposure_dat$id.exposure)
#     }else{dataExistsFor_exposures <- FALSE}
#   }else{dataExistsFor_exposures <- FALSE}
# }else{dataExistsFor_exposures <- FALSE}
# 
# if (exists("outcome_dat")){
#   if(!is.null(outcome_dat)){
#     if(nrow(outcome_dat)>0){
#       dataExistsFor_outcomes <- TRUE
#     }else{dataExistsFor_outcomes <- FALSE}
#   }else{dataExistsFor_outcomes <- FALSE}
# }else{dataExistsFor_outcomes <- FALSE}
# 
# if (exists("dat")){
#   if(!is.null(dat)){
#     if(nrow(dat)>0){
#       dataExistsFor_harmonisation <- TRUE
#     }else{dataExistsFor_harmonisation <- FALSE}
#   }else{dataExistsFor_harmonisation <- FALSE}
# }else{dataExistsFor_harmonisation <- FALSE}
# 
# if (exists("res_all")){
#   if(!is.null(res_all)){
#     if(nrow(res_all)>0){
#       dataExistsFor_MR <- unique(paste(res_all$id.exposure,res_all$id.outcome))
#     }else{
#       dataExistsFor_MR <- FALSE
#       res_all <- data.frame(id.exposure=character(),
#                             id.outcome=character(),
#                             outcome=character(),
#                             exposure=character(),
#                             method=character(),
#                             nsnp=character(),
#                             b=double(),
#                             se=double(),
#                             pval=double(),
#                             analysis_type=character(),
#                             stringsAsFactors=FALSE)
#       }
#   }else{
#     dataExistsFor_MR <- FALSE
#     res_all <- data.frame(id.exposure=character(),
#                           id.outcome=character(),
#                           outcome=character(),
#                           exposure=character(),
#                           method=character(),
#                           nsnp=character(),
#                           b=double(),
#                           se=double(),
#                           pval=double(),
#                           analysis_type=character(),
#                           stringsAsFactors=FALSE)
# }else{
# dataExistsFor_MR <- FALSE
# res_all <- data.frame(id.exposure=character(),
#                       id.outcome=character(),
#                       outcome=character(),
#                       exposure=character(),
#                       method=character(),
#                       nsnp=character(),
#                       b=double(),
#                       se=double(),
#                       pval=double(),
#                       analysis_type=character(),
#                       stringsAsFactors=FALSE)
# }
# 
#These variables need to be initialised first by assigning them NULL, so that they later can be used in iteration without a %variable name% error. Also helps ensure re-running this doesn't duplicate your data.

counter<- -1
res_descriptives<- NULL
res_Q <- NULL
res_tryx <- NULL
res_main_snps <- NULL
res_main_leaveOneOut <- NULL
exposure_dat<-NULL
GWAS_SNPs_inAnalysis <- NULL
res_all <- data.frame(id.exposure=character(),
                      id.outcome=character(),
                      outcome=character(),
                      exposure=character(),
                      method=character(),
                      nsnp=character(),
                      b=double(),
                      se=double(),
                      pval=double(),
                      analysis_type=character(),
                      stringsAsFactors=FALSE)
# 
# ########################
# # Set working directory
# ########################
# #CHANGE TO YOUR DIRCTORY
# setwd('C:/Users/Chris Moreno-Stokoe/OneDrive/Documents/Research/Wellbeing and Insomnia MR Study/data/confirmatory mr')
# 
# ########################
# # Read-in trait name and SNP lists
# ########################
# 
# #import csv file with betterNames trait name masterlist and subset to only those used in the present MR
# #Specify trait IDs here
# betterNames <- read_csv("C:/Users/Chris Moreno-Stokoe/OneDrive/Documents/Research/Wellbeing and Insomnia MR Study/data/confirmatory mr/betterNames.csv")
# betterNames <- subset(betterNames, (betterNames$`In study?` != 0 & betterNames$`In study?` != 'snp reference'))
# 
# #import csv file with betterSNPs snp masterlist from published GWASes of traits in analysis for more reliable SNPs in MR
# betterSNPs <- read_csv("C:/Users/Chris Moreno-Stokoe/OneDrive/Documents/Research/Wellbeing and Insomnia MR Study/data/confirmatory mr/betterSNPs.csv")
# 
# #build list of name:Id conversions
# exposure_names <- cbind('id.exposure' = betterNames$`New Id`, 'exp_name' = betterNames$`Custom name`)
# outcome_names <- cbind('id.outcome' = betterNames$`New Id`, 'out_name' = betterNames$`Custom name`)
# 
# 
# ########################
# ########################
# # Extract instruments
# ########################
# ########################
# 
# 
# #For each trait in analysis, extract instruments according to properties
# for (trait in betterNames$`New Id`){
# 
#   #IF trait has pre-specified SNPs:
#   if (sum(!is.na(select(betterSNPs, trait)))>0){
#     #Conduct liberal GWAS to ensure we identify desired SNPs and collect their info & exposure info
#     exposure_dat_chunk <- extract_instruments(outcomes=trait,
#                                                      force_server = FALSE,
#                                                      p1 = 5e-04,
#                                                      clump = FALSE,
#                                                      p2 = 5e-04,
#                                                      r2 = 0.6,
#                                                      kb = 250)
# 
#     #IF any SNPs found
#     if (!is.null(exposure_dat_chunk)){
# 
#       #Print dataExistsFor_MR message
#       print(paste(nrow(exposure_dat_chunk), ' SNPs found for ', trait),sep="")
# 
#       #Restrict discovered SNPs only to select snps for traits where only specific SNPs are desired
#       selectSNPs <- subset(exposure_dat_chunk, SNP %in% pull(select(betterSNPs, trait)))
#       exposure_dat <- rbind(exposure_dat, selectSNPs)
#       print(paste(nrow(selectSNPs), '/', sum(!is.na(select(betterSNPs, trait))), ' pre-specified SNPs found for ', trait),sep="")
# 
#     }
#   }
#   #ELSE trait does not:
#   else {
#     #Conduct a novel GWAS with conservative values
#     exposure_dat_chunk <- extract_instruments(outcomes=trait,
#                                               force_server = TRUE,
#                                               p1 = 5e-08,
#                                               clump = TRUE,
#                                               p2 = 5e-08,
#                                               r2 = 0.001,
#                                               kb = 10000)
#     exposure_dat <- rbind(exposure_dat, exposure_dat_chunk)
# 
#     #Print dataExistsFor_MR message
#     print(paste(nrow(exposure_dat_chunk), ' SNPs found for ', trait, sep=""))
#   }
# 
# }
# 
# ########################
# # Obtain outcome data
# ########################
# 
# #Obtain outcome data for all SNPs and traits in analysis together
# outcome_dat <- extract_outcome_data(exposure_dat$SNP, betterNames$`New Id`,
#                                     proxies = 1,
#                                     rsq = 0.8,
#                                     align_alleles = 1,
#                                     palindromes = 0,
#                                     maf_threshold = 0.3
# )
# 
# ########################
# # Harmonisation
# ########################
# 
# #Harmonise all data together
# dat <- harmonise_data(
#   exposure_dat = exposure_dat,
#   outcome_dat = outcome_dat,
#   action = 2
# )
# 
# ########################
# # Steiger filtering
# ########################
# 
# #Identify SNPs better associated with the outcome than exposure and add this information to the harmonised dat
# dat_steiger <- steiger_filtering(dat)
# steiger <- data.frame(SNP=dat_steiger$SNP,
#                       exposure=dat_steiger$exposure,
#                       outcome=dat_steiger$outcome,
#                       rsq.exposure=dat_steiger$rsq.exposure,
#                       rsq.outcome=dat_steiger$rsq.outcome,
#                       steiger_dir=dat_steiger$steiger_dir,
#                       steiger_pval=dat_steiger$steiger_pval)
# dat <- merge(dat, steiger, by=c('SNP','outcome','exposure'))
# 
# ########################
# # Tidy data
# ########################
# 
# #Exclude intra-trait data (e.g., BMI with BMI)
# dat <- subset(dat,id.exposure != id.outcome)
# 
# #Give proper working names to traits in analysis (result: no csv-breaking commas, and no empty outcome names)
# dat<-merge(dat,exposure_names,by='id.exposure')
# dat<-merge(dat,outcome_names,by='id.outcome')
# dat$exposure <- dat$exp_name
# dat$outcome <- dat$out_name
# dat <- subset(dat, select=-c(exp_name, out_name))
# 
# #Exclude analyses of the same trait to itself (e.g., BMI to BMI)
# dat <- subset(dat, id.exposure != id.outcome)
# 
# #Ensure strings are NOT factors
# dat <- data.frame(dat, stringsAsFactors=FALSE)
# 

########################
########################
# Conduct MR 
########################
########################

#Run MR & sensitivity tests then save to dataframes

# For each exposure and outcome pair
for (expId in betterNames$`New Id`){
  for (outId in betterNames$`New Id`){
      
      #If running a test within the same trait (e.g., bmi --> bmi) OR if tests already conducted then skip
      if(expId==outId){
        next
      }
    
      ########################
      # Subset data
      ########################
    
      #Subset valid SNPs by exposure, outcome
      datFinal <- subset(dat, id.exposure == expId & id.outcome == outId & mr_keep == TRUE)
      datFinal <- data.frame(unique(datFinal, by = "SNP"), stringsAsFactors=FALSE)
      
      datFinal <- data.frame(unique(datFinal, by = "SNP"), stringsAsFactors=FALSE)
      
      #If no valid SNPs were found, skip analysis
      if(nrow(datFinal) < 1){
        print(paste('Processed. No SNPs available: ', expId, ' on ', outId, sep=""))
        next
      } else {
        print(paste('Processed. Analysing (p08): ', expId, ' on ', outId, sep=""))
      }
      
      ########################
      # Descriptive statistics
      ########################
   
      #Number of SNPs in analysis
      nsnp <- nrow(datFinal)
      
      #Number of SNPs which better predict the outcome than the exposure
      steigerCount <- nrow(subset(datFinal,(steiger_dir == FALSE)))
      
      #Mean, minimum and maximum signal to noise ratios for SNPs in present analysis
      datFinal$fStat <- ((datFinal$beta.exposure)^2/(datFinal$se.exposure)^2)
      f <- c(mean(datFinal$fStat),min(datFinal$fStat),max(datFinal$fStat))
      
      #Classify this analysis for ease of sorting (OPTIONAL)
      
        #Get analysis type from betterNames dataframe
        expType <- subset(betterNames, `New Id` == expId)$`In study?`
        outType <- subset(betterNames, `New Id` == outId)$`In study?`
        
        #Classify by matching analysis types
        if (expType == 'Main' && outType == 'Main'){
          analysis_type <- 'Main'
        } else if (expType == 'Replication' || outType == 'Replication') {
          if (expType == 'Network' || outType == 'Network'){
            analysis_type <- 'None'
          } else {
            analysis_type <- 'Replication'
          }
        } else {
          analysis_type <- 'Network'
        }
      
      ########################
      # Analysis
      ########################
      
      #IVW, Egger, Median & Mode methods
      res <- data.frame(mr(datFinal), stringsAsFactors=FALSE)
        
      #Steiger filtering test
      res_steiger <- mr_steiger2(datFinal$rsq.exposure, datFinal$rsq.outcome, datFinal$samplesize.exposure, datFinal$samplesize.outcome)
        
      #If there are enough SNPs conduct more tests (Egger intercept, single SNP* & MR Tryx*) *=if main/replication
      if(nsnp>1){
        
        #Perform & save egger intercept
        egg <- mr_egger_regression(datFinal$beta.outcome, datFinal$beta.exposure, datFinal$se.exposure, datFinal$se.outcome)
        res_egg <- cbind(id.exposure = datFinal$id.exposure[1], 
                         id.outcome=datFinal$id.outcome[1], 
                         outcome=as.character(datFinal$outcome[1]),
                         exposure=as.character(datFinal$exposure[1]), 
                         method = 'MR Egger intercept', 
                         nsnp=egg$nsnp, #AS NUMBER? OR COMBINE WITH RES AS BEFORE res<-(rbind(res,res_egg))
                         b=egg$b_i,
                         se=egg$se_i,
                         pval=egg$pval_i,
                         analysis_type=analysis_type)
        res_egg <- data.frame(res_egg, stringsAsFactors=FALSE)
        res_all <- rbind(res_all, res_egg)
        
        #Perform & save Cochrane's Q hereogeneity test
        Q <- subset(mr_heterogeneity(datFinal),method!='MR Egger')
        
        #Measurement error I2 regression coeffecient
        i2gx <- Isq(datFinal$beta.exposure, datFinal$se.exposure)
        
        #Test whether sensitivity tests agree in direction with the main IVW/wald estimates
          
          #Get IVW/wald estimate from main analysis
          res_ivw <- subset(res, method == "Inverse variance weighted" | method == "Wald ratio")
          
          #If it is IVW, then convert estimate from a string to number so we can operate on it
          if(res_ivw$method == "Inverse variance weighted"){
            res_ivw_pval=as.numeric(as.character(res_ivw$pval))
          } else {
            res_ivw_pval = res_ivw$pval
          }
          
          #Get sensitivity estimates from main analysis
          res_sensitivity <- subset(res, method != "Inverse variance weighted" & method != "Wald ratio" & method != "MR Egger intercept")

          #Test whether they agree in valence (i.e., if ivw is positive test if sensitivity tests are also positive)
          if(res_ivw$b >= 0){
            res_Disagreement <- sum(res_sensitivity$b < 0)
          }else if(res_ivw$b < 0){
            res_Disagreement <- sum(res_sensitivity$b > 0)
          }
          
      } else {
          Q<-data.frame(Q=NA,Q_df=NA,Q_pval=NA)
          res_Disagreement<-NA
          i2gx <- NA
      }
    
      #If this is the main or replication analyses: MR Tryx and make plots
      if(analysis_type=='Main' || analysis_type=='Replication'){
        
        print('Main/replication analysis. Making plots and performing TRYX')
        
        #Perform & save single SNP analysis
        res_singleSNP <- mr_singlesnp(datFinal)
        res_main_snps <- rbind(res_main_snps, res_singleSNP)
        
        #Perform & save leave one out analysis
        res_leaveOneOut <- mr_leaveoneout(datFinal)
        res_main_leaveOneOut <- rbind(res_main_leaveOneOut, res_leaveOneOut)
        
        #Plot and save plots
          #Make unique filename for each of the plots
          filename=paste(expId,'-',outId,sep="")
          
        #try plots
          #Scatter plot
          try(p1 <- mr_scatter_plot(res, datFinal))
          try(ggsave(p1[[1]], file=paste('cms2020_NetworkMR_Fig_scatter_',filename,'.png',sep=""), width=7, height=7))
          #Forest plot
          try(p2 <- mr_forest_plot(res_singleSNP))
          try(ggsave(p2[[1]], file=paste('cms2020_NetworkMR_Fig_forest_',filename,'.png',sep=""), width=7, height=7))
          #Funnel plot
          try(p3 <- mr_funnel_plot(res_singleSNP))
          try(ggsave(p3[[1]], file=paste('cms2020_NetworkMR_Fig_funnel_',filename,'.png',sep=""), width=7, height=7))
          #Leave-one-out plot
          try(p4 <- mr_leaveoneout_plot(res_leaveOneOut))
          try(ggsave(p4[[1]], file=paste('cms2020_NetworkMR_Fig_leaveOneOut_',filename,'.png',sep=""), width=7, height=7))
          
        #' #MR Tryx
        #' tryxscan <- NULL
        #' tryx.scan <- NULL
        #' tryxanalysis <- NULL
        #' tryx.network <- NULL
        #' 
        #' #'Try' running Radial MR outlier detection so if there is an error (e.g., finding nothing) it doesn't stop the script
        #' try(tryxscan <- tryx.scan(datFinal))
        #' 
        #' #If outliers are found run Tryx analysis
        #' if (!is.null(tryxscan)){
        #'   
        #'   #Select significant candidate trait(s)
        #'   tryxscan <- tryx.sig(tryxscan)
        #'   
        #'   #Plot results
        #'   tryx.network(tryxscan)
        #'   volcano_plot(
        #'     rbind(tryxscan$candidate_exposure_mr, tryxscan$candidate_outcome_mr)
        #'   )
        #'   volcano_plot(tryxscan$candidate_exposure_mr)
        #'   
        #'   #Adjust SNP effects on exposure and outcome given traits
        #'   tryxanalysis <- tryx.analyse(tryxscan)
        #'   tryxanalysis$estimates
        #'   tryxanalysis$plot
        #'   
        #'   #Save plot of tryx results under filename using exposure and outcome names
        #'   ggsave(tryxanalysis$plot, file=paste('cms2020_NetworkMR_Fig_MRTryx_',filename,'.png',sep=""), width=7, height=7)
        #' 
        #'   #Merge tryx results with main results
        #'   res_tryx <- rbind(res_tryx, tryxanalysis$estimates)
          
      }
      
      ########################
      # Save results
      ########################
      
      #Collect descriptives
      descriptives<-cbind(id.exposure = expId,
                          id.outcome = outId,
                          nsnp = nsnp,
                          analysis_type = analysis_type,
                          steiger_nsnpInvalid = steigerCount,
                          steiger_pval = res_steiger$steiger_test,
                          steiger_r2out = res_steiger$r2_out,
                          steiger_r2exp = res_steiger$r2_exp,
                          f_mean = f[1],
                          f_min = f[2],
                          f_max = f[3],
                          nSTestsDisagree = res_Disagreement,
                          I2GX = i2gx,
                          Q = Q$Q,
                          Q_df = Q$Q_df,
                          Q_pval = Q$Q_pval)
      
      #Merge and store results
        #Column bind analysis type for easy sorting
        res$analysis_type <- analysis_type
        datFinal$analysis_type <- analysis_type
        
        #Main results
        res_all <- rbind(res_all, res)
        #Descriptives
        res_descriptives <- rbind(res_descriptives, descriptives)
        #SNPs
        GWAS_SNPs_inAnalysis <- rbind(GWAS_SNPs_inAnalysis, datFinal)
  
  }
  
  #Script dataExistsFor_MR indicator
  counter <- counter+1
  print('########################################################')
  print(paste('Progress: ',(counter / length(betterNames$`New Id`)*100),'%'),sep="")
  print('########################################################')
  
}

########################
########################
# SAVING OUTPUTS
########################
########################

########################
# Make summary statistics
########################
Total_tests <- nrow(res_all)
Total_uniqueSNPs <- nrow(unique(GWAS_SNPs_inAnalysis$SNP))

resultsTable<-data.frame(res_all)
IVWs <- subset(resultsTable, method == "Inverse variance weighted" || method == "Wald ratio")

Total_tests_main <- nrow(subset(IVWs,analysis_type=='Main'))
Total_tests_replication <- nrow(subset(IVWs,analysis_type=='Replication'))
Total_tests_network <- nrow(subset(IVWs,analysis_type=='Network'))



# #Linkage disequilibirum estimates
# 
#   #Get exposure and outcome SNPs
#   snps <- c(select(betterSNPs, expId),select(betterSNPs, outId))
#   b<-select(betterSNPs, expId)
#   c<-select(betterSNPs, outId)
#   d<-c(b,c)
#   a<-as.character(snps)
#   a[a!='NA']
#   
#   betterSNPs[betterSNPs!is.na]
#   e<-NULL
#   f<-betterSNPs[!is.na]
#   
#   for (snp in betterSNPs){
#     e <- c(e,snp)
#   }
#   
#   #Make LD matrix
#   LDmatrix <- ieugwasr::ld_matrix(betterSNPs) #FIX ME! snps[snps!is.na]
#   
#   #Take absolute R values and exclude perfect self correlations
#   LDs <- abs(LDmatrix[LDmatrix!=1])
#   
#   #Make mean, min, max values
#   LD <- c(mean(LDs),min(LDs),max(LDs))

########################
# Save instrument data
########################

#Save the raw SNPs found in the exposure and outcome data (for reference on SNPs unavailable for present analysis)
write.csv(exposure_dat, "cms2020_NetworkMR_GWAS_exposures.csv", row.names=F, quote=F)
write.csv(outcome_dat, "cms2020_NetworkMR_GWAS_outcomes.csv", row.names=F, quote=F)

#Save the SNPs remaining after harmonisation and outcome data search
write.csv(dat, "cms2020_NetworkMR_GWAS_all.csv", row.names=F, quote=F)

#Save SNPs used in each analysis
write.csv(GWAS_SNPs_inAnalysis, "cms2020_NetworkMR_SNPs_inAnalysis.csv", row.names=F, quote=F)

########################
# Save analysis data
########################

#Save results
write.csv(res_all, "cms2020_NetworkMR_results_MR_all.csv", row.names=F, quote=F)
write.csv(subset(res_all,analysis_type=='Main'), "cms2020_NetworkMR_results_MR_main.csv", row.names=F, quote=F)
write.csv(subset(res_all,analysis_type=='Main'||analysis_type=='Replication'), "cms2020_NetworkMR_results_MR_replication.csv", row.names=F, quote=F)
write.csv(subset(res_all,analysis_type=='Main'||analysis_type=='Network'), "cms2020_NetworkMR_results_MR_network.csv", row.names=F, quote=F)

#Save descriptives
write.csv(res_descriptives, "cms2020_NetworkMR_results_descriptives.csv", row.names=F, quote=F)
write.csv(subset(res_descriptives,analysis_type=='Main'||analysis_type=='Replication'), "cms2020_NetworkMR_results_descriptives_replication.csv", row.names=F, quote=F)

#(Due to their different format not all analyses merge into the main results dataframe)

#Save singleSnp
write.csv(res_main_snps, "cms2020_NetworkMR_results_singleSnpResults.csv", row.names=F, quote=F)

#Save leaveOneOut
write.csv(res_main_leaveOneOut, "cms2020_NetworkMR_results_leaveOneOut.csv", row.names=F, quote=F)

#Save mrTryx
#write.csv(res_Tryx, "cms2020_NetworkMR_results_tryx.csv", row.names=F, quote=F)