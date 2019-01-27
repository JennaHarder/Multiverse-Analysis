rm(list=ls())  #clear environment
#Load packages
library(ggplot2)
library(Rmisc)
library(lme4)
library(car)

#for reproducibility
set.seed(162)

#set working directory
setwd("~/QMES/multiverse analysis")

#Read in data
datafrm <- read.csv("alldata3.csv", stringsAsFactors = F, header = T) #dataframe compiling the data from a number of studies
save <- datafrm

#Optional: can reduce time this takes to run by limiting datafrm to the data from a single study (run Line #20). 
  #The script takes about 23 minutes to run if you do this (otherwise, may be more like 1-1.5 hours)
#datafrm <- datafrm[which(datafrm$study=="MastersS2"),]

#Some prep work####

#data structure: trials within participants

#We have two dependent variables:
#Reaction time (in milliseconds)
#Error (binary variable: did the participant make an error on this trial?)

#We have two categorical independent variables: 
#Object = Gun vs. Non-gun 
#Group (race) = Black vs. White.

#We want to effects-code our independent variables (1 vs -1) so that their interaction coefficients represent 
#the overall interaction rather than (e.g.) the change in the effect 
#of "Black" **when Object==Gun** (which is what it would be if we used dummy coding.)
#This will ensure that these interaction terms are telling us what we think they're telling us.
#To make sure R does this appropriately, we need to set these variables as factors (we'll
#do that later) and change the default for how R sets the contrasts on factor levels:
options("contrasts")
options(contrasts = c("contr.helmert", "contr.poly"))

#We'll now create a list of formulae (as string variables) that we'll plug into the  
#analyses in our multiverse program later on:

#fixed effects portions:
anteER <- "error_dummy ~ group + object_effects + group:object_effects " #for error analyses
anteRT <- "RT2 ~ group + object_effects + group:object_effects " #for RT analyses
anteRT2 <- "RT2 ~ group + object_effects " #for model comparisons to get p-values for RT analyses using the lme4 package

#participant random effects portions:
PrandomEffects <- c("", #no random effects for participant
                    "+ (1|s_pnum)", #random intercepts
                    "+ (1+object_effects|s_pnum)", #random slopes for Object
                    "+ (1+group|s_pnum)")  #random slopes for Race
#Notice that none of these allow slopes for both object and group to vary randomly by participant at the same time.
#This is because models with those random effects structures can be difficult for R to specify given the size of 
#these datasets: often don't converge.

#target random effects portions:
TrandomEffects <- c("", #no random effects for target
                    "+ (1|target)", #random intercepts
                    "+ (1+object_effects|target)") #random slopes for Object

#create vector of all the random effects formulas we'll use, 
#including every possible combination of participant and target random effects as listed above:
randomEffects <- vector()
for(p in(PrandomEffects)){
  for(t in(TrandomEffects)){
    randomEffects <- c(randomEffects, paste(p, t))
  }
}
randomEffects <- randomEffects[-1]  #the first item in this vector is blank ("" + ""), so get rid of that

#create vector with names of all the datasets compiled in our dataframe
datasets<-unique(datafrm$study)

#create dataframe with summarized error data, which we'll use to run error rate anova
anovadf <- data.frame(matrix(nrow=1, ncol=5)) #row of 5 NAs to initiate dataframe:  we'll remove this row later
names(anovadf) <- c("s_pnum", "group", "object_effects", "error_rate", "study") #name the columns
for(dset in(datasets)) { #we'll iterate through each dataset
  df <- datafrm[which(datafrm$study==dset),]  #isolate data for the present dataset
  
  rm(temp) #we'll create this dataframe in the next line of code, but since this is a loop, we must first get rid of the old one in iterations 2 through n
  
  #make temporary summary dataset with four rows per each of current dataset's participants:
  temp <- data.frame(matrix(nrow=length(unique(df$s_pnum))*4, ncol=5)) 
  names(temp) <- c("s_pnum", "group", "object_effects", "error_rate", "study")
  temp[,1] <- rep(unique(df$s_pnum), each = 4) #first column: participant numbers
  temp[,2] <- rep(c(-1, 1), each = 2) #second column:  effects-coded levels of race
  temp[,3] <- rep(c(-1, 1), each = 1) #third column:  effects-coded levels of object
  
  #calculate mean error rate for each race:object combination per each participant in dataset:
  for(s in(unique(df$s_pnum))){ #(for each participant)
    for(tr in(c(-1, 1))){ #(for each target race level)
      for(ob in(c(-1, 1))){   #(for each object level)
        temp[which(temp$s_pnum==s & temp$group==tr & temp$object_effects==ob), "error_rate"] <- mean(df[which(df$s_pnum==s & df$group==tr & df$object_effects==ob), "error_dummy"], na.rm=T)
      }
    }
  }
  
  #label which study this was:
  temp$study <- unique(df$study)
  
  #and add it to our larger dataframe (including all the datasets we're using)
  anovadf <- rbind(anovadf, temp)
}
#you will get a warning saying "In rm(temp) : object 'temp' not found." This is not a problem.
anovadf <- anovadf[2:nrow(anovadf),] #remove that first row of NAs we put in earlier


#fix RT data
datafrm$RT <- as.numeric(datafrm$RT) #variable must be numeric
datafrm$RT[which(datafrm$RT<150)] <- NA #responses faster than 150 ms are probably random key-presses: not possible to process and respond this fast
datafrm$RT2 <- ifelse(datafrm$error_dummy==0, datafrm$RT, NA)  #only correct trials
datafrm$RT2 <- log(datafrm$RT2)  #log-transform RT values. Ignore the warning message about "NaNs produced": this is because we had some NA values

#have to adjust a few more things to make R do the ANOVA correctly:
datafrm$group <- as.factor(datafrm$group)
datafrm$object_effects <- as.factor(datafrm$object_effects)

#get rid of dataless rows
datafrm <- datafrm[which(!is.na(datafrm$error_dummy)),]

#create dataframe where the coefficients and p-values we generate will go####

#one row per analysis per dataset
output <- data.frame(matrix(nrow = ((3 + 2*length(randomEffects))*length(datasets)), ncol = 5))  
names(output)<-c("Study", "Analysis", "text", "coefficient", "pval")

#column naming/numbering each of the analyses
analysesRan <- 1:length(randomEffects)
analyses <- c("rtFixed", paste0("rtRandom", analysesRan), 'errorFixed', paste0("errorRandom",analysesRan), "errorANOVA")
output$Analysis<-analyses

#column indicating the random effects structure of each analysis
output$text <- c("<fixed effects only>", randomEffects, "<fixed effects only>", randomEffects, "<n/a>")


#Multiverse####
startTime <- proc.time()
#prep
datasetvec <-  1:length(datasets) #we'll cycle through each dataset
lrf<-length(randomEffects)  #and each of the random effects specifications.
for(dset in(datasetvec)){  #Name each study in the output file
  output$Study[which(is.na(output$Study))[1]:(dset*(lrf*2+3))] <- datasets[dset]
}
#multiverse
tempmat <- sapply(c(datasetvec), FUN = function(dset){
  #tell the user which dataset we're on, so they can track our progress:
  print(paste(dset, datasets[dset], sep=": "))
  #create empty vectors which we'll populate with coefficients and p-values for RT and error analyses:
  outputrt <- vector()  
  outputer <- vector()
  #select the dataset to analyze:
  df<-datafrm[which(datafrm$study==datasets[dset]),]
  
  print("rt analyses") #tell the user we're doing RT analyses now
  
  #Run the reaction time analysis with only fixed effects and get coefficient and p-value
  print("fixed effects") #tell the user which analysis we're on
  #create function to run analysis for this random effects structure:
  fmlaRT <- as.formula(anteRT)
  #Return coefficient and p-value from object:race interaction:
  RTfixedCOEFandPVAL <- c(summary(lm(fmlaRT, data=df))$coefficients[4,1],
                          summary(lm(fmlaRT, data=df))$coefficients[4,4])
  
  #Run the reaction time random effects analyses and get a list of the cofficients and p-values:
  outputrts <- sapply(1:lrf, FUN=function(f){  #cycling through each random effects structure
    print(randomEffects[f]) #tell the user which analysis we're on
    #create formulae using the strings of fixed effects and random effects we created:
    fmlaRT <- as.formula(paste0(anteRT, randomEffects[f]))
    fmlaRT2 <- as.formula(paste0(anteRT2, randomEffects[f]))   
    #create function to run analysis for this random effects structure:
    RTval <- function(fmla1, fmla2){  
      return(tryCatch( #we'll put this in a tryCatch in case we get a warning, which probably means the analysis didn't converge. 
        
        #If we don't get a warning, return object:race interaction coefficient, as well as p-value from comparison of 
        #models with and without object x race interaction:
        expr = c(summary(lmer(fmla1, data=df))$coefficients[4,1], anova((lmer(fmla1, data=df)), (lmer(fmla2, data=df)))[2,8]),
        
        #But if we do get a warning, we can try rerunning the model estimation procedure with more iterations.
        #(could also try using a different optimizer, but here I want to keep procedures consistent across datasets)
        warning=function(w){
          #If adding more iterations doesn't work, we'll need to return a vector of NAs, so to do this we need to have a tryCatch within a tryCatch:
          return(tryCatch(
            
            expr = cbind(summary(
              lmer(fmla1, data=df, control=lmerControl(optCtrl=list(maxfun=3e4))))$coefficients[4,1],
              anova((lmer(fmla1, data=df, control=lmerControl(optCtrl=list(maxfun=3e4)))),
                    (lmer(fmla2, data=df, control=lmerControl(optCtrl=list(maxfun=3e4)))))[2,8]),
            
            warning=function(wn){
              #If we're still getting a warning, print customized warning message:
              print("Regression (RT data) did not converge") 
              #and return NAs for that analysis:
              return(c(NA,NA))  
            }
          ))
        }
      ))
    }
    #run the function we just specified to get coefficient and p-value (or NA, 
    #if there was a warning message):
    vals <- RTval(fmlaRT,fmlaRT2) 
    outputrt <- c(outputrt, vals) #add that p-value to our vector of RT p-values
    return(outputrt)
  }
  )
  
  print("error analyses") #tell the user we're doing error analyses now 
  
  #Run the error analysis with only fixed effects and get coefficient and p-value
  print("fixed effects") #tell the user which analysis we're on
  #create formula to run analysis:
  fmlaER <- as.formula(anteER)
  #Get coefficient and p-value from object:race interaction:
  ERfixedCOEFandPVAL <- c(summary(glm(fmlaER, data = df, family=binomial(link=logit)))$coefficients[4,1],
                          summary(glm(fmlaER, data = df, family=binomial(link=logit)))$coefficients[4,4])
  
  #Run error random effects analyses and make a list of the coefficients and p-values:
  outputers <- sapply(1:lrf, FUN=function(f){ #cycling through each random effects structure
    print(randomEffects[f]) #tell the user which analysis we're on
    #create formula using the strings of fixed effects and random effects we created:
    fmlaER <- as.formula(paste0(anteER, randomEffects[f]))
    #create function to run analysis for this random effects structure:
    reg3 <- function(fmla){
      return(tryCatch(  #we'll put this in a tryCatch in case we get a warning, which probably means the analysis didn't converge. 
        #If no warning, return object:race interaction coefficient and p-value from logistic regression:
        expr = cbind(summary(glmer(fmla, data = df, 
                                   family=binomial(link=logit)))$coefficients[4,1],
                     summary(glmer(fmla, data = df, 
                                   family=binomial(link=logit)))$coefficients[4,4]),
        
        #But if we get a warning, we can try rerunning the model estimation procedure with more iterations.
        #If that doesn't work, though, we'll need to return a vector of NAs, so to do this we need to have a tryCatch within a tryCatch:
        warning=function(w){  
          return(tryCatch(
            expr = cbind(summary(glmer(fmla, data = df, 
                                       family=binomial(link=logit),
                                       lmerControl = glmerControl(optCtrl = list(maxfun=3e4))
            ))$coefficients[4,1],
            summary(glmer(fmla, data = df, 
                          family=binomial(link=logit),
                          lmerControl = glmerControl(optCtrl = list(maxfun=3e4))))$coefficients[4,4]),
            warning=function(wn){
              #if we're still getting a warning, print customized warning message:
              print("Logistic regression (error data) did not converge")
              #and return NAs for that analysis:
              return(c(NA, NA))
            }
          ))
        }
      ))
    }
    #run the function we just specified to get coefficient and p-value (or NAs, 
    #if there was a warning message):
    vals <- reg3(fmlaER)
    #add those values to our vector of error coefficients and p-values
    outputer <- c(outputer, vals) 
    return(outputer)
  }
  )
  #error rate ANOVA:
  anovamod <- lm(error_rate~group*object_effects, data = anovadf[which(anovadf$study==df$study[1]),])
  anovaoutput <- Anova(anovamod, type=3)  #use Type 3 sums of squares
  
  #make vector of coefficients and p-values from all the analyses for this dataset:
  outputvec <- unlist(c(RTfixedCOEFandPVAL, outputrts, ERfixedCOEFandPVAL, outputers, summary(anovamod)$coefficients[4,1], anovaoutput[4,4])) #make vector of coefficients and p-values from all the analyses for this dataset
  #put vector in appropriate format to be added to output file:
  coeffs_pvals <- t(as.matrix(outputvec, nrow=2)) 
  return(coeffs_pvals)
}
)
#the above sapply() function returns a matrix: 
# we want to convert that to a vector and assign it to 
# the p-value column of the output dataframe:
output[,c("coefficient")] <- t(matrix(unlist(tempmat), nrow=2))[,1]
output[,c("pval")] <- t(matrix(unlist(tempmat), nrow=2))[,2]
#create column indicating if object:race interaction was significant in each analysis:
output$sig <- ifelse(output$pval <= 0.05, T, F)
endTime <- proc.time()
print((endTime - startTime)/60)  #third value tells us how long this took us to run

#took about 23 minutes to run 1 dataset
#took about 5 hours to run 19 datasets 