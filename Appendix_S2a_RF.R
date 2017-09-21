##################################################
#
# Random forest analyses with variable selection
# based on Hapfelmeier & Ulm 2013
#
#                    2017.05.12    Masahiro Ryo
#                         masahiroryo@gmail.com
#
##################################################

###############################################################
###     PREPARATION PART
###############################################################
### 1. reading packages (if not installed, automatically install)
package.list <- c("party", "caret","mlr")
tmp.install <- which(lapply(package.list, require, character.only = TRUE)==FALSE)
if(length(tmp.install)>0) install.packages(package.list[tmp.install])
lapply(package.list, require, character.only = TRUE)

# reading the code of Hapfelmeier's Random Forests
source("Appendix_S2b_RF_function.R")
# formula: object of class "formula".
# data: data frame containing the variables.
# nperm: Number of permutation steps used for the permutation test.
# ntree: Number of trees in the Random Forest.
# alpha: Significance level of the permutation test.

### 2. reading the dataset, extracting attribute information
data <- data.frame(read.csv("BDM_Data.csv", na.strings = "NA", header=T, skip=5))


#######################################################################
###     1) preparing an outcome matrix that contains all results and 
###        to be exported as a csv file in the end, and 
###     2) separating predictor variables from the dataset
#######################################################################

### 1. extracting information about predictor and response variables
eliminate <- -c(1:7) # select columns which are not used as predictor variables
nvars     <- ncol(data[eliminate])
ntargets  <- 2 # number of response variables (alpha_EPT, alpha_IBCH...)
varnames  <- colnames(data[eliminate]) # name list for predictor variables
resnames  <- colnames(data)[c(2,4)] # name list for response variables

results   <- data.frame(matrix(NA, ncol=3+nvars, nrow=ntargets))
colnames(results) <- append(c("Response","R2_fit","R2_OOB"),varnames)

     ### [tips] ############## 
     #  1st col. [Response] name of Response variable
     #  2nd col. [RMSE_fit] Root mean squared error for "regression" and accuracy for "classification"
     #  3rd col. [R2_fit]   variance explained with a model (%)
     #           [ _OOB] OOB stands for Out-of-bag (Breiman 1996)
     #                   NOT variance explained (fitting) BUT predictive power (validation)
     #                   This indicates how accurately the built model can predict for a new dataset 


#################################################################################
###     Building a Random Forest model (Breinman 2001) for each response variable
###     For each response variable,
###     1) estimating relative importance score of each predictor variable
###     2) and evaluating the model performance (both fitting and validation)
#################################################################################

### 1. the parameters in RF models 
ntree <- 1000 # the number of decision tree models in a RF algorithm
mtry  <- ceiling(sqrt(nvars)) # the number of predictor variables compared 
                             # at each node in each tree model
alpha <- 0.01    # statistical significance level for Hapfelmeier's RF
              # note that intentionally this is set to 1 to extract all information.
nperm <- 5000  # the number of permutative iterations for estimating p-value (default=400, see Hapfelmeier and Ulm 2013 )

### 2. for loop function for response variables
for (Y in resnames){

  ## 2.1. NA removal from a target response variable (note that not from explanatory variables)
  row.na <- which(is.na(data[,Y]))
  if(length(row.na)>0){data.tmp <- data[-row.na,]}else{data.tmp <- data}
        # note that "data.tmp" is used for modelling (not "data")
  
  ## 2.2. building Hapfelmeier's Random forest model: parameters are ntree (number of trees) and mtry
  formula.all <- as.formula(paste(Y, paste(varnames, collapse=" + "), sep=" ~ "))
  Hapfelmeier.RF <- selection.parallel(formula=formula.all,data=data.tmp,nperm=nperm,ntree=ntree,alpha=1) # to show p.val for all
  var.selected <- which(Hapfelmeier.RF$p.values*nvars<alpha)
  formula.best <- as.formula(paste(Y, paste(varnames[var.selected], collapse=" + "), sep=" ~ "))
  mtry2  <- ceiling(sqrt(length(var.selected))) 
  Best.RF <- cforest(formula=formula.best,data=data.tmp, controls = cforest_unbiased(ntree = ntree, mtry=mtry2))

  if(Y=="alpha_fam"){
    Hapfelmeier.RF.alpha_fam = Hapfelmeier.RF
    Best.RF.alpha_fam = Best.RF
  }else{
    Hapfelmeier.RF.alpha_EPT = Hapfelmeier.RF
    Best.RF.alpha_EPT = Best.RF
  }

  ## 2.3. evaluation
  predictor.varimp <- varimp(Hapfelmeier.RF$forest) # estimating relative importance score of predictor variables 
  eval(parse(text= paste0("obs <- data.tmp$", Y)))     # the observed values of response variable
  pred <- predict(Best.RF)          # the predicted values of response variable
  accuracy.fitting <- postResample(pred, obs)          # evaluating fitting performance
  accuracy.validation <- caret:::cforestStats(Best.RF)  # evaluating predictive performance

  # 2.4. storing the performance and relative importance scores in the outcome matrix
  i_row <-which(Y==resnames)*2-1
  # for relative imporatnce
  results[i_row, 1] <- Y
  results[i_row, c(2)] <- accuracy.fitting[2]
  results[i_row, c(3)] <- accuracy.validation[2]
  results[i_row, -c(1:3)] <- predictor.varimp
  # for p-value
  results[i_row+1,-c(1:3) ] <- Hapfelmeier.RF$p.values*nvars
  write.csv(results,"result.csv", quote=FALSE, row.names=FALSE)  

  
  # 2.5. partial dependence plots for all the variables in Best.RF
  task = makeRegrTask(data = data.tmp[,c(Y, names(var.selected))], target = Y)
  learner = makeLearner(cl="regr.cforest", predict.type = "response", par.vals = list(ntree=1000, mtry=mtry2))
  fit = train(learner, task)
  
  par(mfrow=c(4,5))
  par(mar=c(4,3,1,1))
  n.grids = 30
  for(i in names(var.selected)){
    pdp = generatePartialDependenceData(fit, task,  i, individual=F, resample="bootstrap", gridsize = n.grids)
    if(is.factor(data.tmp[,i])){
      barplot(pdp$data[,1], names.arg= as.character(pdp$data[,2]), xlab=i, ylab="")
    }else{
      plot(pdp$data[,1]~pdp$data[,2], type="l", xlab=i, ylab="")
    }
  }
  # just focus on some of the most important predictors
} 
