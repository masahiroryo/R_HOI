##################################################
#
# Interaction analyses based on Kerry & Okada 2012
# This script needs to be run after running the RF
# within the same environment
#
#                    2017.06.12    Masahiro Ryo
#                         masahiroryo@gmail.com
#
##################################################

###############################################################
###     PREPARATION PART
###############################################################
### 1. reading packages (if not installed, automatically install)
package.list <- c("party", "foreach","doSNOW","mlr")
tmp.install <- which(lapply(package.list, require, character.only = TRUE)==FALSE)
if(length(tmp.install)>0) install.packages(package.list[tmp.install])
lapply(package.list, require, character.only = TRUE)

# reading the function to calculate the relative importance of interactions
source("Appendix_S2d_interaction_function.R")  # the function estimating variable interaction importance

### 2. reading selected variables and best RF model that are gained
###    after running Appendix_S2_a_RF.R
# csv storing a response on 1st col. and predictors on the rest cols.
data <- data.frame(read.csv("data.csv", na.strings = "NA", header=T, skip=1))
# RF: recommended to reduce the number of predictors by variable selection
Best.RF = cforest(alpha_fam~.,data=data, controls = cforest_unbiased(ntree = 1000))
var.list = colnames(Best.RF@data@get("input"))
cmb.3 = combn(var.list,3, simplify = F) # all combinations


###############################################################
###     RUNNING PART
###############################################################
cl<-makeCluster(4) # according to the number of CPU cores
registerDoSNOW(cl) # parallel computing on Windows OS

t<-proc.time()  # stopwatch: start

E.3 = unlist(foreach(i = 1:length(cmb.3), .packages='party') %dopar% {
  intimp(Best.RF, cmb.3[[i]])}) # calculation (main part)

proc.time()-t   # stopwatch: end
stopCluster(cl) # parallel computing: end


###############################################################
###     VISUALIZATION PART (1)
###############################################################
### 1. extracting top 10
top = 10 # check only important ones
high.3 = which(rank(-E.3)<=top)
high.list.3 = as.list(rep(NA,top))

for(i in c(1:top)){
  high.list.3[[i]] <- as.list(rep(NA,2))
  high.list.3[[i]][[1]] = cmb.3[[which(rank(-E.3)==i)]]
  high.list.3[[i]][[2]] = E.3[which(rank(-E.3)==i)]
}

high.list.3 # most important interactions

### 2. frequency distribution of the improtance scores
par(mfrow=c(2,1))
par(mar=c(4,4,1,1))

hist(E.3, col = "#ff00ff40", border = "#ff00ff", breaks = 100, add = TRUE)
hist(E.3, col = "#ff00ff40", border = "#ff00ff", breaks = 100, add = TRUE)


###############################################################
###     VISUALIZATION PART (2) Partial dependence plots
###############################################################
### 1. fitting RF again
task = makeRegrTask(data = data, target = "alpha_fam")
learner = makeLearner(cl="regr.cforest", predict.type = "response", par.vals = list(ntree=1000))
fit = train(learner, task)

### 2. calculating PD scores: e.g. elevation-forest5km-biogeoclimate
ngrid = 30
pd.ind.1 = generatePartialDependenceData(fit, task, c("masl","Woods_prop_5km","biogeo_class"), interaction=T, individual=F, gridsize=ngrid, resample="bootstrap")
result.1 = data.frame(pd.ind.1$data)
ncat.1 = length(summary(result.1[,4])) #number of categorical 
length.1 = nrow(result.1)
result.1 = result.1[-c((1+length.1*(ncat.1/(ncat.1+1))):length.1),]

### 3. preparing axes for visualization
X.1 = result.1[1:ngrid,2]
Y.1 = result.1[seq(1, nrow(result.1)/(ncat.1), by = ngrid),3]
C.1 = names(summary(result.1[,4]))
dupX.1 = which(duplicated(X.1))
dupY.1 = which(duplicated(Y.1))
if(is.factor(Y.1)) Y.1 = as.integer(Y.1)

for(i in c(1:ncat.1)){
  tmp = subset(result.1[,1],result.1[,4]==C.1[i])
  eval(parse(text= paste0(
    "Z.1_",i," = matrix(tmp, nrow=length(X.1),ncol=length(Y.1))"))) 
}
if(length(dupX.1)!=0)X.1 = X.1[-dupX.1]
if(length(dupY.1)!=0)Y.1 = Y.1[-dupY.1]

for(i in c(1:ncat.1)){
  
  eval(parse(text= paste0("
                          if(length(dupX.1)!=0&&length(dupY.1)!=0){
                          Z.1_",i," = Z.1_",i,"[-dupX.1,-dupY.1]
                          }else if(length(dupX.1)!=0){
                          Z.1_",i," = Z.1_",i,"[-dupX.1,]
                          }else if(length(dupY.1)!=0){Z.1_",i," = Z.1_",i,"[,-dupY.1]
                          }", sep="")))
}


### 4. visualization
par(mfrow = c(2,3))
par(mar = c(2, 2, 2, 2))
for(i in c(1,3,4,2,5,6)){
  eval(parse(text= paste0(
    "persp(X.1,Y.1,Z.1_",i,",theta = 45, phi = 20, ticktype='detailed', 
    xlab='masl',ylab='Woods%(5km)',zlab='partial dependence',d=1,
    zlim=range(result.1[,1]), shade=0.1,r=5)"))) 
  title(C.1[i], line = 0)
}
