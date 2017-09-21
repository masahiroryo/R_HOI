##################################################
#
# The R script by Hapfelmeier & Ulm 2013
# modified by Masahiro Ryo for Windows OS
#
#                    2017.05.12    Masahiro Ryo
#                         masahiroryo@gmail.com
#
##################################################

myvarimp <- function(object, mincriterion = 0, conditional = FALSE, threshold = 0.2,
                     pre1.0_0 = conditional, varID) {
  response <- object@responses
    input <- object@data@get("input")
    xnames <- colnames(input)
    inp <- party::initVariableFrame(input, trafo = NULL)
    y <- object@responses@variables[[1]]
    error <- function(x, oob) mean((unlist(x) - y)[oob]^2)
    perror <- rep(0, length(object@ensemble))
    for (b in 1:length(object@ensemble)) {
        tree <- object@ensemble[[b]]
        w <- object@weights[[b]]
        w[w == 0] <- 1
        oob <- object@weights[[b]] == 0
        p <- .Call("R_predict", tree, inp, mincriterion, -1L,
            PACKAGE = "party")
        eoob <- error(p, oob)
        j <- varID
        p <- .Call("R_predict", tree, inp, mincriterion,
                    as.integer(j), PACKAGE = "party")
        perror[b] <- (error(p,oob) - eoob)
    }
    return(MeanDecreaseAccuracy = mean(perror))
}
environment(myvarimp) <- environment(varimp)



# Function for parallel variable selection.
selection.parallel <- function(formula, data, nperm = 400, ntree = 100, alpha = 0.05)
  # formula: object of class "formula".
  # data: data frame containing the variables.
  # nperm: Number of permutation steps used for the permutation test.
  # ntree: Number of trees in the Random Forest.
  # alpha: Significance level of the permutation test.
  {

################################################
# Note that parallel does not work on Windows OS.
  library("foreach")  # to parallelize the main processes
  library("doSNOW")   # to activate foreach in windows
  
  x.names <- all.vars(formula)[-1]
  y.names <- all.vars(formula)[1]
  terms. <- terms.formula(formula)
  x.formula <- attr(terms., "term.labels")
  y.formula <- as.character(attr(terms., "variables"))[2]
  mtry <- ceiling(sqrt(length(x.formula)))
  dat <- subset(data, select = c(y.names, x.names))
  forest <- party::cforest(formula, data = dat, controls = cforest_unbiased(mtry = mtry, ntree = ntree))
  obs.varimp <- varimp(forest)
  perm.mat <- matrix(NA, ncol = length(x.names), nrow = nperm, dimnames = list(1:nperm, x.names))

  cl<-makeCluster(4) #change to your number of CPU cores
  registerDoSNOW(cl)
  
  
  for (j in x.names) {
    cat("\r", "Processing variable ", which(j == x.names), " of ", length(x.names)); flush.console()
    perm.dat <- dat
    perm.mat[, j] <- unlist(foreach(i = 1:nperm, .packages='party',.export="myvarimp") %dopar% {
      perm.dat[, j] <- sample(perm.dat[, j]);
      myvarimp(party::cforest(formula, data = perm.dat, controls = cforest_unbiased(mtry = mtry, ntree = ntree)), varID = which(x.names == j))
    })
   
  }
  stopCluster(cl)
  
  p.vals <- sapply(x.names, function(x) sum(perm.mat[, x] >= obs.varimp[which(x == x.names)]) / nperm)
  p.vals.bonf <- p.vals * length(p.vals)

  if (any(p.vals < alpha)) {
   selection <- names(p.vals)[which(p.vals < alpha)]
   mtry <- ceiling(sqrt(length(selection)))
   forest <- cforest(as.formula(paste(y.formula, ".", sep = " ~ ")), data = subset(dat, select = c(y.names, selection)),
                     controls = cforest_unbiased(mtry = mtry, ntree = ntree))
   p <- p.vals[which(p.vals < alpha)]
   }
  if (any(p.vals.bonf < alpha)) {             
   selection.bonf <- names(p.vals.bonf)[which(p.vals.bonf < alpha)]                         
   mtry <- ceiling(sqrt(length(selection.bonf)))
   forest.bonf <- party::cforest(as.formula(paste(y.formula, ".", sep = " ~ ")), data = subset(dat, select = c(y.names, selection.bonf)),
                          controls = cforest_unbiased(mtry = mtry, ntree = ntree))
   p.bonf <- p.vals.bonf[which(p.vals.bonf < alpha)]
   }
  if (!any(p.vals < alpha)) {
   selection <- c(); forest <- c(); p <- c()
   }
  if (!any(p.vals.bonf < alpha)) {
   selection.bonf <- c(); forest.bonf <- c(); p.bonf <- c()
   }
  Y <- as.numeric(as.character(dat[,y.names]))
  oob.error <- ifelse(length(selection) != 0, mean((Y - as.numeric(as.character(predict(forest, OOB = T))))^2), mean((Y - ifelse(all(Y %in% 0:1), round(mean(Y)), mean(Y)))^2))
  oob.error.bonf <- ifelse(length(selection.bonf) != 0, mean((Y - as.numeric(as.character(predict(forest.bonf, OOB = T))))^2), mean((Y - ifelse(all(Y %in% 0:1), round(mean(Y)), mean(Y)))^2))
  cat("\n", "\r"); flush.console()
  return(list("selection" = selection, "forest" = forest, "oob.error" = oob.error, "p.values" = p,
              "selection.bonf" = selection.bonf, "forest.bonf" = forest.bonf, "oob.error.bonf" = oob.error.bonf, "p.values.bonf" = p.bonf))
}

