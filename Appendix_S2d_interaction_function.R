##################################################
#
# Function: calculating the importance score of interaction
# Interaction analyses based on Kerry & Okada 2012
# extended to 3way interactions
# (as well as 2way interactions)
#
#                    2017.06.12    Masahiro Ryo
#                         masahiroryo@gmail.com
#
##################################################

intimp <- function(object, shuffle) {
  response <- object@responses
  input <- object@data@get("input")
  inp <- initVariableFrame(input, trafo = NULL)
  y <- object@responses@variables[[1]]
  CLASS <- all(response@is_nominal)
  ORDERED <- all(response@is_ordinal)
  if (CLASS) {
    error <- function(x, oob) mean((levels(y)[sapply(x, which.max)] != y)[oob])
  } else {
    if (ORDERED) {
      error <- function(x, oob) mean((sapply(x, which.max) != y)[oob])
    } else {
      error <- function(x, oob) mean((unlist(x) - y)[oob]^2)
    }
  }
  
  perror <- matrix(0, nrow = length(object@ensemble), ncol = 2^length(shuffle))
  
  cmb = list()
  for(i in 1:length(shuffle)) cmb = c(cmb,combn(shuffle, i, simplify = F))
  
  for (b in 1:length(object@ensemble)) {
    tree <- object@ensemble[[b]]
    oob <- object@weights[[b]] == 0
    p <- .Call("R_predict", tree, inp, 0, -1L, PACKAGE = "party")
    eoob <- error(p, oob)

    # set.seed(123)
     rnd = sample(which(oob),length(which(oob)),replace=F)

    for(i in 1:length(cmb)){
      permuted = input
      permuted[oob, cmb[[i]]] = permuted[rnd, cmb[[i]]]
      inp.perm = initVariableFrame(permuted, trafo = NULL)
      p.perm = .Call("R_predict", tree, inp.perm, 0, -1L,  PACKAGE = "party")
      perror[b, i] <- error(p.perm, oob) - eoob
      
    }
  }
  if(length(shuffle)==2) {
    perror[, 4] = perror[, 1] + perror[, 2] - perror[, 3]
  }    
  if(length(shuffle)==3){
    perror[, 8] = perror[, 1] + perror[, 2] + perror[, 3] - 
                  perror[, 4] - perror[, 5] - perror[, 6] + perror[, 7]
  }
  return(MeanDecreaseAccuracy = mean(perror[,ncol(perror)]))
}
