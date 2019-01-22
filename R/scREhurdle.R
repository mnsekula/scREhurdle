#' Hurdle model with random effects for discrete scRNA-seq data
#'
#' Wrapper function used to implement the methodology from the
#' "Detection of differentially expressed genes in discrete-single
#' cell RNA sequencing data using a hurdle model with correlated
#' random effects" manuscript.
#' 
#'
#' @param Y data.frame or matrix where the rows correspond to genes and
#'        columns correspond to cells. Y must contain integers.
#' @param treatGroup vector or factor of length ncol(Y) indicating the 
#'        treatment assigments (control or treatment) of the cells.
#' @param useCDR whether CDR, the proportion of genes expressed per cell, is 
#'        calculated and used in model matrix (X).
#' @param typeRE set type of random effect to use.  Choices are "none", "ind",
#'        or "corr".
#' @param subpop vector or factor of length ncol(Y) indicating the 
#'        cluster/subpopulation assignment of cells. Used ONLY in the 
#'        correlated RE model. Cells in different treatment groups should 
#'        be clustered into separate subpopulations.
#' @param addCovariates data.frame containing values of additional covariates to
#'        add to model matrix (X). It is recommended that continuous values are scaled
#'        to have mean 0 and standard deviation 0.5.
#' @param coefSamps when set to "treatment", only regression coefficients related to 
#'        determining DE genes (i.e., the treatment indicator coefficients) are stored in the 
#'        results, thereby minimizing object size. The "all" option will store all regression 
#'        coefficients. Note: 2 x nrow(Y) x (ncol(X)+1) parameters are stored with "all". 
#' @param parSamps character vector indicating additional model parameters to 
#'        store. Choices include: "omega", "phi", "lambda1", "lambda2", "sigma2".
#'        Additional choices for correlated RE model include: "gamma_t", 
#'        "omega_star".
#' @param adjustMethod p.adjust.method passed to p.adjust() function.
#' @param adapt_engaged passed to vb() function in rstan.
#' @param tol passed to vb() function in rstan.
#' @param eta passed to vb() function in rstan. Eta values as high as 0.4 may work for some datasets. 
#'        Increasing eta should decrease computational time but could cause errors from Stan. Reduce 
#'        eta if model mispecification errors occur.
#' @param output_samples passed to vb() function in rstan.
#' @param stan_seed seed for random number generation passed to vb() function in rstan.
#' @param ... additional parameters to be passed to vb() function in rstan.
#'
#' @return scREhurdle.fit object containing:
#' 
#' \itemize{ 
#'   \item{stanFit}{: an object of stanfit-class.}
#'   \item{inputData}{: list of data input into the Stan model.}
#'   \item{deTab}{: data.frame consisting of beta estimates, statistics, and p-values
#'         associated with identifying DE genes.}
#'   \item{geneTab}{: data.frame providing estimates for all stored model parameters 
#'         corresponding to genes (i.e. beta_L, beta_C, zeta_L, zeta_C, phi).}
#'   \item{cellTab}{: data.frame providing cellular details and stored model parameter 
#'         estimates corresponding to cells (i.e. treatGroup, CDR, subpop, omega).}
#'   \item{hyperEst}{: vector with stored hyperparameter estimates of "lambda1", "lambda2", 
#'         and/or "sigma2" if included in parSamps.}
#' }
#' 
#' @examples
#' \dontrun{
#' data(toyDat)
#' gene.names <- rownames(toyDat$scData)
#' treatment <- substring(colnames(toyDat$scData),1,1)
#' 
#' ## DE analysis with independent random effects
#' exIRE <- scREhurdle(Y=toyDat$scData,treatGroup=treatment,useCDR=TRUE,typeRE = "ind",
#' coefSamps="treatment",eta=0.4,stan_seed=523)
#' 
#' ## DE genes
#' gene.names[which(exIRE$deTab$chisq.padj <= 0.05)]
#' 
#' ## DE analysis with correlated random effects
#' exCRE <- scREhurdle(Y=toyDat$scData,treatGroup=treatment,useCDR=TRUE,typeRE = "corr",
#' subpop=toyDat$subpop, coefSamps="treatment",eta=0.4,stan_seed=523)
#' 
#' ## DE genes
#' gene.names[exCRE$deTab$chisq.padj <= 0.05]
#' }
#'
#' @export
scREhurdle <- function(Y, treatGroup, useCDR = TRUE, typeRE = "ind", subpop = NULL, addCovariates = NULL,
                     coefSamps = c("treatment", "all"), parSamps = NULL, adjustMethod = "BH",
                     adapt_engaged = FALSE, tol = 1e-4, eta = 0.2, output_samples = 1000, stan_seed = 123, ...){
  
  typeRE <- match.arg(typeRE, c("none","ind","corr"))
  if(!is.null(parSamps)){
    parSamps <- match.arg(parSamps, c("omega", "phi", "lambda1", "lambda2", "omega_star", "gamma_t", "sigma2"), 
                          several.ok = TRUE)
    if(typeRE == "ind" & "omega_star" %in% parSamps | typeRE == "ind" & "gamma_t" %in% parSamps){
      stop("omega_star and/or gamma_t cannot be sampled from the indepedent RE model: remove these terms from parSamps")
    }
    if(typeRE == "none" & any(c("omega","omega_star","gamma_t","sigma2") %in% parSamps)){
      stop("only phi, lambda1, lambda2 can be sampled from the no RE model: remove other terms from parSamps")
    }
  }
  coefSamps <- match.arg(coefSamps, c("treatment", "all"))
  
  # Check for count data
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol     
  isCountData <- all(is.wholenumber(Y))
  if(identical(isCountData, FALSE)){
    stop("Y must consist of only integers")
  }
  
  # Check for number of groups
  if(length(unique(treatGroup)) != 2){
    stop("treatGroup must consist of only two groups")
  }
  
  G <- nrow(Y)
  N <- ncol(Y)
  treatGroup <- as.factor(treatGroup)
  X <- X.mod <- model.matrix(~ treatGroup)
  X[,2] <- scale(X[,2], scale = FALSE)
  cell.mat <- data.frame(treatGroup)
  if(!is.null(colnames(Y))){
    rownames(cell.mat) <- colnames(Y)
  }else{
    rownames(cell.mat) <- paste0("Cell",1:N)
  }
  if(isTRUE(useCDR)){
    genesPerCell <- apply(Y, 2, function(x) sum(x > 0))
    CDR <- genesPerCell/nrow(Y)
    cell.mat$CDR <- CDR
    CDR.scale <- scale(CDR)/2
    CDR.mod <- model.matrix(~ CDR.scale - 1)
    X <- cbind(X, CDR.mod)
  }
  
  if(!is.null(addCovariates)){
    for(a in 1:ncol(addCovariates)){
      X.temp <- model.matrix(~ addCovariates[,a] - 1)
      if(is.factor(addCovariates[,a])){
        colnames(X.temp) <- paste0(names(addCovariates)[a], levels(addCovariates[,a]))
      }else{
        colnames(X.temp) <- names(addCovariates)[a]
      }
      X <- cbind(X, X.temp)
      cell.mat <- cbind(cell.mat,addCovariates[,a])
      names(cell.mat)[ncol(cell.mat)] <- names(addCovariates)[a]
    }
  }
  
  pars <- switch(coefSamps,
                 "treatment" = c("beta_L1", "beta_C1"),
                 "all" = c("beta_L", "beta_C"))
  addArgs <- list(...)
  prior.scale <- c(10, rep(2.5, ifelse(typeRE == "none",ncol(X)-1,ncol(X))))
  gene_data <- list(
    G = G,           
    N = N,           
    M = ifelse(typeRE == "none",ncol(X),ncol(X)+1),   
    y = Y,
    X = X,
    prior_s = prior.scale)
  
  if(typeRE == "corr"){
    ## Correlated random effects
    pars <- c(pars, parSamps)
    subpop <- as.factor(subpop)
    cell.mat$subpop <- subpop
    ## Assign gammas to clusters (control clusters are first)
    c.clust <- sort(unique(subpop[which(X.mod[,2]==0)]))
    t.clust <- sort(unique(subpop[which(X.mod[,2]==1)]))
    K0 <- length(c.clust)
    K1 <- length(t.clust)
    gamAssign <- numeric(N)
    for(k in 1:K0){
      c.ind <- c.clust[k]
      gamAssign[which(subpop == c.ind)] <- k
    }
    for(k in 1:K1){
      t.ind <- t.clust[k]
      gamAssign[which(subpop == t.ind)] <- K0+k
    }
    J <- model.matrix(~as.factor(gamAssign) -1)
    
    subpop_data <- list(
      K = ncol(J),
      K0 = K0,
      J = J)
    gene_data <- c(gene_data, subpop_data)
    
    cat("Running CRE Stan model... \n")
    mod <- stanmodels$CRE
    fit <- do.call(rstan::vb,c(list(
      mod,
      data = gene_data,    # named list of data
      adapt_engaged = adapt_engaged,
      eta = eta,
      tol_rel_obj = tol,
      seed = stan_seed,    # seed for Markov Chain
      pars = pars,
      include = TRUE,
      output_samples = output_samples
    ),addArgs))
    
  }else{
    pars <- c(pars, parSamps)
    if(typeRE == "ind"){
      ## Independent random effects
      cat("Running IRE Stan model... \n")
      mod <- stanmodels$IRE
    }else{
      ## No random effects
      cat("Running NRE Stan model... \n")
      mod <- stanmodels$NRE
    }
    fit <- do.call(rstan::vb,c(list(
      mod,
      data = gene_data,
      adapt_engaged = adapt_engaged,
      eta = eta,
      tol_rel_obj = tol,
      seed = stan_seed,
      pars = pars,
      include = TRUE,
      output_samples = output_samples
    ),addArgs))
  }
  
  ########## Sampling & Data Analysis ##########
  fit.dat <- rstan::As.mcmc.list(fit)[[1]]
  geneNames <- rownames(Y)
  
  if(coefSamps == "treatment"){
    betaL1.ind <- grep(pars[1], colnames(fit.dat))
    betaC1.ind <- grep(pars[2], colnames(fit.dat))
    betaL1.samp <- fit.dat[,betaL1.ind]
    betaC1.samp <- fit.dat[,betaC1.ind]
    gene.mat <- matrix(apply(fit.dat[,c(betaL1.ind,betaC1.ind)],2,mean),nrow=nrow(Y))
    rownames(gene.mat) <- geneNames
    colnames(gene.mat) <- c(paste0(colnames(X)[2],"_L"), paste0(colnames(X)[2],"_C"))
    fit.red <- fit.dat[,-c(betaL1.ind,betaC1.ind)]
  }else{
    betaL.ind <- grep(pars[1], colnames(fit.dat))
    betaC.ind <- grep(pars[2], colnames(fit.dat))
    betaL1.samp <- fit.dat[,betaL.ind[(nrow(Y)+1):(2*nrow(Y))]]
    betaC1.samp <- fit.dat[,betaC.ind[(nrow(Y)+1):(2*nrow(Y))]]
    gene.mat <- matrix(apply(fit.dat[,c(betaL.ind,betaC.ind)],2,mean),nrow=nrow(Y))
    rownames(gene.mat) <- geneNames
    if(typeRE == "none"){
      coefNames <- c(gsub("[()]", "", colnames(X)))
    }else{
      coefNames <- c(gsub("[()]", "", colnames(X)),"zeta")
    }
    colnames(gene.mat) <- c(paste0(coefNames,"_L"), paste0(coefNames,"_C"))
    fit.red <- fit.dat[,-c(betaL.ind,betaC.ind)]
  }
  gene.mat <- data.frame(gene.mat)
  treat.mat <- matrix(NA, nrow = nrow(Y), ncol = 9)
  rownames(treat.mat) <- geneNames
  colnames(treat.mat) <- c(paste0(colnames(X)[2],"_L"),"L.Z","L.pval", paste0(colnames(X)[2],"_C"), "C.Z", "C.pval", "chisq", "chisq.pval", "chisq.padj")
  for(i in 1:nrow(Y)){
    treat.mat[i,] <- unlist(WaldBayes(betaL1.samp[,i], betaC1.samp[,i]))
  }
  treat.mat[,9] <- p.adjust(treat.mat[,8], method = adjustMethod)
  treat.mat <- data.frame(treat.mat)
  
  if(!is.null(parSamps)){
    parEst <- apply(fit.red,2,mean)[-ncol(fit.red)]
    if("gamma_t" %in% parSamps){
      cell.mat$gamma_t <- J%*%matrix(parEst[grep("gamma_t",names(parEst))])
      parEst <- parEst[-grep("gamma_t",names(parEst))]
    }
    if("omega_star" %in% parSamps){
      cell.mat$omega_star <- parEst[grep("omega_star",names(parEst))]
      parEst <- parEst[-grep("omega_star",names(parEst))]
    }
    if("omega" %in% parSamps){
      cell.mat$omega <- parEst[grep("omega",names(parEst))]
      parEst <- parEst[-grep("omega",names(parEst))]
    }
    if("phi" %in% parSamps){
      gene.mat$phi <- parEst[grep("phi",names(parEst))]
      parEst <- parEst[-grep("phi",names(parEst))]
    }
    if("sigma2" %in% parSamps){
      if(length(grep("sigma2",names(parEst))) == 1){
        names(parEst)[grep("sigma2",names(parEst))] <- "sigma2"
      }else{
        names(parEst)[grep("sigma2",names(parEst))] <- c("sigma2_*","sigma2_0","sigma2_1")
      }
    }
  }else{
    parEst <- NULL
  }
  res <- list(stanFit = fit, inputData = gene_data, deTab = treat.mat, geneTab = gene.mat, cellTab = cell.mat, hyperEst = parEst)
  class(res) <- "scREhurdle.fit"
  return(res)
}

#' Two-dimensional Wald statistic
#' 
#' Internal function used in scREhurdle used to calculate
#' the 2-D Wald statistic
#' 
#' @param x vector of beta_L1 samples
#' @param y vector of beta_C1 samples
#' @noRd
WaldBayes <- function(x,y){
  est.hat <- matrix(c(mean(x), mean(y)), nrow = 1, ncol = 2)
  n <- length(x) #n = length(x) = length(y)
  x.var.est <- sum((x-est.hat[1])^2)/n #var(x)
  y.var.est <- sum((y-est.hat[2])^2)/n #var(y)
  xy.cov.est <- sum((x-est.hat[1])*(y-est.hat[2]))/n #cov(x,y)
  xy.cov <- matrix(c(x.var.est, xy.cov.est, xy.cov.est, y.var.est), nrow = 2)
  single.cov <- xy.cov[1,2]
  Chi.stat <- est.hat%*%solve(xy.cov)%*%t(est.hat)
  Z <- est.hat/sqrt(c(xy.cov[1,1], xy.cov[2,2]))
  p.vals.sing <- 2*ifelse(Z > 0, pnorm(Z, lower.tail = FALSE), pnorm(Z, lower.tail = TRUE))
  p.val.comb <- pchisq(Chi.stat, 2, lower.tail = FALSE)
  return(list(betaL1 = est.hat[1], betaL1.Z = Z[1], betaL1.pValue = p.vals.sing[1], 
              betaC1 = est.hat[2], betaC1.Z = Z[2], betaC1.pValue = p.vals.sing[2], 
              Chi2.stat = Chi.stat, Chi2.pValue = p.val.comb, adj.Chi2.pValue = NA))
}


#' Prints scREhurdle.fit object
#' @export
#' @noRd
`print.scREhurdle.fit` <-
  function(x, ...){
    print(round(x$deTab[,c(1,4,7,8,9)],6))
  }

#' Summary for scREhurdle.fit object
#' @export
#' @noRd
`summary.scREhurdle.fit` <-
  function(object, ...){
    object.class <- c(class(object$stanFit), "-none-",class(object$deTab),class(object$geneTab),class(object$cellTab),class(object$hyperEst))
    object.length <- c(length(object$stanFit), length(object$inputData), length(object$deTab),length(object$geneTab),length(object$cellTab),length(object$hyperEst))
    object.mode <- c(storage.mode(object$stanFit), storage.mode(object$inputData), storage.mode(object$deTab),storage.mode(object$geneTab),storage.mode(object$cellTab),
                     storage.mode(object$hyperEst))
    object.frame <- data.frame(Class = object.class, Length = object.length, Mode = object.mode)
    rownames(object.frame) <- attributes(object)$names
    print(object.frame)
  }



#' Toy scRNA-seq dataset and corresponding clustering assignment
#'
#' toyDat contains two objects: scDat and subpop.
#' The data.frame called "scDat" contains gene expression counts. 25 cells belong to the control (C) group and 25 cells belong to treatment (T) group as indicated
#' in the column names.
#' The object "subpop" contains the subpopulation clustering assignment for the cells. The control group is clustered into 2 subpopulations and the treatment 
#' group is clustered into 3 supopulations.
#'
#'
#' @format A list containing the following objects:
#'  \itemize{
#'   \item{scDat}{: data.frame with 100 rows (genes) and 50 columns (cells)}
#'   \item{subpop}{: subpopulation clustering assignment for the 50 cells}
#' }
#' 
#' @examples
#' data(toyDat)
#' # Obtain dataset
#' toyDat$scData
#' 
#' # Obtain clustering assignment
#' toyDat$subpop
"toyDat"
