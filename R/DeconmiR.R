#' @title Cell type deconvolution from miRNA expression profiling
#'
#' @description Cell type proportions prediction from miRNA expression profiling
#' data by applying multiple deconvolution methods. Inference proceeds via one
#' of 5 methods (Robust Partial Correlations-RPC, Support Vector Regression-SVR,
#' Constrained Projection-CP, Non-Negative Least Squares-NNLS, Solution of Linear
#' Equations-SLE), as determined by the user.
#'
#' @param data.m miRNA Expression profile: A data matrix with rows labeling the
#' miRNAs (should use miRBase ID as in ref.m) and columns labeling samples
#' (e.g. primary tumor specimens). Missing value is not allowed and all values
#' should be log-transformed.
#'
#' @param ref.m Reference profile: The signature matrix should be a miRNAs (rows)
#' by cell types (columns) matrix containing log-normalized expression
#' values of signature miRNAs.
#'
#' @param method Deconvolution method: Choice of the deconvolution methods to
#' be used: ("RPC", "SVR", "CP", "NNLS", "SLE"). By default is the "RPC".
#'
#' @param maxit The limit of the number of IWLS iterations, only used for "RPC" mode.
#'
#' @param nu.v A vector of several candidate nu values, which is needed for
#' nu-classification, nu-regression, and one-classification in svm, only used
#' for "SVR" mode. The best estimation results among all candidate nu will be
#' automatically returned.
#'
#' @param constraint Choice of either 'inequality' or 'equality' normalization
#' constraint, only used in "CP" mode, By default is 'inequality' (i.e sum of
#' weights adds to a number less or equal than 1), which was implemented in
#' Houseman et al (2012).
#'
#' @return A list with the following entries: estF: a matrix of the estimated
#' fractions, rows are samples and columns are cell types; ref: the reference
#' matrix used; dataREF: the subset of the input data matrix with only the
#' miRNAs defined in the reference matrix.
#'
#' @export
#'
#' @examples
#' data(miRExpCount_LAML.m)
#' data(refMatrix_Blood1.m)
#' miRexpRPM_LAML.m <- getRPM(miRExpCount_LAML.m)
#' miRexpNorm_LAML.m <- logNormalize(miRexpRPM_LAML.m)
#' est.o <- DeconmiR(miRexpNorm_LAML.m, refMatrix_Blood1.m, method = 'RPC')
#' cellfrac.m <- est.o$estF
#'
DeconmiR <- function(data.m, ref.m, method = c("RPC", "SVR", "CP", "NNLS", "SLE"),
  maxit = 50, nu.v = c(0.25, 0.5, 0.75), constraint = c("inequality", "equality")) {

  method <- match.arg(method)
  constraint <- match.arg(constraint)
  if (!method %in% c("RPC", "SVR", "CP", "NNLS", "SLE"))
    stop("Input a valid method!")
  if (method == "RPC") {
    out.o <- DoRPC(data.m, ref.m, maxit)
  } else if (method == "SVR") {
    out.o <- DoSVR(data.m, ref.m, nu.v)
  } else if (method == "CP") {
    if (!constraint %in% c("inequality", "equality")) {
      # make sure constraint must be inequality or equality
      stop("constraint must be inequality or equality when using CP!")
    } else out.o <- DoCP(data.m, ref.m, constraint)
  }else if (method == "NNLS") {
    out.o <- DoNNLS(data.m, ref.m)
  }else if (method == "SLE") {
    out.o <- DoSLE(data.m, ref.m)
  }
  return(out.o)
}

#### SLE: Solution of Linear Equations
DoSLE <- function(data.m, ref.m){
  map.idx <- match(rownames(ref.m),rownames(data.m));
  rep.idx <- which(is.na(map.idx)==FALSE);

  data2.m <- data.m[map.idx[rep.idx],];
  ref2.m <- ref.m[rep.idx,];

  est.m <- matrix(nrow=ncol(data2.m),ncol=ncol(ref2.m));
  colnames(est.m) <- colnames(ref2.m);
  rownames(est.m) <- colnames(data2.m);

  for (s in 1:ncol(data2.m)){
    A <- ref2.m
    b <- data2.m[,s]
    coef.v <- solve(t(A) %*% A) %*% (t(A) %*% b)
    coef.v[which(coef.v < 0)] <- 0
    coefSum <- sum(coef.v)
    coef.v <- coef.v/sum(coef.v)
    est.m[s,] = coef.v
  }
  return(list(estF=est.m,ref=ref2.m,dataREF=data2.m));
}

#' @importFrom nnls nnls
#### NNLS: Non-negative Least Squares
DoNNLS <- function(data.m, ref.m){
  map.idx <- match(rownames(ref.m),rownames(data.m));
  rep.idx <- which(is.na(map.idx)==FALSE);

  data2.m <- data.m[map.idx[rep.idx],];
  ref2.m <- ref.m[rep.idx,];

  est.m <- matrix(nrow=ncol(data2.m),ncol=ncol(ref2.m));
  colnames(est.m) <- colnames(ref2.m);
  rownames(est.m) <- colnames(data2.m);

  for (s in 1:ncol(data2.m)){
    coef.v <- nnls(ref2.m, data2.m[,s])[[1]]
    coefSum <- sum(coef.v)
    coef.v <- coef.v/sum(coef.v)
    est.m[s,] = coef.v
  }
  return(list(estF=est.m,ref=ref2.m,dataREF=data2.m));
}

#' @importFrom MASS rlm
### RPC: Robust partial correlation
DoRPC <- function(avdata.m,ref.m,maxit){
  map.idx <- match(rownames(ref.m),rownames(avdata.m));
  rep.idx <- which(is.na(map.idx)==FALSE);

  data2.m <- avdata.m[map.idx[rep.idx],];
  ref2.m <- ref.m[rep.idx,];

  est.m <- matrix(nrow=ncol(data2.m),ncol=ncol(ref2.m));
  colnames(est.m) <- colnames(ref2.m);
  rownames(est.m) <- colnames(data2.m);
  for(s in 1:ncol(data2.m)){
    rlm.o <- rlm( data2.m[,s] ~ ref2.m, maxit=maxit)
    coef.v <- summary(rlm.o)$coef[2:(ncol(ref2.m)+1),1];
    coef.v[which(coef.v<0)] <- 0;
    total <- sum(coef.v);
    coef.v <- coef.v/total
    est.m[s,] <- coef.v;
  }
  return(list(estF=est.m,ref=ref2.m,dataREF=data2.m));
}

#' @importFrom e1071 svm
### SVR: support vector regression/CIBERSORT
DoSVR<- function(avdata.m,ref.m,nu.v=c(0.25,0.5,0.75)){
  map.idx <- match(rownames(ref.m),rownames(avdata.m));
  rep.idx <- which(is.na(map.idx)==FALSE);

  data2.m <- avdata.m[map.idx[rep.idx],];
  ref2.m <- ref.m[rep.idx,];

  est.lm <- list();
  nui <- 1;
  for(nu in nu.v){
    est.m <- matrix(nrow=ncol(data2.m),ncol=ncol(ref2.m));
    colnames(est.m) <- colnames(ref2.m);
    rownames(est.m) <- colnames(data2.m);
    for(s in 1:ncol(data2.m)){
      svm.o <- svm(x=ref2.m,y=data2.m[,s],scale = TRUE, type="nu-regression", kernel ="linear", nu = nu);
      coef.v <- t(svm.o$coefs) %*% svm.o$SV;
      coef.v[which(coef.v<0)] <- 0;
      total <- sum(coef.v);
      coef.v <- coef.v/total;
      est.m[s,] <- coef.v;
    }
    est.lm[[nui]] <- est.m;
    print(nui);
    nui <- nui+1;
  }

  #### select best nu
  rmse.m <- matrix(NA,nrow=ncol(avdata.m),ncol=length(nu.v));
  for(nui in 1:length(nu.v)){
    reconst.m <- ref2.m %*% t(est.lm[[nui]]);
    for(s in 1:ncol(avdata.m)){
      rmse.m[s,nui] <- sqrt(mean((data2.m[,s] - reconst.m[,s])^2));
    }
    print(nui);
  }
  colnames(rmse.m) <- nu.v;
  nu.idx <- apply(rmse.m,1,which.min);
  estF.m <- est.m;
  for(s in 1:nrow(estF.m)){
    estF.m[s,] <- est.lm[[nu.idx[s]]][s,];
  }
  return(list(estF=estF.m,nu=nu.v[nu.idx],ref=ref2.m,dataREF=data2.m));
}

#' @importFrom quadprog solve.QP
### Houseman CP/quadratic programming
DoCP <- function(beta.m, ref.m, constraint) {
  ### define D matrix
  nCT <- ncol(ref.m)
  D <- 2 * apply(ref.m, 2, function(x) colSums(x * ref.m))

  ### for inequality and equality, coe.v is different
  if (constraint == "inequality") {
    coe.v <- c(-1, 0)
  } else coe.v <- c(1, 1)

  ### define constraints
  A.m <- matrix(0, nrow = nCT, ncol = nCT)
  diag(A.m) <- rep(1, nCT)
  A.m <- cbind(rep(coe.v[1], nCT), A.m)
  b0.v <- c(coe.v[1], rep(0, nCT))

  ### define d-vector and solve for each sample
  nS <- ncol(beta.m)
  westQP.m <- matrix(NA, ncol = ncol(ref.m), nrow = nS)
  colnames(westQP.m) <- colnames(ref.m)
  rownames(westQP.m) <- colnames(beta.m)

  map.idx <- match(rownames(ref.m), rownames(beta.m))
  rep.idx <- which(is.na(map.idx) == FALSE)
  for (s in seq_len(nS)) {
    tmp.v <- beta.m[, s]
    d.v <- as.vector(2 * matrix(tmp.v[map.idx[rep.idx]], nrow = 1) %*% ref.m[rep.idx,
    ])
    qp.o <- solve.QP(D, d.v, A.m, b0.v, meq = coe.v[2])
    westQP.m[s, ] <- qp.o$sol
    message(s)
  }

  return(list(estF = westQP.m, ref = ref.m[rep.idx, ], dataREF = beta.m[map.idx[rep.idx],
  ]))
}



