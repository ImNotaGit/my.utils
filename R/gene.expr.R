## ----functions for processing gene expression data----


get.major.genes <- function(x) {
  # given a character vector of gene identifiers (symbols) in x, remove:
  # genes w/o formal names (LOC, OTTHUMG)
  # miRNA genes (MIR)
  # long non-protein-coding RNA genes (LINC, -AS, -IT, -OT)
  # protein-coding genes of unknown function on the opposite strand (-OS)
  # return a logical vector indicating the genes to be kept
  !grepl("^(---|LOC|OTTHUMG|MIR|LINC)", x) & !grepl("-(OS|AS|IT|OT)[0-9]*$", x)
}


rm.batch.eff <- function(...) {
  # remove the known batch effect in expression data. should pass in arbitrary numbers of gene-by-sample expression matrices. return a list of batch effect-corrected data, in the original order.

  tmp <- list(...)
  combined.dat <- cbind(...)
  batch.id <- rep(1:length(tmp), sapply(tmp, ncol))
  res <- sva::ComBat(combined.dat, batch.id)
  lapply(1:length(tmp), function(i) res[, batch.id==i])
}


prep.data <- function(dat, log="default", norm.method="loess") {
  # prepare expression data: transformation, normalization, etc.
  # dat: either a matrix (of gene-by-sample) or an ExpressionSet object
  # log: whether to log2(x+1) data; "default" automatically determines whether to transform
  # norm.method: either "loess" or "quantile"; "loess" won't work if the data contains NA

  if (class(dat)=="ExpressionSet") {
    # for featureData, fix potential improper column names so that later limma::topTable can use them
    fvarLabels(dat) <- make.names(fvarLabels(dat))
    mat <- exprs(dat)
  } else if (is.matrix(dat)) mat <- dat

  # log2 transform
  if (log=="default") {
    # following the method as in GEO2R, for microarray data
    qx <- as.numeric(quantile(mat, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm=TRUE))
    log <- (qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  }
  if (log) {
    nlt0 <- sum(mat<0, na.rm=TRUE)
    if (nlt0>0) {
      warning(sprintf("Log-transformation error: there are %d negative values in the data, the data may be already on log-scale.\nHere is a summary of the data:\n", nlt0))
      print(summary(as.vector(mat)))
      stop()
    }
    mat <- log2(mat+1)
    cat("log2-transformation performed.\n")
  } else cat("log2-transformation NOT performed.\n")

  # normalization
  if (norm.method=="loess") {
    nna <- sum(is.na(mat))
    if (nna>0) {
      stop(sprintf("Loess normalization error: there are %d NA/NaN's in the data.\n", nna))
    } else mat <- affy::normalize.loess(mat, log.it=FALSE)
  } else if (norm.method=="quantile") {
    mat <- limma::normalizeQuantiles(mat)
  } else cat("Normalization NOT performed.\n")

  # return
  if (class(dat)=="ExpressionSet") {
    exprs(dat) <- mat
    return(dat)
  } else if (is.matrix(dat)) return(mat)
}


de <- function(dat, pheno, model="~.", coef, robust=FALSE, trend=FALSE) {
  # differential expression analysis
  # dat: either a matrix (of gene-by-sample) or an ExpressionSet object
  # pheno: phenotypic data as a data.frame with the same order of samples
  # model: the linear model to use for DE, by default a linear model containing all variables in pheno
  # coef: character, the name of the variable (and its level, if categorical) of interest for which the linear model coefficients to be displayed, e.g. if there's a variable named "gender" with two levels "female" and "male" with "female" being the reference level, then we may use coef="gendermale"
  # robust, trend: parameters for limma eBayes, set to TRUE for log-transformed RNA-seq data

  if (class(dat)=="ExpressionSet") {
    mat <- exprs(dat)
    rownames(mat) <- fData(dat)$Gene.symbol
  } else if (is.matrix(dat)) mat <- dat

  design <- model.matrix(as.formula(model), pheno)
  fit <- limma::lmFit(mat, design)
  fit <- limma::eBayes(fit, robust=robust, trend=trend)
  res <- as.data.table(topTable(fit, coef=coef, number=Inf, genelist=rownames(mat)))
  setnames(res, c("id","log.fc","ave.expr","t","pval","padj","B"))
  res
}

