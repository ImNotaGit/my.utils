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
  # prepare microarray gene expression data: transformation, normalization, etc.
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
      mat[mat<0] <- 0
      warning("There are ",nlt0," negative values in the data, these are set to 0 for log-transformation.")
    }
    mat <- log2(mat+1)
    message("log2-transformation performed.")
  } else message("log2-transformation NOT performed.")

  # normalization
  if (norm.method=="loess") {
    nna <- sum(is.na(mat))
    if (nna>0) {
      mat[is.na(mat)] <- 0
      warning("There are ",nna," NA/NaN's in the data, these are set to 0 for loess normalization.")
    }
    mat <- affy::normalize.loess(mat, log.it=FALSE)
  } else if (norm.method=="quantile") {
    mat <- limma::normalizeQuantiles(mat)
  } else message("Normalization NOT performed.")

  # return
  if (class(dat)=="ExpressionSet") {
    exprs(dat) <- mat
    return(dat)
  } else if (is.matrix(dat)) return(mat)
}


de <- function(dat, pheno, model="~.", coef, robust=FALSE, trend=FALSE, gene.colname=TRUE) {
  # differential expression analysis with limma
  # dat: either a matrix (of gene-by-sample) or an ExpressionSet object
  # pheno: phenotypic data as a data.frame with the same order of samples
  # model: the linear model to use for DE, by default a linear model containing all variables in pheno (w/o interaction terms)
  # coef: character, the name of the variable (and its level, if categorical) of interest for which the linear model coefficients to be displayed, e.g. if there's a variable named "gender" with two levels "female" and "male" with "female" being the reference level, then we may use coef="gendermale"
  # robust, trend: parameters for limma eBayes, set to TRUE for log-transformed RNA-seq data (but for RNA-seq doing DE with read count data using other methods is recommended)
  # gene.colname: the column name for gene symbols in fData(dat) if dat is an ExpressionSet; if TRUE then will try to get gene symbols automatically from fData(dat); if FALSE then will not try to get gene symbols; gene symbols will be added to the returned DE table

  if (class(dat)=="ExpressionSet") {
    mat <- exprs(dat)
    if (!isFALSE(gene.colname)) {
      if (is.character(gene.colname)) {
      	if (!gene.colname %in% names(fData(dat))) stop("`",gene.colname,"` not in fData(dat).") else idx <- gene.colname
      } else {
        idx <- which(tolower(names(fData(dat))) %in% c("gene.symbol","gene_symbol","gene symbol","symbol","orf"))
        if (length(idx)==0) stop("Gene symbol annotation not found in fData(dat), please check.")
        idx <- idx[1]
        message("Used the column `",names(fData(dat))[idx],"` in fData(dat) as gene symbols.")
      }
      gns <- fData(dat)[[idx]]
    } else gns <- rownames(mat)
  } else if (is.matrix(dat)) mat <- dat

  design <- model.matrix(as.formula(model), pheno)
  fit <- limma::lmFit(mat, design)
  fit <- limma::eBayes(fit, robust=robust, trend=trend)
  res <- tryCatch(as.data.table(limma::topTable(fit, coef=coef, number=Inf, genelist=gns)),
                  error=function(e) as.data.table(limma::topTable(fit, coef=coef, number=Inf)))
  setnames(res, "ID", "id")
  setnames(res, "logFC", "log.fc")
  setnames(res, "AveExpr", "ave.expr")
  setnames(res, "P.Value", "pval")
  setnames(res, "adj.P.Val", "padj")
  res
}


rm.low.genes <- function(dat, rm.low.frac.gt=0.5, count.cutoff=10) {
  # remove genes with low count
  # dat: gene-by-sample expression matrix of raw counts

  tot.cnts <- colSums(dat)
  cutoff <- count.cutoff/median(tot.cnts)*1e6
  dat1 <- edgeR::cpm(dat, lib.size=tot.cnts)
  keep <- rowMeans(dat1<cutoff) <= rm.low.frac.gt
  message(sum(keep), " rows (genes/transcripts) remaining.")
  dat[keep, ]
}


get.tmm.log.cpm <- function(dat) {
  # get log2(cpm+1) values with edgeR (TMM-normalized), from raw counts
  # dat: gene-by-sample expression matrix of raw counts; should have low genes already filtered out

  dge <- edgeR::DGEList(counts=dat)
  dge <- edgeR::calcNormFactors(dge)
  edgeR::cpm(dge, log=TRUE, prior.count=1)
}


de.edger <- function(dat, pheno, model="~.", coef, lfc.cutoff=0) {
  # differential expression analysis with edgeR
  # dat: gene-by-sample expression matrix of raw counts; should have low genes already filtered out
  # pheno: phenotypic data as a data.frame with the same order of samples
  # model: the model to use for DE, by default a linear model containing all variables in pheno (w/o interaction terms)
  # coef: character, the name of the variable (and its level, if categorical) of interest for which the linear model coefficients to be displayed, e.g. if there's a variable named "gender" with two levels "female" and "male" with "female" being the reference level, then we may use coef="gendermale"

  dge <- edgeR::DGEList(counts=dat)
  dge <- edgeR::calcNormFactors(dge)
  design <- model.matrix(as.formula(model), pheno)
  dge <- edgeR::estimateDisp(dge, design)

  fit <- edgeR::glmQLFit(dge, design)
  if (lfc.cutoff==0) test.res <- edgeR::glmQLFTest(fit, coef=coef) else test.res <- edgeR::glmTreat(fit, coef=coef, lfc=lfc.cutoff)
  res <- edgeR::topTags(test.res, n=Inf)[[1]]
  res <- cbind(id=row.names(res), as.data.table(res))
  setnames(res, "logFC", "log.fc")
  setnames(res, "logCPM", "log.cpm")
  setnames(res, "PValue", "pval")
  setnames(res, "FDR", "padj")
  res
}

