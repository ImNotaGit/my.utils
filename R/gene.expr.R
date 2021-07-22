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

get.gene.lengths <- function(txdb=NULL, fn=NULL, format=c("auto","gff3","gtf")) {
  # compute length of each gene give a gtf or gff file or a TxDb object
  # txdb: txdb object; fn: file name; format: format of file
  # use summed length of non-overlapping exons per gene (the "union" method)

  if (is.null(txdb)) {
  	if (is.null(fn)) stop("Need to provide either exdb or fn.")
  	format <- match.arg(format)
    txdb <- GenomicFeatures::makeTxDbFromGFF(fn, format=format)
  }
  
  # exons by gene
  exons.by.gene <- GenomicFeatures::exonsBy(txdb, by="gene")
  # summed length of non-overlapping exons per gene
  sum(IRanges::width(GenomicRanges::reduce(exons.by.gene)))
}


rm.batch.eff <- function(...) {
  # remove the known batch effect in expression data. should pass in arbitrary numbers of gene-by-sample expression matrices. return a list of batch effect-corrected data, in the original order.

  tmp <- list(...)
  combined.dat <- cbind(...)
  batch.id <- rep(1:length(tmp), sapply(tmp, ncol))
  res <- sva::ComBat(combined.dat, batch.id)
  lapply(1:length(tmp), function(i) res[, batch.id==i])
}


prep.array <- function(dat, log="default", norm.method="loess") {
  # prepare microarray gene expression data: transformation, normalization, etc.
  # dat: either a matrix (of gene-by-sample) or an ExpressionSet object
  # log: whether to log2(x+1) data; "default" automatically determines whether to transform
  # norm.method: either "loess" or "quantile"; "loess" won't work if the data contains NA, so if any NA will first set these to 0 with warning

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

prep.data <- prep.array


voom <- function(dat, pheno=NULL, model=~., design=NULL, quantile=FALSE, ...) {
  # perform limma::voom normalization for RNAseq data
  # dat: gene-by-sample expression matrix of raw counts; should have low genes already filtered out
  # pheno: phenotypic data as a data.frame with the same order of samples; model: the model to be used for downstream DE; together with pheno, this will be used to generate the design matrix passed to limma::voom
  # or provide design matrix in design; if both pheno and design are NULL, then design will be NULL
  # quantile: whether to apply quantile normalization, if TRUE, will pass normalize.method="quantile" to limma::voom
  # ...: passed to limma::voom
  # return a list(voom, design, genes), where voom is the limma::voom() output, i.e. an EList object, design is the design matrix, and genes is rownames(dat)

  if (!is.null(pheno)) {
    vs <- unique(c(all.vars(model), names(model.frame(model, pheno))))
    vs <- vs[vs!="." & !grepl("\\(|\\)", vs)]
    ccs <- complete.cases(pheno[, vs, with=FALSE])
    if (any(!ccs)) message("Removed ", sum(!ccs), " samples with incomplete (NA) covariate data.")
    dat <- dat[, ccs]
    phe <- pheno[ccs]
    design <- model.matrix(model, phe)
  }
  gs <- rownames(dat)
  
  dat <- edgeR::DGEList(counts=dat)
  dat <- edgeR::calcNormFactors(dat)
  if (quantile) nm <- "quantile" else nm <- "none"
  v <- limma::voom(counts=dat, design=design, normalize.method=nm, ...)
  res <- list(voom=v, design=design, genes=gs)
}

de.limma <- function(dat, pheno=NULL, model=~., design=NULL, coef, robust=FALSE, trend=FALSE, gene.colname=TRUE) {
  # differential expression analysis with limma
  # dat: either a matrix (of gene-by-sample) or an ExpressionSet object, or output from voom()
  # pheno: phenotypic data as a data.frame with the same order of samples; model: the linear model to use for DE, by default a linear model containing all variables in pheno (w/o interaction terms); these two will be used to compute the design matrix
  # or provide design matrix in design, pheno and design cannot both be NULL, unless if dat is output from voom() then will use the design matrix saved in there, and pheno, model and design will be ignored
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
    rownames(mat)[is.na(rownames(mat))] <- ""
    gns[is.na(gns)] <- ""
  } else if (is.matrix(dat)) {
  	mat <- dat
  	rownames(mat)[is.na(rownames(mat))] <- ""
  	gns <- rownames(mat)
  } else if (is.list(dat) && "voom" %in% names(dat) && class(dat$voom)=="EList") {
    # if dat is output from voom()
    mat <- dat$voom
    gns <- dat$genes
    design <- dat$design
    if (!is.null(design)) message("Using the design matrix saved in `dat`, ignoring `pheno` and `model`, if provided.")
  }

  if (is.null(design)) {
    if (is.null(pheno)) stop("Need to provide either `pheno` with `model`, or `design`.")
    vs <- unique(c(all.vars(model), names(model.frame(model, pheno))))
    vs <- vs[vs!="." & !grepl("\\(|\\)", vs)]
    ccs <- complete.cases(pheno[, vs, with=FALSE])
    if (any(!ccs)) message("Removed ", sum(!ccs), " samples with incomplete (NA) covariate data.")
    if (is.matrix(mat)) mat <- mat[, ccs]
    phe <- pheno[ccs]
    design <- model.matrix(model, phe)
  }
  
  fit <- limma::lmFit(mat, design=design)
  fit <- limma::eBayes(fit, robust=robust, trend=trend)
  res <- tryCatch(limma::topTable(fit, coef=coef, number=Inf, genelist=gns),
                  error=function(e) limma::topTable(fit, coef=coef, number=Inf))
  if (!"id" %in% tolower(names(res))) res <- cbind(ID=row.names(res), res)
  res <- as.data.table(res)
  setnames(res, c("ID","logFC","AveExpr","P.Value","adj.P.Val"), c("id","log.fc","ave.expr","pval","padj"), skip_absent=TRUE)
  res
}

de <- de.limma

rm.low.genes <- function(dat, rm.low.frac.gt=0.5, count.cutoff=10) {
  # remove genes with low count
  # dat: gene-by-sample expression matrix of raw counts, or an "eset" i.e. list(expr, pheno, geneid)

  is.eset <- is.list(dat) && all(c("expr","pheno") %in% names(dat))
  if (is.eset) {
    dat0 <- dat
    dat <- dat$expr
  }

  tot.cnts <- colSums(dat)
  cutoff <- count.cutoff/median(tot.cnts)*1e6
  dat1 <- edgeR::cpm(dat, lib.size=tot.cnts)
  keep <- rowMeans(dat1<cutoff) <= rm.low.frac.gt
  message(sum(keep), " rows (genes/transcripts) remaining.")
  if (is.eset) res <- select.genes.eset(dat0, keep) else res <- dat[keep, ]
  res
}


get.tmm.log.cpm <- function(dat, prior.count=1) {
  # get log2(cpm+1) values with edgeR (TMM-normalized), from raw counts
  # dat: gene-by-sample expression matrix of raw counts; should have low genes already filtered out

  dge <- edgeR::DGEList(counts=dat)
  dge <- edgeR::calcNormFactors(dge)
  edgeR::cpm(dge, log=TRUE, prior.count=prior.count)
}


de.edger <- function(dat, pheno=NULL, model=~., design=NULL, coef, lfc.cutoff=0) {
  # differential expression analysis with edgeR
  # dat: gene-by-sample expression matrix of raw counts; should have low genes already filtered out
  # pheno: phenotypic data as a data.frame with the same order of samples
  # model: the model to use for DE, by default a linear model containing all variables in pheno (w/o interaction terms)
  # pheno and model will be used to compute the design matrix; or provide design matrix in design; pheno and design cannot be both NULL
  # coef: character, the name of the variable (and its level, if categorical) of interest for which the linear model coefficients to be displayed, e.g. if there's a variable named "gender" with two levels "female" and "male" with "female" being the reference level, then we may use coef="gendermale"

  if (is.null(design)) {
    if (is.null(pheno)) stop("Need to provide either `pheno` with `model`, or `design`.")
    vs <- unique(c(all.vars(model), names(model.frame(model, pheno))))
    vs <- vs[vs!="." & !grepl("\\(|\\)", vs)]
    ccs <- complete.cases(pheno[, vs, with=FALSE])
    if (any(!ccs)) message("Removed ", sum(!ccs), " samples with incomplete (NA) covariate data.")
    dat <- dat[, ccs]
    phe <- pheno[ccs]
    design <- model.matrix(model, phe)
  }
  
  dge <- edgeR::DGEList(counts=dat)
  dge <- edgeR::calcNormFactors(dge)
  dge <- edgeR::estimateDisp(dge, design)

  fit <- edgeR::glmQLFit(dge, design)
  if (lfc.cutoff==0) test.res <- edgeR::glmQLFTest(fit, coef=coef) else test.res <- edgeR::glmTreat(fit, coef=coef, lfc=lfc.cutoff)
  res <- edgeR::topTags(test.res, n=Inf)[[1]]
  res <- cbind(id=row.names(res), as.data.table(res))
  setnames(res, c("logFC","logCPM","PValue","FDR"), c("log.fc","log.cpm","pval","padj"), skip_absent=TRUE)
  res
}


de.deseq2 <- function(dat, pheno=NULL, model=~., design=NULL, coef, ...) {
  # differential expression analysis with DESeq2
  # dat: gene-by-sample expression matrix of raw counts; should have integer counts (if not then rounded) and have low genes already filtered out
  # pheno: phenotypic data as a data.table with the same order of samples
  # model: the model to use for DE, by default a linear model containing all variables in pheno (w/o interaction terms)
  # pheno and model will be used to compute the design matrix; or provide design matrix in design; pheno and design cannot be both NULL
  # coef: character, the name of the variable (and its level, if categorical) of interest for which the linear model coefficients to be displayed, e.g. if there's a variable named "gender" with two levels "female" and "male" with "female" being the reference level, then we may use coef="gendermale"
  # coef is passed to the name argument (if coef has length 1) or the contrasts argument (if coef has length >1) for DESeq2::results; for the latter, e.g. coef=c("group", "trt", "ctrl") will return results for the 'trt' level compared to 'ctrl' level of the `group` variable
  # ...: passed to DESeq2::results
  
  if (any(!is.wholenumber(dat))) {
    message("DESeq2 only accepts integer read counts. There are non-integer values in `dat` and they have been rounded.")
    dat <- round(dat)
  }
  
  if (is.null(design)) {
    if (is.null(pheno)) stop("Need to provide either `pheno` with `model`, or `design`.")
    pheno <- copy(as.data.table(pheno))
    vs <- unique(c(all.vars(model), names(model.frame(model, pheno))))
    vs <- vs[vs!="." & !grepl("\\(|\\)", vs)]
    ccs <- complete.cases(pheno[, vs, with=FALSE])
    if (any(!ccs)) message("Removed ", sum(!ccs), " samples with incomplete (NA) covariate data.")
    dat <- dat[, ccs]
    phe <- pheno[ccs]
    tmp <- sapply(phe[, vs, with=FALSE], function(x) !is.numeric(x) & !is.factor(x))
    if (any(tmp)) {
      message("These non-numeric variables included in the model are not factors:")
      message(cc(vs[tmp]))
      message("They are converted to factors.")
      phe[, c(vs[tmp]):=lapply(.SD, factor), .SDcols=vs[tmp]]
    }
    design <- model.matrix(model, phe)
  }
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=dat, colData=phe, design=design)
  dds1 <- DESeq2::DESeq(dds)
  if (length(coef)==1) de.res <- DESeq2::results(dds1, name=coef, ...) else de.res <- DESeq2::results(dds1, contrast=coef, ...)
  gid <- rownames(de.res)
  de.res <- as.data.table(de.res)[, .(id=gid, ave.expr=baseMean, log.fc=log2FoldChange, lfc.se=lfcSE, pval=pvalue, padj=padj)]
}
