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


rm.low.genes.sr <- function(dat, rm.low.frac.gt=0.5, low.cutoff=0, assay="RNA", slot=c("counts","data","scale.data"), renormalize=TRUE, ...) {
  # remove genes with low expression from a Seurat object, based on data in the specified assay and slot
  # renormalize: whether to renormalize the assay with Seurat::NormalizeData
  # ...: passed to Seurat::NormalizeData

  if (!requireNamespace("Seurat", quietly=TRUE)) {
    stop("Package \"Seurat\" needed for this function to work.")
  }
  library(Seurat)

  slot <- match.arg(slot)
  keep <- Matrix::rowMeans(GetAssayData(dat, assay=assay, slot=slot)<=low.cutoff) <= rm.low.frac.gt
  message(sum(keep), " rows (genes/transcripts) remaining.")
  res <- subset(dat, features=rownames(dat[[assay]])[keep])
  if (renormalize) res <- NormalizeData(res, ...)
  res
}


get.tmm.log.cpm <- function(dat, prior.count=1) {
  # get log2(cpm+1) values with edgeR (TMM-normalized), from raw counts
  # dat: gene-by-sample expression matrix of raw counts; should have low genes already filtered out

  dge <- edgeR::DGEList(counts=dat)
  dge <- edgeR::calcNormFactors(dge)
  edgeR::cpm(dge, log=TRUE, prior.count=prior.count)
}


de.edger <- function(dat, pheno=NULL, model=~., design=NULL, coef, lfc.cutoff=0, robust=FALSE, ...) {
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

  fit <- edgeR::glmQLFit(dge, design, robust=robust, ...)
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
  } else {
    if (is.null(pheno)) phe <- data.table(idx=1:ncol(dat)) else phe <- pheno
  }
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=dat, colData=phe, design=design)
  dds1 <- DESeq2::DESeq(dds)
  if (length(coef)==1) de.res <- DESeq2::results(dds1, name=coef, ...) else de.res <- DESeq2::results(dds1, contrast=coef, ...)
  gid <- rownames(de.res)
  de.res <- as.data.table(de.res)[, .(id=gid, ave.expr=baseMean, log.fc=log2FoldChange, lfc.se=lfcSE, pval=pvalue, padj=padj)][order(padj, pval)]
}


de.gampoi <- function(dat, pheno, model=~., design, coef, size.factors, pseudobulk=NULL, return.fit=FALSE, ...) {
  # differential expression analysis with glmGamPoi, this may be used for single-cell RNA-seq data
  # dat: gene-by-sample expression matrix of log-normalized expression value, with gene ID/symbol as rownames and sample ID/barcode as colnames
  # pheno: phenotypic data as a data.table with the same order of samples
  # model: the model to use for DE, by default a linear model containing all variables in pheno (w/o interaction terms)
  # pheno and model will be used to compute the design matrix; or provide design matrix in design; pheno and design cannot both be missing; the design matrix should have proper column names
  # coef: a single character passed to glmGamPoi::test_de `contrast`, e.g. can be the name of the variable (and its level, if categorical) of interest for which the linear model coefficients to be displayed, e.g. if there's a variable named "group" with two levels "control" and "treated" with "control" being the reference level, then we may use coef="grouptreated", corresponding to the result of comparing treated to control group; or a list of two contrast vectors, each named by the colnames of the design matrix, the first corresponds to the baseline (e.g. control), the second corresponds to the group of interest (e.g. treated)
  # size.factors: passed to glmGamPoi::glm_gp `size_factors`, if missing use the default "normed_sum"
  # pseudobulk: passed to glmGamPoi::test_de `pseudobulk_by`
  # return.fit: if TRUE, return list(fit=fit, de.res=de.res), where fit is the output from glmGamPoi::glm_gp; otherwise return de.res
  # ...: passed to glmGamPoi::glm_gp
  
  if (!requireNamespace("glmGamPoi", quietly=TRUE)) {
    stop("Package \"glmGamPoi\" needed for this function to work.")
  }

  if (missing(size.factors)) size.factors <- "normed_sum"

  if (missing(design)) {
    if (missing(pheno)) stop("Need to provide either `pheno` with `model`, or `design`.")
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
  
  fit <- glmGamPoi::glm_gp(as.matrix(dat), design=design, size_factors=size.factors, ...)
  de.res <- as.data.table(glmGamPoi::test_de(fit, contrast=coef, pseudobulk_by=pseudobulk))[order(adj_pval, pval), .(id=name, log.fc=lfc, F=f_statistic, pval, padj=adj_pval)]

  if (return.fit) list(fit=fit, de.res=de.res) else de.res
}


de.mast <- function(dat, pheno, model=~., design, cdr=TRUE, coef, lfc.cutoff=0, pos.only=FALSE, lfc.only=FALSE, nc=1L, return.fit=FALSE, ...) {
  # differential expression analysis with MAST, used for e.g. single-cell RNA-seq data
  # dat: gene-by-sample expression matrix (assuming sparse) of log-normalized expression value, with gene ID/symbol as rownames and sample ID/barcode as colnames
  # pheno: phenotypic data as a data.table with the same order of samples
  # model: the model to use for DE, by default a linear model containing all variables in pheno (w/o interaction terms)
  # pheno and model will be used to compute the design matrix; or provide design matrix in design; pheno and design cannot both be missing; the design matrix should have proper column names
  # cdr: whether to include cellular detection rate (i.e. fraction of >0 genes in each cell) as a covariate
  # coef: a single character, the name of the variable (and its level, if categorical) of interest for which the linear model coefficients to be displayed, e.g. if there's a variable named "group" with two levels "control" and "treated" with "control" being the reference level, then we may use coef="grouptreated", corresponding to the result of comparing treated to control group; or a list of two contrast vectors, each named by the colnames of the design matrix, the first corresponds to the baseline (e.g. control), the second corresponds to the group of interest (e.g. treated)
  # lfc.cutoff: a non-negative number, lower cutoff for log fold-change (will return genes whose log fold-change >= this value); if set >0, the P values will no longer be valid
  # pos.only: if TRUE, will return genes with log fold-change >= lfc.cutoff; otherwise, return genes with abs(log fold-change) >= lfc.cutoff
  # *note: I include lfc.cutoff and pos.only arguments here (rather than downstream) because I want to filter genes before doing LR test to save time (if p values are needed)
  # lfc.only: if TRUE, will return only log fold-change without p values
  # nc: number of cores to use for multi-core parallelization
  # return.fit: if TRUE, return list(fit=zlm.fit, summ=zlm.summ, de.res=de.res), where zlm.fit is the output from MAST::zlm, and zlm.summ if the output from summary(zlm.fit); otherwise return de.res
  # ...: passed to MAST::zlm

  if (!requireNamespace("MAST", quietly=TRUE)) {
    stop("Package \"MAST\" (GitHub fork: ImNotaGit/MAST) needed for this function to work.")
  }

  if (missing(design)) {
    if (missing(pheno)) stop("Need to provide either `pheno` with `model`, or `design`.")
    #vs <- unique(c(all.vars(model), names(model.frame(model, pheno))))
    #vs <- vs[vs!="." & !grepl("\\(|\\)", vs)]
    vs <- all.vars(model)
    if ("." %in% vs) vs <- unique(c(setdiff(vs, "."), names(pheno)))
    ccs <- complete.cases(pheno[, vs, with=FALSE])
    if (any(!ccs)) message("Removed ", sum(!ccs), " samples with incomplete (NA) covariate data.")
    dat <- dat[, ccs]
    cdat <- pheno[ccs]
    if (cdr) {
      cdat <- cbind(cdat, cdr_=Matrix::colMeans(dat>0))
      model <- update(model, ~.+cdr_)
    }
  } else {
    if (!missing(pheno)) message("`design` was provided, ignoring `pheno` and `model`.")
    colnames(design) <- make.names(colnames(design))
    if (cdr) design <- cbind(design, cdr_=Matrix::colMeans(dat>0))
    #model <- ~.+0 # doesn't work
    model <- as.formula(sprintf("~ %s + 0", paste(sprintf("`%s`", colnames(design)), collapse=" + ")))
    cdat <- design
  }
  
  sca <- MAST::FromMatrix(exprsArray=as.matrix(dat),
    cData=cbind(data.frame(wellKey=colnames(dat), row.names=colnames(dat)), cdat),
    fData=data.frame(primerid=rownames(dat), row.names=rownames(dat)),
    check_sanity=FALSE)
  
  nc.bak <- options("mc.cores")$mc.cores
  if (nc>1) options(mc.cores=nc)
  zlm.fit <- MAST::zlm(formula=model, sca=sca, parallel=isTRUE(nc>1), ...)
  
  if (is.character(coef)) {
    if (missing(design)) design <- zlm.fit@LMlike@modelMatrix
    coef <- make.names(coef)
    c0 <- setNames(ifelse(colnames(design)==make.names("(Intercept)"), 1, 0), colnames(design))
    c1 <- as.matrix(setNames(ifelse(colnames(design)==coef, 1, 0), colnames(design)))
    colnames(c1) <- coef
    # add intercept to c1
    c1 <- c1 + c0
    args.getlfc <- list(contrast0=c0, contrast1=c1)
    args.lrt <- list(hypothesis=MAST::CoefficientHypothesis(coef))
    coef1 <- coef
  } else if (is.list(coef)) {
    coef <- lapply(coef, function(x) setNames(x, make.names(names(x))))
    c0 <- as.matrix(coef[[1]])
    c1 <- as.matrix(coef[[2]])[rownames(c0), , drop=FALSE]
    colnames(c0) <- colnames(c1) <- "x"
    args.getlfc <- list(contrast0=coef[[1]], contrast1=c1)
    args.lrt <- list(hypothesis=c1-c0)
    coef1 <- "x"
  }
  
  if (lfc.cutoff>0 || pos.only) {
    # pre-filtering genes before doing statistical tests to save time
    if (lfc.cutoff>0 && !lfc.only) warning("Filtering genes based on log fold-change, P values no longer valid.")
    if (pos.only) gidx <- do.call(MAST::getLogFC, c(list(zlmfit=zlm.fit), args.getlfc))[contrast==coef1, logFC>=lfc.cutoff] else gidx <- do.call(MAST::getLogFC, c(list(zlmfit=zlm.fit), args.getlfc))[contrast==coef1, abs(logFC)>=lfc.cutoff]
    gidx[is.na(gidx)] <- FALSE
    if (sum(gidx)==0) {
      warning("No gene passes log fold-change cutoff, NULL returned.")
      return(NULL)
    }
    for (i in slotNames(zlm.fit)) {
      si <- slot(zlm.fit, i)
      if (is.array(si)) {
        if (is.matrix(si)) slot(zlm.fit, i) <- si[gidx, , drop=FALSE] else slot(zlm.fit, i) <- si[, , gidx, drop=FALSE]
      }
    }
    zlm.fit@sca <- zlm.fit@sca[gidx, ]
  }
  
  if (lfc.only) {
    zlm.summ <- MAST::summary(zlm.fit, logFC=do.call(MAST::getLogFC, c(list(zlmfit=zlm.fit), args.getlfc)), doLRT=FALSE)$datatable
    if (nc>1) options(mc.cores=nc.bak) # reset
    # discrete component (logistic)
    de.res <- zlm.summ[contrast==coef1 & component=='D', .(id=primerid, coef.d=coef, ci95lo.d=ci.lo, ci95up.d=ci.hi)]
    # continuous component
    tmp <- zlm.summ[contrast==coef1 & component=='C', .(id=primerid, coef.c=coef, ci95lo.c=ci.lo, ci95up.c=ci.hi)]
    de.res <- merge(de.res, tmp, by="id")
    # logFC
    tmp <- zlm.summ[contrast==coef1 & component=='logFC', .(id=primerid, log.fc=coef, ci95lo.lfc=ci.lo, ci95up.lfc=ci.hi)]
    de.res <- merge(de.res, tmp, by="id")
  } else {
    zlm.summ <- MAST::summary(zlm.fit, logFC=do.call(MAST::getLogFC, c(list(zlmfit=zlm.fit), args.getlfc)), doLRT=setNames(list(do.call(MAST::lrTest, c(list(object=zlm.fit), args.lrt))), coef1), parallel=isTRUE(nc>1))$datatable
    if (nc>1) options(mc.cores=NULL) # reset
    # discrete component (logistic)
    de.res <- zlm.summ[contrast==coef1 & component=='D', .(id=primerid, coef.d=coef, ci95lo.d=ci.lo, ci95up.d=ci.hi, pval.d=`Pr(>Chisq)`)][, padj.d:=p.adjust(pval.d, "BH")]
    # continuous component
    tmp <- zlm.summ[contrast==coef1 & component=='C', .(id=primerid, coef.c=coef, ci95lo.c=ci.lo, ci95up.c=ci.hi, pval.c=`Pr(>Chisq)`)][, padj.c:=p.adjust(pval.c, "BH")]
    de.res <- merge(de.res, tmp, by="id")
    # logFC
    tmp <- merge(
      zlm.summ[contrast==coef1 & component=='logFC', .(id=primerid, log.fc=coef, ci95lo.lfc=ci.lo, ci95up.lfc=ci.hi)],
      zlm.summ[contrast==coef1 & component=='H', .(id=primerid, pval=`Pr(>Chisq)`)][, padj:=p.adjust(pval, "BH")],
      by="id"
    )
    de.res <- merge(de.res, tmp, by="id")
  }
  #de.res <- de.res[order(-abs(log.fc))]
  de.res <- de.res[order(padj, pval)]
  if (return.fit) list(fit=zlm.fit, summ=zlm.summ, de.res=de.res) else de.res
}


make.pseudobulk <- function(mat, mdat, blk, ncells.cutoff=10) {
  # create pseudobulk gene expression data for single-cell RNA-seq
  # mat: gene-by-cell matrix (assuming sparse)
  # mdat: cell meta data, data.frame or other objects coerceable to data.table
  # blk: names of one or more sample/bulk variables (as in column names of mdat)
  # ncells.cutoff: only keep samples/bulks with > this number of cells
  # return: list(mat=mat, mdat=mdat), where mat is the pseudobulk gene-by-sample/bulk matrix, mdat is the corresponding meta data for the pseudobulk samples, the latter is obtained by dropping any variables (columns) with non-unique values within a bulk from the original cell-level mdat

  mdat <- as.data.table(mdat)
  tmp <- sapply(blk, function(i) anyNA(mdat[, i, with=FALSE]))
  if (any(tmp)) warning(sprintf("These bulk variables contain NA! NA will be kept as a separate level:\n%s\n", paste(blk[tmp], collapse=", ")))
  blk <- do.call(paste, c(unname(mdat[, blk, with=FALSE]), sep="_"))
  tmp <- table(blk)>ncells.cutoff
  if (any(!tmp)) {
    idx <- blk %in% names(tmp)[tmp]
    warning(sprintf("These bulks contain <=%d cell, they will be discarded:\n%s\n", ncells.cutoff, paste(names(tmp)[!tmp], collapse=", ")))
    blk <- blk[idx]
    mdat <- mdat[idx]
    mat <- mat[, idx, drop=FALSE]
  }
  blk <- factor(blk)
  tmp <- Matrix::sparseMatrix(i=1:length(blk), j=as.numeric(blk), dims=c(length(blk), length(levels(blk))), dimnames=list(NULL, levels(blk)))
  mat <- as.matrix(mat %*% tmp)
  tmp <- sapply(mdat, function(x) any(colSums(table(x, blk)>0)>1))
  if (any(tmp)) message(sprintf("These variables will be dropped since they have multiple values/levels within a bulk:\n%s\n", paste(names(mdat)[tmp], collapse=", "))
)
  mdat[, blk:=blk]
  mdat <- unique(mdat[, !tmp, with=FALSE])[match(colnames(mat), blk), -"blk"]
  list(mat=mat, mdat=mdat)
}


