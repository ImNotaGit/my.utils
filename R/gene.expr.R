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
  # quantile: whether to apply quantile normalization, if TRUE, will pass normalize.method="quantile" to limma::voom; this will be ignored if `normalize.method` is specified in ...
  # ...: passed to edgeR::calcNormFactors (e.g. `method`) and limma::voom
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
  dat <- pass3dots(edgeR::calcNormFactors.DGEList, dat, ...)

  if ("normalize.method" %in% names(list(...))) {
    v <- pass3dots(limma::voom, counts=dat, design=design, ...)
  } else {
    if (quantile) nm <- "quantile" else nm <- "none"
    v <- pass3dots(limma::voom, counts=dat, design=design, normalize.method=nm, ...)
  }

  list(voom=v, design=design, genes=gs)
}

de.limma <- function(dat, pheno=NULL, model=~., design=NULL, coef, contrast, reduced.model, contr.to.coef=FALSE, gene.colname=TRUE, keep.fit=FALSE, ...) {
  # differential expression analysis with limma
  # pheno: phenotypic data as a data.table with the same order of samples
  # model: the model to use for DE, by default a linear model containing all variables in pheno (w/o interaction terms)
  # design: design matrix for DE
  # if design is NULL, pheno and model will be used to compute the design matrix; otherwise design will be used, and pheno with model will be ignored if provided; need to provide either pheno with model, or design
  # coef, contrast and reduced.model are different ways to specify the model terms for which to perform DE test and to return DE results; provide at most one of these; if more than one is provided, will use reduced.model over contrast over coef; all can be provided as single items as described below or lists of items (named lists recommended) for multiple tests, in which case a corresponding list of DE result tables will be returned; if none of these three is provided, will return results for all coefficients in the model
  # coef: numeric or character vector of model coefficients (corresponding to columns of design matrix); if length>1, the coefficient (logFC) of each and a single P value for joint testing (?) will be returned
  # contrast: numeric contrast vector or matrix (for the latter, one contrast per column), or character vector specifying one or more contrasts in terms of the column names of the design matrix (in which case it will be converted to contrast vector/matrix with limma::makeContrasts); the matrix/multiple contrasts case is handled in the same way as the case of coef with length>1
  # reduced.model: formula of the reduced model (works only if pheno and model are provided), or vector of model coefficients (columns of design matrix) to keep (i.e., the opposite to coef, which specifies the coefficients to drop)
  # contr.to.coef: whether to reform the design matrix with limma::contrastAsCoef such that contrasts become coefficients
  # gene.colname: the column name for gene symbols in fData(dat) if dat is an ExpressionSet; if TRUE then will try to get gene symbols automatically from fData(dat); if FALSE then will not try to get gene symbols; gene symbols will be added to the returned DE table
  # keep.fit: if TRUE, then also return the fitted model in addition to the DE result table(s) as list(fit=fit, de.res=de.res), otherwise return de.res
  # ...: passed to limma::lmFit, limma::eBayes and limma::topTable, e.g. `robust` and `trend` for limma::eBayes, set these to TRUE for log-transformed RNA-seq data (but for RNA-seq doing DE with read count data using other methods can be more recommended)

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

  pars <- .process.de.params(dat=mat, pheno=pheno, model=model, design=design, coef=coef, contrast=contrast, reduced.model=reduced.model, contr.to.coef=contr.to.coef)

  if (!contr.to.coef) {
    fit <- pass3dots(limma::lmFit, pars$dat, design=pars$design, ...)
    if (!is.null(pars$contrast)) {
      pars$coef <- lapply(pars$contrast, function(x) if (is.matrix(x)) 1:ncol(x) else 1)
      tmp <- mapply(function(contrast, coef) {
        fit <- limma::contrasts.fit(fit, contrast)
        fit <- pass3dots(limma::eBayes, fit, ...)
        de.res <- tryCatch(pass3dots(limma::topTable, fit, coef=coef, number=Inf, genelist=gns, ...), error=function(e) pass3dots(limma::topTable, fit, coef=coef, number=Inf, ...))
        list(fit=fit, de.res=de.res)
      }, pars$contrast, pars$coef, SIMPLIFY=FALSE)
      fit <- lapply(tmp, function(x) x$fit)
      de.res <- lapply(tmp, function(x) x$de.res)
    } else if (!is.null(pars$coef)) {
      fit <- pass3dots(limma::eBayes, fit, ...)
      de.res <- lapply(pars$coef, function(x) {
        tryCatch(pass3dots(limma::topTable, fit, coef=x, number=Inf, genelist=gns, ...), error=function(e) pass3dots(limma::topTable, fit, coef=x, number=Inf, ...))
      })
    }
  } else {
    tmp <- mapply(function(design, coef) {
      fit <- pass3dots(limma::lmFit, pars$dat, design=design, ...)
      fit <- pass3dots(limma::eBayes, fit, ...)
      de.res <- tryCatch(pass3dots(limma::topTable, fit, coef=coef, number=Inf, genelist=gns, ...), error=function(e) pass3dots(limma::topTable, fit, coef=coef, number=Inf, ...))
      list(fit=fit, de.res=de.res)
    }, pars$design, pars$coef, SIMPLIFY=FALSE)
    fit <- lapply(tmp, function(x) x$fit)
    de.res <- lapply(tmp, function(x) x$de.res)
  }

  de.res <- lapply(de.res, function(x) {
    res <- as.data.table(x)
    if (!"id" %in% tolower(names(res))) res <- cbind(ID=row.names(res), res)
    setnames(res, c("ID","logFC","AveExpr","P.Value","adj.P.Val"), c("id","log.fc","ave.expr","pval","padj"), skip_absent=TRUE)
    res[order(padj, pval)]
  })

  if (length(de.res)==1) de.res <- de.res[[1]]
  if (keep.fit) list(fit=fit, de.res=de.res) else de.res
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


.process.de.params <- function(dat, pheno, model=~., design, coef, contrast, reduced.model, contr.to.coef=FALSE, make.coef.names=FALSE) {

  if (missing(design) || is.null(design)) {
    if (missing(pheno) || is.null(pheno)) stop("Need to provide either `pheno` with `model`, or `design`.")
    pheno <- as.data.table(pheno)
    vs <- unique(c(all.vars(model), names(model.frame(model, pheno))))
    vs <- vs[vs!="." & !grepl("\\(|\\)", vs)]
    ccs <- complete.cases(pheno[, vs, with=FALSE])
    if (any(!ccs)) {
      dat <- dat[, ccs]
      pheno <- pheno[ccs]
      message("Removed ", sum(!ccs), " samples with incomplete (NA) covariate data.")
    } else pheno <- copy(pheno)
    tmp <- sapply(pheno[, vs, with=FALSE], function(x) !is.numeric(x) & !is.factor(x))
    if (any(tmp)) {
      message("These non-numeric variables included in the model are not factors:")
      message(cc(vs[tmp]))
      message("They are converted to factors.")
      pheno[, c(vs[tmp]):=lapply(.SD, factor), .SDcols=vs[tmp]]
    }
    design <- model.matrix(model, pheno)
    red.mod.fml.ok <- TRUE
  } else {
    if (!(missing(pheno) || is.null(pheno))) {
      warning("Both `pheno` with `model` and `design` are provided, will use `design` and ignore `pheno` with `model`")
    } else pheno <- NULL
    ccs <- rep(TRUE, nrow(design))
    red.mod.fml.ok <- FALSE
  }

  if (!(missing(reduced.model) || is.null(reduced.model))) {
    if (!(missing(contrast) || is.null(contrast)) || !(missing(coef) || is.null(coef))) {
      warning("`reduced.model` is provided, will ignore `contrast` and `coef`.")
      contrast <- coef <- NULL
    }
    if (!is.list(reduced.model)) reduced.model <- list(reduced.model)
    coef <- lapply(reduced.model, function(x) {
      if (class(x)=="formula") {
        if (!red.mod.fml.ok) stop("Not using `pheno` with `model`, `reduced.model` can only be a single or a list of numeric/character vectors.")
        setdiff(colnames(design), colnames(model.matrix(x, pheno)))
      } else if (is.character(x)) {
        setdiff(colnames(design), x)
      } else if (is.numeric(x)) {
        colnames(design)[-x]
      } else stop("Invalid `reduced.model`, should be a single or a list of formulae or numeric/character vectors.")
    })
  } else reduced.model <- NULL

  if (!(missing(contrast) || is.null(contrast))) {
    if (!(missing(coef) || is.null(coef))) {
      warning("`contrast` is provided, will ignore `coef`.")
      coef <- NULL
    }
    if (!is.list(contrast)) contrast <- list(contrast)
    contrast <- lapply(contrast, function(x) {
      if (is.character(x)) makeContrasts(contrasts=x, levels=design, check.names=FALSE) else x # using my copy of makeContrasts
    })
    if (contr.to.coef) {
      message("Reforming design matrix with limma::contrastAsCoef such that contrasts become coefficients.")
      tmp <- lapply(contrast, function(x) {
        tmp <- limma::contrastAsCoef(design, x)
        list(design=tmp$design, coef=tmp$coef)
      })
      design <- lapply(tmp, function(x) x$design)
      coef <- lapply(tmp, function(x) x$coef)
    }
  } else {
    contrast <- NULL
    if (missing(coef) || is.null(coef)) {
      message("No `coef`, `contrast`, or `reduced.model` was provided, will return all coefficients in the model.")
      coef <- colnames(design)
      coef <- setNames(as.list(coef), coef)
    }
  }

  if (!is.null(coef)) {
    names(coef)[names(coef)=="(Intercept)"] <- "Intercept"
    if (make.coef.names) {
      coef <- lapply(coef, function(x) {
        x[x=="(Intercept)"] <- "Intercept"
        make.names(x)
      })
    }
  }

  list(dat=dat, pheno=pheno, model=model, design=design, coef=coef, contrast=contrast, ccs=ccs)
}


.calc.norm.factors.with.ctrl <- function(object, ctrl, lib.size = NULL, method=c("TMM", "RLE"), refColumn = NULL, doWeighting = TRUE, Acutoff = -1e+10) {
  # my copy of edgeR::calcNormFactors.default allowing for specifying of control features/genes, i.e. features/genes that are exprected to remain stable across conditions.
  # ctrl: one or more control features/genes as row indices or row names of object, which is a gene (feature)-by-sample matrix

  if (!requireNamespace("edgeR", quietly=TRUE)) {
    stop("Package \"edgeR\" needed for this function to work.")
  }

  calc.factor.tmm <- function(obs, ctrl, ref, libsize.obs, libsize.ref, doWeighting, Acutoff) {
    # ctrl: index of control feature(s)
    obs <- as.numeric(obs)
    ref <- as.numeric(ref)
    if (is.null(libsize.obs)) nO <- sum(obs) else nO <- libsize.obs
    if (is.null(libsize.ref)) nR <- sum(ref) else nR <- libsize.ref
    logR <- log2((obs/nO)/(ref/nR))
    absE <- (log2(obs/nO) + log2(ref/nR))/2
    v <- (nO - obs)/nO/obs + (nR - ref)/nR/ref
    logR <- logR[ctrl]
    absE <- absE[ctrl]
    v <- v[ctrl]
    fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
    if (sum(fin)==0) stop("No control feature has finite logR with finite absE>Acutoff.")
    logR <- logR[fin]
    absE <- absE[fin]
    v <- v[fin]
    if (max(abs(logR)) < 1e-06) return(1)
    if (doWeighting) f <- sum(logR/v, na.rm = TRUE)/sum(1/v, na.rm = TRUE) else f <- mean(logR, na.rm = TRUE)
    if (is.na(f)) f <- 0
    2^f
  }

  x <- as.matrix(object)
  if (any(is.na(x))) stop("NA counts not permitted")
  if (is.character(ctrl)) ctrl <- match(ctrl, rownames(x))
  if (all(is.na(ctrl))) stop("None of the provided control features is found in data.")
  if (anyNA(ctrl)) message("Some control features not in data.")
  ctrl <- ctrl[!is.na(ctrl)]
  if (any(colSums(x[ctrl, , drop=FALSE])==0)) stop("There exist some samples where all control features have zero counts. Normalization with control features cannot proceed.")
  nsamples <- ncol(x)
  if (is.null(lib.size)) {
    lib.size <- colSums(x)
  } else {
    if (anyNA(lib.size)) stop("NA lib.sizes not permitted")
    if (length(lib.size) != nsamples) {
      if (length(lib.size) > 1L) warning(".calc.norm.factors.tmm.with.ctrl: length(lib.size) doesn't match number of samples", call. = FALSE)
      lib.size <- rep_len(lib.size, nsamples)
    }
  }
  method <- match.arg(method)
  allzero <- .rowSums(x > 0, nrow(x), nsamples) == 0L
  if (any(allzero)) x <- x[!allzero, , drop = FALSE]
  f <- switch(method, TMM = {
    if (is.null(refColumn)) {
      f75 <- suppressWarnings(edgeR:::.calcFactorQuantile(data = x, lib.size = lib.size, p = 0.75))
      if (median(f75) < 1e-20) refColumn <- which.max(colSums(sqrt(x))) else refColumn <- which.min(abs(f75 - mean(f75)))
    }
    f <- rep_len(NA_real_, nsamples)
    for (i in 1:nsamples) f[i] <- calc.factor.tmm(obs = x[, i], ctrl=ctrl, ref = x[, refColumn], libsize.obs = lib.size[i], libsize.ref = lib.size[refColumn], doWeighting = doWeighting, Acutoff = Acutoff)
    f
  }, RLE = {
    gm <- exp(rowMeans(log(x[ctrl, , drop=FALSE])))
    f <- apply(x[ctrl, , drop=FALSE], 2, function(u) median((u/gm)[gm > 0]))
    if (anyNA(f)) stop("There are NA's in normalization factors.")
    f/lib.size
  })
  f <- f/exp(mean(log(f)))
  names(f) <- colnames(x)
  f
}


de.edger <- function(dat, pheno, model=~., design, coef, contrast, reduced.model, contr.to.coef=FALSE, lfc.cutoff=0, keep.fit=FALSE, ctrl.features, norm.factors, ...) {
  # differential expression analysis with edgeR
  # dat: gene-by-sample expression matrix of raw counts; should have low genes already filtered out
  # pheno: phenotypic data as a data.table with the same order of samples
  # model: the model to use for DE, by default a linear model containing all variables in pheno (w/o interaction terms)
  # design: design matrix for DE
  # if design is NULL, pheno and model will be used to compute the design matrix; otherwise design will be used, and pheno with model will be ignored if provided; need to provide either pheno with model, or design
  # coef, contrast and reduced.model are different ways to specify the model terms for which to perform DE test and to return DE results; provide at most one of these; if more than one is provided, will use reduced.model over contrast over coef; all can be provided as single items as described below or lists of items (named lists recommended) for multiple tests, in which case a corresponding list of DE result tables will be returned; if none of these three is provided, will return results for all coefficients in the model
  # coef: numeric or character vector of model coefficients (corresponding to columns of design matrix); if length>1, the coefficient (logFC) of each and a single P value for joint testing (?) will be returned
  # contrast: numeric contrast vector or matrix (for the latter, one contrast per column), or character vector specifying one or more contrasts in terms of the column names of the design matrix (in which case it will be converted to contrast vector/matrix with limma::makeContrasts); the matrix/multiple contrasts case is handled in the same way as the case of coef with length>1
  # reduced.model: formula of the reduced model (works only if pheno and model are provided), or vector of model coefficients (columns of design matrix) to keep (i.e., the opposite to coef, which specifies the coefficients to drop)
  # contr.to.coef: whether to reform the design matrix with limma::contrastAsCoef such that contrasts become coefficients
  # keep.fit: if TRUE, then also return the fitted model in addition to the DE result table(s) as list(fit=fit, de.res=de.res), otherwise return de.res
  # ctrl.features: control features/genes that are expected to remain stable across conditions; if provided (i.e. not missing or NULL), will use my custom function .calc.norm.factors.with.ctrl instead of edge::calcNormFactors, and assign the result to dge@samples$norm.factors; for this, the normalization method (specified via `method` argument in ...) can only be "TMM" (default) or "RLE"
  # norm.factors: custom edgeR normalization factors, a numeric vector, will directly set dge$samples$norm.factors to this if provided; otherwise (i.e. if missing or NULL), will use edgeR::calcNormFactors (or my custom function, if ctrl.features is provided)
  # note: in edgeR, norm.factors is added on to library size (total counts per sample), i.e. norm.factors*library.size is used for normalization
  # note: there are two ways to use only library.size for normalization, one is to specify norm.factors=1, the other is to pass method="none" within ... (see below)
  # ...: any possible additional arguments of calcNormFactors, estimateDisp, glmQLFit, glmQLFTest, and glmTreat, e.g. `method` for calcNormFactors, `trend.method` for estimateDisp, `robust` for estimateDisp and glmQLFit (somewhat different meanings in both but may as well be set identically), `abundance.trend` for glmQLFit, etc.

  if (!requireNamespace("edgeR", quietly=TRUE)) {
    stop("Package \"edgeR\" needed for this function to work.")
  }

  pars <- .process.de.params(dat=dat, pheno=pheno, model=model, design=design, coef=coef, contrast=contrast, reduced.model=reduced.model, contr.to.coef=contr.to.coef)
  dge <- edgeR::DGEList(counts=pars$dat)
  if (missing(norm.factors) || is.null(norm.factors)) {
    if (missing(ctrl.features) || is.null(ctrl.features)) {
      dge <- pass3dots(edgeR::calcNormFactors.DGEList, dge, ...)
    } else {
      m <- list(...)$method
      if (!is.null(m) && !m %in% c("TMM","RLE")) stop("Only \"TMM\" and \"RLE\" are supported normalization methods when `ctrl.features` is provided.")
      norm.factors <- pass3dots(.calc.norm.factors.with.ctrl, dge$counts, ctrl=ctrl.features, lib.size=dge$samples$lib.size, ...)
      dge$samples$norm.factors <- norm.factors
    }
  } else {
    if (length(norm.factors)==1) norm.factors <- rep(norm.factors, ncol(pars$dat)) else norm.factors <- norm.factors[pars$ccs]
    if (is.null(names(norm.factors))) names(norm.factors) <- colnames(pars$dat)
    dge$samples$norm.factors <- norm.factors
  }

  if (!contr.to.coef) {
    dge <- pass3dots(edgeR::estimateDisp.DGEList, dge, pars$design, ...)
    fit <- pass3dots(edgeR::glmQLFit.DGEList, dge, pars$design, ...)
    if (lfc.cutoff==0) {
      if (!is.null(pars$contrast)) {
        de.res <- lapply(pars$contrast, function(x) pass3dots(edgeR::glmQLFTest, fit, contrast=x, ...))
      } else if (!is.null(pars$coef)) {
        de.res <- lapply(pars$coef, function(x) pass3dots(edgeR::glmQLFTest, fit, coef=x, ...))
      }
    } else {
      if (!is.null(pars$contrast)) {
        de.res <- lapply(pars$contrast, function(x) pass3dots(edgeR::glmTreat, fit, contrast=x, lfc=lfc.cutoff, ...))
      } else if (!is.null(pars$coef)) {
        de.res <- lapply(pars$coef, function(x) pass3dots(edgeR::glmTreat, fit, coef=x, lfc=lfc.cutoff, ...))
      }
    }
  } else {
    tmp <- mapply(function(design, coef) {
      dge <- pass3dots(edgeR::estimateDisp.DGEList, dge, design, ...)
      fit <- pass3dots(edgeR::glmQLFit.DGEList, dge, design, ...)
      if (lfc.cutoff==0) de.res <- pass3dots(edgeR::glmQLFTest, fit, coef=coef, ...) else de.res <- pass3dots(edgeR::glmTreat, fit, coef=coef, lfc=lfc.cutoff, ...)
      list(fit=fit, de.res=de.res)
    }, pars$design, pars$coef, SIMPLIFY=FALSE)
    fit <- lapply(tmp, function(x) x$fit)
    de.res <- lapply(tmp, function(x) x$de.res)
  }

  de.res <- lapply(de.res, function(x) {
    tmp <- edgeR::topTags(x, n=Inf)[[1]]
    res <- cbind(id=row.names(tmp), as.data.table(tmp))
    setnames(res, c("logFC","logCPM","PValue","FDR"), c("log.fc","log.cpm","pval","padj"), skip_absent=TRUE)
    setnames(res, stringr::str_replace_all(names(res), "logFC", "log.fc"))
    res[order(padj, pval)]
  })

  if (length(de.res)==1) de.res <- de.res[[1]]
  if (keep.fit) list(fit=fit, de.res=de.res) else de.res
}


de.deseq2 <- function(dat, pheno, model=~., design, coef, contrast, reduced.model, contr.to.coef=FALSE, keep.fit=FALSE, nc=1L, ctrl.features, size.factors, ...) {
  # differential expression analysis with DESeq2
  # dat: gene-by-sample expression matrix of raw counts; should have low genes already filtered out
  # pheno: phenotypic data as a data.table with the same order of samples
  # model: the model to use for DE, by default a linear model containing all variables in pheno (w/o interaction terms)
  # design: design matrix for DE
  # if design is NULL, pheno and model will be used to compute the design matrix; otherwise design will be used, and pheno with model will be ignored if provided; need to provide either pheno with model, or design
  # coef, contrast and reduced.model are different ways to specify the model terms for which to perform DE test and to return DE results; provide at most one of these; if more than one is provided, will use reduced.model over contrast over coef; all can be provided as single items as described below or lists of items (named lists recommended) for multiple tests, in which case a corresponding list of DE result tables will be returned; if none of these three is provided, will return results for all coefficients in the model
  # coef: numeric or character vector of model coefficients (corresponding to columns of design matrix); if length>1, the coefficient (logFC) of each and a single P value for joint testing (?) will be returned
  # contrast: numeric contrast vector or matrix (for the latter, one contrast per column), or character vector specifying one or more contrasts in terms of the column names of the design matrix (in which case it will be converted to contrast vector/matrix with limma::makeContrasts); the matrix/multiple contrasts case is handled in the same way as the case of coef with length>1
  # note: coef is passed to the `name` argument of DESeq2::results, while contrast is passed to the `contrast` argument of DESeq2::results; the character vector mode for contrast of DESeq2::results is not supported here (where, e.g. contrast=c("group", "trt", "ctrl") will return results for the 'trt' level compared to 'ctrl' level of the `group` variable)
  # reduced.model: formula of the reduced model (works only if pheno and model are provided), or vector of model coefficients (columns of design matrix) to keep (i.e., the opposite to coef, which specifies the coefficients to drop)
  # contr.to.coef: whether to reform the design matrix with limma::contrastAsCoef such that contrasts become coefficients
  # keep.fit: if TRUE, then also return the fitted model in addition to the DE result table(s) as list(fit=fit, de.res=de.res), otherwise return de.res
  # nc: number of cores for parallelization
  # ctrl.features: control features/genes that are expected to remain stable across conditions; if provided (i.e. not missing or NULL), will be passed to the `controlGenes` argument of DESeq2::estimateSizeFactors, but unlike `controlGenes`, this can be provided as a character vector of feature/gene symbols (rownames of dat)
  # size.factors: custom size factors for DESeq2, a numeric vector, will directly set sizeFactors(dds) to this if provided; otherwise (i.e. if missing or NULL), will use DESeq2::estimateSizeFactors
  # note: the DESeq2 size factor is directly used for normalization (unlike the edgeR normalization factor, which is multiplied by library size before being used for normalization)
  # ...: passed to DESeq2::DESeq and DESeq2::results

  if (!requireNamespace(c("DESeq2","BiocParallel"), quietly=TRUE)) {
    stop("Packages \"DESeq2\" and \"BiocParallel\" needed for this function to work.")
  }

  if (any(!is.wholenumber(dat))) {
    message("DESeq2 only accepts integer read counts. There are non-integer values in `dat` and they have been rounded.")
    dat <- round(dat)
  }

  pars <- .process.de.params(dat=dat, pheno=pheno, model=model, design=design, coef=coef, contrast=contrast, reduced.model=reduced.model, contr.to.coef=contr.to.coef, make.coef.names=TRUE)
  if (is.null(pars$pheno)) pars$pheno <- data.table(idx=1:ncol(pars$dat))
  if (nc>1) bp <- BiocParallel::MulticoreParam(workers=nc, progressbar=TRUE, RNGseed=0) else bp <- BiocParallel::bpparam()

  if (!(missing(ctrl.features) || is.null(ctrl.features))) {
    if (is.character(ctrl.features)) ctrl.features <- match(ctrl.features, rownames(pars$dat))
    if (all(is.na(ctrl.features))) stop("None of the provided `ctrl.features` is found in data.")
    if (anyNA(ctrl.features)) message("Some `ctrl.features` not in data.")
    ctrl.features <- ctrl.features[!is.na(ctrl.features)]
  }

  if (!(missing(size.factors) || is.null(size.factors))) {
    if (length(size.factors)==1) size.factors <- rep(size.factors, ncol(pars$dat)) else size.factors <- size.factors[pars$ccs]
  }

  if (!contr.to.coef) {
    dds <- DESeq2::DESeqDataSetFromMatrix(countData=pars$dat, colData=pars$pheno, design=pars$design)
    if (missing(size.factors) || is.null(size.factors)) {
      dds <- pass3dots(DESeq2:::estimateSizeFactors.DESeqDataSet, dds, controlGenes=ctrl.features, ...)
    } else sizeFactors(dds) <- size.factors
    fit <- pass3dots(DESeq2::DESeq, dds, parallel=nc>1, BPPARAM=bp, ...)
    if (!is.null(pars$contrast)) {
      de.res <- lapply(pars$contrast, function(x) pass3dots(DESeq2::results, fit, contrast=x, parallel=nc>1, BPPARAM=bp, ...))
    } else if (!is.null(pars$coef)) {
      de.res <- lapply(pars$coef, function(x) pass3dots(DESeq2::results, fit, name=x, parallel=nc>1, BPPARAM=bp, ...))
    }
  } else {
    tmp <- mapply(function(design, coef) {
      dds <- DESeq2::DESeqDataSetFromMatrix(countData=pars$dat, colData=pars$pheno, design=design)
      if (missing(size.factors) || is.null(size.factors)) {
        dds <- pass3dots(DESeq2:::estimateSizeFactors.DESeqDataSet, dds, controlGenes=ctrl.features, ...)
      } else sizeFactors(dds) <- size.factors
      fit <- pass3dots(DESeq2::DESeq, dds, parallel=nc>1, BPPARAM=bp, ...)
      de.res <- pass3dots(DESeq2::results, fit, name=coef, parallel=nc>1, BPPARAM=bp, ...)
      list(fit=fit, de.res=de.res)
    }, pars$design, pars$coef, SIMPLIFY=FALSE)
    fit <- lapply(tmp, function(x) x$fit)
    de.res <- lapply(tmp, function(x) x$de.res)
  }

  de.res <- lapply(de.res, function(x) {
    res <- cbind(id=row.names(x), as.data.table(x))
    setnames(res, c("baseMean","log2FoldChange","lfcSE","pvalue"), c("ave.expr","log.fc","lfc.se","pval"), skip_absent=TRUE)
    res[order(padj, pval)]
  })

  if (length(de.res)==1) de.res <- de.res[[1]]
  if (keep.fit) list(fit=fit, de.res=de.res) else de.res
}


de.glmgampoi <- function(dat, pheno, model=~., design, coef, contrast, reduced.model, contr.to.coef=FALSE, size.factors, pseudobulk=NULL, keep.fit=FALSE, ...) {
  # differential expression analysis with glmGamPoi, this may be used for single-cell RNA-seq data
  # dat: gene-by-sample expression matrix of log-normalized expression value, with gene ID/symbol as rownames and sample ID/barcode as colnames
  # pheno: phenotypic data as a data.table with the same order of samples
  # model: the model to use for DE, by default a linear model containing all variables in pheno (w/o interaction terms)
  # design: design matrix for DE
  # if design is NULL, pheno and model will be used to compute the design matrix; otherwise design will be used, and pheno with model will be ignored if provided; need to provide either pheno with model, or design
  # coef, contrast and reduced.model are different ways to specify the model terms for which to perform DE test and to return DE results; provide at most one of these; if more than one is provided, will use reduced.model over contrast over coef; all can be provided as single items as described below or lists of items (named lists recommended) for multiple tests, in which case a corresponding list of DE result tables will be returned; if none of these three is provided, will return results for all coefficients in the model
  # coef: numeric or character vector of model coefficients (corresponding to columns of design matrix); if length>1, the coefficient (logFC) of each and a single P value for joint testing (?) will be returned
  # contrast: numeric contrast vector or matrix (for the latter, one contrast per column), or character vector specifying one or more contrasts in terms of the column names of the design matrix (in which case it will be converted to contrast vector/matrix with limma::makeContrasts); the matrix/multiple contrasts case is handled in the same way as the case of coef with length>1
  # note: coef is passed to the `name` argument of DESeq2::results, while contrast is passed to the `contrast` argument of DESeq2::results; the character vector mode for contrast of DESeq2::results is not supported here (where, e.g. contrast=c("group", "trt", "ctrl") will return results for the 'trt' level compared to 'ctrl' level of the `group` variable)
  # reduced.model: formula of the reduced model (works only if pheno and model are provided), or vector of model coefficients (columns of design matrix) to keep (i.e., the opposite to coef, which specifies the coefficients to drop)
  # contr.to.coef: whether to reform the design matrix with limma::contrastAsCoef such that contrasts become coefficients
  # size.factors: passed to glmGamPoi::glm_gp `size_factors`, if missing use the default "normed_sum"
  # pseudobulk: passed to glmGamPoi::test_de `pseudobulk_by`
  # keep.fit: if TRUE, then also return the fitted model in addition to the DE result table(s) as list(fit=fit, de.res=de.res), otherwise return de.res
  # nc: number of cores for parallelization
  # ...: passed to glmGamPoi::glm_gp and glmGamPoi::test_de

  if (!requireNamespace("glmGamPoi", quietly=TRUE)) {
    stop("Package \"glmGamPoi\" needed for this function to work.")
  }

  if (missing(size.factors)) size.factors <- "normed_sum"

  pars <- .process.de.params(dat=dat, pheno=pheno, model=model, design=design, coef=coef, contrast=contrast, reduced.model=reduced.model, contr.to.coef=contr.to.coef)
  # glmGamPoi::test_de has a bug if requesting contrast="(Intercept)", so for now I will simply exclude the intercept if need to return all coefficients
  pars$coef <- pars$coef[names(pars$coef)!="Intercept"]

  if (!contr.to.coef) {
    fit <- pass3dots(glmGamPoi::glm_gp, as.matrix(pars$dat), design=pars$design, size_factors=size.factors, ...)
    if (!is.null(pars$contrast)) {
      de.res <- lapply(pars$contrast, function(x) pass3dots(glmGamPoi::test_de, fit, contrast=x, pseudobulk_by=pseudobulk, ...))
    } else if (!is.null(pars$coef)) {
      de.res <- lapply(pars$coef, function(x) pass3dots(glmGamPoi::test_de, fit, contrast=x, pseudobulk_by=pseudobulk, ...))
    }
  } else {
    tmp <- mapply(function(design, coef) {
      fit <- pass3dots(glmGamPoi::glm_gp, as.matrix(pars$dat), design=design, size_factors=size.factors, ...)
      de.res <- pass3dots(glmGamPoi::test_de, fit, contrast=coef, pseudobulk_by=pseudobulk, ...)
      list(fit=fit, de.res=de.res)
    }, pars$design, pars$coef, SIMPLIFY=FALSE)
    fit <- lapply(tmp, function(x) x$fit)
    de.res <- lapply(tmp, function(x) x$de.res)
  }

  de.res <- lapply(de.res, function(x) {
    res <- as.data.table(x)
    setnames(res, c("name","lfc","f_statistic","adj_pval"), c("id","log.fc","F","padj"), skip_absent=TRUE)
    res[order(padj, pval)]
  })

  if (length(de.res)==1) de.res <- de.res[[1]]
  if (keep.fit) list(fit=fit, de.res=de.res) else de.res
}


de.mast <- function(dat, pheno, model=~., design, cdr=TRUE, coef, lfc.cutoff=0, pos.only=FALSE, lfc.only=FALSE, nc=1L, keep.fit=FALSE, ...) {
  # differential expression analysis with MAST, used for e.g. single-cell RNA-seq data
  # dat: gene-by-sample expression matrix (assuming sparse) of log-normalized expression value, with gene ID/symbol as rownames and sample ID/barcode as colnames
  # pheno: phenotypic data as a data.table with the same order of samples
  # model: the model to use for DE, by default a linear model containing all variables in pheno (w/o interaction terms)
  # pheno and model will be used to compute the design matrix; or provide design matrix in design; pheno and design cannot both be missing; the design matrix should have proper column names
  # cdr: whether to include cellular detection rate (i.e. fraction of >0 genes in each cell) as a covariate
  # coef: character vector of model coefficients to test for and to return, e.g. each can be the name of the variable (and its level, if categorical) of interest for which the linear model coefficients to be displayed, e.g. if there's a variable named "group" with two levels "control" and "treated" with "control" being the reference level, then we may use coef="grouptreated", corresponding to the result of comparing treated to control group;
  # or a list of two contrast vectors, each named by the colnames of the design matrix, the first corresponds to the baseline (e.g. control), the second corresponds to the group of interest (e.g. treated)
  # if length of coef is >1, the returned de.res will be a list of DE result tables named by coef; otherwise de.res will be a single table
  # lfc.cutoff: a non-negative number, lower cutoff for log fold-change (will return genes whose log fold-change >= this value); if set >0, the P values will no longer be valid
  # pos.only: if TRUE, will return genes with log fold-change >= lfc.cutoff; otherwise, return genes with abs(log fold-change) >= lfc.cutoff
  # *note: I include lfc.cutoff and pos.only arguments here (rather than downstream) because I want to filter genes before doing LR test to save time (if p values are needed)
  # lfc.only: if TRUE, will return only log fold-change without p values
  # nc: number of cores to use for multi-core parallelization
  # keep.fit: if TRUE, return list(fit=zlm.fit, summ=zlm.summ, de.res=de.res), where zlm.fit is the output from MAST::zlm, and zlm.summ if the output from summary(zlm.fit); otherwise return de.res
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
  if (keep.fit) list(fit=zlm.fit, summ=zlm.summ, de.res=de.res) else de.res
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
  blk <- do.call(paste, c(unname(as.list(mdat[, blk, with=FALSE])), list(sep="_")))
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


