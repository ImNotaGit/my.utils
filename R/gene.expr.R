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


.process.de.params <- function(dat, pheno, model=~., design, coef, contrast, reduced.model, contr.to.coef=FALSE, make.contr=TRUE, make.coef.names=FALSE) {
  # a common helper function to various de.* functions, preprocessing their input to prepare the design matrix and intended DE tests
  # dat: gene-by-sample expression matrix of raw counts; should have low genes already filtered out
  # pheno: phenotypic data as a data.table with the same order of samples
  # model: the model to use for DE, by default a linear model containing all variables in pheno (w/o interaction terms)
  # design: design matrix for DE
  # if design is NULL, pheno and model will be used to compute the design matrix; otherwise design will be used, model will be ignored, pheno will not be used to form the design matrix used for DE if provided, but will be used to figure out terms in the reduced model if reduced.model is provided as formula; need to provide either pheno with model, or design
  # coef, contrast and reduced.model are different ways to specify the model coefficients or comparisons for which to test DE and return results, any combination of these can be provided, with each can be provided as single items or lists of items (named lists recommended) for multiple tests; these specifications will be combined in the order of coef, reduced.model, and contrast (with the orders within each kept the same as the input if multiple items were provided), and list of DE result tables in the corresponding order will be returned; if none of these three is provided, will return results for all coefficients in the model; see below for details;
  # coef: numeric or character vector of model coefficients (corresponding to columns of design matrix); if length>1, the coefficient (logFC) of each and a single P value for joint testing (i.e. an ANOVA/drop()-like test, depending on specific methods) will be returned; for multiple tests, provide a list of such items;
  # contrast: numeric contrast vector or matrix (for the latter, one contrast per column), or character vector specifying one or more contrasts in terms of the column names of the design matrix (e.g. sth like "grpA-grpB", in which case it will be converted to contrast vector/matrix with limma::makeContrasts); the matrix case and character vector of length>1 case is handled in the same way as the case of coef with length>1 (i.e. "joint" testing), which may be useful in some cases, but this requires the multiple contrasts are not collinear and probably require contr.to.coef=TRUE for many methods (not tested); for the common use case of separately checking multiple contrasts, pass them in a list of numeric vectors (instead of a single contrast matrix) or a list of atom character vectors (instead of a single length>1 character vector)
  # reduced.model: formula of the reduced model (works only if pheno is provided), or vector of model coefficients (columns of design matrix) to keep (i.e., the opposite to coef, which specifies the coefficients to drop), or design matrix of the reduced model; for multiple tests, provide a list of such items;
  # contr.to.coef: whether to reform the design matrix with limma::contrastAsCoef such that contrasts become coefficients
  # make.contr: if FALSE, the character contrast input will be left as is without being converted to contrast vectors with limma::makeContrasts; this is for de.glmgampoi, where (for now) `contrast` needs to be character, and needs to set contr.to.coef=FALSE
  # make.coef.names: whether to apply make.names() to model coefficients -- set to TRUE for de.deseq2, which somehow internally changes the names of the model variables (not sure whether also with make.names(), but converting coefs with make.names() works at least for the currently tested cases)
  # will return a list(dat, pheno, model, design, design.contr, coef, coef.contr, contrast, reduced.model, reduced.model.contr, ccs, named); named: if any of the input coef/contrast/reduced.model was named;
  # the returned $dat and $pheno will contain only complete.cases wrt model variables, and $ccs is the index for the complete cases wrt the original input dat/pheno;
  # the returned $coef will correspond to input `coef` + coefficients converted from input `reduced.model`;
  # the returned $contrast will correspond to the input `contrast`; additionally if contr.to.coef is TRUE, input `contrast` will be converted to coefficients and saved in the returned $coef.contr, and in the mean time $design.contr will be the transformed design matrices corresponding to each contr.to.coef conversion
  # the returned $reduced.model are the design matrices of reduced models corresponding to all input `coef` and input `reduced.model`;
  # if contr.to.coef is TRUE, the returned $reduced.model.contr is the design matrices of reduced models corresponding to all coefficients converted from input `contrast`;
  # the $reduced.model and $reduced.model.contr output are for de.glmgampoi to be passed to the `reduced_design` argument of glmGamPoi::test_de, or for de.deseq2 with LRT test, to be passed to the `reduced` argument of DESeq2::DESeq (DESeq2 with LRT has a different workflow and requires different arguments from Wald test)
  # missing items in the returned list will be NULL

  if (!(missing(pheno) || is.null(pheno))) {
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
  }

  if (missing(design) || is.null(design)) {
    if (missing(pheno) || is.null(pheno)) stop("Need to provide either `pheno` with `model`, or `design`.")
    design <- model.matrix(model, pheno)
  } else {
    if (!(missing(pheno) || is.null(pheno))) {
      message("Both `pheno` with `model` and `design` are provided, will use `design`, instead of creating design matrix from `pheno` with `model`")
    }
    ccs <- rep(TRUE, nrow(design))
  }

  named <- FALSE

  # handling `coef`
  if (!(missing(coef) || is.null(coef))) {
    if (!is.list(coef)) coef <- list(coef)
    for (tmp in coef) {
      if (any(!tmp %in% colnames(design))) stop("Invalid `coef`: some are not present in the model or `design`.")
    }
    if (is.null(names(coef))) names(coef) <- sapply(coef, paste, collapse=";") else named <- TRUE
    names(coef)[names(coef)=="(Intercept)"] <- "Intercept"
  } else coef <- list()

  # handing `contrast`
  design.contr <- NULL
  coef.contr <- NULL
  if (!(missing(contrast) || is.null(contrast))) {
    if (!is.list(contrast)) contrast <- list(contrast)
    if (is.null(names(contrast))) names(contrast) <- sapply(1:length(contrast), function(i) if (is.character(contrast[[i]])) contrast[[i]] else paste0("contrast",i)) else named <- TRUE
    contrast <- lapply(contrast, function(x) {
      if (is.character(x)) {
        if (make.contr || contr.to.coef) {
          makeContrasts(contrasts=x, levels=design, check.names=FALSE) # using my copy of makeContrasts
        } else x
      } else if (is.vector(x) && !is.list(x)) {
        t(t(x))
      } else if (is.matrix(x)) {
        x
      } else stop("Invalid `contrast`, should be a single or a list of character/numeric vectors of matrices.")
    })
    if (contr.to.coef) {
      message("Reforming design matrix with limma::contrastAsCoef such that contrasts become coefficients.")
      tmp <- lapply(contrast, function(x) {
        tmp <- limma::contrastAsCoef(design, x)
        list(design=tmp$design, coef=colnames(tmp$design)[tmp$coef])
      })
      design.contr <- lapply(tmp, function(x) x$design)
      coef.contr <- lapply(tmp, function(x) x$coef)
    }
  } else contrast <- NULL

  # handling `reduced.model`: will convert to coefficients
  if (!(missing(reduced.model) || is.null(reduced.model))) {
    if (!is.list(reduced.model)) reduced.model <- list(reduced.model)
    if (is.null(names(reduced.model))) names(reduced.model) <- sapply(1:length(reduced.model), function(i) if (class(reduced.model[[i]])=="formula") deparse(reduced.model[[i]]) else if (is.matrix(reduced.model[[i]])) paste0("reduced_model",i) else paste0("keep:",paste(reduced.model[[i]],collapse=";"))) else named <- TRUE
    tmp <- lapply(reduced.model, function(x) {
      if (class(x)=="formula") {
        if (missing(pheno) || is.null(pheno)) stop("`reduced.model` is provided as formula but `pheno` is missing, cannot form the design matrix for reduced model.")
        setdiff(colnames(design), colnames(model.matrix(x, pheno)))
      } else if (is.matrix(x)) {
        setdiff(colnames(design), colnames(x))
      } else if (is.character(x)) {
        setdiff(colnames(design), x)
      } else if (is.numeric(x)) {
        colnames(design)[-x]
      } else stop("Invalid `reduced.model`, should be a single or a list of formulae or numeric/character vectors.")
    })
    coef <- c(coef, tmp)
  } else reduced.model <- NULL

  if (length(coef)==0) coef <- NULL

  # if none of `coef`, `contrast`, or `reduced.model` was provided, return all coefficients in the model
  if (is.null(coef) && is.null(contrast)) {
    message("No `coef`, `contrast`, or `reduced.model` was provided, will return all coefficients in the model.")
    coef <- colnames(design)
    coef <- setNames(as.list(coef), coef)
    names(coef)[names(coef)=="(Intercept)"] <- "Intercept"
    named <- TRUE
  }

  # the input `reduced.model` if any has been converted to `coef` earlier;
  # now we reassign output `reduced.model` -- this is all the `coef` available converted back to design matrices of reduced models
  # this is used for the `reduced_design` argument of glmGamPoi::test_de, or for de.deseq2 with LRT test, to be passed to the `reduced` argument of DESeq2::DESeq
  if (!is.null(coef)) {
    reduced.model <- mapply(function(des, x) {
      res <- des[, !colnames(des) %in% x, drop=FALSE]
      if ("contrasts" %in% names(attributes(des))) attr(res, "contrasts") <- attr(des, "contrasts")
      res
    }, list(design), coef, SIMPLIFY=FALSE)
    names(reduced.model) <- names(coef)
  }
  # and if contr.to.coef, convert all the coefs from contrasts into design matrices of reduced models
  if (!is.null(design.contr)) {
    reduced.model.contr <- mapply(function(des, x) {
      res <- des[, !colnames(des) %in% x, drop=FALSE]
      if ("contrasts" %in% names(attributes(des))) attr(res, "contrasts") <- attr(des, "contrasts")
      res
    }, design.contr, contrast, SIMPLIFY=FALSE)
  } else reduced.model.contr <- NULL

  # make.coef.names
  if (!is.null(coef) && make.coef.names) {
    coef <- lapply(coef, function(x) {
      x[x=="(Intercept)"] <- "Intercept"
      make.names(x)
    })
  }

  list(dat=dat, pheno=pheno, model=model, design=design, design.contr=design.contr, coef=coef, coef.contr=coef.contr, contrast=contrast, reduced.model=reduced.model, reduced.model.contr=reduced.model.contr, ccs=ccs, named=named)
}


voom <- function(dat, pheno, model=~., design, quantile=FALSE, ...) {
  # perform limma::voom normalization for RNAseq data
  # dat: gene-by-sample expression matrix of raw counts; should have low genes already filtered out
  # pheno: phenotypic data as a data.frame with the same order of samples; model: the model to be used for downstream DE; together with pheno, this will be used to generate the design matrix passed to limma::voom
  # or provide design matrix in design;
  # if design is NULL, pheno and model will be used to compute the design matrix; otherwise design will be used, and pheno with model will be ignored if provided; need to provide either pheno with model, or design
  # quantile: whether to apply quantile normalization, if TRUE, will pass normalize.method="quantile" to limma::voom; this will be ignored if `normalize.method` is specified in ...
  # ...: passed to edgeR::calcNormFactors (e.g. `method`) and limma::voom
  # return a list(voom, design, genes), where voom is the limma::voom() output, i.e. an EList object, design is the design matrix, and genes is rownames(dat)

  if (missing(design) || is.null(design)) {
    if (missing(pheno) || is.null(pheno)) stop("Need to provide either `pheno` with `model`, or `design`.")
    pheno <- as.data.table(pheno)
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

  list(voom=v, design=design, pheno=pheno, genes=gs)
}

de.limma <- function(dat, pheno=NULL, model=~., design=NULL, coef, contrast, reduced.model, contr.to.coef=FALSE, gene.colname=TRUE, keep.fit=FALSE, ...) {
  # differential expression analysis with limma
  # dat: gene-by-sample expression matrix of raw counts; should have low genes already filtered out; or an ExpressionSet, or output from voom
  # pheno: phenotypic data as a data.table with the same order of samples
  # model: the model to use for DE, by default a linear model containing all variables in pheno (w/o interaction terms)
  # design: design matrix for DE
  # if design is NULL, pheno and model will be used to compute the design matrix; otherwise design will be used, model will be ignored, pheno will not be used to form the design matrix used for DE if provided, but will be used to figure out terms in the reduced model if reduced.model is provided as formula; need to provide either pheno with model, or design
  # coef, contrast and reduced.model are different ways to specify the model coefficients or comparisons for which to test DE and return results, any combination of these can be provided, with each can be provided as single items or lists of items (named lists recommended) for multiple tests; these specifications will be combined in the order of coef, reduced.model, and contrast (with the orders within each kept the same as the input if multiple items were provided), and list of DE result tables in the corresponding order will be returned; if none of these three is provided, will return results for all coefficients in the model; see below for details;
  # coef: numeric or character vector of model coefficients (corresponding to columns of design matrix); if length>1, the coefficient (logFC) of each and a single P value for joint testing (i.e. an ANOVA/drop()-like test, depending on specific methods) will be returned; for multiple tests, provide a list of such items;
  # contrast: numeric contrast vector or matrix (for the latter, one contrast per column), or character vector specifying one or more contrasts in terms of the column names of the design matrix (e.g. sth like "grpA-grpB", in which case it will be converted to contrast vector/matrix with limma::makeContrasts); the matrix case and character vector of length>1 case is handled in the same way as the case of coef with length>1 (i.e. "joint" testing), which may be useful in some cases, but this requires the multiple contrasts are not collinear and probably require contr.to.coef=TRUE for many methods (not tested); for the common use case of separately checking multiple contrasts, pass them in a list of numeric vectors (instead of a single contrast matrix) or a list of atom character vectors (instead of a single length>1 character vector)
  # reduced.model: formula of the reduced model (works only if pheno is provided), or vector of model coefficients (columns of design matrix) to keep (i.e., the opposite to coef, which specifies the coefficients to drop), or design matrix of the reduced model; for multiple tests, provide a list of such items;
  # contr.to.coef: whether to reform the design matrix with limma::contrastAsCoef such that contrasts become coefficients
  # gene.colname: the column name for gene symbols in fData(dat) if dat is an ExpressionSet; if TRUE then will try to get gene symbols automatically from fData(dat); if FALSE then will not try to get gene symbols; gene symbols will be added to the returned DE table
  # keep.fit: if TRUE, then also return the fitted model in addition to the DE result table(s) as list(fit=fit, summary=de.res), otherwise return de.res
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
    pheno <- dat$pheno
    if (!is.null(design)) message("Using the design matrix saved in `dat`, ignoring `pheno` and `model`, if provided.")
  }

  pars <- .process.de.params(dat=mat, pheno=pheno, model=model, design=design, coef=coef, contrast=contrast, reduced.model=reduced.model, contr.to.coef=contr.to.coef)

  # handle pars$coef, i.e. input coef + input reduced.model, there will be just a single model fitting
  if (!is.null(pars$coef)) {
    fit0 <- pass3dots(limma::lmFit, pars$dat, design=pars$design, ...)
    fit <- pass3dots(limma::eBayes, fit0, ...)
    de.res <- lapply(pars$coef, function(x) {
      tryCatch(pass3dots(limma::topTable, fit, coef=x, number=Inf, genelist=gns, ...), error=function(e) pass3dots(limma::topTable, fit, coef=x, number=Inf, ...))
    })
    fit <- list(fit)
  } else {
    fit <- list()
    de.res <- list()
  }
  
  # handle contrast, there can be multiple model fittings, one per contrast
  if (!is.null(pars$contrast)) {
    if (!contr.to.coef) {
      if (!exists("fit0")) fit0 <- pass3dots(limma::lmFit, pars$dat, design=pars$design, ...)
      pars$coef <- lapply(pars$contrast, function(x) if (is.matrix(x)) 1:ncol(x) else 1)
      tmp <- mapply(function(contrast, coef) {
        fit <- limma::contrasts.fit(fit0, contrast)
        fit <- pass3dots(limma::eBayes, fit, ...)
        de.res <- tryCatch(pass3dots(limma::topTable, fit, coef=coef, number=Inf, genelist=gns, ...), error=function(e) pass3dots(limma::topTable, fit, coef=coef, number=Inf, ...))
        list(fit=fit, de.res=de.res)
      }, pars$contrast, pars$coef.contr, SIMPLIFY=FALSE)
      fit <- c(fit, lapply(tmp, function(x) x$fit))
      de.res <- c(de.res, lapply(tmp, function(x) x$de.res))
    } else {
      tmp <- mapply(function(design, coef) {
        fit <- pass3dots(limma::lmFit, pars$dat, design=design, ...)
        fit <- pass3dots(limma::eBayes, fit, ...)
        de.res <- tryCatch(pass3dots(limma::topTable, fit, coef=coef, number=Inf, genelist=gns, ...), error=function(e) pass3dots(limma::topTable, fit, coef=coef, number=Inf, ...))
        list(fit=fit, de.res=de.res)
      }, pars$design.contr, pars$coef.contr, SIMPLIFY=FALSE)
      fit <- c(fit, lapply(tmp, function(x) x$fit))
      de.res <- c(de.res, lapply(tmp, function(x) x$de.res))
    }
  }

  de.res <- lapply(de.res, function(x) {
    if (!"id" %in% tolower(names(x))) x <- cbind(ID=row.names(x), x)
    res <- as.data.table(x)
    setnames(res, c("ID","logFC","AveExpr","P.Value","adj.P.Val"), c("id","log.fc","ave.expr","pval","padj"), skip_absent=TRUE)
    res[order(padj, pval)]
  })

  if (!pars$named && length(de.res)==1) de.res <- de.res[[1]]
  if (length(fit)==1) fit <- fit[[1]]
  if (keep.fit) list(fit=fit, summary=de.res) else de.res
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
  # if design is NULL, pheno and model will be used to compute the design matrix; otherwise design will be used, model will be ignored, pheno will not be used to form the design matrix used for DE if provided, but will be used to figure out terms in the reduced model if reduced.model is provided as formula; need to provide either pheno with model, or design
  # coef, contrast and reduced.model are different ways to specify the model coefficients or comparisons for which to test DE and return results, any combination of these can be provided, with each can be provided as single items or lists of items (named lists recommended) for multiple tests; these specifications will be combined in the order of coef, reduced.model, and contrast (with the orders within each kept the same as the input if multiple items were provided), and list of DE result tables in the corresponding order will be returned; if none of these three is provided, will return results for all coefficients in the model; see below for details;
  # coef: numeric or character vector of model coefficients (corresponding to columns of design matrix); if length>1, the coefficient (logFC) of each and a single P value for joint testing (i.e. an ANOVA/drop()-like test, depending on specific methods) will be returned; for multiple tests, provide a list of such items;
  # contrast: numeric contrast vector or matrix (for the latter, one contrast per column), or character vector specifying one or more contrasts in terms of the column names of the design matrix (e.g. sth like "grpA-grpB", in which case it will be converted to contrast vector/matrix with limma::makeContrasts); the matrix case and character vector of length>1 case is handled in the same way as the case of coef with length>1 (i.e. "joint" testing), which may be useful in some cases, but this requires the multiple contrasts are not collinear and probably require contr.to.coef=TRUE for many methods (not tested); for the common use case of separately checking multiple contrasts, pass them in a list of numeric vectors (instead of a single contrast matrix) or a list of atom character vectors (instead of a single length>1 character vector)
  # reduced.model: formula of the reduced model (works only if pheno is provided), or vector of model coefficients (columns of design matrix) to keep (i.e., the opposite to coef, which specifies the coefficients to drop), or design matrix of the reduced model; for multiple tests, provide a list of such items;
  # contr.to.coef: whether to reform the design matrix with limma::contrastAsCoef such that contrasts become coefficients
  # lfc.cutoff: log fold-change cutoff for returning results
  # keep.fit: if TRUE, then also return the fitted model in addition to the DE result table(s) as list(fit=fit, summary=de.res), otherwise return de.res
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

  # handle pars$coef, i.e. input coef + input reduced.model; only one `fit`
  if (!is.null(pars$coef)) {
    dge <- pass3dots(edgeR::estimateDisp.DGEList, dge, pars$design, ...)
    fit <- list(pass3dots(edgeR::glmQLFit.DGEList, dge, pars$design, ...))
    if (lfc.cutoff==0) {
      de.res <- lapply(pars$coef, function(x) pass3dots(edgeR::glmQLFTest, fit[[1]], coef=x, ...))
    } else {
      de.res <- lapply(pars$coef, function(x) pass3dots(edgeR::glmTreat, fit[[1]], coef=x, lfc=lfc.cutoff, ...))
    }
  } else {
    fit <- list()
    de.res <- list()
  }

  # handle contrast
  if (!is.null(pars$contrast)) {
    # if !contr.to.coef, only one `fit`
    if (!contr.to.coef) {
      if (length(fit)==0) {
        dge <- pass3dots(edgeR::estimateDisp.DGEList, dge, pars$design, ...)
        fit <- list(pass3dots(edgeR::glmQLFit.DGEList, dge, pars$design, ...))
      }
      if (lfc.cutoff==0) {
        de.res <- c(de.res, lapply(pars$contrast, function(x) pass3dots(edgeR::glmQLFTest, fit[[1]], contrast=x, ...)))
      } else {
        de.res <- c(de.res, lapply(pars$contrast, function(x) pass3dots(edgeR::glmTreat, fit[[1]], contrast=x, lfc=lfc.cutoff, ...)))
      }
    } else {
    # if contr.to.coef, multiple `fit`, one per contrast
      tmp <- mapply(function(design, coef) {
        dge <- pass3dots(edgeR::estimateDisp.DGEList, dge, design, ...)
        fit <- pass3dots(edgeR::glmQLFit.DGEList, dge, design, ...)
        if (lfc.cutoff==0) de.res <- pass3dots(edgeR::glmQLFTest, fit, coef=coef, ...) else de.res <- pass3dots(edgeR::glmTreat, fit, coef=coef, lfc=lfc.cutoff, ...)
        list(fit=fit, de.res=de.res)
      }, pars$design.contr, pars$coef.contr, SIMPLIFY=FALSE)
      fit <- c(fit, lapply(tmp, function(x) x$fit))
      de.res <- c(de.res, lapply(tmp, function(x) x$de.res))
    }
  }

  de.res <- lapply(de.res, function(x) {
    tmp <- edgeR::topTags(x, n=Inf)[[1]]
    res <- cbind(id=row.names(tmp), as.data.table(tmp))
    setnames(res, c("logFC","logCPM","PValue","FDR"), c("log.fc","log.cpm","pval","padj"), skip_absent=TRUE)
    setnames(res, stringr::str_replace_all(names(res), "logFC", "log.fc"))
    res[order(padj, pval)]
  })

  if (!pars$named && length(de.res)==1) de.res <- de.res[[1]]
  if (length(fit)==1) fit <- fit[[1]]
  if (keep.fit) list(fit=fit, summary=de.res) else de.res
}


de.deseq2 <- function(dat, pheno, model=~., design, coef, contrast, reduced.model, contr.to.coef=FALSE, test=c("LRT","Wald"), keep.fit=FALSE, nc=1L, ctrl.features, size.factors, ...) {
  # differential expression analysis with DESeq2
  # dat: gene-by-sample expression matrix of raw counts; should have low genes already filtered out
  # pheno: phenotypic data as a data.table with the same order of samples
  # model: the model to use for DE, by default a linear model containing all variables in pheno (w/o interaction terms)
  # design: design matrix for DE
  # if design is NULL, pheno and model will be used to compute the design matrix; otherwise design will be used, model will be ignored, pheno will not be used to form the design matrix used for DE if provided, but will be used to figure out terms in the reduced model if reduced.model is provided as formula; need to provide either pheno with model, or design
  # coef, contrast and reduced.model are different ways to specify the model coefficients or comparisons for which to test DE and return results, any combination of these can be provided, with each can be provided as single items or lists of items (named lists recommended) for multiple tests; these specifications will be combined in the order of coef, reduced.model, and contrast (with the orders within each kept the same as the input if multiple items were provided), and list of DE result tables in the corresponding order will be returned; if none of these three is provided, will return results for all coefficients in the model; see below for details;
  # coef: numeric or character vector of model coefficients (corresponding to columns of design matrix); if length>1, the coefficient (logFC) of each and a single P value for joint testing (i.e. an ANOVA/drop()-like test, depending on specific methods) will be returned; for multiple tests, provide a list of such items;
  # contrast: numeric contrast vector or matrix (for the latter, one contrast per column), or character vector specifying one or more contrasts in terms of the column names of the design matrix (e.g. sth like "grpA-grpB", in which case it will be converted to contrast vector/matrix with limma::makeContrasts); the matrix case and character vector of length>1 case is handled in the same way as the case of coef with length>1 (i.e. "joint" testing), which may be useful in some cases, but this requires the multiple contrasts are not collinear and probably require contr.to.coef=TRUE for many methods (not tested); for the common use case of separately checking multiple contrasts, pass them in a list of numeric vectors (instead of a single contrast matrix) or a list of atom character vectors (instead of a single length>1 character vector)
  # reduced.model: formula of the reduced model (works only if pheno is provided), or vector of model coefficients (columns of design matrix) to keep (i.e., the opposite to coef, which specifies the coefficients to drop), or design matrix of the reduced model; for multiple tests, provide a list of such items;
  # note: coef is passed to the `name` argument of DESeq2::results, while contrast is passed to the `contrast` argument of DESeq2::results; the character vector mode for contrast of DESeq2::results is not supported here (where, e.g. contrast=c("group", "trt", "ctrl") will return results for the 'trt' level compared to 'ctrl' level of the `group` variable)
  # contr.to.coef: whether to reform the design matrix with limma::contrastAsCoef such that contrasts become coefficients
  # test: type of test to perform; here I change the default (as in DESeq2::DESeq) from "Wald" to "LRT"
  # keep.fit: if TRUE, then also return the fitted model in addition to the DE result table(s) as list(fit=fit, summary=de.res), otherwise return de.res
  # nc: number of cores for parallelization
  # ctrl.features: control features/genes that are expected to remain stable across conditions; if provided (i.e. not missing or NULL), will be passed to the `controlGenes` argument of DESeq2::estimateSizeFactors, but unlike `controlGenes`, this can be provided as a character vector of feature/gene symbols (rownames of dat)
  # size.factors: custom size factors for DESeq2, a numeric vector, will directly set sizeFactors(dds) to this if provided; otherwise (i.e. if missing or NULL), will use DESeq2::estimateSizeFactors
  # note: the DESeq2 size factor is directly used for normalization (unlike the edgeR normalization factor, which is multiplied by library size before being used for normalization)
  # ...: passed to DESeq2::DESeq and DESeq2::results

  if (!requireNamespace(c("DESeq2","BiocParallel"), quietly=TRUE)) {
    stop("Packages \"DESeq2\" and \"BiocParallel\" needed for this function to work.")
  }

  test <- match.arg(test)
  if (test=="LRT" && !(missing(contrast) || is.null(contrast)) && !contr.to.coef) {
    message("Provided `contrast` with test=\"LRT\", will enforce contr.to.coef=TRUE.")
    contr.to.coef <- TRUE
  }

  if (nc>1) bp <- BiocParallel::MulticoreParam(workers=nc, progressbar=TRUE, RNGseed=0) else bp <- BiocParallel::bpparam()

  if (any(!is.wholenumber(dat))) {
    message("DESeq2 only accepts integer read counts. There are non-integer values in `dat` and they have been rounded.")
    dat <- round(dat)
  }
  pars <- .process.de.params(dat=dat, pheno=pheno, model=model, design=design, coef=coef, contrast=contrast, reduced.model=reduced.model, contr.to.coef=contr.to.coef, make.coef.names=TRUE)
  if (is.null(pars$pheno)) pars$pheno <- data.table(idx=1:ncol(pars$dat))

  if (!(missing(ctrl.features) || is.null(ctrl.features))) {
    if (is.character(ctrl.features)) ctrl.features <- match(ctrl.features, rownames(pars$dat))
    if (all(is.na(ctrl.features))) stop("None of the provided `ctrl.features` is found in data.")
    if (anyNA(ctrl.features)) message("Some `ctrl.features` not in data.")
    ctrl.features <- ctrl.features[!is.na(ctrl.features)]
  } else ctrl.features <- NULL

  if (!(missing(size.factors) || is.null(size.factors))) {
    if (length(size.factors)==1) size.factors <- rep(size.factors, ncol(pars$dat)) else size.factors <- size.factors[pars$ccs]
  } else size.factors <- NULL

  if (!is.null(pars$coef) || (!is.null(pars$contrast) && test=="Wald" && !contr.to.coef)) {
    dds <- DESeq2::DESeqDataSetFromMatrix(countData=pars$dat, colData=pars$pheno, design=pars$design)
    if (is.null(size.factors)) {
      if (is.null(ctrl.features)) {
        dds <- pass3dots(DESeq2:::estimateSizeFactors.DESeqDataSet, dds, ...)
      } else {
        dds <- pass3dots(DESeq2:::estimateSizeFactors.DESeqDataSet, dds, controlGenes=ctrl.features, ...)
      }
    } else sizeFactors(dds) <- size.factors
  }

  if (test=="LRT") {
    # LRT, multiple `fit`, one per test
    # handle pars$reduced.model, i.e. input coef + input reduced.model combined
    if (!is.null(pars$reduced.model)) {
      fit <- lapply(pars$reduced.model, function(x) pass3dots(DESeq2::DESeq, dds, test=test, reduced=x, parallel=nc>1, BPPARAM=bp, ...))
      de.res <- mapply(function(x,i) pass3dots(DESeq2::results, x, name=i, parallel=nc>1, BPPARAM=bp, ...), fit, pars$coef, SIMPLIFY=FALSE)
    } else {
      fit <- list()
      de.res <- list()
    }
    # handle contrast; with LRT contr.to.coef was enforced
    if (!is.null(pars$contrast)) {
      tmp <- mapply(function(design, reduced.model, coef) {
        dds <- DESeq2::DESeqDataSetFromMatrix(countData=pars$dat, colData=pars$pheno, design=design)
        if (is.null(size.factors)) {
          if (is.null(ctrl.features)) {
            dds <- pass3dots(DESeq2:::estimateSizeFactors.DESeqDataSet, dds, ...)
          } else {
            dds <- pass3dots(DESeq2:::estimateSizeFactors.DESeqDataSet, dds, controlGenes=ctrl.features, ...)
          }
        } else sizeFactors(dds) <- size.factors
        fit <- pass3dots(DESeq2::DESeq, dds, test=test, reduced=reduced.model, parallel=nc>1, BPPARAM=bp, ...)
        de.res <- pass3dots(DESeq2::results, fit, name=coef, parallel=nc>1, BPPARAM=bp, ...)
        list(fit=fit, de.res=de.res)
      }, pars$design.contr, pars$reduced.model.contr, pars$coef.contr, SIMPLIFY=FALSE)
      fit <- c(fit, lapply(tmp, function(x) x$fit))
      de.res <- c(de.res, lapply(tmp, function(x) x$de.res))
    }
  } else if (test=="Wald") {
    # handle pars$coef, i.e. input coef + input reduced.model combined; Wald test, just one `fit`
    if (!is.null(pars$coef)) {
      fit <- list(pass3dots(DESeq2::DESeq, dds, test=test, parallel=nc>1, BPPARAM=bp, ...))
      de.res <- lapply(pars$coef, function(x) pass3dots(DESeq2::results, fit[[1]], name=x, parallel=nc>1, BPPARAM=bp, ...))
    } else {
      fit <- list()
      de.res <- list()
    }
    # handle contrast
    if (!is.null(pars$contrast)) {
      if (!contr.to.coef) {
        # just one `fit`
        if (length(fit)==0) fit <- list(pass3dots(DESeq2::DESeq, dds, test=test, parallel=nc>1, BPPARAM=bp, ...))
        de.res <- c(de.res, lapply(pars$contrast, function(x) pass3dots(DESeq2::results, fit[[1]], contrast=x, parallel=nc>1, BPPARAM=bp, ...)))
      } else {
        # multiple `fit`, one per contrast
        tmp <- mapply(function(design, coef) {
          dds <- DESeq2::DESeqDataSetFromMatrix(countData=pars$dat, colData=pars$pheno, design=design)
          if (is.null(size.factors)) {
            if (is.null(ctrl.features)) {
              dds <- pass3dots(DESeq2:::estimateSizeFactors.DESeqDataSet, dds, ...)
            } else {
              dds <- pass3dots(DESeq2:::estimateSizeFactors.DESeqDataSet, dds, controlGenes=ctrl.features, ...)
            }
          } else sizeFactors(dds) <- size.factors
          fit <- pass3dots(DESeq2::DESeq, dds, test=test, parallel=nc>1, BPPARAM=bp, ...)
          de.res <- pass3dots(DESeq2::results, fit, name=coef, parallel=nc>1, BPPARAM=bp, ...)
          list(fit=fit, de.res=de.res)
        }, pars$design.contr, pars$coef.contr, SIMPLIFY=FALSE)
        fit <- c(fit, lapply(tmp, function(x) x$fit))
        de.res <- c(de.res, lapply(tmp, function(x) x$de.res))
      }
    }
  }

  de.res <- lapply(de.res, function(x) {
    res <- cbind(id=row.names(x), as.data.table(x))
    setnames(res, c("baseMean","log2FoldChange","lfcSE","pvalue"), c("ave.expr","log.fc","lfc.se","pval"), skip_absent=TRUE)
    res[order(padj, pval)]
  })

  if (!pars$named && length(de.res)==1) de.res <- de.res[[1]]
  if (length(fit)==1) fit <- fit[[1]]
  if (keep.fit) list(fit=fit, summary=de.res) else de.res
}


de.glmgampoi <- function(dat, pheno, model=~., design, coef, contrast, reduced.model, contr.to.coef=FALSE, size.factors, pseudobulk=NULL, keep.fit=FALSE, ...) {
  # differential expression analysis with glmGamPoi, this may be used for single-cell RNA-seq data
  # dat: gene-by-sample expression matrix of log-normalized expression value, with gene ID/symbol as rownames and sample ID/barcode as colnames
  # pheno: phenotypic data as a data.table with the same order of samples
  # model: the model to use for DE, by default a linear model containing all variables in pheno (w/o interaction terms)
  # design: design matrix for DE
  # if design is NULL, pheno and model will be used to compute the design matrix; otherwise design will be used, model will be ignored, pheno will not be used to form the design matrix used for DE if provided, but will be used to figure out terms in the reduced model if reduced.model is provided as formula; need to provide either pheno with model, or design
  # coef, contrast and reduced.model are different ways to specify the model coefficients or comparisons for which to test DE and return results, any combination of these can be provided, with each can be provided as single items or lists of items (named lists recommended) for multiple tests; these specifications will be combined in the order of coef, reduced.model, and contrast (with the orders within each kept the same as the input if multiple items were provided), and list of DE result tables in the corresponding order will be returned; if none of these three is provided, will return results for all coefficients in the model; see below for details;
  # coef: numeric or character vector of model coefficients (corresponding to columns of design matrix); if length>1, the coefficient (logFC) of each and a single P value for joint testing (i.e. an ANOVA/drop()-like test, depending on specific methods) will be returned; for multiple tests, provide a list of such items;
  # contrast: numeric contrast vector or matrix (for the latter, one contrast per column), or character vector specifying one or more contrasts in terms of the column names of the design matrix (e.g. sth like "grpA-grpB", in which case it will be converted to contrast vector/matrix with limma::makeContrasts); the matrix case and character vector of length>1 case is handled in the same way as the case of coef with length>1 (i.e. "joint" testing), which may be useful in some cases, but this requires the multiple contrasts are not collinear and probably require contr.to.coef=TRUE for many methods (not tested); for the common use case of separately checking multiple contrasts, pass them in a list of numeric vectors (instead of a single contrast matrix) or a list of atom character vectors (instead of a single length>1 character vector)
  # reduced.model: formula of the reduced model (works only if pheno is provided), or vector of model coefficients (columns of design matrix) to keep (i.e., the opposite to coef, which specifies the coefficients to drop), or design matrix of the reduced model; for multiple tests, provide a list of such items;
  # contr.to.coef: whether to reform the design matrix with limma::contrastAsCoef such that contrasts become coefficients
  # size.factors: passed to glmGamPoi::glm_gp `size_factors`, if missing use the default "normed_sum"
  # pseudobulk: passed to glmGamPoi::test_de `pseudobulk_by`
  # keep.fit: if TRUE, then also return the fitted model in addition to the DE result table(s) as list(fit=fit, summary=de.res), otherwise return de.res
  # nc: number of cores for parallelization
  # ...: passed to glmGamPoi::glm_gp and glmGamPoi::test_de

  if (!requireNamespace("glmGamPoi", quietly=TRUE)) {
    stop("Package \"glmGamPoi\" needed for this function to work.")
  }

  if (missing(size.factors)) size.factors <- "normed_sum"

  pars <- .process.de.params(dat=dat, pheno=pheno, model=model, design=design, coef=coef, contrast=contrast, reduced.model=reduced.model, contr.to.coef=contr.to.coef, make.contr=FALSE)
  # so glmGamPoi::test_de has some issue handling coefficients with special characters in their names, including `(Intercept)` and interaction terms like `a:b`, and the solution is simply to quote them with ``
  if (!is.null(pars$coef)) pars$coef <- lapply(pars$coef, function(x) sprintf("`%s`", x))

  # handle pars$coef (or pars$reduced.model), both correspond to input coef + input reduced.model combined
  if (!is.null(pars$coef)) {
    fit <- list(pass3dots(glmGamPoi::glm_gp, as.matrix(pars$dat), design=pars$design, size_factors=size.factors, ...))
    de.res <- mapply(function(coef, reduced.model) {
      if (length(coef)==1) pass3dots(glmGamPoi::test_de, fit[[1]], contrast=coef, pseudobulk_by=pseudobulk, ...)
        else pass3dots(glmGamPoi::test_de, fit[[1]], reduced_design=reduced.model, pseudobulk_by=pseudobulk, ...)
    }, pars$coef, pars$reduced.model, SIMPLIFY=FALSE)
  } else {
    fit <- list()
    de.res <- list()
  }

  # handle contrast
  if (!is.null(pars$contrast)) {
    if (!contr.to.coef) {
      if (length(fit)==0) fit <- list(pass3dots(glmGamPoi::glm_gp, as.matrix(pars$dat), design=pars$design, size_factors=size.factors, ...))
      de.res <- c(de.res, lapply(pars$contrast, function(x) pass3dots(glmGamPoi::test_de, fit[[1]], contrast=x, pseudobulk_by=pseudobulk, ...)))
    } else {
      tmp <- mapply(function(design, coef, reduced.model) {
        fit <- pass3dots(glmGamPoi::glm_gp, as.matrix(pars$dat), design=design, size_factors=size.factors, ...)
        if (length(coef)==1) de.res <- pass3dots(glmGamPoi::test_de, fit, contrast=coef, pseudobulk_by=pseudobulk, ...)
            else de.res <- pass3dots(glmGamPoi::test_de, fit, reduced_design=reduced.model, pseudobulk_by=pseudobulk, ...)
        list(fit=fit, de.res=de.res)
      }, pars$design.contr, pars$coef.contr, pars$reduced.model.contr, SIMPLIFY=FALSE)
      fit <- c(fit, lapply(tmp, function(x) x$fit))
      de.res <- c(de.res, lapply(tmp, function(x) x$de.res))
    }
  }

  de.res <- lapply(de.res, function(x) {
    res <- as.data.table(x)
    setnames(res, c("name","lfc","f_statistic","adj_pval"), c("id","log.fc","F","padj"), skip_absent=TRUE)
    res[order(padj, pval)]
  })

  if (!pars$named && length(de.res)==1) de.res <- de.res[[1]]
  if (length(fit)==1) fit <- fit[[1]]
  if (keep.fit) list(fit=fit, summary=de.res) else de.res
}


de.mast <- function(dat, pheno, model=~., design, cdr=TRUE, thres=FALSE, coef, contrast, reduced.model, lfc.cutoff=0, pos.only=FALSE, lfc.only=FALSE, nc=1L, keep.fit=FALSE, ...) {
  # differential expression analysis with MAST, used for e.g. single-cell RNA-seq data
  # dat: gene-by-sample expression matrix (assuming sparse) of log-normalized expression value, with gene ID/symbol as rownames and sample ID/barcode as colnames
  # pheno: phenotypic data as a data.table with the same order of samples
  # model: the model to use for DE, by default a linear model containing all variables in pheno (w/o interaction terms)
  # pheno and model will be used to compute the design matrix; or provide design matrix in design; pheno and design cannot both be missing; the design matrix should have proper column names
  # to specify a mixed effect model, need to provide pheno and model; design can only be used to specify common fixed effect models
  # cdr: whether to include cellular detection rate (i.e. fraction of >0 genes in each cell) as a covariate; if TRUE, CDR will be computed as the fraction of >0 values in each column of `dat`
  # thres: not sure whether and how this works yet -- whether to perform adaptive thresholding (see MAST vignette "Using MAST with RNASeq"); if TRUE, will apply the automatic thresholds; or provide custom thresholds
  # coef, contrast and reduced.model are different ways to specify the model coefficients or comparisons for which to test DE and return results, any combination of these can be provided, with each can be provided as single items or lists of items (named lists recommended) for multiple tests; these specifications will be combined in the order of coef, reduced.model, and contrast (with the orders within each kept the same as the input if multiple items were provided), and list of DE result tables in the corresponding order will be returned; if none of these three is provided, will return results for all coefficients in the model; see below for details;
  # coef: character vector of model coefficients (corresponding to columns of design matrix); if length>1, will return p values for the LR test for dropping all of the provided coefs, in this case the logFC returned probably won't be very meaningful; for multiple tests, provide a list of such items;
  # contrast: needs to be either a length-2 character vector or a two-column matrix (or a list of such items for multiple tests);
  #   in the length-2 character vector case, each element specify a contrast in terms of the model coefficients (i.e. column names of the design matrix), the first corresonds to the baseline (e.g. control), the second corresponds to the group of interest (e.g. treated);
  #   in the two-column matrix case, each column is a contrast vecter, the first column corresponds to the baseline (e.g. control), the second corresponds to the group of interest (e.g. treated);
  #   will test the contrast #2-#1; two contrasts are needed (in contrast to common glm where a single contrast vector is sufficient) because only then can the logFC be clearly defined (this is special to the two-component MAST model; for the LRT actually just a single contrast vector is sufficient)
  # reduced.model: formula of the reduced model (works only for fixed effect models and if pheno is provided), or vector of model coefficients (columns of design matrix) to keep (i.e., the opposite to coef, which specifies the coefficients to drop), or design matrix of the reduced model; for multiple tests, provide a list of such items;
  # lfc.cutoff: a non-negative number, lower cutoff for log fold-change (will return genes whose log fold-change >= this value); if set >0, the P values will no longer be valid
  # pos.only: if TRUE, will return genes with log fold-change >= lfc.cutoff; otherwise, return genes with abs(log fold-change) >= lfc.cutoff
  # *note: I include lfc.cutoff and pos.only arguments here (rather than downstream) because I want to filter genes before doing LR test to save time (if p values are needed)
  # lfc.only: if TRUE, will return only log fold-change without p values
  # nc: number of cores to use for multi-core parallelization
  # keep.fit: if TRUE, return list(fit=zlm.fit, zlm.summ=zlm.summ, summary=de.res), where zlm.fit is the output from MAST::zlm, and zlm.summ if the output from summary(zlm.fit); otherwise return de.res
  # ...: passed to MAST::zlm

  if (!requireNamespace("MAST", quietly=TRUE)) {
    stop("Package \"MAST\" (GitHub fork: ImNotaGit/MAST) needed for this function to work.")
  }

  if (!(missing(pheno) || is.null(pheno))) {
    pheno <- as.data.table(pheno)
    vs <- all.vars(model)
    if ("." %in% vs) vs <- unique(c(setdiff(vs, "."), names(pheno))) # this is different from the other de.* functions as MAST supports mixed effect models
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
    if (cdr) {
      pheno <- cbind(pheno, cdr_=Matrix::colMeans(dat>0))
      model <- update(model, ~.+cdr_)
    }
    cdat <- pheno
  }

  if (!(missing(design) || is.null(design))) {
    if (!(missing(pheno) || is.null(pheno))) {
      message("Both `pheno` with `model` and `design` are provided, will use `design`, instead of creating design matrix from `pheno` with `model`")
    }
    colnames(design) <- make.names(colnames(design))
    if (cdr) design <- cbind(design, cdr_=Matrix::colMeans(dat>0))
    #model <- ~.+0 # doesn't work
    model <- as.formula(sprintf("~ %s + 0", paste(sprintf("`%s`", colnames(design)), collapse=" + ")))
    cdat <- design
    ccs <- rep(TRUE, nrow(design))
  } else if (missing(pheno) || is.null(pheno)) stop("Need to provide either `pheno` with `model`, or `design`.")

  sca <- MAST::FromMatrix(exprsArray=as.matrix(dat),
    cData=cbind(data.frame(wellKey=colnames(dat), row.names=colnames(dat)), cdat),
    fData=data.frame(primerid=rownames(dat), row.names=rownames(dat)),
    check_sanity=FALSE)

  if (isTRUE(thres)) {
    thres <- MAST::thresholdSCRNACountMatrix(assay(sca), nbins=20, min_per_bin=30)
    MAST::assays(sca, withDimnames=FALSE) <- c(MAST::assays(sca), list(thresh=thres$counts_threshold))
  }
  
  nc.bak <- options("mc.cores")$mc.cores
  if (nc>1) options(mc.cores=nc)
  zlm.fit <- MAST::zlm(formula=model, sca=sca, parallel=isTRUE(nc>1), ...)
  
  if (missing(design) || is.null(design)) design <- zlm.fit@LMlike@modelMatrix

  # if none of `coef`, `contrast`, or `reduced.model` was provided, return all coefficients in the model
  if ((missing(coef) || is.null(coef)) && (missing(contrast) || is.null(contrast)) && (missing(reduced.model) || is.null(reduced.model))) {
    message("No `coef`, `contrast`, or `reduced.model` was provided, will return all coefficients in the model.")
    coef <- colnames(design)
    coef <- setNames(as.list(coef), coef)
  }

  # input `coef`
  if (!(missing(coef) || is.null(coef))) {
    if (!is.list(coef)) coef <- list(coef)
    if (is.null(names(coef))) names(coef) <- sapply(coef, paste, collapse=";") else named <- TRUE
    names(coef)[names(coef)=="(Intercept)"] <- "Intercept"
  } else coef <- list()

  # input `reduced.model`: convert to coefficients
  if (!(missing(reduced.model) || is.null(reduced.model))) {
    if (!is.list(reduced.model)) reduced.model <- list(reduced.model)
    if (is.null(names(reduced.model))) names(reduced.model) <- sapply(1:length(reduced.model), function(i) if (class(reduced.model[[i]])=="formula") deparse(reduced.model[[i]]) else if (is.matrix(reduced.model[[i]])) paste0("reduced_model",i) else paste0("keep:",paste(reduced.model[[i]],collapse=";"))) else named <- TRUE
    tmp <- lapply(reduced.model, function(x) {
      if (class(x)=="formula") {
        if (missing(pheno) || is.null(pheno)) stop("`reduced.model` is provided as formula but `pheno` is missing, cannot form the design matrix for reduced model.")
        setdiff(colnames(design), make.names(colnames(model.matrix(x, pheno)))) # only supports fixed effect models
      } else if (is.matrix(x)) {
        setdiff(colnames(design), make.names(colnames(x)))
      } else if (is.character(x)) {
        setdiff(colnames(design), make.names(x))
      } else if (is.numeric(x)) {
        colnames(design)[-x]
      } else stop("Invalid `reduced.model`, should be a single or a list of formulae or numeric/character vectors.")
    })
    coef <- c(coef, tmp)
  }

  # handling coef from input coef + coef converted from input reduced.model
  if (length(coef)>0) {
    args.summ <- lapply(coef, function(x) {
      x <- make.names(x)
      c0 <- setNames(ifelse(colnames(design)==make.names("(Intercept)"), 1, 0), colnames(design))
      c1 <- as.matrix(setNames(ifelse(colnames(design) %in% c(x, make.names("(Intercept)")), 1, 0), colnames(design)))
      if (length(x)==1) colnames(c1) <- coef1 <- x else colnames(c1) <- coef1 <- "x"
      # some of the logFC probably won't be meaningful but anyway
      args.getlfc <- list(contrast0=c0, contrast1=c1)
      args.lrt <- list(hypothesis=MAST::CoefficientHypothesis(x))
      list(getlfc=args.getlfc, lrt=args.lrt, coef1=coef1)
    })
  } else args.summ <- list()

  # handling `contrast`
  if (!(missing(contrast) || is.null(contrast))) {
    if (!is.list(contrast)) contrast <- list(contrast)
    if (is.null(names(contrast))) names(contrast) <- sapply(1:length(contrast), function(i) if (is.character(contrast[[i]])) contrast[[i]] else paste0("contrast",i)) else named <- TRUE
    tmp <- lapply(contrast, function(x) {
      if (is.character(x) && length(x)==2) {
        x <- makeContrasts(contrasts=x, levels=design, check.names=FALSE) # using my copy of makeContrasts
        c0 <- x[, 1]
        c1 <- x[, 2, drop=FALSE]
      } else if (is.matrix(x)) {
        c0 <- setNames(x[,1], make.names(rownames(x)))[colnames(design)]
        c1 <- as.matrix(setNames(x[,2], make.names(rownames(x)))[colnames(design)])
      } else stop("Invalid `contrast`, should be a single or a list of: either length-2 character vector or two-column matrix specifying two-group comparison.")
      colnames(c1) <- coef1 <- "x"
      args.getlfc <- list(contrast0=c0, contrast1=c1)
      args.lrt <- list(hypothesis=c1-c0)
      list(getlfc=args.getlfc, lrt=args.lrt, coef1=coef1)
    })
    args.summ <- c(args.summ, tmp)
  }

  if (lfc.cutoff>0 || pos.only) {
    # pre-filtering genes before doing statistical tests to save time
    if (lfc.cutoff>0 && !lfc.only) warning("Filtering genes based on log fold-change, P values no longer valid.")
    if (pos.only) gidx <- do.call(MAST::getLogFC, c(list(zlmfit=zlm.fit), args.summ[[1]]$getlfc))[contrast==args.summ[[1]]$coef1, logFC>=lfc.cutoff] else gidx <- do.call(MAST::getLogFC, c(list(zlmfit=zlm.fit), args.summ[[1]]$getlfc))[contrast==args.summ[[1]]$coef1, abs(logFC)>=lfc.cutoff]
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
  
  res <- lapply(args.summ, function(x) {
    if (lfc.only) {
      zlm.summ <- MAST::summary(zlm.fit, logFC=do.call(MAST::getLogFC, c(list(zlmfit=zlm.fit), x$getlfc)), doLRT=FALSE)$datatable
      if (nc>1) options(mc.cores=nc.bak) # reset
      # discrete component (logistic)
      de.res <- zlm.summ[contrast==x$coef1 & component=='D', .(id=primerid, coef.d=coef, ci95lo.d=ci.lo, ci95up.d=ci.hi)]
      # continuous component
      tmp <- zlm.summ[contrast==x$coef1 & component=='C', .(id=primerid, coef.c=coef, ci95lo.c=ci.lo, ci95up.c=ci.hi)]
      de.res <- merge(de.res, tmp, by="id")
      # logFC
      tmp <- zlm.summ[contrast==x$coef1 & component=='logFC', .(id=primerid, log.fc=coef, ci95lo.lfc=ci.lo, ci95up.lfc=ci.hi)]
      de.res <- merge(de.res, tmp, by="id")
    } else {
      zlm.summ <- MAST::summary(zlm.fit, logFC=do.call(MAST::getLogFC, c(list(zlmfit=zlm.fit), x$getlfc)), doLRT=setNames(list(do.call(MAST::lrTest, c(list(object=zlm.fit), x$lrt))), x$coef1), parallel=isTRUE(nc>1))$datatable
      if (nc>1) options(mc.cores=NULL) # reset
      # discrete component (logistic)
      de.res <- zlm.summ[contrast==x$coef1 & component=='D', .(id=primerid, coef.d=coef, ci95lo.d=ci.lo, ci95up.d=ci.hi, pval.d=`Pr(>Chisq)`)][, padj.d:=p.adjust(pval.d, "BH")]
      # continuous component
      tmp <- zlm.summ[contrast==x$coef1 & component=='C', .(id=primerid, coef.c=coef, ci95lo.c=ci.lo, ci95up.c=ci.hi, pval.c=`Pr(>Chisq)`)][, padj.c:=p.adjust(pval.c, "BH")]
      de.res <- merge(de.res, tmp, by="id")
      # logFC
      tmp <- merge(
        zlm.summ[contrast==x$coef1 & component=='logFC', .(id=primerid, log.fc=coef, ci95lo.lfc=ci.lo, ci95up.lfc=ci.hi)],
        zlm.summ[contrast==x$coef1 & component=='H', .(id=primerid, pval=`Pr(>Chisq)`)][, padj:=p.adjust(pval, "BH")],
        by="id"
      )
      de.res <- merge(de.res, tmp, by="id")
    }
    #de.res <- de.res[order(-abs(log.fc))]
    de.res <- de.res[order(padj, pval)]
    list(summ=zlm.summ, de.res=de.res)
  })
  zlm.summ <- lapply(res, function(x) x$summ)
  de.res <- lapply(res, function(x) x$de.res)
  if (length(de.res)==1) {
    de.res <- de.res[[1]]
    zlm.summ <- zlm.summ[[1]]
  }
  if (keep.fit) list(fit=zlm.fit, zlm.summ=zlm.summ, summary=de.res) else de.res
}


make.pseudobulk <- function(mat, mdat, blk, ncells.cutoff=10) {
  # create pseudobulk gene expression data for single-cell RNA-seq
  # mat: gene-by-cell matrix (assuming sparse)
  # mdat: cell meta data, data.frame or other objects coerceable to data.table
  # blk: names of one or more sample/bulk variables (as in column names of mdat)
  # ncells.cutoff: only keep samples/bulks with >= this number of cells
  # return: list(mat=mat, mdat=mdat), where mat is the pseudobulk gene-by-sample/bulk matrix, mdat is the corresponding meta data for the pseudobulk samples, the latter is obtained by dropping any variables (columns) with non-unique values within a bulk from the original cell-level mdat

  mdat <- as.data.table(mdat)
  tmp <- sapply(blk, function(i) anyNA(mdat[, i, with=FALSE]))
  if (any(tmp)) warning(sprintf("These bulk variables contain NA! NA will be kept as a separate level:\n%s\n", paste(blk[tmp], collapse=", ")))
  blk <- do.call(paste, c(unname(as.list(mdat[, blk, with=FALSE])), list(sep="_")))
  tmp <- table(blk)>=ncells.cutoff
  if (any(!tmp)) {
    warning(sprintf("These bulks contain <%d cells, they will be discarded:\n%s\n", ncells.cutoff, paste(names(tmp)[!tmp], collapse=", ")))
    if (all(!tmp)) stop("All bulks contain <%d cells, no sample is left.", ncells.cutoff)
    idx <- blk %in% names(tmp)[tmp]
    blk <- blk[idx]
    mdat <- mdat[idx]
    mat <- mat[, idx, drop=FALSE]
  }
  blk <- factor(blk)
  tmp <- Matrix::sparseMatrix(i=1:length(blk), j=as.numeric(blk), dims=c(length(blk), length(levels(blk))), dimnames=list(NULL, levels(blk)))
  mat <- as.matrix(mat %*% tmp)
  tmp <- sapply(mdat, function(x) any(colSums(table(x, blk, useNA="ifany")>0)>1))
  if (any(tmp)) message(sprintf("These variables will be dropped since they have multiple values/levels within a bulk:\n%s\n", paste(names(mdat)[tmp], collapse=", ")))
  mdat[, blk:=blk]
  mdat <- unique(mdat[, !tmp, with=FALSE])[match(colnames(mat), blk), -"blk"]
  list(mat=mat, mdat=mdat)
}


