## ----functions for pre-processing eset or expression matrix----

# An "eset" is my self-made "expression set" data - a list with $pheno: a data.table of phenotypic data, $expr: a named gene-by-sample gene expression matrix, and $geneid: a data.table of gene id information


select.samples.eset <- function(eset, ids) {
  # select a subset of samples (e.g. patients) from eset. ids should be a vector of sample ids or row indeces (numeric or logic) for eset$pheno.

  # if ids is numeric (row index numbers), no need to convert ids. out of bounds indeces will automatically lead to error
  # if ids is character, assume it's a vector of sample ids
  if (is.character(ids)) {
    # convert ids from symbols to row indeces
    ids <- match(ids, colnames(eset$expr))
    # if all ids are NA, stop
    if (all(is.na(ids))) {
      stop("select.ids.eset: None of the given ids is in the dataset.")
    }
    # if some of the ids NA, give a warning, and select the overlapping ids
    if (any(is.na(ids))) {
      warning("select.ids.eset: Some of the ids are not in the dataset.")
      ids <- ids[!is.na(ids)]
    }
  }
  # if ids is logic, assume it's logic row indeces. out of bounds indeces will automatically lead to error
  if (is.logical(ids)) {
    if (length(ids)<ncol(eset$expr)) stop("select.samples.eset: ids given as logic vector, but with shorter length than the number of samples in the dataset.")
    if (!any(ids)) stop("select.ids.eset: ids given as logic vector, but all its values are FALSE.")
  }
  list(expr=eset$expr[, ids, drop=FALSE], pheno=eset$pheno[ids, ], geneid=eset$geneid)
}

filter.eset <- function(eset, condition) {
  # filter a eset, applying the given condition on eset$pheno.
  e <- substitute(condition)
  ind.logi <- eval(e, eset$pheno, parent.frame())
  # change NA to FALSE
  ind.logi[is.na(ind.logi)] <- FALSE
  # if the result is empty, stop
  if (!any(ind.logi)) stop("filter.eset: The given condition results in empty output.")
  list(expr=eset$expr[, ind.logi, drop=FALSE], pheno=eset$pheno[ind.logi, ], geneid=eset$geneid)
}


rm.pheno.na.eset <- function(eset, pheno.vars=c("surv_status", "surv_days", "gender", "age"), drop.cutoff=0.5, select.vars=FALSE) {
  # remove NAs (i.e. remove samples) in the specified phenotypic variables in eset$pheno. However, if a specified variable has too many NAs (>a fraction=drop.cutoff), the variable will not be used, and if select.cols=TRUE, the variable will not be included in the returned data.

  # avoid changing the original data.table eset$pheno
  this.pheno <- copy(eset$pheno)
  # if any of the vars in pheno.vars are not in names(eset$pheno), give a warning and remove those vars
  vars.avail <- pheno.vars %in% names(this.pheno)
  if (!all(vars.avail)) {
    warning(sprintf("rm.pheno.na.eset: One or more variables in pheno.vars are not present in the data. The variables not present:\n%s.\n", paste(pheno.vars[!vars.avail], collapse=", ")))
    pheno.vars <- pheno.vars[vars.avail]
  }
  # this is specially for my prepared tcga data, where for surv_days (survival time) I regard 0 as NA
  if ("surv_days" %in% pheno.vars) this.pheno[surv_days==0, surv_days:=NA]
  # check if any vars in pheno.vars have too many NAs in their columns
  ncutoff <- nrow(this.pheno) * drop.cutoff
  # logical vector for each var in pheno.vars of whether it has too many NAs
  vars.2many.na <- sapply(pheno.vars, function(x) sum(is.na(this.pheno[[x]])) > ncutoff)
  # print warnings about vars with too many NAs and remove them from pheno.vars
  if (any(vars.2many.na)) {
    warning(sprintf("rm.pheno.na.eset: One or more variables specified in pheno.vars have too many NAs (> number of samples * drop.cutoff), and they will not be used for removing NAs, also they will not be included for select.vars=TRUE:\n%s.\n", paste(pheno.vars[vars.2many.na], collapse=", ")))
    pheno.vars <- pheno.vars[!vars.2many.na]
  }
  # get row indeces for the samples with NAs
  na.row.ind <- unique(unlist(lapply(pheno.vars, function(x) which(is.na(this.pheno[[x]])))))
  # remove NA rows
  this.pheno <- this.pheno[!na.row.ind]
  # if want to select only those columns in pheno.vars, do so
  if (select.vars) this.pheno <- this.pheno[, pheno.vars, with=FALSE]
  list(expr=eset$expr[, -na.row.ind], pheno=this.pheno, geneid=eset$geneid)
}


select.genes.eset <- function(eset, genes) {
  # select a subset of genes from eset. genes should be a vector of gene symbols or row indeces (numeric or logic) for eset$expr.

  # if genes is numeric (row index numbers), no need to convert genes. out of bounds indeces will automatically lead to error
  # if genes is character, assume it's a vector of gene symbols
  if (is.character(genes)) {
    # convert genes from symbols to row indeces
    genes <- match(genes, rownames(eset$expr))
    # if all genes are NA, stop
    if (all(is.na(genes))) {
      stop("select.genes.eset: None of the given genes is in the dataset.")
    }
    # if some of the genes NA, give a warning, and select the overlapping genes
    if (any(is.na(genes))) {
      warning("select.genes.eset: Some of the genes are not in the dataset.")
      genes <- genes[!is.na(genes)]
    }
  }
  # if genes is logic, assume it's logic row indeces. out of bounds indeces will automatically lead to error
  if (is.logical(genes)) {
    if (length(genes)<nrow(eset$expr)) stop("select.genes.eset: Genes given as logic vector, but with shorter length than the number of genes in the dataset.")
    if (!any(genes)) stop("select.genes.eset: Genes given as logic vector, but all its values are FALSE.")
  }
  list(expr=eset$expr[genes, , drop=FALSE], pheno=eset$pheno, geneid=eset$geneid[genes, ])
}

rm.na.genes.eset <- function(eset, rm.frac.cutoff=0.9) {
  # remove from eset those genes that has expression value of NA in >= a certain fraction(=rm.frac.cutoff) of the samples
  nacutoff <- ncol(eset$expr) * rm.frac.cutoff
  ind.logi <- apply(eset$expr, 1, function(x) sum(is.na(x)) < nacutoff)
  select.genes.eset(eset, ind.logi)
}


rm.low.genes.eset <- function(eset, zero.cutoff=0.9, med.cutoff=NULL) {
  # remove from eset those genes that has expression value of 0 in >= a certain fraction(=zero.cutoff) of the samples, and/or those genes with median expression level <= med.cutoff.
  ind.logi.0 <- rep(TRUE, nrow(eset$expr))
  ind.logi.med <- rep(TRUE, nrow(eset$expr))
  if (!is.null(zero.cutoff)) {
    n0cutoff <- ncol(eset$expr) * zero.cutoff
    ind.logi.0 <- apply(eset$expr, 1, function(x) sum(x==0, na.rm=TRUE) < n0cutoff)
  }
  if (!is.null(med.cutoff)) {
    ind.logi.med <- apply(eset$expr, 1, function(x) median(x, na.rm=TRUE) > med.cutoff)
  }
  ind.logi <- ind.logi.0 & ind.logi.med
  select.genes.eset(eset, ind.logi)
}


filter.genes.eset <- function(eset, lb=0, rm.frac.cutoff=0.9) {
  # this may replace the rm.low.genes.eset above, the old one kept for compatibility
  # remove from eset those genes that has expression value <= lb (lower bound) in >= a certain fraction(=rm.frac.cutoff) of the samples
  cutoff <- ncol(eset$expr) * rm.frac.cutoff
  ind.logi <- apply(eset$expr, 1, function(x) sum(x<=lb, na.rm=TRUE) < cutoff)
  select.genes.eset(eset, ind.logi)
}


cbind.eset <- function(eset.x, eset.y, genes="common", gene.id.col="current_id") {
  # combine two esets by sample: assume that these contain non-overlapping samples
  # genes matched by rownames of expr, genes="common" (or "intersect") or "all" (or "union") or "x" or "y"
  # to create $geneid, need to provide the id column of geneid in gene.id.col

  g <- switch(genes,
              all=,
              union=union(rownames(eset.x$expr), rownames(eset.y$expr)),
              common=,
              intersect=intersect(rownames(eset.x$expr), rownames(eset.y$expr)),
              x=rownames(eset.x$expr),
              y=rownames(eset.y$expr))
  gx <- match(g, rownames(eset.x$expr))
  gy <- match(g, rownames(eset.y$expr))
  expr <- cbind(eset.x$expr[gx,], eset.y$expr[gy,])
  rownames(expr) <- g
  tmp <- list(eset.x$pheno, eset.y$pheno)
  names(tmp) <- c(as.character(substitute(eset.x)), as.character(substitute(eset.y)))
  pheno <- rbindlist(tmp, fill=TRUE, idcol="dataset")
  geneid <- rbind(eset.x$geneid, eset.y$geneid, fill=TRUE)[match(g, get(gene.id.col))]
  list(expr=expr, pheno=pheno, geneid=geneid)
}


trans4m.eset <- function(eset, method=inv.norm, by=21) {
  # transform the expression matrix of eset with the specified method and order (i.e. the default 21 means first by column then by row)
  list(expr=trans4m(eset$expr, method=method, by=by), pheno=eset$pheno, geneid=eset$geneid)
}


