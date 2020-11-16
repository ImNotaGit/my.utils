## ----functions for basic math and statistics----


jaccard <- function(a, b) {
  # jaccard index between two sets
  uniqueN(intersect(a, b)) / uniqueN(union(a, b))
}


qrank <- function(vec) {
  # covert a numeric vector to quantile rank. if there's NA, keep it as it is.
  vec.rank <- rank(vec, na.last="keep", ties.method="average")
  vec.rank / (sum(!is.na(vec.rank)) + 1)
}


inv.norm <- function(vec) {
  # normalize a vector to normal distribution
  qnorm(qrank(vec))
}


trans4m <- function(dat, method=inv.norm, by=21) {
  # Transform a vector, or matrix (or matrix-like object).
  # Arguments:
  ## dat, the data to be transformed, a vector or a matrix (or matrix-like object)
  ## method, the method to be used for transformation, a function. By default, method=inv.norm (i.e. quantile normalization)
  ## by, the way of transformation for a matrix, one of 0, 1, 2, 12, 21. By default, by=21, which means first by column and then by row; by=1 for by row only; by=2 for by column only; by=12 first row then column; by=0 treat the matrix as a (1D) vector.

  if (is.vector(dat) || by==0) {
    return(method(dat))
  } else if (by==21) {
    res <- apply(dat, 2, method)
    res <- apply(res, 1, method)
    res <- t(res)
  } else if (by==1) {
    res <- apply(dat, 1, method)
    res <- t(res)
  } else if (by==2) {
    res <- apply(dat, 2, method)
  } else if (by==12) {
    res <- apply(dat, 1, method)
    res <- apply(res, 1, method) # note that here is still "1", not a typo
  } else {
    stop("Invalid value of argument: by. Should be one of: 0, 1, 2, 12, 21.")
  }
  if (is.matrix(dat)) {
    rownames(res) <- rownames(dat)
    colnames(res) <- colnames(dat)
  } else if (is.data.frame(dat)) {
    rownames(res) <- row.names(dat)
    colnames(res) <- names(dat)
  }

  return(res)
}


wilcox <- function(arg1, arg2=NULL, ...) {
  # run Wilcoxon test on two samples, return a named vector with the p value and the rank-biserial correlation for the test. The two samples can be given as two vectors in arg1 and arg2, but can also be provided in the formula form with the formula in arg1, and an optional data.frame-like object in arg2. The ... is to be passed to wilcox.test() as extra arguments and should be named.

  # NOTE: in case of formula and data.frame, although I can pass them directly to wilcox.test, I choose to extract the two samples explicitly, based on the order of the levels of the group factor (again, although this is the same default behavior of wilcox.test).

  # if arg1 and arg2 are both vectors
  if (is.vector(arg1) && is.vector(arg2)) {
    s1 <- arg1
    s2 <- arg2
  } else if (inherits(arg1, "formula")) {
    if (is.data.frame(arg2)) {
      # if arg1 is a formula and arg2 is a data.frame-like object, assume arg1 is sth as simple as y~x
      group <- factor(arg2[[deparse(arg1[[3]])]]) # make sure group is a factor. the order of levels won't change if it is already a factor
      value <- arg2[[deparse(arg1[[2]])]]
    } else {
      # if arg1: formula, but arg2: not data.frame-like, assume arg1 is sth like dat$y~dat$x
      group <- factor(eval(arg1[[3]])) # make sure group is a factor. the order of levels won't change if it is already a factor
      value <- eval(arg1[[2]])
    }
    # the two samples
    s1 <- value[group==levels(group)[1]]
    s2 <- value[group==levels(group)[2]]
  }

  # run wilcox test
  tryCatch({

    wilcox.res <- wilcox.test(s1, s2, ...)
    # p value
    wilcox.p <- wilcox.res$p.value
    # effect size for paired test, use the rank correlation
    more.args <- list(...)
    if (is.null(more.args$paired) || more.args$paired==FALSE) {
      # effect size for unpaired test: use the rank-biserial correlation (cf. wikipedia). this value is between -1 and +1. positive wilcox.r means s2 larger than s1.
      wilcox.r <- unname(1 - 2 * wilcox.res$statistic / (sum(!is.na(s1))*sum(!is.na(s2))))
    } else if (more.args$paired) {
      ranksum.pos <- wilcox.res$statistic
      nr <- sum(s1!=s2, na.rm=TRUE)
      ranksum <- (1+nr)*nr/2
      ranksum.neg <- ranksum - ranksum.pos
      wilcox.r <- unname((ranksum.neg - ranksum.pos) / ranksum)
    }

    # return
    data.table(r.wilcox=wilcox.r, pval=wilcox.p)

  }, error=function(e) {
    data.table(r.wilcox=NA, pval=NA)
  })

}


wilcox3 <- function(arg1, arg2=NULL, arg3=NULL, ...) {
  # run Wilcoxon test on three samples, return a named vector with the p value and the rank-biserial correlation for the test. The samples can be given as three vectors in arg1-3, but can also be provided in the formula form with the formula in arg1, and an optional data.frame-like object in arg2. The ... is to be passed to wilcox.test() as extra arguments and should be named.

  # as in my wilcox(), extract the two samples explicitly, based on the order of the levels of the group factor.
  # if arg1-3 all vectors
  if (is.vector(arg1) && is.vector(arg2) && is.vector(arg3)) {
    s1 <- arg1
    s2 <- arg2
    s3 <- arg3
  } else if (inherits(arg1, "formula")) {
    if (is.data.frame(arg2)) {
      # if arg1 is a formula and arg2 is a data.frame-like object, assume arg1 is sth as simple as y~x
      group <- factor(arg2[[deparse(arg1[[3]])]]) # make sure group is a factor. the order of levels won't change if it is already a factor
      value <- arg2[[deparse(arg1[[2]])]]
    } else {
      # if arg1: formula, but arg2: not data.frame-like, assume arg1 is sth like dat$y~dat$x
      group <- factor(eval(arg1[[3]])) # make sure group is a factor. the order of levels won't change if it is already a factor
      value <- eval(arg1[[2]])
    }
    # the samples
    gpl <- levels(group)
    s1 <- value[group==gpl[1]]
    s2 <- value[group==gpl[2]]
    s3 <- value[group==gpl[3]]
  }

  # run wilcox.test for all pairs, using my wilcox()
  w12 <- wilcox(s1, s2, ...)
  w23 <- wilcox(s2, s3, ...)
  w13 <- wilcox(s1, s3, ...)
  # labels for group comparisons
  if (exists("gpl")) {
    lab <- c(paste(gpl[1],gpl[2],sep=" vs "), paste(gpl[2],gpl[3],sep=" vs "), paste(gpl[1],gpl[3],sep=" vs "))
  } else lab <- c("12", "23", "13")
  res <- as.data.table(rbind(w12, w23, w13))
  res[, compare:=lab]
  setcolorder(res, c("compare", "r.wilcox", "pval"))
  return(res)
}


adjust.pval <- function(pvaldt, method="BH", in.place=TRUE, key="pval") {
  # adjust the p values in the data.table pvaldt using the specified method. modify pvaldt in place by default. it identifies columns of p values by key="pval" in the column name.

  # if don't want to modify in place, make a copy
  if (!in.place) pvaldt <- copy(pvaldt)
  # get the p value columns and do the adjustment
  padj <- pvaldt[, lapply(.SD, p.adjust, method=method), .SDcols=grep(key, names(pvaldt))]
  # rename these adjusted p value columns from "...pval..." to "...padj..."
  setnames(padj, str_replace(names(padj), key, "padj"))
  # add the padj columns to pvaldt
  pvaldt[, (names(padj)):=padj]
}


make.confus.mat <- function(qset, refset, uset, margins=TRUE) {

  # make a confusion matrix of TP/FP/TN/FN given a reference set (or actual positive set, refset) and a query set (or predicted positive set, qset), with the background being the universal set (uset).
  # set margins=TRUE to add margins

  # if qset is empty or NA, return NA
  if (length(qset)==0) {
    warning("In make.confus.mat: qset has zero length. NA returned.\n")
    return(NA)
  } else if (length(qset)==1 && is.na(qset)) {
    warning("In make.confus.mat: qset is NA. NA returned.\n")
    return(NA)
  }
  # make sure uset has unique items
  uset <- unique(uset)
  # logical vector denoting whether each item is in qset or not
  qsetl <- factor(uset %in% qset, levels=c(TRUE, FALSE), labels=c("Positive","Negative"))
  # logical vector denoting whether each item is in refset or not
  refsetl <- factor(uset %in% refset, levels=c(TRUE, FALSE), labels=c("Positive","Negative"))
  # a confusion matrix as a table
  res <- table(`Query/Prediction`=qsetl, `Reference/Actual`=refsetl)
  # add margins
  if (margins) res <- addmargins(res)
  # return
  return(res)
}


confus.mat.quant <- function(..., index="tpr") {

  # calculate specified quantify related to a confusion matrix

  x <- list(...)
  # if the first item of ... (i.e. x[[1]]) is NA, return NA; if a matrix, assume it is the confusion matrix in the standard form as returned by make.confus.mat, otherwise regard ... as the arguments to be passed to make.confus.mat
  if (length(x[[1]])==0) {
    warning("In confus.mat.quant: first argument has zero length. NA returned.\n")
    return(NA)
  } else if (length(x[[1]])==1 && is.na(x[[1]])) {
    warning("In confus.mat.quant: first argument is NA. NA returned.\n")
    return(NA)
  }
  if (is.matrix(x[[1]])) {
    conf <- addmargins(x[[1]][1:2, 1:2])
  } else conf <- make.confus.mat(..., margins=TRUE)

  tp <- conf["Positive", "Positive"]
  fp <- conf["Positive", "Negative"]
  tn <- conf["Negative", "Negative"]
  fn <- conf["Negative", "Positive"]
  actualp <- conf["Sum", "Positive"]
  actualn <- conf["Sum", "Negative"]
  predp <- conf["Positive", "Sum"]
  predn <- conf["Negative", "Sum"]
  tot <- conf["Sum", "Sum"]
  # result
  switch(tolower(index),
         tp = tp,
         fp = fp,
         tn = tn,
         fn = fn,
         actualp = actualp,
         actualn = actualn,
         predp = predp,
         predn = predn,
         tot = tot,

         precision =,
         tdr =,
         ppv = tp / predp,
         npv = tn / predn,

         sensitivity =,
         recall =,
         tpr = tp / actualp,
         fpr = fp / actualn,
         fdr = fp / predp,

         specificity =,
         spc =,
         tnr = tn / actualn,
         fnr = fn / actualp,

         accuracy =,
         acc = (tp + tn) / tot,
         f1 = (2*tp) / (2*tp + fp + fn),
         mcc = (tp*tn - fp*fn) / sqrt(predp*actualp*predn*actualn)
  )
}


enrich.test <- function(qset=NULL, refset=NULL, uset=NULL, confus.mat=NULL, ...) {

  # Fisher's exact test of enrichment of a reference set (or actual positive set, refset) in a query set (or predicted positive set, qset), with the background being the universal set (uset), or from a given confusion matrix as confus.mat.
  # ... passed to fisher.test()

  if (is.null(confus.mat)) {
    conf <- make.confus.mat(qset, refset, uset, margins=FALSE)
  } else conf <- confus.mat[1:2, 1:2]

  # fisher's exact test
  res <- fisher.test(conf, ...)
  res$table <- conf
  res$summary <- data.table(OR=res$estimate, pval=res$p.value)
  return(res)
}


class.enrich.pval <- function(tab.qset, tab.uset, class.name) {
  # given tables for class composition of the query set and the universal set, and one of the class names, calculate the p value for enrichment of that class in the query set.
  # the qset table should have names %in% that of the uset table, and class.name should be one of the names of the qset table.
  tp <- tab.qset[class.name]
  fp <- sum(tab.qset) - tp
  p <- tab.uset[class.name]
  n <- sum(tab.uset) - p
  fn <- p - tp
  tn <- n - fp
  # confusion matrix
  conf <- matrix(c(tp,fn,fp,tn),2)
  enrich.test(confus.mat=conf)$p.value
}


enrich.gsets <- function(fg, gsets, bg, nc=1L, overlap.cutoff=0, padj.cutoff=1.1) {
  # the over-representation type of enrichment test (with Fisher's exact test here.)
  # fg: query genes; gsets: gene sets as a list object; bg: background genes; overlap.cutoff: only select those gene sets with >this value overlap with fg genes; padj.cutoff: fdr threshold.
  # nc: number of cores

  fg1 <- intersect(fg, bg)
  tmp <- sapply(gsets, function(x) sum(x %in% fg1))
  gsets <- gsets[tmp>overlap.cutoff]

  if (length(fg)==0) {
    warning("The number of query genes is zero, NULL returned.\n")
    return(NULL)
  }

  enrich.gset0 <- function(fg, gset, bg) {
    res <- enrich.test(qset=fg, refset=gset, uset=bg, alternative="greater")
    data.table(overlap.size=res$table[1,1], gene.set.size=res$table[1,1]+res$table[2,1], odds.ratio=res$estimate, pval=res$p.value, overlap.genes=list(unique(intersect(intersect(fg, gset),bg))))
  }
  res <- mclapply(gsets, enrich.gset0, fg=fg, bg=bg, mc.cores=nc)
  res <- rbindlist(res, idcol="gene.set")
  if (ncol(res)==0) return(NULL)
  res[, padj:=p.adjust(pval, method="BH")]
  res <- res[order(padj,pval)][padj<padj.cutoff]
  setcolorder(res, c("gene.set","odds.ratio","pval","padj","gene.set.size","overlap.size","overlap.genes"))
  res
}


gsea <- function(dat, gsets, x="log.fc", id="id", seed=1, ...) {
  # GSEA analysis given a data.frame/data.table, x is the column name of the value, id is the column name of the label of x, gsets is a list of gene sets
  # ... is passed to fgsea::fgsea

  x <- dat[[x]]
  names(x) <- dat[[id]]
  set.seed(seed)
  res <- fgsea::fgsea(pathways=gsets, stats=x, nperm=1e4, ...)
  res[order(padj,pval)]
}


run.f <- function(f, dat, model, coef, ..., drop.test="none", drop=coef, keep.fit) {
  # helper function for constructing the various run* functions

  tmp <- list(...)
  args <- c(list(formula=model, data=dat), tmp[names(tmp) %in% names(formals(f))])

  tryCatch({
    fit <- do.call(f, args)
    tmp <- coef(summary(fit))
    res <- data.table(coef=tmp[coef, colnames(tmp) %in% c("Estimate","coef")], pval=tmp[coef, colnames(tmp) %in% c("Pr(>|t|)","Pr(>|z|)")])
    # "anova-like" test if required
    if (drop.test!="none") {
      tmp <- drop1(fit, as.formula(paste("~",drop)), test=drop.test)
      res[, pval:=tmp[drop, names(tmp) %in% c("Pr(>F)","Pr(>Chi)")]]
    }
    if (keep.fit) list(fitted.model=fit, summary.table=res) else res
  }, error=function(e) {
    warning("Error caught by tryCatch, NA returned: ", e, call.=FALSE, immediate.=TRUE)
    if (keep.fit) list(fitted.model=e, summary.table=data.table(coef=NA, pval=NA)) else data.table(coef=NA, pval=NA)
  })
}


run.lm <- function(dat, model = y ~ x*z, coef="x", ..., drop.test=c("none","Chisq","F"), drop=coef, keep.fit=FALSE) {
  # perform a linear regression (wrapper around lm)
  # dat: a data.table containing covariates; model: formula for regression; coef: for which variable/term should the regression coefficient and p value be returned; the default values are intended to be an example
  # ...: additional variables to be used for fitting the model, if they appear in the model formula and are not in dat; additional arguments to lm() can also be provided here, so need to be careful to avoid any name conflicts
  # drop.test: whether to use p value from a test based on nested models (F test or likelihood ratio test) by dropping the variable of interest as given in drop; if not ("none") the p value will be those from summary(lm()), i.e. based on t test
  # keep.fit: if TRUE, will return list(fitted.model, summary.table), else simply return summary.table, which is a data.table containing the coefficient and p value for the variable/term of interest

  drop.test <- match.arg(drop.test)
  run.f(f=lm, dat=dat, model=model, coef=coef, ..., drop.test=drop.test, drop=drop, keep.fit=keep.fit)
}


run.lmer <- function(dat, model = y ~ x + (x|cluster), coef="x", ..., keep.fit=FALSE) {
  # fit a multilevel linear model (wrapper around lmerTest::lmer)
  # dat: a data.table containing covariates; model: formula for model; coef: for which variable/term should the regression coefficient and p value be returned; the default values are intended to be an example
  # ...: additional variables to be used for fitting the model, if they appear in the model formula and are not in dat; additional arguments to lmer() can also be provided here, so need to be careful to avoid any name conflicts
  # keep.fit: if TRUE, will return list(fitted.model, summary.table), else simply return summary.table, which is a data.table containing the coefficient and p value for the variable/term of interest

  library(lmerTest) # the package was imported but not attached; attach it if this function is called
  run.f(f=lmer, dat=dat, model=model, coef=coef, ..., keep.fit=keep.fit)
}


run.glm <- function(dat, model = y ~ x*z, coef="x", family=binomial, ..., drop.test=c("none","Chisq","F","Rao"), drop=coef, keep.fit=FALSE) {
  # fit a generalized linear model (wrapper around glm); family default to binomial for logistic regression
  # dat: a data.table containing covariates; model: formula for regression; coef: for which variable/term should the regression coefficient and p value be returned; the default values are intended to be an example
  # ...: additional variables to be used for fitting the model, if they appear in the model formula and are not in dat; additional arguments to glm() can also be provided here, so need to be careful to avoid any name conflicts
  # drop.test: whether to use p value from a test based on nested models (e.g. likelihood ratio test, i.e. "Chisq") by dropping the variable of interest as given in drop; if not ("none") the p value will be those from summary(glm()), i.e. based on Wald test that may be problematic for small sample size
  # keep.fit: if TRUE, will return list(fitted.model, summary.table), else simply return summary.table, which is a data.table containing the coefficient and p value for the variable/term of interest

  drop.test <- match.arg(drop.test)
  run.f(f=glm, dat=dat, model=model, coef=coef, family=family, ..., drop.test=drop.test, drop=drop, keep.fit=keep.fit)
}


run.clm <- function(dat, model = y ~ x*z, coef="x", ..., drop.test=c("none","Chisq"), drop=coef, keep.fit=FALSE) {
  # fit a cumulative link model (e.g. ordinal logistic regression) (wrapper around ordinal::clm); the default correspond to an ordinal logistic regression
  # dat: a data.table containing covariates; model: formula for regression; coef: for which variable/term should the regression coefficient and p value be returned; the default values are intended to be an example
  # ...: additional variables to be used for fitting the model, if they appear in the model formula and are not in dat; additional arguments to clm() can also be provided here, so need to be careful to avoid any name conflicts
  # drop.test: whether to use p value from a test based on nested models (e.g. likelihood ratio test, i.e. "Chisq") by dropping the variable of interest as given in drop; if not ("none") the p value will be those from summary(clm()), i.e. based on Wald test that may be problematic for small sample size
  # keep.fit: if TRUE, will return list(fitted.model, summary.table), else simply return summary.table, which is a data.table containing the coefficient and p value for the variable/term of interest

  library(ordinal) # the package was imported but not attached; attach it if this function is called
  drop.test <- match.arg(drop.test)
  run.f(f=clm, dat=dat, model=model, coef=coef, ..., drop.test=drop.test, drop=drop, keep.fit=keep.fit)
}


run.cox <- function(dat, model = Surv(surv_days, surv_status) ~ x + age + strata(gender) + cluster(id), coef="x", ..., drop.test=c("none","Chisq"), drop=coef, keep.fit=FALSE) {
  # perform a Cox regression (wrapper around coxph)
  # dat: a data.table containing covariates; model: formula for Cox regression; coef: for which variable/term should the Cox regression coefficient and p value be returned; the default values are intended to be an example
  # ...: additional variables to be used for fitting the Cox model, if they appear in the model formula and are not in dat; additional arguments to coxph() can also be provided here, so need to be careful to avoid any name conflicts
  # drop.test: whether to use p value from a test based on nested models (e.g. likelihood ratio test, i.e. "Chisq") by dropping the variable of interest as given in drop; if not ("none") the p value will be those from summary(coxph()), i.e. based on Wald test that may be problematic for small sample size
  # keep.fit: if TRUE, will return list(fitted.model, summary.table), else simply return summary.table, which is a data.table containing the coefficient and p value for the variable/term of interest

  library(survival) # the package was imported but not attached; attach it if this function is called
  drop.test <- match.arg(drop.test)
  run.f(f=coxph, dat=dat, model=model, coef=coef, ..., drop.test=drop.test, drop=drop, keep.fit=keep.fit)
}

