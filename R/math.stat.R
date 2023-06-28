## ----functions for basic math and statistics----

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

jaccard <- function(a, b) {
  # jaccard index between two sets
  uniqueN(intersect(a, b)) / uniqueN(union(a, b))
}

dmode <- function(x, ...) {
  # get the mode from the density estimation of a vector of continuous values
  d <- density(x, na.rm=TRUE, ...)
  d$x[which.max(d$y)]
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


wilcox <- function(s1, s2, ...) {
  # run Wilcoxon test on two samples given as two vectors in s1 and s2, return a data.table with the rank-biserial correlation and p value for the test.

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
    message(e)
    message("NA's returned.")
    data.table(r.wilcox=NA_real_, pval=NA_real_)
  })

}


wilcox3 <- function(s1, s2, s3, snames=NULL, ...) {
  # run Wilcoxon test on three samples given as three vectors in s1-s3, return a data.table with the rank-biserial correlations and p values for all pairs of tests
  # snames: if not NULL, should be a length-3 vector of the names for s1-s3

  w12 <- wilcox(s1, s2, ...)
  w23 <- wilcox(s2, s3, ...)
  w13 <- wilcox(s1, s3, ...)
  # labels for group comparisons
  if (!is.null(snames)) {
    lab <- c(paste(snames[1],snames[2],sep=" vs "), paste(snames[2],snames[3],sep=" vs "), paste(snames[1],snames[3],sep=" vs "))
  } else lab <- c("1vs2", "2vs3", "1vs3")
  res <- cbind(comparison=lab, rbind(w12, w23, w13))
  res[, padj:=p.adjust(pval, "BH")]
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


get.roc1 <- function(x, pos, neg, x.names=NULL, dir=c(1,-1), ci=FALSE, msg=TRUE, ...) {
  # get ROC data or AUROC -- wrapper around pROC::roc
  # x is a vector of predictor score, where by default larger score corresponds to positive case (dir==1)
  # x needs to be named by the cases, otherwise provide the names in x.names
  # pos: a vector of the names of positive cases
  # neg: a vector of the names of negative cases
  # ci: whether to compute the 95% CI of sensitivity values used for plotting
  # msg: whether to print message on the AUROC value and its 95% CI
  # return a list(roc=<output of pROC::roc>, ci=<output of pROC::ci.se>, auc=<auc value>, auc.ci=<vector of 2: 95%CI of auc>)
  # ... passed to pROC::roc

  if (!is.null(x.names)) names(x) <- x.names
  else if (is.null(names(x))) stop("x needs to be named.")
  pos <- pos[pos %in% names(x)]
  neg <- neg[neg %in% names(x)]
  dir <- match.arg(as.character(dir[1]), c("1","-1"))
  dir <- switch(dir, `1`="<", `-1`=">")
  roc.obj <- pROC::roc(controls=x[neg], cases=x[pos], direction=dir, ci=TRUE, plot=FALSE, ...)
  if (msg) {
    print(roc.obj$auc)
    print(roc.obj$ci)
  }
  if (ci) {
    message("Computing sensitivity CI values...")
    ci.se.obj <- pROC::ci.se(roc.obj, specificities=seq(0,1,0.025))
  } else ci.se.obj <- NULL

  res <- list(roc=roc.obj, ci=ci.se.obj, auc=roc.obj$auc, auc.ci=roc.obj$ci[c(1,3)])
}


get.roc <- function(x, pos, neg, x.names=NULL, curve=TRUE, ...) {
  # get ROC data or AUROC -- wrapper around PRROC::roc.curve
  # x is a vector of predictor score, where by default larger score corresponds to positive case
  # x needs to be named by the cases, otherwise provide the names in x.names
  # pos: a vector of the names of positive cases
  # neg: a vector of the names of negative cases
  # curve: if false, return only AUROC value, otherwise return the object from PRROC::roc.curve(curve=TRUE)
  # ... passed to PRROC::roc.curve

  if (!is.null(x.names)) names(x) <- x.names
  else if (is.null(names(x))) stop("x needs to be named.")
  pos <- pos[pos %in% names(x)]
  neg <- neg[neg %in% names(x)]
  res <- PRROC::roc.curve(scores.class0=x[pos], scores.class1=x[neg], curve=curve, ...)
  if (curve) res else res$auc
}


get.prc <- function(x, pos, neg, x.names=NULL, curve=TRUE, ...) {
  # get precision-recall curve data or AUPRC -- wrapper around PRROC::pr.curve
  # x is a vector of predictor score, where by default larger score corresponds to positive case
  # x needs to be named by the cases, otherwise provide the names in x.names
  # pos: a vector of the names of positive cases
  # neg: a vector of the names of negative cases
  # curve: if false, return only AUPRC value, otherwise return the object from PRROC::pr.curve(curve=TRUE)
  # ... passed to PRROC::pr.curve

  if (!is.null(x.names)) names(x) <- x.names
  else if (is.null(names(x))) stop("x needs to be named.")
  pos <- pos[pos %in% names(x)]
  neg <- neg[neg %in% names(x)]
  res <- PRROC::pr.curve(scores.class0=x[pos], scores.class1=x[neg], curve=curve, ...)
  res$auc <- res$auc.integral
  if (curve) res else res$auc
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


fisher.simple <- function (x)  {
  # simple function for fisher's test for 2x2 table given in x
  # return list(estimate=odds.ratio, p.value),
  # p.value equivalent to fisher.test(x, alternative="greater", or=1)
  # odds.ratio simply approximated by ad/bc in the 2x2 table

  m <- x[1,1]+x[2,1]
  n <- x[1,2]+x[2,2]
  k <- x[1,1]+x[1,2]
  or <- x[1,1]*x[2,2]/x[1,2]/x[2,1]
  pval <- phyper(x[1,1]-1, m, n, k, lower.tail=FALSE)
  list(estimate=or, p.value=pval)
}


enrich.gsets <- function(fg, gsets, bg, nc=1L, overlap.cutoff=0, padj.cutoff=1.1, name="gene", simple=FALSE) {
  # the over-representation type of enrichment test (with Fisher's exact test here.)
  # fg: query genes; gsets: gene sets as a list object; bg: background genes; overlap.cutoff: only select those gene sets with >this value overlap with fg genes; padj.cutoff: fdr threshold.
  # nc: number of cores
  # can be used for sets of other items rather than genes, `name` is used to modify the column names in output table to reflect this
  # simple: use fisher.simple() for p value, and simple approximate the odds ratio as ad/bc in the 2x2 table

  if (length(fg)==0) {
    warning("The size query set is zero, NULL returned.\n")
    return(NULL)
  }

  # helper function
  enrich.gset <- function(fg, gset, bg, simple) {
    mat <- make.confus.mat(qset=fg, refset=gset, uset=bg, margins=FALSE)
    if (simple) res <- fisher.simple(mat) else res <- fisher.test(mat, alternative="greater", conf.int=FALSE)
    data.table(overlap.size=mat[1,1], set.size=mat[1,1]+mat[2,1], odds.ratio=res$estimate, pval=res$p.value, overlap=list(unique(intersect(intersect(fg, gset),bg))))
  }

  fg1 <- intersect(fg, bg)
  tmp <- sapply(gsets, function(x) sum(x %in% fg1))
  gsets <- gsets[tmp>overlap.cutoff]

  res <- mclapply(gsets, enrich.gset, fg=fg, bg=bg, simple=simple, mc.cores=nc)
  res <- rbindlist(res, idcol="set")
  if (ncol(res)==0) return(NULL)
  res[, padj:=p.adjust(pval, method="BH")]
  res <- res[order(padj,pval)][padj<padj.cutoff]
  setcolorder(res, c("set","odds.ratio","pval","padj","set.size","overlap.size","overlap"))
  setnames(res, c("set","set.size","overlap"), c(paste0(name,".set"), paste0(name,".set.size"), paste0("overlap.",name,"s")))
  res
}


enrich.combo.sets <- function(fg1, fg2, refs1, refs2, bg1, bg2, nc=1L, overlap.cutoff=0, padj.cutoff=1.1, simple=TRUE) {
  # fg1 and fg2: vectors of equal length, all the element-wise pairs {(fg1[i],fg2[i])} formed by these is the "foreground" set; the items in a pair should be different (e.g. A-A should not be present)
  # refs1 and refs2: named lists of set annotations to be used on the 1st and 2nd item respectively, i.e. sth like list(set1=c("x1","x2",...), ...)
  # each pair-wise combination of refs1[[i]] and refs2[[j]] will from a "reference" set, which contain all possible ordered pairs formed by different items in refs1[[i]] and refs2[[j]]
  # bg1 and bg2: vectors representing backgrounds for the 1st and 2nd item; each should contain unique items; all pairs {(bg1[i],bg2[j]) | bg1[i]!=bg2[j]} will form the "background" set
  # nc: number of cores
  # overlap.cutoff: only cases where number of overlap between "foreground" and "reference" set > this value will be kept for P value adjustment
  # padj.cutoff: only cases with BH-adjusted P < this value will be returned
  # simple: use fisher.simple() for p value, and simple approximate the odds ratio as ad/bc in the 2x2 table

  # helper function for one pair of ret sets
  enrich.combo.set <- function(fg1, fg2, ref1, ref2, ref1.name, ref2.name, bg1, bg2, overlap.cutoff, simple) {
    r1 <- ref1[ref1 %in% bg1]
    r2 <- ref2[ref2 %in% bg2]
    
    tmp <- fg1 %in% r1 & fg2 %in% r2
    x11 <- sum(tmp)
    if (x11<=overlap.cutoff) return(NULL)
    overlap <- paste0("(",fg1[tmp],",",fg2[tmp],")")
    x12 <- length(f1) - x11
    ref.n <- length(r1)*length(r2)-sum(r1 %in% r2)  # size of ref
    x21 <- ref.n - x11
    bg.n <- length(bg1)*length(bg2)-sum(bg1 %in% bg2)  # size of bg
    x22 <- bg.n - ref.n - x12
    #           | in ref | not in ref | sum
    # in fg     |  x11   |    x12     | length(fg1) (==length(fg2))
    # not in fg |  x21   |    x22     |
    # sum       |  ref.n |            | bg.n
    mat <- matrix(c(x11,x21,x12,x22), 2)
    # fisher.test and summarize result
    if (simple) res <- fisher.simple(mat) else res <- fisher.test(mat, alternative="greater", conf.int=FALSE)
    data.table(ref.set1=ref1.name, ref.set2=ref2.name, overlap.size=mat[1,1], ref.set.size=mat[1,1]+mat[2,1], overlap.pairs=list(overlap), odds.ratio=res$estimate, pval=res$p.value)  
  }

  tmp <- fg1 %in% bg1 & fg2 %in% bg2
  fg1 <- fg1[tmp]
  fg2 <- fg2[tmp]
  # reduce ref set sizes
  refs1 <- refs1[sapply(refs1, function(x) sum(x %in% fg1)>overlap.cutoff)]
  refs2 <- refs2[sapply(refs2, function(x) sum(x %in% fg2)>overlap.cutoff)]

  cl <- makeForkCluster(nc)
  registerDoParallel(cl)

  res <- foreach(ref1=refs1, ref1n=names(refs1), .combine=rbind) %:%
    foreach(ref2=refs2, ref2n=names(refs2), .combine=rbind) %dopar%
      enrich.combo.set(fg1, fg2, ref1, ref2, ref1n, ref2n, bg1, bg2, overlap.cutoff, simple)
  
  stopCluster(cl)
  stopImplicitCluster()

  if (is.null(res)) return(NULL)
  res[, padj:=p.adjust(pval, method="BH")]
  res <- res[order(padj,pval)][padj<padj.cutoff]
  setcolorder(res, c("ref.set1","ref.set2","odds.ratio","pval","padj","ref.set.size","overlap.size","overlap.pairs"))
  res
}


enrich.combo.sets2 <- function(fg1, fg2, refs1, refs2=refs1, bg1, bg2, nc=1L, overlap.cutoff=0, padj.cutoff=1.1, simple=TRUE) {
  # like enrich.combo.sets, but treating item-pairs as unordered, i.e. as sets, i.e. A-B is the same as B-A
  # fg1 and fg2: vectors of equal length, all element-wise unordered pairs {{fg1[i],fg2[i]}} formed by these is the "foreground" set; if sth like A-A is present, it will be removed; if A-B and B-A are both present, only one will be kept
  # refs1 and refs2: named lists of set annotations to be used on the 1st and 2nd item respectively, i.e. sth like list(set1=c("x1","x2",...), ...)
  # each unordered pair-wise combination of refs1[[i]] and refs2[[j]] (i<=j) will from a "reference" set, which contain all possible unordered pairs formed by different items in refs1[[i]] and refs2[[j]]
  # bg1 and bg2: vectors representing backgrounds for the fg1 and fg2 item; each should contain unique items; all unordered pairs {{bg1[i],bg2[j]} | bg1[i]!=bg2[j]} will form the "background" set
  # nc: number of cores
  # overlap.cutoff: only cases where number of overlap between "foreground" and "reference" set > this value will be kept for P value adjustment
  # padj.cutoff: only cases with BH-adjusted P < this value will be returned
  # simple: use fisher.simple() for p value, and simple approximate the odds ratio as ad/bc in the 2x2 table

  # helper function for one pair of ret sets
  enrich.combo.set <- function(fg1, fg2, ref1, ref2, ref1.name, ref2.name, bg1, bg2, overlap.cutoff, simple) {
    if (ref1.name>ref2.name) return(NULL)
    # x11, x12
    tmp <- (fg1 %in% ref1 & fg2 %in% ref2) | (fg2 %in% ref1 & fg1 %in% ref2)
    x11 <- sum(tmp) # (fg1,fg2) was already filtered to be inside (bg1,bg2), so w/o modifying ref1/2 this can correctly compute x11
    if (x11<=overlap.cutoff) return(NULL)
    overlap <- paste0("(",fg1[tmp],",",fg2[tmp],")")
    x12 <- length(fg1) - x11
    # size of ref, x21
    tmp <- as.data.table(rbind(expand.grid(a=ref1, b=ref2), expand.grid(a=ref2, b=ref1)))
    tmp <- unique(tmp[a<b])
    ref.n <- sum((tmp$a %in% bg1 & tmp$b %in% bg2) | (tmp$b %in% bg1 & tmp$a %in% bg2))
    x21 <- ref.n - x11
    # size of bg, x22
    tmp <- sum(bg1 %in% bg2)
    bg.n <- length(bg1)*length(bg2)- (1+tmp)*tmp/2
    x22 <- bg.n - ref.n - x12
    #           | in ref | not in ref | sum
    # in fg     |  x11   |    x12     | length(fg1) (==length(fg2))
    # not in fg |  x21   |    x22     |
    # sum       |  ref.n |            | bg.n
    mat <- matrix(c(x11,x21,x12,x22), 2)
    # fisher.test and summarize result
    if (simple) res <- fisher.simple(mat) else res <- fisher.test(mat, alternative="greater", conf.int=FALSE)
    data.table(ref.set1=ref1.name, ref.set2=ref2.name, overlap.size=mat[1,1], ref.set.size=mat[1,1]+mat[2,1], overlap.pairs=list(overlap), odds.ratio=res$estimate, pval=res$p.value)  
  }

  # preprocess fg1 and fg2
  tmp <- fg1 %in% bg1 & fg2 %in% bg2 # here won't consider fg2 %in% bg1 & fg1 %in% bg2
  tmp <- data.table(a=c(fg1[tmp],fg2[tmp]), b=c(fg2[tmp],fg1[tmp]))
  tmp <- unique(tmp[a<b])
  fg1 <- tmp$a
  fg2 <- tmp$b
  # reduce ref set sizes
  refs1 <- refs1[sapply(refs1, function(x) sum(x %in% c(fg1,fg2))>overlap.cutoff)]
  refs2 <- refs2[sapply(refs2, function(x) sum(x %in% c(fg1,fg2))>overlap.cutoff)]

  cl <- makeForkCluster(nc)
  registerDoParallel(cl)

  res <- foreach(ref1=refs1, ref1n=names(refs1), .combine=rbind) %:%
    foreach(ref2=refs2, ref2n=names(refs2), .combine=rbind) %dopar%
      enrich.combo.set(fg1, fg2, ref1, ref2, ref1n, ref2n, bg1, bg2, overlap.cutoff, simple)
  
  stopCluster(cl)
  stopImplicitCluster()

  if (is.null(res)) return(NULL)
  res[, padj:=p.adjust(pval, method="BH")]
  res <- res[order(padj,pval)][padj<padj.cutoff]
  setcolorder(res, c("ref.set1","ref.set2","odds.ratio","pval","padj","ref.set.size","overlap.size","overlap.pairs"))
  res
}


gsea <- function(dat, gsets, x="log.fc", id="id", seed=1, nc=1L, ...) {
  # GSEA analysis given a data.frame/data.table, x is the column name of the value, id is the column name of the label of x, gsets is a list of gene sets
  # nc: number of cores
  # ... is passed to fgsea::fgsea

  x <- dat[[x]]
  names(x) <- dat[[id]]
  set.seed(seed)
  if (nc>1) {
    set.seed(seed)
    res <- fgsea::fgsea(pathways=gsets, stats=x, BPPARAM=BiocParallel::MulticoreParam(workers=nc, progressbar=TRUE, RNGseed=seed), ...)
  } else {
    set.seed(seed)
    res <- fgsea::fgsea(pathways=gsets, stats=x, ...)
  }
  res[order(padj,pval)]
}


run.wilcox <- function(dat, model, ...) {
  # run Wilcoxon test for 2 or 3 groups, return a data.table with the p value and the rank-biserial correlation for the test.
  # dat: a data.table; model: formula for the wilcoxon test, e.g. y~x where y contains the data values and x defines the groups

  # NOTE: although I can pass dat and model directly to wilcox.test, I choose to extract the groups explicitly, based on the order of the levels of the group factor (again, although this is the same default behavior of wilcox.test).

  group <- factor(dat[[deparse(model[[3]])]]) # make sure group is a factor. the order of levels won't change if it is already a factor
  value <- dat[[deparse(model[[2]])]]
  grps <- levels(group)
  if (length(grps)==2) {
    s1 <- value[group==grps[1]]
    s2 <- value[group==grps[2]]
    res <- wilcox(s1, s2, ...)
  } else if (length(grps)==3) {
    s1 <- value[group==grps[1]]
    s2 <- value[group==grps[2]]
    s3 <- value[group==grps[3]]
    res <- wilcox3(s1, s2, s3, snames=grps, ...)
  } else {
    stop("The number of groups to be compared should be 2 or 3.")
  }
  res
}


run.cor <- function(dat, model, ...) {
  # run cor.test
  # dat: a data.table; model: formula for the test, e.g. y~x where x and y are the two variables to correlate with each other, note that the order doesn't matter; this format is different from that of cor.test for formula

  y <- dat[[deparse(model[[2]])]]
  x <- dat[[deparse(model[[3]])]]
  tmp <- cor.test(x, y, ...)
  res <- data.table(est=tmp$estimate, pval=tmp$p.value)
  setnames(res, "est", if (names(tmp$estimate)=="cor") "r" else names(tmp$estimate))
  res
}


run.f <- function(f, ..., coef, drop.test="none", drop=coef, keep.fit) {
  # helper function for constructing the various run* functions
  # coef: if missing or NULL, will return results for all coefficients in the model

  tmp <- list(...)
  args <- tmp[names(tmp) %in% names(formals(f))]

  tryCatch({
    fit <- do.call(f, args)
    tmp <- coef(summary(fit))
    if (missing(coef) || is.null(coef)) coef <- rownames(tmp)
    if (is.numeric(coef)) coef <- rownames(tmp)[coef]
    if (any(!coef %in% rownames(tmp))) stop("Some `coef` not in model.")
    res <- data.table(coef=tmp[coef, colnames(tmp) %in% c("Estimate","coef")], se=tmp[coef, colnames(tmp) %in% c("Std. Error","se(coef)")], pval=tmp[coef, colnames(tmp) %in% c("Pr(>|t|)","Pr(>|z|)")])
    if (length(coef)>1) res <- cbind(x=coef, res)
    # "anova-like" test if required
    if (drop.test!="none") {
      tmp <- drop1(fit, as.formula(paste("~",drop)), test=drop.test)
      res[, pval:=tmp[drop, names(tmp) %in% c("Pr(>F)","Pr(>Chi)")]]
    }
    if (keep.fit) list(fit=fit, summary=res) else res
  }, error=function(e) {
    warning("Error caught by tryCatch, NA returned: ", e, call.=FALSE, immediate.=TRUE)
    res <- data.table(coef=NA, se=NA, pval=NA)
    if (exists(coef) && is.vector(coef) && length(coef)>1) res <- cbind(x=coef, res)
    if (keep.fit) list(fit=e, summary=res) else res
  })
}


run.lm <- function(dat, model = y ~ x*z, design=NULL, y=NULL, ..., coef="x", drop.test=c("none","Chisq","F"), drop=coef, keep.fit=FALSE) {
  # perform a linear regression (wrapper around lm)
  # dat: a data.table containing covariates; model: formula for regression;
  # coef: for which variable/term should the regression coefficient and p value be returned;
  # the default values of model and coef are intended to be an example
  # alternatively, provide design (design/model matrix) and y (vector of dependent variable values) instead of model and dat; but can also keep dat and provide y as the name of the dependent variable in dat
  # ...: additional variables to be used for fitting the model, if they appear in the model formula and are not in dat; additional arguments to lm() can also be provided here, so need to be careful to avoid any name conflicts
  # drop.test: whether to use p value from a test based on nested models (F test or likelihood ratio test) by dropping the variable of interest as given in drop; if not ("none") the p value will be those from summary(lm()), i.e. based on t test
  # keep.fit: if TRUE, will return list(fit, summary), else simply return summary, which is a data.table containing the coefficient and p value for the variable/term of interest

  drop.test <- match.arg(drop.test)
  if (is.null(design)) {
    run.f(f=lm, data=dat, formula=model, ..., coef=coef, drop.test=drop.test, drop=drop, keep.fit=keep.fit)
  } else {
    colnames(design) <- make.names(colnames(design))
    design <- as.data.frame(design)
    model <- as.formula(sprintf("`_y` ~ %s + 0", paste(sprintf("`%s`", colnames(design)), collapse=" + ")))
    if (is.character(y)) y <- dat[[y]]
    design <- cbind(`_y`=y, design)
    run.f(f=lm, data=design, formula=model, ..., coef=make.names(coef), drop.test=drop.test, drop=make.names(drop), keep.fit=keep.fit)
  }
}


run.lmer <- function(dat, model = y ~ x + (x|cluster), coef="x", ..., keep.fit=FALSE) {
  # fit a multilevel linear model (wrapper around lmerTest::lmer)
  # dat: a data.table containing covariates; model: formula for model; coef: for which variable/term should the regression coefficient and p value be returned; the default values are intended to be an example
  # ...: additional variables to be used for fitting the model, if they appear in the model formula and are not in dat; additional arguments to lmer() can also be provided here, so need to be careful to avoid any name conflicts
  # keep.fit: if TRUE, will return list(fit, summary), else simply return summary, which is a data.table containing the coefficient and p value for the variable/term of interest

  library(lmerTest) # the package was imported but not attached; attach it if this function is called
  run.f(f=lmer, data=dat, formula=model, ..., coef=coef, keep.fit=keep.fit)
}


run.glm <- function(dat, model = y ~ x*z, design=NULL, y=NULL, family=binomial, ..., coef="x", drop.test=c("none","Chisq","F","Rao"), drop=coef, keep.fit=FALSE) {
  # fit a generalized linear model (wrapper around glm); family default to binomial for logistic regression
  # dat: a data.table containing covariates; model: formula for regression;
  # coef: for which variable/term should the regression coefficient and p value be returned;
  # the default values of model and coef are intended to be an example
  # alternatively, provide design (design/model matrix) and y (vector of dependent variable values) instead of model and dat; but can also keep dat and provide y as the name of the dependent variable in dat
  # ...: additional variables to be used for fitting the model, if they appear in the model formula and are not in dat; additional arguments to glm() can also be provided here, so need to be careful to avoid any name conflicts
  # drop.test: whether to use p value from a test based on nested models (e.g. likelihood ratio test, i.e. "Chisq") by dropping the variable of interest as given in drop; if not ("none") the p value will be those from summary(glm()), i.e. based on Wald test that may be problematic for small sample size
  # keep.fit: if TRUE, will return list(fit, summary), else simply return summary, which is a data.table containing the coefficient and p value for the variable/term of interest

  drop.test <- match.arg(drop.test)
  if (is.null(design)) {
    run.f(f=glm, data=dat, formula=model, family=family, ..., coef=coef, drop.test=drop.test, drop=drop, keep.fit=keep.fit)
  } else {
    if (is.character(y)) y <- dat[[y]]
    run.f(f=glm.fit, x=design, y=y, family=family, ..., coef=coef, drop.test=drop.test, drop=drop, keep.fit=keep.fit)
  }
}


run.clm <- function(dat, model = y ~ x*z, coef="x", ..., drop.test=c("none","Chisq"), drop=coef, keep.fit=FALSE) {
  # fit a cumulative link model (e.g. ordinal logistic regression) (wrapper around ordinal::clm); the default correspond to an ordinal logistic regression
  # dat: a data.table containing covariates; model: formula for regression; coef: for which variable/term should the regression coefficient and p value be returned; the default values are intended to be an example
  # ...: additional variables to be used for fitting the model, if they appear in the model formula and are not in dat; additional arguments to clm() can also be provided here, so need to be careful to avoid any name conflicts
  # drop.test: whether to use p value from a test based on nested models (e.g. likelihood ratio test, i.e. "Chisq") by dropping the variable of interest as given in drop; if not ("none") the p value will be those from summary(clm()), i.e. based on Wald test that may be problematic for small sample size
  # keep.fit: if TRUE, will return list(fit, summary), else simply return summary, which is a data.table containing the coefficient and p value for the variable/term of interest

  library(ordinal) # the package was imported but not attached; attach it if this function is called
  drop.test <- match.arg(drop.test)
  run.f(f=clm, data=dat, formula=model, ..., coef=coef, drop.test=drop.test, drop=drop, keep.fit=keep.fit)
}


run.cox <- function(dat, model = Surv(surv_days, surv_status) ~ x + age + strata(gender) + cluster(id), coef="x", ..., drop.test=c("none","Chisq"), drop=coef, keep.fit=FALSE) {
  # perform a Cox regression (wrapper around coxph)
  # dat: a data.table containing covariates; model: formula for Cox regression; coef: for which variable/term should the Cox regression coefficient and p value be returned; the default values are intended to be an example
  # ...: additional variables to be used for fitting the Cox model, if they appear in the model formula and are not in dat; additional arguments to coxph() can also be provided here, so need to be careful to avoid any name conflicts
  # drop.test: whether to use p value from a test based on nested models (e.g. likelihood ratio test, i.e. "Chisq") by dropping the variable of interest as given in drop; if not ("none") the p value will be those from summary(coxph()), i.e. based on Wald test that may be problematic for small sample size
  # keep.fit: if TRUE, will return list(fit, summary), else simply return summary, which is a data.table containing the coefficient and p value for the variable/term of interest

  library(survival) # the package was imported but not attached; attach it if this function is called
  drop.test <- match.arg(drop.test)
  run.f(f=coxph, data=dat, formula=model, ..., coef=coef, drop.test=drop.test, drop=drop, keep.fit=keep.fit)
}

get.survdiff.pval <- function(survdiff) {
  # get p value from survival::survdiff result
  # this function is copied from fastStat::survdiff_p.value

  if (is.matrix(survdiff$obs)) {
    otmp <- apply(survdiff$obs, 1, sum)
    etmp <- apply(survdiff$exp, 1, sum)
  } else {
    otmp <- survdiff$obs
    etmp <- survdiff$exp
  }
  df <- (etmp > 0)
  if (sum(df) < 2) {
    chi <- 0
    return(1)
  } else {
    temp2 <- ((otmp - etmp)[df])[-1]
    vv <- (survdiff$var[df, df])[-1, -1, drop = FALSE]
    chi <- sum(solve(vv, temp2) * temp2)
    survdiff.pvalue = 1 - pchisq(chi, length(temp2))
    return(survdiff.pvalue)
  }
}

run.survdiff <- function(dat, model = Surv(surv_days, surv_status) ~ x + y + strata(z), coef="x=levelx1, y=levely1", ..., keep.fit=FALSE) {
  # perform a survival::survdiff test
  # dat: a data.table containing covariates; model: formula for survdiff;
  # coef: name of one particular row of the survdiff table; the "observed" and "expected" values in this row will be returned;
  # the default values of model and coef are intended to be examples
  # ...: additional variables for survival::survdiff
  # keep.fit: if TRUE, will return list(fit, summary), else simply return summary, which is a data.table containing the "observed" and "expected" values and p value

  library(survival) # the package was imported but not attached; attach it if this function is called

  res <- tryCatch({
    fit <- survdiff(model, dat, ...)
    i <- names(fit$n)==coef
    if (!any(i)) stop("Invalid coef provided! Should be the name of a row of the survdiff summary table.")
    if (is.matrix(fit$obs)) {
      Exp <- sum(fit$exp[i,])
      obs <- sum(fit$obs[i,])
    } else {
      Exp <- fit$exp[i]
      obs <- fit$obs[i]
    }
    p <- get.survdiff.pval(fit)
    res <- data.table(exp=Exp, obs=obs, pval=p)
    if (keep.fit) list(fit=fit, summary=res) else res
  }, error=function(e) {
    warning("Error caught by tryCatch, NA returned: ", e, call.=FALSE, immediate.=TRUE)
    if (keep.fit) list(fit=e, summary=data.table(coef=NA, pval=NA)) else data.table(exp=NA, obs=NA, pval=NA)
  })
}


run.dirichreg <- function(dat, pheno, model, model.type=c("common","alternative"), base, auto.base=TRUE, auto.base.cutoff=0.05, coef, ..., keep.fit=FALSE) {
  # fit a Dirichlet regression model (wrapper around DirichletReg::DirichReg)
  # dat: a matrix of sample-by-compositional variables; will convert to fractions within each row (i.e. divided by row sums)
  # pheno: a data.table/data.frame containing covariates
  # model: formula for regression (do not include left-hand side)
  # model.type: use common or alternative parameterization for Dirichlet regression
  # base: the base/reference category, either character (corresponding to colnames of dat), or numeric (corresponding to column index of dat), or NULL; only has effect if model.type is alternative
  # auto.base: whether to automatically select the base/reference category; only has effect if base is missing or NULL; if FALSE (and base is missing or NULL), will just use the first category (i.e. base=1)
  # the automatic selection of base follows the scCODA method, i.e. will select the category with the smallest dispersion in fractions across all samples, among the categories that have missing or zero values in less than auto.base.cutoff fraction of samples
  # coef: for which variable/term should the regression coefficient and p value be returned; is missing, return all coefficients
  # ...: additional arguments passed to DirichletReg::DirichReg
  # keep.fit: if TRUE, will return list(fit, summary), else simply return summary, which is a data.table containing the coefficient and p value for the variable/term of interest

  if (!requireNamespace("DirichletReg", quietly=TRUE)) {
    stop("Package \"DirichletReg\" needed for this function to work.")
  }

  library(DirichletReg) # the package was not attached; attach it if this function is called
  
  model.type <- match.arg(model.type)
  cs <- colnames(dat)
  dat <- dat/rowSums(dat, na.rm=TRUE)
  if ("sub.comp" %in% names(list(...))) stop("The `sub.comp` argument is provided. This is not supported yet and the base may not be what is specified.")

  if (model.type=="alternative") {
    if (missing(base) || is.null(base)) {
      if (auto.base) {
        disp <- apply(dat, 2, function(x) var(x)/mean(x))
        idx <- colMeans(dat==0)<auto.base.cutoff
        if (sum(idx)==0) stop(sprintf("No category can be selected as base, may need to increase `auto.base.cutoff` (to at least %d.5/%d)", min(colSums(dat==0)), nrow(dat)))
        base <- which(disp==min(disp[idx]))
        message(sprintf("Automatically selected \"%s\" as the base category.", cs[base]))
      } else {
        base <- 1
        message(sprintf("Using \"%s\" as the base category.", cs[base]))
      }
    } else if (is.character(base)) {
      base <- which(cs==base)
      if (length(base)==0) stop("The specified `base` is not among the colnames of `dat`.")
    } else if (is.numeric(base)) {
      message(sprintf("Using \"%s\" as the base category.", cs[base]))
    }
  } else base <- 1 # if model.type=="common", base has no effect

  dat <- as.data.frame(dat)
  dat$Y <- DR_data(dat, base=base)
  dat <- cbind(dat, as.data.frame(pheno))
  model <- as.formula(paste("Y", paste(as.character(model), collapse=" "))) # cannot simply use update.formula, or the resulting formula will give error in DirichReg
  fit <- do.call(DirichReg, c(list(formula=model, data=dat, model=model.type, base=base), list(...))) # have to use do.call otherwise cannot pass in the formula properly due to internal design of DirichReg
  capture.output(summ <- summary(fit), file="/dev/null")
  # extracting model coefficients and p values following DirichletReg::print.summary_DirichletRegModel
  if (model.type=="common") {
    res <- rbindlist(lapply(setNames(seq_along(summ$varnames), nm=summ$varnames), function(i) {
      tmp <- summ$coef.mat[ifelse(i==1, 1, summ$coef.ind[i-1]+1):summ$coef.ind[i], , drop=FALSE]
      cbind(x.var=rownames(tmp), as.data.table(tmp))
    }), idcol="id")
  } else {
    fit$base.var <- cs[base] # keep record of the name of the base category
    printed.var <- 1
    set.size <- summ$n.vars[1]
    res <- data.table()
    x.vars <- rownames(summ$coef.mat)[printed.var:(printed.var+set.size-1)]
    for (i in seq_along(summ$varnames)) {
      if (i==summ$base) {
        res <- rbind(res, data.table(id=summ$varnames[i], x.var=x.vars), fill=TRUE)
      } else {
        tmp <- summ$coef.mat[printed.var:(printed.var+set.size-1), , drop=FALSE]
        res <- rbind(res, cbind(id=summ$varnames[i], x.var=rownames(tmp), as.data.table(tmp)), fill=TRUE)
        printed.var <- printed.var + set.size
      }
    }
    tmp <- summ$coef.mat[printed.var:length(summ$coefficients), , drop=FALSE]
    res <- rbind(res, cbind(id="precision", x.var=rownames(tmp), as.data.table(tmp)), fill=TRUE)
  }
  setnames(res, c("Estimate","Std. Error","z value", "Pr(>|z|)"), c("coef","se","z","pval"))
  res <- lapply(split(res, by="x.var"), function(x) {
    if (model.type=="common") x[, padj:=p.adjust(pval, "BH")] else x[id!="precision", padj:=p.adjust(pval, "BH")]
    x <- x[order(padj, pval, na.last=FALSE), -"x.var"]
    rbind(x[id=="precision"], x[id!="precision"])
  })
  if (!missing(coef)) {
    if (any(!coef %in% names(res))) {
      warning("Invalid `coef`; return results for all coefficients, please double check.")
    } else res <- res[names(res) %in% coef]
  }
  if (length(res)==1) res <- res[[1]]
  if (keep.fit) list(fit=fit, summary=res) else res
}


simple.nested.model.matrix0 <- function(...) {
  # create a proper design matrix for a simple nested design involving two categorical variables
  # need to provide exactly two named arguments, the first will be referred to as A and the second B to facilitate description below, but in the output both variables will be named by their corresponding argument names;
  # A and B should have the same length w/o NA's, both will be treated as categorical, A should have two (or more) levels, B should be nested within A; it doesn't matter if the level coding of B is reused across levels of A (e.g. if within one level of A the B levels are named 1 and 2, then in another level of A the B levels can still be named 1, 2, etc., or it can be named uniquely w/o name-sharing; in both cases it will be understood that all the B levels are distinct from each other.)
  # control-treatment contrast will be used for A and sum contrast will be used for B
  # a common scenario for this is a control-treatment (i.e. bi-level A) experimental design, with the control and treatment arms each containing multiple independent individuals (i.e. multi-level B; the individuals in the control and treatment groups are different), and there are replicated data points for each individual;
  # we want to create a design matrix whose first non-intercept term corresponds to the effect of treatment vs control: mean(mean of each treated individual) - mean(mean of each control individual)
  # this will return the proper design matrix for the specific case as described above, as a replacement for the base R `model.matrix(~A/B, data)`, since the latter in general does not give what is desired.
  # there should be at least one level of A that contains >=2 levels of B for this function to be useful (although otherwise this function still returns the correct result, i.e. the design matrix reduces to ~A)
  
  args <- list(...)
  if (length(args)!=2 || is.null(names(args))) stop("Please provide two named arguments.")
  a <- factor(args[[1]])
  b <- factor(args[[2]])
  an <- names(args)[1]
  bn <- names(args)[2]
  if (length(a)!=length(b)) stop("The provided variables are not of the same length.")
  if (anyNA(a) || anyNA(b)) stop("Please remove NA's in the provided variables.")
  la <- levels(a)
  lb <- lapply(la, function(ai) unique(b[a==ai]))
  names(lb) <- la
  mat <- matrix(rep(1,length(a)), ncol=1)
  colnames(mat) <- "(Intercept)"
  for (ai in la[-1]) {
    xi <- matrix(ifelse(a==ai, 1, 0), ncol=1)
    colnames(xi) <- paste0(an,ai)
    mat <- cbind(mat, xi)
  }
  for (ai in la) {
    for (bi in lb[[ai]][-1]) {
      xi <- matrix(ifelse(a!=ai, 0, ifelse(b==bi, 1, ifelse(b==lb[[ai]][1], -1, 0))), ncol=1)
      colnames(xi) <- paste0(an,ai,":",bn,bi)
      mat <- cbind(mat, xi)
    }
  }
  rownames(mat) <- 1:nrow(mat)
  tmp <- list("contr.treatment", "contr.sum")
  names(tmp) <- c(an, bn)
  attr(mat, "contrasts") <- tmp
  mat
}


simple.nested.model.matrix <- function(dat, a, b) {
  # create a proper design matrix for a simple nested design involving two categorical variables, the first will be referred to as A and the second B to facilitate description below
  # A and B should have the same length w/o NA's, both will be treated as categorical, A should have two (or more) levels, B should be nested within A; it doesn't matter if the level coding of B is reused across levels of A (e.g. if within one level of A the B levels are named 1 and 2, then in another level of A the B levels can still be named 1, 2, etc., or it can be named uniquely w/o name-sharing; in both cases it will be understood that all the B levels are distinct from each other.)
  # control-treatment contrast will be used for A and sum contrast will be used for B
  # a common scenario for this is a control-treatment (i.e. bi-level A) experimental design, with the control and treatment arms each containing multiple independent individuals (i.e. multi-level B; the individuals in the control and treatment groups are different), and there are replicated data points for each individual;
  # we want to create a design matrix whose first non-intercept term corresponds to the effect of treatment vs control: mean(mean of each treated individual) - mean(mean of each control individual)
  # this will return the proper design matrix for the specific case as described above, as a replacement for the base R `model.matrix(~A/B, dat)`, since the latter in general does not give what is desired.
  # there should be at least one level of A that contains >=2 levels of B for this function to be useful (although otherwise this function still returns the correct result, i.e. the design matrix reduces to ~A)
  # dat: a data.frame or data.table containing the two variables A and B described above (but they can have arbitrary names in dat, and in the output both variables will be named by their orignal names)
  # a, b: character, names of the variables within dat, dat[[a]] and dat[[b]] will be A and B respectively

  args <- setNames(list(dat[[a]], dat[[b]]), c(a, b))
  do.call(simple.nested.model.matrix0, args)
}


makeContrasts <- function(..., contrasts=NULL, levels, check.names=TRUE) {
  #	Construct matrix of custom contrasts
  #	Gordon Smyth
  #	30 June 2003.  Last modified 2 April 2010.
  # this is a copy of limma::makeContrasts from version 3.54.1
  # I added the argument check.names to allow for disabling checking for validity of coefficient names, because it may not matter if I use the resulting contrast matrix for other purposes/in other packages;
  # otherwise, if a design matrix from model.matrix is passed to contrasts, interaction terms like a:b would be invalid and cause this function to stop,

	e <- substitute(list(...))
	if(is.factor(levels)) levels <- levels(levels)
	if(!is.character(levels)) levels <- colnames(levels)
	if(levels[1]=="(Intercept)") {
		levels[1] <- "Intercept"
		warning("Renaming (Intercept) to Intercept")
	}
  if (check.names) {
	  notvalid <- (levels != make.names(levels))
	  if(any(notvalid)) stop("The levels must by syntactically valid names in R, see help(make.names).  Non-valid names: ",paste(levels[notvalid],collapse=","))
  }
	n <- length(levels)
	if(n < 1) stop("No levels to construct contrasts from")
	indicator <- function(i,n) {
		out <- rep(0,n)
		out[i] <- 1
		out
	}
	levelsenv <- new.env()
	for (i in 1:n) assign(levels[i], indicator(i,n), pos=levelsenv)

#	Contrasts given as character vector
	if(!is.null(contrasts)) {
		if(length(e)>1) stop("Can't specify both ... and contrasts")
		e <- as.character(contrasts)
		ne <- length(e)
		cm <- matrix(0,n,ne, dimnames=list(Levels=levels,Contrasts=e))
		if(ne==0) return(cm)
		for (j in 1:ne) {
			ej <- parse(text=e[j])
			cm[,j] <- eval(ej, envir=levelsenv)
		}
		return(cm)
	}

#	Contrasts given as list of expressions
	ne <- length(e)
	enames <- names(e)[2:ne]
	easchar <- as.character(e)[2:ne]
	if(is.null(enames))
		cn <- easchar
	else
		cn <- ifelse(enames=="",easchar,enames)
	cm <- matrix(0,n,ne-1, dimnames=list(Levels=levels,Contrasts=cn))
	if(ne < 2) return(cm)
	for (j in 1:(ne-1)) {
		ej <- e[[j+1]]
		if(is.character(ej)) ej <- parse(text=ej)
		ej <- eval(ej, envir=levelsenv)
#		Character variable
		if(!is.numeric(ej)) {
			colnames(cm)[j] <- as.character(ej)[1]
			if(is.character(ej)) ej <- parse(text=ej)
			ej <- eval(ej, envir=levelsenv)
		}
		cm[,j] <- ej
	}
	cm
}


get.fml.terms <- function(fml) {
  # get all terms in the rhs of a formula, returned as a character vector, e.g. apply this on ~a*b or y~a*b gives c("a","b","a:b")
  attr(terms(fml), "term.labels")
}

get.fml.vars <- function(fml) {
  # get all variables in the rhs of a formula, returned as a character vector, e.g. apply this on ~a*b or y~a*b gives c("a","b")
  res <- rownames(attr(terms(fml), "factors"))
  # if lhs is present, remove lhs
  if (length(fml)==3) res <- res[res!=as.character(fml)[2]]
  res
}

drop.vars <- function(fml, x) {
  # drop any terms involving the specified variables `x` (given as a character vector) from the rhs of formula `fml`
  # this does not work for mixed effect model formula
  ftm <- terms(fml)
  tmp <- attr(ftm, "factors")
  idx <- which(colSums(tmp[x, , drop=FALSE])>0)
  formula(drop.terms(ftm, dropx=idx, keep.response=length(fml)==3))
}


