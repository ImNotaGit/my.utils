## ----functions for quick data exploration and visualization----


my.cols <- function(x) {
  # my custom list of easily distinguishable colors for categorical variables with potentially many levels
  # note: order not optimized yet
  cols <- c("#E31A1C", "#295ED4", "#008B00", "#6A3D9A", "#FFFF00", "#FFACFD", "#00FFFF", "#8B4500", "#D9D9D9", "#00FF7F", "#FF1493", "#FFD700", "#7F7F7F", "#66CD00", "#FF7F00", "#FF7D7D", "#ADFF2F", "#AB82FF", "#D2B48C", "#CD853F", "#333333", "#00008B", "#B03060", "#9400D3", "#8B8B00", "#528B8B", "#7EC0EE", "#FFE4C4", "#FF4500", "#CD96CD")
  if (missing(x)) return(cols)
  if (length(x)==1 && is.numeric(x) && is.wholenumber(x)) {
    if (x>length(cols)) stop(sprintf("%d colors are available while you asked for %d", length(cols), x))
    return(cols[1:x])
  } else {
    if (length(x)>length(cols)) stop(sprintf("%d colors are available while you asked for %d.", length(cols), length(x)))
    return(cn1(cols[1:length(x)], x))
  }
}


subchunkify <- function(p, w, h, nm=NULL, ...) {
  # function to create sub-chunks within normal chunks for plots with modified sizes in Rmd documents
  # note: other non-subchunkified plots in the parent chunk may not be produced at the correct position (i.e. messed-up order), so if using subchunkify, then all plots within the parent chunk should be subchunkified;
  # additional chunk options like `message=FALSE`, etc. can be provided in ..., however it seems that in order for this to work, the parent chunk should also have the same options set; subchunks w/o these options can override the parent chunk options
  p.deparsed <- paste0(deparse(function() {p}), collapse="")
  more.args <- deparse(c(...))
  if (more.args=="NULL") more.args <- "" else more.args <- stringr::str_sub(more.args, 3, -2)
  if (is.null(nm)) nm <- sprintf("subchunk_%d", floor(runif(1)*1e6)) else nm <- paste0("subchunk_", nm)
  sub.chunk <- sprintf("\n```{r %s, fig.width=%s, fig.height=%s, echo=FALSE, %s}\n(%s)()\n```\n", nm, w, h, more.args, p.deparsed)
  cat(trimws(knitr::knit(text=knitr::knit_expand(text=sub.chunk), quiet=TRUE)))
}


plot.pca <- function(mat, pc.x=1, pc.y=2, data=NULL, color=NULL, shape=NULL, size=NULL, label=NULL, label.size=3, label.subset=NULL, label.outliers=TRUE, outliers.cutoff=0.01, alpha=0.8, center=TRUE, scale=TRUE, ...) {
  # PCA plot, given a matrix mat of sample-by-variable
  # pc.x and pc.y: PC's to plot on the x any y axes
  # data: data.table with the same number of rows as mat for additional variables of the samples
  # color, shape, size, label: vectors corresponding to the samples (rows of mat) for plotting, or column names in data if data is given
  # label.size: size of label text
  # label.subset: a logical or index vector corresponding to label, or a character vector of a subset of labels; the subset of samples to add label
  # label.outliers: whether or not to specifically label the ourliers in red text (independent of other labeling settings); for the label text, will use value in `label` if given, otherwise use rownames(mat); if rownames(mat)=NULL, will use row indices
  # ourliers.cutoff: the alpha level for determining outliers with a Chi-Squared distribution of the Mahalanobis distances
  # alpha: a single alpha value for the plotting, not used as aes() for plotting
  # center, scale, ...: passed to prcomp()

  if (scale) {
    id <- apply(mat, 2, uniqueN)==1
    s <- sum(id)
    if (s!=0) {
      message(sprintf("removed %d columns of zero variance.", s))
      mat <- mat[, !id]
    }
  }

  res <- prcomp(mat, center=center, scale.=scale, ...)
  tot.var <- sum(res$sdev^2)
  varx <- sprintf("PC %d (%.2f%%)", pc.x, res$sdev[pc.x]^2 /tot.var*100)
  vary <- sprintf("PC %d (%.2f%%)", pc.y, res$sdev[pc.y]^2 /tot.var*100)

  mat.pc <- cbind(x=res$x[, pc.x], y=res$x[, pc.y])
  dat <- data.table(x=res$x[, pc.x], y=res$x[, pc.y])
  if (!is.null(data)) {
    if (nrow(data)!=nrow(dat)) stop("mat and data have different numbers of rows!")
    dat <- cbind(dat, data)
  }
  if (!(length(color)==1 && color %in% names(dat) || is.null(color))) {
    dat[, color:=color]
    color <- "color"
  }
  if (!(length(shape)==1 && shape %in% names(dat) || is.null(shape))) {
    dat[, shape:=shape]
    shape <- "shape"
  }
  if (!(length(size)==1 && size %in% names(dat) || is.null(size))) {
    dat[, size:=size]
    size <- "size"
  }
  if (label.outliers && is.null(label)) {
    if (is.null(rownames(mat))) label <- 1:nrow(mat) else label <- rownames(mat)
  }
  if (!(length(label)==1 && label %in% names(dat) || is.null(label))) {
    dat[, label:=label]
    label <- "label"
  }
  if (!is.null(label.subset)) {
    if (is.logical(label.subset)) label.subset <- which(label.subset)
    if (is.character(label.subset)) label.subset <- which(rownames(mat) %in% label.subset)
  }
  if (label.outliers) {
    md <- mahalanobis(mat.pc, colMeans(mat.pc), cov(mat.pc))
    cutoff <- qchisq(p=1-outliers.cutoff, df=ncol(mat.pc))
    id.outliers <- which(md>cutoff)
    if (!is.null(label.subset)) label.subset <- union(label.subset, id.outliers)
  }

  p <- ggplot(dat, aes(x=x, y=y)) +
    xlab(varx) + ylab(vary) +
    geom_point(aes_string(color=color, shape=shape, size=size), alpha=alpha) +
    theme_classic()
  if (!is.null(label)) {
    if (is.null(label.subset)) dat[, lab.flag:=TRUE] else dat[label.subset, lab.flag:=TRUE]
    dat[, lab.color:="black"]
    if (label.outliers && length(id.outliers)>0) dat[id.outliers, lab.color:="red2"]
    p <- p + geom_text_repel(data=dat[lab.flag==TRUE], aes_string(label=label), color=dat[lab.flag==TRUE, lab.color], size=label.size)
  }
  print(p)

  invisible(list(pca=res, plot.data=dat, plot=p))
}


plot.roc <- function(dat, col="blue4", rev.lgd=FALSE, lgd.tit="theshold", lab=TRUE, lab.size=3.5, lab.posi=c(0.25,0.25)) {
  # plot ROC curve from dat, which is the outpur from get.roc1
  # col: curve color, a single color, or TRUE (varying color by predictor threshold), or a function for transforming the threshold values
  # lgd.tit: title of the legend for color; rev.lgd: whether to reverse legend scale
  # lab: whether to add label of AUROC value and CI; if so, lab.size and lab.posi specify the size and position

  dat.xy <- as.data.table(pROC::coords(dat$roc, "all", transpose=FALSE))[order(-specificity, sensitivity)]

  p <- ggplot(dat.xy) + scale_x_reverse() +
    xlab("Specificity") + ylab("Sensitivity") +
    geom_abline(slope=1, intercept=1, linetype="dashed", alpha=0.7, size=0.2) +
    theme_classic() +
    theme(axis.title.y=element_text(size=14),
      axis.title.x=element_text(size=14),
      axis.text.y=element_text(size=12),
      axis.text.x=element_text(size=12),
      legend.title=element_text(size=12),
      legend.text=element_text(size=12))
  
  if (!is.null(dat$ci)) {
    dat.ci <- data.table(sp=as.numeric(rownames(dat$ci)), se.min=dat$ci[,1], se.max=dat$ci[,3])
    p <- p + geom_ribbon(data=dat.ci, aes(x=sp, ymin=se.min, ymax=se.max), fill="grey50", alpha=0.2)
  }

  if (isTRUE(col) || is.function(col)) {
  	if (is.function(col)) dat.xy[, threshold:=col(threshold)]
  	p <- p + geom_line(aes(x=specificity, y=sensitivity, color=threshold))
  	if (rev.lgd) p <- p + scale_color_viridis_c(name=lgd.tit, direction=-1, guide=guide_colorbar(reverse=TRUE))
  	  else p <- p + scale_color_viridis_c(name=lgd.tit)
  } else {
  	p <- p + geom_line(aes(x=specificity, y=sensitivity), color=col)
  }

  if (lab) {
    lab <- sprintf("AUC=%.3f\n95%% CI:\n(%.3f,%.3f)", dat$auc, dat$auc.ci[1], dat$auc.ci[2])
    p <- p + annotate("text", x=lab.posi[1], y=lab.posi[2], label=lab, size=lab.size)
  }

  return(p)
}


cp.groups <- function(..., ylab="Value", geoms=c("box","violin","jitter"), plab=c(12,23,13), rlab=TRUE, lab.size=4, more.args=list()) {

  # summary grouped data by plotting the groups side-by-side as boxplots (w/ jitter and violin plots), and when there are 2 or 3 groups, print the wilcoxon test p values and r values between each pair of groups in the x axis title.
  # assume the groups of data are given as vectors in ..., or given as a single list of vectors. The first item in ... will be checked, and if it is a list, this single list will be used.
  # geoms: types of figure layers to plot
  # plab: label p value(s) for which pair(s) of comparison; rlab: whether to label rank biserial correlation; lab.size: size of these labels
  # extra arguments to be passed to wilcoxon test should be provided as a list in more.args

  l <- list(...)
  # if the first item in ... (i.e. the first element in l) is a list, assume this is the single object containing the grouped data
  if (is.list(l[[1]])) {l <- l[[1]]}
  # otherwise assume the groups of data are given as vectors in ..., and if they are not named we try to name them by their object name.
  else if (is.null(names(l))) {
    tryCatch({
      tmp <- str_split(deparse(substitute(c(...))), "c\\(|\\, |\\)")[[1]]
      names(l) <- tmp[c(-1, -length(tmp))]
    }, error=function(e) NULL)
  }
  dat <- rbindlist(lapply(l, data.table), idcol="group")
  dat[, group:=factor(group, levels=unique(group))] # specify levels to keep the original order as in ...
  setnames(dat, "V1", "value")
  grps <- levels(dat$group)
  xlabs <- dat[, .(n=.N), by=group][match(grps, group), sprintf("%s\nn=%d", group, n)]

  addm <- function(x, a) (1+a)*max(x[is.finite(x)])+(-a)*min(x[is.finite(x)])
  # if there are only 2 or 3 groups, do wilcoxon test for each pair of groups
  ll <- length(l)
  if (ll==2) {
    tmp <- do.call(wilcox, c(list(value~group, dat), more.args))
    stat <- dat[, .(id=12, x1=1, x2=2, x=1.5, y=addm(value,0.1), p=tmp$pval, r=tmp$r.wilcox)]
  } else if (ll==3) {
    tmp <- do.call(wilcox3, c(list(value~group, dat), more.args))
    stat <- data.table(id=numeric(0), x=numeric(0), y=numeric(0), p=numeric(0), r=numeric(0))
    if (12 %in% plab) stat <- rbind(stat, dat[group %in% levels(group)[1:2], .(id=12, x=1.5, y=addm(value,0.1), p=tmp$pval[1], r=tmp$r.wilcox[1])])
    if (23 %in% plab) stat <- rbind(stat, dat[group %in% levels(group)[2:3], .(id=23, x=2.5, y=addm(value,0.1), p=tmp$pval[2], r=tmp$r.wilcox[2])])
    if (13 %in% plab) {
      if (length(stat$y)==0) a <- 0.1
      else if (max(stat$y)>addm(dat$value,-0.1)) a <- 0.22
      else if (max(stat$y)>addm(dat$value,-0.2)) a <- 0.17
      else a <- 0.1
      stat <- rbind(stat, data.table(id=13, x=2, y=addm(c(stat$y,dat$value),a), p=tmp$pval[3], r=tmp$r.wilcox[3]))
    }
    stat <- merge(data.table(id=c(12,23,13), x1=c(1,2,1), x2=c(2,3,3)), stat, by="id", all=FALSE)
  }
  if (rlab) stat[, lab:=latex2exp::TeX(sprintf("$\\overset{P=%.2g}{r_{rb}=%.2g}$", p, r), output="character")] else stat[, lab:=sprintf("P=%.2g", p)]
  
  # plot summary
  formaty <- function(y) sprintf("%.2g", y)
  p <- ggplot(dat, aes(x=group, y=value)) +
    scale_x_discrete(labels=xlabs) +
    scale_y_continuous(name=ylab, labels=formaty)
  if (any(c("jitter","j") %in% geoms)) {
    if (any(c("violin","v","box","b") %in% geoms)) p <- p + geom_jitter(aes(color=group), size=0.8, width=0.15, height=0.02, alpha=0.4)
    else p <- p +
      geom_jitter(aes(color=group), size=0.8, width=0.2, height=0.02, alpha=0.8) +
      stat_summary(aes(color=group), fun.data=mean_se, geom="pointrange") # plot a line with central dot for mean+/-se
  }
  if (any(c("violin","v") %in% geoms)) {
    if (any(c("box","b") %in% geoms)) p <- p + geom_violin(aes(color=group), scale="width", width=0.7, alpha=0)
    else p <- p + geom_violin(aes(color=group, fill=group), scale="width", width=0.7, alpha=0.3)
  }
  if (any(c("box","b") %in% geoms)) {
    if (any(c("violin","v") %in% geoms)) w <- 0.3 else w <- 0.6
    p <- p + geom_boxplot(aes(color=group), width=w, size=0.8, alpha=0)
  }
  p <- p +
    scale_color_brewer(palette="Set1") +
    scale_fill_brewer(palette="Set1") +
    theme_classic() +
    theme(axis.title.y=element_text(size=15),
      axis.title.x=element_blank(),
      axis.text.y=element_text(size=12),
      axis.text.x=element_text(size=14, hjust=1, angle=35),
      legend.position="none")

  if (ll==2 || ll==3) {
    p <- p + geom_blank(data=dat[, .(group=group[1], y=addm(c(stat$y,value),0.05))], aes(y=y))
    if (rlab) {
      p <- p +
        geom_text(data=stat, aes(x=x, y=y, label=lab), size=lab.size, parse=TRUE) +
        geom_bracket(xmin=stat$x1, xmax=stat$x2, y.position=stat$y, label="", color="grey50")
    } else {
      p <- p + geom_bracket(xmin=stat$x1, xmax=stat$x2, y.position=stat$y, label=stat$lab, label.size=lab.size)
    }
  }

  return(p)
}


mcp.groups <- function(..., ylab="Value", more.args=list()) {

  # a "multiple" version of cp.groups: for the same set of groups, there are multiple comparisons to be done because there are different variables or stratification within the set of groups.
  # same as in cp.groups, summary grouped data by plotting the groups side-by-side as boxplots (w/ jitter and violin plots), and when there are 2 or 3 groups, print the wilcoxon test p values and r values between each pair of groups in the x axis title. The only difference is that the multiple comparisons are further displayed in different panels.
  # assume that each item in ... is a list of vectors/data.frame/data.table, each list-like object contains the data for one group, and the vectors are stratified data (or different variables for that same group). So each item in ... should be a list-like object of the same length (representing the number of strata or variables) and will be compared one-on-one in order.
  # or, ... can be a single list of lists containing all the above data
  # extra arguments to be passed to wilcoxon test should be provided as a list in more.args

  l <- list(...)
  # if the first item in ... (i.e. the first element in l) is a list of lists, assume this is the single object containing the grouped data. we check this by is.list(l[[1]][[1]])
  if (is.list(l[[1]][[1]])) {l <- l[[1]]}
  # otherwise assume the groups of data are given as vectors in ..., and if they are not named we try to name them by their object name.
  else if (is.null(names(l))) {
    tryCatch({
      tmp <- str_split(deparse(substitute(c(...))), "c\\(|\\, |\\)")[[1]]
      names(l) <- tmp[c(-1, -length(tmp))]
    }, error=function(e) NULL)
  }

  # if the lengths of each of ... are not all the same, stop
  if (length(unique(sapply(l, length)))!=1) stop("In mcp.groups: lengths of items provided in ... not all equal.\n")
  # if each of ... (as a list) has different names, stop
  # first check whether each of ... as a list has no names
  if (!is.null(unlist(sapply(l, names)))) {
    if (any(apply(sapply(l, names), 1, function(x) length(unique(x)))!=1)) stop("In mcp.groups: the names of the lists provided in ... don't match exactly.\n")
  }

  dat <- lapply(l, function(onegroup) {
    onegroup.dt <- rbindlist(lapply(onegroup, data.table), idcol="strata")
    onegroup.dt[, strata:=factor(strata, levels=unique(strata))] # specify levels to keep the original order
  })
  dat <- rbindlist(dat, idcol="group")
  dat[, group:=factor(group, levels=unique(group))] # specify levels to keep the original order as in ...
  setnames(dat, "V1", "value")

  datl <- split(dat, dat[, strata])
  # if there are only 2 or 3 groups, do wilcoxon test for each pair of groups, for each stratus
  ll <- length(l)
  if (ll==2) {
    stat.out <- sapply(datl, function(x) {
      stat <- do.call(wilcox, c(list(value~group, x), more.args))
      stat.p <- stat$pval
      stat.r <- stat$r.wilcox
      sprintf("wilcox\np=%.2g\nr=%.2g", stat.p, stat.r)
    })
  } else if (ll==3) {
    stat.out <- sapply(datl, function(x) {
      stat <- do.call(wilcox3, c(list(value~group, x), more.args))
      stat.p <- stat$pval
      stat.r <- stat$r.wilcox
      paste0("wilcox (12,23,13)\np=", paste(sprintf("%.2g", stat.p), collapse="; "), "\nr=", paste(sprintf("%.2g", stat.r), collapse="; "))
    })
  } #else stat.out <- rep("", length(datl))

  if (ll==2 || ll==3) {
    # stat test result data
    stat.out <- data.table(x=mean(1:ll), y=sapply(datl, function(x) x[, 1.12*max(value[is.finite(value)])-0.12*min(value[is.finite(value)])]), strata=factor(names(stat.out), levels=names(stat.out)), s=stat.out) # is.finite will ignore Inf, -Inf, NA and NaN
  }

  # blank data used to adjust axis limits
  blk <- dat[, .(ymax=1.2*max(value[is.finite(value)])-0.2*min(value[is.finite(value)]), group=group), by=strata] # is.finite will ignore Inf, -Inf, NA and NaN

  # plot summary
  p <- ggplot(dat, aes(x=group, y=value)) +
    scale_x_discrete() +
    scale_y_continuous(name=ylab, labels=function(y) sprintf("%.2g", y)) +
    facet_wrap(~strata, scales="free_y") +
    geom_jitter(color="grey", size=1, width=0.15) +
    geom_violin(aes(color=group), scale="width", width=0.6, alpha=0) +
    geom_boxplot(width=0.3, size=0.8, alpha=0) +
    geom_blank(data=blk, aes(y=ymax)) +
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_text(size=15), axis.text.y=element_text(size=12), legend.position="bottom")

  if (ll==2 || ll==3) p <- p + geom_text(data=stat.out, aes(x=x, y=y, label=s), size=4)

  return(p)
}


plot.groups <- function(dat, xvar, yvar, xlab=xvar, ylab=yvar, facet=NULL, geoms=c("box","violin","jitter"), plab=c(12,23,13), rlab=TRUE, lab.size=4, more.args=list()) {
  # cp.groups and mcp.groups, but for the "long format" table as used by ggplot2 by default

  if (!is.factor(dat[[xvar]])) dat[[xvar]] <- factor(dat[[xvar]])

  if (is.null(facet)) {
    ii <- cn(levels(dat[[xvar]]))
    x <- lapply(ii, function(i) dat[[yvar]][dat[[xvar]]==i])
    cp.groups(x, ylab=ylab, geoms=geoms, plab=plab, rlab=rlab, lab.size=lab.size, more.args=more.args)
  } else {
    if (!is.factor(dat[[facet]])) dat[[facet]] <- factor(dat[[facet]])
    ii <- cn(levels(dat[[xvar]]))
    jj <- cn(levels(dat[[facet]]))
    x <- lapply(ii, function(i) {lapply(jj, function(j) dat[[yvar]][dat[[xvar]]==i & dat[[facet]]==j])})
    mcp.groups(x, ylab=ylab, more.args=more.args)
  }
}


plot.xy <- function(x, y, xlab=deparse(substitute(x)), ylab=deparse(substitute(y)), trend="lm", cor.method=c("pearson","spearman","kendall"), density="auto", lab.posi="auto", lab.size=4) {

  q <- qplot(x=x, y=y, xlab=xlab, ylab=ylab) +
    theme_classic() +
    theme(axis.title.y=element_text(size=14),
      axis.title.x=element_text(size=14),
      axis.text.y=element_text(size=12),
      axis.text.x=element_text(size=12))

  if (density=="auto") {
    mat <- table(cut(x,breaks=5), cut(y,breaks=5))
    density <- sum(mat>500)>5
  }
  if (density==TRUE) q <- q + stat_density_2d(aes(color=..level..), geom="polygon", alpha=0) + theme(legend.position="none")

  if ("lm" %in% trend) {
    q <- q + geom_smooth(method=lm, color="blue", size=0.8, fill="blue", alpha=0.2)
    cor.method <- match.arg(cor.method)
    ct <- cor.test(x, y, method=cor.method)
    symb <- switch(cor.method, pearson="r", spearman="rho", kendall="tau")
    p <- ct$p.value
    r <- ct$estimate
    if (p>=2.2e-16) lab <- sprintf("P=%.3g\n%s=%.3f", p, symb, r) else lab <- sprintf("P<2.2e-16\n%s=%.3f", symb, r)
    if (length(lab.posi)==1 && lab.posi=="auto") {
      if (!exists("mat")) mat <- table(cut(x,breaks=5), cut(y,breaks=5))
      tmp <- which(mat==min(mat))
      if (r>0) tmp1 <- match(c(5,10,4,15,9,3,21,22,16,23,17,11,20,14,8,2,24,18,12,6,25,19,13,7,1), tmp)
      else tmp1 <- match(c(25,20,24,15,19,23,1,2,6,3,7,11,10,14,18,22,4,8,12,16,5,9,13,17,21), tmp)
      tmp <- tmp[tmp1[!is.na(tmp1)][1]]
      i <- tmp %% 5
      if (i==0) i <- 5
      j <- (tmp-1) %/% 5 + 1
      lab.posi <- c((1.1-0.2*i)*min(x,na.rm=TRUE)+(0.2*i-0.1)*max(x,na.rm=TRUE), (1.1-0.2*j)*min(y,na.rm=TRUE)+(0.2*j-0.1)*max(y,na.rm=TRUE))
    }
    q <- q + annotate("text", x=lab.posi[1], y=lab.posi[2], label=lab, size=lab.size)
  }
  if ("loess" %in% trend) {
    set.seed(1)
    q <- q + geom_smooth(method=loess, color="red", size=0.8, fill="red", alpha=0.2)
  }
  return(q)
}


plot.pair.corrs <- function(datx, daty, xlab, ylab) {
  lx <- melt(datx, value.name=xlab)
  lx[, variable:=factor(variable, levels=unique(variable))]
  ly <- melt(daty, value.name=ylab)
  dat <- cbind(lx, ly[,!"variable"])

  lmres <- rbindlist(by(dat, dat[,variable], function(d) {
    m <- summary(lm(d[[ylab]]~d[[xlab]]))
    data.table(x=d[, 0.8*min(get(xlab))+0.2*max(get(xlab))], y=d[, 0.8*max(get(ylab))+0.2*min(get(ylab))], txt=sprintf("R2=%.2g\np=%.2g", m$r.squared, m$coefficients["score1","Pr(>|t|)"]))
  }), idcol="variable")
  lmres[, variable:=factor(variable, levels=levels(dat[,variable]))]

  p <- ggplot(dat, aes(x=xlab, y=ylab)) +
    geom_point() +
    facet_wrap(~variable, scales="free") +
    stat_smooth(method="lm") +
    geom_text(data=lmres, aes(x=x, y=y, label=txt), size=4) +
    theme_classic() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_text(size=12), axis.text.y=element_text(size=10), legend.position="bottom")
  return(p)
}


plot.dot <- function(dat, x="odds.ratio", y="gene.set", color="padj", size="overlap.size", xlab=NULL) {
  size1 <- size
  dat <- dat[order(get(x))]
  dat[, c(y):=factor(get(y), levels=get(y))]
  if (is.null(xlab)) xlab <- x
  p <- ggplot(dat, aes(x=get(x), y=get(y))) +
    xlab(xlab) +
    theme_classic() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_text(size=10),
          axis.title.x=element_text(size=12))
  if (!is.null(color) && is.null(size1)) p <- p + geom_point(aes(color=get(color))) + scale_color_continuous(low="red3", high="grey", name=color, guide=guide_colorbar(reverse=TRUE))
  if (!is.null(size1) && is.null(color)) p <- p + geom_point(aes(size=get(size1))) + scale_size_continuous(name=size1)
  if (!is.null(color) && !is.null(size1)) p <- p + geom_point(aes(color=get(color), size=get(size1))) + scale_color_continuous(low="red3", high="grey", name=color, guide=guide_colorbar(reverse=TRUE)) + scale_size_continuous(name=size1)

  return(p)
}


thm <- function(x.tit=NA, x.txt=NA, y.tit=NA, y.txt=NA, tit=NA, face=NA,
	lgd=c(NA,"none","right","bottom","left","top"), lgd.dir=c(NA,"vertical","horizontal"), lgd.box=c(NA,"vertical","horizontal"),
	lgd.tit=NA, lgd.txt=NA, lgd.key=NA, lgd.margin=NA, plt.margin=NA) {
  # shorthand for ggplot2::theme() used for adjusting axes labels, legends, and plot margins
  # NULL means element_blank(), NA means not specified (i.e. default)
  # for arguments corresponding to text elements, can give a single number for text size, or give a named list of parameters, which will be passed to element_text()
  # for lgd.key, give a single number for legend key size in pt
  # for lgd.margin and plt.margin, give a vector of 4 number, representing the margin values in pt for top, right, bottom, and left

  lgd <- match.arg(lgd)
  lgd.dir <- match.arg(lgd.dir)
  lgd.box <- match.arg(lgd.box)

  f <- function(x, u=NULL) {
  	if (is.null(x) || length(x)==0) {
  	  element_blank()
  	} else if (is.numeric(x)) {
  	  if (is.null(u)) element_text(size=x) else unit(x, u)
  	} else if (is.list(x)) {
  	  do.call(element_text, x)
  	} else if (length(x)==1 && is.na(x)) {
  	  NULL
  	} else x
  }

  pars <- list(
  	axis.title.x=f(x.tit),
  	axis.text.x=f(x.txt),
  	axis.title.y=f(y.tit),
    axis.text.y=f(y.txt),
    plot.title=f(tit),
    strip.text=f(face),
    legend.position=f(lgd),
    legend.direction=f(lgd.dir),
    legend.title=f(lgd.tit),
    legend.text=f(lgd.txt),
    legend.key.size=f(lgd.key, "pt"),
    legend.box=f(lgd.box),
    legend.box.margin=f(lgd.margin, "pt"),
    plot.margin=f(plt.margin, "pt")
  )

  pars <- pars[!sapply(pars, is.null)]  
  do.call(theme, pars)
}
