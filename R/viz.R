## ----functions for quick data exploration and visualization----

plot.pca <- function(mat, pc.x=1, pc.y=2, color=NULL, shape=NULL, label=NULL, center=TRUE, scale=TRUE, ...) {
  # PCA plot, given a matrix mat of sample-by-variable
  # pc.x and pc.y: PC's to plot on the x any y axes
  # color, shape, label: vectors corresponding to the samples (rows of mat) for plotting
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

  p <- qplot(res$x[, pc.x], res$x[, pc.y], xlab=varx, ylab=vary,
    color=color, shape=shape) + theme_classic()

  if (!is.null(label)) p <- p + geom_text_repel(label=label)

  return(p)
}


cp.groups <- function(..., ylab="Value", more.args=list()) {

  # summary grouped data by plotting the groups side-by-side as boxplots (w/ jitter and violin plots), and when there are 2 or 3 groups, print the wilcoxon test p values and r values between each pair of groups in the x axis title.
  # assume the groups of data are given as vectors in ..., or given as a single list of vectors. The first item in ... will be checked, and if it is a list, this single list will be used.
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

  # if there are only 2 or 3 groups, do wilcoxon test for each pair of groups
  ll <- length(l)
  if (ll==2) {
    stat <- do.call(wilcox, c(list(value~group, dat), more.args))
    stat.p <- stat["pval.wilcox"]
    stat.r <- stat["r.wilcox"]
    stat.out <- sprintf("wilcox p = %.2g; r = %.2g", stat.p, stat.r)
  } else if (ll==3) {
    stat <- do.call(wilcox3, c(list(value~group, dat), more.args))
    stat.p <- stat[, pval.wilcox]
    stat.r <- stat[, r.wilcox]
    stat.out <- paste0("wilcox p(12,23,13) = ", paste(sprintf("%.2g", stat.p), collapse="; "), "\nr = ", paste(sprintf("%.2g", stat.r), collapse="; "))
  } else stat.out <- "Groups"

  # plot summary
  formaty <- function(y) sprintf("%.2g", y)
  p <- ggplot(dat, aes(x=group, y=value)) +
    scale_x_discrete(name=stat.out) +
    scale_y_continuous(name=ylab, labels=formaty) +
    geom_jitter(color="grey", size=1, width=0.15, height=0.02) +
    geom_violin(color="blue", scale="width", width=0.6, alpha=0) +
    geom_boxplot(width=0.3, size=0.8, alpha=0) +
    theme_classic() +
    theme(axis.title.y=element_text(size=15), axis.title.x=element_text(size=12), axis.text=element_text(size=12))

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
      stat.p <- stat["pval.wilcox"]
      stat.r <- stat["r.wilcox"]
      sprintf("wilcox\np=%.2g\nr=%.2g", stat.p, stat.r)
    })
  } else if (ll==3) {
    stat.out <- sapply(datl, function(x) {
      stat <- do.call(wilcox3, c(list(value~group, x), more.args))
      stat.p <- stat[, pval.wilcox]
      stat.r <- stat[, r.wilcox]
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


plot.groups <- function(dat, xvar, yvar, xlab=xvar, ylab=yvar, facet=NULL) {

  # ensure dat$xvar and dat$facet are factors
  dat[[xvar]] <- factor(dat[[xvar]])
  if (!is.null(facet)) dat[[facet]] <- factor(dat[[facet]])
  # plot
  p <- ggplot(dat, aes(x=get(xvar), y=get(yvar))) +
    scale_x_discrete(name=xlab) +
    scale_y_continuous(name=ylab, labels=function(y) sprintf("%.2g", y)) +
    geom_jitter(color="grey", size=1, width=0.15, height=0.01) +
    theme_classic()
  if (is.null(facet)) {
    p <- p +
      geom_violin(color="blue", scale="width", width=0.6, alpha=0) +
      geom_boxplot(width=0.3, size=0.8, alpha=0) +
      theme(axis.title.y=element_text(size=12), axis.text.y=element_text(size=10), axis.title.x=element_text(size=12), axis.text.x=element_text(size=10, hjust=1, angle=30))
  } else {
    p <- p + facet_wrap(~get(facet), scales="free_y") +
      geom_violin(aes(color=get(xvar)), scale="width", width=0.6, alpha=0) +
      labs(color=xlab) + # this is to change lengend title
      geom_boxplot(width=0.3, size=0.8, alpha=0) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_text(size=12), axis.text.y=element_text(size=10), legend.position="bottom")
  }

  return(p)
}


plot.xy <- function(x, y, trend="lm", cor.method="pearson", xlab=deparse(substitute(x)), ylab=deparse(substitute(y))) {

  q <- qplot(x=x, y=y, xlab=xlab, ylab=ylab) + theme_classic()
  if ("lm" %in% trend) {
    q <- q + geom_smooth(method=lm, color="blue", size=0.8, fill="blue", alpha=0.2)
    ct <- cor.test(x, y, method=cor.method)
    symb <- switch(cor.method, pearson="r", spearman="rho", kendall="tau")
    p <- ct$p.value
    r <- ct$estimate
    if (r>0) {
      q <- q + geom_text(aes(x=0.9*min(x,na.rm=TRUE)+0.1*max(x,na.rm=TRUE), y=0.1*min(y,na.rm=TRUE)+0.9*max(y,na.rm=TRUE), label=sprintf("%s = %.3f\nP = %.3g", symb, r, p)))
    } else {
      q <- q + geom_text(aes(x=0.1*min(x,na.rm=TRUE)+0.9*max(x,na.rm=TRUE), y=0.1*min(y,na.rm=TRUE)+0.9*max(y,na.rm=TRUE), label=sprintf("%s = %.3f\nP = %.3g", symb, r, p)))
    }
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
  dat <- dat[order(get(x))]
  dat[, c(y):=factor(get(y), levels=get(y))]
  if (is.null(xlab)) xlab <- x
  p <- ggplot(dat, aes(x=get(x), y=get(y))) +
    xlab(xlab) +
    theme_classic() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_text(size=10),
          axis.title.x=element_text(size=12))
  if (!is.null(color) && is.null(size)) p <- p + geom_point(aes(color=get(color))) + scale_color_continuous(low="red3", high="grey", name=color, guide=guide_colorbar(reverse=TRUE))
  if (!is.null(size) && is.null(color)) p <- p + geom_point(aes(size=get(size))) + scale_size_continuous(name=size)
  if (!is.null(color) && !is.null(size)) p <- p + geom_point(aes(color=get(color), size=get(size))) + scale_color_continuous(low="red3", high="grey", name=color, guide=guide_colorbar(reverse=TRUE)) + scale_size_continuous(name=size)

  return(p)
}

