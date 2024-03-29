## ----functions for quick data exploration and visualization----


my.cols <- function(x, dup.last=FALSE, na.rm=TRUE, na.color=NULL, no.grey=FALSE) {
  # my custom list of easily distinguishable colors for categorical variables with potentially many levels
  # note: order not optimized yet
  # x: if missing, will return all available colors; or a single number n, return n colors from top of the list;
  #    or a vector, if the vector contains no duplicates, return as many colors from top of the list and named by the vector elements in order;
  #    if the vector contains duplicates, will give colors for sort(table(x), dec=T)
  # dup.last: if not enough colors, whether to duplicate the last color or throw an error
  # na.rm: if x is a vector, whether to exclude NA (if any); when including NA, it will be placed last
  # na.color: if NULL and na.rm=FALSE, will use a color in the list for NA; if not NULL, then specify the color for NA (and will automatically assume na.rm=FALSE)
  # no.grey: exclude grey-ish colors from the list to choose from

  if (no.grey) {
    # this is to avoid confusion with NA, which is by default colored grey
    cols <- c("#E31A1C", "#295ED4", "#008B00", "#6A3D9A", "#FFFF00", "#FFACFD", "#00FFFF", "#8B4500", "#FFE4C4", "#00FF7F", "#FF1493", "#FFD700", "#FF7F00", "#66CD00", "#FF7D7D", "#AB82FF", "#D2B48C", "#ADFF2F", "#CD853F", "#00008B", "#B03060", "#9400D3", "#8B8B00", "#528B8B", "#7EC0EE", "#FF4500", "#CD96CD")
  } else {
    cols <- c("#E31A1C", "#295ED4", "#008B00", "#6A3D9A", "#FFFF00", "#FFACFD", "#00FFFF", "#8B4500", "#FFE4C4", "#00FF7F", "#FF1493", "#FFD700", "#FF7F00", "#DADADA", "#66CD00", "#FF7D7D", "#D2B48C", "#AB82FF", "#667788", "#ADFF2F", "#CD853F", "#333333", "#00008B", "#B03060", "#9400D3", "#8B8B00", "#528B8B", "#7EC0EE", "#FF4500", "#CD96CD")
  }

  if (missing(x)) return(cols)

  get.n.cols <- function(n) {
    if (n>length(cols)) {
      if (dup.last) {
        res <- c(cols, rep(cols[length(cols)], n-length(cols)))
      } else stop(sprintf("%d colors are available while you asked for %d", length(cols), n))
    } else res <- cols[1:n]
    res
  }
  if (length(x)==1 && is.numeric(x) && is.wholenumber(x)) {
    res <- get.n.cols(x)
  } else {
    any.na <- anyNA(x)
    if (any(duplicated(x))) x <- names(sort(table(x), decreasing=TRUE))
    if (!na.rm && is.null(na.color) && any.na) x <- c(x, NA)
    res <- cn1(get.n.cols(length(x)), x)
    if (!is.null(na.color) && any.na) res <- c(res, cn1(na.color, NA))
  }
  res
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


plot.xy <- function(x, y, xlab=NULL, ylab=NULL, dat=NULL, color=NULL, shape=NULL, size=NULL, label=NULL, label.subset=NULL, label.outliers=FALSE, outliers.cutoff=0.01, label.size=3, alpha=0.8, density="auto", diag=FALSE, trend="lm", cor.method=c("pearson","spearman","kendall"), cor.posi="auto", cor.size=4, do.plot=TRUE) {
  # 2-D scatter plot
  # x, y: numeric vectors of the same size, or column names in dat if dat is not NULL and provided as a data.frame or data.table or matrix
  # color, shape, size: vectors in the same order as x/y/rows of dat for plotting (can be a named list to customize legend title), or column names in dat if dat is not NULL
  # label: will label the points if this is anything but NULL; this can be a vector of labels in the same order as x/y/rows of dat, or column names in dat if dat is not NULL, or simply TRUE (in which case will use row names of dat, or row indices if there are no row names)
  # label.subset: specify the subset to add labels; a logical or numeric index vector in the same order as label/rows of dat, or a character vector of a subset of labels; will label all points if NULL
  # label.outliers: if TRUE, will label ourliers in red, independent from other label settings (even if label is NULL); to label only the outliers (and no other points), need to set label.subset to NA
  # ourliers.cutoff: the alpha level for determining outliers with a Chi-Squared distribution of the Mahalanobis distances
  # cor.posi, cor.size: for the text label on correlation
  # do.plot: whether to do plotting, if FALSE, will only return the ggplot object

  if (is.null(dat)) {
    if (is.null(xlab)) xlab <- deparse(substitute(x))
    if (is.null(ylab)) ylab <- deparse(substitute(y))
    dat <- data.table(x=x, y=y)
    x <- "x"
    y <- "y"
    ids <- 1:nrow(dat)
  } else {
    if (is.null(xlab)) xlab <- x
    if (is.null(ylab)) ylab <- y
    ids <- rownames(dat) # works for all types
    if (is.null(ids)) ids <- 1:nrow(dat)
    dat <- copy(as.data.table(dat))
  }
  if (!(length(color)==1 && color %in% names(dat) || is.null(color))) {
    dat[, color:=unlist(color)]
    if (!is.null(names(color))) {
      setnames(dat, "color", names(color))
      color <- names(color)
    } else color <- "color"
  }
  if (!(length(shape)==1 && shape %in% names(dat) || is.null(shape))) {
    dat[, shape:=unlist(shape)]
    if (!is.null(names(shape))) {
      setnames(dat, "shape", names(shape))
      shape <- names(shape)
    } else shape <- "shape"
  }
  if (!(length(size)==1 && size %in% names(dat) || is.null(size))) {
    dat[, size:=unlist(size)]
    if (!is.null(names(size))) {
      setnames(dat, "size", names(size))
      size <- names(size)
    } else size <- "size"
  }

  if (isTRUE(label)) label <- ids
  if (label.outliers && is.null(label)) {
    label <- ids
    if (is.null(label.subset)) label.subset <- NA
  }
  if (!(length(label)==1 && label %in% names(dat) || is.null(label))) {
    dat[, label:=label]
    label <- "label"
  }
  if (!is.null(label.subset)) {
    if (length(label.subset)==1 && label.subset %in% names(dat)) {
      label.subset <- which(dat[[label.subset]])
    } else if (is.character(label.subset)) {
      label.subset <- which(ids %in% label.subset)
    } else if (is.logical(label.subset)) {
      label.subset <- which(label.subset)
    }
  }
  if (label.outliers) {
    tmp <- data.matrix(dat[, c(x, y), with=FALSE])
    md <- mahalanobis(tmp, colMeans(tmp), cov(tmp))
    cutoff <- qchisq(p=1-outliers.cutoff, df=ncol(tmp))
    id.outliers <- which(md>cutoff)
    if (!is.null(label.subset)) label.subset <- union(label.subset, id.outliers)
  }

  p <- ggplot(dat, aes_string(x=x, y=y)) +
    xlab(xlab) + ylab(ylab) +
    geom_point(aes_string(color=color, shape=shape, size=size), alpha=alpha) +
    theme_classic()

  if (density=="auto") {
    dmat <- table(cut(dat[[x]], breaks=5), cut(dat[[y]], breaks=5))
    density <- sum(dmat>500)>5
  }
  if (density==TRUE) p <- p + stat_density_2d(aes(color=..level..), geom="polygon", alpha=0) + theme(legend.position="none")

  if ("loess" %in% trend) {
    set.seed(1)
    p <- p + geom_smooth(method=loess, color="red2", size=0.8, fill="red2", alpha=0.2)
  }
  if (diag) p <- p + geom_abline(slope=1, intercept=0, color="grey")
  if ("lm" %in% trend) {
    p <- p + geom_smooth(method=lm, color="blue", size=0.8, fill="blue", alpha=0.2)
    cor.method <- match.arg(cor.method)
    ct <- cor.test(dat[[x]], dat[[y]], method=cor.method)
    symb <- switch(cor.method, pearson="r", spearman="rho", kendall="tau")
    pval <- ct$p.value
    r <- ct$estimate
    if (pval>=2.2e-16) lab <- sprintf("%s=%.3f\nP=%.3g", symb, r, pval) else lab <- sprintf("%s=%.3f\nP<2.2e-16", symb, r)
    if (length(cor.posi)==1 && cor.posi=="auto") {
      if (!exists("dmat")) dmat <- table(cut(dat[[x]], breaks=5), cut(dat[[y]], breaks=5))
      tmp <- which(dmat==min(dmat))
      if (r>0) tmp1 <- match(c(5,10,4,15,9,3,21,22,16,23,17,11,20,14,8,2,24,18,12,6,25,19,13,7,1), tmp)
      else tmp1 <- match(c(25,20,24,15,19,23,1,2,6,3,7,11,10,14,18,22,4,8,12,16,5,9,13,17,21), tmp)
      tmp <- tmp[tmp1[!is.na(tmp1)][1]]
      i <- tmp %% 5
      if (i==0) i <- 5
      j <- (tmp-1) %/% 5 + 1
      cor.posi <- c((1.1-0.2*i)*min(dat[[x]],na.rm=TRUE)+(0.2*i-0.1)*max(dat[[x]],na.rm=TRUE), (1.1-0.2*j)*min(dat[[y]],na.rm=TRUE)+(0.2*j-0.1)*max(dat[[y]],na.rm=TRUE))
    }
    p <- p + annotate("text", x=cor.posi[1], y=cor.posi[2], label=lab, size=cor.size)
  }

  if (!is.null(label)) {
    if (is.null(label.subset)) dat[, lab.flag:=TRUE] else dat[label.subset, lab.flag:=TRUE]
    dat[, lab.color:="black"]
    if (label.outliers && length(id.outliers)>0) dat[id.outliers, c("outlier", "lab.color"):=list(TRUE, "red2")]
    p <- p + geom_text_repel(data=dat[lab.flag==TRUE], aes_string(label=label), color=dat[lab.flag==TRUE, lab.color], size=label.size)
  }

  if (do.plot) print(p)
  invisible(list(p=p, plot.data=dat))
}


plot.pca <- function(mat, pc.x=1, pc.y=2, max.pc=50, data=NULL, color=NULL, shape=NULL, size=NULL, label=NULL, label.subset=NULL, label.outliers=TRUE, outliers.cutoff=0.01, label.size=3, alpha=0.8, ld.color=NULL, ld.rev=NULL, do.plot="pc", center=TRUE, scale=TRUE, ...) {
  # PCA plot, given a matrix mat (or data.frame, or data.table where the first column is sample name/ID) of sample-by-variable
  # or mat can be the prcomp object, then center, scale and ... will be ignored
  # pc.x and pc.y: PC's to plot on the x any y axes
  # max.pc: max number of PCs to plot on the scree plot
  # data: data.table with the same number of rows as mat for additional variables of the samples
  # color, shape, size, label*, outlier*, alpha: for PC plot with plot.xy
  # ld.color: vector corresponding to the variables (columns of mat) for plotting
  # ld.rev: logical vector corresponding to the variables, whether to plot its loading in the reverse direction; or provide names of the variables to be reversed
  # do.plot: what to plot, a vector of one or more of "scree", "pc", "loading"; or "all" as a shorthand for all plots; set to NULL to disable plotting
  # center, scale, ...: passed to prcomp()

  if (!is.null(do.plot) && any(!do.plot %in% c("scree", "pc", "loading", "all"))) stop('`do.plot` should be among "scree", "pc", "loading", and "all".')
  if (length(do.plot)==1 && do.plot=="all") do.plot <- c("scree", "pc", "loading")
  if ("prcomp" %in% class(mat)) {
    res <- mat
  } else {
    if (length(class(mat))==1 && class(mat)=="data.frame") mat <- data.matrix(mat) else if ("data.table" %in% class(mat)) mat <- dt2mat(mat)
    if (scale) {
      id <- apply(mat, 2, uniqueN)==1
      s <- sum(id)
      if (s!=0) {
        message(sprintf("removed %d columns of zero variance.", s))
        mat <- mat[, !id]
      }
    }
    ids <- rownames(mat)
    res <- prcomp(mat, center=center, scale.=scale, ...)
  }

  tot.var <- sum(res$sdev^2)
  varx <- sprintf("PC %d (%.1f%%)", pc.x, res$sdev[pc.x]^2 /tot.var*100)
  vary <- sprintf("PC %d (%.1f%%)", pc.y, res$sdev[pc.y]^2 /tot.var*100)

  var.dat <- data.table(PC=1:length(res$sdev), `Each PC`=res$sdev^2/tot.var*100)
  elbow <- get.elbow(var.dat)
  var.dat <- var.dat[1:min(max.pc, .N)]
  var.dat[, Cumulative:=cumsum(`Each PC`)]
  var.dat <- melt(var.dat, id.vars="PC", variable.name="what", value.name="y")
  pvar <- ggplot(var.dat, aes(x=PC, y=y, color=what)) +
    scale_x_continuous("PC", breaks=1:length(res$sdev)) + scale_y_continuous("% Total Variance", n.breaks=10) +
    geom_vline(xintercept=elbow, color="grey", linetype="dashed") +
    geom_point() +
    geom_text_repel(data=var.dat[!(PC==1 & what=="Cumulative")], aes(x=PC, y=y, label=sprintf("%.1f",y)), color="grey10", size=label.size) +
    geom_line(aes(group=what)) +
    scale_color_manual(values=c("grey10","darkblue")) +
    theme_classic() +
    theme(panel.grid.major.y=element_line(size=0.4),
          legend.position="none")

  dat <- cbind(x=res$x[, pc.x], y=res$x[, pc.y])
  if (!is.null(data)) {
    if (nrow(data)!=nrow(dat)) stop("mat and data have different numbers of rows!")
    dat <- cbind(dat, as.data.frame(data))
  }
  p <- plot.xy(x="x", y="y", xlab=varx, ylab=vary, dat=dat, color=color, shape=shape, size=size, label=label, label.subset=label.subset, label.outliers=label.outliers, outliers.cutoff=outliers.cutoff, label.size=label.size, alpha=alpha, density=FALSE, diag=FALSE, trend=NULL, do.plot=FALSE)

  circ <- data.table(x=0.5*cos(seq(0,2*pi, len=500)), y=0.5*sin(seq(0,2*pi, len=500)))
  loadings <- data.table(lab=rownames(res$rotation), x=res$rotation[, pc.x], y=res$rotation[, pc.y])
  if (!is.null(ld.color)) {
    if ("prcomp" %in% class(mat) && length(ld.color)!=nrow(loadings)) warning("length of `ld.color` does not match the number of loadings, will skip coloring of loadings; please double check `nrow(mat$rotation)`.")
    loadings[, color:=ld.color[colnames(mat) %in% loadings$lab]]
    ld.color <- "color"
  }
  if (!is.null(ld.rev)) {
    if (is.logical(ld.rev)) {
      if ("prcomp" %in% class(mat) && length(ld.rev)!=nrow(loadings)) warning("length of `ld.rev` does not match the number of loadings, will skip reversed plotting of loadings; please double check `nrow(mat$rotation)`.")
      ld.rev <- ld.rev[colnames(mat) %in% loadings$lab]
    } else if (is.character(ld.rev)) ld.rev <- loadings$lab %in% ld.rev
  } else ld.rev <- FALSE
  loadings[, rev:=ld.rev]
  loadings[rev==TRUE, c("x", "y"):=list(-x, -y)]
  pld <- ggplot()+ xlab(varx) + ylab(vary) +
    geom_path(data=circ, aes(x=x, y=y), color="grey", linetype="dashed") +
    geom_hline(yintercept=0, color="grey", linetype="dashed") +
    geom_vline(xintercept=0, color="grey", linetype="dashed")
  if (loadings[, any(!rev)]) pld <- pld + geom_segment(data=loadings[rev==FALSE], aes_string(x=0, y=0, xend="x", yend="y", color=ld.color), arrow=arrow(length=unit(5,"pt"), type="closed"), lwd=0.7, alpha=0.7)
  if (loadings[, any(rev)]) pld <- pld + geom_segment(data=loadings[rev==TRUE], aes_string(x=0, y=0, xend="x", yend="y", color=ld.color), arrow=arrow(angle=150, length=unit(6,"pt")), lwd=0.7, alpha=0.7)
  pld <- pld + geom_text_repel(data=loadings, aes(x=x, y=y, label=lab), size=label.size) +
    coord_equal() +
    theme_classic()

  if ("scree" %in% do.plot) {
    print(pvar)
  }
  if ("pc" %in% do.plot) {
    print(p$p)
  }
  if ("loading" %in% do.plot) {
    print(pld)
  }

  invisible(list(pca=res, scree.plot.data=var.dat, scree.plot=pvar, elbow=elbow, pc.plot.data=p$plot.data, pc.plot=p$p, loading.plot.data=loadings, loading.plot=pld, outliers=if ("outlier" %in% names(p$plot.data)) p$plot.data[outlier==TRUE, label] else NULL))
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
    #tmp <- do.call(wilcox, c(list(value~group, dat), more.args))
    tmp <- do.call(wilcox, c(list(s1=l[[1]], s2=l[[2]]), more.args))
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
