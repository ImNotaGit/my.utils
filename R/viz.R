## ----functions for quick data exploration and visualization----


my.cols <- function(x, dup.last=FALSE, na.rm=TRUE, na.color=NULL, no.grey=FALSE, mute=FALSE, l=50, c=100) {
  # my custom list of easily distinguishable colors for categorical variables with potentially many levels
  # note: order not optimized yet
  # x: if missing, will return all available colors; or a single number n, return n colors from top of the list;
  #    or a factor, return as many colors from top of the list and named by the levels of the factor in order;
  #    or a vector, if the vector contains no duplicates, return as many colors from top of the list and named by the vector elements in order;
  #    if the vector contains duplicates, will give colors for sort(table(x), dec=T)
  # dup.last: if not enough colors, whether to duplicate the last color or throw an error
  # na.rm: if x is a vector, whether to exclude NA (if any); when including NA, it will be placed last
  # na.color: if NULL and na.rm=FALSE, will use a color in the list for NA; if not NULL, then specify the color for NA (and will automatically assume na.rm=FALSE)
  # no.grey: exclude grey-ish colors from the list to choose from
  # mute: if true will "dampen" the colors with scales::muted, l and c are passed to muted

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
    if (is.factor(x)) x <- levels(x) else if (any(duplicated(x))) x <- names(sort(table(x), decreasing=TRUE)) else x <- x[!is.na(x)]
    if (!na.rm && is.null(na.color) && any.na) x <- c(x, NA)
    res <- cn1(get.n.cols(length(x)), x)
    if (!is.null(na.color) && any.na) res <- c(res, cn1(na.color, NA))
  }
  if (mute) res <- scales::muted(res, l=l, c=c)
  res
}


subchunkify <- function(p, w, h, nm=NULL, ...) {
  # function to create sub-chunks within normal chunks for plots with modified sizes in Rmd documents
  # note: other non-subchunkified plots in the parent chunk may not be produced at the correct position (i.e. messed-up order), so if using subchunkify, then all plots within the parent chunk should be subchunkified;
  # additional chunk options like `message=FALSE`, etc. can be provided in ..., however it seems that in order for this to work, the parent chunk should also have the same options set; subchunks w/o these options can override the parent chunk options
  p.deparsed <- paste0(deparse(function() {p}), collapse="")
  more.args <- deparse(c(...))
  if (more.args=="NULL") more.args <- "" else more.args <- stringr::str_sub(more.args, 3, -2)
  if (is.null(nm)) {
    set.seed(NULL)
    nm <- sprintf("subchunk_%s", floor(runif(1)*1e10))
  } else nm <- paste0("subchunk_", nm)
  sub.chunk <- sprintf("\n```{r %s, fig.width=%s, fig.height=%s, echo=FALSE, %s}\n(%s)()\n```\n", nm, w, h, more.args, p.deparsed)
  cat(trimws(knitr::knit(text=knitr::knit_expand(text=sub.chunk), quiet=TRUE)))
}


plot.xy <- function(x, y, dat=NULL, xlab=NULL, ylab=NULL, color=NULL, shape=NULL, size=NULL, label=NULL, label.subset=NULL, label.outliers=FALSE, outliers.cutoff=0.01, label.size=3, alpha=0.8, density="auto", diag=FALSE, trend="lm", cor.method=c("pearson","spearman","kendall"), cor.pos="auto", cor.size=4, do.plot=TRUE) {
  # 2-D scatter plot
  # x, y: numeric vectors of the same size, or column names in dat if dat is not NULL and provided as a data.frame or data.table or matrix
  # color, shape, size: vectors in the same order as x/y/rows of dat for plotting (can be a named list to customize legend title), or column names in dat if dat is not NULL
  # label: will label the points if this is anything but NULL; this can be a vector of labels in the same order as x/y/rows of dat, or column names in dat if dat is not NULL, or simply TRUE (in which case will use row names of dat, or row indices if there are no row names)
  # label.subset: specify the subset to add labels; a logical or numeric index vector in the same order as label/rows of dat, or a character vector of a subset of labels; will label all points if NULL
  # label.outliers: if TRUE, will label ourliers in red, independent from other label settings (even if label is NULL); to label only the outliers (and no other points), need to set label.subset to NA
  # ourliers.cutoff: the alpha level for determining outliers with a Chi-Squared distribution of the Mahalanobis distances
  # cor.pos, cor.size: for the text label on correlation, its position and size
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
    ct <- cor.test(as.numeric(dat[[x]]), as.numeric(dat[[y]]), method=cor.method)
    symb <- switch(cor.method, pearson="r", spearman="rho", kendall="tau")
    pval <- ct$p.value
    r <- ct$estimate
    if (is.na(pval)) lab <- sprintf("%s=%.3f\nP=NA", symb, r) else if (pval>=2.2e-16) lab <- sprintf("%s=%.3f\nP=%.3g", symb, r, pval) else lab <- sprintf("%s=%.3f\nP<2.2e-16", symb, r)
    if (length(cor.pos)==1 && cor.pos=="auto") {
      if (!exists("dmat")) dmat <- table(cut(dat[[x]], breaks=5), cut(dat[[y]], breaks=5))
      tmp <- which(dmat==min(dmat))
      if (is.na(r) || r>0) tmp1 <- match(c(5,10,4,15,9,3,21,22,16,23,17,11,20,14,8,2,24,18,12,6,25,19,13,7,1), tmp)
      else tmp1 <- match(c(25,20,24,15,19,23,1,2,6,3,7,11,10,14,18,22,4,8,12,16,5,9,13,17,21), tmp)
      tmp <- tmp[tmp1[!is.na(tmp1)][1]]
      i <- tmp %% 5
      if (i==0) i <- 5
      j <- (tmp-1) %/% 5 + 1
      cor.pos.x <- (1.1-0.2*i)*min(as.numeric(dat[[x]]),na.rm=TRUE)+(0.2*i-0.1)*max(as.numeric(dat[[x]]),na.rm=TRUE)
      cor.pos.y <- (1.1-0.2*j)*min(as.numeric(dat[[y]]),na.rm=TRUE)+(0.2*j-0.1)*max(as.numeric(dat[[y]]),na.rm=TRUE)
      if (lubridate::is.Date(dat[[x]])) cor.pos.x <- as.Date(cor.pos.x, origin="1970-01-01")
      if (lubridate::is.Date(dat[[y]])) cor.pos.y <- as.Date(cor.pos.y, origin="1970-01-01")
    } else {
      cor.pos.x <- cor.pos[1]
      cor.pos.y <- cor.pos[2]
    }
    p <- p + annotate("text", x=cor.pos.x, y=cor.pos.y, label=lab, size=cor.size)
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

  dat <- data.frame(x=res$x[, pc.x], y=res$x[, pc.y])
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

  if (isTRUE(label) || (label.outliers && is.null(label)) || !(length(label)==1 && label %in% names(dat))) label <- "label"
  invisible(list(pca=res, scree.plot.data=var.dat, scree.plot=pvar, elbow=elbow, pc.plot.data=p$plot.data, pc.plot=p$p, loading.plot.data=loadings, loading.plot=pld, outliers=if ("outlier" %in% names(p$plot.data)) p$plot.data[outlier==TRUE][[label]] else NULL)) # for outliers need to use [[label]] instead of get(label)
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


xa <- function(x, a) (1-a)*min(x[is.finite(x)]) + a*max(x[is.finite(x)]) # is.finite will ignore Inf, -Inf, NA and NaN


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

  # if there are only 2 or 3 groups, do wilcoxon test for each pair of groups
  ll <- length(l)
  if (ll==2) {
    #tmp <- do.call(wilcox, c(list(value~group, dat), more.args))
    tmp <- do.call(wilcox, c(list(s1=l[[1]], s2=l[[2]]), more.args))
    stat <- dat[, .(id=12, x1=1, x2=2, x=1.5, y=xa(value,1.1), p=tmp$pval, r=tmp$r.wilcox)]
  } else if (ll==3) {
    tmp <- do.call(wilcox3, c(list(value~group, dat), more.args))
    stat <- data.table(id=numeric(0), x=numeric(0), y=numeric(0), p=numeric(0), r=numeric(0))
    if (12 %in% plab) stat <- rbind(stat, dat[group %in% levels(group)[1:2], .(id=12, x=1.5, y=xa(value,1.1), p=tmp$pval[1], r=tmp$r.wilcox[1])])
    if (23 %in% plab) stat <- rbind(stat, dat[group %in% levels(group)[2:3], .(id=23, x=2.5, y=xa(value,1.1), p=tmp$pval[2], r=tmp$r.wilcox[2])])
    if (13 %in% plab) {
      if (length(stat$y)==0) a <- 1.1
      else if (max(stat$y)>xa(dat$value,0.9)) a <- 1.22
      else if (max(stat$y)>xa(dat$value,0.8)) a <- 1.17
      else a <- 1.1
      stat <- rbind(stat, data.table(id=13, x=2, y=xa(c(stat$y,dat$value),a), p=tmp$pval[3], r=tmp$r.wilcox[3]))
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
    p <- p + geom_blank(data=dat[, .(group=group[1], y=xa(c(stat$y,value),1.05))], aes(y=y))
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
    stat.out <- data.table(x=mean(1:ll), y=sapply(datl, function(x) x[, xa(value, 1.12)]), strata=factor(names(stat.out), levels=names(stat.out)), s=stat.out)
  }

  # blank data used to adjust axis limits
  blk <- dat[, .(ymax=xa(value, 1.2), group=group), by=strata]

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


.plot.fracs <- function(dat, xlab=NULL, ylab="Fraction", tit=NULL, facet.tit=NULL, lab=NULL, lab.size=3, lab.frac.ge=0.05, clrs=NULL, mdat=NULL, xclrs=NULL, xord="clust", clust.fun=NULL, dendro=TRUE, dend.scale=1, rs.dend=0.12, ori=c("h", "v"), no.axs=FALSE, no.nlab=FALSE, lgd.pos="bottom", lgd.tit=NULL, lgd.mgn=NULL, rs.lgd=0.1, ret.axs=FALSE, ...) {
  # inner plotting function of plot.fracs for plotting proportion stacked barplots with geom_bar(..., position="fill")
  # facet.tit: if not NULL, will place everything under a single facet with this label
  # no.axs: do not plot X or Y axis depending on plot direction
  # no.nlab: do not add the "Total #" label
  # the above options are used for plotting each group in grouped bar plots (the groups will be combined with cowplot::plot_grid downstream)
  # ret.axs: return only the Y axis; for dendro=TRUE this did not work well and was not used

  dir <- match.arg(ori)

  if (ret.axs) no.axs <- FALSE
  if (length(xord)==1 && xord=="clust" && uniqueN(dat$x)==1) {
    xord <- "keep"
    dendro <- FALSE # there will be an issue if there are multiple sample groups (xgrps in plot.fracs) and the other groups have >1 samples and dendro=TRUE; for now I ignore such cases
  }
  if (length(xord)==1 && xord=="clust") {
    mat <- dt2mat(dcast(dat, x~ygrp, value.var="y"))
    mat[is.na(mat)] <- 0
    if (is.null(clust.fun)) {
      hc <- hclust(dist(mat))
      xlvls <- rownames(mat)[hc$order]
    } else {
      hc <- clust.fun(mat)
      if (class(hc)=="hclust") {
        xlvls <- rownames(mat)[hc$order]
      } else if (is.vector(hc)) {
        xlvls <- hc
      } else stop("`clust.fun` applied to the data matrix did not return an hclust object or a vector of reordered sample names.")
    }
    if (dendro) {
      tmp <- ggdendro::dendro_data(as.dendrogram(hc), type="rectangle")
      dend <- ggplot(ggdendro::segment(tmp)) +
        xlab(xlab) +
        geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
        ggdendro::theme_dendro()
      if (dir=="h") {
        dend <- dend + coord_flip(clip="off") + scale_x_reverse() + scale_y_reverse()
        if (!is.null(xlab)) dend <- dend + theme(axis.title.y=element_text())
        if (ret.axs) {
          dend <- dend + theme(
            axis.line.x=element_line(color="white"),
            axis.ticks.x=element_line(color="white"),
            axis.title.x=element_text(color="white"),
            axis.text.x=element_text(color="white")
          )
        }
      } else if (dir=="v") {
        dend <- dend + coord_cartesian(clip="off") + scale_y_reverse()
        if (!is.null(xlab)) dend <- dend + theme(axis.title.x=element_text(color="black")) # needs to specify color to override theme_dendro(), a bug
        if (ret.axs) {
          dend <- dend + theme(
            axis.line.y=element_line(color="white"),
            axis.ticks.y=element_line(color="white"),
            axis.title.y=element_text(color="white"),
            axis.text.y=element_text(color="white")
          )
        }
      }
    }
  } else if (length(xord)==1 && xord=="keep") {
    if (is.factor(dat$x)) xlvls <- levels(dat$x) else xlvls <- unique(dat$x)
  } else if (length(xord)==1 && xord=="default") {
    xlvls <- levels(factor(dat$x))
  } else if (setequal(xord, as.character(dat$x))) {
    xlvls <- xord
  } else stop("Invalid `xord`.")
  xlvls <- xlvls[xlvls %in% dat$x]
  if (dir=="h") xlvls <- rev(xlvls)

  dat1 <- copy(dat)
  dat1[, x:=factor(x, levels=xlvls)]
  if (is.factor(dat1$ygrp)) ylvls <- rev(levels(dat1$ygrp)) else ylvls <- rev(levels(factor(dat1$ygrp)))
  dat1[, ygrp:=factor(ygrp, levels=ylvls)]
  if (!is.null(facet.tit)) dat1[, fct:=facet.tit]
  if (!is.null(mdat)) mdat <- mdat[match(xlvls, x)]
  if (!is.null(lab)) {
    if (isTRUE(lab)) dat1[, label:=ygrp] else dat1[, label:=get(lab)]
    dat1[y<lab.frac.ge, label:=""]
  }
  p <- ggplot(dat1, aes(x=x, y=y, fill=ygrp)) +
    ggtitle(tit) + xlab(xlab) + scale_y_continuous(name=ylab, breaks=seq(0, 1, 0.1)) +
    { if (is.null(facet.tit)) NULL else facet_wrap(~fct, strip.position=switch(dir, h="right", v="top")) } +
    geom_bar(stat="identity", position="fill", width=0.85, color="grey50", size=0.2) +
    { if ("sep" %in% names(mdat)) geom_vline(xintercept=mdat[, which(sep!=sep[c(2:.N,.N)])+0.5], size=0.2, color="grey50") else NULL } +
    { if (is.null(clrs)) scale_fill_discrete(name=lgd.tit, guide=guide_legend(title.position=if (lgd.pos=="bottom") "bottom" else "top", title.hjust=if (lgd.pos=="bottom") 0.5 else 0, ...)) else scale_fill_manual(name=lgd.tit, values=clrs, guide=guide_legend(title.position=if (lgd.pos=="bottom") "bottom" else "top", title.hjust=if (lgd.pos=="bottom") 0.5 else 0, ...)) } +
    { if (is.null(lab)) NULL else geom_text(aes(label=label), size=lab.size, position=position_fill(vjust=0.5)) } +
    { if ("ntot" %in% names(mdat) && !no.nlab) annotate(geom="text", x=length(xlvls)+1, y=1.05, label="Total #", size=3, angle=switch(dir, h=0, v=90)) else NULL } +
    { if ("ntot" %in% names(mdat)) annotate(geom="text", x=1:length(xlvls), y=1.05, label=mdat$ntot, size=3, angle=switch(dir, h=0, v=90)) else NULL } +
    theme_classic()

  if (dir=="h") {
    p <- p + coord_flip(clip="off") + theme(
      plot.title=if (is.null(tit)) element_blank() else element_text(),
      axis.text.x=if (no.axs) element_blank() else element_text(hjust=1, angle=35),
      axis.title.x=if (no.axs) element_blank() else element_text(),
      axis.line.x=if (no.axs) element_blank() else element_line(),
      axis.ticks.x=if (no.axs) element_blank() else element_line(),
      axis.text.y=if (is.null(xclrs)) element_text() else element_text(color=xclrs[xlvls]),
      axis.title.y=if (is.null(xlab) || dendro) element_blank() else element_text(),
      panel.grid.major.x=element_line(size=0.2, color="grey50"),
      strip.text=element_text(size=9, margin=margin(0,2,0,2)),
      strip.background=element_rect(size=0.6),
      plot.margin=if (no.nlab) margin(0,0,0,5*max(nchar(as.character(xlvls))-(10+10*(0:(length(xlvls)-1))),0)) else margin(10,0,0,5*max(nchar(as.character(xlvls))-(10+10*(0:(length(xlvls)-1))),0)),
      legend.position=lgd.pos,
      legend.title=if (is.null(lgd.tit)) element_blank() else element_text(size=8.5),
      legend.text=element_text(size=8),
      legend.key.size=unit(8,"pt"),
      legend.box.margin=if (lgd.pos=="bottom") margin(c(-10,10,0,0)) else margin()
    )
  } else if (dir=="v") {
    if (dendro) {
      ag <- 90
      hj <- 1
      vj <- 0.5
    } else {
      ag <- 40
      hj <- 1
      vj <- NULL
    }
    p <- p + coord_cartesian(clip="off") + theme(
      plot.title=if (is.null(tit)) element_blank() else element_text(),
      axis.text.y=if (no.axs) element_blank() else element_text(),
      axis.title.y=if (no.axs) element_blank() else element_text(),
      axis.line.y=if (no.axs) element_blank() else element_line(),
      axis.ticks.y=if (no.axs) element_blank() else element_line(),
      axis.text.x=if (is.null(xclrs)) element_text(angle=ag, hjust=hj, vjust=vj) else element_text(angle=ag, hjust=hj, vjust=vj, color=xclrs[xlvls]),
      axis.title.x=if (is.null(xlab) || dendro) element_blank() else element_text(),
      panel.grid.major.y=element_line(size=0.2, color="grey50"),
      strip.text=element_text(size=9, margin=margin(2,0,2,0)),
      strip.background=element_rect(size=0.6),
      plot.margin=if (no.nlab) margin(0,0,0,5*max(nchar(as.character(xlvls))-(10+10*(0:(length(xlvls)-1))),0)) else margin(0,10,0,5*max(nchar(as.character(xlvls))-(10+10*(0:(length(xlvls)-1))),0)),
      legend.position=lgd.pos,
      legend.title=if (is.null(lgd.tit)) element_blank() else element_text(size=8.5),
      legend.text=element_text(size=8),
      legend.key.size=unit(8,"pt"),
      legend.box.margin=if (!is.null(lgd.mgn)) margin(lgd.mgn) else if (lgd.pos=="bottom") margin(c(-10,10,0,0)) else margin()
    )
  }

  if (ret.axs) {
    if (dir=="h") {
      axs <- gtable::gtable_filter(ggplotGrob(p), 'axis-b|xlab-b', trim=FALSE)
      axs <- gtable::gtable_add_padding(axs, unit(c(0,0,6,0), "pt"))
    } else if (dir=="v") {
      axs <- gtable::gtable_filter(ggplotGrob(p), 'axis-l|ylab-l', trim=FALSE)
      axs <- gtable::gtable_add_padding(axs, unit(c(0,0,0,6), "pt"))
    }
    if (dendro) {
      if (dir=="h") {
        daxs <- gtable::gtable_filter(ggplotGrob(dend), 'axis-b|xlab-b', trim=FALSE)
        daxs <- gtable::gtable_add_padding(daxs, unit(c(0,0,6,0), "pt"))
        axs <- arrangeGrob(axs, left=daxs, widths=c(1, 0.15))
      } else if (dir=="v") {
        daxs <- gtable::gtable_filter(ggplotGrob(dend), 'axis-l|ylab-l', trim=FALSE)
        daxs <- gtable::gtable_add_padding(daxs, unit(c(0,1,0,6), "pt"))
        axs <- arrangeGrob(axs, bottom=daxs, heights=c(1, 0.1))
      }
    }
    return(axs)
  }

  if (dendro) {
    if (dir=="h") {
      cowplot::plot_grid(dend, p, nrow=1, align="h", axis="tb", rel_widths=c(rs.dend, 1), scale=c(1, dend.scale))
    } else if (dir=="v") {
      if (lgd.pos=="bottom") {
        lgd <- cowplot::get_legend(p)
        p <- p + theme(legend.position="none")
        cowplot::plot_grid(p, dend, lgd, ncol=1, align="v", axis="lr", rel_heights=c(1, rs.dend, rs.lgd), scale=c(1, dend.scale, 1))
      } else {
        cowplot::plot_grid(p, dend, ncol=1, align="v", axis="lr", rel_heights=c(1, rs.dend), scale=c(1, dend.scale))
      }
    }
  } else p
}

plot.fracs <- function(dat, mode=c("count", "frac"), xlab, ylab="Fraction", tit=NULL, xvar=NULL, ygrp=NULL, yvar=NULL, lab=NULL, lab.size=3, lab.frac.ge=0.05, mdat=NULL, mdat.xvar=NULL, ntot=NULL, xgrp=NULL, xcol=NULL, xsep=NULL, xord="clust", clust.fun=NULL, dendro=TRUE, dend.scale=1, rs.dend=0.12, yord=c("default", "keep"), lgd.tit=NULL, lgd.pos=c("bottom", "right", "none"), rs.xgrp=NULL, rs.lgd=NULL, ori=c("h", "v"), ...) {
  # function for visualizing fraction data with proportion stacked bar plots (i.e. with geom_bar(..., position="fill"))
  # dat: a variable-by-sample matrix or a data.table in the long format
  # mode: whether the data values are counts or fractions
  # xvar, ygrp, yvar: if dat is a long data.table, these are required and are the column names in dat corresponding to the sample, variable, and variable value (i.e. count or fraction values) respectively
  # lab: if not NULL, labels will be placed on each bar when fraction>=lab.frac.ge: if TRUE, the labels will be the variable names; if dat is a long data.table, can also be a column name and the corresponding values in this column will be used as the labels
  # mdat: optional, a data.table of sample-level meta data, one row per sample
  # mdat.xvar: the column name in mdat for sample (if NULL will use xvar, and if xvar is NULL will assume the first column)
  # ntot: if dat is fraction, can provide total count per sample data here; if dat is count, will calculate total count from dat and this will be ignored; provided either as a vector (preferentially named by sample name) or a column name in mdat (if the latter is provided)
  # xgrp: sample group, if not NULL, will plot bar plots for each group and then combine all plots with cowplot::plot_grid; provided similarly as ntot; this can be a factor and its levels determine the order of the plots
  # xcol: variable used to color the sample names (i.e. axis texts); provided similarly as xgrp; if a factor the levels will affect color mapping
  # xsep: variable used to group the sample within each plot by separating lines; provided similarly as ntot
  # xord: sample order in each plot, "clust" or "default" or "keep" or a vector of custom order, "clust" for clustering (by default hclust on dist on default parameters, can also provide clust.fun, a function that takes the data matrix and returns the ordered sample names after clustering), "default" for default converting-to-factor behavior, "keep" for keeping the original order in dat
  # dendro: if TRUE and xord="clust", show dendrogram; will be ignored if xord!="clust"
  # dend.scale: scaling of dendrogram, used to manually ensure alignment between the dendrogram and main plot (unfortunately I did not find a general way of perfect alignment automatically); if xgrp is not NULL, may provide a vector in the same length and order as the sample groups
  # rs.dend: size of dendrogram relative to main plot, provide as a single number
  # yord: ordering of the variables, "default" for roughly decreasing order by the median fraction across sample, "keep" for keeping original order in dat
  # lgd.tit, lgd.pos: legend title and position
  # ori: plot orientation (horizontal or vertical; if "h", the bars will be horizontal, i.e. each sample is a row, Y-axis is sample; vice versa for "v"); note: initially I named this "dir" but there seems to be a conflict with ... passed to guide_legend()
  # rs.xgrp: for adjusting relative sizes of the sample groups; if needed provide as a vector with the same length as the number of groups (xgrp) and in the corresponding order (in terms of the final order of xgrp); in terms of corresponding to the plots, the order is left to right or top to bottom; the values are adjustments to default, i.e. should usually be around 0 +- ~0.? to < ~5
  # rs.lgd: size of legend relative to main plot, provide as needed as a single number
  # ...: passed to scale_fill_?(guide=guide_legend(...))

  mode <- match.arg(mode)
  yord <- match.arg(yord)
  dir <- match.arg(ori)
  lgd.pos <- match.arg(lgd.pos)
  if (length(xord)>1 || xord!="clust") dendro <- FALSE

  if (is.matrix(dat)) {
    if ("table" %in% class(dat)) class(dat) <- "matrix"
    if (mode=="count") {
      ntot <- colSums(dat)
      dat <- dat/rep(colSums(dat), each=nrow(dat))
    }
    dat <- melt(cbind(ygrp=rownames(dat), as.data.table(dat)), id.vars="ygrp", variable.name="x", value.name="y")
    if (missing(xlab)) xlab <- NULL
  } else {
    dat <- dat[, .(x=get(xvar), ygrp=get(ygrp), y=get(yvar))]
    if (mode=="count") {
      ntot <- dat[, .(n=sum(y, na.rm=TRUE)), by=x][, setNames(n, x)]
      dat <- dat[, .(ygrp, y=y/sum(y, na.rm=TRUE)), by=x]
    }
    if (missing(xlab)) xlab <- xvar
  }

  tmp <- dat[, .(ymed=median(y, na.rm=TRUE), ymax=max(y, na.rm=TRUE)), by=ygrp][order(-ymed)]
  if (yord=="keep") {
    if (is.factor(dat$ygrp)) ylvls <- levels(dat$ygrp) else ylvls <- unique(dat$ygrp)
  } else if (yord=="default") ylvls <- as.character(tmp$ygrp)
  ncl <- length(my.cols())
  if (length(ylvls)<=ncl) {
    clrs <- my.cols(ylvls)
  } else {
    if (yord=="default") {
      tmp1 <- tmp[, which(ymax>0.05)]
      tmp2 <- tmp[, which(ymax<=0.05)]
      idxs <- 1:(ncl-1)
      n1 <- sum(tmp1>=ncl)
      n2 <- sum(tmp2<ncl)
      if (n1>0 && n2>0) {
        if (n1>=n2) idxs[idxs %in% tmp2[tmp2<ncl]] <- tmp1[tmp1>=ncl][1:n2]
          else idxs[idxs %in% rev(tmp2[tmp2<ncl])[1:n1]] <- tmp1[tmp1>=ncl]
      }
      clrs <- my.cols(tmp[idxs, ygrp])
      clrs <- c(clrs, setNames(rep(my.cols()[ncl], tmp[-idxs, .N]), tmp[-idxs, ygrp]))
      ylvls <- names(clrs)
    } else if (yord=="keep") {
      clrs <- my.cols(ylvls[1:ncl])
      clrs <- c(clrs, setNames(rep(my.cols()[ncl], length(ylvls)-ncl), ylvls[-1:-ncl]))
    }
  }
  dat[, ygrp:=factor(ygrp, levels=ylvls)]

  if (!is.null(mdat)) {
    mdat <- copy(mdat)
    if (!is.null(mdat.xvar)) tmp <- mdat.xvar else if (!is.null(xvar)) tmp <- xvar else tmp <- names(mdat)[1]
    setnames(mdat, tmp, "x")
    if (mode=="count") mdat[, ntot:=ntot[x]] else if (!is.null(ntot)) setnames(mdat, ntot, "ntot")
    if (!is.null(xgrp)) {
      setnames(mdat, xgrp, "grp")
      if (!is.factor(mdat$grp)) mdat[, grp:=factor(grp)]
    }
    if (!is.null(xcol)) {
      setnames(mdat, xcol, "col")
      if (!is.factor(mdat$col)) mdat[, col:=factor(col)]
      xclrs <- mdat[, setNames(my.cols(levels(col))[col], x)]
    } else xclrs <- NULL
    if (!is.null(xsep)) setnames(mdat, xsep, "sep")
  } else {
    mdat <- data.table(x=unique(dat$x))
    if (mode=="count") {
      mdat[, ntot:=ntot[x]]
    } else if (!is.null(ntot)) {
      if (is.null(names(ntot))) mdat[, ntot:=ntot] else mdat[, ntot:=ntot[x]]
    }
    if (!is.null(xgrp)) {
      if (is.null(names(xgrp))) mdat[, grp:=xgrp] else mdat[, grp:=xgrp[x]]
      if (!is.factor(mdat$grp)) mdat[, grp:=factor(grp)]
    }
    if (!is.null(xcol)) {
      if (is.null(names(xcol))) mdat[, col:=xcol] else mdat[, col:=xcol[x]]
      if (!is.factor(mdat$col)) mdat[, col:=factor(col)]
      xclrs <- mdat[, setNames(my.cols(levels(col))[col], x)]
    } else xclrs <- NULL
    if (!is.null(xsep)) {
      if (is.null(names(xsep))) mdat[, sep:=xsep] else mdat[, sep:=xsep[x]]
    }
  }
  if (!any(c("ntot", "grp", "col", "sep") %in% names(mdat))) mdat <- NULL

  if (!is.null(xgrp)) {
    grps <- levels(mdat$grp)
    if (length(dend.scale)==1) dend.scale <- rep(dend.scale, length(grps))
    if (is.null(names(dend.scale))) names(dend.scale) <- grps
    p.list <- lapply(grps, function(i) {
      mdat1 <- mdat[grp==i]
      dat1 <- dat[x %in% mdat1$x]
      # if w/o dendrogram, I first plot each group w/o Y-axis and later add a common Y-axis, this will make it easier to set a useable default rs.xgrp
      # but if with dendrogram, I cannot get the above approach to work so I include Y-axis in the first/last plot; as a result manual adjustment of rs.xgrp is usually necessary
      if (dendro) {
        .plot.fracs(dat1, xlab=NULL, ylab=ylab, facet.tit=i, lab=lab, lab.size=lab.size, lab.frac.ge=lab.frac.ge, clrs=clrs, mdat=mdat1, xclrs=xclrs, xord=xord, clust.fun=clust.fun, dendro=dendro, dend.scale=dend.scale[i], rs.dend=rs.dend, ori=dir, no.axs=switch(dir, h=i!=grps[length(grps)], v=i!=grps[1]), no.nlab=switch(dir, h=i!=grps[1], v=i!=grps[length(grps)]), lgd.pos="none")
      } else .plot.fracs(dat1, xlab=NULL, ylab=ylab, facet.tit=i, lab=lab, lab.size=lab.size, lab.frac.ge=lab.frac.ge, clrs=clrs, mdat=mdat1, xclrs=xclrs, xord=xord, clust.fun=clust.fun, dendro=dendro, dend.scale=dend.scale[i], rs.dend=rs.dend, ori=dir, no.axs=TRUE, no.nlab=switch(dir, h=i!=grps[1], v=i!=grps[length(grps)]), lgd.pos="none")
    })
    rel <- mdat[, .(n=.N), by=grp][order(grp), n]
    rel[rel %in% 1:2] <- rel[rel %in% 1:2]+0.4
    rel[rel %in% 3:4] <- rel[rel %in% 3:4]+0.2
    if (!is.null(rs.xgrp)) rel <- rel+rs.xgrp
    if (lgd.pos=="bottom") {
      lgd.mgn <- switch(dir, h=c(-10,10,0,0), v=c(0,10,0,0))
    } else lgd.mgn <- NULL
    lgd <- cowplot::get_legend(.plot.fracs(dat, facet.tit="", clrs=clrs, mdat=mdat, dendro=FALSE, ori=dir, lgd.tit=lgd.tit, lgd.pos=lgd.pos, lgd.mgn=lgd.mgn, ...))
    if (dir=="h") {
      rel[1] <- rel[1]+1
      if (dendro) {
        tmp <- 1+0.0075*min(sum(rel), 50)/max(rel[length(rel)]/sum(rel), 0.1)
        rel[length(rel)] <- rel[length(rel)]+tmp
      }
      p <- cowplot::plot_grid(plotlist=p.list, ncol=1, align="v", axis="lr", rel_heights=rel)
      if (!dendro) {
        mdat1 <- mdat[grp==grps[length(grps)]]
        dat1 <- dat[x %in% mdat1$x]
        y.axs <- .plot.fracs(dat1, xlab=NULL, ylab=ylab, facet.tit=grps[length(grps)], lab=lab, lab.size=lab.size, lab.frac.ge=lab.frac.ge, clrs=clrs, mdat=mdat1, xclrs=xclrs, xord=xord, clust.fun=clust.fun, dendro=dendro, dend.scale=dend.scale[grps[length(grps)]], rs.dend=rs.dend, ori=dir, no.nlab=TRUE, lgd.pos="none", ret.axs=TRUE)
        p <- gridExtra::arrangeGrob(p, bottom=y.axs)
      }
      if (!is.null(xlab)) {
        x.tit <- grid::textGrob(xlab, gp=grid::gpar(fontsize=10), rot=90)
        p <- gridExtra::arrangeGrob(p, left=x.tit)
      }
      if (lgd.pos=="bottom") {
        if (is.null(rs.lgd)) rs.lgd <- c(0.2+min(0.25*uniqueN(dat$x)+0.45, 2), 0.2)
        p <- cowplot::plot_grid(p, lgd, ncol=1, rel_heights=rs.lgd)
      } else if (lgd.pos=="right") {
        if (is.null(rs.lgd)) rs.lgd <- c(5, 0.1+0.1*max(nchar(ylvls)))
        p <- cowplot::plot_grid(p, lgd, nrow=1, rel_widths=rs.lgd)
      }
    } else if (dir=="v") {
      rel[length(rel)] <- rel[length(rel)]+1
      if (dendro) {
        tmp <- 1+0.0075*min(sum(rel), 50)/max(rel[1]/sum(rel), 0.1)
        rel[1] <- rel[1]+tmp
      }
      p <- cowplot::plot_grid(plotlist=p.list, nrow=1, align="h", axis="tb", rel_widths=rel)
      if (!dendro) {
        mdat1 <- mdat[grp==grps[1]]
        dat1 <- dat[x %in% mdat1$x]
        y.axs <- .plot.fracs(dat1, xlab=NULL, ylab=ylab, facet.tit=grps[1], lab=lab, lab.size=lab.size, lab.frac.ge=lab.frac.ge, clrs=clrs, mdat=mdat1, xclrs=xclrs, xord=xord, clust.fun=clust.fun, dendro=dendro, dend.scale=dend.scale[grps[1]], rs.dend=rs.dend, ori=dir, no.nlab=TRUE, lgd.pos="none", ret.axs=TRUE)
        p <- gridExtra::arrangeGrob(p, left=y.axs)
      }
      if (!is.null(xlab)) {
        x.tit <- grid::textGrob(xlab, gp=grid::gpar(fontsize=10))
        p <- gridExtra::arrangeGrob(p, bottom=x.tit)
      }
      if (lgd.pos=="bottom") {
        if (is.null(rs.lgd)) rs.lgd <- 0.1
        p <- cowplot::plot_grid(p, lgd, ncol=1, rel_heights=c(1, rs.lgd))
      } else if (lgd.pos=="right") {
        if (is.null(rs.lgd)) rs.lgd <- c(0.2+min(0.25*uniqueN(dat$x)+0.45, 2), 0.1+0.1*max(nchar(ylvls)))
        p <- cowplot::plot_grid(p, lgd, nrow=1, rel_widths=rs.lgd)
      }
    }
    if (!is.null(tit)) {
      p.tit <- grid::textGrob(tit, gp=grid::gpar(fontsize=12))
      p <- cowplot::plot_grid(gridExtra::arrangeGrob(p, top=p.tit))
    }
  } else {
    if (is.null(rs.lgd)) rs.lgd <- 0.1
    p <- .plot.fracs(dat, xlab=xlab, ylab=ylab, tit=tit, lab=lab, lab.size=lab.size, lab.frac.ge=lab.frac.ge, clrs=clrs, mdat=mdat, xord=xord, clust.fun=clust.fun, dendro=dendro, dend.scale=dend.scale, rs.dend=rs.dend, ori=dir, lgd.tit=lgd.tit, lgd.pos=lgd.pos, rs.lgd=rs.lgd, ...)
  }
  p
}


sc.dotplotly <- function(dat, gns=NULL, mdat, grp, blk=NULL, std=TRUE, exp=TRUE, f1=mean, t1=NULL, f2=mean, t2=NULL, expr.cutoff=0, ncells.cutoff=3, gene.anno=NULL, grp.anno=NULL, gene.txt=NULL, grp.txt=NULL, grp.name="cluster", xlab=str_to_title(grp.name), flip=FALSE, lo.args=list(), ...) {
  # interactive dotplot with heatmaply for single-cell gene expression data
  # dat, gns, mdat, grp, blk, exp, f1, t1, f2, t2, expr.cutoff, ncells.cutoff: passed to `summ.expr.by.grp`
  # std: whether to plot the standardized (i.e. scaled) expression values across groups; if TRUE and `t1` and `t2` are not specified, will automatically set `t1` or `t2` to `scale` as appropriate depending on whether `blk` is given; if at least one of `t1` and `t2` is given, the transformation will be determined by `t1` and `t2`
  # independently, `std` will also have effect on the dot color scheme (3-color or 2-color) and label_names
  # gene.anno: data.frame or data.table, first column should contain gene symbols/names, each subsequent column is a variable used to annotate the genes
  # grp.anno: similar to gene.anno, annotation for the groups
  # gene.txt: a named vector of common hovertext for each gene, named by the gene symbols
  # grp.txt: similar to gene.txt, a named vector of common hovertext for each group, named by the group levels
  # grp.name: what the groups represent, e.g. "cluster"; this is used in the hovertext
  # flip: by default will plot gene-by-group dot plot, and xlab is for group; if flip=TRUE will flip the X and Y axes (and xlab will become ylab)
  # lo.args: pased to plotly::layout
  # ...: passed to heatmaply

  for (pkg in c("heatmaply", "plotly")) {
    if (!requireNamespace(pkg, quietly=TRUE)) {
      stop(sprintf("Package \"%s\" needed for this function to work.", pkg))
    }
  }

  if (is.null(t1) && is.null(t2) && std) {
    if (is.null(blk)) t2 <- scale else t1 <- scale
  }
  tmp <- summ.expr.by.grp(dat, gns, mdat, grp, blk, exp, f1, t1, f2, t2, pct=TRUE, expr.cutoff, ncells.cutoff, ret.no.t=TRUE, ret.grp.sizes=TRUE)
  # todo: use ret.no.t=TRUE to also plot a non-scaled plot (provide this as an option)
  avg <- tmp$avg
  #avg[is.na(avg)] <- 0
  pct <- tmp$pct
  gns <- rownames(avg)
  if (!is.null(t1) || !is.null(t2)) txt.mat <- sprintf("log expr: %.3g\n", tmp$avg.no.trans) else txt.mat <- ""
  txt.mat <- paste0(txt.mat, sprintf("# of cells in %s: %d\n", grp.name, rep(tmp$grp.ncells[colnames(avg)], each=nrow(avg))))
  if (!is.null(blk)) {
    txt.mat <- paste0(txt.mat, sprintf("# of samples in %s: %d\n", grp.name, rep(tmp$grp.nblks[colnames(avg)], each=nrow(avg))))
  }
  if (!is.null(gene.txt)) txt.mat <- paste0(txt.mat, paste0(rep(gene.txt[gns], ncol(avg)), "\n"))
  if (!is.null(grp.txt)) txt.mat <- paste0(txt.mat, paste0(rep(grp.txt[colnames(avg)], each=nrow(avg)), "\n"))
  dim(txt.mat) <- dim(avg)
  if (flip) {
    avg <- t(avg)
    pct <- t(pct)
    txt.mat <- t(txt.mat)
  }
  nr <- nrow(avg)
  nch.r <- max(nchar(rownames(avg)), na.rm=TRUE)
  nc <- ncol(avg)
  nch.c <- max(nchar(colnames(avg)), na.rm=TRUE)
  if (std) f.color <- ggplot2::scale_color_gradient2(low="deepskyblue3", mid="grey90", high="indianred3", midpoint=0) else f.color <- ggplot2::scale_color_gradient(low="grey90", high="deepskyblue3") # for dot plot, need to use scale_color rather than scale_fill (the latter for normal heatmap)
  args <- list(
    x=avg,
    node_type="scatter",
    point_size_mat=pct,
    label_names=if (flip) c(grp.name, "gene", ifelse(std, "std expr", "log expr")) else c("gene", grp.name, ifelse(std, "std expr", "log expr")),
    label_format_fun=function(...) sprintf("%.3g", ...),
    point_size_name="% expr",
    custom_hovertext=txt.mat,
    scale_fill_gradient_fun=f.color,
    ...
  )
  if (flip) args$ylab <- xlab else args$xlab <- xlab
  if (!"column_text_angle" %in% names(args)) args$column_text_angle <- 40
  if (!"grid_size" %in% names(args)) args$grid_size <- 0.05
  if (!"branches_lwd" %in% names(args)) args$branches_lwd <- 0.3
  if (!is.null(gene.anno)) {
    gene.anno <- dt2mat(gene.anno, data=FALSE)
    if (is.null(dim(gene.anno))) gene.anno <- cbind(anno=gene.anno)
    if (flip) {
      args$col_side_colors <- gene.anno[match(gns, rownames(gene.anno)), , drop=FALSE]
      args$col_side_colors[is.na(args$col_side_colors)] <- "NA"
      args$col_side_palette <- colorspace::rainbow_hcl
    } else {
      args$row_side_colors <- gene.anno[match(gns, rownames(gene.anno)), , drop=FALSE]
      args$row_side_colors[is.na(args$row_side_colors)] <- "NA"
      args$row_side_palette <- colorspace::rainbow_hcl
    }
  }
  if (!is.null(grp.anno)) {
    grp.anno <- dt2mat(grp.anno, data=FALSE)
    if (flip) {
      args$row_side_colors <- grp.anno[match(rownames(avg), rownames(grp.anno)), , drop=FALSE]
      args$row_side_palette <- colorspace::rainbow_hcl
    } else {
      args$col_side_colors <- grp.anno[match(colnames(avg), rownames(grp.anno)), , drop=FALSE]
      args$col_side_palette <- colorspace::rainbow_hcl
    }
  }
  if (flip) {
    if (!is.null(grp.anno)) nc1 <- ncol(grp.anno) else nc1 <- 0
    if (!is.null(gene.anno)) nr1 <- ncol(gene.anno) else nr1 <- 0
  } else {
    if (!is.null(gene.anno)) nc1 <- ncol(gene.anno) else nc1 <- 0
    if (!is.null(grp.anno)) nr1 <- ncol(grp.anno) else nr1 <- 0
  }
  if (!"subplot_widths" %in% names(args)) {
    args$subplot_widths <- c(nc/(nc+0.7*nc1+1), 0.7*nc1/(nc+0.7*nc1+1), 1/(nc+0.7*nc1+1))
    args$subplot_widths <- args$subplot_widths[args$subplot_widths!=0]
  }
  if (!"subplot_heights" %in% names(args)) {
    args$subplot_heights <- c(1.5/(nr+nr1+1.5), nr1/(nr+nr1+1.5), nr/(nr+nr1+1.5))
    args$subplot_heights <- args$subplot_heights[args$subplot_heights!=0]
  }
  if (is.null(lo.args$width)) lo.args$width <- 50*(nc+0.7*nc1+1)+10*nch.r+150
  if (is.null(lo.args$height)) lo.args$height <- 27*(nr+nr1+1.5)+10*nch.c+50
  if (is.null(lo.args$autosize)) lo.args$autosize <- TRUE

  p <- do.call(heatmaply::heatmaply, args)
  lo.args$p <- p
  do.call(plotly::layout, lo.args)
}


#####  functions and helper functions for heatmap plotting with ComplexHeatmap #####

.xcut <- function(x, ys, x.seps="s", from0=FALSE, sym=FALSE, m3d=TRUE, ys.m3d=NULL, trans=NULL) {
  # helper function to map values in x to a scale of y values given in `ys`; x is a numeric vector or matrix, ys is a vector representing the scale
  # length(ys) should be >=2, the first and last values are mapped to min(x) and max(x) if m3d=FALSE, otherwise to median(x)+-3*mad(x); the "internal" values of `ys` are mapped to x values given in `x.seps`
  # x.seps can also be "q", in which case quantiles of x will be mapped to `ys`, or "s" in which case equally spaced cut points from min-max or median+-3mad will be mapped to `ys`
  # trans: a function to transform x and x.seps as a first step before mapping
  # from0: if TRUE then all x (after transformation) should >=0
  # sym: whether the result mapping should be made symmetric, only effective if from0=FALSE and the length of ys is an odd number
  # ys.m3d: only effective if m3d=TRUE, the y scale values used for out-of-m3d range x values

  if (is.numeric(x.seps) && length(x.seps)!=length(ys)-2) stop("`x.seps` should have a length of length(ys)-2.")
  q <- is.character(x.seps) && x.seps=="q"
  if (!is.null(trans)) {
    if (!is.function(trans)) stop("`trans` should be a function.")
    x <- trans(x)
    if (is.numeric(x.seps)) x.sep <- trans(x.seps)
  }
  if (from0 && any(x<0, na.rm=TRUE)) stop("`from0` is true but the data contains negative values.")
  x1 <- x[is.finite(x)] # will remove na's as well
  ma <- max(x1)
  mi <- min(x1)
  if (m3d) {
    m <- median(x1)
    d3 <- 3*mad(x1)
    hi <- m+d3
    lo <- m-d3
    if (hi>ma) hi <- ma else hi <- max(x1[x1<hi])
    if (lo<mi) lo <- mi else lo <- min(x1[x1>lo])
    if (from0) {
      if (q) xs <- c(quantile(c(0,x1[x1<=hi]), seq(0, 1, len=length(ys))), ma) else xs <- c(seq(0, hi, len=length(ys)), ma)
      if (!is.null(ys.m3d)) ys <- c(ys, ys.m3d) else ys <- c(ys, ys[length(ys)])
    } else {
      if (q) xs <- c(mi, quantile(x1[x1>=lo & x1<=hi], seq(0, 1, len=length(ys))), ma) else xs <- c(mi, seq(lo, hi, len=length(ys)), ma)
      if (!is.null(ys.m3d)) {
        if (length(ys.m3d)!=2) stop("`ys.m3d` should have a length of 2.")
        ys <- c(ys.m3d[1], ys, ys.m3d[2])
      } else ys <- c(ys[1], ys, ys[length(ys)])
    }
  } else {
    if (q) xs <- quantile(x1, seq(0, 1, len=length(ys))) else xs <- seq(mi, ma, len=length(ys))
    if (from0) xs[1] <- 0
  }
  if (is.numeric(x.seps)) {
    rng <- xs[c(1, length(xs))]
    if (m3d) {
      if (from0) xs[2:(length(xs)-2)] <- x.seps else xs[3:(length(xs)-2)] <- x.seps
    } else xs[2:(length(xs)-1)] <- x.seps
    xs <- sort(xs)
    idx <- xs<rng[1] | xs>rng[2]
    if (any(idx)) {
      xs <- xs[!idx]
      ys <- ys[!idx]
    }
    if (!from0 && sym && length(xs)%%2==1) {
      i <- length(xs) %/% 2
      mid <- xs[i+1]
      if ((mid-xs[1])<(xs[length(xs)]-mid)) xs[1:i] <- 2*mid-xs[length(xs):(i+2)] else xs[(i+2):length(xs)] <- 2*mid-xs[i:1]
    }
  }
  if (xs[1]==xs[2]) {
    xs <- xs[-1]
    ys <- ys[-1]
  }
  if (xs[length(xs)]==xs[length(xs)-1]) {
    xs <- xs[-length(xs)]
    ys <- ys[-length(ys)]
  }
  list(x=xs, y=ys)
}

.fcolor <- function(x, cols, pal=1, x.seps="s", from0=FALSE, sym=FALSE, m3d=TRUE, cols.m3d=NULL, na.col="grey50", trans=NULL, ...) {
  # helper function to map values in x to a color scale in `cols`; x is a numeric vector or matrix, cols is a vector representing the scale
  # if length(cols)==1, then cols can be 2 or 3 corresponding to preset colors (2 is based on 2 colors, and 3 is based on 3 colors for bi-directional data with mid point 0)
  # pal: 1 or 2, preset palette for when cols=2 or 3, otherwise will be ignored
  # other arguments passed to .xcut; ... passed to circlize::colorRamp2

  args <- list(...)
  if (length(cols)==1) {
    if (cols==2) {
      if (from0) {
        if (pal==1) cols <- c("grey90", "red2") else cols <- c("grey90", "indianred3")
        if (m3d) cols.m3d <- "darkred"
      } else {
        if (pal==1) cols <- c("blue2", "orange") else cols <- c("dodgerblue4", "goldenrod1")
        if (m3d) {
          if (pal==1) cols.m3d <- c("darkblue", "yellow") else cols.m3d <- c("midnightblue", "yellow")
        }
      }
    } else if (cols==3) {
      x.seps <- 0
      from0 <- FALSE
      sym <- TRUE
      if (pal==1) cols <- c("blue2", "grey90", "red2") else cols <- c("deepskyblue3", "grey90", "indianred3")
      if (m3d) {
        if (pal==1) cols.m3d <- c("darkblue", "darkred") else cols.m3d <- c("royalblue4", "darkred")
      }
    }
  }
  m <- .xcut(x, cols, x.seps, from0, sym, m3d, cols.m3d, trans)
  if (uniqueN(m$x)==1) {
    # the results need to be a function similar to that returned by circlize::colorRamp2
    f <- function(i) {
      res <- rep(colorRampPalette(m$y)(3)[2], length(i))
      res[is.na(i)] <- na.col
      res
    }
    attr(f, "breaks") <- m$x[1]*c(0.99, 1, 1.01)
    attr(f, "colors") <- rep(colorRampPalette(m$y)(3)[2], 3)
  } else {
    # instead of directly return the function generated by circlize::colorRamp2, need to handle trans and NA (for some reason the na_col in Heatmap does not work?)
    f <- function(i) {
      f1 <- circlize::colorRamp2(m$x, m$y, ...)
      if (is.null(trans)) res <- f1(i) else res <- f1(trans(i))
      res[is.na(i)] <- na.col
      res
    }
    attr(f, "breaks") <- m$x
    attr(f, "colors") <- m$y
  }
  attr(f, "transparency") <- if ("transparency" %in% names(args)) args$transparency else 0
  attr(f, "space") <- if ("space" %in% names(args)) args$space else "LAB"
  f
}

.fcolor1 <- function(x, cols=NULL, mute=FALSE, ...) {
  # helper function to generate colors for a vector x
  if (is.numeric(x)) {
    if (is.null(cols)) {
      if (any(x>0, na.rm=TRUE) && any(x<0, na.rm=TRUE)) .fcolor(x, cols=3, ...) else .fcolor(x, cols=2, ...)
    } else .fcolor(x, cols, ...)
  } else {
    # unlike na_col in Heatmap, na_col does work in HeatmapAnnotation so handling NA here is unnecessary, but I keep it anyway
    res <- my.cols(x, na.color="white", mute=mute)
    names(res)[is.na(names(res))] <- "NA" # need to do this otherwise HeatmapAnnotation will give error
    res
  }
}

.fsize <- function(mat, sizes=c(0.1, 0.55), x.seps="s", from0=TRUE, sym=FALSE, m3d=FALSE, sizes.m3d=NULL, trans=NULL) {
  # helper function to map values in mat to a scale of dot sizes (radius) for ComplexHeatmap; mat is a numeric matrix, sizes is a vector representing the scale
  m <- .xcut(mat, sizes, x.seps, from0, sym, m3d, sizes.m3d, trans)
  xs <- m$x
  sizes <- m$y
  if (!is.null(trans)) mat <- trans(mat)

  function(i=NULL, j=NULL, x=NULL) {
    if (is.null(x)) x <- mat[i, j] else if (!is.null(trans)) x <- trans(x)
    if (xs[1]==xs[length(xs)]) {
      if (is.infinite(x) && x>0) {
        s <- sizes[length(sizes)]
      } else if (is.infinite(x) && x<0) {
        s <- sizes[1]
      } else s <- (sizes[1] + sizes[length(sizes)])/2
    } else s <- approx(xs, sizes, x, rule=2, ties="ordered")$y
    if (is.na(s)) NA else s*unit(3, "mm")
  }
}

plot.hm <- function(mat, smat=NULL, sp=FALSE, xlab=NULL, ylab=NULL, x.anno=NULL, y.anno=NULL, name=" ", sname=" ", cols=3, pal=1, seps="s", from0=FALSE, sym=TRUE, m3d=TRUE, cols.m3d=NULL, na.col="grey50", trans=NULL, sizes=c(0.1, 0.55), s.seps="s", s.from0=TRUE, s.sym=FALSE, s.m3d=FALSE, sizes.m3d=NULL, s.trans=NULL, colf.anno=.fcolor1, cellf=NULL, lab.pos="bl", clust="xy", dend.pos="tl", anno.pos="bl", lgd.pos=c("b", "r"), lgd.ori=c("default", "h", "v"), anno.lgd.pos=c("b", "r"), anno.lgd.ori=c("default", "h", "v"), lgd.key.nmax=5, pack.lgd=c("default", "h", "v"), merge.lgd=NULL, ...) {
  # plot heatmap (or dot plot) with ComplexHeatmap
  # mat: for heatmap color; smat: for dot size, if provided will do dot plot; will assume that mat and smat have the same dimensions and row/column orders
  # sp: if TRUE, will assume that smat contains adjusted P values; can provide significance cutoff in s.seps, e.g. s.seps=0.05
  # name and sname are legend labels for mat and smat respectively
  # x.anno and y.anno: data.frame or data.table, first column should contain IDs (i.e. col/rownames of mat), each subsequent column is a variable used to annotate the columns/rows
  # cols to trans passed to .fcolor for heatmap color
  # sizes to s.trans passed to .fsize for dot size
  # colf.anno: a function for generating annotation colors
  # cellf: pass a custom cellf replacing that defined in this function, the env of cellf will be changed to the runtime env of this function, so can refer to mat, smat, rf, colf, etc.
  # clust: if to cluster neither columns nor rows, set to ""

  lgd.pos <- match.arg(lgd.pos)
  lgd.ori <- match.arg(lgd.ori)
  anno.lgd.pos <- match.arg(anno.lgd.pos)
  anno.lgd.ori <- match.arg(anno.lgd.ori)
  pack.lgd <- match.arg(pack.lgd)
  if (is.null(merge.lgd)) merge.lgd <- lgd.pos==anno.lgd.pos
  if (lgd.ori=="default") lgd.ori <- switch(lgd.pos, b="horizontal", r="vertical") else lgd.ori <- switch(lgd.ori, h="horizontal", v="vertical")
  if (anno.lgd.ori=="default") anno.lgd.ori <- switch(anno.lgd.pos, b="horizontal", r="vertical") else anno.lgd.ori <- switch(anno.lgd.ori, h="horizontal", v="vertical")
  lgd.key.nr <- switch(anno.lgd.pos, b=NULL, r=lgd.key.nmax)
  lgd.key.nc <- switch(anno.lgd.pos, b=lgd.key.nmax, r=NULL)
  pos.map <- c(t="top", r="right", b="bottom", l="left")
  lgd.pos <- pos.map[lgd.pos]
  anno.lgd.pos <- pos.map[anno.lgd.pos]

  colf <- .fcolor(mat, cols, pal, seps, from0, sym, m3d, cols.m3d, na.col, trans)

  if (!is.null(smat)) {
    if (sp) {
      smat[smat<2.2e-16] <- 2.2e-16
      smat <- -log10(smat)
      if (sname==" ") sname <- expression(P[adj])
      if (is.character(s.seps)) {
        p.cut <- 0.05
        s.seps <- -log10(0.05)
        sizes <- c(0.1, 0.2, 0.55)
      } else {
        p.cut <- s.seps[1]
        s.seps <- -log10(s.seps)
        if (length(sizes)==2) sizes <- c(0.1, 0.2, 0.55)
      }
    }
    rf <- .fsize(smat, sizes, s.seps, from0=TRUE, sym=FALSE, m3d=FALSE, sizes.m3d, trans=NULL)
    if (is.null(cellf)) {
      if (sp) {
        cellf <- function(j, i, x, y, width, height, fill) {
          grid.rect(x=x, y=y, width=width, height=height, gp=gpar(col=NA, fill=NA))
          grid.circle(x=x, y=y, r=rf(i, j), gp=gpar(fill=colf(mat[i, j]), col=if (smat[i, j]>s.seps[1]) "black" else NA))
        }
      } else {
        cellf <- function(j, i, x, y, width, height, fill) {
          grid.rect(x=x, y=y, width=width, height=height, gp=gpar(col=NA, fill=NA))
          grid.circle(x=x, y=y, r=rf(i, j), gp=gpar(fill=colf(mat[i, j]), col=NA))
        }
      }
    } else {
      environment(cellf) <- environment()
    }
  } else cellf <- NULL

  # temporary workaround for clustering involving Inf
  mat1 <- mat
  if (clust!="" && any(is.infinite(mat))) {
    if (any(is.infinite(mat[mat>0]))) {
      tmp <- 10*max(mat[is.finite(mat)])
      if (tmp==0) tmp <- 10*max(abs(mat[is.finite(mat)]))
      if (tmp==0) tmp <- 1
      mat1[is.infinite(mat) & mat>0] <- tmp
    }
    if (any(is.infinite(mat[mat<0]))) {
      tmp <- 10*min(mat[is.finite(mat)])
      if (tmp==0) tmp <- 10*max(abs(mat[is.finite(mat)]))
      if (tmp==0) tmp <- -1
      mat1[is.infinite(mat) & mat<0] <- tmp
    }
  }

  if (!is.null(y.anno)) {
    y.anno <- y.anno[match(rownames(mat), y.anno[[1]]), -1, drop=FALSE] # drop=FALSE compatible with both data.table and data.frame
    row.cols <- lapply(y.anno, colf.anno)
    row.ha <- HeatmapAnnotation(
      which="row", df=y.anno, col=row.cols, na_col="white",
      simple_anno_size=unit(2.5,"mm"), gap=unit(2,"points"),
      show_annotation_name=TRUE, annotation_name_side=if (is.null(x.anno)) "bottom" else "top", annotation_name_gp=gpar(fontsize=9, fontface="plain", srt=40), annotation_name_rot=40,
      annotation_legend_param=list(direction=anno.lgd.ori, nrow=lgd.key.nr, ncol=lgd.key.nc, by_row=anno.lgd.ori=="horizontal", grid_width=unit(3, "mm"), grid_height=unit(3, "mm"), title_gp=gpar(fontsize=9, fontface="plain"), labels_gp=gpar(fontsize=8, fontface="plain"))
    )
  } else row.ha <- NULL

  if (!is.null(x.anno)) {
    x.anno <- x.anno[match(colnames(mat), x.anno[[1]]), -1, drop=FALSE] # drop=FALSE compatible with both data.table and data.frame
    col.cols <- lapply(x.anno, colf.anno)
    col.ha <- HeatmapAnnotation(
      which="column", df=x.anno, col=col.cols, na_col="white",
      simple_anno_size=unit(2.5,"mm"), gap=unit(2,"points"),
      show_annotation_name=TRUE, annotation_name_side="left", annotation_name_gp=gpar(fontsize=9, fontface="plain"),
      annotation_legend_param=list(direction=anno.lgd.ori, nrow=lgd.key.nr, ncol=lgd.key.nc, by_row=anno.lgd.ori=="horizontal", grid_width=unit(3, "mm"), grid_height=unit(3, "mm"), title_gp=gpar(fontsize=9, fontface="plain"), labels_gp=gpar(fontsize=8, fontface="plain"))
    )
  } else col.ha <- NULL

  hm <- Heatmap(mat1,
    col=colf,
    na_col="grey50",
    name=name,
    width=ncol(mat)*unit(if (any(grepl("\n", colnames(mat)))) 7 else 5, "mm"),
    height=nrow(mat)*unit(5, "mm"),
    row_title=ylab,
    row_title_gp=gpar(fontsize=10),
    row_title_side=pos.map[str_extract(lab.pos, "[lr]")],
    row_names_gp=gpar(fontsize=9),
    row_names_side=pos.map[str_extract(lab.pos, "[lr]")],
    cluster_rows=grepl("y", clust),
    row_dend_side=if (grepl("y", clust)) pos.map[str_extract(dend.pos, "[lr]")] else "left",
    row_dend_width=unit(7, "mm"),
    row_dend_gp=gpar(lwd=0.2),
    column_title=xlab,
    column_title_gp=gpar(fontsize=10),
    column_title_side=pos.map[str_extract(lab.pos, "[tb]")],
    column_names_gp=gpar(fontsize=9),
    column_names_side=pos.map[str_extract(lab.pos, "[tb]")],
    column_names_rot=40,
    cluster_columns=grepl("x", clust),
    column_dend_side=if (grepl("x", clust)) pos.map[str_extract(dend.pos, "[tb]")] else "top",
    column_dend_height=unit(7, "mm"),
    column_dend_gp=gpar(lwd=0.2),
    border=TRUE,
    border_gp=gpar(lwd=0.2),
    rect_gp=if (is.null(smat)) gpar(col=NA) else gpar(type="none"),
    cell_fun=cellf,
    left_annotation=if (grepl("l", anno.pos)) row.ha else NULL,
    right_annotation=if (grepl("r", anno.pos)) row.ha else NULL,
    top_annotation=if (grepl("t", anno.pos)) col.ha else NULL,
    bottom_annotation=if (grepl("b", anno.pos)) col.ha else NULL,
    show_heatmap_legend=is.null(smat),
    heatmap_legend_param=list(direction=lgd.ori, grid_width=unit(3, "mm"), grid_height=unit(3, "mm"), title_gp=gpar(fontsize=9, fontface="plain"), labels_gp=gpar(fontsize=8, fontface="plain")),
    ...
  )

  if (!is.null(smat)) {
    if (pack.lgd=="default") pack.lgd <- switch(lgd.pos, bottom="horizontal", right="vertical") else pack.lgd <- switch(pack.lgd, h="horizontal", v="vertical")
    lgd.c <- Legend(col_fun=colf, title=name, title_gp=gpar(fontsize=9, fontface="plain"), labels_gp=gpar(fontsize=8, fontface="plain"), direction=lgd.ori)
    if (sp) {
      tmp <- seq(0, max(smat, na.rm=TRUE), len=4)[-1]
      if (all(tmp<s.seps[1])) lbs <- tmp else if (all(tmp>=s.seps[1])) lbs <- c(-log10(0.2), tmp) else lbs <- sort(c(s.seps, tmp))
      grs <- lapply(lbs, function(i) function(x, y, w, h) grid.circle(x, y, r=rf(x=i), gp=gpar(fill="grey", col=if (i>=s.seps[1]) "black" else NA)))
      lbs <- ifelse(lbs==s.seps[1], paste0("<", p.cut), ifelse(lbs==-log10(0.2), paste0(">", p.cut), ifelse(lbs==-log10(2.2e-16), "<2.2e-16", sprintf("%.1g", 10^(-lbs)))))
      tmp <- !duplicated(lbs)
      lbs <- lbs[tmp]
      grs <- grs[tmp]
    } else {
      if (is.character(s.seps)) {
        if (s.seps=="s") {
          if (s.from0) lbs <- seq(0, max(smat[is.finite(smat)]), len=4) else lbs <- seq(min(smat[is.finite(smat)]), max(smat[is.finite(smat)]), len=4)
        } else {
          lbs <- quantile(smat[is.finite(smat)], seq(0, 1, len=4))
          if (s.from0) lbs[1] <- 0
        }
      } else {
        if (s.from0) lbs <- c(0, s.seps, max(smat[is.finite(smat)])) else lbs <- c(min(smat[is.finite(smat)]), s.seps, max(smat[is.finite(smat)]))
      }
      lbs <- unique(lbs)
      grs <- lapply(lbs, function(i) function(x, y, w, h) grid.circle(x, y, r=rf(x=i), gp=gpar(fill="grey", col=NA)))
      lbs <- sprintf("%.1f", lbs)
    }
    lgd.s <- Legend(title=sname, title_gp=gpar(fontsize=9, fontface="plain"), labels=lbs, labels_gp=gpar(fontsize=8, fontface="plain"), graphics=grs, nrow=switch(lgd.ori, horizontal=1, vertical=NULL), ncol=switch(lgd.ori, horizontal=NULL, vertical=1), by_row=lgd.ori=="horizontal", direction=lgd.ori)
    lgd <- packLegend(lgd.c, lgd.s, direction=pack.lgd)
    draw(hm, column_title=NULL, heatmap_legend_side=anno.lgd.pos, merge_legend=merge.lgd, annotation_legend_list=lgd, annotation_legend_side=lgd.pos)
  } else draw(hm, heatmap_legend_side=lgd.pos, merge_legend=merge.lgd, annotation_legend_side=anno.lgd.pos, legend_grouping="original")
}


sc.dotplot <- function(dat, gns=NULL, mdat, grp, blk=NULL, std=TRUE, exp=TRUE, f1=mean, t1=NULL, f2=mean, t2=NULL, expr.cutoff=0, ncells.cutoff=3, gene.anno=NULL, markers=NULL, flag=FALSE, grp.map=NULL, markers.wl=NULL, grp.anno=NULL, xlab, flip=TRUE, cols=NULL, pal=1, m3d=FALSE, ...) {
  # dotplot with ComplexHeatmap for single-cell gene expression data
  # dat, gns, mdat, grp, blk, exp, f1, t1, f2, t2, expr.cutoff, ncells.cutoff: passed to `summ.expr.by.grp`
  # if gns is NULL and `markers` is provided, will automatically set gns to unique(unlist(markers))
  # std: whether to plot the standardized (i.e. scaled) expression values across groups; if TRUE and `t1` and `t2` are not specified, will automatically set `t1` or `t2` to `scale` as appropriate depending on whether `blk` is given; if at least one of `t1` and `t2` is given, the transformation will be determined by `t1` and `t2`
  # independently, `std` will also have effect on the dot color scheme (3-color or 2-color) and label_names
  # gene.anno: a data.frame or data.table of gene annotation, first column should contain gene symbols/names
  # grp.anno: similar to gene.anno, group annotation, first column should contain group names
  # markers: a list of marker gene annotation, e.g. list(CD8T=c("CD8A", "CD8B")), if provided will add to gene.anno
  # flag: if TRUE and `markers` is provided, will highlight the cases where the gene expression in a group is not as expected (either unexpected high expression or unexpected low expression)
  # grp.map: a named list or vector mapping each group name (in grp) to the class (e.g. cell type) name(s) of `markers` (i.e. names(markers)); group names not in `grp.map` will be kept as is before trying to match names(markers); if NULL, will assume that the group names already match names(markers)
  # markers.wl: a list of whitelist markers to avoid flagging, e.g. list(B="CD4"), these markers are allowed to be expressed in the corresponding classes/cell types; if provided will automatically set flag=TRUE
  # if markers.wl is not NULL and a class/cell type is in `markers` but not in `markers.wl`, the class/cell type won't be included for flagging of unexpected high expression
  # xlab: if missing will set to `grp`, to not show xlab set it to NULL
  # flip: when FALSE will plot gene-by-group dot plot, and xlab is for group; if flip=TRUE (default) will flip the X and Y axes (and xlab will become ylab)
  # ...: passed to plot.hm

  if (is.null(t1) && is.null(t2) && std) {
    if (is.null(blk)) t2 <- scale else t1 <- scale
  }
  if (is.null(gns) && !is.null(markers)) gns <- unique(unlist(markers))
  if (missing(xlab)) xlab <- grp
  tmp <- summ.expr.by.grp(dat, gns, mdat, grp, blk, exp, f1, t1, f2, t2, pct=TRUE, expr.cutoff, ncells.cutoff, ret.no.t=FALSE, ret.grp.sizes=TRUE)
  # todo: use ret.no.t=TRUE to also plot a non-scaled plot (provide this as an option)
  avg <- tmp$avg
  pct <- tmp$pct
  # remove any genes whose value is NaN in all groups
  idx <- rowSums(!is.na(avg))>0
  if (any(!idx)) message(sprintf("These genes have NA values across all groups, they will be excluded:\n%s", paste(rownames(avg)[!idx], collapse=", ")))
  avg <- avg[idx, , drop=FALSE]
  pct <- pct[idx, , drop=FALSE]
  # add group sizes (cell and block/sample numbers) to group names
  grps <- colnames(avg)
  if (!is.null(blk)) {
    grps <- setNames(sprintf("%s (n=%d;%d)", grps, tmp$grp.ncells[grps], tmp$grp.nblks[grps]), grps)
  } else {
    grps <- setNames(sprintf("%s (n=%d)", grps, tmp$grp.ncells[grps]), grps)
  }
  colnames(avg) <- colnames(pct) <- grps
  if (!is.null(grp.anno)) grp.anno[[1]] <- grps[grp.anno[[1]]]
  # markers
  cellf <- NULL
  if (!is.null(markers)) {
    tmp <- rbindlist(lapply(markers, function(x) data.table(gene=x)), idcol="class")
    tmp <- tmp[, .(`marker of`=paste(sort(unique(class)), collapse=",")), by=gene]
    if (!is.null(gene.anno)) {
      names(gene.anno)[1] <- "gene"
      gene.anno <- merge(gene.anno, tmp, by="gene", all=TRUE)
    } else gene.anno <- tmp[, .(gene, `marker of`)]

    if (flag || !is.null(markers.wl)) {
      cutoffs <- apply(avg, 1, function(x) {
        x <- x[!is.na(x)]
        if (length(x)<=2) return(c(lo=Inf, hi=-Inf))
        set.seed(1)
        km <- kmeans(x, 2, nstart=5)
        hi <- min(x[km$cluster==which.max(km$centers)])
        lo <- max(x[km$cluster==which.min(km$centers)])
        c(lo=lo, hi=hi)
      })
      gns <- rownames(avg)
      grp.map <- c(as.list(grp.map), as.list(cn(setdiff(names(grps), names(grp.map)))))
      fmat <- sapply(setNames(names(grps), grps), function(g) {
        x <- setNames(rep(0, length(gns)), gns)
        cts <- grp.map[[g]]
        if (!any(cts %in% names(markers)) || (!is.null(markers.wl) && !any(cts %in% names(markers.wl)))) return(x)
        hi.exp <- gns %in% unlist(markers[cts])
        lo.obs <- avg[, grps[g]]<cutoffs["hi",]
        x[hi.exp & lo.obs] <- 1
        if (!is.null(markers.wl) && !any(cts %in% names(markers.wl))) return(x)
        lo.exp <- gns %in% setdiff(unlist(markers[!names(markers) %in% cts]), unlist(markers.wl[cts]))
        hi.obs <- avg[, grps[g]]>cutoffs["lo",]
        x[lo.exp & hi.obs] <- -1
        x
      })
      if (flip) fmat <- t(fmat)
      assign(".tmp.cellf.env", new.env(), envir=.GlobalEnv)
      .tmp.cellf.env$fmat <- fmat
      on.exit(rm(.tmp.cellf.env, envir=.GlobalEnv))
      cellf <- function(j, i, x, y, width, height, fill) {
        fmat <- .tmp.cellf.env$fmat
        grid.rect(x=x, y=y, width=width, height=height, gp=gpar(col=NA, fill=NA))
        grid.circle(x=x, y=y, r=rf(i, j), gp=gpar(fill=colf(mat[i, j]), col=if (fmat[i, j]==1) "red2" else if (fmat[i, j]==-1) "blue2" else NA))
      }
    }
  }

  if (flip) {
    avg <- t(avg)
    pct <- t(pct)
    x.anno <- gene.anno
    y.anno <- grp.anno
    ylab <- xlab
    xlab <- NULL
    clus <- "y"
  } else {
    x.anno <- grp.anno
    y.anno <- gene.anno
    ylab <- NULL
    clus <- "x"
  }
  if (is.null(cols)) {
    if (std) cols <- 3 else cols <- 2
  }

  tryCatch({
    plot.hm(avg, pct, sp=FALSE, xlab, ylab, x.anno, y.anno, name=ifelse(std, "std expr", "log expr"), sname="% expr", cols=cols, pal=pal, m3d=m3d, cellf=cellf, ...)
  }, error=function(e) {
    # if error, could be that clustering failed due to NA's
    message("Will try replacing all NA's with 0. Alternatively, try pass clust='' to disable clustering of rows and columns of the heatmap.")
    avg[is.na(avg)] <- 0
    pct[is.na(pct)] <- 0
    plot.hm(avg, pct, sp=FALSE, xlab, ylab, x.anno, y.anno, name=ifelse(std, "std expr", "log expr"), sname="% expr", cols=cols, pal=pal, m3d=m3d, cellf=cellf, ...)
  })
}

