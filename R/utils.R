## ----some utility functions for common small tasks----

env <- Sys.getenv
fp <- file.path
carg <- function(x=TRUE) commandArgs(trailingOnly=x)

co <- function(..., c="") paste(c(...), collapse=c)
cc <- function(...) paste(c(...), collapse=", ")

cn <- function(...) {
  # create a named vector with name being the same as vector content
  res <- c(...)
  names(res) <- res
  res
}

cn1 <- function(x, names.x) {
  names(x) <- names.x
  x
}

unq <- function(x) unique(unlist(x))

factor1 <- function(x) factor(x, levels=x)


rmt <- function(port=4321, ...) {
  # start rmote server
  rmote::start_rmote(port=port, ...)
}

rmt0 <- rmote::stop_rmote

pd <- rmote::plot_done

hgd <- function(port=4321, token=FALSE, ...) {
  # start httpgd server
  httpgd::hgd(port=port, token=token, ...)
}

hgd0 <- httpgd::hgd_close


hh <- function(x, nr=5, nc=nr) {
  # check the first nr rows and nc columns of a matrix-like object
  if (nr>nrow(x)) nr <- nrow(x)
  if (nc>ncol(x)) nc <- ncol(x)
  x[1:nr, 1:nc]
}


dt2mat <- function(x, rn=1, rm=1) {
  # convert a data.table to matrix, using the rn column as rownames, and remove columns as specified by rm
  # rn and rm can be indices or column names
  rns <- x[[rn]]
  res <- data.matrix(x[, -rm, with=FALSE])
  rownames(res) <- rns
  res
}


write.tab <- function(x, file, cn=TRUE, rn=FALSE, sep="\t", quote="auto", append=FALSE) {
  # write table with changed defaults
  tryCatch(fwrite(x, file=file, append=append, quote=quote, sep=sep, na="NA", col.names=cn, row.names=rn),
  	error=function(e) write.table(x, file=file, append=append, quote=ifelse(quote=="auto",FALSE,quote), sep=sep, na="NA", col.names=cn, row.names=rn))
}


rgx.esc <- function(x) {
  # convert string x to a fully-escaped regex
  str_replace_all(x, "(\\W)", "\\\\\\1")
}


im <- function(x, y, sep=" ; ", grep=FALSE) {
  # extension of %in%
  # x and y are both vectors and can contain multiple sep-separated items in each of their elements, im(x,y) is similar to x %in% y, but is TRUE as long as at least one item of the element of x is in the collection of all items of y (when grep=FALSE), or any item of the element of x matches any substring of any element of y (when grep=TRUE).

  if (grep) {
    #x.regex <- str_replace_all(x, sep, "|")
    x.regex <- str_split(x, sep)
    x.regex <- lapply(x.regex, function(x) paste(regex.escape(x), collapse="|"))
    y.combined <- paste(y, collapse=" ; ")
    res <- sapply(x.regex, function(i) grepl(i, y.combined, ignore.case=TRUE))
  } else {
    x.items <- str_split(x, sep)
    y.all.items <- unique(unlist(str_split(y, sep)))
    res <- sapply(x.items, function(i) any(i %in% y.all.items))
    res <- unname(res)
  }
  return(res)

}


mmatch <- function(x, y, return="l", simplify=TRUE, cplx=FALSE, sep=" ; ", grep=FALSE) {
  # multiple match
  # similar to match(x,y), but match all occurences of each element of x in y
  # return="l" for returning a list, or "dt" for returning a data.table, or "v" for vector
  # when x has length 1 and simplify=TRUE, return a vector
  # when return a vector and simplify=TRUE, return the unique y indeces

  ### Besides when cplx=TRUE, both x and y can contain multiple sep-separated items in each of their elements, and mmatch(x,y) gives a match as long as any item of the element of x matches any item of the element of y (when grep=FALSE), or any item of the element of x matches any substring of any element of y (when grep=TRUE).

  res <- lapply(x, function(xi) {
    if (cplx) {
      ind.y <- which(im(y, xi, sep=sep, grep=grep))
    } else ind.y <- which(y==xi)
    if (length(ind.y)==0) return(NA)
    return(ind.y)
  })

  if (length(x)==1 && simplify) return(unlist(res))
  if (return=="v") {
    res <- unlist(res)
    if (simplify) res <- unique(res)
  }
  if (return=="dt") {
    res <- lapply(res, function(i) data.table(ind.y=i))
    res <- rbindlist(res, idcol="ind.x")
  }

  return(res)

}


replace.na <- function(DT, to=0, col.mode="numeric", in.place=TRUE) {
  # change all the NAs in the columns of the type col.mode in a data.table into the value of to, by default, change all numeric NAs to 0s.
  if (!in.place) DT <- copy(DT)
  for (j in seq_along(DT)) {
    tmp <- DT[[j]]
    if (mode(tmp)==col.mode) set(DT, i=which(is.na(tmp)), j, to)
  }
  return(DT)
}


convert.gene.id <- function(x, from=c("ensembl.gene","ensembl.tx","ensembl.prot","refseq.nm","refseq.np","entrez","uniprot","hgnc","mgi","symbol",
  "affy_hc_g110","affy_hg_focus","affy_hg_u133a","affy_hg_u133a_2","affy_hg_u133b","affy_hg_u133_plus_2","affy_hg_u95a","affy_hg_u95av2","affy_hg_u95b","affy_hg_u95c","affy_hg_u95d","affy_hg_u95e","affy_hta_2_0","affy_huex_1_0_st_v2","affy_hugenefl","affy_hugene_1_0_st_v1","affy_hugene_2_0_st_v1","affy_primeview","affy_u133_x3p","agilent_cgh_44b","agilent_gpl6848","agilent_sureprint_g3_ge_8x60k","agilent_sureprint_g3_ge_8x60k_v2","agilent_wholegenome","agilent_wholegenome_4x44k_v1","agilent_wholegenome_4x44k_v2","codelink_codelink","illumina_humanht_12_v3","illumina_humanht_12_v4","illumina_humanref_8_v3","illumina_humanwg_6_v1","illumina_humanwg_6_v2","illumina_humanwg_6_v3","phalanx_onearray",
  "affy_mg_u74a","affy_mg_u74av2","affy_mg_u74b","affy_mg_u74bv2","affy_mg_u74c","affy_mg_u74cv2","affy_moe430a","affy_moe430b","affy_moex_1_0_st_v1","affy_mogene_1_0_st_v1","affy_mogene_2_1_st_v1","affy_mouse430a_2","affy_mouse430_2","affy_mu11ksuba","affy_mu11ksubb","illumina_mouseref_8","illumina_mousewg_6_v1","illumina_mousewg_6_v2"),
  to=c("symbol","ensembl.gene","ensembl.tx","ensembl.prot","refseq.nm","refseq.np","entrez","uniprot","hgnc","mgi",
  "affy_hc_g110","affy_hg_focus","affy_hg_u133a","affy_hg_u133a_2","affy_hg_u133b","affy_hg_u133_plus_2","affy_hg_u95a","affy_hg_u95av2","affy_hg_u95b","affy_hg_u95c","affy_hg_u95d","affy_hg_u95e","affy_hta_2_0","affy_huex_1_0_st_v2","affy_hugenefl","affy_hugene_1_0_st_v1","affy_hugene_2_0_st_v1","affy_primeview","affy_u133_x3p","agilent_cgh_44b","agilent_gpl6848","agilent_sureprint_g3_ge_8x60k","agilent_sureprint_g3_ge_8x60k_v2","agilent_wholegenome","agilent_wholegenome_4x44k_v1","agilent_wholegenome_4x44k_v2","codelink_codelink","illumina_humanht_12_v3","illumina_humanht_12_v4","illumina_humanref_8_v3","illumina_humanwg_6_v1","illumina_humanwg_6_v2","illumina_humanwg_6_v3","phalanx_onearray",
  "affy_mg_u74a","affy_mg_u74av2","affy_mg_u74b","affy_mg_u74bv2","affy_mg_u74c","affy_mg_u74cv2","affy_moe430a","affy_moe430b","affy_moex_1_0_st_v1","affy_mogene_1_0_st_v1","affy_mogene_2_1_st_v1","affy_mouse430a_2","affy_mouse430_2","affy_mu11ksuba","affy_mu11ksubb","illumina_mouseref_8","illumina_mousewg_6_v1","illumina_mousewg_6_v2"),
  species=c("hs","mm"), use.biomart=TRUE, keep=TRUE) {
  # convert gene ids, default is from some id to gene symbol
  # x: vector of ids
  # use.biomart: for now only TRUE is implemented
  # return a data.table with two columns "from" and "to"; "from" is in the order of x (if keep=TRUE, otherwise from will only contain the items with successful matches and in no particular order); "to" is a vector for 1:1 mapping or a list of 1:many mapping; unmapped items will be NA (if keep=TRUE, otherwise these won't be included)

  #from <- match.arg(from)
  #to <- match.arg(to)
  #species <- match.arg(species)
  from <- from[1] # occasionally I needed to use sth else, but still would like to provide the list of commonly used options for information; thus this way
  to <- to[1]
  species <- species[1]

  if (use.biomart) {
    from <- switch(from,
      ensembl.gene="ensembl_gene_id",
      ensembl.tx="ensembl_transcript_id",
      ensembl.prot="ensembl_peptide_id",
      refseq.nm="refseq_mrna",
      refseq.np="refseq_peptide",
      entrez="entrezgene_id",
      uniprot="uniprotswissprot",
      hgnc="hgnc_id",
      mgi="mgi_id",
      from
    )
    to <- switch(to,
      ensembl.gene="ensembl_gene_id",
      ensembl.tx="ensembl_transcript_id",
      ensembl.prot="ensembl_peptide_id",
      refseq.nm="refseq_mrna",
      refseq.np="refseq_peptide",
      entrez="entrezgene_id",
      uniprot="uniprotswissprot",
      hgnc="hgnc_id",
      mgi="mgi_id",
      to
    )
    if (from=="symbol") from <- switch(species, hs="hgnc_symbol", mm="mgi_symbol")
    if (to=="symbol") to <- switch(species, hs="hgnc_symbol", mm="mgi_symbol")
    db <- switch(species, hs="hsapiens_gene_ensembl", mm="mmusculus_gene_ensembl", species)
    #mart <- biomaRt::useMart(biomart="ensembl", dataset=db)
    mart <- biomaRt::useEnsembl(biomart="genes", dataset=db)
    mapp <- as.data.table(biomaRt::getBM(attributes=c(from, to), filters=from, values=x, mart=mart))
    setnames(mapp, c("from","to"))
  } else {
    stop("Not implemented yet.")
  }

  mapp <- mapp[to!="" & !is.na(to)]
  # if unique map, simplify `to` to a vector, otherwise keep as a list
  n <- mapp[, .(n=uniqueN(to)), by=from][, unique(n)]
  if (length(n)==1 && n==1) {
  	message("Mapping is unique, the `to` column in the returned table is a vector.")
  	if (keep) mapp <- mapp[match(x, from)][, from:=x]
  } else {
  	message("Mapping is not unique, the `to` column in the returned table is a list.")
  	if (keep) {
  	  mapp <- rbind(mapp, data.table(from=unique(setdiff(x, mapp$from)), to=NA))
  	  mapp <- mapp[, .(to=list(unique(to))), by=from][match(x, from)]
  	} else mapp <- mapp[, .(to=list(unique(to))), by=from]	
  }
  mapp
}


convert.gene.id2 <- function(x, from=c("symbol","ensembl.gene","ensembl.tx","ensembl.prot","refseq.nm","refseq.np","entrez","uniprot","hgnc","mgi",
  "affy_hc_g110","affy_hg_focus","affy_hg_u133a","affy_hg_u133a_2","affy_hg_u133b","affy_hg_u133_plus_2","affy_hg_u95a","affy_hg_u95av2","affy_hg_u95b","affy_hg_u95c","affy_hg_u95d","affy_hg_u95e","affy_hta_2_0","affy_huex_1_0_st_v2","affy_hugenefl","affy_hugene_1_0_st_v1","affy_hugene_2_0_st_v1","affy_primeview","affy_u133_x3p","agilent_cgh_44b","agilent_gpl6848","agilent_sureprint_g3_ge_8x60k","agilent_sureprint_g3_ge_8x60k_v2","agilent_wholegenome","agilent_wholegenome_4x44k_v1","agilent_wholegenome_4x44k_v2","codelink_codelink","illumina_humanht_12_v3","illumina_humanht_12_v4","illumina_humanref_8_v3","illumina_humanwg_6_v1","illumina_humanwg_6_v2","illumina_humanwg_6_v3","phalanx_onearray",
  "affy_mg_u74a","affy_mg_u74av2","affy_mg_u74b","affy_mg_u74bv2","affy_mg_u74c","affy_mg_u74cv2","affy_moe430a","affy_moe430b","affy_moex_1_0_st_v1","affy_mogene_1_0_st_v1","affy_mogene_2_1_st_v1","affy_mouse430a_2","affy_mouse430_2","affy_mu11ksuba","affy_mu11ksubb","illumina_mouseref_8","illumina_mousewg_6_v1","illumina_mousewg_6_v2"),
  to=c("symbol","ensembl.gene","ensembl.tx","ensembl.prot","refseq.nm","refseq.np","entrez","uniprot","hgnc","mgi",
  "affy_hc_g110","affy_hg_focus","affy_hg_u133a","affy_hg_u133a_2","affy_hg_u133b","affy_hg_u133_plus_2","affy_hg_u95a","affy_hg_u95av2","affy_hg_u95b","affy_hg_u95c","affy_hg_u95d","affy_hg_u95e","affy_hta_2_0","affy_huex_1_0_st_v2","affy_hugenefl","affy_hugene_1_0_st_v1","affy_hugene_2_0_st_v1","affy_primeview","affy_u133_x3p","agilent_cgh_44b","agilent_gpl6848","agilent_sureprint_g3_ge_8x60k","agilent_sureprint_g3_ge_8x60k_v2","agilent_wholegenome","agilent_wholegenome_4x44k_v1","agilent_wholegenome_4x44k_v2","codelink_codelink","illumina_humanht_12_v3","illumina_humanht_12_v4","illumina_humanref_8_v3","illumina_humanwg_6_v1","illumina_humanwg_6_v2","illumina_humanwg_6_v3","phalanx_onearray",
  "affy_mg_u74a","affy_mg_u74av2","affy_mg_u74b","affy_mg_u74bv2","affy_mg_u74c","affy_mg_u74cv2","affy_moe430a","affy_moe430b","affy_moex_1_0_st_v1","affy_mogene_1_0_st_v1","affy_mogene_2_1_st_v1","affy_mouse430a_2","affy_mouse430_2","affy_mu11ksuba","affy_mu11ksubb","illumina_mouseref_8","illumina_mousewg_6_v1","illumina_mousewg_6_v2"),
  from.sp=c("mm","hs"), to.sp=c("hs","mm"), keep=TRUE) {
  # convert gene ids across species, for now only between human and mice; default is mice gene symbol to human gene symbol
  # x: vector of ids
  # return a data.table with two columns "from" and "to"; "from" is in the order of x (if keep=TRUE, otherwise from will only contain the items with successful matches and in no particular order); "to" is a vector for 1:1 mapping or a list of 1:many mapping; unmapped items will be NA (if keep=TRUE, otherwise these won't be included)

  #from <- match.arg(from)
  #to <- match.arg(to)
  #from.sp <- match.arg(from.sp)
  #to.sp <- match.arg(to.sp)
  from <- from[1] # occasionally I needed to use sth else, but still would like to provide the list of commonly used options for information; thus this way
  to <- to[1]
  from.sp <- from.sp[1]
  to.sp <- to.sp[1]

  from <- switch(from,
    ensembl.gene="ensembl_gene_id",
    ensembl.tx="ensembl_transcript_id",
    ensembl.prot="ensembl_peptide_id",
    refseq.nm="refseq_mrna",
    refseq.np="refseq_peptide",
    entrez="entrezgene_id",
    uniprot="uniprotswissprot",
    hgnc="hgnc_id",
    mgi="mgi_id",
    from
  )
  to <- switch(to,
    ensembl.gene="ensembl_gene_id",
    ensembl.tx="ensembl_transcript_id",
    ensembl.prot="ensembl_peptide_id",
    refseq.nm="refseq_mrna",
    refseq.np="refseq_peptide",
    entrez="entrezgene_id",
    uniprot="uniprotswissprot",
    hgnc="hgnc_id",
    mgi="mgi_id",
    to
  )
  if (from=="symbol") from <- switch(from.sp, hs="hgnc_symbol", mm="mgi_symbol")
  if (to=="symbol") to <- switch(to.sp, hs="hgnc_symbol", mm="mgi_symbol")
  from.db <- switch(from.sp, hs="hsapiens_gene_ensembl", mm="mmusculus_gene_ensembl", from.sp)
  to.db <- switch(to.sp, hs="hsapiens_gene_ensembl", mm="mmusculus_gene_ensembl", to.sp)
  #from.mart <- biomaRt::useMart(biomart="ensembl", dataset=from.db)
  from.mart <- biomaRt::useEnsembl(biomart="genes", dataset=from.db)
  #to.mart <- biomaRt::useMart(biomart="ensembl", dataset=to.db)
  to.mart <- biomaRt::useEnsembl(biomart="genes", dataset=to.db)
  mapp <- as.data.table(biomaRt::getLDS(attributes=from, filters=from, values=x, mart=from.mart, attributesL=to, martL=to.mart))
  setnames(mapp, c("from","to"))

  mapp <- mapp[to!="" & !is.na(to)]
  # if unique map, simplify `to` to a vector, otherwise keep as a list
  n <- mapp[, .(n=uniqueN(to)), by=from][, unique(n)]
  if (length(n)==1 && n==1) {
  	message("Mapping is unique, the `to` column in the returned table is a vector.")
  	if (keep) mapp <- mapp[match(x, from)][, from:=x]
  } else {
  	message("Mapping is not unique, the `to` column in the returned table is a list.")
  	if (keep) {
  	  mapp <- rbind(mapp, data.table(from=unique(setdiff(x, mapp$from)), to=NA))
  	  mapp <- mapp[, .(to=list(unique(to))), by=from][match(x, from)]
  	} else mapp <- mapp[, .(to=list(unique(to))), by=from]	
  }
  mapp
}


convert.gset.species <- function(gs, from="hs", to="mm") {
  # convert gene sets containing gene symbols from on species to another, by default from human to mice
  # gs: gene sets in a list

  gns <- unique(unlist(gs))
  from1 <- from
  to1 <- to
  mapp <- suppressMessages(convert.gene.id2(gns, from.sp=from1, to.sp=to1))
  lapply(gs, function(x) {
    res <- mapp[from %in% x, unique(unlist(to))]
    res[!is.na(res)]
  })
}


qlapply <- function(f, ..., pkgs=c(), njobs, mem=16, config=list(P="short"), fail_on_error=FALSE) {
  # paralleled lapply of function f with clustermq::Q
  # first item in ... should be the variable to iterate over, named by the corresponding argument name of f; the names of the returned list will be equal to names(list(...)[[1]])
  # other items in ... should be either other constant arguments passed to f or variables to be exported to f (i.e. those not defined within f)
  # mem is memory per job in G; could also specify in config as list(mem=xx, ...), in which case the mem argument will be ignored
  # the names "P" and "mem" used in config correspond to my current template in use (data/.clustermq.sge.tmpl)

  l <- list(...)
  x <- l[1]
  l <- l[-1]
  fargs <- formalArgs(f)
  const <- l[names(l) %in% fargs]
  export <- l[!names(l) %in% fargs]
  if (missing(njobs)) njobs <- min(length(x[[1]]), 100)
  if (!"mem" %in% names(config)) config$mem <- mem # the name "mem" used in config correspond to my current template in use (data/.clustermq.sge.tmpl)
  args <- c(x, list(fun=f, const=const, export=export, pkgs=pkgs, n_jobs=njobs, template=config, fail_on_error=fail_on_error))
  res <- do.call(clustermq::Q, args=args)
  names(res) <- names(x[[1]])
  res
}


apply1 <- function(dat, pheno, f, ..., var.name="x", nc=1L) {
  # dat is a gene(feature)-by-sample matrix; pheno is a data.table of sample-level covariates; this function will apply a function f (one of the run.* functions) for certain statistical testing to each gene/feature
  # each gene/feature will be cbind-ed with pheno with column name var.name, and the resulting data will be passed to the "dat" argument of f; provide additional arguments to f in ...
  # this function will rbind all f outputs, adding an "id" column of gene/feature name (using the rownames of dat or dat$expr), add BH-adjusted p value with adjust.pval() then order by padj as output
  # nc: number of cores

  args <- list(...)
  args$dat <- copy(pheno)

  cl <- makeCluster(nc, type="FORK")
  res <- rbindlist(parApply(cl, dat, 1, function(x) {
    args$dat[[var.name]] <- x
    do.call(f, args)
  }), idcol="id")
  stopCluster(cl)

  res <- adjust.pval(res)
  res[order(padj, pval)]
}


