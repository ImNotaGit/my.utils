## ----some utility functions for common small tasks----


cn <- function(...) {
  res <- c(...)
  names(res) <- res
  res
}


hh <- function(x, nr=5, nc=nr) {
  # check the first nr rows and nc columns of a matrix-like object
  if (nr>nrow(x)) nr <- nrow(x)
  if (nc>ncol(x)) nc <- ncol(x)
  x[1:nr, 1:nc]
}


write.tab <- function(x, file, append=FALSE, col.names=TRUE) {
  # write a tsv file w/o quotation and w/o row names, but by default with col names
  fwrite(x, file=file, append=append, quote=FALSE, sep="\t", na="NA", col.names=col.names)
}


regex.escape <- function(x) {
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

