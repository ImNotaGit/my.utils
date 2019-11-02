## ----functions for working with multi-type cancer data----


parse.cancer.type <- function(type) {
  # type should be an atom vector, and it will be parsed into a vector of TCGA cancer type names as below
  if (length(type)!=1) stop("In parse.cancer.type: the length of type should be 1.\n")
  all.tcga.types <- c("ACC", "BLCA", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "PAAD", "PCPG", "READ", "SARC", "SKCM", "STAD", "THCA", "THYM", "UVM", "BRCA", "CESC", "OV", "UCEC", "UCS", "PRAD", "TGCT")
  res <- switch(tolower(type),
                adrenocortical =,
                adrenal = "ACC",
                `urinary bladder` = "BLCA",
                breast = "BRCA",
                cervical =,
                endocervical = "CESC",
                `gall bladder` =,
                cholangia =,
                cholangial = "CHOL",
                colon = "COAD",
                coadread =,
                crc =,
                colorectum =,
                colorectal = c("COAD", "READ"),
                lymphoma = "DLBC",
                esophagus =,
                esophageal = "ESCA",
                glioblastoma = "GBM",
                gbmlgg =,
                brain =,
                glioma = c("GBM", "LGG"),
                `head and neck` = "HNSC",
                kipan =,
                kidney = c("KICH", "KIRC", "KIRP"),
                leukemia =,
                blood =,
                aml = "LAML",
                liver =,
                hcc =,
                hepatocellular = "LIHC",
                lung = c("LUAD", "LUSC"),
                mesothelioma = "MESO",
                ovary =,
                ovarian = "OV",
                pancreas =,
                pancreatic = "PAAD",
                pheochromocytoma =,
                paraganglioma = "PCPG",
                nerve =,
                nervous =,
                `nervous system` = c("GBM", "LGG", "PCPG"),
                pc =,
                prostate = "PRAD",
                rectum =,
                rectal = "READ",
                sarcoma = "SARC",
                skin =,
                melanoma = c("SKCM", "UVM"),
                stomach = "STAD",
                stes =,
                `stomach and esophageal` = c("ESCA", "STAD"),
                testicle =,
                testicular = "TGCT",
                thyroid = "THCA",
                thymus =,
                thymoma = "THYM",
                uterus =,
                uterine = c("UCEC", "UCS"),
                m =,
                male = c("PRAD", "TGCT"),
                f =,
                female = c("BRCA", "CESC", "OV", "UCEC", "UCS"),
                mf =,
                monosex = c("BRCA", "CESC", "OV", "UCEC", "UCS", "PRAD", "TGCT"),
                bothsex =,
                unisex = c("ACC", "BLCA", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "PAAD", "PCPG", "READ", "SARC", "SKCM", "STAD", "THCA", "THYM", "UVM"),
                all =,
                pan =,
                pancancer = all.tcga.types,
                toupper(type)
  )
  if (length(res)==1) {
    if (!res %in% all.tcga.types) stop("In parse.cancer.type: invalid type.\n")
  }
  return(res)
}


select.cancers <- function(type) {

  # return a named vector of selected types of cancers. type="male" for cancers only in males; similarly, "female" (note: BRCA is regarded as a female cancer); "mf" for male+female; "unisex" for cancers that can happen in both sexes; "simple" for male+female+unisex; "aggregate" for aggregate data of multiple types; "unisex.w.aggr" for unisex+aggregate; "all" for simple+aggregate; or give a vector of specific cancer code(s).

  # different categories of cancers
  cancer.unisex <- c("ACC", "BLCA", "CHOL", "COAD", "DLBC", "ESCA", "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", "LUAD", "LUSC", "MESO", "PAAD", "PCPG", "READ", "SARC", "SKCM", "STAD", "THCA", "THYM", "UVM")
  cancer.male <- c("PRAD", "TGCT")
  cancer.female <- c("BRCA", "CESC", "OV", "UCEC", "UCS")
  cancer.aggregate.nopan <- c("COADREAD", "GBMLGG", "KIPAN", "STES")
  cancer.aggregate <- c("COADREAD", "GBMLGG", "KIPAN", "STES", "unisex", "PAN")
  all.valid.types <- c(cancer.unisex, cancer.male, cancer.female, cancer.aggregate)

  # if length(type)>1, we assume it's provided as a vector of simple cancer type names, so return itself
  if (length(type)>1) {
    if (!all(type %in% all.valid.types)) stop("In select.cancers: invalid type.\n")
    names(type) <- type
    return(type)
  }

  # switch only works for variables of length==1
  res <- switch(type,
                male = cancer.male,
                female = cancer.female,
                mf = c(cancer.male, cancer.female),
                unisex = cancer.unisex,
                simple = c(cancer.male, cancer.female, cancer.unisex),
                aggregate = cancer.aggregate,
                unisex.w.aggr = c(cancer.unisex, cancer.aggregate),
                all.but.pan = c(cancer.male, cancer.female, cancer.unisex, cancer.aggregate.nopan),
                all = c(cancer.male, cancer.female, cancer.unisex, cancer.aggregate),
                type)

  if (!(res %in% all.valid.types)) stop("In select.cancers: invalid type.\n")
  names(res) <- res
  return(res)

}


cancers.lapply <- function(cancers, f, ...) {
  # lapply() f() on cancers with try-catch, then remove empty (NA or NULL) elements
  res <- lapply(select.cancers(cancers),
                function(cancer) {
                  tryout <- try(f(cancer, ...))
                  if (class(tryout)=="try-error") {
                    cat("Error while processing cancer type: ", cancer, "\n")
                    return(NULL)
                  }
                  return(tryout)
                }
  )
  res[!sapply(res, function(x) identical(x, NA) || identical(x, NULL))]
}

