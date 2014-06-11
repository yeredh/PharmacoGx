########################
## Benjamin Haibe-Kains
## All rights Reserved
## October 23, 2013
########################

#################################################
## Get fRMA normalized CGP data from InSilicoDB
##
## 
#################################################


###############
## TODO cell_id ==?? cellid???
###############




`getCGP` <- 
function (gene=TRUE, tmpdir="tmp", delete.tmpdir=FALSE, cosmic.annotation=FALSE, cosmic.version="v68", 
  replicates=c("last", "first", "all", "mean", "median"), verbose=FALSE, downloadMethod="wget") {

  replicates <- match.arg(replicates)
  
  ## create directories for temporary files
  if(!file.exists(tmpdir)) { dir.create(tmpdir, showWarnings=FALSE, recursive=TRUE) }
  
  badchars <- "[\xb5]|[]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
  
  if (verbose) { message("Downloading the genomic data of the Cancer Genome Project from InSilicoDB") }
  inSilicoDb::InSilicoLogin(login="bhaibeka@gmail.com", password="747779bec8a754b91076d6cc1f700831")
  # inSilicoDb::getCurationInfo(dataset="ISDB12210")
  platf <- inSilicoDb::getPlatforms(dataset="ISDB12210")
  eset <- inSilicoDb::getDatasets(dataset="ISDB12210", norm="FRMA", curation="24802", features="PROBE")
  inSilicoDb::InSilicoLogout()
  
  ## only one platform, may be subject to change
  platf <- platf[[1]]
  eset <- eset[[1]]
  
  colnames(Biobase::pData(eset)) <- gsub(badchars, "_", colnames(Biobase::pData(eset)))
  
  ## gene centric expression
  if (gene) {
    if (verbose) { message("Gene centric data") }
    eset <- MetaGx::probeGeneMapping(eset=eset, platform="GPL96", method="jetset")
  }
  
  ## replicated experiments
  pheno <- Biobase::pData(eset)
  switch(replicates,
    "first" = {
      if (verbose) { message("First experiment of each replicate is kept") }
      iix <- order(pheno[ , "file_day"], pheno[ , "file_hour"], decreasing=FALSE, na.last=TRUE)
      ix <- rownames(pheno)[iix][!duplicated(pheno[iix, "cell_id"])]
      Biobase::exprs(eset) <- Biobase::exprs(eset)[ , ix, drop=FALSE]
      Biobase::pData(eset) <- Biobase::pData(eset)[ix, , drop=FALSE]
    },
    "last" = {
      if (verbose) { message("Last experiment of each replicate is kept") }
      iix <- order(pheno[ , "file_day"], pheno[ , "file_hour"], decreasing=TRUE, na.last=TRUE)
      ix <- rownames(pheno)[iix][!duplicated(pheno[iix, "cell_id"])]
      Biobase::exprs(eset) <- Biobase::exprs(eset)[ , ix, drop=FALSE]
      Biobase::pData(eset) <- Biobase::pData(eset)[ix, , drop=FALSE]
    },
    "mean" = {
      if (verbose) { message("Mean of replicates is computed") }
      stop(sprintf("Method replicates %s not implemented yet", replicates))
    },
    "median" = {
      if (verbose) { message("Median of replicates is computed") }
      stop(sprintf("Method replicates %s not implemented yet", replicates))
      ix <- pheno[duplicated(pheno[ , "cell_id"]), "cell_id"]
      nn <- NULL
      ex <- matrix(NA, nrow=nrow(Biobase::fData(eset)), ncol=sum(!duplicated(pheno[ , "cell_id"]), na.rm=TRUE), dimnames=list(rownames(Biobase::fData(eset)), nn))
      for (i in 1:length(ix)) {
        iix <- !is.na(pheno[ , "cell_id"]) & pheno[ , "cell_id"] == ix[i]
        xx <- apply(Biobase::exprs(eset)[ , iix], 1, get(replicates), na.rm=TRUE)
      }
    },
    "all" = {
      if (verbose) { message("All replicates are kept in the expression set") }
    }  
  )
  
  ## replace sample names by cell line names
  nn <- genefu::rename.duplicate(Biobase::pData(eset)[ , "cell_id"], sep="_")$new.x
  rownames(Biobase::pData(eset)) <- colnames(Biobase::exprs(eset)) <- nn
  
  
  ## ftp for new cgp data
  ftpdir <- "ftp://ftp.ebi.ac.uk//pub/databases/microarray/data/experiment/MTAB/E-MTAB-783/"
  
  
  
  ## get drug sensivity data
  if (verbose) { message("Download and format drug sensitivity data") }
  tmpfiles <- NULL

  ## download drug sensitivity (release 2)
  myfn <- file.path(tmpdir, "cgp_drug_sensitivity.csv")
  if(!file.exists(myfn)) {
    if (verbose) { message("Download drug sensitivity measurements") }
    dwl.status <- download.file(url="ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-5.0/gdsc_manova_input_w5.csv", destfile=myfn, method=downloadMethod)
    if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
    tmpfiles <- c(tmpfiles, myfn)
  }

  ## download drug concentration (release 2)
  myfn <- file.path(tmpdir, "cgp_drug_concentration.csv")
  if(!file.exists(myfn)) {
    if (verbose) { message("Download screening drug concentrations") }
    dwl.status <- download.file(url="ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-5.0/gdsc_compounds_conc_w5.csv", destfile=myfn, method=downloadMethod)
    if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
    tmpfiles <- c(tmpfiles, myfn)
  }

  ## download cell line annotations and COSMIC IDs
  myfn <- file.path(tmpdir, "celline_annotations.RData")
  if(!file.exists(myfn)) {
    if (verbose) { message("Download cell lines annotations") }
    ## annotations from GDSC (Genomics of Drug Sensitivity in Cancer)
    dwl.status <- download.file(url="ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/release-5.0/gdsc_cell_lines_w5.csv", destfile=file.path(tmpdir, "cgp_celline_collection.csv"), method=downloadMethod)
    tmpfiles <- c(tmpfiles, file.path(tmpdir, "cgp_celline_collection.csv"))
    if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }
    celline.gdsc <- read.csv(file=file.path(tmpdir, "cgp_celline_collection.csv"), stringsAsFactors=FALSE)
    celline.gdsc[celline.gdsc == "" | celline.gdsc == " " | celline.gdsc == "  "] <- NA
    celline.gdsc <- celline.gdsc[!is.na(celline.gdsc[ , "CELL_LINE_NAME"]), , drop=FALSE]
    dupln <- unique(celline.gdsc[ , "CELL_LINE_NAME"][duplicated(celline.gdsc[ , "CELL_LINE_NAME"])])
    celline.gdsc <- celline.gdsc[!duplicated(celline.gdsc[ , "CELL_LINE_NAME"]), , drop=FALSE]
    celline.gdsc[ , "COSMIC_ID"] <- as.character(celline.gdsc[ , "COSMIC_ID"])
    rownames(celline.gdsc) <- celline.gdsc[ , "CELL_LINE_NAME"]
    ## annotations from COSMIC
    if (cosmic.annotation) {
      
      #       dwl.status <- download.file(url=sprintf("ftp://ftp.sanger.ac.uk/pub/CGP/cell_lines_project/data_export/CosmicCellLineProject_%s.tsv.gz", cosmic.version), destfile=file.path(tmpdir, sprintf("CosmicCellLineProject_%s.tsv.gz", cosmic.version)), method=downloadMethod)
      
      dwl.status <- getCosmic(em="bhk.labgroup@gmail.com", passw="pharmacogenomics", directory=file.path(tmpdir, sprintf("CosmicCellLineProject_%s.tsv.gz", cosmic.version)))
      if(dwl.status != 0) { stop("Download failed, please rerun the pipeline! It may be that there is a new version of the file CosmicCellLineProject, please look at ftp://ftp.sanger.ac.uk/pub/CGP/cosmic/data_export/ and update the script accordingly ...") }
      ## untar
      res <- R.utils::gunzip(filename=file.path(tmpdir, sprintf("CosmicCellLineProject_%s.tsv.gz", cosmic.version)), overwrite=TRUE)
      tmpfiles <- c(tmpfiles, file.path(tmpdir, sprintf("CosmicCellLineProject_%s.tsv", cosmic.version)))
      celline.cosmic <- read.csv(file=file.path(tmpdir, sprintf("CosmicCellLineProject_%s.tsv", cosmic.version)), sep="\t", stringsAsFactors=FALSE)
      cosmic.celline <- cosmic.celline[complete.cases(cosmic.celline[ , c("Sample.name", "Sample.source")]) & cosmic.celline[ , "Sample.source"] == "cell-line", , drop=FALSE]
      cosmic.celline[cosmic.celline == "NS" | cosmic.celline == "" | cosmic.celline == " " | cosmic.celline == "  "] <- NA
      ## remove cell line with no name
      celline.cosmic <- celline.cosmic[!is.na(celline.cosmic[ , "Sample.name"]), , drop=FALSE]
      ## merge the gene targets
      dupln <- unique(celline.cosmic[ , "Sample.name"][duplicated(celline.cosmic[ , "Sample.name"])])
      tt <- celline.cosmic
      iix.rm <- NULL
      for(i in 1:length(dupln)) {
        duplix <- celline.cosmic[ ,"Sample.name"] == dupln[i]
        iix <- sort((which(duplix)), decreasing=FALSE)[1]
        iix.rm <- c(iix.rm, setdiff(which(duplix), iix))
        tt[iix, "Gene.name"] <- paste(celline.cosmic[duplix, "Gene.name"], collapse="///")
        tt[iix, "UniProt.ID"] <- paste(celline.cosmic[duplix, "UniProt.ID"], collapse="///")
        tt[iix, "Zygosity"] <- paste(celline.cosmic[duplix, "Zygosity"], collapse="///")
        tt[iix, "CDS_MUT_SYNTAX"] <- paste(celline.cosmic[duplix, "CDS_MUT_SYNTAX"], collapse="///")
        tt[iix, "AA_MUT_SYNTAX"] <- paste(celline.cosmic[duplix, "AA_MUT_SYNTAX"], collapse="///")
        tt[iix, "NCBI36.genome.position"] <- paste(celline.cosmic[duplix, "NCBI36.genome.position"], collapse="///")
        tt[iix, "GRCh37.genome.position"] <- paste(celline.cosmic[duplix, "GRCh37.genome.position"], collapse="///")
      }
      tt <- tt[-iix.rm, , drop=FALSE]
      rownames(tt) <- tt[ , "Sample.name"]
      celline.cosmic <- tt
    } else {
      celline.cosmic <- cbind("Sample.name"=celline.gdsc[ , "CELL_LINE_NAME"], "ID_sample"=celline.gdsc[ , "COSMIC_ID"])
      rownames(celline.cosmic) <- celline.cosmic[ , "Sample.name"]
    }
    ## merge GDSC and COSMIC annotations through COSMIC_ID
    # iix <- which(!is.na(celline.gdsc[ , "COSMIC_ID"]) & !is.element(celline.gdsc[ , "COSMIC_ID"], celline.cosmic[ , "ID_sample"]))
    iix <- which(complete.cases(celline.gdsc[ , c("CELL_LINE_NAME", "COSMIC_ID")]) & !is.element(celline.gdsc[ , "COSMIC_ID"], celline.cosmic[ , "ID_sample"]) & !is.element(celline.gdsc[ , "CELL_LINE_NAME"], celline.cosmic[ , "Sample.name"]))
    if (length(iix) == 0) {
      tt <- data.frame(matrix(NA, nrow=nrow(celline.cosmic) + length(iix), ncol=ncol(celline.cosmic), dimnames=list(c(rownames(celline.cosmic), rownames(celline.gdsc)[iix]), colnames(celline.cosmic))))
      tt[rownames(celline.cosmic), ] <- celline.cosmic
      tt[rownames(celline.gdsc)[iix], "Sample.name"] <- celline.gdsc[iix, "CELL_LINE_NAME"]
      tt[rownames(celline.gdsc)[iix], "ID_sample"] <- celline.gdsc[iix, "COSMIC_ID"]
      celline <- tt
      colnames(celline)[match(c("Sample.name", "ID_sample"), colnames(celline))] <- c("CELL_LINE_NAME", "COSMIC_ID")
    } else {
      if (cosmic.annotation) {
        celline <- celline.cosmic
        colnames(celline)[match(c("Sample.name", "ID_sample"), colnames(celline))] <- c("CELL_LINE_NAME", "COSMIC_ID")
      } else {
        celline <- celline.gdsc
      }
    }
    save(list=c("celline.cosmic", "celline.gdsc", "celline"), compress=TRUE, file=myfn)
  } else {
    load(myfn)
  }

  ## download drug information
  message("Download drug information")
  myfn <- file.path(tmpdir, "cgp_drug_information.csv")
  if (!file.exists(myfn)) {
    # dwl.status <- download.file(url="http://www.cancerrxgene.org/action/ExportJsonTable/CSV", destfile=file.path(path.drug, "dwl", "export-Automatically_generated_table_data.csv"), quiet=TRUE)
    # if(dwl.status != 0) { stop("Download failed, please rerun the pipeline!") }  
    tables <- XML::readHTMLTable("http://www.cancerrxgene.org/translation/Drug")
    drugs <- tables[1][[1]]
    write.csv(drugs, row.names=FALSE, file=myfn)
  }
  myfn <- file.path(tmpdir, "nature11005-s2.zip")
  if(!file.exists(myfn)) {
    if (verbose) { message("Download nature supplementary information") }
    dwl.status <- download.file(url="http://www.nature.com/nature/journal/v483/n7391/extref/nature11005-s2.zip", destfile=myfn, method=downloadMethod)
    if(dwl.status != 0) { stop("Download failed, please rerun the script!") }
  }
  ff <- as.character(unzip(zipfile=file.path(tmpdir, "nature11005-s2.zip"), list=TRUE)[1, 1])
  unzip(zipfile=file.path(tmpdir, "nature11005-s2.zip"), exdir=file.path(tmpdir))
  tmpfiles <- c(tmpfiles, file.path(tmpdir, "nature_supplementary_information.xls"))
  
  ## phenotype for the drugs
  if (verbose) { message("Read drug sensitivity measurements") }
  myfn <- file.path(tmpdir, "cgp_drug_sensitivity.RData")
  if(!file.exists(myfn)) {
    drugpheno <- read.csv(file.path(tmpdir, "cgp_drug_sensitivity.csv"), stringsAsFactors=FALSE)
    drugpheno[drugpheno == "" | drugpheno == " "] <- NA
    save(list="drugpheno", compress=TRUE, file=myfn)
    tmpfiles <- c(tmpfiles, myfn)
  } else { load(myfn) }
  ## format column names
  coln2 <- unlist(drugpheno[1, ,drop=TRUE])
  coln2[coln2 == ""] <- NA
  #drugpheno <- drugpheno[-1, ,drop=FALSE]
  drugpheno <- drugpheno[!is.na(drugpheno[ , "Cell.Line"]), ,drop=FALSE]
  coln <- colnames(drugpheno)
  coln2[is.na(coln2)] <- coln[is.na(coln2)]
  coln2 <- genefu::rename.duplicate(x=coln2, sep="_dupl")$new.x
  myx <- sapply(sapply(strsplit(coln2, "_"), function(x) { return(x[[1]]) }), Hmisc::all.is.numeric)
  coln2[myx] <- paste("drugid", gsub(pattern=badchars, replacement="_", x=toupper(coln2[myx])), sep="_")
  colnames(drugpheno) <- coln2
  ## drug identifiers and names
  dn <- toupper(gsub(badchars, "", sapply(strsplit(coln, "_"), function(x) { return(x[[1]]) })))
  ## manual curation for drug names starting with a figure
  dn[!is.na(dn) & dn == "X17AAG"] <- "17AAG"
  dn[!is.na(dn) & dn == "X681640"] <- "681640"
  did <- sapply(strsplit(coln2, "_"), function(x) { if(x[[1]] == "drugid") { return(x[[2]]) } else { return(NA) } })
  drugnid <- cbind("drug.name"=dn, "drug_id"=did)[!is.na(did) & !duplicated(did), ]
  rownames(drugnid) <- paste("drugid", drugnid[ , "drug_id"], sep="_")

  ## cell line identifiers
  dupln <- duplicated(drugpheno[ ,"Cell.Line"])
  if(sum(dupln) > 1) { warning("some cell lines are duplicated, only the first instance is kept") }
  drugpheno <- drugpheno[!dupln, , drop=FALSE]
  
  if(any(!is.element(drugpheno[ ,"Cell.Line"], celline[ , "CELL_LINE_NAME"]))) { warning("Some cell line with drug sensitivity data have no annotations") }
  celln <- as.character(drugpheno[ ,"Cell.Line"])
  drugpheno <- data.frame("cell_id"=celln, drugpheno, stringsAsFactors=FALSE)
  rownames(drugpheno) <- celln

  ## get mutational data, i.e., protein coding variants
  ## Genetic mutation data for cancer genes. Includes MSI status (1 = unstable and 0 = stable) and gene-fusions. A binary code 'x::y' description is used for each gene where 'x' identifies a coding variant and 'y' indicates copy number information from SNP6.0 data. For gene fusions, cell lines are identified as fusion not-detected (0) or the identified fusion is given. The following abbreviations are used: not analysed (na), not detected or wild-type (wt), no copy number information (nci).
  ## we assume that AKT2 and WT1 are the first and last genes in the file
  #rangeg <- which(colnames(drugpheno) == "AKT2"):which(colnames(drugpheno) == "WT1")
  
  rangeg <- which(colnames(drugpheno) == "AKT2"):which(colnames(drugpheno) == "MLL_AFF1")
  
  mutation <- as.matrix(drugpheno[ , rangeg, drop=FALSE])
  mutation <- apply(X=mutation, MARGIN=c(1, 2), FUN=function(x) {
    x <- unlist(strsplit(x, split="::"))
    if(length(x) == 2) {
      ####
      # TODO Ask benjamin
      # What is going on here?
      ####
      if(!is.na(x[[1]]) && (x[[1]] == "na" || x[[1]] == "p.?" || x[[1]] == "p.0?")) {
        x <- NA
      } else {
        x <- x[[1]]
      }
    } else { x <- NA }
    return(x)
  })
  
  ## url for cell line collection
  celline <- data.frame("cell_id"=as.character(celline[ , "CELL_LINE_NAME"]), celline, stringsAsFactors=FALSE)
  ## add url based on COSMIC IDs
  uurl <- paste("http://cancer.sanger.ac.uk/cell_lines/sample/overview?id=", celline[ , "COSMIC_ID"], sep="")
  uurl[is.na(celline[ , "COSMIC_ID"])] <- NA
  celline <- data.frame("cell_id"=celline[ , "cell_id"], "link"=uurl, celline[ , !is.element(colnames(celline), "cell_id")], stringsAsFactors=FALSE)

  ## make sure that cell_id are not factors
  celline[, "cell_id"] <- as.character(celline[, "cell_id"])
  Biobase::pData(eset)[, "cell_id"] <- as.character(Biobase::pData(eset)[, "cell_id"])
  drugpheno[, "cell_id"] <- as.character(drugpheno[, "cell_id"])
  
  #####
  ## TODO What is going on here?
  #####
  
  ## union of all cell line with data
  cellnall <- sort(unique(c(row.names(Biobase::pData(eset)), rownames(mutation), as.character(drugpheno[ ,"cell_id"]))))
  
  ## update sampleinfo
  dd <- data.frame(matrix(NA, ncol=ncol(Biobase::pData(eset)), nrow=length(cellnall), dimnames=list(cellnall, colnames(Biobase::pData(eset)))), check.names=FALSE)
  dd[rownames(Biobase::pData(eset)), colnames(Biobase::pData(eset))] <- Biobase::pData(eset)
  Biobase::pData(eset) <- dd
  
  ## update gene expressions
  dd <- matrix(NA, nrow=nrow(Biobase::exprs(eset)), ncol=length(cellnall), dimnames=list(rownames(Biobase::exprs(eset)), cellnall))
  dd[rownames(Biobase::exprs(eset)), colnames(Biobase::exprs(eset))] <- Biobase::exprs(eset)
  Biobase::exprs(eset) <- dd
  
  ## update drugpheno
  dd <- data.frame(matrix(NA, ncol=ncol(drugpheno), nrow=length(cellnall)))
  rownames(dd) <- cellnall
  colnames(dd) <- colnames(drugpheno)
  newlev <- sapply(drugpheno, levels)
  newlev$cell_id <- cellnall
  ## genefu::setcolclass.df is NOT in the Bioconductor package, but in the package on github (https://github.com/bhaibeka/genefu/blob/master/R/setcolclass.df.R)
  dd <- genefu::setcolclass.df(df=dd, colclass=sapply(drugpheno, class), factor.levels=newlev)
  dd[rownames(drugpheno),colnames(drugpheno)] <- drugpheno
  dd[ ,"cell_id"] <- cellnall
  drugpheno <- dd

  ## update mutation
  dd <- matrix(NA, ncol=ncol(mutation), nrow=length(cellnall), dimnames=list(cellnall, colnames(mutation)))
  dd[rownames(mutation), colnames(mutation)] <- mutation
  mutation <- t(dd)
  
  ## update celline
  dd <- data.frame(matrix(NA, ncol=ncol(celline), nrow=length(cellnall), dimnames=list(cellnall, colnames(celline))), check.names=FALSE)
  iix <- intersect(rownames(celline), cellnall)
  dd[iix, colnames(celline)] <- celline[iix, , drop=FALSE]
  celline <- dd
  celline[ , "cell_id"] <- celline[ , "CELL_LINE_NAME"] <- rownames(celline)
  ## annotate cell lines with curated tissue type
  tissue.type <- read.csv(file.path(system.file("extdata", package="PharmacoGx"), "cell_line_collection_all.csv"), stringsAsFactors=FALSE)
  rownames(tissue.type) <- tissue.type[ , 1]
  celline <- cbind("tissue.type"=tissue.type[match(celline[ , "cell_id"], tissue.type[ , "cell_id"]), "tissue.type"], celline)


  # ## reproducibility between different screening sites
 #    ## camptothecin was screened at MGH (drug id 195) and WTSI (drug id 1003)
 #    ## data only available in the supplementary infomration of the Nature website
 #    myfn2 <- file.path(saveres, "nature_supplinfo_drugpheno_cgp.RData")
 #    if(!file.exists(myfn2)) {
 #      drugpheno.nature <- gdata::read.xls(xls=file.path(path.drug, "nature_supplementary_information.xls"), sheet=2)
 #      drugpheno.nature[drugpheno.nature == "" | drugpheno.nature == " "] <- NA
 #      save(list="drugpheno.nature", compress=TRUE, file=myfn2)
 #    } else { load(myfn2) }
 #    ## format column names
 #    coln2 <- gsub(" ", "", sapply(drugpheno.nature[1,], as.character))
 #    coln2[coln2 == ""] <- NA
 #    drugpheno.nature <- drugpheno.nature[-1, ,drop=FALSE]
 #    coln <- colnames(drugpheno.nature)
 #    coln2[is.na(coln2)] <- coln[is.na(coln2)]
 #    coln2 <- genefu::rename.duplicate(x=coln2, sep="_dupl")$new.x
 #    myx <- sapply(sapply(strsplit(coln2, "_"), function(x) { return(x[[1]]) }), Hmisc::all.is.numeric)
 #    coln2[myx] <- paste("drugid", gsub(pattern=badchars, replacement="_", x=toupper(coln2[myx])), sep="_")
 #    colnames(drugpheno.nature) <- coln2
 #    myx <- sapply(strsplit(colnames(drugpheno.nature), "_"), function(x) { return(all(x[c(length(x)-1, length(x))] == c("IC", "50"))) })
 #    ic50 <- drugpheno.nature[ , myx, drop=FALSE]
 #    nn <- dimnames(ic50)
 #    nn[[2]] <- gsub("_IC_50", "", nn[[2]])
 #    ic50 <- apply(ic50, 2, as.numeric)
 #    dimnames(ic50) <- nn
 #    ic50 <- exp(ic50) / 10^6
 #    ## camptothecin
 #    pdf(file.path(saveres, "cgp_camptothecin_mgh_wtsi_paper.pdf"))
 #    yy <- -log10(ic50[ , "drugid_195", drop=FALSE])
 #    xx <- -log10(ic50[ , "drugid_1003", drop=FALSE])
 #    ccix <- complete.cases(xx, yy)
 #    nnn <- sum(ccix)
 #    cc <- cor.test(x=xx, y=yy, method="spearman", use="complete.obs", alternative="greater")
 #    cci <- spearmanCI(x=cc$estimate, n=sum(ccix))
 #    par(mar=c(4, 4, 3, 1) + 0.1)
 #    llim <- round(range(c(xx, yy), na.rm=TRUE) * 10) / 10
 #    myScatterPlot(x=xx, y=yy, xlab="-log10 IC50 (WTSI)", ylab="-log10 IC50 (MGH)", main="CAMPTOTHECIN", pch=16, method="transparent", transparency=0.75)
 #    legend(x=par("usr")[1], y=par("usr")[4], xjust=0.075, yjust=0.85, bty="n", legend=sprintf("Rs=%.3g, p=%.1E, n=%i", cc$estimate, cc$p.value, nnn), text.font=2)
 #    dev.off()
 # 
 # 
 # 
 # 





  # ## drug information
#   if (verbose) { message("Read drug information") }
#   myfn2 <- file.path(tmpdir, "nature_supplinfo_druginfo_cgp.RData")
#   if(!file.exists(myfn2)) {
#     druginfo <- gdata::read.xls(xls=file.path(tmpdir, "Supplementary_data_final_Apr6.xlsx"), sheet=4)
#     druginfo[druginfo == "" | druginfo == " "] <- NA
#     save(list="druginfo", compress=TRUE, file=myfn2)
#   } else { load(myfn2) }
#   druginfo <- data.frame("drug_id"=gsub(pattern =badchars, replacement="", x=toupper(druginfo[ ,"Drug.ID"])), druginfo, stringsAsFactors=FALSE)
#   rownames(druginfo) <- druginfo[ ,"drug_id"] <- paste("drugid", as.character(druginfo[ ,"drug_id"]), sep="_")
# 
#   ## drug concentration
#   if (verbose) { message("Read drug concentration") }
#   drugconc <- read.csv(file.path(tmpdir, "cgp_drug_concentration.csv"), stringsAsFactors=FALSE)
#   drugconc[drugconc == "" | drugconc == " "] <- NA
#   drugconc <- data.frame("drug.name"=toupper(gsub(badchars, "", drugconc[ ,"Compound.Name"])), drugconc, stringsAsFactors=FALSE)
#   if(all(!is.element(drugconc[ , "drug.name"], drugnid[ , "drug.name"]))) { stop("Screening concentration for drugs ithout identifiers!") }
#   rownames(drugconc) <- rownames(drugnid)[match(drugconc[ , "drug.name"], drugnid[ , "drug.name"])]
#   drugconc <- data.frame("drug_id"=rownames(drugconc), drugconc, stringsAsFactors=FALSE)
# 
#   ## combine all drugs
#   dix <- sort(unique(c(rownames(druginfo), rownames(drugconc), paste("drugid", sapply(strsplit(colnames(drugpheno)[grep("^drugid_", colnames(drugpheno))], "_"), function(x) { return(x[[2]]) }), sep="_"))))
#   ## update druginfo
#   druginfo2 <- data.frame(matrix(NA, nrow=length(dix), ncol=ncol(druginfo), dimnames=list(dix, colnames(druginfo))))
#   newlev <- sapply(druginfo, levels)
#   newlev$drug_id <- sapply(strsplit(dix, split="_"), function(x) { return(x[2]) })
#   druginfo2 <- genefu::setcolclass.df(df=druginfo2, colclass=sapply(druginfo, class), factor.levels=newlev)
#   druginfo2[match(rownames(druginfo), dix), colnames(druginfo)] <- druginfo
#   druginfo2[ , "drug_id"] <- newlev$drug_id
#   druginfo <- druginfo2
#   ## update drugconc
#   drugconc2 <- data.frame(matrix(NA, nrow=length(dix), ncol=ncol(drugconc), dimnames=list(dix, colnames(drugconc))))
#   newlev <- sapply(drugconc, levels)
#   newlev$drug_id <- sapply(strsplit(dix, split="_"), function(x) { return(x[2]) })
#   drugconc2 <- genefu::setcolclass.df(df=drugconc2, colclass=sapply(drugconc, class), factor.levels=newlev)
#   drugconc2[match(rownames(drugconc), dix), colnames(drugconc)] <- drugconc
#   drugconc2[ , "drug_id"] <- newlev$drug_id
#   drugconc <- drugconc2
# 
#   ## report concentrations per cell line and per drug
#   drugconc2 <- data.frame(matrix(NA, nrow=nrow(drugconc) * length(cellnall), ncol=6, dimnames=list(paste(rep(cellnall, times=nrow(drugconc)), rep(rownames(drugconc), each=length(cellnall)), sep="..."), c("cell_id", "drug_id", "drug_name", "nbr_conc_tested", "min_Dose_uM", "max_Dose_uM"))))
#   drugconc2[ , "cell_id"] <- rep(cellnall, times=nrow(drugconc))
#   drugconc2[ , "drug_id"] <- rep(rownames(drugconc), each=length(cellnall))
#   drugconc2[ , "drug_name"] <- rep(as.character(drugconc[ ,"drug.name"]), each=length(cellnall))
#   ## as mentioned in the supplementary information of Garnett et al., a single cell line is used on each plate and treated with 28 different drugs over a 9-pt, 256-fold concentration range
#   drugconc2[ , "nbr_conc_tested"] <- 9
#   drugconc2[ , "min_Dose_uM"] <- rep(drugconc[ , "Min.Concentration.micromolar."], each=length(cellnall))
#   drugconc2[ , "max_Dose_uM"] <- rep(drugconc[ , "Max.Concentration.micromolar."], each=length(cellnall))
#   drugconc <- drugconc2
#   drugconc[, "cell_id"] <- as.character(drugconc[, "cell_id"])

## drug information
  if (verbose) { message("Read drug information") }
  druginfo <- read.csv(file.path(tmpdir, "cgp_drug_information.csv"))
  druginfo[!is.na(druginfo) & (druginfo == " " | druginfo == " ")] <- NA
  druginfo <- data.frame("drug.name"=toupper(gsub(badchars, "", druginfo[ ,"Name"])), druginfo)
  myx <- match(druginfo[ , "drug.name"], drugnid[ , "drug.name"])
  if (any(is.na(myx))) { stop ("Some drugs have missing annotations") }
  ## correct ambiguity for AZD6482: drugid_156 corresponds to the first occurence of AZD6482 while drugid_1066 corresponds to the second
  ## table(!is.na(drugpheno[ , "drugid_156_AUC"]))
  ## table(!is.na(drugpheno[ , "drugid_1066_AUC"]))
  myx[druginfo[ , "drug.name"] == "AZD6482"][2] <- which(drugnid[ , "drug.name"] == "AZD6482")[2]
  druginfo <- data.frame("drugid"=rownames(drugnid)[myx], drugnid[myx, , drop=FALSE], druginfo)
  rownames(druginfo) <- as.character(druginfo[ ,"drugid"])
  ## complement drug infomration with the supplementary infomration from the Nature website
  myfn2 <- file.path(tmpdir, "nature_supplinfo_druginfo_cgp.RData")
  if(!file.exists(myfn2)) {
    druginfo.nature <- gdata::read.xls(xls=file.path(tmpdir, "nature_supplementary_information.xls"), sheet=4)
    druginfo.nature[druginfo.nature == "" | druginfo.nature == " "] <- NA
    save(list="druginfo.nature", compress=TRUE, file=myfn2)
  } else { load(myfn2) }
  rownames(druginfo.nature) <- paste("drugid", druginfo.nature[ , "Drug.ID"], sep="_")
  druginfo <- data.frame(druginfo, druginfo.nature[rownames(druginfo), c("Brand.name", "Site.of.screening", "Drug.type", "Drug.class.I", "Drug.class.II", "Target.family", "Effector.pathway.biological.process", "Clinical.trials", "Source")])

  ## drug concentration
  message("Read drug concentration")
  drugconc <- read.csv(file.path(tmpdir, "cgp_drug_concentration.csv"))
  drugconc[!is.na(drugconc) & (drugconc == "" | drugconc == " ")] <- NA
  drugconc <- data.frame("drug.name"=toupper(gsub(badchars, "", drugconc[ ,"Compound.Name"])), drugconc)
  if(all(!is.element(drugconc[ , "drug.name"], drugnid[ , "drug.name"]))) { stop("Screening concentration for drugs without identifiers!") }
  myx <- match(drugconc[ , "drug.name"], drugnid[ , "drug.name"])
  ## correct ambiguity for AZD6482: drugid_156 corresponds to the first occurence of AZD6482 while drugid_1066 corresponds to the second
  myx[drugconc[ , "drug.name"] == "AZD6482"][2] <- which(drugnid[ , "drug.name"] == "AZD6482")[2]
  rownames(drugconc) <- rownames(drugnid)[myx]
  drugconc <- data.frame("drugid"=rownames(drugconc), drugconc)

  ## combine all drugs
  dix <- sort(unique(c(rownames(druginfo), rownames(drugconc), paste("drugid", sapply(strsplit(colnames(drugpheno)[grep("^drugid_", colnames(drugpheno))], "_"), function(x) { return(x[[2]]) }), sep="_"))))
  ## update druginfo
  druginfo2 <- data.frame(matrix(NA, nrow=length(dix), ncol=ncol(druginfo), dimnames=list(dix, colnames(druginfo))))
  newlev <- sapply(druginfo, levels)
  newlev$drugid <- sapply(strsplit(dix, split="_"), function(x) { return(x[2]) })
  druginfo2 <- genefu::setcolclass.df(df=druginfo2, colclass=sapply(druginfo, class), factor.levels=newlev)
  druginfo2[match(rownames(druginfo), dix), colnames(druginfo)] <- druginfo
  druginfo2[ , "drugid"] <- newlev$drugid
  druginfo <- druginfo2
  ## update drugconc
  drugconc2 <- data.frame(matrix(NA, nrow=length(dix), ncol=ncol(drugconc), dimnames=list(dix, colnames(drugconc))))
  newlev <- sapply(drugconc, levels)
  newlev$drugid <- sapply(strsplit(dix, split="_"), function(x) { return(x[2]) })
  drugconc2 <- genefu::setcolclass.df(df=drugconc2, colclass=sapply(drugconc, class), factor.levels=newlev)
  drugconc2[match(rownames(drugconc), dix), colnames(drugconc)] <- drugconc
  drugconc2[ , "drugid"] <- newlev$drugid
  drugconc <- drugconc2

  ## report concentrations per cell line and per drug
  drugconc2 <- data.frame(matrix(NA, nrow=nrow(drugconc) * length(cellnall), ncol=6, dimnames=list(paste(rep(rownames(drugconc), times=length(cellnall)), rep(cellnall, each=nrow(drugconc)), sep="..."), c("cellid", "drugid", "drug.name", "nbr.conc.tested", "min.Dose.uM", "max.Dose.uM"))))
  drugconc2[ , "cellid"] <- rep(cellnall, times=nrow(drugconc))
  drugconc2[ , "drugid"] <- rep(rownames(drugconc), each=length(cellnall))
  drugconc2[ , "drug.name"] <- rep(as.character(drugconc[ ,"drug.name"]), each=length(cellnall))
  ## as mentioned in the supplementary information of Garnett et al., a single cell line is used on each plate and treated with 28 different drugs over a 9-pt, 256-fold concentration range
  drugconc2[ , "nbr.conc.tested"] <- 9
  drugconc2[ , "min.Dose.uM"] <- rep(drugconc[ , "Min.Concentration.micromolar."], each=length(cellnall))
  drugconc2[ , "max.Dose.uM"] <- rep(drugconc[ , "Max.Concentration.micromolar."], each=length(cellnall))
  drugconc <- drugconc2
  drugconc[, "cell_id"] <- as.character(drugconc[, "cell_id"])
  




  ## IC50 in micro molar
  if (verbose) { message("Extracting IC50 values") }

  myx <- sapply(strsplit(colnames(drugpheno), "_"), function(x) { return(all(x[c(length(x)-1, length(x))] == c("IC", "50"))) })
  ic50 <- drugpheno[ ,myx,drop=FALSE]
  nn <- dimnames(ic50)
  nn[[2]] <- gsub("_IC_50", "", nn[[2]])
  ic50 <- apply(ic50, 2, as.numeric)
  dimnames(ic50) <- nn
  ic50 <- exp(ic50)

  ## sensitivity calling using waterfall plot
  ic50.call <- NULL
  for(i in 1:ncol(ic50)) {
    ic50.call <- cbind(ic50.call, callingWaterfall(x=ic50[ ,i], type="IC50", intermediate.fold=c(4, 1.2, 1.2), cor.min.linear=0.95, plot=FALSE, name=sprintf("%s (CGP)", colnames(ic50)[i])))
  }
  dimnames(ic50.call) <- dimnames(ic50)

  ## activity area
  if (verbose) { message("Extracting AUC values") }

  ## continuous values
  myx <- sapply(strsplit(colnames(drugpheno), "_"), function(x) { return(all(x[c(length(x))] == c("AUC"))) })
  auc <- drugpheno[ , myx, drop=FALSE]
  nn <- dimnames(auc)
  nn[[2]] <- gsub("_AUC", "", nn[[2]])
  auc <- apply(auc, 2, as.numeric)
  ## AUC for sensitivity
  auc <- 1 - auc
  dimnames(auc) <- nn
  
  ## sensitivity calling using waterfall plot
  auc.call <- NULL
  for(i in 1:ncol(auc)) {
    auc.call <- cbind(auc.call, callingWaterfall(x=auc[ ,i], type="AUC", intermediate.fold=c(4, 1.2, 1.2), cor.min.linear=0.95, plot=FALSE, name=sprintf("%s (CGP)", colnames(auc)[i])))
  }
  dimnames(auc.call) <- dimnames(auc)
   
  if (delete.tmpdir){
    ## delete temporary files
    sapply(tmpfiles, function (x) { unlink(x, recursive=TRUE) })
  }

  return (list("expression"=eset, "mutation"=mutation, "cell.line"=celline, "drug"=list("info"=druginfo, "concentration"=drugconc, "duration"=NULL, "ic50"=ic50, "ic50.call"=ic50.call, "auc"=auc, "auc.call"=auc.call)))
}

## End
