#' Count full length reads
#'
#' @param alignments: alignments file from minimap2
#' @param linksDatabase: created with prepare_links_database
#'
#' @return
#' @export
#'
#' @examples
countLinks <- function(alignments, linksDatabase) {

  # make single nt starts
  startsAlignemnts <-  prepareForCountStarts(alignments, 1)
  startAlignments <-
    readTSSassignment(startsAlignemnts, linksDatabase$TSSCoordinate.base)
  startAlignments <-
    startAlignments %>% as.data.frame(.) %>% dplyr::select(name, promoter_id)
  # make single nt ends
  endsAlignemnts <-  prepareForCountEnds(alignments, 1)
  endsAlignemnts <-
    readTESassignment(endsAlignemnts, linksDatabase$TESCoordinate.base)
  endsAlignemnts <-
    endsAlignemnts %>% as.data.frame(.) %>% dplyr::select(name, tes_id)
  # make pairs
  pairsTested <-
    left_join(endsAlignemnts, startAlignments, by = "name") %>% dplyr::filter(gsub("\\:.*", "", tes_id) == gsub("\\:.*", "", promoter_id))
  pairsTested <-
    left_join(
      pairsTested,
      linksDatabase$pairDataBase %>% dplyr::distinct(tes_id, .keep_all = TRUE) %>% dplyr::select(gene_id, tes_id),
      by = "tes_id"
    )  %>%
    mutate(pairs_id = paste0(
      gene_id,
      ":",
      gsub(".*:", "", pairsTested$promoter_id),
      ":",
      gsub(".*:", "", pairsTested$tes_id)
    ))
  # summarize only reads that expand both pairs
  countedPairs <- pairsTested %>% group_by(pairs_id) %>% tally()
  countedPairsFinal <-
    left_join(
      countedPairs,
      pairsTested %>% dplyr::distinct(pairs_id, .keep_all = TRUE) %>% dplyr::select(gene_id, pairs_id),
      by = "pairs_id"
    )
  # normalize expression in CPMs
  dd <- edgeR::DGEList(counts = countedPairsFinal$n)
  dge <- edgeR::calcNormFactors(dd)
  countedPairsFinal$pairs_cpm <- as.numeric (edgeR::cpm(dge))
  res <- list()
  res$countedPairsFinal <- countedPairsFinal
  res$pairsTested <- pairsTested
  return(res)
}
#' prepare counts trimming the reads to their most 5' into a 1 nt window for counting
#'
#' @param x alignments: alignments file from minimap2
#' @param window window for trimming the reads to their most 5'/3'
#'
#' @return
#' @export
#'
#' @examples
prepareForCountStarts <- function(x, window) {
  alignments <- x
  pos <- alignments[alignments@strand == "+",]
  neg <- alignments[alignments@strand == "-",]
  GenomicRanges::end(pos) <- GenomicRanges::start(pos) + window
  GenomicRanges::start(neg) <- GenomicRanges::end(neg) - window
  shortstarts <- c(pos, neg)
  return(shortstarts)
}
#' prepare counts trimming the reads to their most 3' into a 1 nt window for counting
#'
#' @param x  alignments: alignments file from minimap2
#' @param window window for trimming the reads to their most 5'/3'
#' @return
#' @export
#'
#' @examples
prepareForCountEnds <- function(x, window) {
  alignments <- x
  pos <- alignments[alignments@strand == "+",]
  neg <- alignments[alignments@strand == "-",]
  GenomicRanges::start(pos) <- GenomicRanges::end(pos) - window
  GenomicRanges::end(neg) <- GenomicRanges::start(neg) + window
  shortstarts <- c(pos, neg)
  return(shortstarts)
}
#' Title
#'
#' @param startsAlignemnts  trimmed reads  output from prepareForCountStarts
#' @param TSSCoordinate.base reference database for TSS
#'
#' @return
#' @export
#'
#' @examples
readTSSassignment <- function(startsAlignemnts, TSSCoordinate.base) {
  tssDb <- TSSCoordinate.base
  ovlps <- GenomicRanges::findOverlaps(startsAlignemnts , tssDb , maxgap = 50)
  promoterIds <- tssDb[subjectHits(ovlps),]$count
  promoterStarts <- GenomicRanges::start(tssDb[subjectHits(ovlps),])
  startsAlignemnts2 <- startsAlignemnts[queryHits(ovlps),]
  startsAlignemnts2$promoter_id <- promoterIds
  startsAlignemnts2$promoterStarts <- promoterStarts
  startsAlignemnts2$dist2Assignment <-
    abs(GenomicRanges::start(startsAlignemnts2) - startsAlignemnts2$promoterStarts)
  startsAlignemnts2 <-
    startsAlignemnts2[order(startsAlignemnts2$dist2Assignment, decreasing = FALSE),]
  startsAlignemnts2 <-
    startsAlignemnts2[!duplicated(startsAlignemnts2$name),]
  return(startsAlignemnts2)
}
#' Title
#'
#' @param endsAlignements trimmed reads  output from prepareForCountEnds
#' @param TESCoordinate.base reference database for TES
#'
#' @return
#' @export
#'
#' @examples
readTESassignment <- function(endsAlignements, TESCoordinate.base) {
  tesDb <- TESCoordinate.base
  ovlps <- findOverlaps(endsAlignements , tesDb , maxgap = 150)
  tesIds <- tesDb[subjectHits(ovlps),]$count
  endStarts <- GenomicRanges::start(tesDb[subjectHits(ovlps),])
  endsAlignemnts2 <- endsAlignements[queryHits(ovlps),]
  endsAlignemnts2$tes_id <- tesIds
  endsAlignemnts2$endStarts <- endStarts
  endsAlignemnts2$dist2Assignment <-
    abs(GenomicRanges::start(endsAlignemnts2) - endsAlignemnts2$endStarts)
  endsAlignemnts2 <-
    endsAlignemnts2[order(endsAlignemnts2$dist2Assignment, decreasing = FALSE),]
  endsAlignemnts2 <-
    endsAlignemnts2[!duplicated(endsAlignemnts2$name),]
  return(endsAlignemnts2)
}

#' Identify full length reads for downstream analysis
#'
#' @param annotPath : reference annotaiton path
#' @param alignmentsFile : bam file alignments object
#' @param tss.ntwindow : window for the database assignment for full length at thes TSS
#' @param tes.ntwindow : window for the database assignment for full length at thes TES
#'
#' @return
#' @export
#'
#' @examples
get_sailor_full_lengths <- function(annotPath, alignmentsFile, tss.ntwindow, tes.ntwindow) {
  annot <- rtracklayer::import.gff(annotPath)
  linksDatabase <- prepareLinksDatabase(annot, tss.window=tss.ntwindow, tes.window=tes.ntwindow)
  alignmentsFile <- GenomicRanges::GRanges(alignmentsFile)
  alignmentsFile$name <- names(alignmentsFile)
  names(alignmentsFile) <- NULL
  countsLinks <- countLinks(alignmentsFile, linksDatabase)
  return(countsLinks)
  }



