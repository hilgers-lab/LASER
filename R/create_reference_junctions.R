#' Load short read sjb data
#'
#' @param sjbPath
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @return
#' @export
#'
#' @examples
read_sjout <- function(sjbPath) {
  junctionFile <- readr::read_table(sjbPath, col_names = FALSE)
  colnames(junctionFile) <-
    c(
      "seqnames",
      "start",
      "end",
      "strand",
      "intronMotif",
      "comp.ref",
      "read.count",
      "multimappers",
      "overhang"
    )

  junctionFile <- junctionFile %>%
    dplyr::mutate(
      juncID = paste0(.data$seqnames, ":", .data$start, "-", .data$end),
      strand = dplyr::case_when(
        .data$strand == 0 ~ "*",
        .data$strand == 1 ~ "+",
        TRUE ~ "-"
      )
    )
  return(junctionFile)
}
#' Extract junctions from annotation
#'
#' @param refAnnot
#' @importFrom rlang .data
#' @return
#' @export
#'
#' @examples
make_annotation_txf <- function(refAnnot){
  txdb <- GenomicFeatures::makeTxDbFromGRanges(refAnnot)
  annot_junctions <- SGSeq::convertToTxFeatures(txdb)
  annot_junctions <- annot_junctions[SGSeq::type(annot_junctions) == "J"]
  GenomicRanges::start(annot_junctions) <- GenomicRanges::start(annot_junctions) + 1
  GenomicRanges::end(annot_junctions) <- GenomicRanges::end(annot_junctions) - 1
  return(annot_junctions)
}

#' Annotate junctions
#'
#' @param AnnotJunctions
#' @param refAnnot
#' @param linksDatabase
#' @importFrom rlang .data
#' @return
#' @export
#'
#' @examples
get_junctions_from_annot <- function(AnnotJunctions,linksDatabase,genesMultipleJunctions){
  # overlap junctions to windows to filter
  message("Removing junctions overlapping TSS and 3'end windows")
  junctionsInPromoters <- GenomicRanges::findOverlaps(AnnotJunctions, linksDatabase$promoterDependentWindow)
  juncNotIn5prime <- AnnotJunctions[-S4Vectors::queryHits(junctionsInPromoters)]
  junctionsIn3ends <- GenomicRanges::findOverlaps(AnnotJunctions, linksDatabase$tesDependentWindow)
  juncNotIn3prime <- AnnotJunctions[-S4Vectors::queryHits(junctionsIn3ends)]
  junctionsToWork <- unique( c(juncNotIn5prime,juncNotIn3prime))
  # keep only genes with more than one junction combinations
  message("Removing genes with more")
  # remove single junction genes
  message("Removing single CDS isoform genes")
  junctionsToWork <- junctionsToWork %>%
    as.data.frame(.) %>% tidyr::unnest(., geneName) %>%
    dplyr::filter(geneName %in% genesMultipleJunctions$geneName) %>%
    GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  junctionsToWork <- unique(junctionsToWork)
  junctionsToWork$txName <- NULL
  junctionsToWork$juncID <-
    paste0( GenomicRanges::seqnames(junctionsToWork),
           ":",
           GenomicRanges::start(junctionsToWork),
           "-",
           GenomicRanges::end(junctionsToWork))
  GenomicRanges::mcols(junctionsToWork) <-  GenomicRanges::mcols(junctionsToWork)[c(2,3)]
  junctionRef <- junctionsToWork
  junctionRef$gene_id <- junctionRef$geneName
  junctionRef$geneName <- NULL
  message("Reference junctions ready")
  return(junctionRef)
}

#' Get junctions from short read sequencing data
#'
#' @param pathSJ.out
#' @param min.counts
#' @param AnnotJunctions
#' @param refAnnot
#' @param linksDatabase
#' @importFrom rlang .data
#' @return
#' @export
#'
#' @examples
get_junctions_from_short_reads <- function(pathSJ.out, min.counts,refAnnot,linksDatabase,genesMultipleJunctions){
  message("Load splice junction files")
  junctionCounts <- read_sjout(pathSJ.out)
  junctionCounts <- junctionCounts %>%
    dplyr::filter(.data$read.count > min.counts) %>%
    dplyr::select(1, 2, 3, 4, 7, 10) %>%
    GenomicRanges::makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  # overlap junctions to windows to filter
  junctionsInPromoters <- GenomicRanges::findOverlaps(junctionCounts, linksDatabase$promoterDependentWindow)
  juncNotIn5prime <- junctionCounts[-S4Vectors::queryHits(junctionsInPromoters)]
  junctionsIn3ends <- GenomicRanges::findOverlaps(junctionCounts, linksDatabase$tesDependentWindow)
  juncNotIn3prime <- junctionCounts[-S4Vectors::queryHits(junctionsIn3ends)]
  shortReadJunctions <- unique( c(juncNotIn5prime,juncNotIn3prime))
  # annotate junctions to genes
  geneBins <- unique(unlist( range( split(refAnnot, f=refAnnot$gene_id) )))
  hits2Genes <- GenomicRanges::findOverlaps(  shortReadJunctions, geneBins, type = "within")
  shortReadJunctions <- shortReadJunctions[S4Vectors::queryHits(hits2Genes),]
  shortReadJunctions$gene_id <- names(geneBins[S4Vectors::subjectHits(hits2Genes),])
  GenomicRanges::mcols(shortReadJunctions)<- GenomicRanges::mcols(shortReadJunctions)[c("gene_id", "juncID")]
  shortReadJunctions <-  unique(shortReadJunctions)
  # keep only genes with more than one junction combinations
  message("single isoform genes")
  shortReadJunctions <- shortReadJunctions[shortReadJunctions$gene_id %in% genesMultipleJunctions$geneName,]
  return(shortReadJunctions)
}
#' Title
#'
#' @param pathSJ.out
#' @param min.jcounts
#' @param AnnotJunctions
#' @param refAnnot
#' @param type
#' @importFrom rlang .data
#' @return
#' @export
#'
#' @examples
create_reference_junctions <- function(pathSJ.out, min.jcounts, refAnnot, type){
  message("Extracting junctions from annotation")
  AnnotJunctions <- make_annotation_txf(refAnnot)
  message("Preparing Annotation of 5'and 3' ends")
  linksDatabase <- prepareLinksDatabase(refAnnot,
                                              tss.window = 50,
                                              tes.window = 150)
  # window of exons dependent on ATSS usage
  linksDatabase$promoterDependentWindow <-
    unlist( range(
      split(
        linksDatabase$TSSCoordinate.base,
        f = linksDatabase$TSSCoordinate.base$value.group_name
      )
    ))

  # window of exons dependent on alternative 3'end usage
  linksDatabase$tesDependentWindow <-
    unlist( range(
      split(
        linksDatabase$TESCoordinate.base,
        f = linksDatabase$TESCoordinate.base$value.group_name
      )
    ))

  # get genes with more 1 different junction
  genesMultipleJunctions <- AnnotJunctions %>% as.data.frame(.) %>%
    tidyr::unnest(., txName) %>% tidyr::unnest(., geneName) %>%
    mutate(JunID = paste0(seqnames, start, end)) %>%
    dplyr::group_by(txName) %>%
    dplyr::mutate(junName = paste0(JunID, collapse = "|")) %>%
    dplyr::distinct(junName, .keep_all = TRUE) %>%
    dplyr::group_by(geneName) %>% dplyr::tally() %>% dplyr::filter(n>1)

  if(type=="short"){
    message("Extracting junctions from short reads")
    shortReadJunctions <- get_junctions_from_short_reads(pathSJ.out,
                                                         min.counts=min.jcounts,
                                                         refAnnot,
                                                         linksDatabase,
                                                         genesMultipleJunctions)
    message("Extracting junctions from reference annotation")
    annotationJunctions <- get_junctions_from_annot(AnnotJunctions, linksDatabase,genesMultipleJunctions)
    reference_junctions <- c(shortReadJunctions, annotationJunctions)
    return(reference_junctions)
  }else{
    annotation_junctions <- get_junctions_from_annot(AnnotJunctions, refAnnot)
    return(annotation_junctions)
  }
}
