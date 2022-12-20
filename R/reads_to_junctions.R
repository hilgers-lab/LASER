#' junctions to reads assignments
#'
#' @param bamPath: path to bed file
#' @param refJunAnnot: path to reference junctions in genomic ranges format
#'
#'
#' @return
#' @export
#'
#' @examples
read_to_junctions <- function(bamPath,refJunAnnot, pathAnnot){
  bamAlignments <- GenomicAlignments::readGAlignments(bamPath, use.names = TRUE)
  # filter reads by full length
  fullLengths <- get_sailor_full_lengths(pathAnnot, bamAlignments, tss.ntwindow=50, tes.ntwindow=150)
  bamAlignments <- bamAlignments[names(bamAlignments) %in% fullLengths$pairsTested$name,]
  read_junctions <- GenomicAlignments::junctions(bamAlignments, use.mcols = TRUE)
  reference_junction_data <- read_refjunction_correction(read_junctions, refJunAnnot, 10, 10)
  bam_recovered_junctions <-  make_junction_database(read_junctions)
  reads_to_junctions <- inner_join(bam_recovered_junctions, reference_junction_data, by = "raw_juncID")
  read.features <- fullLengths$pairsTested %>% group_by(name) %>% distinct(name, .keep_all= TRUE)%>% dplyr::rename(read_id=name)
  reads_to_junctions <- left_join(reads_to_junctions, read.features, by="read_id")
  return(reads_to_junctions)
}
#' Title
#'
#' @param bamJunctions
#' @param shortReadJunctions
#' @param dist.5ss
#' @param dist.3ss
#' @importFrom rlang .data
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges findOverlaps
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @return
#' @export
#'
#' @examples
read_refjunction_correction <- function(bamJunctions,
                                        shortReadJunctions,
                                        dist.5ss,
                                        dist.3ss) {
  # Take bam file and summarize all junctions
  tmpRef <- unique(unlist(bamJunctions))
  tmpRef$raw_juncID <-
    paste0(seqnames(tmpRef), ":",
           start(tmpRef), "-",
           end(tmpRef))
  refJunctions <- tmpRef
  # Assign found junctions in bam to reference junctions
  hitsovlps <- findOverlaps(refJunctions, shortReadJunctions)
  junctionsInShort <- refJunctions[queryHits(hitsovlps), ]
  cat("Assigned to junctions reads: " , length(unique(junctionsInShort)))
  junctionsInShort$start_short <-
    start(shortReadJunctions[subjectHits(hitsovlps),])
  junctionsInShort$end_short <-
    end(shortReadJunctions[subjectHits(hitsovlps),])
  junctionsInShort$strand_short <-
    strand(shortReadJunctions[subjectHits(hitsovlps),])
  junctionsInShort$dist_start <-
    abs(start(shortReadJunctions[subjectHits(hitsovlps),]) -  start(junctionsInShort))
  junctionsInShort$dist_end <-
    abs(end(shortReadJunctions[subjectHits(hitsovlps),]) -  end(junctionsInShort))
  junctionsShortAssignments <-
    junctionsInShort %>% as.data.frame(., row.names = NULL) %>%
    arrange(dist_end + dist_start) %>%
    dplyr::mutate (totalDist = dist_end + dist_start) %>%
    dplyr::distinct(raw_juncID, .keep_all = TRUE)
  # filter by distance and assign new id
  cat("Junctions before filtering " ,
      nrow(junctionsShortAssignments))
  junctionsShortAssignments <-
    junctionsShortAssignments %>%
    dplyr::filter(dist_start <= dist.5ss &
                    dist_end <= dist.3ss) %>%
    dplyr::mutate(new_junID = paste0(seqnames, ":",
                                     start_short, "-",
                                     end_short))
  cat("Junctions after filtering " ,
      nrow(junctionsShortAssignments))
  return(junctionsShortAssignments)
}

#' merge junctions to reads
#'
#' @param bamJunctions
#'
#' @return
#' @export
#'
#' @examples
make_junction_database <- function(bamJunctions) {
  junBam <- unlist(bamJunctions)
  junBam$raw_juncID <-
    paste0(seqnames(junBam),":",
           start(junBam),"-",
           end(junBam))
  junBam$read_id <- names(unlist(bamJunctions))
  # add later
  #junBam[junBam$read_id %in% fullLength,]
  junBam <-
    junBam %>% as.data.frame(., row.names = NULL) %>% dplyr::select(raw_juncID, read_id)
  return(junBam)
}




