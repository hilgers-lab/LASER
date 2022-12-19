#' Compute couplings exon junctions and 3'ends
#'
#' @param junctions: junction data with to junction-reads assignments
#' @param read_ids: read_ids for full length reads
#' @param pathJunctionRef: reference database used.
#'
#' @return
#' @export
#'
#' @examples
compute_exon_3end_couplings <- function(junctions,read_ids,pathJunctionRef){
  juncSums <- readr::read_tsv(junctions) %>% dplyr::select(new_junID, read_id)
  fullLengths <- readr::read_tsv(read_ids)
  junctionRef<- rtracklayer::import.gff(pathJunctionRef)
  # keep only junctions in database
  juncSums <- juncSums %>% filter(new_junID %in% junctionRef$juncID)
  readsToJunctions <- left_join(fullLengths,
                     juncSums %>% dplyr::rename(name = read_id),
                     by = "name")
  # summary junction combinations found.add less memory consuming id
  tmpRefjun <- readsToJunctions %>%
    distinct(new_junID, .keep_all = TRUE) %>%
    arrange(new_junID) %>%
    group_by(gene_id) %>%
    mutate(jun_id = paste0("J", sequence(n()))) %>%
    ungroup()
  readsToJunctions <- left_join(
    readsToJunctions,
    tmpRefjun %>% dplyr::select(new_junID,jun_id),
    by =  "new_junID"
  )
  readsToJunctionsCollapse <- readsToJunctions %>%
    group_by(name) %>%
    mutate(juncComb = paste0(jun_id, collapse = ",")) %>% mutate(juncComb=paste0(gene_id, ":", juncComb))
  junCombCounts <- readsToJunctionsCollapse %>%
    distinct(name , .keep_all = TRUE) %>%
    group_by(juncComb) %>%
    tally()
  # keep only genes with more than 2 junction combinations with at least 2 reads each
  message("working with genes with at least 2 junctions expressed with more than 2 read counts")
  junCombCounts <-
    left_join(
      junCombCounts,
      readsToJunctionsCollapse %>%
        ungroup() %>%
        distinct(juncComb, .keep_all = TRUE) %>%
        dplyr::select(juncComb, gene_id),
      by = "juncComb"
    )
  genesWithMoreThan1JunExpressed <- junCombCounts %>%
    filter(n >= 2) %>%
    group_by(gene_id) %>%
    tally() %>%
    filter(n > 1)
  message("genesToWork:", nrow(genesWithMoreThan1JunExpressed))
  # lets start with tss
  # make pairs ids for junction-tss combinations
  testSub <-
    readsToJunctions %>%
    filter(gene_id %in% genesWithMoreThan1JunExpressed$gene_id) %>%
    mutate(
      tssJunID = paste0(tes_id, ":", new_junID),
      isoID = paste0(promoter_id, ":", new_junID, ":", tes_id)
    )
  cat("test")
  # lets make small database to annotate the counts
  simpleT <-
    testSub %>%
    dplyr::select(tes_id, promoter_id, new_junID, tssJunID, isoID) %>%
    distinct(isoID, .keep_all = TRUE)
  ### summarize junctions
  # counts per pair
  tssJunCounts <- testSub %>%
    group_by(tssJunID) %>%
    tally()
  tssJunCounts <-
    left_join(
      tssJunCounts,
      simpleT %>% distinct(tssJunID, .keep_all = TRUE) %>% dplyr::select(tssJunID, tes_id, new_junID),
      by = "tssJunID"
    ) %>% mutate(gene_id = gsub("\\:.*", "", .data$tes_id) )
  annotPairsExp <- tssJunCounts
  # make table of junction fractions
  annotPairsExp <-
    annotPairsExp %>% dplyr::group_by(new_junID) %>% dplyr::mutate(jun_sum = sum(n)) %>%
    group_by(tes_id) %>% dplyr::mutate(start_sum = sum(n)) %>% dplyr::group_by(tssJunID) %>% dplyr::mutate(pairs_sum = sum(n)) %>%
    group_by(gene_id) %>% mutate(geneMean = sum(n))
  # FilterGenes with more than one 3 ends
  genesATSS <-
    annotPairsExp %>% dplyr::distinct(tes_id, .keep_all = TRUE) %>% dplyr::group_by(gene_id) %>% dplyr::tally() %>% dplyr::filter(n >1) %>% dplyr::pull(gene_id)
  annotPairsExp <- annotPairsExp %>% dplyr::filter(gene_id %in% genesATSS)
  # Keep only genes with mutiple junction combinations
  genesMJC <-
    annotPairsExp %>% dplyr::distinct(new_junID, .keep_all = TRUE) %>% dplyr::group_by(gene_id) %>% tally() %>% filter(n >1) %>% dplyr::pull(gene_id)
  annotPairsExp <- annotPairsExp %>% filter(gene_id %in% genesMJC)
  ### Testing
  ##########################################################################################
  perGeneList <- split(annotPairsExp, f = annotPairsExp$gene_id)
  couplingsmatrix <- lapply(perGeneList, function(x) {
    x1 <-
      x %>% maditr::dcast(new_junID ~ tes_id, value.var = "n") %>% as.data.frame(.)
    x1 <- na.omit(x1)
    if (nrow(x1) > 0) {
      rownames(x1) <- x1$new_junID
      x1$new_junID <- NULL
      return(x1)
    } else{
      return(x1)
    }
  })
  couplingsmatrix <- Filter(function(x)
    nrow(x) > 1, couplingsmatrix)
  couplingsmatrix <-
    lapply(couplingsmatrix, function(x) {
      x[is.na(x)] <- 0
      x
    })
  couplingsmatrix <- lapply(couplingsmatrix, function(x) {
    x + 1
  })
  ## fixes remove single junction genes.
  perGeneChisqTest <- lapply(couplingsmatrix, function(x) {
    # MONTE CARLO SIMULATION P-VALUE
    x2 <- chisq.test(x, simulate.p.value = TRUE)
    x1 <- x2$p.value
    return(x1)
  })

  perGeneSquareRes <- lapply(couplingsmatrix, function(x) {
    x2 <- chisq.test(x, simulate.p.value = TRUE)
    # RESIDUALS from chisq.test
    x1 <- x2$statistic
    return(x1)
  })

  makeDFPerGene <- function(x){
    td <- chisq.test(x)
    tdo <- as.data.frame(td$observed) %>%
      dplyr::mutate(new_junID=rownames(.))
    tdo <- reshape2::melt( tdo) %>%
      mutate(pairs_id = paste0(new_junID,":",variable)) %>%
      dplyr::select(pairs_id, value) %>%
      dplyr::rename(observed=value)
    tde <- as.data.frame(td$expected) %>%
      dplyr::mutate(new_junID=rownames(.))
    tde <- reshape2::melt( tde) %>%
      mutate(pairs_id = paste0(new_junID,":",variable)) %>%
      dplyr::select(pairs_id, value) %>%
      dplyr::rename(expected=value)
    tdr <- as.data.frame(td$residuals) %>%
      dplyr::mutate(new_junID=rownames(.))
    tdr <- reshape2::melt( tdr) %>%
      dplyr::mutate(pairs_id = paste0(new_junID,":",variable),
             gene_id=gsub("\\:.*","",new_junID) ) %>%
      dplyr::select(pairs_id, value, gene_id) %>%
      dplyr::rename(residuals=value)
    d<- left_join(tdo, tde, by="pairs_id")
    d<- left_join(d, tdr, by="pairs_id")
    return(d)
  }
  x2 <-  lapply( couplingsmatrix, makeDFPerGene)
  x3 <- do.call(rbind,x2)
  longTable <- x3
  tt2 <- test2 %>%
    dplyr::select(new_junID,jun_id) %>%
    distinct(new_junID, .keep_all = TRUE) %>%
    dplyr::rename(gene_id=new_junID)
  longTable <- left_join(longTable,  tt2, by="gene_id")
  affectedGenes <- as.data.frame(do.call(rbind, perGeneChisqTest)) %>%
    dplyr::rename(p.value.chisq = 1) %>%
    dplyr::mutate(gene_id = rownames(.))
  affectedGenesXsquare <- as.data.frame(do.call(rbind, perGeneSquareRes)) %>%
    dplyr::rename(sqrtRes = 1) %>%
    dplyr::mutate(gene_id = rownames(.))
  affectedGenes<- left_join(affectedGenes, affectedGenesXsquare, by="gene_id")
  affectedGenes$adj.pvalue <- p.adjust(affectedGenes$p.value.chisq)
  resultTest  <- list()
  resultTest$couplingsPerGene  <- affectedGenes
  resultTest$couplingsPerJunction  <- longTable
 return(resultTest)
}

