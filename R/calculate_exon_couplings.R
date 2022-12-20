#' Title
#'
#' @param junctions
#' @param read_ids
#' @param pathJunctionRef
#'
#' @return
#' @export
#'
#' @examples
calculate_exon_couplings <- function(junctions,reference_junctions){
  couplingsResult <- list()
  res_couplings_tss <- compute_exon_tss_couplings(junctions,reference_junctions)
  res_couplings_tes <- compute_exon_3end_couplings(junctions,reference_junctions)
  couplingsResult$TSS <- res_couplings_tss
  couplingsResult$TES <- res_couplings_tes
  return(couplingsResult)
}
