#annot_path <- system.file('data/dm6.annot.gtf.gz', package = 'SaiLoR')
annot_path <- system.file("exdata/dm6.annot.gtf.gz", package="SaiLoR")
junction_path <- system.file("exdata/short_read_junctions.SJ.out.tab", package = 'SaiLoR')


#annot_path <- "../data/dm6.annot.gtf.gz"
ref_annot <- rtracklayer::import.gff(annot_path)
# test junction reference generation
#junction_path <- "../data/short_read_junctions.SJ.out.tab"
test_that('create_reference_junctions correctly extracts reference', {
  expect_error(create_reference_junctions( junction_path, min.jcounts = 2 , ref_annot, type="short"),regexp=NA)
})

test_that('create_reference_junctions returns expected output genomic ranges',{

  ## Windows OS throws warnings
  suppressWarnings(
    reference_junctions <- create_reference_junctions( junction_path, min.jcounts = 2 , ref_annot, type="short")
  )
  expect_s4_class(reference_junctions, 'GRanges')
})


