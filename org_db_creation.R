
#--- read gene info file

gene_info <- readr::read_delim(file = "data_files/FungiDB-63_CneoformansH99_GO.txt", delim = "\t", comment = "!", col_names = FALSE) %>%
  dplyr::select(c(X2,X3,X10)) %>%
  dplyr::distinct()

colnames(gene_info) <- c("GID","SYMBOL", "GENENAME")

#--- read gff

cnGFF <- GenomicFeatures::makeTxDbFromGFF(file = "~/Documents/CDK7_project/RNASeq/AGRF_CAGRF220611046_HGHV5DRX2/GENOME/FungiDB-59_CneoformansH99.gff")

cn_chr <- GenomicFeatures::genes(cnGFF) %>%
  data.frame() %>%
  dplyr::distinct() %>%
  dplyr::select(c(gene_id, seqnames))

colnames(cn_chr) <- c("GID","CHROMOSOME")

#--- read GO table

GO_table <- readr::read_delim(file = "data_files/Cneoh99_go_data.txt",comment = "!", col_names = T)

GO_data <- GO_table %>%
  dplyr::select(c("GID", "GO_ID", "GO_EVIDENCE_CODE")) %>%
  dplyr::distinct()

colnames(GO_data) <- c("GID", "GO", "EVIDENCE")

term_to_gene <- GO_table %>% dplyr::select(c("GO_ID", "GID"))

term_to_name <- GO_table %>% dplyr::select(c("GO_ID","GO_TERM_NAME"))


AnnotationForge::makeOrgPackage(gene_info=gene_info, chromosome=cn_chr, go=GO_data,
                                version="0.1",
                                maintainer="pooja sethiya <pooja.sethiya@sydney.edu.au>",
                                author="pooja sethiya",
                                outputDir = ".",
                                tax_id="235443",
                                genus="Cryptococcus",
                                species="neoformansH99v64",
                                goTable="go")


install.packages("./org.CneoformansH99v64.eg.db", repos=NULL, type="source")
library("org.CneoformansH99v64.eg.db")
