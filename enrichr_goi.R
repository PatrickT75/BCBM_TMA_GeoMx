library(enrichR)

websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}

dbs <- c("MSigDB_Hallmark_2020")

## find overlaps
# overlap <- list(Tumor.Brain=v1, Tumor.Breast=v2, Immune.Breast=v3, Immune.Brain=v4)
# ggvenn(overlap)
# ggsave("figs/venn/tumor-brain-vs-breast-venn-diagram.png")


## find genes of interest
# genes_of_interest <- setdiff(v3, union(v4, union(v2, v1)))

do_enrichr <- function(genelist){
if (websiteLive) {
  ###### define list of genes
  
  # ## some overlap of genes
  # enriched <- enrichr(genes_of_interest, dbs)
  
  # ## volcano plot significant genes 
  # enriched <- enrichr(rownames(tumor_lme %>% filter(diffexpressed == "UP")), dbs)
  
  ## genes from a sample-heatmap cluster
  
  genelist <- genelist
  enriched <- enrichr(genelist, dbs)
  
  for(db in dbs){
    print(plotEnrich(enriched[[db]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value"))
    # ggsave("figs/heatmap-enrichr/nopdx-enrichr/immune-brain-tnbc-vs-hrpos-cluster4.png")
    # write.csv(enriched[[db]], paste("tables/venn/all-tumor-genes-", db, ".csv", sep=""), row.names=FALSE)
  }
}
}
