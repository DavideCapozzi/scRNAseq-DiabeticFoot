# Ultra-robust R script: Match genes with manual aliases + BioMart + HGNC

library(tidyverse)
library(openxlsx)
library(stringr)
library(biomaRt)
library(httr)
library(jsonlite)

# == ===================================================================
# 1. INPUT FILES
# =====================================================================
deg_file <- "~/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/6_differential_expression_analysis/results/DE_healed_vs_not_healed_ALL_clusters.xlsx"
output_file <- "~/Library/CloudStorage/GoogleDrive-federica.dannunzio@uniroma1.it/Drive condivisi/sc-FEDE_DAVIDE/01_second_new_analysis/6_differential_expression_analysis/results/DEGs_Genes_of_Interest.xlsx"

# =====================================================================
# 2. GENE LISTS
# =====================================================================
visit1_genes <- c("IL-15RA", "IL-12B", "TNF", "CCL3", "CCL4", "MCP-4", "IL-18R1")
visit2_genes <- c("SPARCL1", "FCN2", "VASN", "CHL1", "CCL18", "SAA4", "F11", "SELL", "F7", "LTBP2", "CNDP1", "CRTAC1", "ICAM3")
all_genes <- unique(c(visit1_genes, visit2_genes))

# =====================================================================
# 3. MANUAL ALIAS LIST (curated and reliable)
# =====================================================================
manual_alias <- list(
  "IL15RA"   = c("IL15RA", "IL-15RA", "CD215"),
  "IL12B"    = c("IL12B", "IL-12B", "NFSF2", "CLMF"),
  "TNF"      = c("TNF", "TNFA", "TNFSF2"),
  "CCL3"     = c("CCL3", "MIP1A", "MIP-1A"),
  "CCL4"     = c("CCL4", "MIP1B", "MIP-1B"),
  "CCL13"    = c("MCP-4", "MCP4", "CCL13", "NCC-1"),
  "IL18R1"   = c("IL18R1", "IL-18R1", "IL18RA", "CD218A", "IL1RRP"),
  "SPARCL1"  = c("SPARCL1", "SC1", "MAST9", "HEVIN"),
  "FCN2"     = c("FCN2", "FCNL", "EBP-37", "P35", "Ficolin-2"),
  "VASN"     = c("VASN", "Vasorin", "SLITL2", "ATIA"),
  "CHL1"     = c("CHL1", "L1CAM2", "CALL", "NHBP"),
  "CCL18"    = c("CCL18", "PARC", "MIP-4", "MIP4", "DCCK1"),
  "SAA4"     = c("SAA4", "C-SAA", "CSAA", "A-SAA4"),
  "F11"      = c("F11", "FXI", "PTA"),
  "SELL"     = c("SELL", "CD62L", "LAM1", "LECAM1", "L-Selectin"),
  "F7"       = c("F7", "FVII", "SPCA"),
  "LTBP2"    = c("LTBP2"),
  "CNDP1"    = c("CNDP1"),
  "CRTAC1"   = c("CRTAC1"),
  "ICAM3"    = c("ICAM3", "CD50")
)

# =====================================================================
# 4. Retrieve synonyms from BioMart
# =====================================================================
get_biomart_synonyms <- function(genes) {
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  attrs <- c("hgnc_symbol", "external_synonym")
  res <- getBM(attributes = attrs, filters = "hgnc_symbol", values = genes, mart = ensembl)
  res <- res %>% filter(hgnc_symbol != "") %>% distinct()
  res_list <- list()
  for (g in genes) {
    syns <- res %>% filter(hgnc_symbol == g) %>% pull(external_synonym) %>% na.omit()
    syn_list <- unique(c(g, syns))
    res_list[[g]] <- syn_list
  }
  return(res_list)
}

# =====================================================================
# 5. Retrieve synonyms from HGNC REST API (robust)
# =====================================================================
get_hgnc_synonyms <- function(genes) {
  base_url <- "https://rest.genenames.org/fetch/symbol/"
  syn_list <- list()
  for (g in genes) {
    url <- paste0(base_url, g)
    resp <- tryCatch(GET(url, add_headers(`Accept` = 'application/json')), error = function(e) NULL)
    if (!is.null(resp) && status_code(resp) == 200) {
      data <- content(resp, as = "text", encoding = "UTF-8")
      json <- fromJSON(data, simplifyDataFrame = FALSE)
      if (!is.null(json$response$docs) && length(json$response$docs) > 0) {
        doc <- json$response$docs[[1]]
        if (!is.null(doc$alias_symbol) && length(doc$alias_symbol) > 0) {
          syn_list[[g]] <- unique(c(g, doc$alias_symbol))
        } else { syn_list[[g]] <- g }
      } else { syn_list[[g]] <- g }
    } else { syn_list[[g]] <- g }
  }
  return(syn_list)
}

# =====================================================================
# 6. Combine manual + BioMart + HGNC synonyms
# =====================================================================
biomart_syns <- get_biomart_synonyms(all_genes)
hgnc_syns <- get_hgnc_synonyms(all_genes)

alias_map <- list()
for (g in all_genes) {
  alias_map[[g]] <- unique(c(manual_alias[[g]], biomart_syns[[g]], hgnc_syns[[g]]))
}

alias_df <- stack(alias_map) %>% rename(Alias = values, Gene = ind) %>% mutate(Alias_clean = toupper(gsub("[^A-Z0-9]", "", Alias)))

# =====================================================================
# 7. Load DEG matrix
# =====================================================================
de <- read.xlsx(deg_file, sheet = 1)
if (!"gene" %in% tolower(names(de))) {
  guess <- which(sapply(de, function(x) any(grepl("^[A-Za-z0-9_-]+", head(x)))))
  names(de)[guess[1]] <- "gene"
}
de$gene_clean <- toupper(gsub("[^A-Z0-9]", "", de$gene))

# =====================================================================
# 8. Matching function
# =====================================================================
match_genes <- function(gene_list) {
  results <- list()
  for (g in gene_list) {
    aliases_clean <- alias_df %>% filter(Gene == g) %>% pull(Alias_clean)
    idx_exact <- which(de$gene_clean %in% aliases_clean)
    idx_regex <- unique(unlist(lapply(aliases_clean, function(al) grep(al, de$gene_clean))))
    all_idx <- unique(c(idx_exact, idx_regex))
    if (length(all_idx) > 0) {
      temp <- de[all_idx, ]
      temp$Target <- g
      temp$Match <- "Found"
      results[[g]] <- temp
    } else {
      missing <- tibble(Target = g, Match = "Not Found")
      results[[g]] <- missing
    }
  }
  bind_rows(results)
}

# =====================================================================
# 9. Run matching for VISIT1 + VISIT2
# =====================================================================
res_visit1 <- match_genes(visit1_genes)
res_visit2 <- match_genes(visit2_genes)

# =====================================================================
# 10. Export Excel with multiple sheets
# =====================================================================
wb <- createWorkbook()
addWorksheet(wb, "VISIT1_Found")
writeData(wb, "VISIT1_Found", res_visit1 %>% filter(Match == "Found"))
addWorksheet(wb, "VISIT1_Missing")
writeData(wb, "VISIT1_Missing", res_visit1 %>% filter(Match == "Not Found"))
addWorksheet(wb, "VISIT2_Found")
writeData(wb, "VISIT2_Found", res_visit2 %>% filter(Match == "Found"))
addWorksheet(wb, "VISIT2_Missing")
writeData(wb, "VISIT2_Missing", res_visit2 %>% filter(Match == "Not Found"))
saveWorkbook(wb, output_file, overwrite = TRUE)
cat("\nUltra-robust matching complete! File saved to:\n", output_file, "\n")
