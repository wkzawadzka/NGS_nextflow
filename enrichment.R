# -----------------------------
# LOAD LIBRARIES
# -----------------------------
library(dplyr)
library(readr)
library(clusterProfiler)
library(org.Gg.eg.db)
library(ggplot2)
library(enrichplot)
library(flextable)
library(officer)
library(purrr)

# -----------------------------
# GLOBAL VARIABLES
# -----------------------------
font.size <- 10
today <- Sys.Date()
thin_black_border <- fp_border(color = "black", width = 1)

# -----------------------------
# FUNCTION: CREATE FLExtable
# -----------------------------
make_flextable <- function(df, fontname = "Calibri", fontsize_header = 9, fontsize_body = 9, fontsize_gene = 8) {
  df %>%
    flextable::flextable() %>%
    bold(bold = TRUE, part = "header") %>%
    fontsize(size = fontsize_header, part = "header") %>%
    fontsize(size = fontsize_body, part = "body") %>%
    flextable::font(fontname = fontname, part = "all") %>%
    fontsize(j = "Gene IDs", size = fontsize_gene) %>%
    hline_top(border = thin_black_border, part = "all") %>%
    hline_bottom(border = thin_black_border, part = "all") %>%
    colformat_char(na_str = "-") %>%
    colformat_date(na_str = "-") %>%
    colformat_datetime(na_str = "-") %>%
    colformat_double(na_str = "-") %>%
    colformat_int(na_str = "-") %>%
    colformat_lgl(na_str = "-") %>%
    colformat_num(na_str = "-")
}

# -----------------------------
# INPUT FILE
# -----------------------------
file <- "GATK_VEP_Q04KhDtKbq340t7l.txt"
title <- "gatk"

# -----------------------------
# READ AND PROCESS SNP DATA
# -----------------------------
snps <- read_tsv(file) %>%
  dplyr::rename(
    `Input` = `#Uploaded_variation`,
    `Gene symbol` = SYMBOL,
    `Clinical significance` = CLIN_SIG,
    Impact = IMPACT,
    Biotype = BIOTYPE,
    `Existing variation` = Existing_variation
  ) %>%
  filter(!(`Gene symbol` == "-")) %>%
  # summarize per ENSG (just write after comma if multiple lines)
  group_by(`Gene symbol`, `Input`, `Gene`) %>%
  summarise(
    across(
      everything(),
      ~ paste(unique(na.omit(.x)), collapse = ", "),
      .names = "{.col}"
    ),
    .groups = "drop"
  )
snps_ensembl <- unique(snps$Gene)

# -----------------------------
# MAP ENSEMBL TO ENTREZ
# -----------------------------
entrez_ids <- mapIds(
  org.Gg.eg.db,
  keys = snps_ensembl,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
) %>% na.omit() %>% unique()

# -----------------------------
# FUNCTION: ENRICHMENT ANALYSIS
# -----------------------------
run_enrichGO <- function(entrez_ids, OrgDb = org.Gg.eg.db) {
  ego <- enrichGO(gene = entrez_ids, OrgDb = OrgDb, keyType = "ENTREZID", ont = "BP")
  egox <- setReadable(ego, OrgDb, "ENTREZID")
  return(egox)
}

ego <- run_enrichGO(entrez_ids)
ego_results <- as.data.frame(ego)

# -----------------------------
# FORMAT RESULTS FOR REPORT
# -----------------------------
format_enrichment_table <- function(enrich_df) {
  enrich_df %>%
    dplyr::select(
      -c(zScore, pvalue, qvalue, Count, RichFactor, GeneRatio, BgRatio)
    ) %>%
    mutate(
      geneID = gsub("/", ", ", geneID),
      p.adjust = round(p.adjust, 5),
      FoldEnrichment = round(FoldEnrichment, 3)
    ) %>%
    dplyr::rename(
      `Fold Enrichment` = FoldEnrichment,
      `Gene IDs` = geneID,
      `Adjusted p-value` = p.adjust,
      `GO ID` = ID,
      `Biological Process` = Description
    ) %>%
    arrange(desc(`Fold Enrichment`)) %>%
    flextable() %>%
    fontsize(size = 9, part = "header") %>%
    fontsize(size = 9) %>%
    flextable::font(fontname = "Calibri") %>%
    fontsize(j = "Gene IDs", size = 8) %>%
    bold(j = 2) %>%
    flextable::width(j = 1:5, c(0.7, 1.4, 0.9, 0.9, 4))
}

tab_ego <- format_enrichment_table(ego_results)

# -----------------------------
# CREATE DOCX REPORT
# -----------------------------
doc <- read_docx("szablon.docx") %>%
  body_add_par(paste("Date:", today)) %>%
  body_add_par("") %>%
  body_add_par("Methods", style = "heading 1") %>%
  body_add_par("") %>%
  body_add_par(
    "To identify biological pathways overrepresented among the SNPs, pathway enrichment analysis was performed using the Gene Ontology Biological Process (GO BP) and Kyoto Encyclopedia of Genes and Genomes (KEGG) database. Gene annotations from VEP output were used, and Ensembl gene identifiers were extracted. For KEGG analysis, Ensembl IDs were mapped to Entrez gene identifiers using the org.Gg.eg.db package.",
    style = "Normal"
  ) %>%
  body_add_par("") %>%
  body_add_par("GO BP and KEGG over-representation analysis was carried out with clusterProfiler. Statistical significance of enriched pathways was determined with a hypergeometric test, and multiple testing correction was performed using Benjamini-Hochberg.", style = "Normal") %>%
  body_add_break() %>%
  body_add_par("Enrichment Analysis", style = "heading 1") %>%
  body_add_par("Table 1 GO Biological Process enrichment of SNPs", style = "table title") %>%
  body_add_par(paste0(nrow(ego_results), " significantly enriched GO BP pathways found.")) %>%
  body_add_par("") %>%
  body_add_flextable(tab_ego) %>%
  body_add_par("")

# -----------------------------
# PLOTTING FUNCTION
# -----------------------------
plot_figs <- function(enrich_obj, type = "GO Biological Process", typecm = "GO_BP") {
  
  # dotplot
  p <- dotplot(enrich_obj, showCategory = 10, title = type)
  png_file <- paste0(title, "_", typecm, "_dotplot.png")
  ggsave(png_file, plot = p, width = 10, height = 10, dpi = 300)
  
  doc <<- doc %>%
    body_add_par(paste0("Figure: ", typecm, " enrichment dotplot"), style = "table title") %>%
    body_add_img(src = png_file, width = 10, height = 10, style = "centered")
  
  # upsetplot
  p2 <- upsetplot(enrich_obj)
  png_file2 <- paste0(title, "_", typecm, "_upsetplot.png")
  ggsave(png_file2, plot = p2, width = 10, height = 6, dpi = 300)
  
  doc <<- doc %>%
    body_add_par(paste0("Figure: ", type, " enrichment upsetplot"), style = "table title") %>%
    body_add_img(src = png_file2, width = 10, height = 6, style = "centered")
  
  # cnetplot
  p3 <- cnetplot(enrich_obj, showCategory = 10, circular = TRUE, colorEdge = TRUE)
  png_file3 <- paste0(title, "_", typecm, "_cnet.png")
  ggsave(png_file3, plot = p3, width = 10, height = 8, dpi = 300)
  
  doc <<- doc %>%
    body_add_par(paste0("Figure: ", type, " enrichment cnetplot"), style = "table title") %>%
    body_add_img(src = png_file3, width = 10, height = 8, style = "centered")
}

plot_figs(ego, type = "GO Biological Process", typecm = "GO_BP")

# -----------------------------
# KEGG ANALYSIS
# -----------------------------
ekegg <- enrichKEGG(entrez_ids, organism = "gga")
ekeggx <- setReadable(ekegg, 'org.Gg.eg.db', 'ENTREZID')
ekegg_results <- as.data.frame(ekeggx)

format_ekegg_table <- function(df) {
  df %>%
    dplyr::select(
      -c(category, subcategory, zScore, pvalue, qvalue, Count, RichFactor, GeneRatio, BgRatio)
    ) %>%
    mutate(
      geneID = gsub("/", ", ", geneID),
      p.adjust = round(p.adjust, 5),
      FoldEnrichment = round(FoldEnrichment, 3)
    ) %>%
    dplyr::rename(
      `Fold Enrichment` = FoldEnrichment,
      `Gene IDs` = geneID,
      `Adjusted p-value` = p.adjust,
      `KEGG ID` = ID,
      `Biological Process` = Description
    ) %>%
    arrange(desc(`Fold Enrichment`)) %>%
    flextable() %>%
    fontsize(size = 9, part = "header") %>%
    fontsize(size = 9) %>%
    flextable::font(fontname = "Calibri") %>%
    fontsize(j = "Gene IDs", size = 8) %>%
    bold(j = 2) %>%
    flextable::width(j = 1:5, c(0.7, 1.4, 0.9, 0.9, 4))
}

tab_ekegg <- format_ekegg_table(ekegg_results)

if (nrow(ekegg_results) > 0) {
  doc <- doc %>%
    body_add_par("Table 2 KEGG enrichment", style = "table title") %>%
    body_add_flextable(tab_ekegg)
} else {
  doc <- doc %>%
    body_add_par("Table 2 KEGG enrichment", style = "table title") %>%
    body_add_par("") %>%
    body_add_par("No significantly enriched KEGG pathways were identified.")
}

if (nrow(ekegg_results) > 1) {
  plot_figs(ekegg, type = "KEGG", typecm = "EKEGG")
}

# -----------------------------
# SAVE DOCX
# -----------------------------
print(doc, target = paste0(title, "_NGS_report_", today, ".docx"))
