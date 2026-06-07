suppressPackageStartupMessages({
  library(biomaRt)
})

pick_attr <- function(mart, candidates) {
  attrs <- listAttributes(mart)$name
  hit <- candidates[candidates %in% attrs]
  if (length(hit) == 0) stop("Missing attribute(s): ", paste(candidates, collapse = ", "), call. = FALSE)
  hit[[1]]
}

pick_filter <- function(mart, candidates) {
  flt <- listFilters(mart)$name
  hit <- candidates[candidates %in% flt]
  if (length(hit) == 0) stop("Missing filter(s): ", paste(candidates, collapse = ", "), call. = FALSE)
  hit[[1]]
}

strand_to_symbol <- function(x) {
  ifelse(is.na(x), NA_character_, ifelse(x == 1, "+", ifelse(x == -1, "-", "*")))
}

coord_string <- function(chr, start, end, strand_sym) {
  paste0(chr, ":", start, "-", end, ":", strand_sym)
}

overlaps <- function(a_start, a_end, b_start, b_end) {
  !(is.na(a_start) | is.na(a_end) | is.na(b_start) | is.na(b_end)) &
    (a_end >= b_start) & (a_start <= b_end)
}

safe_getBM <- function(tag, ...) {
  message("getBM: ", tag)
  out <- tryCatch(
    getBM(...),
    error = function(e) {
      stop("getBM failed at [", tag, "]: ", conditionMessage(e), call. = FALSE)
    }
  )
  out
}

fallback_gene_query_map <- function(gene_symbols) {
  fallback_aliases <- list(
    KIR2DL5 = c("KIR2DL5A", "KIR2DL5B"),
    SLAMF3 = "LY9",
    SLAMF5 = "CD84"
  )

  requested_gene_id <- character(0)
  query_symbol <- character(0)

  for (gene_symbol in gene_symbols) {
    aliases <- fallback_aliases[[gene_symbol]]
    if (is.null(aliases)) {
      aliases <- gene_symbol
    }
    requested_gene_id <- c(requested_gene_id, rep(gene_symbol, length(aliases)))
    query_symbol <- c(query_symbol, aliases)
  }

  out <- data.frame(
    requested_gene_id = requested_gene_id,
    query_symbol = query_symbol,
    stringsAsFactors = FALSE
  )
  out <- unique(out)
  out
}

build_biomart_tables_from_queries <- function(mart, gene_query_df, label_prefix = "biomart") {
  if (!is.data.frame(gene_query_df) || nrow(gene_query_df) < 1L) {
    stop("gene_query_df must be a non-empty data.frame.", call. = FALSE)
  }

  A_SYMBOL <- pick_attr(mart, c("hgnc_symbol"))
  F_SYMBOL <- pick_filter(mart, c("hgnc_symbol"))

  A_GENE_CHR <- pick_attr(mart, c("chromosome_name"))
  A_GENE_START <- pick_attr(mart, c("start_position"))
  A_GENE_END <- pick_attr(mart, c("end_position"))
  A_GENE_STRAND <- pick_attr(mart, c("strand"))

  A_TX_ID <- pick_attr(mart, c("ensembl_transcript_id"))
  F_TX_ID <- pick_filter(mart, c("ensembl_transcript_id"))
  A_TX_CHR <- pick_attr(mart, c("chromosome_name"))
  A_TX_START <- pick_attr(mart, c("transcript_start"))
  A_TX_END <- pick_attr(mart, c("transcript_end"))
  A_TX_STRAND <- pick_attr(mart, c("strand"))

  A_EXON_ID <- pick_attr(mart, c("ensembl_exon_id"))
  A_EXON_RANK <- pick_attr(mart, c("rank"))
  A_EXON_START <- pick_attr(mart, c("exon_chrom_start", "exon_genomic_start"))
  A_EXON_END <- pick_attr(mart, c("exon_chrom_end", "exon_genomic_end"))
  A_EXON_CHR <- A_GENE_CHR
  A_EXON_STRAND <- A_GENE_STRAND

  A_GCOD_START <- pick_attr(mart, c("genomic_coding_start"))
  A_GCOD_END <- pick_attr(mart, c("genomic_coding_end"))

  query_symbols <- unique(gene_query_df$query_symbol)

  gene_raw <- safe_getBM(
    tag = paste0(label_prefix, "_gene_raw"),
    attributes = c(A_SYMBOL, A_GENE_CHR, A_GENE_START, A_GENE_END, A_GENE_STRAND),
    filters = F_SYMBOL,
    values = query_symbols,
    mart = mart
  )
  if (nrow(gene_raw) < 1L) {
    return(NULL)
  }

  gene_raw <- merge(
    gene_raw,
    gene_query_df,
    by.x = A_SYMBOL,
    by.y = "query_symbol",
    all.x = FALSE,
    all.y = FALSE
  )
  gene_raw$strand_sym <- strand_to_symbol(gene_raw[[A_GENE_STRAND]])
  gene_tbl <- data.frame(
    gene_id = gene_raw$requested_gene_id,
    gene_coord = coord_string(gene_raw[[A_GENE_CHR]], gene_raw[[A_GENE_START]], gene_raw[[A_GENE_END]], gene_raw$strand_sym),
    stringsAsFactors = FALSE
  )
  gene_tbl <- gene_tbl[!duplicated(gene_tbl[, c("gene_id", "gene_coord")]), , drop = FALSE]

  iso_tx_raw <- safe_getBM(
    tag = paste0(label_prefix, "_iso_tx_raw"),
    attributes = c(A_SYMBOL, A_TX_ID, A_TX_CHR, A_TX_START, A_TX_END, A_TX_STRAND),
    filters = F_SYMBOL,
    values = query_symbols,
    mart = mart
  )
  if (nrow(iso_tx_raw) < 1L) {
    return(list(
      gene_tbl = gene_tbl,
      exon_tbl = gene_tbl[0, c("gene_id"), drop = FALSE],
      iso_tbl = gene_tbl[0, c("gene_id"), drop = FALSE],
      mapped_exons_to_isoforms = gene_tbl[0, c("gene_id"), drop = FALSE],
      found_requested_gene_ids = unique(gene_raw$requested_gene_id)
    ))
  }

  iso_tx_raw <- merge(
    iso_tx_raw,
    gene_query_df,
    by.x = A_SYMBOL,
    by.y = "query_symbol",
    all.x = FALSE,
    all.y = FALSE
  )
  iso_tx_raw$strand_sym <- strand_to_symbol(iso_tx_raw[[A_TX_STRAND]])
  iso_tx_raw$isoform_coord <- coord_string(
    iso_tx_raw[[A_TX_CHR]],
    iso_tx_raw[[A_TX_START]],
    iso_tx_raw[[A_TX_END]],
    iso_tx_raw$strand_sym
  )

  tx_ids <- unique(iso_tx_raw[[A_TX_ID]])
  tx_ids <- tx_ids[!is.na(tx_ids) & tx_ids != ""]
  if (length(tx_ids) < 1L) {
    return(list(
      gene_tbl = gene_tbl,
      exon_tbl = gene_tbl[0, c("gene_id"), drop = FALSE],
      iso_tbl = gene_tbl[0, c("gene_id"), drop = FALSE],
      mapped_exons_to_isoforms = gene_tbl[0, c("gene_id"), drop = FALSE],
      found_requested_gene_ids = unique(gene_raw$requested_gene_id)
    ))
  }

  map_exon_raw <- safe_getBM(
    tag = paste0(label_prefix, "_map_exon_raw"),
    attributes = c(
      A_TX_ID,
      A_EXON_ID, A_EXON_RANK,
      A_EXON_CHR, A_EXON_START, A_EXON_END, A_EXON_STRAND,
      A_GCOD_START, A_GCOD_END
    ),
    filters = F_TX_ID,
    values = tx_ids,
    mart = mart
  )

  map_exon_raw <- merge(
    map_exon_raw,
    iso_tx_raw[, c("requested_gene_id", A_TX_ID, "strand_sym")],
    by = A_TX_ID,
    all.x = TRUE
  )

  map_exon_raw$exon_coord <- coord_string(
    map_exon_raw[[A_EXON_CHR]],
    map_exon_raw[[A_EXON_START]],
    map_exon_raw[[A_EXON_END]],
    strand_to_symbol(map_exon_raw[[A_EXON_STRAND]])
  )

  exon_overlaps_cds <- overlaps(
    a_start = map_exon_raw[[A_EXON_START]],
    a_end = map_exon_raw[[A_EXON_END]],
    b_start = map_exon_raw[[A_GCOD_START]],
    b_end = map_exon_raw[[A_GCOD_END]]
  )

  mapped_exons_to_isoforms <- data.frame(
    gene_id = map_exon_raw$requested_gene_id,
    isoform_id = map_exon_raw[[A_TX_ID]],
    exon_id = map_exon_raw[[A_EXON_ID]],
    exon_rank = suppressWarnings(as.integer(map_exon_raw[[A_EXON_RANK]])),
    exon_coord = map_exon_raw$exon_coord,
    exon_overlaps_cds = exon_overlaps_cds,
    stringsAsFactors = FALSE
  )
  mapped_exons_to_isoforms <- mapped_exons_to_isoforms[
    !is.na(mapped_exons_to_isoforms$gene_id) & mapped_exons_to_isoforms$gene_id != "" &
      !is.na(mapped_exons_to_isoforms$isoform_id) & mapped_exons_to_isoforms$isoform_id != "" &
      !is.na(mapped_exons_to_isoforms$exon_id) & mapped_exons_to_isoforms$exon_id != "" &
      !is.na(mapped_exons_to_isoforms$exon_coord) & mapped_exons_to_isoforms$exon_coord != "",
    ,
    drop = FALSE
  ]
  mapped_exons_to_isoforms <- mapped_exons_to_isoforms[
    !duplicated(mapped_exons_to_isoforms[, c("gene_id", "isoform_id", "exon_id", "exon_rank", "exon_coord")]),
    ,
    drop = FALSE
  ]

  exon_tbl <- mapped_exons_to_isoforms[, c("gene_id", "exon_id", "exon_coord"), drop = FALSE]
  exon_tbl <- exon_tbl[!duplicated(exon_tbl[, c("gene_id", "exon_coord")]), , drop = FALSE]

  cds_start <- tapply(
    map_exon_raw[[A_GCOD_START]],
    map_exon_raw[[A_TX_ID]],
    function(x) if (all(is.na(x))) NA_integer_ else min(x, na.rm = TRUE)
  )
  cds_end <- tapply(
    map_exon_raw[[A_GCOD_END]],
    map_exon_raw[[A_TX_ID]],
    function(x) if (all(is.na(x))) NA_integer_ else max(x, na.rm = TRUE)
  )
  cds_span <- data.frame(
    isoform_id = names(cds_start),
    cds_start = as.integer(cds_start),
    cds_end = as.integer(cds_end),
    stringsAsFactors = FALSE
  )

  iso_raw2 <- merge(
    iso_tx_raw,
    cds_span,
    by.x = A_TX_ID,
    by.y = "isoform_id",
    all.x = TRUE
  )
  cds_coord2 <- ifelse(
    is.na(iso_raw2$cds_start) | is.na(iso_raw2$cds_end),
    NA_character_,
    coord_string(iso_raw2[[A_TX_CHR]], iso_raw2$cds_start, iso_raw2$cds_end, iso_raw2$strand_sym)
  )

  iso_tbl <- data.frame(
    gene_id = iso_raw2$requested_gene_id,
    isoform_id = iso_raw2[[A_TX_ID]],
    cds_coord = cds_coord2,
    isoform_coord = iso_raw2$isoform_coord,
    stringsAsFactors = FALSE
  )
  iso_tbl <- iso_tbl[
    !is.na(iso_tbl$gene_id) & iso_tbl$gene_id != "" &
      !is.na(iso_tbl$isoform_id) & iso_tbl$isoform_id != "",
    ,
    drop = FALSE
  ]
  iso_tbl <- iso_tbl[!duplicated(iso_tbl[, c("gene_id", "isoform_id", "cds_coord", "isoform_coord")]), , drop = FALSE]

  list(
    gene_tbl = gene_tbl,
    exon_tbl = exon_tbl,
    iso_tbl = iso_tbl,
    mapped_exons_to_isoforms = mapped_exons_to_isoforms,
    found_requested_gene_ids = unique(c(gene_tbl$gene_id, exon_tbl$gene_id, iso_tbl$gene_id))
  )
}

export_biomart_tables_hg19 <- function(
    genes_csv_path = "config/genes.csv",
    out_dir = "intermediate/biomart"
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  genes_df <- read.csv(genes_csv_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"gene_id" %in% colnames(genes_df)) {
    stop("Expected column 'gene_id' in ", genes_csv_path, " (HGNC symbols).", call. = FALSE)
  }
  gene_symbols <- trimws(as.character(genes_df$gene_id))
  gene_symbols <- unique(gene_symbols[!is.na(gene_symbols) & gene_symbols != ""])
  if (length(gene_symbols) == 0) stop("No gene symbols found in genes.csv column gene_id.", call. = FALSE)
  
  primary_query_df <- data.frame(
    requested_gene_id = gene_symbols,
    query_symbol = gene_symbols,
    stringsAsFactors = FALSE
  )
  grch37_mart <- useMart(
    biomart = "ensembl",
    dataset = "hsapiens_gene_ensembl",
    host = "https://grch37.ensembl.org"
  )
  primary_tables <- build_biomart_tables_from_queries(
    mart = grch37_mart,
    gene_query_df = primary_query_df,
    label_prefix = "grch37"
  )
  if (is.null(primary_tables)) {
    stop("GRCh37 BioMart query returned no rows for the requested genes.", call. = FALSE)
  }

  gene_tbl <- primary_tables$gene_tbl
  exon_tbl <- primary_tables$exon_tbl
  iso_tbl <- primary_tables$iso_tbl
  mapped_exons_to_isoforms <- primary_tables$mapped_exons_to_isoforms

  returned_symbols <- unique(c(gene_tbl$gene_id, iso_tbl$gene_id, exon_tbl$gene_id))
  missing_symbols <- setdiff(gene_symbols, returned_symbols)

  supplemental_sources <- character(0)
  supplemented_gene_symbols <- character(0)
  if (length(missing_symbols) > 0L) {
    fallback_query_df <- fallback_gene_query_map(missing_symbols)
    current_mart <- useMart(
      biomart = "ensembl",
      dataset = "hsapiens_gene_ensembl",
      host = "https://www.ensembl.org"
    )
    fallback_tables <- build_biomart_tables_from_queries(
      mart = current_mart,
      gene_query_df = fallback_query_df,
      label_prefix = "current_ensembl_fallback"
    )

    if (!is.null(fallback_tables)) {
      supplemented_gene_symbols <- intersect(missing_symbols, fallback_tables$found_requested_gene_ids)
      if (nrow(fallback_tables$gene_tbl) > 0L) {
        gene_tbl <- rbind(gene_tbl, fallback_tables$gene_tbl)
        gene_tbl <- gene_tbl[!duplicated(gene_tbl[, c("gene_id", "gene_coord")]), , drop = FALSE]
      }
      if (nrow(fallback_tables$exon_tbl) > 0L) {
        exon_tbl <- rbind(exon_tbl, fallback_tables$exon_tbl)
        exon_tbl <- exon_tbl[!duplicated(exon_tbl[, c("gene_id", "exon_coord")]), , drop = FALSE]
      }
      if (nrow(fallback_tables$iso_tbl) > 0L) {
        iso_tbl <- rbind(iso_tbl, fallback_tables$iso_tbl)
        iso_tbl <- iso_tbl[!duplicated(iso_tbl[, c("gene_id", "isoform_id", "cds_coord", "isoform_coord")]), , drop = FALSE]
      }
      if (nrow(fallback_tables$mapped_exons_to_isoforms) > 0L) {
        mapped_exons_to_isoforms <- rbind(mapped_exons_to_isoforms, fallback_tables$mapped_exons_to_isoforms)
        mapped_exons_to_isoforms <- mapped_exons_to_isoforms[
          !duplicated(mapped_exons_to_isoforms[, c("gene_id", "isoform_id", "exon_id", "exon_rank", "exon_coord")]),
          ,
          drop = FALSE
        ]
      }
      if (length(supplemented_gene_symbols) > 0L) {
        supplemental_sources <- c(
          supplemental_sources,
          "HGNC REST search endpoints for approved-symbol resolution",
          "Current Ensembl BioMart (https://www.ensembl.org) for missing-gene transcript/exon annotation"
        )
      }
    }

    returned_symbols <- unique(c(gene_tbl$gene_id, iso_tbl$gene_id, exon_tbl$gene_id))
    missing_symbols <- setdiff(gene_symbols, returned_symbols)
  }
  missing_tbl <- data.frame(gene_id = missing_symbols, stringsAsFactors = FALSE)

  NK_genes_info <- gene_tbl
  NK_exons_info <- exon_tbl
  NK_isoform_info <- iso_tbl

  save(NK_genes_info, file = file.path(out_dir, "NK_genes_info.RDA"))
  save(NK_exons_info, file = file.path(out_dir, "NK_exons_info.RDA"))
  save(NK_isoform_info, file = file.path(out_dir, "NK_isoform_info.RDA"))
  save(mapped_exons_to_isoforms, file = file.path(out_dir, "mapped_exons_to_isoforms.RDA"))

  invisible(list(
    gene = NK_genes_info,
    exon = NK_exons_info,
    isoform = NK_isoform_info,
    exon_isoform_map = mapped_exons_to_isoforms,
    missing = missing_tbl,
    supplemented = data.frame(gene_id = supplemented_gene_symbols, stringsAsFactors = FALSE),
    supplemental_sources = unique(supplemental_sources)
  ))
}
