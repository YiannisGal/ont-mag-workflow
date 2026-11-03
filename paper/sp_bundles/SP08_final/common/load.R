# Shared loaders for SP08: per-dataset count table, two-pass taxonomy, sample->depth (from metadata files, never column order)
suppressMessages(library(data.table))
B  <- "${DOCS_DIR}"; E <- file.path(B, "Etoliko_138_2"); O <- file.path(B, "Etoliko_new")
S  <- file.path(B, "Metagenome_analysis_PhD/workspace_2/99_thesis_audit/SP08_final")
QZ <- file.path(S, "inputs/qza")
DEPTHS <- sprintf("%d-%d", seq(0, 65, 5), seq(5, 70, 5)); MID <- seq(2.5, 67.5, 5)
dlab <- function(x) { p <- as.integer(strsplit(trimws(x), "[-_]")[[1]]); sprintf("%d-%d", p[1], p[2]) }
rtax <- function(p) { t <- fread(p, sep = "\t", header = TRUE, quote = ""); setNames(t[[2]], t[[1]]) }
rtab <- function(p) { m <- as.matrix(read.delim(p, row.names = 1, check.names = FALSE, comment.char = "", skip = if (grepl("^# Constructed", readLines(p, 1))) 1 else 0)); m }
bmeta <- local({ m <- read.delim(file.path(E, "Bacteria/out_3/qplot_input_fixed/Bacteria Metadata.txt"), check.names = FALSE)
                 setNames(sapply(m$Depth, dlab), sub("BA$", "", m$SampleID)) })          # YGnn -> depth
ymap <- function(ids, suf) { k <- sub(paste0(suf, "$"), "", ids); d <- bmeta[k]; if (anyNA(d)) stop("no depth for ", paste(ids[is.na(d)], collapse = ",")); setNames(d, ids) }
ometa <- local({ m <- read.delim(file.path(O, "Metadata_minion.tsv")); setNames(sapply(m$Depth, dlab), trimws(m$ClusterID)) })
DS <- list(
  "V3-V4" = list(tab = file.path(E, "Bacteria/out_3/otutab_trim.txt"),
                 p1 = file.path(QZ, "Bacteria_out_3_taxonomy_blast_silva_138_2.tsv"), p2 = file.path(QZ, "Bacteria_out_3_taxonomy_blast_silva_138_2_75.tsv"),
                 p1lab = "SILVA 138.2, 95%", p2lab = "SILVA 138.2, 75%", meta = function(s) ymap(s, "BA"), fa = file.path(E, "Bacteria/out_3/otus_trimmed.fa")),
  "V4-V5" = list(tab = file.path(E, "Archaea/out/otutab_trim.txt"),
                 p1 = file.path(QZ, "Archaea_out_taxonomy_blast_silva.tsv"), p2 = file.path(QZ, "Archaea_out_taxonomy_blast_silva_138_2_75.tsv"),
                 p1lab = "SILVA (RESCRIPt, pre-138.2), 95% [no 138.2 95% run exists]", p2lab = "SILVA 138.2, 75%", meta = function(s) ymap(s, "AR"), fa = file.path(E, "Archaea/out/otus_trimmed.fa")),
  "ITS2"  = list(tab = file.path(E, "Eukaryotes/out/otutab.txt"),
                 p1 = file.path(E, "Eukaryotes/out/taxonomy_new/taxonomy.tsv"), p2 = file.path(E, "Eukaryotes/out/taxonomy_blast_75/taxonomy_75.tsv"),
                 p1lab = "UNITE_new, 95% [no 97% run exists]", p2lab = "UNITE, 75%", meta = function(s) ymap(s, "EY"), fa = file.path(E, "Eukaryotes/out/otus_trimmed.fa")),
  "ONT14" = list(tab = file.path(O, "Results_new/final_output/cluster_tab_norm.txt"),
                 p1 = file.path(O, "Results_new/taxonomy_blast_138_2_95/taxonomy.tsv"), p2 = file.path(O, "Results_new/taxonomy_blast/taxonomy.tsv"),
                 p1lab = "SILVA 138.2, 95%", p2lab = "SILVA 138, 90% (only second pass on file)", meta = function(s) { d <- ometa[s]; if (anyNA(d)) stop("ONT no depth"); setNames(d, s) },
                 fa = file.path(O, "Results_new/final_output/cluster_seq_trimmed.fa")))
hasphy <- function(t) grepl("p__[A-Za-z]", t) & !grepl("p__(unidentified|Incertae)", t)
load_ds <- function(n) {
  d <- DS[[n]]; m <- rtab(d$tab); t1 <- rtax(d$p1); t2 <- rtax(d$p2)
  keep <- intersect(rownames(m), names(t1)); m <- m[keep, , drop = FALSE]
  t1 <- t1[keep]; t2 <- t2[keep]; t2[is.na(t2)] <- "Unassigned"
  fin <- ifelse(t1 != "Unassigned" & !is.na(t1), t1, t2)       # pass 1 if assigned, else pass 2
  plast <- grepl("o__Chloroplast", t1) | grepl("o__Chloroplast", t2)
  dep <- d$meta(colnames(m))
  list(name = n, m = m, t1 = t1, t2 = t2, fin = fin, plast = plast, dep = dep, n_tab_rows = nrow(rtab(d$tab)))
}
by_depth <- function(x, w = rep(TRUE, nrow(x$m))) {   # summed counts per depth, plastids removed
  m <- x$m[!x$plast & w, , drop = FALSE]; sapply(DEPTHS, function(k) rowSums(m[, x$dep == k, drop = FALSE]))
}
