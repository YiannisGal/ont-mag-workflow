source(file.path(Sys.getenv("S"), "common/load.R")); suppressMessages(library(Kendall)); set.seed(42)
out <- file.path(S, "C3_unassigned_depth"); res <- list(); qc <- c()
DPANN <- "p__(Nanoarchaeota|Micrarchaeota|Aenigmarchaeota|Altiarchaeota|Iainarchaeota|Diapherotrites|Huberarchaeota|Nanohaloarchaeota|Woesearchaeota|Pacearchaeota|Undinarchaeota|Parvarchaeota|DPANN)"
permsp <- function(x, y, B = 9999) { r <- cor(x, y, method = "spearman"); rp <- replicate(B, cor(sample(x), y, method = "spearman")); c(r, (sum(abs(rp) >= abs(r) - 1e-12) + 1) / (B + 1)) }
prof <- list()
for (n in names(DS)) {
  x <- load_ds(n); D <- by_depth(x); tot <- colSums(D)
  sets <- list(`phylum-unassigned` = !hasphy(x$fin), uOTU = x$t1 == "Unassigned" | is.na(x$t1))
  if (n == "V4-V5") sets$`DPANN share` <- grepl(DPANN, x$fin)
  np <- sum(!x$plast)
  qc <- c(qc, sprintf("%s: table rows %d, OTUs with taxonomy %d, plastid OTUs excluded %d, libraries %d, depths %d (per-depth n: %s), per-depth totals min %d max %d; pass1=%s, pass2=%s",
     n, x$n_tab_rows, nrow(x$m), sum(x$plast), ncol(x$m), length(unique(x$dep)), paste(table(factor(x$dep, DEPTHS)), collapse = ","), min(tot), max(tot), DS[[n]]$p1lab, DS[[n]]$p2lab))
  for (s in names(sets)) {
    w <- sets[[s]][!x$plast]; pct <- 100 * colSums(D[w, , drop = FALSE]) / tot
    overall <- 100 * sum(D[w, ]) / sum(D); notu <- sum(w)
    sp <- permsp(MID, pct); ex <- suppressWarnings(cor.test(MID, pct, method = "spearman", exact = TRUE)$p.value); mk <- MannKendall(pct)
    prof[[paste(n, s)]] <- data.table(dataset = n, measure = s, depth = DEPTHS, mid = MID, pct = pct)
    res[[length(res) + 1]] <- data.table(dataset = n, measure = s, n_depths = 14, n_OTUs_in_set = notu, n_OTUs_total = np, overall_read_pct = overall,
      min_pct = min(pct), min_depth = DEPTHS[which.min(pct)], max_pct = max(pct), max_depth = DEPTHS[which.max(pct)],
      spearman_rho = sp[1], p_perm9999 = sp[2], p_exact = ex, MK_tau = as.numeric(mk$tau), MK_p = as.numeric(mk$sl))
  }
}
r <- rbindlist(res); r[, p_BH := p.adjust(p_perm9999, "BH")]; r[, MK_p_BH := p.adjust(MK_p, "BH")]
fwrite(r, file.path(out, "C3_unassigned_vs_depth.tsv"), sep = "\t"); fwrite(rbindlist(prof), file.path(out, "C3_per_depth_profiles.tsv"), sep = "\t")
writeLines(qc, file.path(out, "QC_inputs.txt")); print(r, digits = 3); cat(qc, sep = "\n")
