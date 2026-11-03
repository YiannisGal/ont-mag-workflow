# SP07 B3: community (Hellinger rel. abundance per depth, Bray-Curtis) ~ sediment geochemistry (units-corrected table07).
# Extractions summed per depth, then relative abundance (as SP03/SP04). Archaea/ITS2 depths from the bacterial sample with the same YG number.
suppressMessages({library(vegan); library(data.table)}); set.seed(42); PERM <- 9999
B <- "${DOCS_DIR}"; E <- file.path(B, "Etoliko_138_2")
S <- file.path(B, "Metagenome_analysis_PhD/workspace_2/99_thesis_audit/SP07_final"); out <- file.path(S, "B3_geochem_stats")
DEPTHS <- sprintf("%d-%d", seq(0, 65, 5), seq(5, 70, 5)); MID <- setNames(seq(2.5, 67.5, 5), DEPTHS)
dlab <- function(x) { p <- as.integer(strsplit(trimws(x), "[-_–]")[[1]]); sprintf("%d-%d", p[1], p[2]) }
rtab <- function(p) as.matrix(read.delim(p, row.names = 1, check.names = FALSE, comment.char = "", skip = if (grepl("^# Constructed", readLines(p, 1))) 1 else 0))
bm <- read.delim(file.path(E, "Bacteria/out_3/qplot_input_fixed/Bacteria Metadata.txt"), check.names = FALSE)
yg <- setNames(sapply(bm$Depth, dlab), sub("BA$", "", bm$SampleID))
om <- read.delim(file.path(E, "ONT/final_output/Metadata_minion.txt")); ont <- setNames(sapply(om$Depth, dlab), trimws(om$ClusterID))
DS <- list("Illumina bacteria V3-V4" = list(p = file.path(E, "Bacteria/out_3/otutab_norm.txt"), d = function(s) yg[sub("BA$", "", s)]),
           "Illumina archaea V4-V5"  = list(p = file.path(E, "Archaea/out/otutab_norm.txt"),    d = function(s) yg[sub("AR$", "", s)]),
           "ITS2"                    = list(p = file.path(E, "Eukaryotes/out/otutab_norm.txt"), d = function(s) yg[sub("EY$", "", s)]),
           "ONT bacteria"            = list(p = file.path(E, "ONT/final_output/cluster_tab_trimmed.txt"), d = function(s) ont[s]))
md5 <- sapply(c(sapply(DS, `[[`, "p"), file.path(S, "inputs/table07_sediment_geochemistry_units_corrected.tsv")), tools::md5sum)
fwrite(data.table(md5 = md5, file = names(md5)), file.path(out, "inputs_md5.tsv"), sep = "\t", col.names = FALSE)

g <- read.delim(file.path(S, "inputs/table07_sediment_geochemistry_units_corrected.tsv"), check.names = FALSE)
rownames(g) <- sapply(g[[1]], dlab); stopifnot(identical(rownames(g), DEPTHS))
pick <- function(p) { j <- grep(p, names(g)); stopifnot(length(j) == 1); g[[j]] }
chem <- data.frame(LOI550 = pick("^Organic carbon \\(LOI550"), LOI950 = pick("^LOI950"), Polysaccharides = pick("polysaccharides"),
                   Proteins = pick("proteins"), Phenols = pick("phenols"), NH4N = pick("^NH"), TKN = pick("TKN"), Sulfate = pick("Sulfate"),
                   TP_rank = rank(pick("phosphorus"), ties.method = "average"), row.names = DEPTHS)   # TP relative only -> ranks; NO3 excluded
vif_drop <- function(env) { cur <- names(env); ord <- c()
  repeat { if (length(cur) < 2) break
    v <- sapply(cur, function(k) { r2 <- summary(lm(env[[k]] ~ ., data = env[setdiff(cur, k)]))$r.squared; if (r2 >= 1 - 1e-12) Inf else 1 / (1 - r2) })
    if (max(v) <= 10) break
    dk <- names(which.max(v)); ord <- c(ord, sprintf("%s(%s)", dk, ifelse(is.finite(max(v)), sprintf("%.1f", max(v)), "Inf"))); cur <- setdiff(cur, dk) }
  list(keep = cur, order = ord) }
safe <- function(expr) tryCatch(expr, error = function(e) NULL)
res <- list(); vp <- list(); mg <- list(); ef <- list(); qc <- c(); diag <- c()
for (n in names(DS)) {
  m <- rtab(DS[[n]]$p); dep <- DS[[n]]$d(colnames(m)); if (anyNA(dep)) stop(n, ": samples without depth: ", paste(colnames(m)[is.na(dep)], collapse = ","))
  D <- t(sapply(DEPTHS, function(k) rowSums(m[, dep == k, drop = FALSE]))); D <- D[rowSums(D) > 0, colSums(D) > 0, drop = FALSE]
  dd <- rownames(D); Y <- decostand(sweep(D, 1, rowSums(D), "/"), "hellinger")
  env <- as.data.frame(scale(chem[dd, ])); v <- vif_drop(env); cur <- v$keep
  envR <- env[cur]; envR$depth <- MID[dd]; bc <- vegdist(Y, "bray")
  full <- capscale(as.formula(paste("Y ~", paste(cur, collapse = " + "))), data = envR, distance = "bray")
  an <- safe(anova(full, permutations = PERM)); pt <- safe(permutest(full, permutations = PERM))
  m0 <- capscale(Y ~ 1, data = envR, distance = "bray")
  sel <- safe(ordiR2step(m0, scope = formula(full), R2scope = TRUE, permutations = PERM, trace = FALSE))
  sv <- if (is.null(sel)) character(0) else attr(terms(sel), "term.labels")
  marg <- safe(anova(full, by = "margin", permutations = PERM))
  e <- envfit(scores(m0, display = "sites", choices = 1:2), envR, permutations = PERM)
  partv <- if (length(sv)) sv else cur
  part <- safe(capscale(as.formula(paste("Y ~", paste(partv, collapse = " + "), "+ Condition(depth)")), data = envR, distance = "bray"))
  pa <- if (is.null(part)) NULL else safe(anova(part, permutations = PERM))
  dpm <- capscale(Y ~ depth, data = envR, distance = "bray"); da <- anova(dpm, permutations = PERM)
  man <- mantel(bc, dist(MID[dd]), method = "pearson", permutations = PERM)
  V <- safe(varpart(bc, ~ depth, as.formula(paste("~", paste(partv, collapse = " + "))), data = envR))
  num <- function(x) if (is.null(x) || length(x) == 0 || is.na(x[1])) NA_real_ else x[1]
  res[[n]] <- data.table(dataset = n, n = nrow(Y), depths = paste(dd, collapse = ","), retained_VIF = paste(cur, collapse = ","),
     dropped_VIF_order = if (length(v$order)) paste(v$order, collapse = " > ") else "none",
     selected_ordiR2step = if (length(sv)) paste(sv, collapse = ",") else "none",
     full_adjR2 = num(RsquareAdj(full)$adj.r.squared), full_F = num(an$F), full_p = num(an$`Pr(>F)`),
     sel_adjR2 = if (length(sv)) RsquareAdj(sel)$adj.r.squared else NA, partial_terms = paste(partv, collapse = ","),
     partial_adjR2 = if (is.null(part)) NA else num(RsquareAdj(part)$adj.r.squared), partial_F = num(pa$F), partial_p = num(pa$`Pr(>F)`),
     depth_only_adjR2 = RsquareAdj(dpm)$adj.r.squared, depth_only_F = da$F[1], depth_only_p = da$`Pr(>F)`[1], mantel_r_depth = man$statistic, mantel_p = man$signif)
  if (!is.null(V)) { fr <- V$part$indfract
    vp[[n]] <- data.table(dataset = n, chemistry_set = paste(partv, collapse = ","), depth_only = fr$Adj.R.square[1], geochem_only = fr$Adj.R.square[2], shared = fr$Adj.R.square[3], residual = fr$Adj.R.square[4]) }
  if (!is.null(marg)) mg[[n]] <- data.table(dataset = n, term = rownames(marg), F = marg$F, p = marg$`Pr(>F)`)[term != "Residual"]
  ef[[n]] <- data.table(dataset = n, variable = names(e$vectors$r), r2 = e$vectors$r, p = e$vectors$pvals)
  diag <- c(diag, sprintf("%s: permutest full model F = %s, p = %s; residual df = %s; full-model constrained eigenvalues = %s", n,
     ifelse(is.null(pt), "not estimable", sprintf("%.3f", pt$F.0)), ifelse(is.null(an), "NA", sprintf("%.4f", an$`Pr(>F)`[1])),
     ifelse(is.null(an), "0", an$Df[2]), paste(sprintf("%.3f", head(full$CCA$eig, 5)), collapse = ",")))
  if (n == "Illumina bacteria V3-V4") qc <- c(qc, sprintf("%s: depth-only Mantel r = %.3f (p = %.4f); ms v2 Mantel r ~ 0.60 -> %s", n, man$statistic, man$signif,
     ifelse(man$statistic > 0 && abs(man$statistic - 0.60) <= 0.15, "PASS (same direction, similar size)", "FAIL")))
  else qc <- c(qc, sprintf("%s: depth-only Mantel r = %.3f (p = %.4f) [info]", n, man$statistic, man$signif))
}
r <- rbindlist(res); r[, full_p_BH := p.adjust(full_p, "BH")]; r[, partial_p_BH := p.adjust(partial_p, "BH")]
fwrite(r, file.path(out, "B3_dbRDA_results.tsv"), sep = "\t"); fwrite(rbindlist(vp), file.path(out, "B3_varpart.tsv"), sep = "\t")
fwrite(rbindlist(mg), file.path(out, "B3_marginal_tests.tsv"), sep = "\t"); fwrite(rbindlist(ef), file.path(out, "B3_envfit.tsv"), sep = "\t")
fwrite(data.table(depth = DEPTHS, scale(chem)), file.path(out, "B3_predictors_standardised.tsv"), sep = "\t")
writeLines(c(qc, "", "db-RDA residual diagnostics (permutest on the full model):", diag), file.path(out, "QC.txt"))
print(r, digits = 3); print(rbindlist(vp), digits = 3); cat(qc, diag, sep = "\n")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
