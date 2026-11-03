# C2: community (Hellinger rel. abund. per depth, Bray-Curtis) ~ ranked sediment chemistry, db-RDA family, n = 14 per dataset
source(file.path(Sys.getenv("S"), "common/load.R")); suppressMessages(library(vegan)); set.seed(42); PERM <- 9999
out <- file.path(S, "C2_geochem_stats")
g <- read.delim(file.path(S, "inputs/table07_sediment_geochemistry_units_corrected.tsv"), check.names = FALSE)
g$depth <- sapply(gsub("–", "-", g[[1]]), dlab); stopifnot(identical(g$depth, DEPTHS))
pick <- function(p) { j <- grep(p, names(g)); stopifnot(length(j) == 1); g[[j]] }
chem <- data.frame(LOI550 = pick("^Organic carbon \\(LOI550"), LOI950 = pick("^LOI950"), Polysaccharides = pick("polysaccharides"), Proteins = pick("proteins"),
                   Phenols = pick("phenols"), NH4N = pick("^NH"), TKN = pick("TKN"), Sulfate = pick("Sulfate"), row.names = DEPTHS)   # TP, NO3 excluded
env <- as.data.frame(scale(apply(chem, 2, rank, ties.method = "average"))); rownames(env) <- DEPTHS
# VIF drop, one at a time
vord <- c(); cur <- names(env)
repeat { v <- sapply(cur, function(k) 1 / (1 - summary(lm(env[[k]] ~ ., data = env[setdiff(cur, k)]))$r.squared))
  cat("VIF:", paste(sprintf("%s=%.2f", names(v), v), collapse = " "), "\n"); if (max(v) <= 10) break
  dk <- names(which.max(v)); vord <- c(vord, sprintf("%s(%.1f)", dk, max(v))); cur <- setdiff(cur, dk) }
envR <- env[cur]; envR$depth <- MID
cat("dropped (order):", if (length(vord)) vord else "none", "\n")
res <- list(); vp <- list(); mg <- list(); ef <- list(); qc <- c()
mantel_ms <- c("V3-V4" = 0.60, "ONT14" = 0.43)
for (n in names(DS)) {
  x <- load_ds(n); D <- t(by_depth(x)); D <- D[, colSums(D) > 0]; Y <- decostand(sweep(D, 1, rowSums(D), "/"), "hellinger")
  bc <- vegdist(Y, "bray"); fchem <- as.formula(paste("Y ~", paste(cur, collapse = " + ")))
  full <- capscale(fchem, data = envR, distance = "bray"); pt <- permutest(full, permutations = PERM); an <- anova(full, permutations = PERM)
  m0 <- capscale(Y ~ 1, data = envR, distance = "bray")
  sel <- ordiR2step(m0, scope = formula(full), R2scope = TRUE, permutations = PERM, trace = FALSE)
  sv <- attr(terms(sel), "term.labels"); if (!length(sv)) sv <- character(0)
  marg <- anova(full, by = "margin", permutations = PERM)
  e <- envfit(scores(m0, display = "sites", choices = 1:2), envR, permutations = PERM)
  partv <- if (length(sv)) sv else cur
  part <- capscale(as.formula(paste("Y ~", paste(partv, collapse = " + "), "+ Condition(depth)")), data = envR, distance = "bray")
  pa <- anova(part, permutations = PERM)
  dep <- capscale(Y ~ depth, data = envR, distance = "bray"); da <- anova(dep, permutations = PERM)
  man <- mantel(bc, dist(MID), method = "pearson", permutations = PERM)
  V <- varpart(bc, ~ depth, as.formula(paste("~", paste(partv, collapse = " + "))), data = envR)
  fr <- V$part$indfract
  res[[n]] <- data.table(dataset = n, n = nrow(Y), retained_VIF = paste(cur, collapse = ","), dropped_VIF_order = if (length(vord)) paste(vord, collapse = " > ") else "none",
     selected_ordiR2step = if (length(sv)) paste(sv, collapse = ",") else "none", full_adjR2 = RsquareAdj(full)$adj.r.squared, full_F = an$F[1], full_p = an$`Pr(>F)`[1],
     permutest_F = pt$F.0, sel_adjR2 = if (length(sv)) RsquareAdj(sel)$adj.r.squared else NA, partial_terms = paste(partv, collapse = ","),
     partial_adjR2 = RsquareAdj(part)$adj.r.squared, partial_F = pa$F[1], partial_p = pa$`Pr(>F)`[1],
     depth_only_adjR2 = RsquareAdj(dep)$adj.r.squared, depth_only_F = da$F[1], depth_only_p = da$`Pr(>F)`[1], mantel_r_depth = man$statistic, mantel_p = man$signif)
  vp[[n]] <- data.table(dataset = n, chemistry_set = paste(partv, collapse = ","), depth_only = fr$Adj.R.square[1], chemistry_only = fr$Adj.R.square[2], shared = fr$Adj.R.square[3], residual = fr$Adj.R.square[4])
  mg[[n]] <- data.table(dataset = n, term = rownames(marg), F = marg$F, p = marg$`Pr(>F)`)[term != "Residual"]
  ef[[n]] <- data.table(dataset = n, variable = names(e$vectors$r), r2 = e$vectors$r, p = e$vectors$pvals)
  qc <- c(qc, sprintf("%s: depth-only Mantel r = %.3f (p = %.4f); ms v2 = %s -> %s; depth-only db-RDA adj R2 = %.3f; permutest full F = %.3f",
     n, man$statistic, man$signif, ifelse(n %in% names(mantel_ms), mantel_ms[n], "not given"),
     ifelse(n %in% names(mantel_ms), ifelse(sign(man$statistic) == 1 && abs(man$statistic - mantel_ms[n]) <= 0.15, "PASS (direction and approx. size)", "FAIL"), "INFO"), RsquareAdj(dep)$adj.r.squared, pt$F.0))
}
r <- rbindlist(res); r[, full_p_BH := p.adjust(full_p, "BH")]; r[, partial_p_BH := p.adjust(partial_p, "BH")]
fwrite(r, file.path(out, "C2_dbRDA_results.tsv"), sep = "\t"); fwrite(rbindlist(vp), file.path(out, "C2_varpart.tsv"), sep = "\t")
fwrite(rbindlist(mg), file.path(out, "C2_marginal_tests.tsv"), sep = "\t"); fwrite(rbindlist(ef), file.path(out, "C2_envfit.tsv"), sep = "\t")
fwrite(data.table(depth = DEPTHS, envR), file.path(out, "C2_predictors_ranked_standardised.tsv"), sep = "\t")
writeLines(c(paste("VIF dropped (order):", if (length(vord)) paste(vord, collapse = " > ") else "none"), qc), file.path(out, "QC.txt"))
print(r, digits = 3); print(rbindlist(vp), digits = 3); cat(qc, sep = "\n")
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
