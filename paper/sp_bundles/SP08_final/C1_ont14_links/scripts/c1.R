# C1: link classes (SP03 rules, coverage = alignment over the shorter sequence), ONT14 depth profiles, SP03 comparison, Fig 3.50 input
source(file.path(Sys.getenv("S"), "common/load.R")); out <- file.path(S, "C1_ont14_links"); P3 <- file.path(S, "../SP03_16S_links")
x <- load_ds("ONT14"); m <- x$m[!x$plast, ]; tot <- colSums(m)
D <- by_depth(x); Dra <- sweep(D, 2, colSums(D), "/")
ra_all <- rowSums(m) / sum(m); rnk <- setNames(rank(-ra_all, ties.method = "min"), names(ra_all))
peak <- function(o) DEPTHS[which.max(Dra[o, ])]
h <- fread(file.path(out, "hits/genomes_vs_ont14.tsv"), header = FALSE, col.names = c("q","t","id","alnlen","mism","opens","qcov","tcov","ql","tl"))
h[, genome := sub("\\|.*", "", q)]; h[, copy := sapply(strsplit(q, "\\|"), `[`, 2)]
h[, cov_shorter := ifelse(t == "*", 0, pmin(100, 100 * alnlen / pmin(ql, tl)))]
cls <- function(id, cv) ifelse(id >= 99 & cv >= 95, "exact", ifelse(id >= 97, "same_genus_level", ifelse(id >= 94.5, "related", "no_link")))
h[, class := ifelse(t == "*", "no_hit_>=90%", cls(id, cov_shorter))]
h[, plastid := t %in% names(x$plast)[x$plast]]
h[, taxon := ifelse(t == "*", "", x$fin[t])]
fwrite(h, file.path(out, "C1_all_hits_ont14.tsv"), sep = "\t")
# depth profiles of every linked cluster
lk <- unique(h[class %in% c("exact","same_genus_level","related") & !plastid, t])
idx <- split(colnames(m), x$dep[colnames(m)])
pr <- rbindlist(lapply(lk, function(o) rbindlist(lapply(DEPTHS, function(k) { s <- idx[[k]]; per <- m[o, s] / tot[s]
  data.table(cluster = o, depth_cm = k, n_libraries = length(s), rel_abund_depth = Dra[o, k], min_library = min(per), max_library = max(per), mean_library = mean(per), overall_rank = rnk[o], overall_rel_abund = ra_all[o]) }))))
fwrite(pr, file.path(out, "C1_depth_profiles_ont14.tsv"), sep = "\t")
# best link per genome (highest identity, then coverage)
g <- fread(file.path(P3, "rrna/genomes.tsv")); ord <- c(exact = 1, same_genus_level = 2, related = 3, no_link = 4, "no_hit_>=90%" = 5)
h[, full_copy := ql >= 1200]   # SP03 "|full" tag is unreliable (74-bp fragments are tagged full)
h[, cord := ord[class]]
# best link: full-length copies first (partial barrnap copies only if no full copy hits), then strongest class, then identity
best <- h[!plastid | t == "*"][order(genome, -full_copy, cord, -id, -cov_shorter)][, .SD[1], by = genome]
old <- fread(file.path(P3, "links.tsv"))[dataset == "ONT_full"]; setnames(old, c("OTU","identity","OTU_cov"), c("t","id","tcov"))
old[, idn := suppressWarnings(as.numeric(id))]; oldb <- old[!is.na(idn)][order(genome, -idn)][, .SD[1], by = genome]
op <- fread(file.path(P3, "depth_profiles.tsv"))[dataset == "ONT_full"][, .(old_peak = depth_cm[which.max(rel_abund_mean)]), by = OTU]
rows <- rbindlist(lapply(g$genome, function(gg) {
  b <- best[genome == gg]; o <- oldb[genome == gg]; nc <- g[genome == gg, n_16S]
  linked <- nrow(b) && b$class %in% c("exact","same_genus_level","related")
  ocl <- if (nrow(o)) o$link_class else if (nc == 0) "unlinkable (no 16S)" else "no_hit_>=90%"
  olinked <- ocl %in% c("exact","same_genus_level","related")
  data.table(genome = gg, source = g[genome == gg, source], taxon = g[genome == gg, taxon], n_16S = nc,
    ONT14_cluster = if (nrow(b) && b$t != "*") b$t else "none", identity = if (nrow(b) && b$t != "*") sprintf("%.1f", b$id) else "none",
    aln_len = if (nrow(b) && b$t != "*") b$alnlen else "none", cov_shorter = if (nrow(b) && b$t != "*") sprintf("%.1f", b$cov_shorter) else "none",
    class = if (nc == 0) "unlinkable (no 16S)" else if (nrow(b)) b$class else "none",
    rank = if (linked) rnk[b$t] else "none", peak_depth_14 = if (linked) peak(b$t) else "none",
    ONT14_taxon = if (nrow(b) && b$t != "*") x$fin[b$t] else "none",
    SP03_cluster = if (nrow(o)) o$t else "none", SP03_identity = if (nrow(o)) o$id else "none", SP03_class = ocl,
    SP03_peak_depth = if (olinked) op[OTU == o$t, old_peak] else "none") }))
rows[, changed := ifelse(class != SP03_class | peak_depth_14 != SP03_peak_depth, paste0("yes: ", ifelse(class != SP03_class, sprintf("class %s->%s", SP03_class, class), ""),
     ifelse(class != SP03_class & peak_depth_14 != SP03_peak_depth, "; ", ""), ifelse(peak_depth_14 != SP03_peak_depth, sprintf("peak %s->%s", SP03_peak_depth, peak_depth_14), "")), "no")]
fwrite(rows, file.path(out, "C1_ont14_links.tsv"), sep = "\t")
# Fig 3.50 input: genome x depth, all classes, Illumina from SP03 (unchanged) + new ONT14
il <- fread(file.path(P3, "links.tsv"))[dataset != "ONT_full" & link_class %in% c("exact","same_genus_level","related")]
il[, idn := as.numeric(identity)]; ilb <- il[order(genome, dataset, -idn)][, .SD[1], by = .(genome, dataset)]
ip <- fread(file.path(P3, "depth_profiles.tsv"))
fig <- rbind(
  rbindlist(lapply(seq_len(nrow(ilb)), function(i) { r <- ilb[i]; p <- ip[dataset == r$dataset & OTU == r$OTU]
    data.table(genome = r$genome, dataset = r$dataset, feature = r$OTU, identity = r$identity, link_class = r$link_class, depth_cm = p$depth_cm, rel_abund = p$rel_abund_mean) })),
  rbindlist(lapply(which(rows$class %in% c("exact","same_genus_level","related")), function(i) { r <- rows[i]
    data.table(genome = r$genome, dataset = "ONT14", feature = r$ONT14_cluster, identity = r$identity, link_class = r$class, depth_cm = DEPTHS, rel_abund = Dra[r$ONT14_cluster, ]) })))
fwrite(dcast(fig, genome + dataset + feature + identity + link_class ~ factor(depth_cm, DEPTHS), value.var = "rel_abund"), file.path(out, "C1_Fig3_50_input.tsv"), sep = "\t")
# QC
qc <- c(); ctl <- fread(file.path(out, "hits/control_ArcOtu22_vs_ont14.tsv"), header = FALSE)
qc <- c(qc, sprintf("[%s] Arc_Otu22 full-length gene vs ONT14: %s (must be no hit, no error)", ifelse(all(ctl$V2 == "*"), "PASS", "FAIL"), ifelse(all(ctl$V2 == "*"), "no hit >=90%", paste(ctl$V2, ctl$V3))))
sv <- fread(file.path(out, "hits/self_vs_self.tsv"), header = FALSE); nq <- length(unique(h$q)); ns <- length(unique(sv[V1 == V2 & V3 == 100, V1]))
qc <- c(qc, sprintf("[%s] self-identity: %d of %d genome 16S copies (all isolates + MAGs) hit themselves at 100%%", ifelse(ns == nq, "PASS", "FAIL"), ns, nq))
nl <- ncol(x$m); nd <- length(unique(x$dep))
qc <- c(qc, sprintf("[%s] ONT14: %d libraries, %d depths (Results_new/final_output/cluster_tab_norm.txt; 2 of the 70 Metadata_minion barcodes absent: %s)", ifelse(nl == 68 & nd == 14, "PASS", "FAIL"), nl, nd,
   paste(setdiff(names(ometa), colnames(x$m)), collapse = ",")))
meth <- g[grepl("^Methano|^MTP4", taxon), genome]
top05 <- names(sort(D[, "0-5"], decreasing = TRUE))[1:5]
bad <- h[t %in% top05 & genome %in% meth & class %in% c("exact","same_genus_level")]
mp <- rows[genome %in% meth & peak_depth_14 == "0-5"]
qc <- c(qc, sprintf("[%s] 0-5 cm: no top-5 ONT14 cluster at 0-5 cm links (>=97%%) to a methanogen genome (%s); methanogen-linked clusters peaking at 0-5 cm: %s",
   ifelse(nrow(bad) == 0, "PASS", "FAIL"), paste(top05, x$fin[top05], sep = "=", collapse = " | "), if (nrow(mp)) paste(mp$genome, mp$ONT14_cluster, collapse = ",") else "none"))
writeLines(qc, file.path(out, "QC.txt")); cat(qc, sep = "\n")
print(rows[, .(genome, n_16S, ONT14_cluster, identity, cov_shorter, class, rank, peak_depth_14, SP03_cluster, SP03_class, SP03_peak_depth, changed)])
