# C2 — community vs ranked sediment chemistry (H1, second half)
Outputs: `C2_dbRDA_results.tsv`, `C2_marginal_tests.tsv`, `C2_envfit.tsv`, `C2_varpart.tsv`, `C2_predictors_ranked_standardised.tsv`, `QC.txt`, `sessionInfo.txt`. Script `scripts/c2.R`.
Response: per-depth summed counts → relative abundance → Hellinger; Bray–Curtis; n = 14. Predictors: ranks of LOI550, LOI950, polysaccharides, proteins, phenols, NH4⁺-N, TKN, sulfate (TP, NO3⁻ excluded), standardised. VIF: **polysaccharides dropped (VIF 11.4)**; all others ≤ 3.4. Depth = midpoint. vegan 2.7.5, 9,999 permutations, seed 42. Partial model and varpart use the ordiR2step-selected chemistry. Same datasets/tables as C3.

## QC
Depth-only Mantel r: V3–V4 0.633 (ms v2 0.60) **PASS**; ONT14 0.476 (ms v2 0.43) **PASS**; V4–V5 0.567, ITS2 0.368 (no reference). permutest F of full models: 1.45, 1.74, 1.22, 1.19.

## Result
Full chemistry model significant only for archaea (V4–V5, adj R² 0.29, p_BH 0.046). **Beyond depth** (Condition(depth)): archaea TKN + NH4⁺-N adj R² 0.22, p_BH 0.0008; ONT14 TKN 0.11, p_BH 0.014; V3–V4 nothing (phenols' effect is wholly shared with depth); ITS2 phenols 0.04, p_BH 0.078. Associations only — one core, bulk dried sediment, ranked values.
