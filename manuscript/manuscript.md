Decomposing the transportability of immune-state biomarkers across checkpoint immunotherapy

## Abstract
Immune-state biomarkers are proposed to predict benefit from checkpoint blockade, yet their performance can change across cohorts. We tested whether transport failure could be decomposed across predictor construction, measurement and clinical outcome. A 322-gene immune-state representation was frozen in IMvigor210 and transported without re-estimation to independent cohorts. External AUCs were 0.504 (95% CI 0.354–0.658) and 0.660 (0.404–0.862), with both intervals spanning 0.5. Calibration showed distinct transport behaviour from ranking: recalibration improved log-loss and Brier score in some cohorts but did not restore discrimination. Across 21 prespecified representation and preprocessing configurations, no configuration established reproducible external transport; correction of a data-mapping error altered the representation but did not restore transport. Joint efficacy–toxicity modelling was not identifiable in the only cohort containing both endpoints. These findings support evaluating biomarker transportability as a property of the predictor–measurement–outcome system rather than of a feature list alone.

Checkpoint blockade produces heterogeneous tumour control and immune-related toxicity across patients. Immune-state biomarkers have therefore been proposed to capture host–tumour interactions that precede treatment outcomes. Yet the same signature can change performance when applied across diseases, cohorts, assays or clinical definitions. This pattern is usually described as lack of reproducibility, but that description collapses several distinct failure modes.¹–⁵

A predictive claim is not transported as a list of features. It requires a defined biological construct, a measurement process, a predictor and a clinical estimand. If the construct changes, the measurement process observes it differently, or the outcome definition changes, the same numerical score need not represent the same predictive quantity. A loss of external discrimination therefore does not by itself establish absence of biology, just as successful discrimination does not establish calibrated risk.

Checkpoint immunotherapy provides a stringent setting in which to test this distinction because available datasets span tumour transcriptomes, peripheral blood, single-cell profiles and immune-repertoire measurements, while outcomes range from response and durable benefit to organ-specific immune toxicity.³–⁶ These measurements are biologically related but are not automatically statistically interchangeable.

Here we ask whether transport failure of an immune-state biomarker can be decomposed into measurable components of predictor construction, measurement, outcome definition and target-population shift. We freeze a baseline immune-state representation before external testing, distinguish discrimination from calibration, compare the frozen predictor with simpler transport strategies, and perturb the representation in a controlled manner. Joint efficacy–toxicity prediction is treated as a separate identifiability question rather than inferred from unrelated small cohorts.

Our aim is not to declare immune-state biology absent. Rather, we seek a stricter empirical criterion for when a biomarker claim has actually crossed a disease, assay or clinical boundary. This framework treats apparent transport failure as a testable prediction-system problem and separates what can be established by the current data from what requires paired and prospective validation.

## Results

## A frozen immune-state representation does not establish reproducible external discrimination

We defined the primary transport estimand as the out-of-sample performance of a frozen baseline representation S0 derived in IMvigor210, without parameter re-estimation in the target cohort. The source representation comprised 322 genes; after a data-provenance audit recovered the gene-row identity of the anchor matrix, 320 of these genes were measurable in the anchor and formed the cross-cohort projection panel. Under the primary frozen transport analysis, external AUCs were 0.504 (95% CI 0.354–0.658) in GSE165252 and 0.660 (0.404–0.862) in GSE78220; both intervals spanned 0.5 (Fig. 2a). Source-swapped leave-one-cohort-out analyses were retained as sensitivity transports rather than as the primary estimand and gave cohort-specific AUCs spanning 0.391–0.660 (Fig. 2b). A descriptive pooled AUC of 0.615 (95% CI 0.515–0.715; I² = 0.387; k = 3) was not used to license transport because only three cohorts contributed and the target cohorts differed substantially in size (Extended Data Fig. 3).

The within-cohort reference was higher in the source cohort than in either small external cohort, but the within-cohort and external estimates were imprecise. We therefore interpret the primary result as failure to establish reproducible external discrimination under the frozen estimand, rather than as evidence that the underlying immune biology is absent. Source-swapped analyses are reported only to characterise sensitivity to the choice of source and do not redefine the primary scientific question.

## Discrimination and calibration show distinct transport behaviour

External discrimination was accompanied by calibration behaviour that could not be reduced to ranking performance (Fig. 3). In IMvigor210, held-out predictions were concentrated around the source event rate; the 5th–95th percentile of predicted probabilities was 0.310–0.474 (median 0.394), whereas the observed event rate was 0.228. Calibration slopes on the logit scale in the source-swapped sensitivity framework were 2.93 in IMvigor210, −0.06 in GSE165252 and 1.87 in GSE78220. Supervised recalibration improved log-loss from 0.583 to 0.511 and Brier score from 0.196 to 0.168 in IMvigor210, but worsened log-loss from 0.688 to 0.832 in GSE78220, without restoring discrimination (Fig. 3c). Recalibration was unstable in GSE165252, with a slope CI of −6.78 to 3.13.

A prevalence-only intercept model provided probability accuracy comparable to or better than the frozen predictor in the target cohorts. This separation of ranking from probability accuracy shows why an external AUC alone is insufficient to define a transport claim: a predictor can retain some ranking information while losing calibration, or calibration can be altered by outcome prevalence and case mix without requiring a change in the underlying biological association.

## Representation sensitivity constrains, but does not identify, the cause of transport failure

We evaluated 21 non-degenerate configurations spanning gene-set choice, normalization, representation dimension (k = 1, 2, 3, 5) and regularization (Fig. 4a). External AUC ranged from 0.418 to 0.602 in GSE165252 and from 0.513 to 0.679 in GSE78220. No tested configuration produced an external confidence interval excluding 0.5 in the beneficial direction under either percentile or bias-corrected-and-accelerated bootstrap intervals (Fig. 4b). The largest variability arose from construct-level choices (mean variance 0.0028), compared with 0.0002 for normalization, 0.0011 for dimension and 0 for regularization. A selection-adjusted permutation null over the same search showed that the observed maximum external AUC did not exceed the 95th percentile of the permutation max-AUC distribution in either target cohort (0.602 observed versus 0.644 permuted for GSE165252; 0.680 versus 0.744 for GSE78220; both P ≥ 0.39; Fig. 4c), so no tested configuration separated from chance once selection across the 21 candidates was accounted for.

A provenance audit also corrected an earlier data-mapping error that had generated an apparent 70-gene all-zero pattern under positional alignment. After gene-row identity recovery, the final 320-gene panel contained zero all-zero genes, and the first principal component explained 26.5% of variance. Restricting the representation to 50- or 100-gene subsets did not restore external transport. The appropriate interpretation is therefore not that one biological mechanism has been disproved, but that the tested representation choices and the corrected mapping do not account for the observed transport boundary.

## Tested transport strategies do not consistently recover the frozen signal

We compared prevalence-only prediction, the frozen immune-state predictor, supervised recalibration, hierarchical modelling and covariance-alignment/location-scale adaptation, using identical held-out patients where applicable (Fig. 5). In the source-swapped sensitivity evaluation, none of the 24 model-by-target combinations met the prespecified dual licensing convention of a beneficial AUC interval excluding 0.5 together with a calibration slope within 0.10 of 1. The closest source-swapped result was obtained by training in GSE78220 and testing in IMvigor210, with AUC 0.640 (95% CI 0.560–0.710) and a calibration slope of 0.43; an AUC-only rule would classify this result differently from the dual criterion. Two further controls bounded the interpretation. The published 18-gene T-cell-inflamed gene expression profile, transported under the identical frozen estimand, met neither criterion in either external target (AUC 0.546 [0.403–0.697] in GSE165252 and 0.571 [0.318–0.788] in GSE78220; Extended Data Fig. 1), and precision simulations at the observed sample sizes and event rates showed that even a true AUC of 0.75–0.80 excludes 0.5 in only a minority of simulated replications at n = 25 and n = 71 (Extended Data Fig. 9). The observed refusals are therefore compatible with limited target-cohort precision and cannot by themselves distinguish weak underlying biology from measurement non-equivalence.

A clinical-only baseline using the harmonizable covariates sex and advanced stage reached an external AUC of 0.487 (0.305–0.698) in GSE78220 and collapsed to a constant prediction in GSE165252, where patient-level covariates were not deposited. In GSE186143, TCR diversity alone produced an AUC of 0.570 versus 0.544 for the combined immune-state and clinical model; log-loss was 0.696 versus 0.767 and Brier score 0.250 versus 0.271. These comparisons support the narrower conclusion that the tested immune-state transport strategy did not consistently add predictive information beyond simpler comparators in the available external datasets.

## Joint efficacy–toxicity prediction remains an identifiability boundary

The only public cohort containing both efficacy and severe-toxicity endpoints was GSE186143 (n = 47). Its joint occupancy was sparse: 10 patients had durable clinical benefit and severe toxicity, 17 had benefit without severe toxicity, 5 had severe toxicity without durable benefit and 15 had neither outcome (Fig. 6a). Fisher’s exact test gave an odds ratio of 1.76 (P = 0.53). A four-category joint model was not reliably identifiable: log-loss was 1.345 versus 1.356 for an independent-model comparator, the joint model won 12 of 25 cross-validation folds, and a sign-permutation test gave P = 0.986 (Fig. 6b,c). The smallest training cell contained four patients. Toxicity-only cohorts were therefore treated as stress tests rather than as evidence for a shared immune-state mechanism.

## Discussion

Our analysis shows that transportability of an immune-state biomarker cannot be inferred from discrimination alone. In independent checkpoint-immunotherapy cohorts, a frozen immune-state representation did not establish reproducible external discrimination, while calibration displayed partially distinct behaviour. Recalibration improved probability accuracy in some settings without restoring ranking, indicating that a failed transport claim is not a single statistical event.

The controlled representation analysis further shows why post-hoc biomarker repair can be misleading. Correcting the earlier gene-row mapping error removed the apparent 70-gene all-zero pattern, but did not restore external transport. Likewise, varying gene subsets, normalization, dimensionality and regularization did not yield a consistently transportable predictor. These observations do not identify one universal biological cause; instead, they constrain the explanation space. They also demonstrate that provenance verification is a prerequisite for any transport claim.

The central conceptual implication is that transportability should be evaluated as a property of a predictor–measurement–outcome system. A predictor may retain biological relevance while losing calibration as outcome prevalence or case mix changes. Conversely, a model may show discrimination while its probability estimates cease to represent the same clinical risk. Measurement changes can add distortions that are invisible to AUC-only validation. The response-prevalence differences across the three efficacy cohorts and their differing clinical composition therefore belong to the transport problem itself, not merely to the Discussion of limitations.

The present analysis does not establish a shared immune-state mechanism linking efficacy and toxicity. The only cohort with both endpoints is too small for a stable joint model, and juxtaposing unrelated toxicity cohorts cannot substitute for a same-patient design. The next decisive experiment is therefore a paired bridge study followed by independent prospective validation: the same baseline construct should be measured with two platforms in the same patients, treatment-induced change should be sampled longitudinally, and response and severe toxicity should be defined prospectively. This design would allow measurement error, baseline state and treatment-induced displacement to be estimated jointly rather than inferred from separate datasets.

Several limitations define the boundary of the conclusions. The efficacy analysis is anchored in one frozen source context; public cohorts differ in disease, regimen, sampling and assay; the adaptation methods tested here are not exhaustive; and the joint efficacy–toxicity cohort is underpowered. Most importantly, the present data cannot distinguish weak underlying biology from construct or measurement non-equivalence with certainty. We therefore do not claim that immune-state biomarkers are intrinsically non-transportable.

The practical consequence is a stricter standard for a transport claim. A biomarker should be considered transportable only when its construct, measurement and clinical estimand remain sufficiently equivalent, or when changes in these components are explicitly modelled and then validated in an independent population. The appropriate next step is not to assume a sharp separation between efficacy and toxicity, but to allow continuous shared-state and boundary models to compete in adequately powered data. The resulting framework is falsifiable: a well-powered, frozen predictor that preserves both discrimination and calibration would refute a broad non-transportability claim; conversely, failure after controlled measurement harmonization and prospective validation would strengthen the case that the limitation is intrinsic to the biological predictor itself.

## Methods

## Study design and datasets

We analysed public checkpoint-immunotherapy datasets spanning tumour RNA, peripheral blood, single-cell profiles and immune-repertoire measurements. The main efficacy datasets comprised IMvigor210 (n = 298 labelled of 348), GSE165252 (n = 71) and GSE78220 (n = 25). Additional datasets were used as stress tests: GSE166181 (n = 20), GSE186143 (n = 47 for efficacy; 53 for toxicity), GSE189125 (n = 13), GSE216329 (n = 24 patients; 48 samples) and GSE253720 (n = 12). GSE287541 was withheld from the primary transport analysis because its estimand was a within-patient time × group interaction rather than a baseline cross-cohort outcome. Cohort provenance, endpoint definitions and accession-level metadata are stored in the cohort ledger.

## Gene-row identity recovery

The IMvigor210 expression matrix contained 31,085 genes × 348 samples with anonymised row labels. The row order was found to correspond to the alphabetical order of gene symbols in the fData annotation table, rather than to the identifier order assumed in the earlier pipeline. The audit was triggered by discordant sex-anchor profiles and was verified using sex-linked expression anchors, cross-cohort gene-mean rank correlations and file-level provenance hashes. The recovered mapping was exported as a machine-readable row key. Under the corrected mapping, 320 of the 322 source genes were measurable in the anchor dataset, with no all-zero genes in the final 320-gene panel. The previously reported 70-gene all-zero pattern was therefore treated as a data-mapping error and is not used as biological evidence.

## Primary transport estimand

The primary estimand was the out-of-sample performance of P(Y = 1 | S0), where S0 is the baseline immune-state representation learned in IMvigor210 and then frozen. The 322-gene source representation was transported using the 320-gene cross-cohort intersection; target labels were not used to re-estimate the primary predictor. Source-swapped leave-one-cohort-out analyses were treated as sensitivity analyses because changing the source cohort changes the predictor and therefore changes the estimand.

## Immune-state representation

The frozen representation used centered log-ratio-transformed expression followed by principal-component representation, with k = 2 components in the primary analysis. Under the recovered gene identity, PC1 explained 26.5% of variance and PC2 7.3%. The complete preprocessing and gene lists are stored in the analysis manifest.

## Discrimination and calibration

AUC was calculated on held-out patients. Confidence intervals were estimated using 2,000-replicate bootstrap procedures, with the bias-corrected-and-accelerated interval used as the primary interval and percentile intervals reported as sensitivity analyses. Calibration-in-the-large, calibration slope, Brier score and log-loss were calculated from the same held-out predictions. Supervised recalibration was treated as target adaptation and was separated from primary external validation.

## Representation and measurement sensitivity

We evaluated 21 non-degenerate configurations varying gene-set choice, normalization, representation dimension (k = 1, 2, 3, 5) and regularization. The corrected 320-gene mapping was used as the source of truth. Measurable-only and zero-handling analyses were re-specified after gene-row identity recovery. All configurations were preserved in machine-readable form for audit and reproduction.

## Model competition

Tested strategies included prevalence-only prediction, frozen immune-state prediction, supervised recalibration, hierarchical modelling, covariance alignment/location-scale adaptation and a clinical-only baseline using the maximal harmonizable covariate set common to IMvigor210 and GSE78220 (sex and advanced stage). These strategies are interpreted as tested comparators rather than an exhaustive benchmark of all domain-adaptation methods.

## Joint efficacy–toxicity analysis

Efficacy and severe toxicity were jointly evaluated in GSE186143 using the four outcome cells defined by durable clinical benefit and severe toxicity. Fisher’s exact test, held-out predictive metrics and sign-permutation analysis were used. No mechanistic conclusion was drawn from the small joint cohort.

## Reproducibility and analysis freeze

The analysis manifest specifies sample and patient identifiers, exclusions, transformations, model configurations, random seeds, software versions, endpoint definitions and validation freezes. The corrected IMvigor210 row key and provenance checks are versioned with the analysis package; the final source-to-figure mapping is documented in the manifest.

## Data availability

All analysed datasets are publicly available through their source repositories. The corrected IMvigor210 row key and provenance documentation are released with the analysis package.

## Code availability

The analysis code, manifest and environment specification will be provided in a versioned public repository; the permanent repository identifier will be included in the final submission package.

## References

1. 1. Ayers M, et al. IFN-γ-related mRNA profile predicts clinical response to PD-1 blockade. J Clin Invest. 2017;127:2930–2940. https://doi.org/10.1172/JCI91190

2. 2. Cristescu R, et al. Pan-tumor genomic biomarkers for PD-1 checkpoint blockade-based immunotherapy. Science. 2018;362:eaar3599. https://pubmed.ncbi.nlm.nih.gov/30309915/

3. 3. Lozano AX, et al. T cell characteristics associated with toxicity to immune checkpoint blockade in patients with melanoma. Nat Med. 2022;28:353–362. https://doi.org/10.1038/s41591-021-01623-z

4. 4. Sung C, et al. Integrative analysis of risk factors for immune-related adverse events of checkpoint blockade therapy in cancer. Nat Cancer. 2023;4:844–859. https://doi.org/10.1038/s43018-023-00572-5

5. 5. Dimitriou F, et al. A targetable type III immune response with increase of IL-17A expressing CD4+ T cells is associated with immunotherapy-induced toxicity in melanoma. Nat Cancer. 2024;5:1390–1408. https://doi.org/10.1038/s43018-024-00810-4

6. 6. Ma S-C, et al. CRTAM inhibition mitigates toxicity of immune checkpoint inhibitors without antitumor efficacy trade-off. Nat Cancer. 2026;7:650–666. https://doi.org/10.1038/s43018-026-01135-0

## Figure legends

Figure 1 | Transportability is a property of a prediction system. Schematic defining the transport object as construct → measurement → predictor → clinical estimand → target population. Distinguish primary frozen transport from source-swapped sensitivity analyses. The figure should make clear that transportability is assessed separately through discrimination and calibration, and that failure can arise at multiple system levels.

Figure 2 | A frozen immune-state representation does not establish reproducible external discrimination. a, Primary frozen IMvigor210-to-target AUCs with 95% CIs for GSE165252 and GSE78220, with AUC = 0.5 reference. b, Source-swapped sensitivity AUCs shown separately, explicitly labelled as sensitivity analyses. Avoid pooled diamonds in the main panel; move the descriptive pooled estimate to Extended Data.

Figure 3 | Discrimination and calibration show distinct transport behaviour. a, ROC curves for primary frozen predictions. b, Calibration plots or calibration-in-the-large/slope summaries using the identical held-out predictions. c, Before/after supervised recalibration for log-loss and Brier score, with AUC shown separately to demonstrate that probability-scale improvement does not restore ranking.

Figure 4 | Representation and provenance perturbation constrain the explanation space. a, Pre-specified 21-configuration design matrix. b, External AUC distributions with bootstrap intervals and a 0.5 reference; display the full search rather than a selected best model. c, Corrected gene-row mapping and zero-all-zero-gene result, with the prior 70-gene pattern explicitly labelled as a data-mapping error. Any selection-adjusted permutation null should be added here after P0 calculation. d, Selection-adjusted permutation null for the 21-configuration search (2,000 anchor-label permutations, paired across configurations): observed maximum external AUC versus the permutation max-AUC distribution, with the 95th percentile marked.

Figure 5 | Tested transport strategies do not consistently recover the frozen signal. a, Model-by-target AUC and calibration metrics using the same held-out patients. b, Probability-scale metrics for prevalence-only, frozen and adapted strategies. c, Clinical-only comparison and source-swapped near-miss example. Do not rank methods as winners; show the multidimensional trade-off.

Figure 6 | Joint efficacy–toxicity prediction remains an identifiability boundary. a, Four-cell occupancy in GSE186143. b, Cross-validated joint versus independent-model log-loss/Brier summary. c, Sign-permutation null for fold-level differences. Keep the figure framed as an identifiability limit and design implication, not a mechanism-discovery figure.

## Extended Data legends

Extended Data Fig. 1 | Positive-control biomarker transport. Primary positive control using the published 18-gene T-cell-inflamed GEP after completion of the P0 analysis. Six-gene IFN-γ panel is sensitivity only.

Extended Data Fig. 2 | Cohort provenance and estimand ledger. Disease, regimen, sample type, assay, endpoint, event prevalence, available covariates and primary/sensitivity role for every cohort.

Extended Data Fig. 3 | Sensitivity of pooled external summaries. Forest plot of cohort-specific AUCs, heterogeneity and descriptive pooled estimates; explicitly note k = 3 and differing cohort sizes.

Extended Data Fig. 4 | Full 21-configuration representation audit. Complete machine-readable configuration matrix with AUC, Brier, log-loss, calibration slope and bootstrap intervals.

Extended Data Fig. 5 | Gene-row mapping audit. Sex-anchor diagnostics, corrected row key, cross-cohort gene-mean rank correlations and provenance checks supporting the 320-gene source of truth.

Extended Data Fig. 6 | Threshold sensitivity of the transport criterion. AUC-only, calibration-only and joint decision conventions using slope tolerances of ±0.05, ±0.10 and ±0.20 plus calibration-in-the-large.

Extended Data Fig. 7 | Case-mix and outcome-shift audit. Response prevalence, stage composition, regimen and harmonizable covariates across efficacy cohorts, with standardized differences where definitions permit.

Extended Data Fig. 8 | Toxicity-only stress tests. Results from GSE189125, GSE216329 and GSE253720 with explicit estimand and denominator limitations.

Extended Data Fig. 9 | External-validation precision and simulation. Operating characteristics for the observed cohort sizes and event rates; show the relationship between true AUC/calibration and probability of excluding AUC = 0.5 under the specified analysis plan.

Extended Data Fig. 10 | Reproducibility manifest. Full analysis provenance, frozen seeds, software environment, source-file hashes, source-to-figure mapping and repository identifiers.

## Supplementary Information

Supplementary Note 1 | Cohort provenance and endpoint harmonisation Full cohort ledger, inclusion/exclusion logic, assay and sampling metadata, endpoint definitions and disease/regimen boundaries.

Supplementary Note 2 | Statistical estimands and transport criteria Primary frozen estimand, sensitivity estimands, discrimination/calibration definitions and rationale for reporting continuous transport metrics alongside the dual licensing convention.

Supplementary Note 3 | Representation audit Complete 21-configuration grid, gene lists, transformations, regularisation settings and all held-out predictions.

Supplementary Note 4 | Gene-row identity recovery Detailed audit trail for the IMvigor210 mapping correction, provenance hashes, sex-anchor tests and correction of the former 70-gene data-mapping artefact.

Supplementary Note 5 | Model competition Full model specifications, covariance-alignment implementation, recalibration procedure and clinical-only comparator.

Supplementary Note 6 | Joint efficacy–toxicity analysis Cell counts, identifiability assessment, permutation procedure and sensitivity analyses.

Supplementary Tables S1–S8 S1 cohort ledger; S2 gene-row key; S3 primary and sensitivity AUCs; S4 calibration metrics; S5 representation configurations; S6 model competition; S7 case-mix audit; S8 simulation and power/precision scenarios.

Source Data One source-data file per main and Extended Data figure where feasible, with immutable row/column definitions and analysis-manifest linkage.
