Analysis code
Plasma proteomic waist circumference, visceral fat and cardiometabolic risk at the same measured waist and
body-mass index: a UK Biobank cohort study

WHAT IS HERE
  Every analysis, table and data figure in the manuscript and the Supplementary Material is produced by one of
  these scripts. They read UK Biobank individual-level data, which cannot be redistributed (application 347405), so they are
  provided for inspection and for re-use on an approved UK Biobank extract rather than for direct execution. The
  aggregate output is in the source data (Source_data.xlsx, one sheet per table; the same tables as CSV files in
  Source_data.zip). Paths at the top of each script point to the analysis server and are the only thing that needs
  changing to run them elsewhere.

HOW TO READ THEM
  The numbering follows the order in which the analyses were written, not the order of the manuscript. The list below
  groups the scripts by the part of the paper they produce. Several scripts compute more than the paper reports; the
  header of each script lists everything it computes.

  Construction of proWC (Results, first section; Figure 2; Supplementary Figures 1 and 2; Supplementary Tables 1-3;
  Supplementary Methods and Results, section 1)
    figures/Figure2__proWC_construction.R  the out-of-fold predictions analysed throughout (lasso_WC.csv: nested
                                           ten-fold cross-validation, lambda chosen inside each training fold), the
                                           anchoring regression and proWCdelta, the geographic hold-out (lasso_test.csv:
                                           model trained in England, applied to Scotland and Wales), and Figure 2
    01_oof_prowc.R                         the same out-of-fold derivation repeated as a check (model A), with the
                                           anchoring estimated out of fold (model B), the sex-specific score (model C)
                                           and the age-and-sex-only reference (model D)
    02_export_lasso_coefficients.R         the deposited model, refitted to the full analysis set at lambda = 0.0146
                                           (Supplementary Table 1, reproduced to within 1e-12), and Supplementary Table 2 (its notes
                                           say which accuracy figures follow arithmetically from the others)
    88_figure2_model_performance.R         the accuracy figures printed in Figure 2A-D
    34_reduced_protein_scores.R            reduced scores from the first k proteins to enter the LASSO path
                                           (Supplementary Figure 1F; Supplementary Table 3)
    46_imputation_medians.R                the per-protein medians used for missing values
    52_missingness_block_refit.R           refit without the participants missing the block of 1,459 proteins
    56_batch_holdout_and_complete_panel_model.R   leave-one-batch-out accuracy; the complete-panel model
    66_fold_internal_imputation.R          the whole cross-validation repeated with the medians computed inside each
                                           training fold (run one fold per process, then "combine")
    65c_figure2_panels_source_data.R       source data for Figure 2F-I
    make_SupplementaryFigure1_model_development.R   Supplementary Figure 1
    figures/SupplementaryFigure2__sex_age_BMI_stratification.R   Supplementary Figure 2
    86_protein_names_all.R                 full names of all 2,920 proteins (UniProt, curated for the proteins named in
                                           the paper), filled into the protein tables of the source data

  Visceral fat (Results, second section; Figure 3; Supplementary Tables 4 and 5; Supplementary Methods, section 3, and
  Supplementary Results, section 2)
    74_abdominal_mri.R                     abdominal MRI: linkage checks against waist and DXA, timing, visceral and
                                           subcutaneous fat per SD of proWCdelta, participants of stable size
    75_mri_stress_test.R                   measurement error (visceral share of the extra fat), both visits' size fixed,
                                           form of the adjustment, confounding, sex, selection into imaging
    78_supplementary_table4_mri.R          Supplementary Tables 4 and 5, formatted from the output of the scripts here and of 90
                                           (protein names from 86_protein_names_all.R)
    08_bodycomp_anchoring.R                bioimpedance body composition per SD of proWCdelta on the day of the blood sample
    26_repeat_waist_interval.R, 29_waist_change_later_intervals.R, 87_waist_change_by_quintile.R
                                           change in measured waist after baseline (29 gives the slopes quoted in the
                                           text; 87 the means across quintiles of proWCdelta)
    79_mri_checks.R                        further checks: linkage against the waist measured at the imaging
                                           visit, the size of both visits held fixed, the visceral share within each sex,
                                           partial R2, routine blood measurements, attendance weighting, and the
                                           association of every protein with visceral and subcutaneous fat
    80_protein_axis_checks.R, 81_protein_axis_by_sex.R
                                           the protein profile of proWCdelta against those of visceral and subcutaneous
                                           fat, overall and within each sex
    83_protein_axis_comparator.R           the protein profile of a larger measured waist at the same BMI, as a comparator,
                                           and the split-half reliability of the profiles of visceral and subcutaneous fat
    90_same_visit_mri_and_selection.R      proteins measured at the imaging visits: a model restricted to the imaging-visit
                                           panel, fitted at baseline without those participants, related to MRI fat at the
                                           same visit and compared with baseline proteins, also with the same covariates for
                                           both (section D); the stability of the discordance between the two imaging visits,
                                           and the Olink plates of the two samples; and the imaging and disease estimates without the
                                           participants selected into the proteomics project by consortium members or from
                                           its repeat-imaging study (it first reproduces the published estimates and stops
                                           if any differs); its rows are also in Supplementary Table 4 (78)
    91_tape_error_bound.R                  the correlation of the waist between the two imaging visits at fixed BMI, age
                                           and sex in everyone measured at both visits, which sets the largest correlation
                                           of the discordance that tape error alone could produce
    figures/Figure3__visceral_fat.R        Figure 3

  Cardiometabolic risk (Results, third section; Table 2; Figure 4; Supplementary Methods, sections 2 and 4-7, and
  Supplementary Results, sections 3-5)
    57_primary_incident_table2.R           Table 2: participants free of each endpoint by five sources; standardised
                                           risks at the conditional 10th and 90th percentiles
    76_standardised_risk_ci.R              bootstrap confidence intervals for those risks, their difference and ratio
    figures/Figure4__risk_at_same_size.R   Figure 4
    40_multisource_prevalent.R             the multi-source definition of disease present at baseline
    39_continuous_sensitivity.R, 25_continuous_proWCdelta_adjusted.R
                                           Cox, Fine-Gray, the sex-specific score (derived in 01_oof_prowc.R) and
                                           code-set sensitivity analyses
    68b_five_year_landmark.R               the five-year landmark (replaces section 3 of 68_biochem_nohba1c_splines_landmark.R)
    33_proWCdelta_beyond_fat_and_clinical_markers.R, 59_comparisons_free_of_prevalent.R
                                           adjustment for bioimpedance body fat and routine clinical measurements
    37_kidney_function_sensitivity.R       adjustment for creatinine and cystatin C
    44_confounding_technical_holdout.R, 44b_injury_outcomes.R
                                           lifestyle, technical factors, medication; outcomes without a cardiometabolic
                                           pathway; the disease associations in the geographic hold-out
    28_protein_missingness_sensitivity.R, 47_further_robustness.R, 50_additional_analyses.R,
    55_outcomes_error_targets.R, 61_hf_plate_sparse_sex_mortality.R, 64_reliability_missingness_egfr_overlap.R,
    65_twoprotein_holdout20_batch_bodysize.R, 68_biochem_nohba1c_splines_landmark.R
                                           regression calibration, repeat-measurement reliability, Olink plate, sparser
                                           model, men and women, protein missingness and processing batch, hospital
                                           contact, normoglycaemic participants, form of the body-size adjustment
    53_absolute_risk_calibration.R         standardised risks at the unconditional percentiles of proWCdelta
    47_further_robustness.R, part B        incremental value by cross-validated Cox models (Supplementary Results,
                                           section 5)
    41_incremental_cox.R                   the same comparison in a single 80:20 split (source data only)
    54_covariate_imputation_diagnostics.R, 73_observed_covariates.R
                                           the filled covariate values, and the models refitted with the recorded values

  Comparison with a proteomic BMI (Results, fifth section; Supplementary Table 7; Supplementary Methods, section 1,
  and Supplementary Results, section 7)
    38a_oof_comparator_score.R             a proteomic BMI built exactly as proWC (the same proteins with age and sex,
                                           nested ten-fold cross-validation, out-of-fold predictions, anchoring on the
                                           measured value); run as "Rscript 38a_oof_comparator_score.R BMI" (its other
                                           option, body fat percentage, is not used in the paper)
    89_proteomic_bmi_comparison.R          proWCdelta and the proteomic BMI discordance per SD at fixed WC, BMI, age and
                                           sex, alone and in one model, for abdominal fat on MRI and for incident
                                           disease; it first reproduces the published estimates for proWCdelta and stops
                                           if any differs

  Table 1
    77_table1_characteristics.R            baseline characteristics of the cohort and of the imaging participants

  External comparison in the Guangzhou Nutrition and Health Study (Results, fourth section; Supplementary Methods,
  section 9, and Supplementary Results, section 6)
    These scripts read the published supplementary tables of Cai et al. Cell Reports Medicine 2023;4:101172, together
    with the UK Biobank data used elsewhere here.
    G1_gnhs_describe.R                     first description of the seven published tables
    G2_gnhs_analysis_set.R                 baseline analysis sets of the two cohorts; per-protein associations
    G3_cross_platform_concordance.R        agreement of the protein associations between the two studies
    G4_gnhs_discordance_risk.R             the threshold-flag analysis in the two larger cohorts
    G5_gnhs_continuous_waist.R             the full procedure repeated where the waist is recorded in centimetres
    G6_sample_size_control.R               UK Biobank refitted at the same sample size and protein overlap
    G7_gnhs_metabolic_profile.R            the metabolic measurements at the same measured body size
    G8_accuracy_decomposition.R            why the accuracy differs: sample size, protein panel, spread of the waist
    G9_direct_transport.R                  weights fitted here, applied unchanged to the Chinese measurements
    84_supplementary_table6_gnhs.R         Supplementary Table 6, formatted from the Guangzhou source tables (no new analysis)

  Source data and figure post-processing
    85_source_data_workbook.R              the source data as one Excel workbook (Source_data.xlsx)
    figures/postprocess_Figure2E_axis_ticks.ps, figures/postprocess_Figure2E_axis_titles.ps,
    figures/postprocess_uppercase_panel_tags.pl
                                           Ghostscript and Perl steps applied to the Figure 2 PDF: the tick labels and
                                           axis titles of panel E and upper-case panel tags (no plotted value changes)
  Figure 1 is a schematic drawn in PowerPoint and involves no analysis; the tools that assembled the Supplementary
  Material and the reporting checklists are not included.

SOFTWARE
  sessionInfo.txt gives the R version and the versions of the packages used.

