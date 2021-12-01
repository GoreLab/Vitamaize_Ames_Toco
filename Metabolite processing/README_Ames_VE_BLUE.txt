Folder "1_initial_file_processing_diagnostics_VE_20190903" includes pre-BLUE diagnostics.
Folder "2_VE_blue_functions_boxcox_20190903" includes BLUE calculation analyses.
Folder "3_boxcox_trans_for_hmp_gwas" includes Box-cox transformation analyses.

General comments should be kept in mind:
(1) vita traits were processed using source ID rather than pedigree as unique identifier.
(2) For d2017, 10 lines have PedID duplication issues (measured twice) and repeated measurements from Ames17-BV-19L were removed.
(3) Two mistakes were corrected for Event.Name column based on Laura's comment: This planting error means that B73 was really planted in Range 5, Pass 74, Tier 4 in 2017; H22w was really planted in Range 5 Pass 75 Tier 4 that year.
(4) Using New_pass and New_range infomation to calculate BLUE.
(5) Using new wrapper (s1_version) with finalized model file.
(6) B73 was coded as the same across years.
(7) g.T did not converged using wald.asreml(), and update.asreml() was used to fix this bug.