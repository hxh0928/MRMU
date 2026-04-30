test_that("packaged urate-CAD example data are available", {
  data("urate_cad_mrmu_iv", package = "MRMU")
  data("urate_cad_mrmu_background", package = "MRMU")

  expect_equal(urate_cad_mrmu_iv$exposure, "Biomarker_Urate")
  expect_equal(urate_cad_mrmu_iv$confounders, c("Metabolic_SBP", "Metabolic_DBP"))
  expect_equal(urate_cad_mrmu_iv$outcome, "CAD_UKB")
  expect_equal(nrow(urate_cad_mrmu_iv$b.exp), 696)
  expect_equal(ncol(urate_cad_mrmu_iv$b.exp), 3)
  expect_equal(unique(urate_cad_mrmu_iv$Threshold), 5e-5)

  expected_traits <- c("Biomarker_Urate", "Metabolic_SBP", "Metabolic_DBP", "CAD_UKB")
  expect_equal(rownames(urate_cad_mrmu_background$C), expected_traits)
  expect_equal(colnames(urate_cad_mrmu_background$Omega), expected_traits)
})
