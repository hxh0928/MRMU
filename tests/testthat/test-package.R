test_that("package exports primary functions", {
  expect_true(is.function(MRMU))
  expect_true(is.function(MRMU_Input))
  expect_true(is.function(MV_harmonise))
  expect_true(is.function(est_MRMU_ldsc_paras))
})
