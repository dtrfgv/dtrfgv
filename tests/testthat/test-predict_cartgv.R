context("test-predict_cartgv")

test_that("predict_cartgv with phoneme data", {
  data("DataPhoneme_test")
  res_new  <- predict_cartgv(new, tree, carts, coups)$hat.Y
  expect_equal(length(res_new), 311)
})
