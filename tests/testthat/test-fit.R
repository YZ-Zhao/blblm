test_that("multiplication works", {
  fit <- blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
  expect_identical(coef(fit), coef.blblm(fit))
  expect_identical(confint(fit, c("wt", "hp")), confint.blblm(fit, c("wt", "hp")))
  expect_identical(sigma(fit, confidence = TRUE), sigma.blblm(fit, confidence = TRUE))
  expect_identical(predict(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170)), confidence = TRUE), predict.blblm(fit, data.frame(wt = c(2.5, 3), hp = c(150, 170)), confidence = TRUE))
})

