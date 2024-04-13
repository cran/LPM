PWN <-
function (x1,frvol,R,p, irr)
{x = PW(x1,frvol,R)
  out <- matrix(0, nrow(x), ncol(x) + 1)
  for (i in 1:ncol(x)) out[, i] <- sort(x[, i])
  for (j in 1:nrow(x)) out[j, ncol(x) + 1] <- (100 * (2 *
    j - 1))/(2 * nrow(x))
  k <- 1
  while (out[k, ncol(x) + 1] < p) k <- k + 1
  out1 <- out[k, ]
  out2 <- t(out1) * 10
  out3 <- (out2 * 1000)/(30 * irr * 3600)
  result <- list()
  result$MatOrd <- out
  result$Values <- out2[, 1:12]
  result$Flow <- out3[, 1:12]
  return(result)
}
