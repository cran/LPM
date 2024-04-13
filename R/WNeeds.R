WNeeds <-
function (x, frvol, R)
{ erb=F
  R1 = frvol * R * 10
  cat("Maximum water reserve <- ", R1, "\n")
  out = MSeason(FIdrico2(x, 12, R, R1), 12)
  out = out*10
}
