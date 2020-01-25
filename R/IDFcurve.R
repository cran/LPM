IDFcurve <- function (rain, g, s, tc, stvalue1 = 1, stvalue2 = fre, fre, Tr = 200,
                      MP = F, Trplot = F)
{ options(warn = -1)
  rain <- as.matrix(rain)
  Y = (nrow(rain)*fre)/(365*24)
  cat("Years = ", Y, "\n")
  Ev1 <- eventi(rain, g, s)
  Ev2 <- eventi2(Ev1)
  Ev3 <- Ev2
  for (k in 1:g) Ev3[,k]=(Ev2[,k])/(fre*k)
  nr <- nrow(Ev2)
  cat("Number of events = ", nr, "\n")
  ne <- nr/Y
  cat("Number of events per Year = ", ne, "\n")
  Tr <- Tr*ne
  vb1 <- vb(nr - 1)
  ris <- nlm(f3, c(vb1, stvalue1), x3 = Ev2, fre = fre, iterlim = 1000)
  if (MP)
    ris <- nlm(f2, c(ris$estimate, stvalue2), x3 = Ev2, fre = fre,
              iterlim = 1000)
  mu <- mean(ris$estimate[1:nr]) - 0.45006 * sd(ris$estimate[1:nr])
  sigma <- sd(ris$estimate[1:nr])/1.2825
  a <- qgumbel(1 - 1/Tr, sigma = sigma, mu = mu)
  pr <- pgumbel(ris$estimate[1:nr], sigma = sigma, mu = mu)
  tr <- 1/(1-pr)
  tr <- tr/ne
  if (MP) {
    b <- ris$estimate[nr + 2]
    m <- ris$estimate[nr + 1]
    min <- ris$minimum
  }
  else {
    b <- 0
    m <- ris$estimate[nr + 1]
    min <- ris$minimum
  }
  aH <- Hreg(Ev2, ris$estimate[1:nr], m, b, fre)
  iH <- Ireg(Ev3, ris$estimate[1:nr], m, b, fre)
  h <- a * tc/(b + tc)^m
  i <- a/(b + tc)^m
  out <- list(par = c(a, m, h, i, min), Curve = iH)
    ts.plot(ts(t(Ev3[1:10, ]), start = fre, end = fre * g, frequency = 1/fre),
          type = "p",  ylab="I[mm/h]", xlab="t[h]")
  for (w in 1:10) lines(x = seq(fre, fre * g, fre), y = iH[w, ], col = "red")
    if (Trplot)
  legend("topright", paste("Tr.plot =", round(tr[1:10],2), "\n"), text.col = "red")
    if (MP)
    cat("m = ", m, "\n")
  else cat("n = ", 1 - m, "\n")
  cat("b = ", b, "\n")
  cat("h(tc) = ", h, "\n")
  cat("i(tc) = ", i, "\n")
  cat("Offset =", min, "\n")
  out
}



















