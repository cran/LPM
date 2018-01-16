IDFcurve2 <- function (rain, tc, stvalue1 = 1, stvalue2 = 1, t, Tr = 200, MP = F, Trplot = F) 
{ options(warn = -1)  
  rain <- annualmax(rain)
    rain <- as.matrix(rain)
    nr <- nrow(rain)
    vb1 <- vb(nr - 1)
    Ev3 <- rain
    for (k in 1:ncol(rain)) Ev3[,k]=(Ev3[,k])/(t[k])
    ris <- nlm(f4,c(vb1,1),x3=rain,t,iterlim=1000)
 if (MP) 
    ris <- nlm(f5, c(ris$estimate, stvalue2),x3=rain,t,iterlim=1000)
    mu <- mean(ris$estimate[1:nr]) - 0.45006 * sd(ris$estimate[1:nr])
    sigma <- sd(ris$estimate[1:nr])/1.2825
    a <- qGumbel(1 - 1/Tr, sigma = sigma, mu = mu)
    pr <- pGumbel(ris$estimate[1:nr], sigma = sigma, mu = mu)
    tr <- 1/(1-pr)
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
        aH <- Hreg1(rain, ris$estimate[1:nr], m, b, t)
        iH <- Ireg1(rain, ris$estimate[1:nr], m, b, t)
ts.plot(ts(t(Ev3[1:10, ]), start = 1, end = length(t), frequency = 1),type = "p",  ylab="I[mm/h]",lty = t, xlab="t")
    for (w in 1:10) lines(x = seq(1,length(t),1) , y = iH[w, 
        ], col = "red")
if (Trplot)
legend("topright", paste("Tr.plot =", round(tr[1:10],2), "\n"), text.col = "red")
    h <- a * tc/(b + tc)^m
    i <- a/(b + tc)^m
    out <- list(par = c(a, m, h, i, min), I=Ev3, Curve=iH)
    cat("a(Tr) = ", a, "\n")
     if (MP) 
        cat("m = ", m, "\n")
    else cat("n = ", 1 - m, "\n")
    cat("b = ", b, "\n")
    cat("h(tc) = ", h, "\n")
    cat("i(tc) = ", i, "\n")
    cat("Offset =", min, "\n")
    out
}









