#' @keywords internal
Thornthwaite <-
structure(function (x, fre)
{
    Et <- vector()
    for (j in 1:(length(x)/fre)) {
        for (i in 1:fre) {
            I <- Index1(x, fre, j)
            n <- length(x)
            a <- 0.016 * I + 0.5
            c <- 1.6000000000000001 * (10/I)^a
            Et[(j - 1) * fre + i] <- c * x[(j - 1) * fre + i]^a
        }
    }
    Et * 10
}, source = c("function(x, fre)", "{", "\tEt <- vector()", "\tfor(j in 1:(length(x)/fre)) {",
"\t\tfor(i in 1:fre) {", "\t\t\tI <- Index1(x, fre, j)", "\t\t\tn <- length(x)",
"\t\t\ta <- 0.016 * I + 0.5", "\t\t\tc <- 1.6 * (10/I)^a", "\t\t\tEt[(j - 1) * fre + i] <- c * x[(j - 1) * fre + i]^a",
"\t\t}", "\t}", "\tEt * 10", "}"))
#' @keywords internal
PW <-
function (x, frvol, R, type = "Mi")
{
  fre <- 12
  R1 <- frvol * R * 10
  cat("Maximum water reserve <- ", R1, "\n")
  n <- length(x)
  out1 <- MSeason(FIdrico2(x, fre, R, R1, type),
    fre)
    out <- out1
}
#' @keywords internal
Pe1 <-
structure(function (x, fre)
{
    pe <- vector()
    a1 <- MSeason(x, fre)
    for (j in 1:(length(x)/fre)) {
        for (i in 1:fre) {
            if (mean(a1[, i]) < 10)
                R <- 0.94999999999999996
            if (mean(a1[, i]) < 50 && mean(a1[, i]) > 10)
                R <- 0.90000000000000002
            if (mean(a1[, i]) < 100 && mean(a1[, i]) > 50)
                R <- 0.80000000000000004
            if (mean(a1[, i]) < 150 && mean(a1[, i]) > 100)
                R <- 0.69999999999999996
            if (mean(a1[, i]) < 200 && mean(a1[, i]) > 150)
                R <- 0.59999999999999998
            if (mean(a1[, i]) > 200)
                R <- 0.5
            pe[(j - 1) * fre + i] <- R * x[(j - 1) * fre + i]
        }
    }
    pe
}, source = c("function(x, fre)", "{", "\tpe <- vector()", "\ta1 <- MSeason(x, fre)",
"\tfor(j in 1:(length(x)/fre)) {", "\t\tfor(i in 1:fre) {", "\t\t\tif(mean(a1[, i]) < 10)",
"\t\t\t\tR <- 0.95", "\t\t\tif(mean(a1[, i]) < 50 && mean(a1[, i]) > 10)",
"\t\t\t\tR <- 0.9", "\t\t\tif(mean(a1[, i]) < 100 && mean(a1[, i]) > 50)",
"\t\t\t\tR <- 0.8", "\t\t\tif(mean(a1[, i]) < 150 && mean(a1[, i]) > 100)",
"\t\t\t\tR <- 0.7", "\t\t\tif(mean(a1[, i]) < 200 && mean(a1[, i]) > 150)",
"\t\t\t\tR <- 0.6", "\t\t\tif(mean(a1[, i]) > 200)", "\t\t\t\tR <- 0.5",
"\t\t\tpe[(j - 1) * fre + i] <- R * x[(j - 1) * fre + i]", "\t\t}",
"\t}", "\tpe", "}"))
#' @keywords internal
MSeason <-
structure(function (x, fre)
{
    out <- vector()
    appo <- matrix(0, length(x)/fre, fre)
    for (j in 1:fre) {
        for (i in 1:(length(x)/fre)) appo[i, j] <- x[j + (i -
            1) * fre]
        out[j] <- mean(appo[, j])
    }
    appo
}, source = c("function(x, fre)", "{", "\tout <- vector()", "\tappo <- matrix(0, length(x)/fre, fre)",
"\tfor(j in 1:fre) {", "\t\tfor(i in 1:(length(x)/fre))", "\t\t\tappo[i, j] <- x[j + (i - 1) * fre]",
"\t\tout[j] <- mean(appo[, j])", "\t}", "\tappo", "}"))
#' @keywords internal
Index1 <-
structure(function (ls, fre, o)
{
    I1 <- 0
    for (j in 1:fre) I1 <- I1 + (ls[(o - 1) * fre + j]/5)^1.514
    I1
}, source = c("function(ls, fre, o)", "{", "\tI1 <- 0", "\tfor(j in 1:fre)",
"\t\tI1 <- I1 + (ls[(o - 1) * fre + j]/5)^1.514", "\tI1", "}"
))
#' @keywords internal
FIdrico2 <-
function (y, fre, R, R1, type = "Mi")
{ erb=F
  pe <- vector()
  Fi <- vector()
  fi <- vector()
  ET <- Thornthwaite(y[, 2], fre)
  pe <- Pe1(y[, 1], fre)
  if (erb)
    df <- Df(y[, 3], R, type)
  for (k in 1:nrow(y)) {
    if (erb)
      fi[k] <- (1 - df[k]) * ET[k] - pe[k]
    else fi[k] <- ET[k] - pe[k]
  }
  if (erb == F) {
    for (j in 1:(nrow(y) - 1)) {
      if (fi[j] > 0)
        Fi[j] <- fi[j]
      if (fi[j] <= 0) {
        Fi[j] <- 0
        if (fi[j] < -R1)
          fi[j + 1] <- fi[j + 1] - R1
        else fi[j + 1] <- fi[j + 1] + fi[j]
      }
    }
  }
  if (erb) {
    for (j in 1:(nrow(y) - 1)) {
      if (fi[j] > 0)
        Fi[j] <- fi[j]
      if (fi[j] <= 0) {
        Fi[j] <- 0
        if (fi[j] < -R1)
          fi[j + 1] <- fi[j + 1] - R1 + df[j + 1] *
            ET[j + 1]
        else fi[j + 1] <- fi[j + 1] + fi[j]
      }
    }
  }
  if (fi[nrow(y)] < 0)
    Fi[nrow(y)] <- 0
  else Fi[nrow(y)] <- fi[nrow(y)]
  Fi
}
#' @keywords internal
Df <-
structure(function (x, R, type = "Mi")
{
    if (type == "Mi") {
        n <- length(x)
        Y <- vector()
        for (i in 1:n) {
            if ((R == 40) && (x[i] > 50))
                Y[i] <- -47.823309999999999 + 7.2017680000000004 *
                  x[i] - 0.105087 * x[i]^2 + 0.00040200000000000001 *
                  x[i]^3
            if ((R == 40) && (x[i] <= 50))
                Y[i] <- 100
            if ((R == 40) && (x[i] > 115))
                Y[i] <- 0
            if ((R == 60) && (x[i] > 60))
                Y[i] <- 16.139963999999999 + 3.0600170000000002 *
                  x[i] - 0.032861000000000001 * x[i]^2 + 8.5000000000000006e-05 *
                  x[i]^3
            if ((R == 60) && (x[i] <= 60))
                Y[i] <- 100
            if ((R == 60) && (x[i] > 175))
                Y[i] <- 0
            if ((R == 80) && (x[i] > 80))
                Y[i] <- -46.986013999999997 + 3.7565460000000002 *
                  x[i] - 0.028961000000000001 * x[i]^2 + 5.8999999999999998e-05 *
                  x[i]^3
            if ((R == 80) && (x[i] <= 80))
                Y[i] <- 100
            if ((R == 80) && (x[i] > 230))
                Y[i] <- 0
            if (Y[i] < 0)
                Y[i] <- 0
        }
    }
    if (type == "Sb") {
        n <- length(x)
        Y <- vector()
        for (i in 1:n) {
            if ((R == 40) && (x[i] > 50))
                Y[i] <- -74.014476000000002 + 8.4906310000000005 *
                  x[i] - 0.124338 * x[i]^2 + 0.000484 * x[i]^3
            if ((R == 40) && (x[i] <= 50))
                Y[i] <- 100
            if ((R == 40) && (x[i] > 105))
                Y[i] <- 0
            if ((R == 60) && (x[i] > 60))
                Y[i] <- 2.5763669999999999 + 3.5583939999999998 *
                  x[i] - 0.038304999999999999 * x[i]^2 + 0.000101 *
                  x[i]^3
            if ((R == 60) && (x[i] <= 60))
                Y[i] <- 100
            if ((R == 60) && (x[i] > 170))
                Y[i] <- 0
            if ((R == 80) && (x[i] > 80))
                Y[i] <- -76.410022999999995 + 4.4588479999999997 *
                  x[i] - 0.034051999999999999 * x[i]^2 + 6.9999999999999994e-05 *
                  x[i]^3
            if ((R == 80) && (x[i] <= 80))
                Y[i] <- 100
            if ((R == 80) && (x[i] > 230))
                Y[i] <- 0
            if (Y[i] < 0)
                Y[i] <- 0
        }
    }
    if (type == "Ar") {
        n <- length(x)
        Y <- vector()
        for (i in 1:n) {
            if ((R == 40) && (x[i] > 50))
                Y[i] <- 21.908625000000001 + 4.8351009999999999 *
                  x[i] - 0.082493999999999998 * x[i]^2 + 0.00033700000000000001 *
                  x[i]^3
            if ((R == 40) && (x[i] <= 50))
                Y[i] <- 100
            if ((R == 40) && (x[i] > 110))
                Y[i] <- 0
            if ((R == 60) && (x[i] > 60))
                Y[i] <- 83.104050999999998 + 1.3151740000000001 *
                  x[i] - 0.020656999999999998 * x[i]^2 + 6.0000000000000002e-05 *
                  x[i]^3
            if ((R == 60) && (x[i] <= 60))
                Y[i] <- 100
            if ((R == 60) && (x[i] > 175))
                Y[i] <- 0
            if ((R == 80) && (x[i] > 80))
                Y[i] <- 106.304895 + 0.60338999999999998 * x[i] -
                  0.010817 * x[i]^2 + 2.6999999999999999e-05 *
                  x[i]^3
            if ((R == 80) && (x[i] <= 80))
                Y[i] <- 100
            if ((R == 80) && (x[i] > 220))
                Y[i] <- 0
            if (Y[i] < 0)
                Y[i] <- 0
        }
    }
    for (j in 1:n) if (Y[j] > 100)
        Y[j] = 100
    Y/100
}, source = c("function(x, R, type = \"Mi\")", "{", "\tif(type == \"Mi\") {",
"\t\tn <- length(x)", "\t\tY <- vector()", "\t\tfor(i in 1:n) {",
"\t\t\tif((R == 40) && (x[i] > 50))", "\t\t\t\tY[i] <- -47.82331 + 7.201768 * x[i] - 0.105087 * ",
"\t\t\t\t  x[i]^2 + 0.000402 * x[i]^3", "\t\t\tif((R == 40) && (x[i] <= 50))",
"\t\t\t\tY[i] <- 100", "\t\t\tif((R == 40) && (x[i] > 115))",
"\t\t\t\tY[i] <- 0", "\t\t\tif((R == 60) && (x[i] > 60))", "\t\t\t\tY[i] <- 16.139964 + 3.060017 * x[i] - 0.032861 * ",
"\t\t\t\t  x[i]^2 + 8.5e-005 * x[i]^3", "\t\t\tif((R == 60) && (x[i] <= 60))",
"\t\t\t\tY[i] <- 100", "\t\t\tif((R == 60) && (x[i] > 175))",
"\t\t\t\tY[i] <- 0", "\t\t\tif((R == 80) && (x[i] > 80))", "\t\t\t\tY[i] <- -46.986014 + 3.756546 * x[i] - 0.028961 *",
"\t\t\t\t  x[i]^2 + 5.9e-005 * x[i]^3", "\t\t\tif((R == 80) && (x[i] <= 80))",
"\t\t\t\tY[i] <- 100", "\t\t\tif((R == 80) && (x[i] > 230))",
"\t\t\t\tY[i] <- 0", "\t\t\tif(Y[i] < 0)", "\t\t\t\tY[i] <- 0",
"\t\t}}", "\t\t", "\t", "\tif(type == \"Sb\") {", "\t\tn <- length(x)",
"\t\tY <- vector()", "\t\tfor(i in 1:n) {", "\t\t\tif((R == 40) && (x[i] > 50))",
"\t\t\t\tY[i] <- -74.014476 + 8.490631 * x[i] - 0.124338 *",
"\t\t\t\t  x[i]^2 + 0.000484 * x[i]^3", "\t\t\tif((R == 40) && (x[i] <= 50))",
"\t\t\t\tY[i] <- 100", "\t\t\tif((R == 40) && (x[i] > 105))",
"\t\t\t\tY[i] <- 0", "\t\t\tif((R == 60) && (x[i] > 60))", "\t\t\t\tY[i] <- 2.576367 + 3.558394 * x[i] - 0.038305 * ",
"\t\t\t\t  x[i]^2 + 0.000101 * x[i]^3", "\t\t\tif((R == 60) && (x[i] <= 60))",
"\t\t\t\tY[i] <- 100", "\t\t\tif((R == 60) && (x[i] > 170))",
"\t\t\t\tY[i] <- 0", "\t\t\tif((R == 80) && (x[i] > 80))", "\t\t\t\tY[i] <- -76.410023 + 4.458848 * x[i] - 0.034052 *",
"\t\t\t\t  x[i]^2 + 7e-005 * x[i]^3", "\t\t\tif((R == 80) && (x[i] <= 80))",
"\t\t\t\tY[i] <- 100", "\t\t\tif((R == 80) && (x[i] > 230))",
"\t\t\t\tY[i] <- 0", "\t\t\tif(Y[i] < 0)", "\t\t\t\tY[i] <- 0",
"\t\t}}", "\t\t", "\t", "\tif(type == \"Ar\") {", "\t\tn <- length(x)",
"\t\tY <- vector()", "\t\tfor(i in 1:n) {", "\t\t\tif((R == 40) && (x[i] > 50))",
"\t\t\t\tY[i] <- 21.908625 + 4.835101 * x[i] - 0.082494 * ",
"\t\t\t\t  x[i]^2 + 0.000337 * x[i]^3", "\t\t\tif((R == 40) && (x[i] <= 50))",
"\t\t\t\tY[i] <- 100", "\t\t\tif((R == 40) && (x[i] > 110))",
"\t\t\t\tY[i] <- 0", "\t\t\tif((R == 60) && (x[i] > 60))", "\t\t\t\tY[i] <- 83.104051 + 1.315174 * x[i] - 0.020657 * ",
"\t\t\t\t  x[i]^2 + 6e-005 * x[i]^3", "\t\t\tif((R == 60) && (x[i] <= 60))",
"\t\t\t\tY[i] <- 100", "\t\t\tif((R == 60) && (x[i] > 175))",
"\t\t\t\tY[i] <- 0", "\t\t\tif((R == 80) && (x[i] > 80))", "\t\t\t\tY[i] <- 106.304895 + 0.60339 * x[i] - 0.010817 * ",
"\t\t\t\t  x[i]^2 + 2.7e-005 * x[i]^3", "\t\t\tif((R == 80) && (x[i] <= 80))",
"\t\t\t\tY[i] <- 100", "\t\t\tif((R == 80) && (x[i] > 220))",
"\t\t\t\tY[i] <- 0", "\t\t\tif(Y[i] < 0)", "\t\t\t\tY[i] <- 0",
"\t\t}}", "for (j in 1:n) ", "if (Y[j]>100) Y[j]=100\t", "\t\tY/100",
"\t", "}"))
#' @keywords internal
contaleva2 <-
structure(function (x, a, ser)
{
    y <- Levaneg(a)
    c1 <- Contazeri(x)
    appo2 <- Contazeri(Levaneg(a))
    i <- 0
    while (appo2 > c1) {
        c <- 1 * i
        appo <- a + c
        appo2 <- Contazeri(Levaneg(appo))
        cfinb <- c
        cfina <- c - 1
        i <- i + 1
    }
    cat("series", ser, "\n")
    cat("constant fraquency  about", cfina, "\n")
    cat("and ", cfinb, "\n")
    sappo2 <- Contazeri(Levaneg(a + cfina))
    i <- 0
    while (sappo2 > c1) {
        c <- cfina + 0.10000000000000001 * i
        sappo <- a + c
        sappo2 <- Contazeri(Levaneg(sappo))
        c2finb <- c
        c2fina <- c - 0.10000000000000001
        i <- i + 1
    }
    cat("constant frequency about", c2fina, "\n")
    cat("and", c2finb, "\n")
    stappo2 <- Contazeri(Levaneg(a + c2fina))
    i <- 0
    while (stappo2 > c1) {
        c <- c2fina + 0.01 * i
        stappo <- a + c
        stappo2 <- Contazeri(Levaneg(stappo))
        c3finb <- c
        c3fina <- c + 0.01
        i <- i + 1
    }
    cfinale <- (c3fina + c3finb)/2
    cat("final constant frequency:", cfinale, "\n")
    cfinale <- ifelse((cfinale > 0), cfinale, 0)
    cfinale
}, source = c("function (x, a, ser) ", "{", "    y <- Levaneg(a)",
"    c1 <- Contazeri(x)", "    appo2  <-  Contazeri(Levaneg(a))",
"    i  <-  0", "    while (appo2 > c1) {", "        c <- 1 * i",
"        appo <- a + c", "        appo2 <- Contazeri(Levaneg(appo))",
"        cfinb <- c", "        cfina <- c - 1", "        i  <-  i + 1",
"    }", "    cat(\"series\", ser, \"\\n\")", "    cat(\"constant fraquency  about\", cfina, \"\\n\")",
"    cat(\"and \", cfinb, \"\\n\")", "    sappo2  <-  Contazeri(Levaneg(a + cfina))",
"    i  <-  0", "    while (sappo2 > c1) {", "        c <- cfina + 0.1 * i",
"        sappo <- a + c", "        sappo2 <- Contazeri(Levaneg(sappo))",
"        c2finb <- c", "        c2fina <- c - 0.1", "        i  <-  i + 1",
"    }", "    cat(\"constant frequency about\", c2fina, \"\\n\")",
"    cat(\"and\", c2finb, \"\\n\")", "    stappo2  <-  Contazeri(Levaneg(a + c2fina))",
"    i  <-  0", "    while (stappo2 > c1) {", "        c <- c2fina + 0.01 * i",
"        stappo <- a + c", "        stappo2 <- Contazeri(Levaneg(stappo))",
"        c3finb <- c", "        c3fina <- c + 0.01", "        i  <-  i + 1",
"    }", "    cfinale <- (c3fina + c3finb)/2", "    cat(\"final constant frequency:\", cfinale, \"\\n\")",
"    cfinale  <-  ifelse((cfinale > 0), cfinale, 0)", "    cfinale",
"}"))
#' @keywords internal
sync <-
  function (x)
  {
    k <- nrow(x[[1]])
    h <- length(x)
    out1 <- matrix(0, h * nrow(x[[1]]), ncol(x[[1]]))
    out <- vector()
    for (i in 1:h) out[((k^2) * (i - 1) + 1):(i * (k^2))] <- vec(x[[i]])
    R <- id(k * (k * h))
    for (l in 1:(h * k^2)) R[l, l] <- out[l]
    R
  }
#' @keywords internal
vec <-
  function (a)
  {
    out <- vector()
    k <- nrow(a)
    for (i in 1:k) {
      for (j in 1:ncol(a)) out[i + (j - 1) * k] <- a[i, j]
    }
    out
  }
#' @keywords internal
unvec <-
  function (a, k, p)
  {
    out <- matrix(0, k, k * p)
    for (i in 1:k) {
      for (j in 1:(k * p)) out[i, j] <- a[i + (j - 1) * k]
    }
    out
  }
#' @keywords internal
id <-
function (k)
{
    mat <- matrix(0, k, k)
    for (i in 1:k) {
        for (j in 1:k) mat[i, j] <- ifelse(i == j, 1, 0)
    }
    mat
}
#' @keywords internal
modarma <-
function (x, p, q, passo, n1, n, simul, lsign, graph = F)
{
    m <- arima(x, order = c(p, 0, q), method = "ML")
    stdev = vector()
    for (e in 1:(p + q)) stdev[e] = sqrt(m$var.coef[e, e])
    mAR = vector()
    if (p > 0)
        for (i in 1:p) mAR[i] = m$coef[i]
    else mAR = 0
    mMA = vector()
    if (q > 0)
        for (j in 1:q) mMA[j] = m$coef[p + j]
    else mMA = 0
    m1 <- list(series = c("x1$des"), model = list(order = c(p,
        0, q), ar = c(mAR), ma = c(-mMA)))
    m2 = m$residuals
    if (p == 0) {
        E <- Mainf1(c(m1$model$ma), c(0), n1)
    }
    if (q == 0) {
        E <- Mainf1(c(0), c(m1$model$ar), n1)
    }
    if (p != 0 && q != 0) {
        E <- Mainf1(c(m1$model$ma), c(m1$model$ar), n1)
    }
    np <- p + q
    lung <- length(m2)
    residui <- m2[(p + 1):lung]
    lung <- length(residui)
    lung1 <- lung + 1
    se <- var(residui)
    sm <- var(x)
    Test <- PortBIC(residui, np, se, sm, passo, plot = graph,
        lsign = lsign)
    if (simul) {
        sim1 <- list()
        for (i in 1:n) {
            if (p == 1)
                st <- t(startvalue(var(x), 1, 1))
            else st <- t(startvalue(sigmaZ(x, p), 1, p))
            sima <- sample(residui)
            simb <- sample(residui)
            sim1[[i]] <- filter(c(sima, simb), E, sides = 1)[(lung +
                1):(lung * 2)]
            sim1[[i]] <- c(st, sim1[[i]])
        }
    }
    results <- list()
    results$residui <- residui
    results$BIC <- Test
    results$para <- m1$model
    results$stdev <- stdev
    if (simul)
        results$simulazioni <- sim1
    return(results)
}
#' @keywords internal
modfarma <-
function (serie, pp, qq, passo, n, n1, simul, lsign, graph = F)
{
    a1 <- fracdiff(serie, nar = pp, nma = qq, drange = c(0, 0.5))
    stdev = a1$stderror.dpq
    pietra <- list(parameter = list(d = a1$d, AR = c(a1$ar),
        MA = c(a1$ma)))
    a2 <- farima.err(serie, pietra)
    a2 <- a2[(pp + 1):(length(serie))]
    a3 <- PortBIC(a2, pp + qq + 1, var(a2), var(serie), passo = passo,
        plot = graph, lsign = lsign)
    if (simul) {
        if (pp == 0 && qq == 0)
            E <- Mainf3(c(0), c(0), a1$d, n1)
        if (pp == 0 && qq != 0)
            E <- Mainf3(c(0), c(a1$ma), a1$d, n1)
        if (qq == 0 && pp != 0)
            E <- Mainf3(c(a1$ar), c(0), a1$d, n1)
        if (pp != 0 && qq != 0)
            E <- Mainf3(c(a1$ar), c(a1$ma), a1$d, n1)
        lung <- length(a2)
        lung1 <- lung + 1
        sim1 <- list()
        for (i in 1:n) {
            sima <- sample(a2)
            simb <- sample(a2)
            sim1[[i]] <- filter(c(sima, simb), E, sides = 1)[(lung1 -
                pp):((lung) * 2)]
        }
    }
    results <- list()
    results <- pietra
    results$residui <- a2
    results$BIC <- a3
    results$stdev <- stdev
    if (simul)
        results$simulazioni <- sim1
    return(results)
}
#' @keywords internal
Stagiona <-
function (x, xsim, fre, outer, s.win, svar, sim, par)
{
    x1 <- ts(x, frequency = fre)
    x2 <- stlstu(x1, outer = outer, s.window = s.win, svar.window = svar,
        var = T)
    compm <- x2$seas
    compv <- x2$var
    x3 <- prova1(x, compm, fre)
    x4 <- x3 - mean(x)
    x5 <- prova4(x4, compv, fre)
    if (sim == 1) {
        npr <- length(xsim) + par
        xsim0 <- numeric(npr)
        xsim0[1:(npr - par)] <- xsim
        xsim1 <- prova2(xsim0, compv, fre)
        xsim2 <- prova7(xsim1, compm, fre)
        xsim3 <- xsim2 + mean(x)
    }
    results <- list()
    results$des <- x5
    if (sim == 1)
        results$stag <- xsim3
    results$M <- compm
    results$V <- compv
    return(results)
}
#' @keywords internal
arfimainf <-
function (phi, theta, d, n)
{
    tot <- numeric(n)
    dinf <- coeffd2(170, d)
    arinf <- Arinf(phi, theta, n)
    tot <- Prodpoly1(dinf, arinf)
    tot[1:n]
}
#' @keywords internal
Arinf <-
function (phi, theta, n)
{
    p <- length(phi)
    q <- length(theta)
    if (q == 0) {
        q <- 1
        theta <- numeric(q)
    }
    pi <- numeric(n + q + 1)
    pi[q + 1] <- -1
    for (j in 1:n) {
        appo <- 0
        for (k in 1:q) {
            appo <- appo + theta[k] * pi[j - k + q + 1]
        }
        if (j <= p) {
            pi[j + q + 1] <- phi[j] + appo
        }
        else {
            pi[j + q + 1] <- appo
        }
    }
    -pi[(q + 1):(n + q + 1)]
}
#' @keywords internal
farima.err <-
function (data, pietra)
{
    pr10546 <- filter.frac(data, (pietra$parameter$d))
    pr31 <- list()
    if (is.null(pietra$parameter$AR) == T & is.null(pietra$parameter$MA) ==
        T) {
        return(pr10546)
    }
    assign("pr10546", pr10546)
    pr31$series = pr10546
    if (is.null(pietra$parameter$AR) == T & is.null(pietra$parameter$MA) ==
        F) {
        pr31$model$order <- c(0, 0, (length(pietra$parameter$MA)))
        pr31$model$ma <- pietra$parameter$MA
    }
    if (is.null(pietra$parameter$AR) == F & is.null(pietra$parameter$MA) ==
        T) {
        pr31$model$order <- c((length(pietra$parameter$AR)),
            0, 0)
        pr31$model$ar <- pietra$parameter$AR
    }
    if (is.null(pietra$parameter$AR) == F & is.null(pietra$parameter$MA) ==
        F) {
        pr31$model$order <- c((length(pietra$parameter$AR)),
            0, (length(pietra$parameter$MA)))
        pr31$model$ar <- pietra$parameter$AR
        pr31$model$ma <- pietra$parameter$MA
    }
    pr31$model$ndiff <- 0
    pr32 <- arima(pr31$series, order = pr31$model$order)
    remove("pr10546")
    return(pr32$residuals)
}
#' @keywords internal
Mainf1 <-
function (theta, phi, n1)
{
    q <- length(theta)
    p <- length(phi)
    psi <- numeric(n1 + p + 1)
    psi[p + 1] <- -1
    for (j in 1:n1) {
        appo <- 0
        for (k in 1:p) {
            appo <- appo + phi[k] * psi[j - k + p + 1]
        }
        if (j <= q) {
            psi[j + p + 1] <- theta[j] + appo
        }
        else {
            psi[j + p + 1] <- appo
        }
    }
    -psi[(p + 1):(n1 + p + 1)]
}
#' @keywords internal
Mainf3 <-
function (phi, theta, d, n)
{
    appo1 <- numeric()
    arfimainf <- arfimainf(phi, theta, d, n)
    appo <- length(arfimainf)
    appo1 <- -arfimainf[2:appo]
    Mainf2 <- Mainf1(c(0), appo1, n - 1)
    Mainf2 <- Mainf2[1:n]
    Mainf2
}
#' @keywords internal
prova1 <-
function (x, a, frequency)
{
    n <- length(x)
    N <- n/frequency
    a1 <- numeric(n)
    a2 <- numeric(n)
    for (i in 0:(N - 1)) {
        for (j in 1:frequency) {
            a1[j + i * frequency] <- a[j]
        }
    }
    a2 <- x - a1
    a2
}
#' @keywords internal
prova2 <-
function (x, a, frequency)
{
    n <- length(x)
    x1 <- x
    N <- n/frequency
    a1 <- numeric(n)
    a2 <- numeric(n)
    for (i in 0:(N - 1)) {
        for (j in 1:frequency) {
            a1[j + i * frequency] <- a[j]
        }
    }
    a2 <- x1 * a1
    a2
}
#' @keywords internal
prova4 <-
function (x, a, frequency)
{
    n <- length(x)
    N <- n/frequency
    a1 <- numeric(n)
    a2 <- numeric(n)
    for (i in 0:(N - 1)) {
        for (j in 1:frequency) {
            a1[j + i * frequency] <- a[j]
        }
    }
    a2 <- x/a1
    a2
}
#' @keywords internal
prova7 <-
function (x, a, frequency)
{
    n <- length(x)
    N <- n/frequency
    a1 <- numeric(n)
    a2 <- numeric(n)
    for (i in 0:(N - 1)) {
        for (j in 1:frequency) {
            a1[j + i * frequency] <- a[j]
        }
    }
    a2 <- x + a1
    a2
}
#' @keywords internal
coeffd2 <-
function (n, d)
{
    co <- numeric(n + 1)
    co[1] <- 1
    for (i in 1:n) {
        co[i + 1] <- (((-1)^i) * gamma(d + 1))/(gamma(i + 1) *
            gamma(d - i + 1))
    }
    co
}
#' @keywords internal
filter.frac <-
function (xinput, d)
{
    n <- length(xinput)
    b <- c()
    b <- bd1(n, d)
    x <- xinput
    x <- c(x * 0, x)
    if (n <= (length(b)))
        r <- filter(x, b[1:n], sides = 1)[(n + 1):(2 * n)]
    if (n > (length(b)))
        r <- filter(x, b, sides = 1)[(n + 1):(2 * n)]
    return(r)
}
#' @keywords internal
bd1 <-
  function (n, d)
  {
    if (d == 0) {
      result <- c(1, rep(0, n - 1))
    }
    else {
      result <- c(1, gamma(1:50 - d)/(gamma(1:50 + 1) * gamma(-d)))
      result <- c(result, (51:(n - 1))^(-d - 1)/gamma(-d))
    }
    drop(cbind(result))
  }

#' @keywords internal
Prodpoly1 <-
function (a, b)
{
    p <- length(a)
    q <- length(b)
    a <- c(a, numeric(q - p))
    p <- length(a)
    c <- numeric(p + q - 1)
    b <- b[q:1]
    d <- matrix(0, p, q)
    d <- outer(a, b, "*")
    c <- somdiag(d)
    c
}
#' @keywords internal
stlstu <-
function (data, outer = 0, s.window = 1, svar.window = 1, var = F)
{
    if (is.ts(data) == F) {
        cat("Errore. I dati non sono una serie temporale", fill = T)
        return()
    }
    data1 = data - mean(data)
    pr1 = cycle(data1)
    pr2 = tapply(data1, pr1, mean)
    pr2 = pr2 - mean(pr2)
    pr12 = pr2
    if (outer > 0) {
        for (i in 1:outer) {
            pr3 = rep(pr2, len = length(data1))
            pr4 = data1 - pr3
            pr5 = abs(pr4)
            pr6 = 6 * median(pr5)
            pr6v = 36 * median(pr5)
            pr7 = pr5/pr6
            pr7v = pr5/pr6v
            pr7[pr7 > 1] = 1
            pr7v[pr7v > 1] = 1
            pr8 = (1 - pr7^2)^2
            pr8v = (1 - pr7v^2)^2
            pr9 = data1 * pr8
            pr9v = data1 * pr8v
            pr10 = tapply(pr8, pr1, sum)
            pr10v = tapply(pr8v, pr1, sum)
            pr11 = tapply(pr9, pr1, sum)
            pr12 = pr11/pr10
            pr12 = pr12 - mean(pr12)
            pr2 = pr12
        }
    }
    if (s.window != 1) {
        pr12 = rep(pr12, 3)
        loess.control(trace.hat = "approximate")
        pr12 = loess.smooth(seq(1, length(pr12)), pr12, evaluation = length(pr12),
            span = s.window/length(pr12), )$y
        pr12 = pr12[(length(pr12)/3 + 1):(length(pr12) * 2/3)]
    }
    if (var == T) {
        pv1 = tapply(data1, pr1, var)
        if (outer > 0) {
            pv2 = tapply(pr9v, pr1, mean)
            pv2 = rep(pv2, len = length(data1))
            pv2 = (pr9v - pv2)^2
            pv3 = tapply(pv2, pr1, sum)
            pv4 = pv3/pr10v
            pv11 = (length(data1)/frequency(data1))
            pv1 = pv4 * pv11/(pv11 - 1)
        }
    }
    if (svar.window != 1) {
        pv1 = rep(pv1, 3)
        loess.control(trace.hat = "approximate")
        pv1 = loess.smooth(seq(1, length(pv1)), pv1, evaluation = length(pv1),
            span = svar.window/length(pv1), )$y
        pv1 = pv1[(length(pv1)/3 + 1):(length(pv1) * 2/3)]
    }
    pr13 = list()
    pr13$seas = pr12
    if (outer > 0)
        pr13$we = pr8
    if (var == T)
        pr13$var = sqrt(pv1)
    return(pr13)
}
#' @keywords internal
multides <-
function (x, smean, svar, fre, outer)
{
    k <- ncol(x)
    out <- matrix(0, nrow(x), k)
    for (i in 1:k) {
        stag <- Stagiona(x[, i], x[, i], fre, outer, smean,
            svar, 0, 0)
        out[, i] <- stag$des
    }
    out
}
#' @keywords internal
multipro <-
function (x, y)
{
    k <- ncol(x)
    output <- matrix(0, nrow(x), k)
    for (i in 1:k) output[, i] <- rain.adapt(x[, i], y[, i],
        i)
    output
}
#' @keywords internal
multistag <-
function (xsim, x, smean, svar, fre, outer)
{
    k <- ncol(x)
    out <- matrix(0, nrow(x), k)
    for (i in 1:k) {
        stag <- Stagiona(x[, i], xsim[, i], fre, outer, smean,
            svar, 1, 0)
        out[, i] <- stag$stag
    }
    out
}
#' @keywords internal
permres <-
function (res)
{
    pres <- matrix(0, nrow(res), ncol(res))
    for (t in 1:ncol(res)) {
        pres[, t] <- sample(res[, t])
    }
    pres
}
#' @keywords internal
PortBIC <-
function (x, np, sigmaR, sigmaS, passo = passo, lsign = 0.95,
    plot = F)
{
    n <- length(x)
    a <- numeric(passo)
    if (plot) {
        dev.new()
        plot(acf(x, passo), main = "ACF Residuals")
    }
    a1 <- acf(x, passo, plot = F)
    a <- a1$acf[2:(passo + 1)]
    a2 <- a^2
    Q <- n * sum(a2)
    npar <- np
    g <- passo - npar
    Tchi <- qchisq(lsign, g)
    BIC <- (n - npar) * log((n * sigmaR)/(n - npar)) + npar *
        log((n * (sigmaS - sigmaR))/npar)
    results <- list()
    results$Q <- Q
    results$DistrChi <- Tchi
    results$BIC <- BIC
    return(results)
}
#' @keywords internal
Portmanteau <-
function (res, lag)
{
    ac <- acf(res, lag.max = lag, type = "covariance", plot = F)$acf
    a0 <- solve(ac[1, , ])
    P <- 0
    for (i in 2:(lag + 1)) {
        P <- P + sum(diag(t(ac[i, , ]) %*% a0 %*% ac[i, , ] %*%
            a0))
    }
    dim(res)[1] * P
}
#' @keywords internal
restrict <-
function (x, y, n, prob)
{
    prob = (1 - (1 - prob)/2)
    k <- nrow(x[[1]])
    l <- length(x)
    subset <- list()
    for (p in 1:l) {
        subset[[p]] <- matrix(0, k, k)
        for (i in 1:k) {
            for (j in 1:k) if (l == 1)
                subset[[p]][i, j] <- ifelse(abs(x[[p]][i, j]/y[i,
                  j]) >= qt(prob, (n - k * l - 1)), 1, 0)
            else subset[[p]][i, j] <- ifelse(abs(x[[p]][i, j]/y[[p]][i,
                j]) >= qt(prob, (n - k * l - 1)), 1, 0)
        }
    }
    subset
}
#' @keywords internal
sigmaZ <-
function (x, p)
{
    s <- list()
    u <- list()
    for (i in 0:(p - 1)) s[[i + p]] <- acf(x, lag.max = p, type = "covariance", plot = F)$acf[i + 1, , ]
    for (j in 0:(p - 1)) s[[j + 1]] <- t(s[[2 * p - j - 1]])
    for (g in p:1) {
        u[[g]] <- s[[g]]
        for (e in (g + 1):(g + p - 1)) u[[g]] <- cbind(u[[g]],
            s[[e]])
    }
    y <- u[[p]]
    if (p > 1) {
        for (f in (p - 1):1) y <- rbind(y, u[[f]])
    }
    if (p == 1)
        y <- acf(x, lag.max = p, type = "covariance", plot = F)$acf[1, , ]
    y
}
#' @keywords internal
simvar <-
function (coeff, sigZ, residui, n, b, p, nsim, dipperm)
{
    w <- list()
    for (y in 1:nsim) {
        st <- startvalue(sigZ, b, p)
        sim <- matrix(0, b, n)
        appo <- matrix(0, b, n)
        res <- permres(residui)
        if (dipperm == 1) {
            res <- dippermres(residui)
        }
        for (g in 1:p) {
            sim[, g] <- st[, (p - g + 1)]
        }
        for (j in (p + 1):n) {
            for (i in 1:length(coeff)) {
                appo[, j] <- appo[, j] + coeff[[i]] %*% sim[,
                  j - i]
            }
            sim[, j] <- appo[, j] + t(res)[, j - p]
        }
        w[[y]] <- t(sim)
    }
    w
}
#' @keywords internal
somdiag <-
function (x)
{
    a <- ncol(x)
    n <- a * 2 - 1
    c <- numeric(n)
    c[1] <- x[1, a]
    c[n] <- x[a, 1]
    nu <- (2 * a - 4)/2 + 1
    for (i in 1:(n - 2)) {
        if (i < (nu)) {
            c[i + 1] <- sum(diag(x[1:(i + 1), (a - i):a]))
        }
        if (i == (nu)) {
            c[i + 1] <- sum(diag(x))
        }
        if (i > (nu)) {
            c[i + 1] <- sum(diag(x[(1 + i - nu):a, 1:(a - i +
                nu)]))
        }
    }
    c
}
#' @keywords internal
startvalue <-
function (Corrmatrix, k, p)
{
    appo <- vector()
    norm1 <- matrix(0, k, p)
    h <- t(chol(Corrmatrix)) %*% rnorm(k * p)
    for (i in 1:p) {
        for (j in 1:k) {
            appo[j] <- h[k * (i - 1) + j, 1]
        }
        norm1[, i] <- appo
    }
    norm1
}
#' @keywords internal
varcoeff <-
function (k, p, x1)
{
    varc <- list(k)
    appo <- matrix(0, k, k)
    for (j in 1:p) {
        for (i in 1:k) {
            if (k == 1)
                appo[i, ] <- x1$ar[j]
            else appo[i, ] <- x1$ar[j, , i]
        }
        varc[[j]] <- t(appo)
    }
    varc
}
#' @keywords internal
zero <-
function (k)
{
    mat <- matrix(0, k, k)
    for (i in 1:k) {
        for (j in 1:k) mat[i, j] <- ifelse(i == j, 0, 0)
    }
    mat
}
#' @keywords internal
Corrmatrix <-
function (x, p, lag)
{
    y <- acf(x, lag.max = p, plot = F)$acf[lag + 1, , ]
    y
}
#' @keywords internal
dippermres <-
function (res)
{
    h <- sample(time(res))
    pdres <- matrix(0, nrow(res), ncol(res))
    for (i in 1:nrow(res)) pdres[h[i], ] <- res[i, ]
    pdres
}
#' @keywords internal
count <-
function (x)
{
    c <- 0
    for (i in 1:nrow(x)) c <- c + ifelse(x[i, i] == 1, 1, 0)
    c
}
#' @keywords internal
contaleva <-
function (x, a, ser)
{
    y <- Levaneg(a)
    c1 <- Contazeri(x)
    appo2 = 0
    i = 0
    while (appo2 < c1) {
        c <- 0.1 * i
        appo <- y - c
        appo2 <- Contazeri(Levaneg(appo))
        cfinb <- c
        cfina <- c - 0.1
        i = i + 1
    }
    cat("series", ser, "\n")
    cat("constant fraquency  about", cfina, "\n")
    cat("and ", cfinb, "\n")
    sappo2 = Contazeri(Levaneg(y - cfina))
    i = 0
    while (sappo2 < c1) {
        c <- cfina + 0.01 * i
        sappo <- y - c
        sappo2 <- Contazeri(Levaneg(sappo))
        c2finb <- c
        c2fina <- c - 0.01
        i = i + 1
    }
    cat("constant frequency about", c2fina, "\n")
    cat("and", c2finb, "\n")
    stappo2 = Contazeri(Levaneg(y - c2fina))
    i = 0
    while (stappo2 < c1) {
        c <- c2fina + 0.001 * i
        stappo <- y - c
        stappo2 <- Contazeri(Levaneg(stappo))
        c3finb <- c
        c3fina <- c - 0.001
        i = i + 1
    }
    cfinale <- (c3fina + c3finb)/2
    cat("final constant frequency:", cfinale, "\n")
    cfinale = ifelse((cfinale > 0), cfinale, 0)
    cfinale
}
#' @keywords internal
Levaneg <-
function (x)
{
    n <- length(x)
    for (i in 1:n) {
        if (x[i] < 0)
            x[i] <- 0
    }
    x
}
#' @keywords internal
aggiusta <-
function (x, c)
{
    n <- length(x)
    for (i in 1:n) {
        if (x[i] != 0)
            x[i] <- x[i] + c
    }
    x
}
#' @keywords internal
aggiustam <-
function (x, a, ser)
{
    n <- length(a)
    arid <- numeric()
    x1 <- sum(x)
    a1 <- sum(a)
    conta <- 0
    for (i in 1:n) {
        if (a[i] != 0) {
            conta <- conta + 1
            arid[conta] <- a[i]
        }
    }
    if (a1 < x1) {
        i = 0
        x4 = 0
        while (x4 < x1) {
            c <- 0.1 * i
            x3 <- arid + c
            x4 <- sum(x3)
            cfinb <- c
            cfina <- c - 0.1
            i = i + 1
        }
        cat("series:", ser, "\n")
        cat("constant volume about", cfina, "\n")
        cat("and", cfinb, "\n")
        ii = 0
        x4a = sum(arid + cfina)
        while (x4a < x1) {
            cn <- 1e-06 * ii + cfina
            x3a <- arid + cn
            x4a <- sum(x3a)
            cnfinb <- cn
            cnfina <- cn - 1e-06
            ii = ii + 1
        }
        cvfin <- ((cnfinb + cnfina))/2
    }
    if (a1 > x1) {
        i = 0
        x4 = x1 + 1
        while (x4 > x1) {
            c <- 0.1 * i
            x3 <- arid - c
            x4 <- sum(x3)
            cfinb <- c
            cfina <- c - 0.1
            i = i + 1
        }
        cat("series:", ser, "\n")
        cat("costant volume about", cfina, "\n")
        cat("and", cfinb, "\n")
        ii = 0
        x4a = sum(arid + cfina)
        while (x4a > x1) {
            cn <- 1e-06 * ii + cfina
            x3a <- arid - cn
            x4a <- sum(x3a)
            cnfina <- cn
            cnfinb <- cn + 1e-06
            ii = ii + 1
        }
        cvfin <- ((cnfinb + cnfina))/2
    }
    cat("final constant volume:", cvfin, "\n")
    cvfin
}
#'@keywords internal
Contazeri <-
function (x)
{
    conta <- 0
    for (i in 1:length(x)) {
        if (x[i] == 0)
            conta <- conta + 1
    }
    conta1 <- conta/length(x)
    conta1
}
#' @keywords internal
vb <-
function(s)
{out=c(1)
for (j in 1:s)
  out = c(out,1)
out
}
#'@keywords internal
Hreg <-
function (ev, a, m, b,fre)
{
  out = matrix(0, nrow(ev), ncol(ev))
  for (i in 1:ncol(ev)) out[, i] = a * (fre * i)/(b + fre * i) ^m
  out
}
#'@keywords internal
f3 <-
function (e, x3,fre)
{
  e = as.matrix(e)
  r = nrow(e)
  a = e[1:(r - 1)]
  b = 0
  m = e[r]

  x2 <- Hreg(x3, a, m, b, fre)
  out = f1(x3, x2)
  out
}
#'@keywords internal
f2 <-
function (e, x3,fre)
{
  e <- as.matrix(e)
  r <- nrow(e)
  a <- e[1:(r - 2)]
  m <- e[r-1]
  b <- e[r]

  x2 <- Hreg(x3, a, m, b,fre)
  out <- f1(x3, x2)
  out
}
#'@keywords internal
f1 <-
function (x1, x2)
{
  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)
  v <- 0
  r <- nrow(x1)
  c <- ncol(x1)
  {
    for (j in 1:c) {
      for (i in 1:r) v <- v + (x1[i, j] - x2[i, j])^2
    }
  }
  v
}
#'@keywords internal
eventi2 <-
function (x)
{
  cont = 0
  x=as.matrix(x)
  n = nrow(x)
  out2 = matrix(0, n, ncol(x))
  for (i in 1:n) {
    if (x[i, 1] > 0) {
      cont = cont + 1
      out2[cont, ] = x[i, ]
    }
  }
  out3 <- out2[1:cont, ]
  out4 <- out3
  for (j in 1:ncol(out3))
    out4[,j]=sort(out3[,j],decreasing = T )
  out4
}
#'@keywords internal
eventi <-
function (x, g, s)
{
  n <- nrow(as.matrix(x))
  out <- matrix(0, n, g)
  k = 1
  while (k < (n - g)) {
    if (x[k] > s) {
      out[k, 1] <- x[k]
      for (j in 1:(g - 1)) out[k, j + 1] <- out[k, j] +
          x[k + j]
      k = k + g
    }
    else {
      out[k, 1] <- 0
      k = k + 1
    }
  }
  out
}
#'@keywords internal
Ireg1 <-
  function (ev, a, m, b, t)
  {
    out = matrix(0, nrow(ev), ncol(ev))
    for (i in 1:length(t)) out[, i] = a /(b + t[i])^m
    out
  }
#'@keywords internal
Hreg1 <-
  function (ev, a, m, b, t)
  {
    out = matrix(0, nrow(ev), ncol(ev))
    for (i in 1:length(t)) out[, i] = a * (t[i])/(b + t[i])^m
    out
  }
#'@keywords internal
f4 <-
  function (e, x3, t)
  {
    e = as.matrix(e)
    r = nrow(e)
    a = e[1:(r - 1)]
    b = 0
    m = e[r]
    x2 <- Hreg1(x3, a, m, b, t)
    out = f1(x3, x2)
    out
  }
#'@keywords internal
f5 <-
  function (e, x3, t)
  {
    e <- as.matrix(e)
    r <- nrow(e)
    a <- e[1:(r - 2)]
    m <- e[r - 1]
    b <- e[r]
    x2 <- Hreg1(x3, a, m, b, t)
    out <- f1(x3, x2)
    out
  }
#'@keywords internal
Ireg <-
  function (ev, a, m, b, fre)
  {
    out = matrix(0, nrow(ev), ncol(ev))
    for (i in 1:ncol(ev)) out[, i] = a /(b + fre *
                                           i)^m
    out
  }
#'@keywords internal
annualmax <-
  function (x)
{for (k in 1:ncol(x)) x[,k]=sort(x[,k],decreasing=T)
for (j in 1:nrow(x))
{for (i in 1:(ncol(x)-1))
{if (x[j,i]>x[j,i+1])
{cat("warning","\n")
  x[j,i+1]=x[j,i]}
}
}
x}
