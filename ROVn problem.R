# Load libraries
library(tidyverse)
# Define functions
# ROV <- function(K, r, tau, s, a){
#   d <- (a - 1) * log(K) / (s * sqrt(tau))
#   rov <- exp(-(r - s^2 / 2) * tau) * K^a * pnorm(d + s * sqrt(tau)) - exp(-r * tau) * K * pnorm(d)
#   return(rov)
# }
# ROVstar <- function(Kstar, r, tau, s, a, b){
#   dStar <- (log(b) + (a - 1) * log(Kstar)) / (s * sqrt(tau))
#   rovStar <- exp(-r * tau) * (1 / a - 1) * Kstar * pnorm(dStar)
#   return(rovStar)
# }
# rootFn <- function(Kstar, r, tau, s, a, b){
#   dStar <- (log(b) + (a - 1) * log(Kstar)) / (s * sqrt(tau))
#   dROVdK <- exp(-(r - s^2 / 2) * tau) * a * b * Kstar^a * pnorm(dStar + s * sqrt(tau)) - exp(-r * tau) * Kstar * pnorm(dStar)
#   return(dROVdK)
# }
get_di <- function(K_im1, f_im1, tau_i, m, s, eta_im10 = 1){
  d_i <- (log(f_im1 / K_im1) + (m - eta_im10^2 * s^2 / 2) * tau_i) / (s * eta_im10 * sqrt(tau_i))
  return(d_i)
}
ROVi <- function(K_im1, f_im1, tau_i, m, s, eta_im10 = 1){
  d_i <- get_di(K_im1, f_im1, tau_i, m, s, eta_im10)
  u <- d_i + eta_im10 * s * sqrt(tau_i)
  rov <- f_im1 * pnorm(u) - exp(-r * tau_i) * K_im1 * pnorm(d_i)
  probSuc <- pnorm(d_i)
  return(c(rov, probSuc))
}
get_dfidfim1 <- function(K_im1, f_im1, f_i, tau_i, m, s, eta_im10 = 1){
  d_i <- get_di(K_im1, f_im1, tau_i, m, s, eta_im10)
  u <- d_i + eta_im10 * s * sqrt(tau_i)
  if(eta_im10 == 1){dEtaim10dfim1 <- 0}else{dEtaim10dfim1 <- -eta_im10 / f_i}
  dfidfim1 <- pnorm(u) + dEtaim10dfim1 * s * sqrt(tau_i) * f_im1 * dnorm(u)
  return(dfidfim1)
}
ROVn <- function(X0, Kvec, tauVec, m, s){
  f_im1 <- X0; eta_im10 <- 1
  n <- length(Kvec)
  fiOut <- c(); etaOut <- c(); probSucVec <- c()
  for(i in 1:n){
    tau_i <- tauVec[i]; K_im1 <- Kvec[i]
    outROVi <- ROVi(K_im1, f_im1, tau_i, m, s, eta_im10)
    f_i  <- outROVi[1]; probSuc <- outROVi[2]
    dfidfim1 <- get_dfidfim1(K_im1, f_im1, f_i, tau_i, m, s, eta_im10)
    eta_i0 <- f_im1 / f_i * dfidfim1 * eta_im10
    fiOut[i] <- f_i; etaOut[i] <- eta_i0; probSucVec[i] <- probSuc
    f_im1 <- f_i; eta_im10 <- eta_i0
  }
  return(data.frame(f_i = fiOut, eta_i0 = etaOut, probSuc = probSucVec))
}
#===============================================================================
# Define discount rate
rYrly_discrete <- 0.10 #0.035
rQtly_discrete <- (1 + rYrly_discrete)^(1 / 4) - 1
rYrly <- round(log(1 + rYrly_discrete), 3)
r <- round(log(1 + rQtly_discrete), 3)
# Stage durations - last stage first
# periodDurations <- c(3, 4, 2, 3) * 4
# tauVec <- rev(cumsum(periodDurations))
Tvec <- c(3, 0.85, 1.7, 3.4) * 4
tauVec <- rev(cumsum(Tvec))
# Calculate project NPV(Tn) i.e. X0
t <- seq(0, 50, length.out = 50)
alpha <- 0.7; A <- 1
NPVTn <- exp(-rYrly * tauVec[1]) * A * (-expint::gammainc(1 + alpha, rYrly * t) / rYrly^(1 + alpha) + expint::gammainc(1 + alpha, 0) / rYrly^(1 + alpha))
dNVT0 <- A * t^alpha
dNPVT0 <- exp(-rYrly * t) * dNVT0
plot(t, dNVT0)
plot(t, dNPVT0)
plot(t, NPVTn)
X0 <- NPVTn[25] * 10^6
#-------------------------------------------------------------------------------
# Define project NPV volatility:
cv <- 2
s <- cv * r * sqrt(tauVec[1])
# Or derive it from elicited info:
# Define confidence interval of interest, usually 95% (z = 1.96)
# zBound <- 1.96
# # Max and min log return elicited from experts/stakeholders
# # 1) Elicit upper and lower pctg error in project NPV
# # 2) Interpret as upper and lower bounds on 95% conf interval of arithmetic return
# # 3) Convert arithmetic return to log return by formula log ret = log(1 + arith ret)
# yMaxA <- 1.15#1.85
# yMinA <- -0.8#-0.45
# yMax <- log(1 + yMaxA)
# yMin <- log(1 + yMinA)
# # bTerm <- -(2 * yMin + 8 / 3 * zBound^2)
# # yMax <- 1 / 2 * (-bTerm - sqrt(bTerm^2 + 4))
# # Back out s and m
# s <- (yMax - yMin) / (2 * zBound * sqrt(tauVec[1]))
# s / (r * sqrt(tauVec[1]))
#-------------------------------------------------------------------------------
stage1Invest <- 530000
Kvec <- c(350000, 181000, 312000, 336000) # Last stage first
#exp(tauVec[1] * r) * X0
m <- r
dfROVn <- ROVn(X0, Kvec, tauVec, m, s)
f_n <- dfROVn$f_i[nrow(dfROVn)]
f_n - stage1Invest
# NPV
X0 - sum(exp(-r * tauVec) * Kvec) - stage1Invest
#===============================================================================
# "Manual" check step by step
eta_00 <- 1; f0 <- X0
outROVi <- ROVi(Kvec[1], f0, tauVec[1], m, s, eta_00) # Value of option to launch
f1 <- outROVi[1]; probSuc1 <- outROVi[2]
df1df0 <- get_dfidfim1(Kvec[1], X0, f1, tauVec[1], m, s, eta_00)
#pnorm(get_di(Kvec[1], X0, tauVec[1], m, s, eta_00) + eta_00 * s * sqrt(tauVec[1]))
eta_10 <- X0 / f1 * df1df0 * eta_00
outROVi <- ROVi(Kvec[2], f1, tauVec[2], m, s, eta_10) # Value of option on stage 4 research
f2 <- outROVi[1]; probSuc2 <- outROVi[2]
df2df1 <- get_dfidfim1(Kvec[2], f1, f2, tauVec[2], m, s, eta_10)
eta_20 <- f1 / f2 * df2df1 * eta_10
outROVi <- ROVi(Kvec[3], f2, tauVec[3], m, s, eta_20) # Value of option on stage 3 research
f3 <- outROVi[1]; probSuc3 <- outROVi[2]
df3df2 <- get_dfidfim1(Kvec[3], f2, f3, tauVec[3], m, s, eta_20)
eta_30 <- f2 / f3 * df3df2 * eta_20
outROVi <- ROVi(Kvec[4], f3, tauVec[4], m, s, eta_30) # Value of option on stage 2 research
f4 <- outROVi[1]; probSuc4 <- outROVi[2]
df4df3 <- get_dfidfim1(Kvec[4], f3, f4, tauVec[4], m, s, eta_30)
eta_40 <- f3 / f4 * df4df3 * eta_30
#--------------------------------------------------------------------------
dfROVnCheck <- data.frame(f_i = c(f1, f2, f3, f4), eta_i0 = c(eta_10, eta_20, eta_30, eta_40), probSuc = c(probSuc1, probSuc2, probSuc3, probSuc4))
dfROVn - dfROVnCheck
#===============================================================================

























#==========================================================================
#==========================================================================
#==========================================================================
# Scraps
# b <- 3
# a <- 0.75
# thisInt <- c(0, 10^7)
# K0star <- rootSolve::uniroot.all(rootFn, interval = thisInt,
#                                  lower = min(thisInt), upper = max(thisInt),
#                                  r = r, tau = tauVec[1], s = s, a = a, b = b)
# K0star <- K0star[2]
#X0 <- b * K0star^a
# K0 <- 100
#Kvec <- c(1, 0.1, 0.5, 0.4) * K0
# Evaluate f, fStar, and compare, solve for Kstar
f0 <- ROV(K0, X0, r, tau1, m, s)
f0StarEnv <- ROVstar(K0, r, tau1, s, a, b)
thisInt <- c(0, 10^4)
K0star <- rootSolve::uniroot.all(rootFn, interval = thisInt,
                       lower = min(thisInt), upper = max(thisInt),
                       r = r, tau = tau1, s = s, a = a, b = b)
Kstar <- Kstar[2]
dfPlot <- data.frame(K0, f0, f0StarEnv) %>% gather(type, val, f0:f0StarEnv)
gg <- ggplot(dfPlot, aes( x = K0, y = val, group = type, color = type))
gg <- gg + geom_line()
gg <- gg + geom_vline(xintercept = Kstar)
gg
# Evaluate g
X0star <- b * Kstar^a
f1Star <- ROV(Kstar, X0star, r, tau1, m, s)
Kg <- Kstar * 0.10; tau2 <- 4 * 4; m <- r
stdev <- sqrt(2 * r)
f2Star <- f1Star * exp(-r * tau2)
g0 <- ROV(Kg, f2Star, r, tau2, m, stdev)
# Evaluate h
Kh <- Kstar * 0.05; tau3 <- 3 * 4; m <- r
stdev <- sqrt(2 * r)
h0 <- ROV(Kh, g0, r, tau3, m, stdev)
plot(K0, f0)
plot(Kg, g0)
plot(Kh, h0)

#





