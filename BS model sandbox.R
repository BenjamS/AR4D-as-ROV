library(tidyverse)
library(patchwork)

label_size <- 2.5
smallLabel_size <- 2
title_size <- 8
subtitle_size <- 7
legendText_size <- 6
axisText_size <- 6
axisTitle_size <- 7
facetTitle_size <- 7
#======================================================================
OVfun <- function(X0, K, tau, mm, s, r, output = "OV"){
  
  d2 <- (log(X0 / K) + (mm - s^2 / 2) * tau) / (s * sqrt(tau))
  d1 <- d2 + s * sqrt(tau)
  N1 <- pnorm(d1)
  N2 <- pnorm(d2)
  
  OV <- exp((mm - r) * tau) * X0 * N1 - exp(-r * tau) * K * N2
  
  if(output == "OV"){
    out <- OV
  }
  
  if(output == "N2"){
    out <- N2
  }
  
  if(output == "N1"){
    out <- N1
  }
  
  return(out)
}
#=====================================================================
NewOV <- function(X0, K = NULL, s, r){
  u <- 1 - 2 * r / s^2
  Kstar <- X0 * ((2 * r / s^2) / (1 - 2 * r / s^2))^(1 / (2 * r / s^2 - 1))
  if(is.null(K)){
    K <- Kstar
  }
  OV <- -K^(2 * r / s^2) / u * (X0^u - K^u)
  outVec <- c(OV, Kstar)
  return(outVec)
}
#=====================================================================
r <- 0.005
tau <- 15
cv <- seq(.7, 2, length.out = 3)
cvTxt <- c("low", "medium", "high")
s <- round(cv * r * sqrt(tau), 4)
#s <- round(cv * r, 4)
XoK <- seq(0.1, 3, length.out = 50)
X0 <- 1
K <- X0 / XoK
NPV <- X0 - exp(-r * tau) * K
this_s <- s[3]

u <- 1 - 2 * r / this_s^2
OV <- exp(-r * tau) * K^(1 - u) / (2 * u) * (X0^u - K^u)
OV2 <- OVfun(X0, K, tau, r, this_s, r, output = "OV")
#d <- (log(X0 / K) + (r - this_s^2 / 2) * tau) / (this_s * sqrt(tau))
#fStar <- exp(-r * tau) * K * (exp(-this_s^2 / 2 * tau) - pnorm(d))
#dif <- OV2 - fStar
aLpha <- 1.14
bEta <- 1
a0 <- (aLpha * K^bEta - K) / (2 *(aLpha * K^(bEta - 1) - 1 - bEta * log(K) + log(aLpha))) #K / (2 * (1 + log(n) / (n - 1)))
c0 <- 1 / 2 + a0 / K
g0 <- a0 * log(K) - c0 * K
fStar <- c0 * X0 - exp(-r * tau) * a0 * (log(X0) + (r - this_s^2 / 2) * tau) + exp(-r * tau) * g0
#fStar <- c0 * (X0 - exp(-r * tau) * K) + exp(-r * tau) * a0 * (log(K / X0) - (r - this_s^2 / 2) * tau)

dfPlot <- data.frame(OV2, fStar, NPV, XoK) %>% gather(Type, Val, OV2:NPV)
gg <- ggplot(dfPlot, aes(x = XoK, y = Val, group = Type, color = Type))
gg <- gg + geom_line()
gg <- gg + coord_cartesian(xlim = c(0, 2.5), ylim = c(-.1, 0.75))
gg <- gg + geom_hline(yintercept = 0, color = "red", linetype = "dashed")
gg <- gg + geom_vline(xintercept = 1, color = "red", linetype = "dashed")
gg
#



tau <- 1:35
d <- (log(X0 / K) + (r - this_s^2 / 2) * tau) / (this_s * sqrt(tau))
fStar <- exp(-r * tau) * K * (exp(-this_s^2 / 2 * tau) - pnorm(d))













outNew <- NewOV(X0, K = NULL, this_s, r)
#print(outNew)
newROV <- outNew[1]
Kstar <- outNew[2]


















ns <- length(cv)
list_df <- list()
for(i in 1:ns){
  this_s <- s[i]
  ROV <- OVfun(X0, K, tau, r, this_s, r, output = "OV")
  Phi2 <- OVfun(X0, K, tau, r, this_s, r, output = "N2")
  df <- data.frame(xx = X0 / K, Type = paste("ROV", cvTxt[i], "risk"), Phi2, Value = ROV)
  outNew <- NewOV(X0, K = NULL, this_s, r)
  #print(outNew)
  newROV <- outNew[1]
  Kstar <- outNew[2]
  dfNew <- data.frame(xx = X0 / K, Type = paste("New ROV", cvTxt[i], "risk"), Phi2 = NA, Value = newROV)
  list_df[[i]] <- as.data.frame(rbind(df, dfNew))
  
}
df_npv <- data.frame(xx = X0 / K, Type = "NPV", Phi2 = NA, Value = NPV)
list_df[[ ns + 1]] <- df_npv

df_plot <- as.data.frame(do.call(rbind, list_df))
colnames(df_plot)[1] <- "x(0)/K"
# df_plot1 <- subset(df_plot, cv == cv[ns])
# df_plot1$cv <- NULL
# df_plot1 <- df_plot1 %>% gather(type, Value, ROV:NPV)
# df_plot2 <- subset(df_plot, cv != cv[ns])

#n <- length(unique(df_plot1$type))
n <- ns + 1
bag_of_colors <- randomcoloR::distinctColorPalette(k = 2 * n)
color_vec <- sample(bag_of_colors, n)

# gg <- ggplot()
# gg <- gg + geom_line(data = df_plot2, aes(x = `x(0)/K`, y = ROV,
#                                          group = cv),
#                      color = color_vec[2], lwd = 1)
# gg <- gg + geom_line(data = df_plot1, aes(x = `x(0)/K`, y = Value,
#                                           group = type, color = type),
#                      lwd = 1)
#----------------------------------------------------------------------------
# ROV vs. NPV plot
gg <- ggplot(df_plot, aes(x = `x(0)/K`, y = Value,
                          group = Type, color = Type))
gg <- gg + geom_line(lwd = 1)
gg <- gg + scale_color_manual(values = color_vec)
#gg <- gg + labs(y = "Value")
gg <- gg + geom_hline(yintercept = 2, color = "red", linetype = "dashed")
gg <- gg + geom_hline(yintercept = 0)
gg <- gg + theme_bw()
gg <- gg + theme(axis.text = element_blank(),
                 axis.title.y = element_text(size = axisTitle_size),
                 axis.title.x = element_blank(),
                 legend.title = element_blank(),
                 legend.position = "top",
                 legend.text = element_text(size = legendText_size))
gg <- gg + guides(color = guide_legend(nrow = 2, byrow = T))
gg1 <- gg













aa = seq(0, 1, length.out = 20)
r = 0.05
t = 0
v <- r^(-aa - 1) * gamma(aa + 1)

plot(aa, v)



















Ta <- 29
A <- 5 * 7 * 2.5 * 2.54 * 10^-3
k <- 0.025
t <- seq(0, 240, length.out = 50)
T0 <- 9
Tt <- Ta + exp(A * k * t) * (T0 - Ta)
plot(t, Tt)




