install.packages("devtools")
require(devtools)
install_github("braverock/quantstrat")
install.packages("FinancialInstrument") 
install.packages("PerformanceAnalytics")
library(zoo)

# Plot fnc
TwoVarPlot<-function(plotdata,title="",leg1="",leg2=""){
  oldpar<-par()
  par(mar=c(3,4,3,4))
  plot(plotdata[,1],main=title,xlab="",ylab="",col="red",las=1)
  axis(2,col.axis="red",col="red",las=1)
  par(new=TRUE)
  plot(plotdata[,2],xaxt="n", yaxt="n",col="blue",ylab="",xlab="")
  axis(4,col.axis="blue",col="blue",las=1)
  legend("topleft", c(leg1,leg2), text.col=c("red","blue"),cex=0.8,bty="n")
  options(warn=-1)
  try(suppressWarnings(par(oldpar)),silent = T)
  options(warn=0)
}

# Load Data
load("Database.RData")

# Calculate Momentum Factor (12-month return)
# Momentum = Pt/P(t-12) - 1
Close.Data <- Close[,-c(1:4)] # Remove BIST indices

# Remove rows where more than 80% of stocks have NA
na_threshold <- 0.8 * ncol(Close.Data)
too_many_nas <- apply(is.na(Close.Data), 1, sum) > na_threshold
Close.Data.Clean <- Close.Data[!too_many_nas, ]

# Calculate 12-month momentum
Momentum.Factor <- Close.Data.Clean/lag(Close.Data.Clean, -252) - 1  
# % Ratio of the current price (P_t) to the price 12 months ago (P_{t-12})
colnames(Momentum.Factor) <- paste0("MOM_", colnames(Close.Data.Clean))

# Use monthly data - Get month-end dates from the Close data directly
Month.Ends <- xts::endpoints(Close, on="months")

# Get monthly factor values
Month.Ends.Dates <- time(Close)[Month.Ends[-1]]
# Filter momentum factor to monthly data
Factor.Data <- Momentum.Factor[Month.Ends.Dates]

# Remove any NA dates and ensure proper date range
Month.Ends.Dates <- Month.Ends.Dates[!is.na(Month.Ends.Dates)]


head(Factor.Data, 5)

# Plot data
#png("momentum_factor_distribution.png", width=800, height=600)
hist(c(coredata(Factor.Data)), 1000, main="Momentum Factor Distribution")
#dev.off()

# Winsorize Outliers
dat <- coredata(Factor.Data) # convert to matrix
dat[dat < -1] <- -1  # Cap at -100% (cannot lose more than 100%)
dat[dat > 5] <- 5    # Cap at 500% (reasonable upper bound for momentum)
Factor.Data <- zoo(dat, time(Factor.Data))

# Standardize
Z.Score <- function(x) (x-mean(x, na.rm=T))/sd(x, na.rm=T)

Factor.Data.Std <- t(apply(Factor.Data, 1, Z.Score)) # gives matrix
Factor.Data.Std <- zoo(Factor.Data.Std, time(Factor.Data))

# Plot standardized data
#png("standardized_momentum_distribution.png", width=800, height=600)
hist(c(coredata(Factor.Data.Std)), 1000, main="Standardized Momentum Factor Distribution")
#dev.off()

# Get prices and returns
Prices <- window(na.locf(Close[,-c(1:4)]), Month.Ends.Dates)
Rets <- Prices/lag(Prices, -1, na.pad=T) - 1

Bench <- window(na.locf(Close[,"XU100"]), Month.Ends.Dates)
Bench.Rets <- Bench/lag(Bench, -1, na.pad=T) - 1

# Align all datasets to the same time periods
# Use %in% to avoid numeric conversion
factor_dates <- time(Factor.Data.Std)
rets_dates <- time(Rets)
common_dates <- factor_dates[factor_dates %in% rets_dates]

Factor.Data.Std.Aligned <- Factor.Data.Std[common_dates]
Rets.Aligned <- Rets[common_dates]

# Recalculate IC with proper alignment
IC <- vapply(2:min(nrow(Rets.Aligned), nrow(Factor.Data.Std.Aligned) + 1),
             function(i) cor(as.numeric(Factor.Data.Std.Aligned[i-1,]),
                             as.numeric(Rets.Aligned[i,]),
                             use="pairwise.complete.obs"),
             3.54)
IC <- zoo(IC, time(Rets.Aligned)[2:(length(IC)+1)])

# Plot IC
#png("momentum_IC.png", width=800, height=600)
plot(IC, type="h", main="Momentum Factor Info Coefficient", las=1, xlab="", ylab="")
lines(rollmeanr(IC, 12), col=2)
abline(h=mean(IC,na.rm=TRUE), col=4)
title(sub=paste("Avg IC = ", round(100*mean(IC,na.rm = TRUE), 1), "%", sep=""))
#dev.off()

# Calculate Rank IC with aligned data
IC.Rank <- vapply(2:min(nrow(Rets.Aligned), nrow(Factor.Data.Std.Aligned) + 1),
                  function(i) cor(
                    rank(as.numeric(Factor.Data.Std.Aligned[i-1,]), na.last="keep"),
                    as.numeric(Rets.Aligned[i,]),
                    use="pairwise.complete.obs"),
                  3.54)
IC.Rank <- zoo(IC.Rank, time(Rets.Aligned)[2:(length(IC.Rank)+1)])

# Plot Rank IC
#png("momentum_rank_IC.png", width=800, height=600)
plot(IC.Rank, type="h", main="Momentum Factor Rank Info Coefficient",
     las=1, xlab="", ylab="")
lines(rollmeanr(IC.Rank, 12), col=2)
abline(h=mean(IC.Rank,na.rm=TRUE), col=4)
title(sub=paste("Avg Rank IC = ", round(100*mean(IC.Rank,na.rm = TRUE), 1),
                "%", sep=""))
#dev.off()

# Construct 5 portfolios based on momentum
suppressMessages(library(ggplot2))

Find.Portfolios <- function(x) {
  valid_data <- x[!is.na(x)]
  if(length(valid_data) < 5) {
    # If less than 5 valid observations, return NAs
    result <- rep(NA, length(x))
  } else {
    # Create portfolios only for valid data, then map back
    result <- rep(NA, length(x))
    result[!is.na(x)] <- as.numeric(cut_number(valid_data, n = 5))
  }
  return(result)
}

# Use aligned factor data for portfolio construction
Port.IDs <- t(apply(Factor.Data.Std.Aligned, 1, Find.Portfolios))
Port.IDs <- zoo(Port.IDs, time(Factor.Data.Std.Aligned))

head(Port.IDs, 5) # Stock i belongs to which portfolio

# Portfolio Construction Function
Construct.Port <- function(ID){
  Weights <- (Port.IDs == ID) # Weights
  Norm.Weights <- t(apply(Weights, 1, function(x) x/sum(x, na.rm=T)))
  Norm.Weights <- zoo(Norm.Weights, time(Weights)) # Normalized
  
  Norm.Weights <- lag(Norm.Weights, -1) # Lag weights
  
  # Use aligned returns and handle dimension matching
  # Ensure dimensions match by using only overlapping periods
  overlap_dates <- time(Norm.Weights)[time(Norm.Weights) %in% time(Rets.Aligned)]
  Norm.Weights.Overlap <- Norm.Weights[overlap_dates]
  Rets.Overlap <- Rets.Aligned[overlap_dates]
  
  # Calculate portfolio returns
  Port.Rets <- apply(Norm.Weights.Overlap * Rets.Overlap, 1, sum, na.rm=T)
  Port.Rets <- zoo(Port.Rets, overlap_dates)        
  Port <- cumprod(1 + Port.Rets)
  Port
}

# Construct quintile portfolios (P1=Low Momentum, P5=High Momentum)
P1 <- Construct.Port(1)  # Losers
P2 <- Construct.Port(2)
P3 <- Construct.Port(3)
P4 <- Construct.Port(4)
P5 <- Construct.Port(5)  # Winners

All.Portfolios <- merge(P1, P2, P3, P4, P5)

#png("momentum_portfolios.png", width=1000, height=600)
plot(All.Portfolios, screens=1, col=1:5,
     xlab="", ylab="", las=1, main="Momentum Factor Portfolios")
legend("topleft", paste("Port ", 1:5, " (", c("Losers", "Low", "Mid", "High", "Winners"), ")", sep=""), 
       text.col=1:5, bty="n")
#dev.off()

#png("momentum_portfolios_log.png", width=1000, height=600)
plot(All.Portfolios, screens=1, col=1:5,
     xlab="", ylab="", las=1, log="y", main="Momentum Factor Portfolios (Log Scale)")
legend("topleft", paste("Port ", 1:5, " (", c("Losers", "Low", "Mid", "High", "Winners"), ")", sep=""), 
       text.col=1:5, bty="n")
#dev.off()

# Fnc for CAGR return
CAGR <- function(myzoo){
  (tail(myzoo, 1)/as.numeric(head(myzoo, 1)))^
    (365/as.numeric(diff(range(time(myzoo))))) - 1
}

All.CAGR <- cbind(CAGR(All.Portfolios[,1]),
                  CAGR(All.Portfolios[,2]),
                  CAGR(All.Portfolios[,3]),
                  CAGR(All.Portfolios[,4]),
                  CAGR(All.Portfolios[,5]))

#png("CAGR_by_portfolio.png", width=800, height=600)
barplot(as.numeric(All.CAGR), names.arg=paste("Port ", 1:5, sep=""),
        main="CAGR by Momentum Portfolio", las=1)
#dev.off()

All.Vol <- apply(All.Portfolios, 2,
                 function(x) sqrt(12)*sd(diff(log(x))))

#png("volatility_by_portfolio.png", width=800, height=600)
barplot(All.Vol, names.arg=paste("Port ", 1:5, sep=""),
        main="Annualized Volatility by Momentum Portfolio", las=1)
#dev.off()

# Fnc for max drawdown
MaxDD <- function(myzoo){
  min(myzoo/TTR::runMax(myzoo, n=1, cumulative=T) - 1, na.rm=T)  
}

All.MDD <- apply(All.Portfolios, 2, MaxDD)
#png("max_drawdown_by_portfolio.png", width=800, height=600)
barplot(All.MDD, names.arg=paste("Port ", 1:5, sep=""),
        main="Max Drawdowns by Momentum Portfolio", las=1)
#dev.off()

# Construct Momentum Factor Mimicking Portfolio (FMP)
# Winners minus Losers (P5 - P1)
FMP.Rets<-(P5/lag(P5,-1)-1)-(P1/lag(P1,-1)-1)
FMP <- cumprod(1 + FMP.Rets)

#png("FMP_returns.png", width=800, height=600)
barplot(FMP.Rets, main="Momentum FMP Returns (Winners - Losers)")
#dev.off()

#png("FMP_cumulative.png", width=800, height=600)
plot(FMP, main="Momentum Factor Mimicking Portfolio (Winners - Losers)")
#dev.off()

# t-test for momentum factor
print("T-test for Momentum Factor (Portfolio Approach):")
print(t.test(FMP.Rets, alternative="greater"))

# Cross-sectional regression functions
Get.Slope <- function(i){
  dat <- cbind(as.numeric(Rets.Aligned[i,]),
               as.numeric(Factor.Data.Std.Aligned[i-1,]))
  dat <- na.omit(dat)
  
  if(nrow(dat) < 5) {  # Need at least 5 observations for regression
    return(NA)
  }
  
  fit <- lm(dat[,1] ~ dat[,2])
  summary(fit)$coefficients[2,1]
}

# Obtain factor returns through cross sectional regressions
FMP.CS.Rets <- vapply(2:min(nrow(Rets.Aligned), nrow(Factor.Data.Std.Aligned) + 1), Get.Slope, 3.54)
FMP.CS.Rets <- zoo(FMP.CS.Rets, time(Rets.Aligned)[2:(length(FMP.CS.Rets)+1)])

# FMP Cross-Sectional
FMP.CS <- cumprod(1 + FMP.CS.Rets)

#png("FMP_CS_returns.png", width=800, height=600)
barplot(FMP.CS.Rets, main="Momentum FMP Returns (Cross-Sectional)")
#dev.off()

#png("FMP_CS_cumulative.png", width=800, height=600)
plot(FMP.CS, main="Momentum FMP (Cross-Sectional Regression)")
#dev.off()

# t-test for cross-sectional momentum factor
print("T-test for Momentum Factor (Cross-Sectional Approach):")
print(t.test(FMP.CS.Rets, alternative="greater"))

# FMP weights function
Get.FMP.Weights <- function(i){
  dat <- cbind(as.numeric(Rets.Aligned[i,]),
               as.numeric(Factor.Data.Std.Aligned[i-1,]))
  non.nas <- which(apply(dat, 1, function(x) (!is.na(x[1])) & (!is.na(x[2]))))
  dat <- na.omit(dat)
  
  if(nrow(dat) < 5) {
    return(rep(0, ncol(Rets.Aligned)))
  }
  
  XX <- cbind(rep(1, dim(dat)[1]), dat[,2])
  WW <- t(solve(t(XX) %*% XX) %*% t(XX))[,2]
  
  W <- rep(0, ncol(Rets.Aligned))
  W[non.nas] <- WW
  W
}

# Obtain factor returns through cross sectional regressions
FMP.W <- t(vapply(2:min(nrow(Rets.Aligned), nrow(Factor.Data.Std.Aligned) + 1), 
                  Get.FMP.Weights, rep(3.54, ncol(Rets.Aligned))))
FMP.W <- zoo(FMP.W, time(Rets.Aligned)[2:(nrow(FMP.W)+1)])

#png("FMP_weights_last_period.png", width=1200, height=600)
barplot(as.numeric(tail(FMP.W, 1)),
        main="Momentum FMP Weights for Last Period", las=2,
        names.arg=colnames(Rets.Aligned), cex.names=0.4)
#dev.off()

# Align FMP series for comparison
common_fmp_dates <- time(FMP.Rets)[time(FMP.Rets) %in% time(FMP.CS.Rets)]
FMP.Rets.Aligned <- FMP.Rets[common_fmp_dates]
FMP.CS.Rets.Aligned <- FMP.CS.Rets[common_fmp_dates]

#png("FMP_comparison.png", width=1000, height=600)
TwoVarPlot(merge(FMP.Rets.Aligned, FMP.CS.Rets.Aligned), "Momentum FMPs",
           "FMP (Portfolio Approach)",
           "FMP (CS Regression Approach)")
#dev.off()

# Volatility and correlation of FMPs
cat("Volatility of Portfolio FMP:", sqrt(12)*sd(FMP.Rets, na.rm=TRUE), "\n")
cat("Volatility of CS Regression FMP:", sqrt(12)*sd(FMP.CS.Rets, na.rm=TRUE), "\n")
cat("Correlation between FMPs:", cor(FMP.Rets.Aligned, FMP.CS.Rets.Aligned, use="complete.obs"), "\n")

# Define TRY.Rets outside the conditional (initialize as NULL)
TRY.Rets <- NULL

# Use USDTRY as macro factor (if available)
if("Turkey" %in% colnames(FX)) {
  TRY <- window(FX[,"Turkey"], Month.Ends.Dates)
  TRY.Rets <<- TRY/lag(TRY, -1, na.pad=T) - 1
  
  # Use P1-P5 as base assets
  Base.Asset.Rets <- All.Portfolios
  
  # Fnc to compute weights
  Get.TS.Weights <- function(i){
    dat <- na.omit(merge(TRY.Rets, Base.Asset.Rets))
    if(i <= 35 || nrow(dat) < (i-35)) return(rep(NA, 5))
    
    dat <- window(dat, time(TRY.Rets)[(i-35):i])
    
    if(nrow(dat) < 10) return(rep(NA, 5))
    
    fit <- lm(dat[,1] ~ dat[,-1])
    W <- summary(fit)$coefficients[2:6,1]
    W/sum(W)
  }
  
  valid_indices <- 37:min(nrow(Rets.Aligned), length(time(TRY.Rets)))
  FMP.TS.W <- t(vapply(valid_indices, Get.TS.Weights, rep(3.54, 5)))
  FMP.TS.W <- zoo(FMP.TS.W, time(Rets.Aligned)[valid_indices])
  
  #png("FMP_TS_weights.png", width=800, height=600)
  barplot(as.numeric(tail(FMP.TS.W, 1)), names.arg=paste("Port ", 1:5, sep=""),
          las=2, main="Momentum FMP Weights (Time Series)")
  #dev.off()
}

# Now you can use TRY.Rets outside the conditional
if(!is.null(TRY.Rets)) {
  # 1. Check correlation between portfolios and macro factor
  macro_correlations <- cor(merge(All.Portfolios, TRY.Rets), use="complete.obs")
  print("Correlations with TRY:")
  print(macro_correlations[6, 1:5])
} else {
  print("Turkey FX data not available - cannot compute correlations")
}

# Principal component analysis
PCA <- princomp(coredata(diff(log(All.Portfolios))))
summary(PCA)

#png("PCA_scree_plot.png", width=800, height=600)
plot(PCA, main="Eigenvalue Scree Plot")
#dev.off()

loadings(PCA)

PCA.FMP <- loadings(PCA)[,1]
PCA.FMP <- PCA.FMP/sum(PCA.FMP)

#png("PCA_FMP_weights.png", width=800, height=600)
barplot(PCA.FMP, main="PCA Momentum FMP Weights",
        names.arg=paste("Port ", 1:5, sep=""), las=2)
#dev.off()

PCA.Factor.Rets <- zoo(PCA$scores, time(diff(log(All.Portfolios))))
head(PCA.Factor.Rets)

dat <- na.omit(merge(PCA.Factor.Rets, Bench.Rets))
if(ncol(dat) >= 6) {
  correlation <- cor(dat)[6,1]
  cat("Correlation between First PCA Factor and Market:", correlation, "\n")
  
  #png("PCA_vs_market.png", width=1000, height=600)
  TwoVarPlot(log(merge(cumprod(1+dat[,1]), cumprod(1+dat[,6]))),
             "PCA Momentum Factor Model", "log(First PCA Factor)", "log(Market Factor)")
  #dev.off()
}

# Summary Statistics
cat("\n=== MOMENTUM FACTOR ANALYSIS SUMMARY ===\n")
cat("Average IC:", round(100*mean(IC, na.rm=TRUE), 2), "%\n")
cat("Average Rank IC:", round(100*mean(IC.Rank, na.rm=TRUE), 2), "%\n")
cat("Portfolio FMP Annual Return:", round(100*mean(FMP.Rets, na.rm=TRUE)*12, 2), "%\n")
cat("Cross-Sectional FMP Annual Return:", round(100*mean(FMP.CS.Rets, na.rm=TRUE)*12, 2), "%\n")
cat("FMP Correlation:", round(cor(FMP.Rets.Aligned, FMP.CS.Rets.Aligned, use="complete.obs"), 3), "\n")

# Portfolio Performance Summary
cat("\n=== PORTFOLIO PERFORMANCE ===\n")

# Fixed portfolio statistics printing
portfolio_names <- c("Losers", "Low", "Mid", "High", "Winners")

for(i in 1:5) {
  cat("Portfolio", i, "(", portfolio_names[i], "):\n")
  
  # For All.CAGR: it's a zoo object with matrix structure, access by column
  cat("  CAGR:", round(100*All.CAGR[1,i], 2), "%\n")
  
  # For All.Vol and All.MDD: they're vectors, access by index
  cat("  Volatility:", round(100*All.Vol[i], 2), "%\n")
  cat("  Max Drawdown:", round(100*All.MDD[i], 2), "%\n\n")
}

# Alternative approach using as.numeric() to convert zoo to vector
cat("\n=== Alternative approach ===\n")
CAGR_values <- as.numeric(All.CAGR[1,])  # Extract the row as a vector

for(i in 1:5) {
  cat("Portfolio", i, "(", portfolio_names[i], "):\n")
  cat("  CAGR:", round(100*CAGR_values[i], 2), "%\n")
  cat("  Volatility:", round(100*All.Vol[i], 2), "%\n")
  cat("  Max Drawdown:", round(100*All.MDD[i], 2), "%\n\n")
}

cat("All plots saved as PNG files in your working directory!\n")

# Additional analysis to understand CAGR vs FMP weight divergence

# 1. Check correlation between portfolios and macro factor
macro_correlations <- cor(merge(All.Portfolios, TRY.Rets), use="complete.obs")
print("Correlations with TRY:")
print(macro_correlations[6, 1:5])






# Fix rolling correlations with proper data alignment
library(zoo)

# Debug: Check the date ranges and alignment issues
print("P1 date range:")
print(paste("Start:", start(P1), "End:", end(P1)))
print("TRY.Rets date range:")
print(paste("Start:", start(TRY.Rets), "End:", end(TRY.Rets)))

# Method 1: Use proper merge with all=FALSE to avoid overlapping index issues
if(exists("P1") && exists("TRY.Rets") && !is.null(TRY.Rets)) {
  
  # Clean merge - only keep dates that exist in both series
  merged_data_P1 <- merge(P1, TRY.Rets, all=FALSE)
  print("Clean merged P1-TRY data:")
  print(paste("Observations:", nrow(merged_data_P1)))
  print("Date range:")
  print(paste("Start:", start(merged_data_P1), "End:", end(merged_data_P1)))
  print("First few rows:")
  print(head(merged_data_P1))
  print("Structure:")
  print(str(merged_data_P1))
  
  # Remove any rows with NA values
  merged_data_P1_clean <- na.omit(merged_data_P1)
  print(paste("After removing NAs:", nrow(merged_data_P1_clean), "observations"))
  
  if(nrow(merged_data_P1_clean) >= 36) {
    # Calculate rolling correlation using the cleaned data
    rolling_cor_P1 <- rollapply(merged_data_P1_clean, 
                                width=36, 
                                FUN=function(x) {
                                  # x should be a matrix with 2 columns
                                  if(is.matrix(x) && ncol(x) == 2 && nrow(x) >= 10) {
                                    cor(x[,1], x[,2], use="complete.obs")
                                  } else {
                                    NA
                                  }
                                }, 
                                by.column=FALSE, 
                                fill=NA, 
                                align="right")
    
    print("P1 rolling correlation calculated successfully")
    print(paste("Rolling correlation observations:", length(rolling_cor_P1)))
  } else {
    print("Not enough overlapping data for P1 rolling correlation")
    rolling_cor_P1 <- NULL
  }
  
  # Same process for P5
  if(exists("P5")) {
    merged_data_P5 <- merge(P5, TRY.Rets, all=FALSE)
    merged_data_P5_clean <- na.omit(merged_data_P5)
    print(paste("P5 clean observations:", nrow(merged_data_P5_clean)))
    
    if(nrow(merged_data_P5_clean) >= 36) {
      rolling_cor_P5 <- rollapply(merged_data_P5_clean, 
                                  width=36, 
                                  FUN=function(x) {
                                    if(is.matrix(x) && ncol(x) == 2 && nrow(x) >= 10) {
                                      cor(x[,1], x[,2], use="complete.obs")
                                    } else {
                                      NA
                                    }
                                  }, 
                                  by.column=FALSE, 
                                  fill=NA, 
                                  align="right")
      
      print("P5 rolling correlation calculated successfully")
    } else {
      print("Not enough overlapping data for P5 rolling correlation")
      rolling_cor_P5 <- NULL
    }
  }
  
  # Plot rolling correlations if both exist
  if(!is.null(rolling_cor_P1) && !is.null(rolling_cor_P5)) {
    #png("rolling_correlations.png", width=1000, height=600)
    
    # Get the range for y-axis
    y_range <- range(c(as.numeric(rolling_cor_P1), as.numeric(rolling_cor_P5)), na.rm=TRUE)
    
    # Plot P1
    plot(rolling_cor_P1, 
         main="36-Month Rolling Correlation with TRY Returns", 
         ylim=y_range,
         col="blue", 
         lwd=2, 
         ylab="Correlation",
         xlab="Date")
    
    # Add P5 line
    lines(rolling_cor_P5, col="red", lwd=2)
    
    # Add legend
    legend("topright", 
           legend=c("P1 (Losers)", "P5 (Winners)"), 
           col=c("blue", "red"), 
           lwd=2)
    
    # Add horizontal line at zero
    abline(h=0, lty=2, col="gray")
    
    #dev.off()
    
    print("Rolling correlation plot saved successfully")
    
    # Print some summary statistics
    cat("\nRolling Correlation Summary:\n")
    cat("P1 (Losers) - Mean:", round(mean(rolling_cor_P1, na.rm=TRUE), 3), 
        "Min:", round(min(rolling_cor_P1, na.rm=TRUE), 3),
        "Max:", round(max(rolling_cor_P1, na.rm=TRUE), 3), "\n")
    cat("P5 (Winners) - Mean:", round(mean(rolling_cor_P5, na.rm=TRUE), 3), 
        "Min:", round(min(rolling_cor_P5, na.rm=TRUE), 3),
        "Max:", round(max(rolling_cor_P5, na.rm=TRUE), 3), "\n")
    
  } else {
    print("Could not create rolling correlations - insufficient data")
  }
  
} else {
  print("Required data not available for rolling correlation analysis")
}

# Alternative simpler approach for basic correlations
print("\n=== Static Correlations ===")
if(exists("All.Portfolios") && exists("TRY.Rets")) {
  static_correlations <- cor(merge(All.Portfolios, TRY.Rets, all=FALSE), use="complete.obs")
  print("Static correlations with TRY:")
  print(round(static_correlations[6, 1:5], 3))
}

# Plot rolling correlations
#png("rolling_correlations.png", width=1000, height=600)
plot(rolling_cor_P1, main="36-Month Rolling Correlation with TRY", 
     ylim=range(c(rolling_cor_P1, rolling_cor_P5), na.rm=TRUE),
     col="blue", lwd=2)
lines(rolling_cor_P5, col="red", lwd=2)
legend("topright", c("P1 (Losers)", "P5 (Winners)"), col=c("blue", "red"), lwd=2)
#dev.off()

# 3. Check beta relationships
beta_analysis <- function(portfolio_returns, macro_returns) {
  merged_data <- na.omit(merge(portfolio_returns, macro_returns))
  if(nrow(merged_data) < 10) return(NA)
  
  model <- lm(merged_data[,1] ~ merged_data[,2])
  return(list(beta = coef(model)[2], 
              r_squared = summary(model)$r.squared,
              t_stat = summary(model)$coefficients[2,3]))
}

# Calculate betas for each portfolio
portfolio_betas <- sapply(1:5, function(i) {
  port_rets <- diff(log(All.Portfolios[,i]))
  try_rets <- TRY.Rets[time(TRY.Rets) %in% time(port_rets)]
  beta_analysis(port_rets, try_rets)
})

print("Beta Analysis Results:")
print(portfolio_betas)

# 4. Examine factor loadings over time
factor_loadings <- function(window_size = 36) {
  loadings_matrix <- matrix(NA, nrow = length(time(All.Portfolios)) - window_size + 1, ncol = 5)
  
  for(i in window_size:length(time(All.Portfolios))) {
    window_data <- All.Portfolios[(i-window_size+1):i,]
    window_try <- TRY.Rets[time(TRY.Rets) %in% time(window_data)]
    
    if(length(window_try) > 10) {
      try_rets_aligned <- window_try[time(window_try) %in% time(diff(log(window_data)))]
      port_rets_aligned <- diff(log(window_data))[time(diff(log(window_data))) %in% time(try_rets_aligned)]
      
      if(nrow(port_rets_aligned) > 5) {
        model <- lm(try_rets_aligned ~ port_rets_aligned)
        loadings_matrix[i-window_size+1, ] <- coef(model)[-1]
      }
    }
  }
  
  loadings_zoo <- zoo(loadings_matrix, time(All.Portfolios)[window_size:length(time(All.Portfolios))])
  return(loadings_zoo)
}

rolling_loadings <- factor_loadings()

#png("rolling_factor_loadings.png", width=1200, height=800)
plot(rolling_loadings, screens=1, col=1:5, main="Rolling Factor Loadings (36-month window)")
legend("topright", paste("Port", 1:5), col=1:5, lty=1)
#dev.off()

# 5. Compare FMP performance vs individual portfolios
fmp_performance <- function(weights, portfolios) {
  # Create FMP using given weights
  weighted_returns <- sweep(diff(log(portfolios)), 2, weights, "*")
  fmp_returns <- apply(weighted_returns, 1, sum, na.rm=TRUE)
  return(zoo(fmp_returns, time(diff(log(portfolios)))))
}

# Equal-weighted FMP (traditional approach)
equal_weights <- c(-1, 0, 0, 0, 1)  # Short losers, long winners
equal_fmp <- fmp_performance(equal_weights, All.Portfolios)

# Time-series optimized FMP (your approach)
ts_weights <- as.numeric(tail(FMP.TS.W, 1))  # Use latest weights
ts_fmp <- fmp_performance(ts_weights, All.Portfolios)

# Compare performance
#png("fmp_comparison_strategies.png", width=1000, height=1200)
plot(cumprod(1 + ts_fmp), main="FMP Strategy Comparison", ylab="Cumulative Return", xlab="Date", col="blue", lwd=2)
lines(cumprod(1 + equal_fmp),  
      col="red", lwd=2)
legend("topleft", c("Portfolio Approach (P5-P1)", "Time-Series Optimized"), 
       col=c("red", "blue"), lwd=2, bty="n")
#dev.off()

# Performance metrics comparison
cat("Equal Weight FMP - CAGR:", 
    (tail(cumprod(1 + equal_fmp), 1)^(252/length(equal_fmp)) - 1) * 100, "%\n")
cat("Equal Weight FMP - Volatility:", sqrt(252) * sd(equal_fmp, na.rm=TRUE) * 100, "%\n")

cat("TS Optimized FMP - CAGR:", 
    (tail(cumprod(1 + ts_fmp), 1)^(252/length(ts_fmp)) - 1) * 100, "%\n")
cat("TS Optimized FMP - Volatility:", sqrt(252) * sd(ts_fmp, na.rm=TRUE) * 100, "%\n")


# Her iki serinin kümülatif getirisini 1 ile başlat
FMP_cum <- cumprod(1 + FMP.Rets)
FMP_CS_cum <- cumprod(1 + FMP.CS.Rets)

# Her iki seriyi aynı grafikte, aynı eksende çiz
plot(FMP_cum, type="l", col="red", lwd=2, 
     main="FMPs", ylab="Cumulative Return", xlab="Date")
lines(FMP_CS_cum, col="blue", lwd=2)
legend("topleft", legend=c("FMP (Portfolio Approach (P5-P1))", "FMP (CS Regression Approach)"),
       col=c("red","blue"), lwd=2, bty="n")


# === Summary Statistics for P1–P5 and Long/Short Strategy (P5 - P1) ============

# Get monthly returns for P1–P5
Port.Rets <- lapply(1:5, function(i) All.Portfolios[,i] / lag(All.Portfolios[,i], -1) - 1)
names(Port.Rets) <- paste0("P", 1:5)

# Calculate stats
AvgRet <- sapply(Port.Rets, function(r) mean(r, na.rm = TRUE))
AnnRet <- sapply(1:5, function(i) CAGR(All.Portfolios[,i]))
StdDev <- sapply(Port.Rets, function(r) sqrt(12) * sd(r, na.rm = TRUE))
# Align returns and benchmark by index
PctOutPerf <- sapply(Port.Rets, function(r) {
  # Merge to keep only overlapping dates
  aligned <- merge(r, Bench.Rets, join = "inner")
  mean(aligned[,1] > aligned[,2], na.rm = TRUE) * 100
})

# Combine into data.frame
PortStats <- data.frame(
  Portfolio = paste0("P", 1:5),
  AvgRet = round(AvgRet * 100, 2),
  AnnRet = round(AnnRet * 100, 2),
  StdDev = round(StdDev * 100, 2),
  PctOutPerf = round(PctOutPerf, 1)
)

cat("\n===== Portfolio Statistics (P1 to P5) =====\n")
print(PortStats)

# === Long/Short Strategy (P5 - P1) =============================================

LongShort.Rets <- Port.Rets[[5]] - Port.Rets[[1]]
LongShort.CAGR <- CAGR(cumprod(1 + LongShort.Rets))
LongShort.StdDev <- sqrt(12) * sd(LongShort.Rets, na.rm = TRUE)
LongShort.TStat <- t.test(LongShort.Rets, alternative = "greater")$statistic
LongShort.PctOut <- mean(LongShort.Rets > 0, na.rm = TRUE) * 100

cat("\n===== Long/Short Strategy (P5 - P1) =====\n")
cat(sprintf("Avg Ret     : %.2f%%\n", mean(LongShort.Rets, na.rm = TRUE) * 100))
cat(sprintf("Ann Ret     : %.2f%%\n", LongShort.CAGR * 100))
cat(sprintf("Std Dev     : %.2f%%\n", LongShort.StdDev * 100))
cat(sprintf("T-Statistic : %.2f\n", LongShort.TStat))
cat(sprintf("%% OutPerf   : %.1f%%\n", LongShort.PctOut))


save(image)

