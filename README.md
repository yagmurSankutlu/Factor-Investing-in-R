# 📈 Factor Investing in R

This project provides R scripts to see the performances of common investment factor strategies on the Borsa Istanbul (BIST100) using historical stock data.

### Methodology Overview

Each script follows a standard backtesting process. Raw factor data is calculated for each stock, cleaned by winsorizing outliers, and then standardized using a monthly cross-sectional Z-score. Stocks are then sorted into five quintile portfolios based on their factor scores. The performance of these portfolios is analyzed, and a long-short Factor-Mimicking Portfolio (FMP) is constructed to measure the overall factor premium.

---

## 🚀 Getting Started

### 1. Install Packages
Run this command once in your R console to install the necessary packages.

```R
install.packages(c("zoo", "xts", "quantstrat", "quantmod", "roll", "dplyr", "ggplot2", "TTR", "devtools"))


For the P/E and ROE factors, use Database.RData.
load("Database (1).RData")
source("PE_Ratio.R") 
# or
source("Roe_Factor.R")

For Momentum, Trend, Volatility, and the Composite factor analysis, use Database.RData.
load("Database.RData")
source("MomentumFactor.R") 
```

Running a script will generate performance plots and summary statistics in the R console.

P.S : This is a group project if you want to reach the codes of other factors, you can use this link: https://drive.google.com/drive/folders/1I2-mFxhYvuiWoXQ3TzZHxtmDSBoHCzG1?usp=sharing
