library(vrtest)
library(tseries)
library(pracma)
library(trend)
library(DescTools)
library(lmtest)
library(tidyverse)

data <- read.csv("R/data/btc_prices.csv")
data$ticker <- rep("BTC", nrow(data))
data$date <- as.Date(data$date)
data$return <- c(NA, diff(log(data$close)))
prices <- data$close
returns <- data$return[2:nrow(data)]

# Box-Pierce and Ljung-Box tests - Null is no autocorrelation
# Box.test(returns, lag = log(nrow(data)), type = "Box-Pierce", fitdf = 0)
Box.test(returns, lag = floor(log(nrow(data))), type = "Ljung-Box", fitdf = 0)

# Auto.Q test - automatic order of lags and robust to conditional heteroscedasticity
Auto.Q(returns)

# Dickey-Fuller - Null is unit root (non-stationarity)
adf.test(returns)

# KPSS - Null is stationarity
kpss.test(returns, null = "Trend", lshort = TRUE)

# Variance ratio manual - Null is random walk
variance_ratio_test(prices, 8)

# Variance Ratio test - Null is random walk
Auto.VR(returns)$stat
pnorm(Auto.VR(returns)$stat, lower.tail = FALSE)

# Runs test - Null is i.i.d.
runs.test(factor(na.omit(data$return) > mean(data$return, na.rm = TRUE)))

# Bartels test - Null is i.i.d.
# bartels.test(na.omit(data$return)) # Wrong
BartelsRankTest(na.omit(data$return), alternative = "two.sided", method = "normal")

# Hurst exponent
hurstexp(returns)


# Variance Ratio test - Null is random walk
variance_ratio_test <- function(prices, q) {
    log_prices <- log(prices)
    nqp1 <- length(prices)
    n <- floor((nqp1 - 1) / q)

    # Compute the variance estimators
    mu_hat <- 1 / (n * q) * sum(diff(log_prices[1:(n * q + 1)]))
    m <- q * (n * q - q + 1) * (1 - q / (n * q))

    sigma2_a <- 1 / (n * q - 1) * sum((log_prices[2:(n * q + 1)] - log_prices[1:(n * q)] - mu_hat)^2)
    sigma2_c <- 1 / (m) * sum((log_prices[(q + 1):(n * q + 1)] - log_prices[1:(n * q + 1 - q)] - mu_hat * q)^2)

    # Compute the variance ratio
    VR_q <- sigma2_c / sigma2_a

    delta_k <- sapply(1:(q - 1), function(k) {
        num <- sum((log_prices[(k + 2):(n * q + 1)] - log_prices[(k + 1):(n * q)] - mu_hat)^2 * (log_prices[2:(n * q + 1 - k)] - log_prices[1:(n * q - k)] - mu_hat)^2)
        denom <- sum((log_prices[2:(n * q + 1)] - log_prices[1:(n * q)] - mu_hat)^2)
        (n * q * num) / denom
    })

    # Compute asymptotic variance theta(q)
    theta_q <- 4 * sum(((1:(q - 1)) / q)^2 * delta_k)

    # Compute the standard test statistic
    psi_q <- sqrt(n * q) * (VR_q - 1) / sqrt(theta_q)

    # Return results
    list(
        VR_q = VR_q,
        theta_q = theta_q,
        test_statistic = psi_q,
        p_value = 2 * (1 - pnorm(abs(psi_q)))
    )
}

# Bartels RVN test for randomness
n <- length(returns)
ranks <- rank(returns) # Compute ranks of the observations
R_bar <- mean(ranks) # Mean of the ranks

numerator <- sum((ranks[-n] - ranks[-1])^2) # Sum of squared rank differences
denominator <- sum((ranks - R_bar)^2) # Sum of squared rank deviations

RVN <- numerator / denominator # Compute the test statistic

if (n <= 100) {
    p_q <- (5 * n * (n + 1) * (n - 1)^2) / (2 * (n - 2) * (5 * n^2 - 2 * n - 9))
    beta_p_value <- pbeta(RVN / 4, p_q, p_q) # p-value from beta distribution
} else {
    mean_rvn <- 2
    var_rvn <- 4 / n
    z_score <- (RVN - mean_rvn) / sqrt(var_rvn)
    beta_p_value <- pnorm(z_score, lower.tail = F) # p-value from normal distribution
}

list(RVN = RVN, p_value = beta_p_value)






data <- read.csv("R/data/btc_prices.csv")


data_before <- data[data$date >= as.Date("2024-10-06") & data$date <= as.Date("2024-11-05"), ]
data_during <- data[data$date >= as.Date("2024-11-06") & data$date <= as.Date("2024-12-05"), ]
data_after <- data[data$date >= as.Date("2024-12-06") & data$date <= as.Date("2025-01-05"), ]



returns <- diff(log(data_before$close))
adf.test(returns)
returns <- diff(log(data_during$close))
adf.test(returns)
returns <- diff(log(data_after$close))
adf.test(returns)

variance_ratio_test(data_after$close, 3)

variance_ratio_test(subset(data, date >= as.Date("2024-12-06") & date <= as.Date("2025-01-05"))$close, 3)
