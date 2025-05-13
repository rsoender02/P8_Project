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

variance_ratio_test(prices, floor(log(length(prices))))

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







data_before <- data[data$date >= as.Date("2024-10-06") & data$date <= as.Date("2024-11-05"), ]
data_during <- data[data$date >= as.Date("2024-11-06") & data$date <= as.Date("2024-12-05"), ]
data_after <- data[data$date >= as.Date("2024-12-06") & data$date <= as.Date("2025-01-05"), ]









# Loop through each year
for (i in 2015:2024) {
    year_str <- as.character(i)

    # Subset data for the year
    data_list[[year_str]] <- subset(data, date >= as.Date(paste0(i, "-01-01")) & date <= as.Date(paste0(i, "-12-31")))

    if (nrow(data_list[[year_str]]) > 1) {
        data_list[[year_str]]$return <- c(NA, diff(log(data_list[[year_str]]$close)))
        returns <- na.omit(data_list[[year_str]]$return)

        # Perform Tests
        test_results[[year_str]] <- list(
            Ljung_Box = Box.test(returns, lag = min(floor(log(length(returns))), length(returns) - 1), type = "Ljung-Box", fitdf = 0),
            Auto_Q = tryCatch(Auto.Q(returns), error = function(e) NA),
            ADF = adf.test(returns),
            KPSS = kpss.test(returns, null = "Trend", lshort = TRUE),
            Variance_Ratio = variance_ratio_test(data_list[[year_str]]$close, 5),
            Runs_Test = runs.test(factor(returns > mean(returns))),
            Bartels_Test = BartelsRankTest(returns, alternative = "two.sided", method = "normal"),
            Hurst_Exponent = tryCatch(hurstexp(returns), error = function(e) NA),
            Auto_VR_Test = tryCatch(
                {
                    vr_stat <- Auto.VR(returns)$stat
                    list(stat = vr_stat, p_value = pnorm(vr_stat, lower.tail = FALSE))
                },
                error = function(e) list(stat = NA, p_value = NA)
            )
        )
    } else {
        test_results[[year_str]] <- NA # If not enough data, store NA
    }
}


# Initialize an empty list to store p-values
p_values_list <- list()

# Loop through each year and extract p-values
for (i in 2015:2024) {
    year_str <- as.character(i)

    if (!is.null(test_results[[year_str]]) && is.list(test_results[[year_str]])) {
        p_values_list[[year_str]] <- c(
            Ljung_Box = test_results[[year_str]]$Ljung_Box$p.value,
            Auto_Q = test_results[[year_str]]$Auto_Q$Pvalue,
            ADF = test_results[[year_str]]$ADF$p.value,
            KPSS = test_results[[year_str]]$KPSS$p.value,
            Variance_Ratio = test_results[[year_str]]$Variance_Ratio$p_value,
            Runs_Test = test_results[[year_str]]$Runs_Test$p.value,
            Bartels_Test = test_results[[year_str]]$Bartels_Test$p.value,
            Hurst_Exponent = test_results[[year_str]]$Hurst_Exponent$Hs,
            Auto_VR_Test = test_results[[year_str]]$Auto_VR_Test$p_value
        )
    } else {
        # If no data for the year, store NA
        p_values_list[[year_str]] <- rep(NA, 9)
    }
}

# Convert the list to a dataframe
p_values_df <- as.data.frame(do.call(rbind, p_values_list))

# Set row names as years
rownames(p_values_df) <- names(p_values_list)

# Print the final dataframe
print(p_values_df["2020", ])

write.csv(p_values_df, file = "bitcoin_p_values_yearly.csv", row.names = TRUE)

p_values_df




# Initialize an empty list to store test statistics
test_statistic_list <- list()

# Loop through each year and extract test statistics
for (i in 2015:2024) {
    year_str <- as.character(i)

    if (!is.null(test_results[[year_str]]) && is.list(test_results[[year_str]])) {
        test_statistic_list[[year_str]] <- c(
            Ljung_Box = test_results[[year_str]]$Ljung_Box$statistic,
            Auto_Q = test_results[[year_str]]$Auto_Q$Stat,
            ADF = test_results[[year_str]]$ADF$statistic,
            KPSS = test_results[[year_str]]$KPSS$statistic,
            Variance_Ratio = test_results[[year_str]]$Variance_Ratio$test_statistic,
            Runs_Test = test_results[[year_str]]$Runs_Test$statistic,
            Bartels_Test = test_results[[year_str]]$Bartels_Test$statistic,
            Hurst_Exponent = test_results[[year_str]]$Hurst_Exponent$Hs,
            Auto_VR_Test = test_results[[year_str]]$Auto_VR_Test$stat
        )
    } else {
        # If no data for the year, store NA
        test_statistic_list[[year_str]] <- rep(NA, 9)
    }
}

# Convert the list to a dataframe
test_statistic_df <- as.data.frame(do.call(rbind, test_statistic_list))

# Set row names as years
rownames(test_statistic_df) <- names(test_statistic_list)

# Print the final dataframe
print(test_statistic_df["2020", ])

write.csv(test_statistic_df, file = "bitcoin_test_statistics_yearly.csv", row.names = TRUE)


test_statistic_df



# Trump election

# Define custom date ranges
date_ranges <- list(
    "2024-11-06_to_2024-12-05" = c(as.Date("2024-10-06"), as.Date("2024-11-05")),
    "2024-12-06_to_2025-01-05" = c(as.Date("2024-11-06"), as.Date("2024-12-05")),
    "2025-01-06_to_2025-02-05" = c(as.Date("2024-12-06"), as.Date("2025-01-05"))
)


# Initialize result lists
test_results <- list()
p_values_list <- list()
test_statistic_list <- list()

# Loop over the three custom date ranges
for (name in names(date_ranges)) {
    start_date <- date_ranges[[name]][1]
    end_date <- date_ranges[[name]][2]

    subset_data <- subset(data, date >= start_date & date <= end_date)

    if (nrow(subset_data) > 1) {
        subset_data$return <- c(NA, diff(log(subset_data$close)))
        returns <- na.omit(subset_data$return)

        test_results[[name]] <- list(
            Ljung_Box = Box.test(returns, lag = min(floor(log(length(returns))), length(returns) - 1), type = "Ljung-Box", fitdf = 0),
            Auto_Q = tryCatch(Auto.Q(returns), error = function(e) NA),
            ADF = adf.test(returns),
            KPSS = kpss.test(returns, null = "Trend", lshort = TRUE),
            Variance_Ratio = variance_ratio_test(subset_data$close, 3),
            Runs_Test = runs.test(factor(returns > mean(returns))),
            Bartels_Test = BartelsRankTest(returns, alternative = "two.sided", method = "normal"),
            Hurst_Exponent = tryCatch(hurstexp(returns), error = function(e) NA),
            Auto_VR_Test = tryCatch(
                {
                    vr_stat <- Auto.VR(returns)$stat
                    list(stat = vr_stat, p_value = pnorm(vr_stat, lower.tail = FALSE))
                },
                error = function(e) list(stat = NA, p_value = NA)
            )
        )

        # Extract p-values
        p_values_list[[name]] <- c(
            Ljung_Box = test_results[[name]]$Ljung_Box$p.value,
            Auto_Q = test_results[[name]]$Auto_Q$Pvalue,
            ADF = test_results[[name]]$ADF$p.value,
            KPSS = test_results[[name]]$KPSS$p.value,
            Variance_Ratio = test_results[[name]]$Variance_Ratio$p_value,
            Runs_Test = test_results[[name]]$Runs_Test$p.value,
            Bartels_Test = test_results[[name]]$Bartels_Test$p.value,
            Hurst_Exponent = test_results[[name]]$Hurst_Exponent$Hs,
            Auto_VR_Test = test_results[[name]]$Auto_VR_Test$p_value
        )

        # Extract test statistics
        test_statistic_list[[name]] <- c(
            Ljung_Box = test_results[[name]]$Ljung_Box$statistic,
            Auto_Q = test_results[[name]]$Auto_Q$Stat,
            ADF = test_results[[name]]$ADF$statistic,
            KPSS = test_results[[name]]$KPSS$statistic,
            Variance_Ratio = test_results[[name]]$Variance_Ratio$test_statistic,
            Runs_Test = test_results[[name]]$Runs_Test$statistic,
            Bartels_Test = test_results[[name]]$Bartels_Test$statistic,
            Hurst_Exponent = test_results[[name]]$Hurst_Exponent$Hs,
            Auto_VR_Test = test_results[[name]]$Auto_VR_Test$stat
        )
    } else {
        p_values_list[[name]] <- rep(NA, 9)
        test_statistic_list[[name]] <- rep(NA, 9)
    }
}

# Convert to dataframes
p_values_df <- as.data.frame(do.call(rbind, p_values_list))
test_statistic_df <- as.data.frame(do.call(rbind, test_statistic_list))

# Write to CSV
write.csv(p_values_df, "bitcoin_p_values_election.csv", row.names = TRUE)
write.csv(test_statistic_df, file = "bitcoin_test_statistics_election.csv", row.names = TRUE)

# Show result
print(p_values_df)
print(test_statistic_df)
