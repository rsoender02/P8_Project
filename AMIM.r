  library(quantmod)
  library(tseries)
  library(lmtest)

  # Define parameters
  start_date <- as.Date("2014-01-01")
  end_date <- as.Date("2025-03-31")
  symbol <- "BTC-USD"
  output_file <- "X:/Credit Risk/MN/P8/btc_prices.csv"

  # Get historical data and check for errors
  data <- try(getSymbols(symbol, src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE))
  if (!inherits(data, "try-error")) {
    # Extract relevant data and save to CSV
    bitcoin_data <- data.frame(
      date = index(data),
      open = as.numeric(data[, 1]),
      high = as.numeric(data[, 2]),
      low = as.numeric(data[, 3]),
      close = as.numeric(data[, 4]),
      adjusted_close = as.numeric(Ad(data)),
      volume = as.numeric(data[, 5])
    )
    write.csv(bitcoin_data, output_file, row.names = FALSE)
    cat("Data saved to", output_file, "\n")
  } else {
    cat("Failed to download data for", symbol, "\n")
  }


  MIM.roll <- function(data.table, identity.col, Date.col, rollWindow, return.col, min.obs, max.lag, a, force) {
    data <- data.table::copy(data.table)
    data.table::setDT(data)
    data.table::setnames(data,
      old = c(return.col),
      new = c("RET")
    ) ## change column name
    data.table::setorderv(data, c(identity.col, Date.col)) # oder columns

    pb <- txtProgressBar(
      min = 0, style = 3,
      max = nrow(data) - rollWindow * data.table::uniqueN(data[, c(identity.col), with = F])
    ) ### SET progress Bar

    dt.dates <- data[, list(date.join = seq(base::as.Date(get(Date.col), format = "%Y-%m%-%d"),
      by = "-1 day", len = rollWindow
    )),
    by = c(identity.col, Date.col)
    ] ## create an extra table for join

    MIM <- data.table::merge.data.table(dt.dates, data,
      by.x = c(identity.col, "date.join"),
      by.y = c(identity.col, Date.col), all.x = T
    )
    remove(dt.dates)
    MIM <- na.omit(MIM)
    MIM <- MIM[,
      {
        if (length(RET) < min.obs | sd(RET, na.rm = T) == 0) {
          MIM <- as.numeric(NA)
          N <- as.numeric(NA)
        } else {
          f <- ar(RET, aic = TRUE, order.max = max.lag)
          if (length(f$ar) == 0) {
            setTxtProgressBar(pb, .GRP) ### Update the progress
            if (force == FALSE) {
              MIM <- 0.0
              N <- 0.0
            }
           else if (force == TRUE) {
            # If the AIC point zero lag then estimate an AR(1), the coeficient is not significant anyway
            # This will avoid zero in MIM and AMIM
            f <- ar(RET, aic = F, order.max = 1)
            g <- solve(t(chol(f$asy.var.coef))) %*% f$ar ### (P^-1)*beta  or L^-1'Beta. Chol give upper matrix U => t(u)=L
            g <- abs(g) * a ## absolute and scale
            MIM <- sum(g, na.rm = T) / (1 + sum(g, na.rm = T))
            N <- as.numeric(f$order)
          } } else {
            g <- solve(t(chol(f$asy.var.coef))) %*% f$ar ### (P^-1)*beta  or L^-1'Beta. Chol give upper matrix U => t(u)=L
            g <- abs(g) * a ## absolute and scale
            MIM <- sum(g, na.rm = T) / (1 + sum(g, na.rm = T))
            N <- as.numeric(f$order)
          }
        }
        setTxtProgressBar(pb, .GRP) ### Update the progress bar
        list(MIM = MIM, N = N)
      },
      by = c(identity.col, Date.col)
    ]
    return(MIM)
  }

  AMIM.roll <- function(data.table, identity.col, Date.col, rollWindow, return.col, min.obs, max.lag, force) {
    MIM. <- MIM.roll(
      data.table = data.table, identity.col = identity.col, Date.col = Date.col,
      rollWindow = rollWindow, return.col = return.col, min.obs = min.obs,
      max.lag = max.lag, a = 1, force
    ) # compute MIM and N. Force a=1 like in the Tran & Leivrik (2019) paper

    MIM. <- data.table::setDT(MIM.)

    CI <- data.table::copy(AMIM::CI) ## copy the CI data from AMIM package
    data.table::setDT(CI)

    MIM. <- data.table::merge.data.table(x = MIM., y = CI[a == 1, .(N, CI)], all.x = T, by.x = "N", by.y = "N", sort = F)

    MIM.$AMIM <- (MIM.$MIM - MIM.$CI) / (1 - MIM.$CI)

    return(MIM.)
  }

  library(ggplot2)
  data <- read.csv("X:/Credit Risk/MN/P8/btc_prices.csv")
  data$ticker <- rep("BTC", nrow(data))
  data$date <- as.Date(data$date)

data$return <- c(NA, diff(data$close) / head(data$close, -1)) #simple return
data$return <- c(NA, diff(log(data$close))) # log return

#Entire period
  AMIM <- AMIM.roll(
    data.table = data, identity.col = "ticker", rollWindow = 3848,
    Date.col = "date", return.col = "return", min.obs = 3848, max.lag = 3000, force = FALSE)

  AMIM[3848] # AMIM 0

AMIM

# Plot Bitcoin closing prices
price_plot <- ggplot(data, aes(x = date, y = close)) +
  geom_line(color = "#211A52") +
  labs(
    x = "Time",
    y = "Price (USD)"
  ) +
  theme_minimal()

ggsave("btc_price_plot.png", plot = price_plot, width = 10, height = 4, dpi = 300)

# Plot Bitcoin returns
return_plot <- ggplot(data, aes(x = date, y = return)) +
  geom_line(color = "#211A52") +
  labs(
    x = "Time",
    y = "Daily Log Return"
  ) +
  theme_minimal()

ggsave("log_btc_return_plot.png", plot = return_plot, width = 10, height = 4, dpi = 300)

#365 rolling window
  AMIM <- AMIM.roll(
    data.table = data, identity.col = "ticker", rollWindow = 365,
    Date.col = "date", return.col = "return", min.obs = 364, max.lag = 30, force=TRUE
  ) 

unique(AMIM$N)

  k <- 30
  AMIM$moving_avg <- filter(AMIM$AMIM, rep(1/k, k), sides = 1)

AMIM <- na.omit(AMIM)

  p <- ggplot(AMIM, aes(x = date, y = N)) +
    geom_line(color = "#211A52") +
    labs(
      x = expression(t),
    y = ("Number of Lags")
    ) +
    theme_minimal()+
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
  ggsave("number_of_lags_used_365_no_force.png", plot = p, width = 10, height = 4, dpi = 300)

write.csv(AMIM, "X:/Credit Risk/MN/P8/amim_365_data_force.csv", row.names = FALSE)

  p <- ggplot(AMIM, aes(x = date, y = moving_avg)) +
    geom_line(color = "#211A52") +
    labs(
      x = expression(t),
    y = expression(italic(AMIM)[italic(t)*","*italic(365)] ~ "— 30-day MA")
    ) +
    theme_minimal()+
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
  ggsave("amim_moving_average_plot_force.png", plot = p, width = 10, height = 4, dpi = 300)

  p <- ggplot(AMIM, aes(x = date, y = AMIM)) +
    geom_line(color = "#211A52") +
    labs(
      x = "Time",
    y = expression(italic(AMIM)[italic(t)*","*italic(w)])
    ) +
    theme_minimal()+
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

  ggsave("amim_plot_force_log_365.png", plot = p, width = 10, height = 4, dpi = 300)

  p <- ggplot(AMIM, aes(x = date, y = MIM)) +
    geom_line(color = "#211A52") +
    labs(
      x = "Time",
    y = expression(italic(AMIM)[italic(t)*","*italic(w)])
    ) +
    theme_minimal()+
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
  ggsave("mim_plot_force.png", plot = p, width = 10, height = 4, dpi = 300)


#60 rolling window
data2021_2025 <- data[data$date >= "2021-01-01" & data$date <= "2025-03-31", ]

AMIM <- AMIM.roll(
    data.table = data2021_2025, identity.col = "ticker", rollWindow = 60,
    Date.col = "date", return.col = "return", min.obs = 59, max.lag = 10, force=TRUE
  )

unique(AMIM$N)

library(tidyr)
library(dplyr)

AMIM$btc_price <- log(data2021_2025$close[data2021_2025$date %in% AMIM$date])
AMIM <- na.omit(AMIM)

# Get min and max for each series
amim_min <- min(AMIM$AMIM, na.rm = TRUE)
amim_max <- max(AMIM$AMIM, na.rm = TRUE)

btc_min <- min(AMIM$btc_price, na.rm = TRUE)
btc_max <- max(AMIM$btc_price, na.rm = TRUE)

# Scale BTC linearly to match AMIM range
btc_scale <- (amim_max - amim_min) / (btc_max - btc_min)
btc_shift <- amim_min - btc_min * btc_scale
AMIM$btc_rescaled <- AMIM$btc_price * btc_scale + btc_shift

p <- ggplot(AMIM, aes(x = date)) +
  geom_line(aes(y = AMIM), color = "#211A52", size = 0.6) +
  geom_line(aes(y = btc_rescaled), color = "#CC445B", size = 0.6) +
  scale_y_continuous(
    name = expression(italic(AMIM)[italic(t)*","*italic(60)]),
    sec.axis = sec_axis(
      trans = ~ (. - btc_shift) / btc_scale,
      name = "Bitcoin Price (Log USD)"
    )
  ) +
  labs(x = "t") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(color = "#211A52"),
    axis.title.y.right = element_text(color = "#CC445B"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave("log_return_amim_plot_no_force_60_rolling_10_max_lag.png", plot = p, width = 10, height = 4, dpi = 300)


#14 rolling window
data2021_2025 <- data[data$date >= "2021-01-01" & data$date <= "2025-03-31", ]

AMIM <- AMIM.roll(
    data.table = data2021_2025, identity.col = "ticker", rollWindow = 60,
    Date.col = "date", return.col = "return", min.obs = 59, max.lag = 15, force=TRUE
  )

unique(AMIM$N)

library(tidyr)
library(dplyr)

AMIM$btc_price <- log(data2021_2025$close[data2021_2025$date %in% AMIM$date])
AMIM <- na.omit(AMIM)

  p <- ggplot(AMIM, aes(x = date, y = N)) +
    geom_line(color = "#211A52") +
    labs(
      x = expression(t),
    y = ("Number of Lags")
    ) +
    theme_minimal()+
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
  ggsave("number_of_lags_used_60_force.png", plot = p, width = 10, height = 4, dpi = 300)

write.csv(AMIM, "X:/Credit Risk/MN/P8/amim_365_data_force.csv", row.names = FALSE)

# Get min and max for each series
amim_min <- min(AMIM$AMIM, na.rm = TRUE)
amim_max <- max(AMIM$AMIM, na.rm = TRUE)

btc_min <- min(AMIM$btc_price, na.rm = TRUE)
btc_max <- max(AMIM$btc_price, na.rm = TRUE)

# Scale BTC linearly to match AMIM range
btc_scale <- (amim_max - amim_min) / (btc_max - btc_min)
btc_shift <- amim_min - btc_min * btc_scale
AMIM$btc_rescaled <- AMIM$btc_price * btc_scale + btc_shift

p <- ggplot(AMIM, aes(x = date)) +
  geom_line(aes(y = AMIM), color = "#211A52", size = 0.6) +
  geom_line(aes(y = btc_rescaled), color = "#CC445B", size = 0.6) +
  scale_y_continuous(
    name = expression(italic(AMIM)[italic(t)*","*italic(7)]),
    sec.axis = sec_axis(
      trans = ~ (. - btc_shift) / btc_scale,
      name = "Bitcoin Price (Log USD)"
    )
  ) +
  labs(x = "t") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(color = "#211A52"),
    axis.title.y.right = element_text(color = "#CC445B"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave("log_return_amim_plot_force_7_rolling_7_max_lag.png", plot = p, width = 10, height = 4, dpi = 300)

data2024_long_boi <- data[data$date >= "2024-09-08" & data$date <= "2024-11-07", ]
data2024_oct <- data[data$date >= "2024-10-06" & data$date <= "2024-11-06", ]
data2024_nov <- data[data$date >= "2024-11-06" & data$date <= "2024-12-06", ]
data2024_dec <- data[data$date >= "2024-12-06" & data$date <= "2025-01-06", ]

dim(data2024_long_boi)
dim(data2024_nov)

#data2024_oct
  AMIM_long_boi <- AMIM.roll(
    data.table = data2024_long_boi, identity.col = "ticker", rollWindow = 60,
    Date.col = "date", return.col = "return", min.obs = 60, max.lag = 59, force = TRUE)

  AMIM_long_boi[60]

#data2024_oct
  AMIM_oct <- AMIM.roll(
    data.table = data2024_oct, identity.col = "ticker", rollWindow = 31,
    Date.col = "date", return.col = "return", min.obs = 31, max.lag = 30, force = TRUE)

  AMIM_oct[31]

  #data2024_nov
  AMIM_nov <- AMIM.roll(
    data.table = data2024_nov, identity.col = "ticker", rollWindow = 30,
    Date.col = "date", return.col = "return", min.obs = 30, max.lag = 29, force = TRUE)

  AMIM_nov[30] 

  #data2024_dec
  AMIM_dec <- AMIM.roll(
    data.table = data2024_dec, identity.col = "ticker", rollWindow = 31,
    Date.col = "date", return.col = "return", min.obs = 31, max.lag = 30, force = TRUE)

  AMIM_oct[31]
  AMIM_nov[30] 
  AMIM_dec[31] 


dim(data2024_oct)

dim(data_before_covid)
# BEFORE COVID
data_before_covid <- data[data$date >= "2019-11-01" & data$date <= "2020-01-31", ]
  AMIM_before_covid <- AMIM.roll(
    data.table = data_before_covid, identity.col = "ticker", rollWindow = 91,
    Date.col = "date", return.col = "return", min.obs = 91, max.lag = 90, force = TRUE)

  AMIM_before_covid[91] 
dim(data_during_covid)
# DURING COVID
data_during_covid <- data[data$date >= "2020-03-01" & data$date <= "2020-05-01", ]
  AMIM_during_covid <- AMIM.roll(
    data.table = data_during_covid, identity.col = "ticker", rollWindow = 61,
    Date.col = "date", return.col = "return", min.obs = 61, max.lag = 60, force = TRUE)
AMIM_during_covid[62] 
# AFTER COVID
data_after_covid <- data[data$date >= "2021-03-01" & data$date <= "2021-05-01", ]
AMIM_after <- AMIM.roll(
  data.table = data_after_covid, identity.col = "ticker", rollWindow = nrow(unique(data_after_covid["date"])),
  Date.col = "date", return.col = "return", min.obs = nrow(unique(data_after_covid["date"])),
  max.lag = nrow(unique(data_after_covid["date"])) - 1, force = TRUE)
AMIM_after[nrow(AMIM_after)]

