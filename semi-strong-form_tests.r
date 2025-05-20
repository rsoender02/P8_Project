library(vrtest)
library(tseries)
library(pracma)
library(trend)
library(DescTools)
library(lmtest)
library(tidyverse)
library(MASS)

# Use either of (1), (2), (3), or (4)

# (1) Constant-mean-return model
data <- read.csv("path to file.csv")
data$date <- as.Date(data$date)
data$return <- c(NA, diff(log(data$close)))
event_date <- as.Date("2024-11-06")

# (2) SP500
data <- read.csv("path to file.csv")
data$date <- as.Date(data$date)
data_market <- read.csv("path to file")
data_market$date <- as.Date(data_market$date)
data <- data %>% semi_join(data_market, by = "date")
data$return <- c(NA, diff(log(data$close)))
data_market$return <- c(NA, diff(log(data_market$close)))
event_date <- as.Date("2024-11-06")

# (3) CRYPTOCAP_TOTAL2
data <- read.csv("path to file")
data$date <- as.Date(data$date)
data_market <- read.csv("path to file")
names(data_market)[names(data_market) == "time"] <- "date"
data_market$date <- as.Date(data_market$date)
data$return <- c(NA, diff(log(data$close)))
data_market$return <- c(NA, diff(log(data_market$close)))
event_date <- as.Date("2024-11-06")

# (4) CAPM model - use either SP500 or CRYPTOCAP_TOTAL2
data_market <- read.csv("path to file")
data_market <- read.csv("path to file")

data <- read.csv("path to file")
data$date <- as.Date(data$date)
risk_free <- 0.04039 / 252
data_market$date <- as.Date(data_market$date)
data <- data %>% semi_join(data_market, by = "date")
data$return <- c(NA, diff(log(data$close)) - risk_free)
data_market$return <- c(NA, diff(log(data_market$close)) - risk_free)
event_date <- as.Date("2024-11-06")




# Constant-mean-return model
data_ref <- data %>%
  filter(date >= "2024-07-05" & date <= "2024-10-04")

data_ev <- data %>%
  filter(date >= "2024-11-05" & date <= "2024-11-07")

theta_hat <- lm(data_ref$return ~ 1)$coefficients
X <- rep(1, length(data_ref$return))
epsilon_hat <- data_ref$return - X * theta_hat
sigma_hat_squared <- sum(epsilon_hat^2) / (length(epsilon_hat) - 2)

X_ev <- rep(1, length(data_ev$return))
epsilon_hat_ev <- data_ev$return - X_ev * theta_hat
V <- sigma_hat_squared * (diag(length(X_ev)) + X_ev %*% solve(t(X_ev) %*% X_ev) %*% t(X_ev))

gamma <- rep(1, nrow(V))
CAR_hat <- t(gamma) %*% epsilon_hat_ev
variance_CAR_hat <- t(gamma) %*% V %*% gamma
SCAR_hat <- CAR_hat / sqrt(variance_CAR_hat)

p_value <- 2 * (1 - pnorm(abs(SCAR_hat)))
p_value

# CAR plot, constant-mean-return model
data_ev <- data %>%
  filter(date >= "2024-11-05" & date <= "2024-11-14")

theta_hat <- lm(data_ref$return ~ 1)$coefficients
X <- rep(1, length(data_ref$return))
epsilon_hat <- data_ref$return - X * theta_hat
sigma_hat_squared <- sum(epsilon_hat^2) / (length(epsilon_hat) - 2)

X_ev <- rep(1, length(data_ev$return))
epsilon_hat_ev <- data.frame(date = data_ev$date, CAR = cumsum(data_ev$return - X_ev * theta_hat))

CAR_plot <- ggplot(epsilon_hat_ev, aes(x = date, y = CAR * 100)) +
  geom_line(color = "#211A52", size = 1.2) +
  geom_vline(xintercept = as.numeric(event_date), linetype = "dashed", color = "#cc445b") +
  labs(x = "Date", y = "Cumulative abnormal return (%)") +
  scale_x_date(date_labels = "%d %b") +
  theme_minimal()

ggsave("CAR_plot.png", plot = CAR_plot, width = 10, height = 4, dpi = 300)




# Right-extended event window CAR - constant-mean-return model
ref_start <- as.Date("2024-07-05")
ref_end <- as.Date("2024-10-04")
data_ref <- data %>% filter(date >= ref_start & date <= ref_end)

# Fit market model
fit <- lm(data_ref$return ~ 1)
theta_hat <- coef(fit)
epsilon_ref <- resid(fit)
sigma2_hat <- sum(epsilon_ref^2) / (nrow(data_ref) - length(theta_hat))

# 2. Define ev‐start and maximum ev‐end
ev_start <- as.Date("2024-11-05")
max_end <- as.Date("2024-11-14")

# 3. Build sequence of end‐dates and prepare results container
end_dates <- seq(ev_start, max_end, by = "day")
results <- data.frame(
  end_date = end_dates,
  CAR      = NA_real_,
  SCAR     = NA_real_,
  p_value  = NA_real_
)

# 4. Loop over each candidate ev window
for (i in (seq_along(end_dates))) {
  # Subset data for this window
  this_end <- end_dates[i]
  data_ev <- data %>% filter(date >= ev_start & date <= this_end)

  # Design matrix and residuals for the ev window
  X_ev <- rep(1, length(data_ev$return))
  eps_ev <- data_ev$return - X_ev * theta_hat

  # Variance–covariance of eps_ev
  X <- rep(1, length(data_ref$return))
  V <- sigma2_hat * (diag(length(X_ev)) +
    X_ev %*% solve(t(X) %*% X) %*% t(X_ev))

  # CAR and SCAR
  gamma <- rep(1, nrow(V))
  CAR_hat <- as.numeric(t(gamma) %*% eps_ev)
  var_CAR <- as.numeric(t(gamma) %*% V %*% gamma)
  SCAR_hat <- CAR_hat / sqrt(var_CAR)

  # Two‐sided p‐value
  p_val <- 2 * (1 - pnorm(abs(SCAR_hat)))

  # Store
  results$CAR[i] <- CAR_hat
  results$SCAR[i] <- SCAR_hat
  results$p_value[i] <- p_val
}

# 5. Inspect the table of results
print(results)
write.csv(results, file = "bitcoin_eventstudy_results_constant_mean_return_right_extended.csv", row.names = FALSE)



# Rolling event window CAR - constant-mean-return model
# 1) Estimate θ̂ and σ̂² on the reference window
ref_start <- as.Date("2024-07-05")
ref_end <- as.Date("2024-10-04")

data_ref <- data %>% filter(date >= ref_start & date <= ref_end)

fit <- lm(data_ref$return ~ 1)

theta_hat <- coef(fit)
eps_ref <- resid(fit)
sigma2_hat <- sum(eps_ref^2) / (nrow(data_ref) - length(theta_hat))

# 2) Rolling-window parameters
ev_start <- as.Date("2024-11-05")
ev_max <- as.Date("2024-11-14")
L <- 1 # window length in days
step_dates <- seq(ev_start, ev_max - (L - 1), by = "day")

# Prepare storage
results <- data.frame(
  window_start = step_dates,
  window_end   = step_dates + (L - 1),
  CAR          = numeric(length(step_dates)),
  SCAR         = numeric(length(step_dates)),
  p_value      = numeric(length(step_dates))
)

# 3) Loop over each rolling window
for (i in seq_along(step_dates)) {
  ws <- step_dates[i]
  we <- results$window_end[i]

  ev <- data %>% filter(date >= ws & date <= we)
  n <- nrow(ev)

  # residuals on this window
  X_ev <- rep(1, length(ev$return))
  eps_ev <- as.vector(ev$return - X_ev * theta_hat)

  # invertible XtX or fallback
  X <- rep(1, length(data_ref$return))
  XtX <- t(X) %*% X
  invX <- tryCatch(
    solve(XtX),
    error = function(e) ginv(XtX)
  )

  # covariance of epsilons
  V <- sigma2_hat * (diag(n) + X_ev %*% invX %*% t(X_ev))

  # scalar CAR, SCAR, and p
  CAR_val <- sum(eps_ev)
  varCAR <- sum(V) # since γ=1
  SCAR_val <- CAR_val / sqrt(varCAR)
  p_val <- 2 * (1 - pnorm(abs(SCAR_val)))

  results$CAR[i] <- CAR_val
  results$SCAR[i] <- SCAR_val
  results$p_value[i] <- p_val
}

print(results)
write.csv(results, file = "bitcoin_eventstudy_results_constant_mean_return1day.csv", row.names = FALSE)



# Market model - use either SP500, CRYPTOCAP_TOTAL2 or CAPM for either
data <- data %>%
  filter(date >= "2024-07-05" & date <= "2025-03-31")
data_market <- data_market %>%
  filter(date >= "2024-07-05" & date <= "2025-03-31")

data_ref <- data %>%
  filter(date >= "2024-07-05" & date <= "2024-10-04")
data_market_ref <- data_market %>%
  filter(date >= "2024-07-05" & date <= "2024-10-04")

data_ev <- data %>%
  filter(date >= "2024-11-05" & date <= "2024-11-07")
data_market_ev <- data_market %>%
  filter(date >= "2024-11-05" & date <= "2024-11-07")

theta_hat <- lm(data_ref$return ~ data_market_ref$return)$coefficients
X <- cbind(1, data_market_ref$return)
epsilon_hat <- data_ref$return - X %*% theta_hat
sigma_hat_squared <- sum(epsilon_hat^2) / (length(epsilon_hat) - 2)

X_ev <- cbind(1, data_market_ev$return)
epsilon_hat_ev <- data_ev$return - X_ev %*% theta_hat
V <- sigma_hat_squared * (diag(nrow(X_ev)) + X_ev %*% solve(t(X_ev) %*% X_ev) %*% t(X_ev))

gamma <- rep(1, nrow(V))
CAR_hat <- t(gamma) %*% epsilon_hat_ev
variance_CAR_hat <- t(gamma) %*% V %*% gamma
SCAR_hat <- CAR_hat / sqrt(variance_CAR_hat)

p_value <- 2 * (1 - pnorm(abs(SCAR_hat)))
p_value


# CAR plot, market model
data_ev <- data %>%
  filter(date >= "2024-11-05" & date <= "2024-11-14")
data_market_ev <- data_market %>%
  filter(date >= "2024-11-05" & date <= "2024-11-14")

theta_hat <- lm(data_ref$return ~ data_market_ref$return)$coefficients
X <- cbind(1, data_market_ref$return)
epsilon_hat <- data_ref$return - X %*% theta_hat
sigma_hat_squared <- sum(epsilon_hat^2) / (length(epsilon_hat) - 2)

X_ev <- cbind(1, data_market_ev$return)
epsilon_hat_ev <- data.frame(date = data_ev$date, CAR = cumsum(data_ev$return - X_ev %*% theta_hat))

CAR_plot <- ggplot(epsilon_hat_ev, aes(x = date, y = CAR * 100)) +
  geom_line(color = "#211A52", size = 1.2) +
  geom_vline(xintercept = as.numeric(event_date), linetype = "dashed", color = "#cc445b") +
  labs(x = "Date", y = "Cumulative abnormal return (%)") +
  scale_x_date(date_labels = "%d %b") +
  theme_minimal()

ggsave("CAR_plot.png", plot = CAR_plot, width = 10, height = 4, dpi = 300)




# Right-extended event window CAR - market model
ref_start <- as.Date("2024-07-05")
ref_end <- as.Date("2024-10-04")
data_ref <- data %>% filter(date >= ref_start & date <= ref_end)
data_market_ref <- data_market %>% filter(date >= ref_start & date <= ref_end)

# Fit market model
fit <- lm(data_ref$return ~ data_market_ref$return)
theta_hat <- coef(fit)
epsilon_ref <- resid(fit)
sigma2_hat <- sum(epsilon_ref^2) / (nrow(data_ref) - length(theta_hat))

# 2. Define ev‐start and maximum ev‐end
ev_start <- as.Date("2024-11-05")
max_end <- as.Date("2024-11-14")

# 3. Build sequence of end‐dates and prepare results container
end_dates <- seq(ev_start, max_end, by = "day")
# Use the following line if using S&P 500
end_dates <- c(as.Date("2024-11-05"), as.Date("2024-11-06"), as.Date("2024-11-07"), as.Date("2024-11-08"), as.Date("2024-11-11"), as.Date("2024-11-12"), as.Date("2024-11-13"), as.Date("2024-11-14"))
results <- data.frame(
  end_date = end_dates,
  CAR      = NA_real_,
  SCAR     = NA_real_,
  p_value  = NA_real_
)

# 4. Loop over each candidate ev window
for (i in (seq_along(end_dates))) {
  # Subset data for this window
  this_end <- end_dates[i]
  data_ev <- data %>% filter(date >= ev_start & date <= this_end)
  market_ev <- data_market %>% filter(date >= ev_start & date <= this_end)

  # Design matrix and residuals for the ev window
  X_ev <- cbind(1, market_ev$return)
  eps_ev <- data_ev$return - X_ev %*% theta_hat
  X <- cbind(1, data_market_ref$return)

  # Variance–covariance of eps_ev
  if (this_end == ev_start) {
    V <- sigma2_hat * (diag(nrow(X_ev)))
  } else {
    V <- sigma2_hat * (diag(nrow(X_ev)) +
      X_ev %*% solve(t(X) %*% X) %*% t(X_ev))
  }

  # CAR and SCAR
  gamma <- rep(1, nrow(V))
  CAR_hat <- (t(gamma) %*% eps_ev)
  var_CAR <- as.numeric(t(gamma) %*% V %*% gamma)
  SCAR_hat <- CAR_hat / sqrt(var_CAR)

  # Two‐sided p‐value
  p_val <- 2 * (1 - pnorm(abs(SCAR_hat)))

  # Store
  results$CAR[i] <- CAR_hat
  results$SCAR[i] <- SCAR_hat
  results$p_value[i] <- p_val
}

# 5. Inspect the table of results
print(results)
write.csv(results, file = "bitcoin_eventstudy_results_market_right_extended.csv", row.names = FALSE)




# Rolling event window CAR - market model
# 1) Estimate θ̂ and σ̂² on the reference window
ref_start <- as.Date("2024-07-05")
ref_end <- as.Date("2024-10-04")

data_ref <- data %>% filter(date >= ref_start & date <= ref_end)
mkt_ref <- data_market %>% filter(date >= ref_start & date <= ref_end)

fit <- lm(data_ref$return ~ mkt_ref$return)

theta_hat <- coef(fit)
eps_ref <- resid(fit)
sigma2_hat <- sum(eps_ref^2) / (nrow(data_ref) - length(theta_hat))

# 2) Rolling-window parameters
ev_start <- as.Date("2024-11-05")
ev_max <- as.Date("2024-11-14")
L <- 1 # window length in days
step_dates <- seq(ev_start, ev_max - (L - 1), by = "day")
# Use the following line if using S&P 500
step_dates <- c(as.Date("2024-11-05"), as.Date("2024-11-06"), as.Date("2024-11-07"), as.Date("2024-11-08"), as.Date("2024-11-11"), as.Date("2024-11-12"), as.Date("2024-11-13"), as.Date("2024-11-14"))

# Prepare storage
results <- data.frame(
  window_start = step_dates,
  window_end   = step_dates + (L - 1),
  CAR          = numeric(length(step_dates)),
  SCAR         = numeric(length(step_dates)),
  p_value      = numeric(length(step_dates))
)

# 3) Loop over each rolling window
for (i in seq_along(step_dates)) {
  ws <- step_dates[i]
  we <- results$window_end[i]

  ev <- data %>% filter(date >= ws & date <= we)
  mkt_ev <- data_market %>% filter(date >= ws & date <= we)
  n <- nrow(ev)

  # residuals on this window
  X_ev <- cbind(1, mkt_ev$return)
  eps_ev <- as.vector(ev$return - X_ev %*% theta_hat)

  # invertible XtX or fallback
  X <- cbind(1, mkt_ref$return)
  XtX <- t(X) %*% X
  invX <- tryCatch(
    solve(XtX),
    error = function(e) ginv(XtX)
  )

  # covariance of epsilons
  V <- sigma2_hat * (diag(n) + X_ev %*% invX %*% t(X_ev))

  # scalar CAR, SCAR, and p
  CAR_val <- sum(eps_ev)
  varCAR <- sum(V) # since γ=1
  SCAR_val <- CAR_val / sqrt(varCAR)
  p_val <- 2 * (1 - pnorm(abs(SCAR_val)))

  results$CAR[i] <- CAR_val
  results$SCAR[i] <- SCAR_val
  results$p_value[i] <- p_val
}
print(results)
write.csv(results, file = "bitcoin_eventstudy_results_market1day.csv", row.names = FALSE)






# Plot of windows
data_windows <- data %>%
  filter(date >= "2024-06-01" & date <= "2024-12-31")

ref_start <- as.Date("2024-07-05")
ref_end <- as.Date("2024-10-04")
ev_start <- as.Date("2024-11-05")
ev_end <- as.Date("2024-11-14")

event_windows <- ggplot(data_windows, aes(x = date, y = close)) +
  # Shaded areas
  annotate("rect",
    xmin = ref_start, xmax = ref_end,
    ymin = -Inf, ymax = Inf, alpha = 0.3, fill = "grey70"
  ) +
  annotate("rect",
    xmin = ev_start, xmax = ev_end,
    ymin = -Inf, ymax = Inf, alpha = 0.3, fill = "#BB5B17"
  ) +
  # Labels (centered and abbreviated)
  annotate("text",
    x = as.Date("2024-08-20"), y = Inf, label = "Estimation",
    vjust = 2, size = 3.5, fontface = "italic", color = "grey30"
  ) +
  annotate("text",
    x = as.Date("2024-11-09"), y = Inf, label = "Event",
    vjust = 2, size = 3.5, fontface = "italic", color = "#BB5B17"
  ) +
  # Main plot elements
  geom_line(color = "#211A52", size = 1.2) +
  geom_vline(xintercept = as.numeric(event_date), linetype = "dashed", color = "#cc445b") +
  labs(x = "Date", y = "Price (USD)") +
  theme_minimal() +
  scale_x_date(
    date_labels = "%b %Y",
    date_breaks = "1 month",
    limits = range(data_windows$date),
    expand = c(0.01, 0.01)
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

event_windows
ggsave("event_windows.png", plot = event_windows, width = 10, height = 4, dpi = 300)
