# Разность доходов и расходов на работника, млн рублей; Краснодарский край, 2023.
d <- read.csv('data/krasnodar_2023_regression.csv', stringsAsFactors = FALSE)
d$profit_per_worker <- (d$revenue - d$expenditure) / d$employees_count / 1e6
vars <- c('profit_per_worker', 'ln_dist_krasnodar', 'ln_dist_sochi', 'ln_dist_novorossiysk', 'ln_dist_nearest',
          'auto', 'small', 'ln_settl_pop', 'settlement_id')
d <- d[complete.cases(d[vars]) & is.finite(d$ln_settl_pop), ]
stopifnot(nrow(d) > 0, all(d$employees_count > 0), !anyDuplicated(d$tin))
d$population_group <- cut(d$settl_population, c(0, 10000, 100000, Inf),
                          right = FALSE, labels = c('under_10k', '10k_to_100k', '100k_plus'))
stopifnot(!anyNA(d$population_group))

# CR1: поправка на число кластеров и регрессоров; t-распределение с G-1 степенями свободы.
cluster_table <- function(model, cluster) {
  X <- model.matrix(model)
  u <- residuals(model)
  n <- nrow(X); k <- ncol(X); g <- length(unique(cluster))
  stopifnot(length(cluster) == n, g > 1, model$rank == k)
  scores <- rowsum(X * u, group = cluster, reorder = FALSE)
  stopifnot(max(abs(colSums(scores))) < 1e-7)
  bread <- solve(crossprod(X))
  V <- (g / (g - 1)) * ((n - 1) / (n - k)) * bread %*% crossprod(scores) %*% bread
  se <- sqrt(diag(V)); b <- coef(model); critical <- qt(.975, g - 1)
  data.frame(term = names(b), estimate = unname(b), se_cluster = unname(se),
             conf_low = unname(b - critical * se), conf_high = unname(b + critical * se),
             p_value = unname(2 * pt(abs(b / se), g - 1, lower.tail = FALSE)),
             row.names = NULL)
}

q <- quantile(d$profit_per_worker, c(.01, .99))
trimmed <- d[d$profit_per_worker >= q[1] & d$profit_per_worker <= q[2], ]
models <- list(
  krasnodar = lm(profit_per_worker ~ ln_dist_krasnodar, data = d),
  controls = lm(profit_per_worker ~ ln_dist_krasnodar + auto + small + ln_settl_pop, data = d),
  interactions = lm(profit_per_worker ~ ln_dist_krasnodar * (auto + small) + ln_settl_pop, data = d),
  nearest = lm(profit_per_worker ~ ln_dist_nearest + auto + small + ln_settl_pop, data = d),
  trimmed = lm(profit_per_worker ~ ln_dist_krasnodar + auto + small + ln_settl_pop, data = trimmed),
  population_groups = lm(profit_per_worker ~ ln_dist_krasnodar + auto + small + population_group, data = d),
  population_groups_trimmed = lm(profit_per_worker ~ ln_dist_krasnodar + auto + small + population_group, data = trimmed)
)
data_sets <- list(d, d, d, d, trimmed, d, trimmed)
names(data_sets) <- names(models)
for (center in c('krasnodar', 'sochi', 'novorossiysk', 'nearest')) {
  formula <- as.formula(paste('profit_per_worker ~', paste0('ln_dist_', center),
                              '* (auto + small) + ln_settl_pop'))
  for (sample in c('all', 'trimmed')) {
    name <- paste(center, sample, 'interactions', sep = '_')
    data_sets[[name]] <- if (sample == 'all') d else trimmed
    models[[name]] <- lm(formula, data = data_sets[[name]])
  }
}
# Проверка влияния фирмы с максимальным значением показателя.
data_sets$sochi_without_max <- d[-which.max(d$profit_per_worker), ]
models$sochi_without_max <- lm(profit_per_worker ~ ln_dist_sochi * (auto + small) + ln_settl_pop,
                               data = data_sets$sochi_without_max)
tables <- Map(function(m, data) cluster_table(m, data$settlement_id), models, data_sets)
result <- do.call(rbind, Map(function(tab, name) data.frame(model = name, tab), tables, names(tables)))
rownames(result) <- NULL
stats <- data.frame(model = names(models), n = sapply(models, nobs),
                    clusters = sapply(data_sets, function(x) length(unique(x$settlement_id))),
                    r_squared = sapply(models, function(x) summary(x)$r.squared), row.names = NULL)
write.csv(result, 'results/profit_regression_coefficients.csv', row.names = FALSE)
write.csv(stats, 'results/profit_regression_summary.csv', row.names = FALSE)
print(stats, row.names = FALSE)
print(result, row.names = FALSE, digits = 4)
cat('\nГраницы 1–99% для проверки устойчивости:', q, 'млн руб.\n')

# Каждый центр оценивается отдельно: расстояния до центров взаимосвязаны.
center_results <- list()
for (center in c('krasnodar', 'sochi', 'novorossiysk', 'nearest')) {
  distance <- paste0('ln_dist_', center)
  formula <- as.formula(paste('profit_per_worker ~', distance, '* (auto + small) + ln_settl_pop'))
  for (sample in c('all', 'trimmed')) {
    data <- if (sample == 'all') d else trimmed
    m <- models[[paste(center, sample, 'interactions', sep = '_')]]
    X <- model.matrix(m); u <- residuals(m)
    n <- nrow(X); k <- ncol(X); g <- length(unique(data$settlement_id))
    stopifnot(m$rank == k, g > 1)
    scores <- rowsum(X * u, data$settlement_id)
    bread <- solve(crossprod(X))
    V <- g / (g - 1) * (n - 1) / (n - k) * bread %*% crossprod(scores) %*% bread
    for (group in c('retail', 'auto', 'difference')) {
      L <- setNames(rep(0, k), names(coef(m)))
      if (group != 'difference') L[distance] <- 1
      if (group != 'retail') L[paste0(distance, ':auto')] <- 1
      estimate <- sum(L * coef(m)); se <- sqrt(drop(t(L) %*% V %*% L))
      center_results[[length(center_results) + 1]] <- data.frame(
        center = center, sample = sample, group = group, n = n, clusters = g,
        estimate = estimate, conf_low = estimate - qt(.975, g - 1) * se,
        conf_high = estimate + qt(.975, g - 1) * se,
        p_value = 2 * pt(abs(estimate / se), g - 1, lower.tail = FALSE))
    }
  }
}
center_results <- do.call(rbind, center_results)
write.csv(center_results, 'results/profit_regression_centers.csv', row.names = FALSE)
print(center_results, row.names = FALSE, digits = 4)
