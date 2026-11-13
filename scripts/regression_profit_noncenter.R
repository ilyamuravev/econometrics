# Регрессии только по фирмам ВНЕ центров (не в Краснодаре, Сочи и Новороссийске).
# Разность доходов и расходов на работника, млн руб.; Краснодарский край, 2023.
#
# Зачем: в полной выборке 443 из 832 фирм зарегистрированы в самих центрах, где
# расстояние равно 0. Тогда коэффициент при log(1 + км) во многом сравнивает
# «центр против остальных», а не показывает, как результат меняется по мере
# удаления. Здесь центры исключены, и связь оценивается только по изменению
# расстояния между нецентральными пунктами. Для сравнения та же связь оценивается
# на полной выборке с индикаторами центров: они «забирают» скачок между центром
# и периферией, и коэффициент расстояния тоже опирается только на нецентральные пункты.
#
# Спецификации, зависимая переменная и кластерные ошибки — как в
# scripts/regression_profit_road.R; каждая модель — с прямым и дорожным расстоянием.
#
# Запуск из корня проекта (нужен только базовый R):
#   source("scripts/regression_profit_noncenter.R")
#
# Результаты (results/):
#   profit_noncenter_sample.csv       — выборки, число фирм и пунктов
#   profit_noncenter_coefficients.csv — все коэффициенты, кластерные ошибки, 95%-е интервалы, p
#   profit_noncenter_summary.csv      — N, число кластеров, R², исключённая фирма
#   profit_noncenter_centers.csv      — коэффициент расстояния для розницы, автосервиса и их разница
#   profit_noncenter_comparison.csv   — полная выборка и выборка вне центров рядом
# и график figures/11_noncenter_vs_full.png.

d0 <- read.csv('data/krasnodar_2023_regression.csv', stringsAsFactors = FALSE)
d0$profit_per_worker <- (d0$revenue - d0$expenditure) / d0$employees_count / 1e6
centers <- c('krasnodar', 'sochi', 'novorossiysk', 'nearest')
types <- c(straight = 'dist', road = 'road')


# ---- 1. Выборки ---------------------------------------------

vars <- c('profit_per_worker', paste0('ln_dist_', centers), paste0('ln_road_', centers),
          'auto', 'small', 'ln_settl_pop', 'settlement_id')
d <- d0[complete.cases(d0[vars]) & is.finite(d0$ln_settl_pop), ]   # общая выборка, как в дорожных моделях
stopifnot(!anyDuplicated(d$tin), all(d$employees_count > 0))

# Индикаторы регистрации в самом центре (расстояние 0 и по прямой, и по дорогам)
d$in_krasnodar    <- as.integer(d$dist_krasnodar == 0)
d$in_sochi        <- as.integer(d$dist_sochi == 0)
d$in_novorossiysk <- as.integer(d$dist_novorossiysk == 0)
stopifnot(all(d$in_krasnodar + d$in_sochi + d$in_novorossiysk == as.integer(d$dist_nearest == 0)),
          all((d$road_nearest == 0) == (d$dist_nearest == 0)))

nc <- d[d$dist_nearest > 0, ]                                         # фирмы вне центров
q <- quantile(nc$profit_per_worker, c(.01, .99))                      # 1% хвостов — внутри этой выборки
nc_trimmed <- nc[nc$profit_per_worker >= q[1] & nc$profit_per_worker <= q[2], ]

n_settl <- function(x) length(unique(x$settlement_id))
sample_log <- data.frame(
  step = c('common_sample', 'excluded_in_krasnodar', 'excluded_in_sochi', 'excluded_in_novorossiysk',
           'noncenter_sample', 'noncenter_trimmed'),
  firms = c(nrow(d), sum(d$in_krasnodar), sum(d$in_sochi), sum(d$in_novorossiysk), nrow(nc), nrow(nc_trimmed)),
  settlements = c(n_settl(d), 1, 1, 1, n_settl(nc), n_settl(nc_trimmed)),
  note = c('832 фирмы, как в regression_profit_road.R', '', '', '',
           'фирмы, зарегистрированные не в Краснодаре, Сочи и Новороссийске',
           sprintf('profit_per_worker в [%.3f; %.3f] млн руб. (1-й и 99-й процентили выборки вне центров)', q[1], q[2])))
print(sample_log, row.names = FALSE)


# ---- 2. Кластерные ошибки и линейные комбинации (как в regression_profit_road.R) ----

cluster_vcov <- function(model, cluster) {
  X <- model.matrix(model); u <- residuals(model)
  n <- nrow(X); k <- ncol(X); g <- length(unique(cluster))
  stopifnot(length(cluster) == n, g > 1, model$rank == k)
  scores <- rowsum(X * u, group = cluster, reorder = FALSE)
  bread <- solve(crossprod(X))
  (g / (g - 1)) * ((n - 1) / (n - k)) * bread %*% crossprod(scores) %*% bread
}
lincom <- function(model, V, L, g) {
  estimate <- sum(L * coef(model)); se <- sqrt(drop(t(L) %*% V %*% L)); critical <- qt(.975, g - 1)
  data.frame(estimate = estimate, se_cluster = se,
             conf_low = estimate - critical * se, conf_high = estimate + critical * se,
             p_value = 2 * pt(abs(estimate / se), g - 1, lower.tail = FALSE))
}


# ---- 3. Спецификации ----------------------------------------

center_dummies <- 'in_krasnodar + in_sochi + in_novorossiysk'
specs <- data.frame(
  model  = c('krasnodar', 'controls', 'interactions', 'nearest', 'trimmed',
             'full_center_dummies', 'full_center_dummies_interactions'),
  center = c('krasnodar', 'krasnodar', 'krasnodar', 'nearest', 'krasnodar', 'krasnodar', 'krasnodar'),
  rhs    = c('X', 'X + auto + small + ln_settl_pop', 'X * (auto + small) + ln_settl_pop',
             'X + auto + small + ln_settl_pop', 'X + auto + small + ln_settl_pop',
             paste('X + auto + small + ln_settl_pop +', center_dummies),
             paste('X * (auto + small) + ln_settl_pop +', center_dummies)),
  sample = c('noncenter', 'noncenter', 'noncenter', 'noncenter', 'noncenter_trimmed', 'full', 'full'),
  stringsAsFactors = FALSE)
for (center in centers) for (smp in c('noncenter', 'noncenter_trimmed', 'noncenter_no_top_cook')) {
  specs <- rbind(specs, data.frame(model = paste(center, smp, 'interactions', sep = '_'), center = center,
                                   rhs = 'X * (auto + small) + ln_settl_pop', sample = smp,
                                   stringsAsFactors = FALSE))
}
make_formula <- function(rhs, type, center) {
  as.formula(paste('profit_per_worker ~', gsub('X', paste0('ln_', type, '_', center), rhs, fixed = TRUE)))
}

# Самая влиятельная фирма (наибольшее расстояние Кука в любой из двух версий модели)
# исключается из обеих версий, чтобы выборки совпадали
top_cook <- list()
for (center in centers) {
  cook <- sapply(types, function(type)
    cooks.distance(lm(make_formula('X * (auto + small) + ln_settl_pop', type, center), data = nc)))
  i <- which.max(apply(cook, 1, max))
  top_cook[[center]] <- list(row = i, label = paste0(nc$tin[i], ' (', nc$settlement[i], ', ',
                                                     nc$industry[i], ', ', round(nc$profit_per_worker[i], 1), ' млн руб.)'))
}


# ---- 4. Оценивание ------------------------------------------

coef_rows <- list(); summary_rows <- list(); center_rows <- list()
for (i in seq_len(nrow(specs))) {
  sp <- specs[i, ]
  data <- switch(sp$sample, full = d, noncenter = nc, noncenter_trimmed = nc_trimmed,
                 noncenter_no_top_cook = nc[-top_cook[[sp$center]]$row, ])
  dropped <- if (sp$sample == 'noncenter_no_top_cook') top_cook[[sp$center]]$label else ''
  g <- n_settl(data)
  for (distance in names(types)) {
    type <- types[[distance]]
    m <- lm(make_formula(sp$rhs, type, sp$center), data = data)
    V <- cluster_vcov(m, data$settlement_id)
    k <- length(coef(m))
    for (j in seq_len(k)) {
      coef_rows[[length(coef_rows) + 1]] <- data.frame(
        distance = distance, model = sp$model, sample = sp$sample,
        term = sub(paste0('ln_', type, '_'), 'ln_distance_', names(coef(m))[j]),
        n = nobs(m), clusters = g, lincom(m, V, replace(rep(0, k), j, 1), g))
    }
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      distance = distance, model = sp$model, center = sp$center, sample = sp$sample,
      n = nobs(m), clusters = g, r_squared = summary(m)$r.squared, dropped_firm = dropped)
    if (grepl('interactions$', sp$model) && sp$model != 'interactions') {
      x <- paste0('ln_', type, '_', sp$center)
      for (group in c('retail', 'auto', 'difference')) {
        L <- setNames(rep(0, k), names(coef(m)))
        if (group != 'difference') L[x] <- 1
        if (group != 'retail') L[paste0(x, ':auto')] <- 1
        center_rows[[length(center_rows) + 1]] <- data.frame(
          distance = distance, model = sp$model, center = sp$center, sample = sp$sample, group = group,
          n = nobs(m), clusters = g, lincom(m, V, L, g))
      }
    }
  }
}
coefficients <- do.call(rbind, coef_rows)
summary_table <- do.call(rbind, summary_rows)
center_table <- do.call(rbind, center_rows)


# ---- 5. Сравнение с полной выборкой --------------------------

# Полная выборка — из results/profit_road_comparison.csv (scripts/regression_profit_road.R)
full <- read.csv('results/profit_road_comparison.csv', stringsAsFactors = FALSE)
pick_full <- function(model, group) full[full$model == model & full$group == group, ]
pick_nc <- function(model, group, term = NULL) {
  if (is.null(term)) {
    center_table[center_table$model == model & center_table$group == group, ]
  } else {
    x <- coefficients[coefficients$model == model & coefficients$term == term, ]
    x$group <- group; x
  }
}
rows <- list(
  list('Краснодар, с контролями', pick_full('controls', 'all_firms'), pick_nc('controls', 'all_firms', 'ln_distance_krasnodar')),
  list('Ближайший центр, с контролями', pick_full('nearest', 'all_firms'), pick_nc('nearest', 'all_firms', 'ln_distance_nearest')),
  list('Краснодар: розница', pick_full('krasnodar_all_interactions', 'retail'), pick_nc('krasnodar_noncenter_interactions', 'retail')),
  list('Краснодар: автосервис', pick_full('krasnodar_all_interactions', 'auto'), pick_nc('krasnodar_noncenter_interactions', 'auto')),
  list('Краснодар: автосервис − розница', pick_full('krasnodar_all_interactions', 'difference'), pick_nc('krasnodar_noncenter_interactions', 'difference')),
  list('Сочи: розница', pick_full('sochi_all_interactions', 'retail'), pick_nc('sochi_noncenter_interactions', 'retail')),
  list('Сочи: автосервис − розница', pick_full('sochi_all_interactions', 'difference'), pick_nc('sochi_noncenter_interactions', 'difference')),
  list('Новороссийск: розница', pick_full('novorossiysk_all_interactions', 'retail'), pick_nc('novorossiysk_noncenter_interactions', 'retail')),
  list('Ближайший центр: автосервис − розница', pick_full('nearest_all_interactions', 'difference'), pick_nc('nearest_noncenter_interactions', 'difference')))
comparison <- do.call(rbind, lapply(rows, function(r) {
  f <- r[[2]]; nc_est <- r[[3]]
  do.call(rbind, lapply(names(types), function(distance) {
    x <- nc_est[nc_est$distance == distance, ]
    stopifnot(nrow(f) == 1, nrow(x) == 1)
    data.frame(estimate_name = r[[1]], distance = distance,
               full_estimate = f[[paste0(distance, '_estimate')]],
               full_conf_low = f[[paste0(distance, '_conf_low')]],
               full_conf_high = f[[paste0(distance, '_conf_high')]],
               full_p_value = f[[paste0(distance, '_p_value')]], full_n = f$n,
               noncenter_estimate = x$estimate, noncenter_conf_low = x$conf_low,
               noncenter_conf_high = x$conf_high, noncenter_p_value = x$p_value,
               noncenter_n = x$n, noncenter_clusters = x$clusters)
  }))
}))
# Модель на полной выборке с индикаторами центров — для сверки с выборкой вне центров
dummies <- coefficients[coefficients$model == 'full_center_dummies' & coefficients$term == 'ln_distance_krasnodar', ]


# ---- 6. Сохранение и вывод ----------------------------------

write.csv(sample_log, 'results/profit_noncenter_sample.csv', row.names = FALSE)
write.csv(coefficients, 'results/profit_noncenter_coefficients.csv', row.names = FALSE)
write.csv(summary_table, 'results/profit_noncenter_summary.csv', row.names = FALSE)
write.csv(center_table, 'results/profit_noncenter_centers.csv', row.names = FALSE)
write.csv(comparison, 'results/profit_noncenter_comparison.csv', row.names = FALSE)

cat('\nСамые влиятельные фирмы вне центров (исключаются в noncenter_no_top_cook):\n')
for (center in centers) cat(sprintf('  %-13s %s\n', center, top_cook[[center]]$label))
cat('\nПолная выборка против выборки вне центров, млн руб. на работника на единицу log(1 + км):\n')
show <- comparison
show[sapply(show, is.double)] <- lapply(show[sapply(show, is.double)], round, 3)
print(show, row.names = FALSE)
cat('\nКраснодар, с контролями, полная выборка + индикаторы центров:\n')
print(dummies[, c('distance', 'n', 'estimate', 'conf_low', 'conf_high', 'p_value')], row.names = FALSE, digits = 3)


# ---- 7. График: полная выборка и выборка вне центров (расстояние по дорогам) ----

ink <- '#0b0b0b'; ink2 <- '#52514e'
p <- comparison[comparison$distance == 'road', ]
p$y <- rev(seq_len(nrow(p)))
colors <- c(full = '#2a78d6', noncenter = '#eb6834')     # проверены на различимость
shapes <- c(full = 19, noncenter = 15)
off <- 0.16

png('figures/11_noncenter_vs_full.png', width = 1900, height = 1150, res = 180)
par(mar = c(4.3, 1, 1, 1), oma = c(0, 15.5, 3.6, 0), col.axis = ink2, fg = ink2)
plot(NA, xlim = range(c(p$full_conf_low, p$full_conf_high, p$noncenter_conf_low, p$noncenter_conf_high, 0)),
     ylim = c(0.4, nrow(p) + 0.6), yaxt = 'n', ylab = '',
     xlab = 'Коэффициент при log(1 + км по дорогам), млн руб. на работника')
abline(v = 0, col = '#9A9993', lty = 2)
segments(p$full_conf_low, p$y + off, p$full_conf_high, p$y + off, col = colors['full'], lwd = 2)
points(p$full_estimate, p$y + off, pch = shapes['full'], col = colors['full'], cex = 0.95)
segments(p$noncenter_conf_low, p$y - off, p$noncenter_conf_high, p$y - off, col = colors['noncenter'], lwd = 2)
points(p$noncenter_estimate, p$y - off, pch = shapes['noncenter'], col = colors['noncenter'], cex = 0.95)
axis(2, at = p$y, labels = p$estimate_name, las = 1, tick = FALSE, outer = TRUE, cex.axis = 0.85, col.axis = ink)
legend('topright', legend = c(sprintf('все фирмы (%d)', p$full_n[1]),
                                sprintf('вне центров (%d)', p$noncenter_n[1])),
       col = colors, pch = shapes, lwd = 2, bty = 'n', cex = 0.85, text.col = ink)
mtext('Связь с расстоянием: все фирмы и фирмы вне трёх центров',
      side = 3, line = 2, cex = 1, col = ink, font = 2, outer = TRUE)
mtext('Вне центров — фирмы не из Краснодара, Сочи и Новороссийска; 95%-е интервалы',
      side = 3, line = 0.7, cex = 0.8, col = ink2, outer = TRUE)
dev.off()
