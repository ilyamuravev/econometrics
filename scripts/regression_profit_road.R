# Регрессии с расстоянием ПО ДОРОГАМ и сравнение с расстоянием по прямой.
# Разность доходов и расходов на работника, млн руб.; Краснодарский край, 2023.
#
# Зависимая переменная и спецификации — те же, что в scripts/regression_profit.R.
# Каждая модель оценивается дважды: с прямыми (ln_dist_*) и с дорожными
# (ln_road_*) расстояниями на ОДНОЙ И ТОЙ ЖЕ выборке, поэтому различия между
# версиями вызваны только способом измерения расстояния. Базовый скрипт
# regression_profit.R и его результаты не меняются.
#
# Запуск из корня проекта (нужен только базовый R):
#   source("scripts/regression_profit_road.R")
#
# Результаты (results/):
#   profit_road_sample.csv       — выборки, число фирм и пунктов, причины исключений
#   profit_road_coefficients.csv — все коэффициенты, кластерные ошибки, 95%-е интервалы, p
#   profit_road_summary.csv      — N, число кластеров, R², исключённая фирма
#   profit_road_centers.csv      — коэффициент расстояния для розницы, автосервиса и их разница
#   profit_road_comparison.csv   — прямое и дорожное расстояние рядом (для слайдов)
# и график figures/08_road_vs_straight_coefficients.png.

d0 <- read.csv('data/krasnodar_2023_regression.csv', stringsAsFactors = FALSE)
d0$profit_per_worker <- (d0$revenue - d0$expenditure) / d0$employees_count / 1e6
centers <- c('krasnodar', 'sochi', 'novorossiysk', 'nearest')
types <- c(straight = 'dist', road = 'road')   # префиксы столбцов ln_dist_* и ln_road_*


# ---- 1. Проверка дорожных расстояний ------------------------

# Столбцы road_* добавляет scripts/road_distances.R из кэша data/road_distances.csv.
for (center in centers) {
  road <- d0[[paste0('road_', center)]]
  straight <- d0[[paste0('dist_', center)]]
  stopifnot(!is.null(road), !anyNA(road), all(road >= 0),
            max(abs(d0[[paste0('ln_road_', center)]] - log1p(road))) < 1e-9,
            all(road[straight == 0] == 0))   # фирмы в самом центре: 0 и по прямой, и по дорогам
  cat(sprintf('%-13s фирм в центре: %3d | маршрутов короче 0.95 прямой (>5 км): %d\n',
              center, sum(straight == 0), sum(straight > 5 & road < 0.95 * straight)))
}


# ---- 2. Общая выборка ---------------------------------------

vars <- c('profit_per_worker', paste0('ln_dist_', centers), paste0('ln_road_', centers),
          'auto', 'small', 'ln_settl_pop', 'settlement_id')
bad_pop <- !is.finite(d0$ln_settl_pop)               # население 0 -> логарифм не определён
bad_other <- !complete.cases(d0[vars]) & !bad_pop
d <- d0[!bad_pop & !bad_other, ]
stopifnot(nrow(d) > 0, all(d$employees_count > 0), !anyDuplicated(d$tin))
d$population_group <- cut(d$settl_population, c(0, 10000, 100000, Inf),
                          right = FALSE, labels = c('under_10k', '10k_to_100k', '100k_plus'))
stopifnot(!anyNA(d$population_group))

# Проверка чувствительности: без 1% крайних значений зависимой переменной с каждой стороны
q <- quantile(d$profit_per_worker, c(.01, .99))
trimmed <- d[d$profit_per_worker >= q[1] & d$profit_per_worker <= q[2], ]

n_settl <- function(x) length(unique(x$settlement_id))
sample_log <- data.frame(
  step = c('prepared_table', 'excluded_zero_population', 'excluded_missing_values',
           'common_sample', 'trimmed_sample'),
  firms = c(nrow(d0), sum(bad_pop), sum(bad_other), nrow(d), nrow(trimmed)),
  settlements = c(n_settl(d0), length(unique(d0$settlement_id[bad_pop])),
                  length(unique(d0$settlement_id[bad_other])), n_settl(d), n_settl(trimmed)),
  note = c('data/krasnodar_2023_regression.csv',
           paste(unique(d0$settlement[bad_pop]), collapse = '; '),
           paste(unique(d0$settlement[bad_other]), collapse = '; '),
           'одна и та же для прямых и дорожных моделей',
           sprintf('profit_per_worker в [%.3f; %.3f] млн руб. (1-й и 99-й процентили)', q[1], q[2])))
print(sample_log, row.names = FALSE)


# ---- 3. Кластерные стандартные ошибки (как в regression_profit.R) ----

# CR1: поправка на число кластеров и регрессоров; t-распределение с G-1 степенями свободы.
cluster_vcov <- function(model, cluster) {
  X <- model.matrix(model)
  u <- residuals(model)
  n <- nrow(X); k <- ncol(X); g <- length(unique(cluster))
  stopifnot(length(cluster) == n, g > 1, model$rank == k)
  scores <- rowsum(X * u, group = cluster, reorder = FALSE)
  stopifnot(max(abs(colSums(scores))) < 1e-7)
  bread <- solve(crossprod(X))
  (g / (g - 1)) * ((n - 1) / (n - k)) * bread %*% crossprod(scores) %*% bread
}

# Линейная комбинация коэффициентов L'b с интервалом и p-значением
lincom <- function(model, V, L, g) {
  estimate <- sum(L * coef(model)); se <- sqrt(drop(t(L) %*% V %*% L))
  critical <- qt(.975, g - 1)
  data.frame(estimate = estimate, se_cluster = se,
             conf_low = estimate - critical * se, conf_high = estimate + critical * se,
             p_value = 2 * pt(abs(estimate / se), g - 1, lower.tail = FALSE))
}


# ---- 4. Спецификации ----------------------------------------

# X — расстояние до указанного центра: ln_dist_<центр> или ln_road_<центр>
specs <- data.frame(
  model  = c('krasnodar', 'controls', 'interactions', 'nearest', 'trimmed',
             'population_groups', 'population_groups_trimmed'),
  center = c('krasnodar', 'krasnodar', 'krasnodar', 'nearest', 'krasnodar', 'krasnodar', 'krasnodar'),
  rhs    = c('X', 'X + auto + small + ln_settl_pop', 'X * (auto + small) + ln_settl_pop',
             'X + auto + small + ln_settl_pop', 'X + auto + small + ln_settl_pop',
             'X + auto + small + population_group', 'X + auto + small + population_group'),
  sample = c('all', 'all', 'all', 'all', 'trimmed', 'all', 'trimmed'),
  stringsAsFactors = FALSE)
# Каждый центр отдельно: расстояния до центров взаимосвязаны
for (center in centers) for (smp in c('all', 'trimmed', 'no_top_cook')) {
  specs <- rbind(specs, data.frame(model = paste(center, smp, 'interactions', sep = '_'),
                                   center = center, rhs = 'X * (auto + small) + ln_settl_pop',
                                   sample = smp, stringsAsFactors = FALSE))
}

make_formula <- function(rhs, type, center) {
  as.formula(paste('profit_per_worker ~', gsub('X', paste0('ln_', type, '_', center), rhs, fixed = TRUE)))
}

# Проверка влияния: без самой влиятельной фирмы (наибольшее расстояние Кука).
# Фирма выбирается по максимуму из двух версий модели и исключается из ОБЕИХ,
# чтобы прямая и дорожная оценки остались на одной выборке.
top_cook <- list()
for (center in centers) {
  cook <- sapply(types, function(type)
    cooks.distance(lm(make_formula('X * (auto + small) + ln_settl_pop', type, center), data = d)))
  i <- which.max(apply(cook, 1, max))
  top_cook[[center]] <- list(row = i, tin = d$tin[i], settlement = d$settlement[i],
                             same_in_both = all(apply(cook, 2, which.max) == i))
}


# ---- 5. Оценивание ------------------------------------------

coef_rows <- list(); summary_rows <- list(); center_rows <- list()
for (i in seq_len(nrow(specs))) {
  sp <- specs[i, ]
  data <- switch(sp$sample, all = d, trimmed = trimmed, no_top_cook = d[-top_cook[[sp$center]]$row, ])
  dropped <- if (sp$sample == 'no_top_cook') {
    paste0(top_cook[[sp$center]]$tin, ' (', top_cook[[sp$center]]$settlement, ')')
  } else ''
  g <- n_settl(data)
  for (distance in names(types)) {
    type <- types[[distance]]
    m <- lm(make_formula(sp$rhs, type, sp$center), data = data)
    V <- cluster_vcov(m, data$settlement_id)
    k <- length(coef(m))

    for (j in seq_len(k)) {                        # все коэффициенты модели
      L <- replace(rep(0, k), j, 1)
      coef_rows[[length(coef_rows) + 1]] <- data.frame(
        distance = distance, model = sp$model, sample = sp$sample,
        term = sub(paste0('ln_', type, '_'), 'ln_distance_', names(coef(m))[j]),
        lincom(m, V, L, g))
    }
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      distance = distance, model = sp$model, center = sp$center, sample = sp$sample,
      n = nobs(m), clusters = g, r_squared = summary(m)$r.squared, dropped_firm = dropped)

    # Отраслевые коэффициенты расстояния (для микропредприятий):
    # розница = b_dist, автосервис = b_dist + b_dist:auto, разница = b_dist:auto
    if (grepl('interactions$', sp$model) && sp$model != 'interactions') {
      x <- paste0('ln_', type, '_', sp$center)
      for (group in c('retail', 'auto', 'difference')) {
        L <- setNames(rep(0, k), names(coef(m)))
        if (group != 'difference') L[x] <- 1
        if (group != 'retail') L[paste0(x, ':auto')] <- 1
        center_rows[[length(center_rows) + 1]] <- data.frame(
          distance = distance, center = sp$center, sample = sp$sample, group = group,
          n = nobs(m), clusters = g, lincom(m, V, L, g))
      }
    }
  }
}
coefficients <- do.call(rbind, coef_rows)
summary_table <- do.call(rbind, summary_rows)
center_table <- do.call(rbind, center_rows)


# ---- 6. Таблица сравнения для слайдов -----------------------

# Изменение показателя при росте расстояния с 10 до 100 км:
# коэффициент x log((1 + 100) / (1 + 10)), млн руб. на работника
to_100 <- log(101 / 11)
key <- rbind(
  data.frame(block = 'krasnodar_models',
             coefficients[coefficients$term == 'ln_distance_krasnodar' &
                          !grepl('interactions', coefficients$model), ],
             center = 'krasnodar', group = 'all_firms')[, c('block', 'distance', 'model', 'center', 'sample', 'group',
                                                            'estimate', 'conf_low', 'conf_high', 'p_value')],
  data.frame(block = 'nearest_model',
             coefficients[coefficients$model == 'nearest' & coefficients$term == 'ln_distance_nearest', ],
             center = 'nearest', group = 'all_firms')[, c('block', 'distance', 'model', 'center', 'sample', 'group',
                                                          'estimate', 'conf_low', 'conf_high', 'p_value')],
  data.frame(block = 'industry_by_center', model = paste(center_table$center, center_table$sample, 'interactions', sep = '_'),
             center_table)[, c('block', 'distance', 'model', 'center', 'sample', 'group',
                               'estimate', 'conf_low', 'conf_high', 'p_value')])
id <- c('block', 'model', 'center', 'sample', 'group')
comparison <- merge(
  setNames(key[key$distance == 'straight', -2], c(id, paste0('straight_', c('estimate', 'conf_low', 'conf_high', 'p_value')))),
  setNames(key[key$distance == 'road', -2], c(id, paste0('road_', c('estimate', 'conf_low', 'conf_high', 'p_value')))),
  by = id, sort = FALSE)
comparison$straight_effect_10_100km <- comparison$straight_estimate * to_100
comparison$road_effect_10_100km <- comparison$road_estimate * to_100
comparison <- merge(comparison, unique(summary_table[summary_table$distance == 'road',
                                                     c('model', 'n', 'clusters', 'dropped_firm')]),
                    by = 'model', sort = FALSE)
comparison <- comparison[order(match(comparison$block, unique(key$block)),
                               match(comparison$center, centers),
                               match(comparison$sample, c('all', 'trimmed', 'no_top_cook')),
                               match(comparison$group, c('all_firms', 'retail', 'auto', 'difference'))),
                         c(id, 'n', 'clusters', grep('^straight_', names(comparison), value = TRUE),
                           grep('^road_', names(comparison), value = TRUE), 'dropped_firm')]
rownames(comparison) <- NULL


# ---- 7. Сохранение ------------------------------------------

write.csv(sample_log, 'results/profit_road_sample.csv', row.names = FALSE)
write.csv(coefficients, 'results/profit_road_coefficients.csv', row.names = FALSE)
write.csv(summary_table, 'results/profit_road_summary.csv', row.names = FALSE)
write.csv(center_table, 'results/profit_road_centers.csv', row.names = FALSE)
write.csv(comparison, 'results/profit_road_comparison.csv', row.names = FALSE)

cat('\nСамые влиятельные фирмы (исключаются в выборке no_top_cook):\n')
for (center in centers) {
  tc <- top_cook[[center]]
  cat(sprintf('  %-13s ИНН %s, %s; одна и та же в обеих версиях: %s\n',
              center, tc$tin, tc$settlement, tc$same_in_both))
}
cat('\nСравнение: коэффициент при log(1 + расстояние), млн руб. на работника\n')
show <- comparison
num <- sapply(show, is.numeric) & !names(show) %in% c('n', 'clusters')
show[num] <- lapply(show[num], round, 3)
print(show[, c('model', 'sample', 'group', 'n', 'straight_estimate', 'straight_conf_low', 'straight_conf_high',
               'road_estimate', 'road_conf_low', 'road_conf_high')], row.names = FALSE)


# ---- 8. График: прямое и дорожное расстояние рядом -----------

plot_data <- center_table[center_table$sample %in% c('all', 'trimmed'), ]
center_names <- c(krasnodar = 'Краснодар', sochi = 'Сочи', novorossiysk = 'Новороссийск',
                  nearest = 'Ближайший центр')
group_names <- c(retail = 'розница', auto = 'автосервис', difference = 'разница')
rows <- expand.grid(group = names(group_names), center = centers, stringsAsFactors = FALSE)
rows$label <- paste0(center_names[rows$center], ': ', group_names[rows$group])
rows$y <- nrow(rows) + 3 * 0.8 - (seq_len(nrow(rows)) - 1) - rep(0:3, each = 3) * 0.8
colors <- c(straight = '#2a78d6', road = '#eb6834')   # проверено на различимость при дальтонизме
shapes <- c(straight = 19, road = 15)                  # круг / квадрат — различие не только цветом
offset <- c(straight = 0.17, road = -0.17)
sample_titles <- c(all = 'Полная выборка', trimmed = 'Без 1% крайних значений с каждой стороны')

png('figures/08_road_vs_straight_coefficients.png', width = 2000, height = 1150, res = 170)
par(mfrow = c(1, 2), mar = c(4.2, 1, 3, 1), oma = c(0, 13, 2.2, 0), col.axis = '#52514e', fg = '#52514e')
for (smp in c('all', 'trimmed')) {
  p <- plot_data[plot_data$sample == smp, ]
  p$y <- rows$y[match(paste(p$center, p$group), paste(rows$center, rows$group))] + offset[p$distance]
  plot(NA, xlim = range(c(p$conf_low, p$conf_high, 0)), ylim = range(rows$y) + c(-0.6, 0.6),
       yaxt = 'n', ylab = '', xlab = 'Коэффициент при log(1 + км), млн руб. на работника',
       main = sample_titles[smp], cex.main = 1, font.main = 1)
  abline(v = 0, col = '#9A9993', lty = 2)
  segments(p$conf_low, p$y, p$conf_high, p$y, col = colors[p$distance], lwd = 2)
  points(p$estimate, p$y, pch = shapes[p$distance], col = colors[p$distance], cex = 0.95)
  if (smp == 'all') axis(2, at = rows$y, labels = rows$label, las = 1, tick = FALSE, outer = TRUE, cex.axis = 0.85)
}
mtext('Связь расстояния с разностью доходов и расходов на работника: микропредприятия, 95%-е интервалы',
      outer = TRUE, line = 0.8, cex = 1)
legend('bottomright', legend = c('по прямой', 'по дорогам'), col = colors, lwd = 2, pch = shapes,
       bty = 'n', cex = 0.85, inset = c(0, -0.01))
dev.off()


# ---- 9. Независимая проверка кластерных ошибок (если установлен sandwich) ----

if (requireNamespace('sandwich', quietly = TRUE)) {
  m <- lm(make_formula('X * (auto + small) + ln_settl_pop', 'road', 'krasnodar'), data = d)
  diff_se <- max(abs(sqrt(diag(cluster_vcov(m, d$settlement_id))) -
                     sqrt(diag(sandwich::vcovCL(m, cluster = d$settlement_id, type = 'HC1')))))
  cat('\nРасхождение кластерных ошибок с sandwich::vcovCL (HC1):', signif(diff_se, 3), '\n')
}
