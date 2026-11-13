# Наглядное сравнение моделей с расстоянием по прямой и по дорогам.
# Читает результаты scripts/regression_profit_road.R, модели заново не оценивает.
#
# Запуск из корня проекта (нужен только базовый R):
#   source("scripts/plot_road_vs_straight.R")
#
# Результаты:
#   figures/09_road_vs_straight_scatter.png — все 42 оценки: по прямой против по дорогам
#   figures/10_road_vs_straight_table.png   — таблица главных моделей для слайда
#   results/profit_road_comparison_slide.csv — та же таблица текстом (например, для knitr::kable)

cmp <- read.csv('results/profit_road_comparison.csv', stringsAsFactors = FALSE)
stopifnot(nrow(cmp) > 0)

ink <- '#0b0b0b'; ink2 <- '#52514e'; grid_col <- '#d6d5d0'
num <- function(x, digits = 3) {                       # 0,094 и −0,038 — запятая и типографский минус
  x[abs(x) < 0.5 * 10^-digits] <- 0                     # без «−0,000»
  s <- formatC(x, format = 'f', digits = digits, decimal.mark = ',')
  sub('^-', '−', s)
}
excludes_zero <- function(low, high) low > 0 | high < 0

# Сводные числа для подписей
shift <- cmp$road_estimate - cmp$straight_estimate
half_width <- (cmp$straight_conf_high - cmp$straight_conf_low) / 2
changed <- excludes_zero(cmp$straight_conf_low, cmp$straight_conf_high) !=
           excludes_zero(cmp$road_conf_low, cmp$road_conf_high)
cat(sprintf('Оценок: %d | сдвиг в среднем %.4f, максимум %.4f | медиана сдвига %.0f%% полуширины интервала | выводов изменилось: %d\n',
            nrow(cmp), mean(abs(shift)), max(abs(shift)), 100 * median(abs(shift) / half_width), sum(changed)))


# ---- 1. Диаграмма рассеяния: по прямой против по дорогам ----

samples <- c(all = 'полная выборка', trimmed = 'без 1% хвостов', no_top_cook = 'без влиятельной фирмы')
colors  <- c(all = '#2a78d6', trimmed = '#eb6834', no_top_cook = '#1baf7a')  # проверены на различимость
shapes  <- c(all = 19, trimmed = 17, no_top_cook = 15)                       # форма дублирует цвет
lim <- range(c(cmp$straight_estimate, cmp$road_estimate)) + c(-0.02, 0.02)

png('figures/09_road_vs_straight_scatter.png', width = 1700, height = 1500, res = 190)
par(mar = c(4.5, 4.8, 4.6, 1.2), col.axis = ink2, fg = ink2, las = 1)
plot(NA, xlim = lim, ylim = lim, asp = 1,
     xlab = 'Оценка с расстоянием по прямой', ylab = 'Оценка с расстоянием по дорогам')
abline(h = 0, v = 0, col = grid_col)
abline(0, 1, col = ink2, lty = 2)                       # линия равенства: оценки совпадают
for (s in names(samples)) {
  p <- cmp[cmp$sample == s, ]
  points(p$straight_estimate, p$road_estimate, pch = shapes[s], col = colors[s], cex = 1.15)
}
# Подписи к двум наибольшим сдвигам (текст — нейтральным цветом)
top <- order(-abs(shift))[1:2]
text(cmp$straight_estimate[top], cmp$road_estimate[top], pos = 4, offset = 0.6, cex = 0.75, col = ink2,
     labels = c('Сочи: автосервис,\nполная выборка', 'Сочи: автосервис − розница,\nполная выборка'))
text(lim[1], lim[2], adj = c(0, 1), cex = 0.8, col = ink2,
     labels = sprintf('Пунктир — линия равенства.\nСдвиг оценок: в среднем %s, максимум %s млн руб.\nИзменились выводы о значимости: %d из %d',
                      num(mean(abs(shift))), num(max(abs(shift))), sum(changed), nrow(cmp)))
legend('bottomright', legend = samples, col = colors, pch = shapes, bty = 'n', cex = 0.85, text.col = ink)
mtext('Коэффициенты при log(1 + км) почти не меняются при переходе к дорогам',
      side = 3, line = 2.4, cex = 1.05, col = ink, font = 2)
mtext(sprintf('Ключевые оценки (%d шт.): модели с контролями и отраслевые коэффициенты по 4 центрам; млн руб. на работника', nrow(cmp)),
      side = 3, line = 1, cex = 0.8, col = ink2)
dev.off()


# ---- 2. Таблица главных моделей -----------------------------

rows <- data.frame(
  model = c('controls', 'trimmed', 'nearest',
            'krasnodar_all_interactions', 'krasnodar_trimmed_interactions', 'krasnodar_no_top_cook_interactions',
            'sochi_all_interactions', 'sochi_no_top_cook_interactions',
            'sochi_all_interactions', 'sochi_trimmed_interactions', 'sochi_no_top_cook_interactions'),
  group = c('all_firms', 'all_firms', 'all_firms', 'difference', 'difference', 'difference',
            'difference', 'difference', 'retail', 'retail', 'retail'),
  label = c('Краснодар, с контролями', '', 'Ближайший центр, с контролями',
            'Краснодар: автосервис − розница', '', '',
            'Сочи: автосервис − розница', '',
            'Сочи: розница', '', ''),
  stringsAsFactors = FALSE)
tab <- merge(rows, cmp, by = c('model', 'group'), sort = FALSE)
tab <- tab[match(paste(rows$model, rows$group), paste(tab$model, tab$group)), ]
stopifnot(nrow(tab) == nrow(rows), !anyNA(tab$straight_estimate))

ci <- function(e, l, h) paste0(num(e), '  [', num(l), '; ', num(h), ']')
yes_no <- function(l, h) ifelse(excludes_zero(l, h), 'да', 'нет')
slide <- data.frame(
  model    = tab$label,
  sample   = unname(samples[tab$sample]),
  straight = ci(tab$straight_estimate, tab$straight_conf_low, tab$straight_conf_high),
  road     = ci(tab$road_estimate, tab$road_conf_low, tab$road_conf_high),
  shift    = num(tab$road_estimate - tab$straight_estimate),
  excludes_zero = paste(yes_no(tab$straight_conf_low, tab$straight_conf_high), '→',
                        yes_no(tab$road_conf_low, tab$road_conf_high)),
  stringsAsFactors = FALSE)
write.csv(setNames(slide, c('Модель', 'Выборка', 'По прямой [95% ДИ]', 'По дорогам [95% ДИ]',
                            'Сдвиг', 'Интервал исключает 0: прямое → дороги')),
          'results/profit_road_comparison_slide.csv', row.names = FALSE)

n_rows <- nrow(slide)
x <- c(0.010, 0.300, 0.455, 0.655, 0.835, 0.990)        # левые края колонок; последняя — правый край
header <- c('Модель', 'Выборка', 'По прямой [95% ДИ]', 'По дорогам [95% ДИ]', 'Сдвиг', 'Интервал\nисключает 0')
y_top <- n_rows + 1

png('figures/10_road_vs_straight_table.png', width = 2300, height = 1150, res = 190)
par(mar = c(0, 0, 0, 0))
plot.new(); plot.window(xlim = c(0, 1), ylim = c(-1.6, n_rows + 3))
text(0.010, n_rows + 2.5, 'Расстояние по прямой и по дорогам: главные модели',
     adj = c(0, 0.5), font = 2, cex = 1.15, col = ink)
text(0.010, n_rows + 1.85, 'Коэффициент при log(1 + км), млн руб. на работника; отраслевые оценки — для микропредприятий',
     adj = c(0, 0.5), cex = 0.85, col = ink2)
for (i in seq_len(n_rows)) if (i %% 2 == 0) {
  rect(0, y_top - i - 0.5, 1, y_top - i + 0.5, col = '#f2f1ee', border = NA)
}
segments(0, y_top - 0.5, 1, y_top - 0.5, col = ink2)
for (i in which(slide$model != '')[-1]) segments(0, y_top - i + 0.5, 1, y_top - i + 0.5, col = grid_col)
text(x[1:4], y_top, header[1:4], adj = c(0, 0.5), font = 2, cex = 0.82, col = ink)
text(x[5] + 0.02, y_top, header[5], adj = c(0, 0.5), font = 2, cex = 0.82, col = ink)
text(x[6], y_top, header[6], adj = c(1, 0.5), font = 2, cex = 0.82, col = ink)
for (i in seq_len(n_rows)) {
  y <- y_top - i
  text(x[1], y, slide$model[i], adj = c(0, 0.5), cex = 0.82, col = ink)
  text(x[2], y, slide$sample[i], adj = c(0, 0.5), cex = 0.82, col = ink2)
  text(x[3], y, slide$straight[i], adj = c(0, 0.5), cex = 0.82, col = ink, family = 'mono')
  text(x[4], y, slide$road[i], adj = c(0, 0.5), cex = 0.82, col = ink, family = 'mono')
  text(x[5] + 0.02, y, slide$shift[i], adj = c(0, 0.5), cex = 0.82, col = ink, family = 'mono')
  text(x[6], y, slide$excludes_zero[i], adj = c(1, 0.5), cex = 0.82, col = ink)
}
segments(0, 0.5, 1, 0.5, col = ink2)
text(0.010, -0.35, adj = c(0, 0.5), cex = 0.78, col = ink2,
     labels = sprintf('Во всех %d ключевых оценках сдвиг в среднем %s, максимум %s; вывод о том, исключает ли интервал ноль, не изменился ни разу.',
                      nrow(cmp), num(mean(abs(shift))), num(max(abs(shift)))))
text(0.010, -1.05, adj = c(0, 0.5), cex = 0.78, col = ink2,
     labels = 'Выборки: полная — 832 фирмы; без 1% хвостов — 814; без влиятельной фирмы — 831 (без автосервиса в Сочи с 55 млн руб. на работника).')
dev.off()

print(slide, row.names = FALSE)
