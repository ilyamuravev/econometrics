# Разведочный анализ (только то, что нужно по заданию):
#   1) аномальные значения  2) матрица корреляций и проверка связей между признаками
#   3) кластеры            4) карты
# Признаки: доходы (доход на одного работника), количество сотрудников, размер населённого пункта.
# Краснодарский край, 2023, микро- и малые предприятия: продуктовая розница и автосервис.
#
# Запуск в RStudio: Session -> Set Working Directory -> To Source File Location, затем Source.
# Файл krasnodar_all.parquet должен лежать рядом (создаётся автоматически из ../panel.csv, если его нет).
#
# Пакеты (один раз):
# options(timeout = 600)
# install.packages(c("duckdb", "DBI", "dplyr", "ggplot2", "corrplot", "leaflet",
#                    "htmlwidgets", "cluster", "tidyr", "gridExtra", "viridisLite"), type = "binary")

library(DBI)
library(duckdb)
library(dplyr)
library(ggplot2)
library(corrplot)
library(leaflet)
library(htmlwidgets)
library(cluster)

PANEL    <- "../panel.csv"
PARQUET  <- "krasnodar_all.parquet"
USE_DIST <- FALSE   # TRUE — добавить расстояние до Краснодара в корреляции и кластеры
KRD      <- c(lat = 45.0355, lon = 38.9753)
set.seed(0)

# ---------- 1. данные ----------
con <- dbConnect(duckdb())
if (!file.exists(PARQUET)) {
  dbExecute(con, sprintf("
    COPY (SELECT * FROM read_csv('%s', delim=';', header=true, all_varchar=true)
          WHERE region_code='23') TO '%s' (FORMAT parquet)", PANEL, PARQUET))
}
df <- dbGetQuery(con, sprintf("
  SELECT tin, category, activity_code_main, settlement, oktmo,
         try_cast(lat AS double) lat, try_cast(lon AS double) lon,
         try_cast(revenue AS double) revenue, try_cast(employees_count AS double) employees
  FROM '%s'
  WHERE year='2023' AND category IN ('1','2')
    AND (activity_code_main LIKE '45.20%%' OR activity_code_main LIKE '47.11%%' OR activity_code_main LIKE '47.2%%')",
  PARQUET))
# населения в данных нет: размер пункта = число зарегистрированных там фирм (все ОКВЭД)
size <- dbGetQuery(con, sprintf(
  "SELECT oktmo, count(distinct tin) AS n_firms_settl FROM '%s' WHERE year='2023' GROUP BY 1", PARQUET))
dbDisconnect(con, shutdown = TRUE)

haversine <- function(lat, lon, lat0 = KRD["lat"], lon0 = KRD["lon"]) {
  p <- pi / 180
  a <- sin((lat - lat0) * p / 2)^2 + cos(lat * p) * cos(lat0 * p) * sin((lon - lon0) * p / 2)^2
  12742 * asin(sqrt(a))
}

df <- df |>
  left_join(size, by = "oktmo") |>
  mutate(sector = ifelse(startsWith(activity_code_main, "45.20"), "Автосервис", "Продуктовая розница"),
         dist_km = haversine(lat, lon)) |>
  filter(!is.na(dist_km), dist_km <= 400) |>
  distinct(tin, .keep_all = TRUE)               # у части фирм несколько адресов — одна строка на ИНН

d <- df |>
  filter(!is.na(revenue), revenue > 0, !is.na(employees), employees > 0) |>
  mutate(rev_per_emp = revenue / employees,
         ln_rpe = log(rev_per_emp), ln_emp = log(employees),
         ln_firms = log(n_firms_settl), ln_dist = log1p(dist_km))

feat <- c(ln_rpe = "ln доход на работника", ln_emp = "ln сотрудники",
          ln_firms = "ln размер пункта (число фирм)")
if (USE_DIST) feat <- c(feat, ln_dist = "ln расстояние до Краснодара")

cat("фирм в выборке:", nrow(d), "\n"); print(count(d, sector))
print(summary(d[, c("revenue", "employees", "rev_per_emp", "n_firms_settl")]))

# ---------- 2. аномальные значения ----------
tukey <- function(s, k) { q <- quantile(s, c(.25, .75)); i <- q[2] - q[1]; s < q[1] - k * i | s > q[2] + k * i }

d <- d |>
  group_by(sector) |>
  mutate(anomaly = tukey(ln_rpe, 1.5),          # граница Тьюки по доходу на работника, внутри сектора
         extreme = tukey(ln_rpe, 3)) |>         # более строгая граница
  ungroup()

cat("\nАномалии по доходу на работника (Тьюки 1,5*IQR; extreme — 3*IQR):\n")
print(d |> group_by(sector) |> summarise(n = n(), anomaly = sum(anomaly), extreme = sum(extreme)))
cat("\nСамые высокие значения:\n")
print(d |> arrange(desc(rev_per_emp)) |> select(tin, settlement, employees, revenue, rev_per_emp) |> head(8))
cat("\nСамые низкие значения:\n")
print(d |> arrange(rev_per_emp) |> select(tin, settlement, employees, revenue, rev_per_emp) |> head(8))
write.csv(d |> filter(anomaly), "anomalies_me.csv", row.names = FALSE, fileEncoding = "UTF-8")

long <- d |>
  transmute(sector, `Доход на работника, руб` = rev_per_emp, `Сотрудники` = employees,
            `Размер пункта (число фирм)` = n_firms_settl) |>
  tidyr::pivot_longer(-sector, names_to = "var", values_to = "value")
p1 <- ggplot(long, aes(sector, value)) + geom_boxplot(outlier.size = .6) + scale_y_log10() +
  facet_wrap(~var, scales = "free_y") + labs(x = NULL, y = NULL, title = "Распределения и выбросы (лог-шкала)") +
  theme_minimal(base_size = 11)
ggsave("01_outliers_me.png", p1, width = 12, height = 4.5, dpi = 130); print(p1)

# ---------- 3. корреляции ----------
clean <- d |> filter(!anomaly)
groups <- c(list(`Все` = clean), split(clean, clean$sector))

png("02_correlation_matrices_me.png", width = 2000, height = 1300, res = 130)
par(mfrow = c(2, 3))
for (m in c("pearson", "spearman")) for (nm in names(groups)) {
  g <- groups[[nm]]
  cm <- cor(g[, names(feat)], method = m); dimnames(cm) <- list(feat, feat)
  corrplot(cm, method = "color", addCoef.col = "black", number.cex = .9, tl.col = "black", tl.cex = .75,
           col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
           title = sprintf("%s — %s (n=%d)", m, nm, nrow(g)), mar = c(0, 0, 3, 0))
}
dev.off()

cat("\nПроверка корреляций (Спирмен, значимость по парам):\n")
for (nm in names(groups)) {
  g <- groups[[nm]]; cat("\n--", nm, "(n =", nrow(g), ")\n")
  tab <- do.call(rbind, lapply(combn(names(feat), 2, simplify = FALSE), function(p) {
    ct <- cor.test(g[[p[1]]], g[[p[2]]], method = "spearman", exact = FALSE)
    data.frame(признак_1 = feat[[p[1]]], признак_2 = feat[[p[2]]], rho = round(unname(ct$estimate), 3),
               p = signif(ct$p.value, 3), значима_5 = ct$p.value < .05)
  }))
  print(tab, row.names = FALSE)
}

# ---------- 4. кластеры ----------
X <- scale(clean[, names(feat)])
dm <- dist(X)
sil <- sapply(2:7, function(k) mean(silhouette(kmeans(X, centers = k, nstart = 25)$cluster, dm)[, "sil_width"]))
names(sil) <- 2:7
cat("\nSilhouette по k:\n"); print(round(sil, 3))
K <- as.integer(names(which.max(sil))); cat("выбрано k =", K, "\n")
clean$cluster <- factor(kmeans(X, centers = K, nstart = 25)$cluster)

prof <- clean |> group_by(cluster) |>
  summarise(n = n(), доход_на_работника = median(rev_per_emp), сотрудники = median(employees),
            размер_пункта = median(n_firms_settl), доля_розницы = round(mean(sector == "Продуктовая розница"), 2))
print(as.data.frame(prof))
write.csv(prof, "cluster_profile_me.csv", row.names = FALSE, fileEncoding = "UTF-8")

pc1 <- ggplot(clean, aes(ln_emp, ln_rpe, colour = cluster)) + geom_jitter(width = .08, size = 1) +
  labs(x = "ln сотрудники", y = "ln доход на работника", title = "Кластеры (k-means)") + theme_minimal()
pc2 <- ggplot(clean, aes(ln_firms, ln_rpe, colour = cluster)) + geom_point(size = 1) +
  labs(x = "ln размер пункта (число фирм)", y = "ln доход на работника", title = "Кластеры и размер пункта") + theme_minimal()
g4 <- gridExtra::arrangeGrob(pc1, pc2, ncol = 2)
ggsave("03_clusters_me.png", g4, width = 13, height = 5.5, dpi = 130)
grid::grid.newpage(); grid::grid.draw(g4)

# ---------- 5. карты ----------
agg <- d |> group_by(oktmo, settlement, sector) |>
  summarise(lat = first(lat), lon = first(lon), n = n(), rpe = median(rev_per_emp), .groups = "drop")
p5 <- ggplot(agg, aes(lon, lat, colour = log10(rpe), size = n)) + geom_point(alpha = .8) +
  scale_colour_viridis_c(name = "log10 медианного\nдохода на работника") +
  scale_size_continuous(range = c(1.5, 8), name = "фирм") +
  facet_wrap(~sector) + coord_quickmap() + theme_minimal(base_size = 11)
ggsave("04_map_income_me.png", p5, width = 15, height = 6.5, dpi = 130); print(p5)

pal <- colorFactor("Set1", clean$cluster)
m <- leaflet() |>
  addProviderTiles(providers$CartoDB.Positron) |>
  setView(KRD[["lon"]], KRD[["lat"]], zoom = 7) |>
  addCircleMarkers(data = clean, lng = ~lon, lat = ~lat, radius = 3, color = ~pal(cluster), weight = .5, fillOpacity = .8,
                   label = ~sprintf("%s (%s): доход/раб. %.2f млн; сотрудники %.0f; кластер %s",
                                    settlement, sector, rev_per_emp / 1e6, employees, cluster)) |>
  addCircleMarkers(data = filter(d, anomaly), lng = ~lon, lat = ~lat, radius = 5, color = "black", fillOpacity = 1,
                   label = ~sprintf("АНОМАЛИЯ %s: %.2f млн/раб., сотрудники %.0f", settlement, rev_per_emp / 1e6, employees))
saveWidget(m, file = normalizePath("05_map_interactive_me.html", mustWork = FALSE), selfcontained = FALSE)

cat("\nготово, файлы с суффиксом _me в", getwd(), "\n")
