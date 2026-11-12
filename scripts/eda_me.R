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
#                    "htmlwidgets", "cluster", "gridExtra", "viridisLite", "sf", "ggspatial"), type = "binary")

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

# гистограммы в обычных единицах; верхний 1% по доходу и сотрудникам не показан (там самые длинные хвосты)
graphics.off()   # чистое окно графиков, чтобы первым не появлялся пустой лист
q99 <- function(x) unname(quantile(x, .99, na.rm = TRUE))
lab <- scales::label_number(big.mark = " ", decimal.mark = ",")
cut_inc <- q99(d$rev_per_emp); cut_emp <- q99(d$employees)

h1 <- ggplot(filter(d, rev_per_emp <= cut_inc), aes(rev_per_emp / 1e6, fill = sector)) +
  geom_histogram(bins = 40, alpha = .6, position = "identity") +
  scale_x_continuous(labels = lab) + scale_y_continuous(labels = lab) +
  labs(x = "млн руб. на работника в год", y = "число фирм", fill = NULL, title = "Доход на работника",
       subtitle = sprintf("без верхнего 1%% (больше %.1f млн руб.)", cut_inc / 1e6)) +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
h2 <- ggplot(filter(d, employees <= cut_emp), aes(employees, fill = sector)) +
  geom_histogram(breaks = seq(0.5, ceiling(cut_emp) + 0.5, by = 1), alpha = .6, position = "identity") +
  scale_x_continuous(labels = lab) + scale_y_continuous(labels = lab) +
  labs(x = "сотрудников в фирме", y = "число фирм", title = "Число сотрудников",
       subtitle = sprintf("без верхнего 1%% (больше %.0f чел.)", cut_emp)) +
  theme_minimal(base_size = 11) + theme(legend.position = "none")
h3 <- ggplot(mutate(d, size_class = cut(n_firms_settl, breaks = c(0, 500, 2000, 10000, 50000, Inf),
                                        labels = c("до 500", "500–2 000", "2 000–10 000", "10 000–50 000", "больше 50 000"))),
             aes(size_class, fill = sector)) +
  geom_bar(position = "dodge", alpha = .85) +
  scale_y_continuous(labels = lab) +
  labs(x = "фирм, зарегистрированных в пункте", y = "число фирм в выборке", title = "Размер населённого пункта") +
  theme_minimal(base_size = 11) + theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
p1 <- gridExtra::arrangeGrob(h1, h2, h3, ncol = 3)
ggsave("01_distributions_me.png", p1, width = 14, height = 5, dpi = 130)
grid::grid.draw(p1)

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

# кластеры считались на логарифмах, но на графиках показаны обычные единицы
cl <- clean |> filter(employees <= cut_emp, rev_per_emp <= cut_inc)
cl_note <- "без верхнего 1% по доходу и по сотрудникам"
pc1 <- ggplot(cl, aes(employees, rev_per_emp / 1e6, colour = cluster)) +
  geom_jitter(width = .25, height = 0, size = 1, alpha = .7) +
  scale_x_continuous(labels = lab) + scale_y_continuous(labels = lab) +
  labs(x = "сотрудников в фирме", y = "доход на работника, млн руб.", title = "Кластеры (k-means)", subtitle = cl_note) +
  theme_minimal()
pc2 <- ggplot(cl, aes(n_firms_settl / 1000, rev_per_emp / 1e6, colour = cluster)) +
  geom_point(size = 1, alpha = .7) +
  scale_x_continuous(labels = lab) + scale_y_continuous(labels = lab) +
  labs(x = "фирм в населённом пункте, тыс.", y = "доход на работника, млн руб.",
       title = "Кластеры и размер пункта", subtitle = cl_note) +
  theme_minimal()
g4 <- gridExtra::arrangeGrob(pc1, pc2, ncol = 2)
ggsave("03_clusters_me.png", g4, width = 13, height = 5.5, dpi = 130)
grid::grid.newpage(); grid::grid.draw(g4)

# ---------- 5. карты ----------
agg <- d |> group_by(oktmo, settlement, sector) |>
  summarise(lat = first(lat), lon = first(lon), n = n(), rpe = median(rev_per_emp), .groups = "drop")
# карта края: точки пунктов поверх подложки (нужен интернет для загрузки плиток карты;
# если загрузить не удалось, рисуется карта без подложки)
agg_sf <- sf::st_as_sf(agg, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
# граница Краснодарского края (GADM) лежит рядом со скриптом; без неё карта рисуется без контура
border_file <- Filter(file.exists, c("krasnodar_border.geojson", "data/krasnodar_border.geojson", "../data/krasnodar_border.geojson"))
border <- if (length(border_file)) sf::st_read(border_file[1], quiet = TRUE) else NULL
make_map <- function(tiles) {
  ggplot(agg_sf) +
    (if (tiles) ggspatial::annotation_map_tile(type = "cartolight", zoom = 8, progress = "none")) +
    (if (!is.null(border)) geom_sf(data = border, fill = if (tiles) NA else "#eef1f4", colour = "#3d4651", linewidth = .6)) +
    geom_sf(aes(colour = rpe / 1e6, size = n), alpha = .85) +
    scale_colour_viridis_c(name = "Медианный доход\nна работника, млн руб.", trans = "log10",
                           breaks = c(.1, .3, 1, 3, 10, 30), labels = c("0,1", "0,3", "1", "3", "10", "30")) +
    scale_size_continuous(range = c(1.5, 8), name = "фирм в выборке") +
    coord_sf(xlim = c(36.5, 42.0), ylim = c(43.3, 47.0), expand = FALSE) +
    facet_wrap(~sector) + theme_minimal(base_size = 11)
}
p5 <- make_map(TRUE)
ok <- tryCatch({ ggsave("04_map_income_me.png", p5, width = 15, height = 6.5, dpi = 130); print(p5); TRUE },
               error = function(e) { message("Подложка карты недоступна (", conditionMessage(e), "), рисую без неё"); FALSE })
if (!ok) { p5 <- make_map(FALSE); ggsave("04_map_income_me.png", p5, width = 15, height = 6.5, dpi = 130); print(p5) }

pal <- colorFactor("Set1", clean$cluster)
add_border <- function(map) {
  if (is.null(border)) map else addPolygons(map, data = border, fill = FALSE, color = "#333333", weight = 2)
}
m <- leaflet() |>
  addProviderTiles(providers$CartoDB.Positron) |>
  setView(KRD[["lon"]], KRD[["lat"]], zoom = 7) |>
  add_border() |>
  addCircleMarkers(data = clean, lng = ~lon, lat = ~lat, radius = 3, color = ~pal(cluster), weight = .5, fillOpacity = .8,
                   label = ~sprintf("%s (%s): доход/раб. %.2f млн; сотрудники %.0f; кластер %s",
                                    settlement, sector, rev_per_emp / 1e6, employees, cluster)) |>
  addCircleMarkers(data = filter(d, anomaly), lng = ~lon, lat = ~lat, radius = 5, color = "black", fillOpacity = 1,
                   label = ~sprintf("АНОМАЛИЯ %s: %.2f млн/раб., сотрудники %.0f", settlement, rev_per_emp / 1e6, employees))
saveWidget(m, file = normalizePath("05_map_interactive_me.html", mustWork = FALSE), selfcontained = FALSE)

cat("\nготово, файлы с суффиксом _me в", getwd(), "\n")
