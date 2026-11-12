# ============================================================
# Расстояния ПО ДОРОГАМ до центров края (Краснодар, Сочи, Новороссийск)
# ============================================================
#
# Запускать ПОСЛЕ prepare_data_git.R (нужен готовый датасет).
#
# Что делает скрипт:
#   1. Берёт итоговый датасет data/krasnodar_2023_regression.rds
#   2. Для каждого населённого пункта выборки считает расстояние по
#      автодорогам (км) и время в пути (мин) до трёх центров через сервис
#      маршрутизации OSRM (работает на данных OpenStreetMap)
#   3. Сохраняет ответ сервиса в data/road_distances.csv. При следующих
#      запусках расстояния берутся из этого файла: интернет не нужен,
#      результат не меняется от запуска к запуску
#   4. Добавляет в датасет столбцы road_* и ln_road_* и перезаписывает
#      data/krasnodar_2023_regression.csv и .rds
#
# Считаем по населённым пунктам, а не по фирмам: у всех фирм одного пункта
# одинаковые координаты, поэтому запросов нужно в разы меньше.
#
# Источник: OSRM, https://project-osrm.org; дорожный граф © участники
# OpenStreetMap, лицензия ODbL. Публичный демо-сервер — только для
# некоммерческого использования и не чаще 1 запроса в секунду.
#
# Пересчитать расстояния заново: удалите data/road_distances.csv.


# ---- 0. Подготовка ------------------------------------------

# Пакеты ставятся только если их ещё нет (jsonlite входит в tidyverse)
pkgs <- c("dplyr", "readr", "jsonlite")
new_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

library(dplyr)
library(readr)
library(jsonlite)

# Рабочая папка — корень проекта (откройте krasnodar_sme.Rproj).
# Если скрипт запущен из папки scripts/, поднимаемся на уровень выше.
if (!dir.exists("data") && dir.exists("../data")) setwd("..")
getwd()

data_file   <- "data/krasnodar_2023_regression.rds"
cache_file  <- "data/road_distances.csv"
osrm_server <- "https://router.project-osrm.org"
# запасной сервер, если основной недоступен:
# osrm_server <- "https://routing.openstreetmap.de/routed-car"
chunk_size  <- 40     # пунктов в одном запросе (сервер ограничивает размер таблицы)

options(HTTPUserAgent = "HSE econometrics student project (R)")

reg <- readRDS(data_file)


# ---- 1. Центры и населённые пункты --------------------------

# Центры — те же, что для прямых расстояний: берём столбцы dist_<центр>
# и точку, где расстояние до центра равно 0 (фирмы в самом центре).
center_slugs <- setdiff(sub("^dist_", "", grep("^dist_", names(reg), value = TRUE)),
                        "nearest")
centers <- data.frame(slug = center_slugs, lat = NA_real_, lon = NA_real_)
for (i in seq_len(nrow(centers))) {
  in_center <- reg[[paste0("dist_", centers$slug[i])]] == 0
  centers$lat[i] <- reg$lat[in_center][1]
  centers$lon[i] <- reg$lon[in_center][1]
}
if (anyNA(centers)) stop("Не удалось найти координаты центра: ",
                         paste(centers$slug[is.na(centers$lat)], collapse = ", "))
centers

# Уникальные населённые пункты выборки
settl <- reg %>%
  distinct(settlement_id, settlement, lat, lon)
nrow(settl)


# ---- 2. Запрос к OSRM (только для пунктов, которых нет в кэше) ----

# Один запрос к сервису table: от пачки пунктов (src) до центров (dst).
# Возвращает матрицы расстояний (м) и времени (с): строки — пункты,
# столбцы — центры. NA — маршрут не найден.
osrm_table <- function(src, dst) {
  coords <- paste(sprintf("%.6f,%.6f", c(dst$lon, src$lon), c(dst$lat, src$lat)),
                  collapse = ";")                       # OSRM ждёт "долгота,широта"
  n_dst <- nrow(dst)
  query <- paste0(osrm_server, "/table/v1/driving/", coords,
                  "?sources=", paste(n_dst + seq_len(nrow(src)) - 1, collapse = ";"),
                  "&destinations=", paste(seq_len(n_dst) - 1, collapse = ";"),
                  "&annotations=duration,distance")
  con <- url(query)
  txt <- tryCatch(readLines(con, warn = FALSE),
                  error = function(e) {
                    close(con)
                    stop("Сервис OSRM (", osrm_server, ") недоступен. Проверьте интернет ",
                         "или раскомментируйте запасной сервер в начале скрипта.", call. = FALSE)
                  })
  close(con)
  res <- fromJSON(paste(txt, collapse = ""))
  if (res$code != "Ok") {
    stop("OSRM вернул ошибку ", res$code, ". Если TooBig — уменьшите chunk_size.",
         call. = FALSE)
  }
  snap <- res$sources$distance                          # на сколько метров точку сдвинули до дороги
  if (is.null(snap)) snap <- rep(NA_real_, nrow(src))
  list(distance = matrix(res$distances, nrow = nrow(src)),
       duration = matrix(res$durations, nrow = nrow(src)),
       snap     = snap)
}

road <- if (file.exists(cache_file)) {
  read_csv(cache_file, col_types = cols(settlement_id = col_character(),
                                        settlement    = col_character(),
                                        calc_date     = col_character()))
} else {
  NULL
}
todo <- if (is.null(road)) settl else filter(settl, !settlement_id %in% road$settlement_id)
cat("Пунктов в кэше:", ifelse(is.null(road), 0, nrow(road)),
    "| нужно запросить:", nrow(todo), "\n")

if (nrow(todo) > 0) {
  chunks <- split(todo, ceiling(seq_len(nrow(todo)) / chunk_size))
  parts <- list()
  for (k in seq_along(chunks)) {
    ch  <- chunks[[k]]
    tab <- osrm_table(ch, centers)
    part <- data.frame(settlement_id = ch$settlement_id, settlement = ch$settlement,
                       lat = ch$lat, lon = ch$lon, snap_m = round(tab$snap))
    for (j in seq_len(nrow(centers))) {
      part[[paste0("road_", centers$slug[j])]] <- round(tab$distance[, j] / 1000, 2)  # км
      part[[paste0("time_", centers$slug[j])]] <- round(tab$duration[, j] / 60, 1)    # мин
    }
    parts[[k]] <- part
    cat("Запрос", k, "из", length(chunks), "— готово\n")
    Sys.sleep(1.5)                                      # правило сервиса: ≤ 1 запроса в секунду
  }
  road <- bind_rows(road, bind_rows(parts) %>% mutate(calc_date = as.character(Sys.Date())))
  write_csv(road, cache_file)
  cat("Сохранено в", cache_file, "\n")
}


# ---- 3. Проверки --------------------------------------------

road_cols <- paste0("road_", centers$slug)

# (а) Пары без маршрута
sum(is.na(road[, road_cols]))

# (б) Точки, которые пришлось далеко сдвигать до ближайшей дороги (> 2 км):
# координаты пункта могут стоять в поле — тогда расстояние чуть неточное
road %>% filter(snap_m > 2000) %>% select(settlement, snap_m)

# (в) Дорога не может быть заметно короче прямой. Коэффициент извилистости
# = по дорогам / по прямой (для пунктов дальше 5 км от центра)
check <- reg %>%
  distinct(settlement_id, dist_krasnodar) %>%
  left_join(select(road, settlement_id, settlement, road_krasnodar), by = "settlement_id") %>%
  filter(dist_krasnodar > 5) %>%
  mutate(ratio = road_krasnodar / dist_krasnodar)

summary(check$ratio)
check %>% filter(ratio < 0.95)                                     # подозрительные, должно быть пусто
check %>% arrange(desc(ratio)) %>% head(10)                         # самые извилистые маршруты
cor(log(1 + check$dist_krasnodar), log(1 + check$road_krasnodar))   # согласованность с прямым


# ---- 4. Добавляем столбцы в датасет -------------------------

reg <- reg %>%
  select(-starts_with("road_"), -starts_with("ln_road_")) %>%      # если скрипт уже запускали
  left_join(select(road, settlement_id, all_of(road_cols)), by = "settlement_id")

for (s in centers$slug) {
  reg[[paste0("ln_road_", s)]] <- log(1 + reg[[paste0("road_", s)]])
}
# ближайший центр по дорогам (если до какого-то центра маршрута нет — по остальным)
reg$road_nearest    <- apply(as.matrix(reg[, road_cols]), 1, function(x) {
  if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
})
reg$ln_road_nearest <- log(1 + reg$road_nearest)

# Ставим новые столбцы сразу после прямых расстояний
reg <- reg %>%
  relocate(starts_with("road_"), starts_with("ln_road_"), .after = ln_dist_nearest)

write_csv(reg, "data/krasnodar_2023_regression.csv")
saveRDS(reg, "data/krasnodar_2023_regression.rds")

# Краткая проверка результата
dim(reg)
summary(select(reg, dist_krasnodar, road_krasnodar, dist_nearest, road_nearest))
