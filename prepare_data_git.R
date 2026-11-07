# Проект по эконометрике: подготовка данных
#
# Отличия от prepare_data.R:
#   - нет setwd(): откройте krasnodar_sme.Rproj, папка проекта станет
#     рабочей директорией автоматически, у каждого участника своя
#   - пакеты ставятся, только если их ещё нет
#   - вырезка края хранится сжатой: data/panel_kda_2023.csv.gz (~12 МБ
#     вместо ~90 МБ), она лежит в репозитории, поэтому panel.csv
#     (16,7 ГБ) для запуска НЕ нужен
# «Удалённость от крупных центров и доходы малого бизнеса»
# Краснодарский край, 2023 год
#
# Что делает скрипт:
#   1. Достаёт из panel.csv (реестр МСП) строки Краснодарского края за 2023 год
#   2. Чистит данные: типы переменных, коды ОКТМО, повторы ИНН внутри года
#   3. Находит топ-3 города края по численности населения и по экономике
#   4. Отбирает выборку для регрессии по коду ОКВЭД (основной вид деятельности):
#        47.11 — продуктовая розница (магазины, где преобладают продукты питания)
#        45.20 — техническое обслуживание и ремонт автомобилей (автосервис)
#   5. Создаёт переменные расстояния до трёх центров
#   6. Мерджит численность населения населённого пункта (справочник ИНИД №160)
#   7. Сохраняет готовый датасет и таблицы для слайдов
#
# Какие файлы должны лежать в рабочей папке:
#   panel.csv        — «Реестр субъектов МСП за 2016–2024 гг.» (ИНИД, набор №205)
#                      (нужен, только если нет data/panel_kda_2023.csv.gz)
#   settlements_12032021/data.csv — «Населенные пункты России: численность
#                      населения и географические координаты» (ИНИД, набор №160)
#
# Результаты пишутся в папку data/:
#   panel_kda_2023.csv.gz        — строки края за 2023 год (промежуточный файл, сжатый)
#   top3_cities.csv              — рейтинг городов по населению и доходам МСП
#   cleaning_log.csv             — сколько наблюдений осталось после каждого шага
#   krasnodar_2023_regression.*  — итоговый датасет для регрессии


# 0. Подготовка

# Пакеты ставятся только если их ещё нет
pkgs <- c("dplyr", "readr", "stringr")
new_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new_pkgs) > 0) install.packages(new_pkgs)

library(dplyr)
library(readr)
library(stringr)

# Рабочая директория: откройте krasnodar_sme.Rproj — RStudio сам сделает
# папку проекта рабочей. setwd() не нужен, пути ниже относительные.
getwd()

YEAR         <- 2023                     # год наблюдений
REGION_ISO   <- "RU-KDA"                 # Краснодарский край в panel.csv
REGION_NAME  <- "Краснодарский"          # как край называется в справочнике №160

panel_file   <- "panel.csv"
settl_file   <- "settlements_12032021/data.csv"   # справочник населённых пунктов (ИНИД №160)
dir.create("data", showWarnings = FALSE)
extract_file <- "data/panel_kda_2023.csv.gz"   # сжатая вырезка (лежит в репозитории)
old_extract  <- "data/panel_kda_2023.csv"      # несжатая вырезка от prepare_data.R


# Вспомогательные функции

# В panel.csv ОКТМО записан как число, поэтому ведущий ноль потерян:
# у Краснодара должно быть 03701000001, а в файле "3701000001".
# Функция возвращает ноль (11 знаков у населённого пункта, 8 у муниципалитета).
fix_oktmo <- function(x) {
  x <- str_trim(as.character(x))
  case_when(
    is.na(x) | x == "" ~ NA_character_,
    nchar(x) == 10     ~ str_pad(x, 11, pad = "0"),
    nchar(x) == 7      ~ str_pad(x, 8,  pad = "0"),
    TRUE               ~ x
  )
}

# Расстояние по прямой между двумя точками на сфере (формула гаверсинусов), км
haversine_km <- function(lat1, lon1, lat2, lon2) {
  to_rad <- pi / 180
  dlat <- (lat2 - lat1) * to_rad
  dlon <- (lon2 - lon1) * to_rad
  a <- sin(dlat / 2)^2 + cos(lat1 * to_rad) * cos(lat2 * to_rad) * sin(dlon / 2)^2
  2 * 6371 * asin(pmin(1, sqrt(a)))
}

# Ключ для сравнения названий: нижний регистр, ё -> е, без лишних пробелов
name_key <- function(x) {
  str_squish(str_replace_all(str_to_lower(x), "ё", "е"))
}


# 1. Извлечение Краснодарского края из panel.csv

# Файл слишком большой, чтобы загрузить его в память целиком.
# Читаем его кусками по 500 тыс. строк и оставляем только строки,
# где есть "RU-KDA" и год 2023. Это грубый текстовый фильтр, точная
# проверка по столбцам будет в шаге 2..
#
# Порядок такой:
#   - есть data/panel_kda_2023.csv.gz -> шаг пропускается (обычный случай
#     после клонирования репозитория);
#   - есть несжатая data/panel_kda_2023.csv (от prepare_data.R) -> просто
#     сжимаем её, 16 ГБ заново не читаем;
#   - иначе читаем panel.csv (5–15 минут).

if (!file.exists(extract_file) && file.exists(old_extract)) {
  con_out <- gzfile(extract_file, open = "w")
  writeLines(readLines(old_extract), con_out, useBytes = TRUE)
  close(con_out)
  cat("Сжата готовая вырезка:", old_extract, "->", extract_file, "\n")
}

if (!file.exists(extract_file) && !file.exists(panel_file)) {
  stop("Нет ни ", extract_file, ", ни ", panel_file, ". Скачайте набор ИНИД №205 ",
       "(https://data.rcsi.science/data-catalog/datasets/205/) и положите panel.csv ",
       "в папку проекта.")
}

if (!file.exists(extract_file)) {
  t_start <- Sys.time()
  con <- file(panel_file, open = "r")
  out <- gzfile(extract_file, open = "w")   # пишем сразу в сжатом виде
  writeLines(readLines(con, n = 1), out, useBytes = TRUE)   # строка заголовка

  pattern_region <- paste0('"', REGION_ISO, '"')   # "RU-KDA"
  pattern_year   <- paste0('";', YEAR, ';')        # ";2023;  (год — второй столбец)
  n_read <- 0
  n_kept <- 0

  repeat {
    lines <- readLines(con, n = 500000)
    if (length(lines) == 0) break
    keep <- grepl(pattern_region, lines, fixed = TRUE, useBytes = TRUE) &
            grepl(pattern_year,   lines, fixed = TRUE, useBytes = TRUE)
    writeLines(lines[keep], out, useBytes = TRUE)
    n_read <- n_read + length(lines)
    n_kept <- n_kept + sum(keep)
    if (n_read %% 5e6 == 0) {
      cat("Прочитано строк:", format(n_read, big.mark = " "),
          "| оставлено:", format(n_kept, big.mark = " "), "\n")
    }
  }
  close(con)
  close(out)
  cat("Готово. Всего прочитано:", n_read, "оставлено:", n_kept, "\n")
  print(Sys.time() - t_start)
}


# 2. Загрузка и общая очистка

# Все столбцы читаем как текст, типы задаём сами
kda_raw <- read_delim(extract_file, delim = ";",
                      col_types = cols(.default = col_character()),
                      na = c("", "NA"))
dim(kda_raw)

kda_all <- kda_raw %>%
  # точная проверка после текстового фильтра
  filter(region_iso_code == REGION_ISO, year == as.character(YEAR)) %>%
  mutate(
    year            = as.integer(year),
    kind            = as.integer(kind),        # 1 — юрлицо, 2 — ИП, 3 — самозанятый
    category        = as.integer(category),    # 1 — микро, 2 — малое, 3 — среднее
    lat             = as.numeric(lat),
    lon             = as.numeric(lon),
    revenue         = as.numeric(revenue),     # доходы за год, руб.
    expenditure     = as.numeric(expenditure), # расходы за год, руб.
    employees_count = as.numeric(employees_count),
    start_date      = as.Date(start_date),
    end_date        = as.Date(end_date),
    oktmo           = fix_oktmo(oktmo),
    settlement      = str_squish(settlement),
    activity_code_main = str_trim(activity_code_main)
  )

# Один ИНН может встречаться в году несколько раз, если в течение года менялись
# его сведения (адрес, категория, ОКВЭД). Оставляем последнее состояние за год —
# запись с самой поздней датой окончания (при равенстве — с самым поздним началом).
n_rows_before_dedup <- nrow(kda_all)
kda <- kda_all %>%
  arrange(tin, desc(end_date), desc(start_date)) %>%
  distinct(tin, .keep_all = TRUE)

cat("Строк до удаления повторов:", n_rows_before_dedup,
    "| уникальных ИНН:", nrow(kda), "\n")

# Финансовые показатели есть только у юрлиц (у ИП и самозанятых их нет)
kda %>%
  group_by(kind) %>%
  summarise(n = n(),
            share_with_revenue = mean(!is.na(revenue)),
            share_with_employees = mean(!is.na(employees_count)))


# 3. Топ-3 города края по населению и по экономике

if (!file.exists(settl_file)) {
  stop("Не найден файл ", settl_file, ". Скачайте набор ИНИД №160 ",
       "«Населенные пункты России» (https://data.rcsi.science/data-catalog/datasets/160/) ",
       "и положите его в рабочую папку под этим именем.")
}

# Разделитель в справочнике определяем по первой строке (; , или табуляция)
first_line <- readLines(settl_file, n = 1, encoding = "UTF-8")
delims <- c(";", ",", "\t")
settl_delim <- delims[which.max(str_count(first_line, fixed(delims)))]

settl_raw <- read_delim(settl_file, delim = settl_delim,
                        col_types = cols(.default = col_character()))

need_cols <- c("region", "settlement", "type", "population",
               "latitude_dd", "longitude_dd", "oktmo")
if (!all(need_cols %in% names(settl_raw))) {
  stop("В справочнике нет столбцов: ",
       paste(setdiff(need_cols, names(settl_raw)), collapse = ", "),
       ". Есть столбцы: ", paste(names(settl_raw), collapse = ", "))
}

settl <- settl_raw %>%
  filter(str_detect(region, REGION_NAME)) %>%
  transmute(settlement_160 = str_squish(settlement),
            type_160       = type,
            population     = as.numeric(population),
            lat_160        = as.numeric(latitude_dd),
            lon_160        = as.numeric(longitude_dd),
            oktmo          = fix_oktmo(oktmo),
            key            = name_key(settlement))
nrow(settl)   # число населённых пунктов края в справочнике

# (а) По численности населения — из справочника №160
top_pop <- settl %>%
  arrange(desc(population)) %>%
  mutate(rank_pop = row_number()) %>%
  select(rank_pop, city = settlement_160, type_160, oktmo, population) %>%
  head(10)
top_pop

# (б) По экономике. ВВП (ВРП) Росстат считает только для регионов целиком,
# по городам официальной оценки нет. Поэтому используем прокси из тех же данных:
# суммарные доходы МСП-юрлиц, зарегистрированных в населённом пункте, за 2023 год,
# и число субъектов МСП.
# Группируем по ОКТМО и названию: в panel.csv у некоторых сёл стоит код
# соседнего города (например, у с. Цемдолина — код Новороссийска).
top_econ <- kda %>%
  filter(!is.na(oktmo)) %>%
  group_by(oktmo, city = settlement) %>%
  summarise(n_sme       = n(),                                   # все МСП (ЮЛ+ИП)
            revenue_bln = sum(revenue[kind == 1], na.rm = TRUE) / 1e9,
            .groups = "drop") %>%
  arrange(desc(revenue_bln)) %>%
  mutate(rank_econ = row_number()) %>%
  select(rank_econ, city, oktmo, revenue_bln, n_sme) %>%
  head(10)
top_econ

# Сводная таблица по городам из обоих рейтингов
top_cities <- full_join(select(top_pop, city, oktmo, rank_pop, population),
                        select(top_econ, city_econ = city, oktmo, rank_econ,
                               revenue_bln, n_sme),
                        by = "oktmo") %>%
  mutate(city = coalesce(city, city_econ)) %>%
  select(-city_econ) %>%
  filter(rank_pop <= 5 | rank_econ <= 5) %>%
  arrange(rank_pop, rank_econ)
top_cities

# Центры — города из топ-3 хотя бы одного рейтинга
# (если рейтинги совпадают, центров ровно три)
center_oktmo <- union(top_pop$oktmo[1:3], top_econ$oktmo[1:3])

# Короткие латинские имена для названий переменных
slugs <- c("Краснодар" = "krasnodar", "Сочи" = "sochi",
           "Новороссийск" = "novorossiysk", "Армавир" = "armavir",
           "Анапа" = "anapa", "Геленджик" = "gelendzhik")

# Координаты центров берём так же, как у фирм в panel.csv: тогда у фирм,
# зарегистрированных в самом центре, расстояние будет ровно 0.
# Для сравнения рядом — координаты из справочника №160.
centers <- kda %>%
  filter(oktmo %in% center_oktmo, !is.na(lat)) %>%
  count(oktmo, settlement, lat, lon) %>%
  group_by(oktmo) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(key = name_key(settlement)) %>%
  left_join(settl %>%
              distinct(oktmo, key, .keep_all = TRUE) %>%
              select(oktmo, key, population, lat_160, lon_160),
            by = c("oktmo", "key")) %>%
  mutate(slug = ifelse(settlement %in% names(slugs), slugs[settlement],
                       paste0("center", row_number())),
         gap_vs_160_km = haversine_km(lat, lon, lat_160, lon_160)) %>%
  arrange(desc(population)) %>%
  select(settlement, slug, oktmo, population, lat, lon, gap_vs_160_km)
centers

write_excel_csv(top_cities, "data/top3_cities.csv")


# 4. Выборка для регрессии
# Отбор по шагам; число наблюдений после каждого шага пишем в журнал
s0 <- kda
s1 <- filter(s0, kind == 1)                                    # юридические лица
s2 <- filter(s1, category %in% c(1, 2))                        # микро и малые
# ОКВЭД с подкодами: 47.11, 47.11.1 (замороженные продукты), 47.11.2 (незамороженные),
# 47.11.3 (большой ассортимент с преобладанием продуктов); 45.20, 45.20.1 (легковые),
# 45.20.2 (прочие ТС), 45.20.3 (мойка, полировка), 45.20.4 (техпомощь на дорогах)
s3 <- filter(s2, str_detect(activity_code_main, "^47\\.11") |  # продуктовая розница
                 str_detect(activity_code_main, "^45\\.20"))   # ТО и ремонт авто
s4 <- filter(s3, !is.na(revenue), revenue > 0)                 # положительные доходы
s5 <- filter(s4, !is.na(employees_count), employees_count > 0) # есть работники
s6 <- filter(s5, !is.na(lat), !is.na(lon),                     # корректные координаты
                 between(lat, 43, 47.5), between(lon, 36, 42)) # (рамка края)

cleaning_log <- tibble(
  step = c("0. Все МСП края за 2023 г. (уникальные ИНН)",
           "1. Только юридические лица",
           "2. Микро- и малые предприятия",
           "3. ОКВЭД 47.11* или 45.20*",
           "4. Доходы > 0",
           "5. Численность работников > 0",
           "6. Корректные координаты"),
  n = c(nrow(s0), nrow(s1), nrow(s2), nrow(s3), nrow(s4), nrow(s5), nrow(s6))
)
cleaning_log

reg <- s6 %>%
  mutate(
    industry  = factor(ifelse(str_detect(activity_code_main, "^47\\.11"), "retail", "auto"),
                       levels = c("retail", "auto")),   # база — продуктовая розница
    auto      = as.integer(industry == "auto"),         # 1 — автосервис
    size      = factor(category, levels = c(1, 2), labels = c("micro", "small")),
    small     = as.integer(category == 2),              # 1 — малое предприятие
    rev_per_worker = revenue / employees_count,         # доходы на работника, руб.
    log_rpw   = log(rev_per_worker),                    # зависимая переменная
    log_emp   = log(employees_count)
  )

# Флаг выбросов: доходы на работника вне 1–99 перцентилей своей отрасли.
# Наблюдения НЕ удаляем — флаг нужен для проверки устойчивости.
reg <- reg %>%
  group_by(industry) %>%
  mutate(outlier_1_99 = as.integer(
    rev_per_worker < quantile(rev_per_worker, 0.01) |
    rev_per_worker > quantile(rev_per_worker, 0.99))) %>%
  ungroup()


# 5. Переменные расстояния до центров

# Для каждого центра: расстояние по прямой (км) и log(1 + расстояние)
for (i in seq_len(nrow(centers))) {
  d_name <- paste0("dist_", centers$slug[i])
  reg[[d_name]] <- haversine_km(reg$lat, reg$lon, centers$lat[i], centers$lon[i])
  reg[[paste0("ln_", d_name)]] <- log(1 + reg[[d_name]])
}

# Расстояние до ближайшего из центров и какой это центр
dist_cols <- paste0("dist_", centers$slug)
dist_mat  <- as.matrix(reg[, dist_cols])
reg$dist_nearest    <- apply(dist_mat, 1, min)
reg$nearest_center  <- centers$settlement[apply(dist_mat, 1, which.min)]
reg$ln_dist_nearest <- log(1 + reg$dist_nearest)
reg$in_center       <- as.integer(reg$dist_nearest == 0)   # фирма в самом центре

# Проверка: расстояния между самими центрами (Краснодар — Сочи ≈ 170 км по прямой)
round(outer(seq_len(nrow(centers)), seq_len(nrow(centers)), function(i, j)
  haversine_km(centers$lat[i], centers$lon[i], centers$lat[j], centers$lon[j])))


# 6. Мердж с численностью населения населённого пункта

# Шаг 1: по коду ОКТМО И названию. Одного ОКТМО мало: в panel.csv у некоторых
# сёл стоит код соседнего города (с. Цемдолина -> код Новороссийска), и без
# проверки названия селу досталось бы население всего города.
reg <- reg %>% mutate(key = name_key(settlement))

settl_by_oktmo <- settl %>%
  filter(!is.na(oktmo)) %>%
  group_by(oktmo, key) %>%
  filter(n() == 1) %>%                 # пары, которые в справочнике не повторяются
  ungroup() %>%
  select(oktmo, key, settl_pop_oktmo = population, settlement_160_oktmo = settlement_160)

reg <- reg %>%
  left_join(settl_by_oktmo, by = c("oktmo", "key"))

# Шаг 2: для оставшихся — по названию населённого пункта, а среди одноимённых
# выбираем ближайший по координатам (не дальше 10 км)
unmatched <- reg %>%
  filter(is.na(settl_pop_oktmo)) %>%
  distinct(settlement, lat, lon) %>%
  mutate(key = name_key(settlement))

by_name <- unmatched %>%
  inner_join(select(settl, key, population, lat_160, lon_160, settlement_160),
             by = "key", relationship = "many-to-many") %>%
  mutate(d = haversine_km(lat, lon, lat_160, lon_160)) %>%
  filter(d <= 10) %>%
  group_by(settlement, lat, lon) %>%
  slice_min(d, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(settlement, lat, lon, settl_pop_name = population,
         settlement_160_name = settlement_160)

reg <- reg %>%
  left_join(by_name, by = c("settlement", "lat", "lon")) %>%
  mutate(
    settl_population = coalesce(settl_pop_oktmo, settl_pop_name),
    settlement_160   = coalesce(settlement_160_oktmo, settlement_160_name),
    pop_match        = case_when(!is.na(settl_pop_oktmo) ~ "oktmo+name",
                                 !is.na(settl_pop_name)  ~ "name+coords",
                                 TRUE                    ~ "not matched"),
    ln_settl_pop     = log(settl_population)
  ) %>%
  select(-key, -settl_pop_oktmo, -settl_pop_name,
         -settlement_160_oktmo, -settlement_160_name)

table(reg$pop_match)
# Какие пункты не нашлись (если есть — проверить названия вручную)
reg %>% filter(pop_match == "not matched") %>% count(settlement, oktmo, sort = TRUE)


# 7. Итоговый датасет

# Идентификатор населённого пункта для кластеризации стандартных ошибок.
# Только ОКТМО не годится: у некоторых сёл код соседнего города.
reg <- reg %>% mutate(settlement_id = paste(oktmo, settlement, sep = "_"))

reg_final <- reg %>%
  select(tin, org_short_name, industry, auto, activity_code_main,
         category, size, small,
         settlement, settlement_type, oktmo, settlement_id, lat, lon,
         revenue, expenditure, employees_count,
         rev_per_worker, log_rpw, log_emp, outlier_1_99,
         starts_with("dist_"), starts_with("ln_dist_"),
         nearest_center, in_center,
         settl_population, ln_settl_pop, pop_match) %>%
  arrange(industry, dist_nearest)

write_csv(reg_final, "data/krasnodar_2023_regression.csv")
saveRDS(reg_final, "data/krasnodar_2023_regression.rds")
write_excel_csv(cleaning_log, "data/cleaning_log.csv")

# Краткая проверка результата
dim(reg_final)
table(reg_final$industry, reg_final$size)
n_distinct(reg_final$settlement_id)   # число населённых пунктов (кластеры для ст. ошибок)
summary(select(reg_final, rev_per_worker, employees_count,
               dist_krasnodar, dist_nearest, settl_population))
