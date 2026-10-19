# =============================================================================
# Hospitalizations for ambulatory care sensitive conditions (ACSC) in Ecuador,
# 2000-2024. Complete analysis script.
#
# Martínez-Abarca AL, Henríquez-Trujillo MC, Bonilla E, Henríquez-Trujillo AR.
# Repository: https://github.com/rhenriquezitg/acsc_ecuador
# Archive: https://doi.org/10.5281/zenodo.23366479
#
# The script reproduces every number, table and figure of the article and of
# its supplementary material from the aggregated data in data/. Run it from the
# repository root:
#     Rscript acsc_ecuador_analysis.R
# Run time is a few minutes. Labels use UTF-8 characters: on Linux or macOS run R in a
# UTF-8 locale (e.g. LANG=C.UTF-8); R >= 4.2 on Windows uses UTF-8 by default.
#
# Sections
#   0. Settings and packages
#   1. ACSC definition (20 groups) and provinces
#   2. Functions: classification, rates, joinpoint, segmented regression, spatial
#   3. Optional: rebuild data/ from the INEC microdata (needs the inectab package)
#   4. Load the aggregated data
#   5. Descriptive tables
#   6. National series and joinpoint regression
#   7. Segmented regression: 2008 health reform and COVID-19
#   8. Provinces: rates, variation and spatial autocorrelation
#   9. Age, sex and leading diagnoses
#  10. Sensitivity analyses
#  10b. Hospital bed supply by province, 2004-2024
#  11. Figures (article and supplementary material)
#
# Data sources
#   - Hospital discharges (Egresos Hospitalarios / Camas y Egresos Hospitalarios),
#     Instituto Nacional de Estadistica y Censos (INEC), 2000-2024, public and
#     private establishments. Open microdata at https://www.ecuadorencifras.gob.ec
#   - Population: INEC projections by province, sex and single year of age
#     (revision 2024), provinces in current boundaries.
#   - Hospital beds: INEC, Registro Estadistico de Camas y Egresos Hospitalarios, historical
#     series 2024 (Tabulados_series_ECEH_2024.xlsx, series 1.4.2, 1.5.2 and 1.5.3).
#   - Provincial boundaries: ADM1 2024 (ecu_adm_2024).
# =============================================================================


# ---- 0. Settings and packages -----------------------------------------------

REBUILD_FROM_MICRODATA <- FALSE              # TRUE only with the inectab cache (section 3)
CACHE_DIR <- "C:/Datasets/ecuador/cache"     # inectab parquet cache of INEC microdata
SHAPEFILE <- "C:/Datasets/ecuador/ecu_adm_2024/ecu_adm_adm1_2024.shp"

DATA <- "data"; RES <- "results"; FIG <- "figures"
for (d in c(DATA, RES, FIG)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

YEARS <- 2000:2024
NPERM <- 4499            # permutations of the joinpoint test (999 for single groups)
SEED  <- 2026

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2); library(patchwork)
  library(sf); library(spdep); library(sandwich)
})
# Also used: ragg (figure devices), scales.

write_res <- function(d, name) readr::write_csv(as.data.frame(d), file.path(RES, paste0(name, ".csv")), na = "")


# ---- 1. ACSC definition and provinces ---------------------------------------
# Brazilian list (Alfradique et al. 2009) as adapted for Ecuador. Category is
# used for aggregated analyses. Main diagnosis only.

GROUPS <- tibble::tribble(
  ~group, ~name,                                                                ~icd10,                                                       ~category,
   1L, "Anemia",                                                                "D50",                                                        "Other acute",
   2L, "Angina de pecho",                                                       "I20",                                                        "Chronic",
   3L, "Asma",                                                                  "J45-J46",                                                    "Chronic",
   4L, "Deficiencias nutricionales",                                            "E40-E46, E50-E64",                                           "Other acute",
   5L, "Diabetes mellitus",                                                     "E10-E14",                                                    "Chronic",
   6L, "Enfermedad inflamatoria de \u00f3rganos p\u00e9lvicos femeninos",      "N70-N73, N75-N76",                                           "Other acute",
   7L, "Enfermedades cerebrovasculares",                                        "I63-I67, I69, G45-G46",                                      "Chronic",
   8L, "Enfermedades de v\u00edas a\u00e9reas inferiores",                      "J20-J21, J40-J44, J47",                                      "Chronic",
   9L, "Enfermedades prevenibles por inmunizaci\u00f3n",                        "A17.0, A19, A33-A37, A95, B05-B06, B16, B26, G00.0",         "Vaccine-preventable",
  10L, "Enfermedades relacionadas con la atenci\u00f3n prenatal y del parto",   "A50, O23, P35.0",                                            "Infectious",
  11L, "Epilepsia",                                                             "G40-G41",                                                    "Other acute",
  12L, "Gastroenteritis infecciosas y complicaciones",                          "A00-A09, E86",                                               "Infectious",
  13L, "Hipertensi\u00f3n arterial",                                            "I10-I11",                                                    "Chronic",
  14L, "Infecciones de garganta, nariz y o\u00eddos",                           "H66, J00-J03, J06, J31",                                     "Infectious",
  15L, "Infecciones de la piel y tejido celular subcut\u00e1neo",               "A46, L01-L04, L08",                                          "Infectious",
  16L, "Infecciones renales y de v\u00edas urinarias",                          "N10-N12, N30, N34, N39.0",                                   "Infectious",
  17L, "Insuficiencia card\u00edaca",                                           "I50, J81",                                                   "Chronic",
  18L, "Neumon\u00edas bacterianas",                                            "J13-J14, J15.3-J15.4, J15.8-J15.9, J18.1",                   "Vaccine-preventable",
  19L, "\u00dalcera gastrointestinal",                                          "K25-K28, K92.0-K92.2",                                       "Other acute",
  20L, "Otras condiciones evitables",                                           "A15-A16, A17.1-A17.9, A18, A51-A53, B50-B54, B77, I00-I02", "Vaccine-preventable"
)
GROUPS$label <- sprintf("%02d %s", GROUPS$group, GROUPS$name)
GROUPS_EN <- c("Anemia", "Angina pectoris", "Asthma", "Nutritional deficiencies", "Diabetes mellitus",
               "Female pelvic inflammatory disease", "Cerebrovascular diseases", "Lower respiratory tract diseases",
               "Immunization-preventable diseases", "Prenatal and childbirth care conditions", "Epilepsy",
               "Infectious gastroenteritis", "Hypertension", "Ear, nose and throat infections",
               "Skin and subcutaneous infections", "Kidney and urinary tract infections", "Heart failure",
               "Bacterial pneumonia", "Gastrointestinal ulcer", "Other avoidable conditions")

REGION_TAB <- tibble::tibble(
  province = sprintf("%02d", 1:24),
  province_name = c("Azuay", "Bol\u00edvar", "Ca\u00f1ar", "Carchi", "Cotopaxi", "Chimborazo", "El Oro",
                    "Esmeraldas", "Guayas", "Imbabura", "Loja", "Los R\u00edos", "Manab\u00ed", "Morona Santiago",
                    "Napo", "Pastaza", "Pichincha", "Tungurahua", "Zamora Chinchipe", "Gal\u00e1pagos",
                    "Sucumb\u00edos", "Orellana", "Santo Domingo de los Ts\u00e1chilas", "Santa Elena"),
  region = c("Sierra", "Sierra", "Sierra", "Sierra", "Sierra", "Sierra", "Costa", "Costa", "Costa", "Sierra",
             "Sierra", "Costa", "Costa", "Amazon\u00eda", "Amazon\u00eda", "Amazon\u00eda", "Sierra", "Sierra",
             "Amazon\u00eda", "Gal\u00e1pagos", "Amazon\u00eda", "Amazon\u00eda", "Costa", "Costa"))

PROVS <- REGION_TAB$province
SEX_LAB <- c("1" = "Hombre", "2" = "Mujer")
CATEGORIES <- c("Chronic", "Infectious", "Vaccine-preventable", "Other acute")
REGIONS <- c("Costa", "Sierra", "Amazon\u00eda", "Gal\u00e1pagos")
DELIVERY <- "O80-O84"     # deliveries: removed from the denominator in the "excluding deliveries" share
GROUP_DELIVERY <- 99L     # group code given to delivery records (non-ACSC)

# Three-character categories of which only some subcategories are ACSC (e.g. J18 when
# only J18.1 counts). Records coded with the bare category cannot be classified and are
# counted as non-ACSC; the record log reports them.
partial_categories <- function(groups = GROUPS) {
  tok <- trimws(unlist(strsplit(groups$icd10, "[,-]")))
  unique(substr(tok[grepl("\\.", tok)], 1, 3))
}


# ---- 2. Functions -----------------------------------------------------------

# 2.1 Classification (used only when rebuilding from microdata) ---------------

classify_codes <- function(codes, groups = GROUPS) {
  codes <- unique(codes[!is.na(codes)])
  hit <- vapply(groups$icd10, function(s) inectab::icd_match(codes, s), logical(length(codes)))
  if (is.null(dim(hit))) hit <- matrix(hit, nrow = length(codes))
  multi <- rowSums(hit) > 1
  if (any(multi)) stop("Codes in more than one ACSC group: ", paste(codes[multi], collapse = ", "))
  grp <- ifelse(rowSums(hit) == 0, 0L, groups$group[max.col(hit, ties.method = "first")])
  stats::setNames(as.integer(grp), codes)
}

# Aggregates of INEC egresos by year, province (residence and occurrence), sex,
# single year of age and ACSC group, plus monthly national counts and ICD-10 detail.
build_aggregates <- function(cache_dir, years = YEARS, geo_scheme = "dpa_2024", groups = GROUPS) {
  ds <- inectab::open_cache(cache_dir)
  partial <- partial_categories(groups)
  annual <- list(); monthly <- list(); icd <- list(); logs <- list()
  for (y in years) {
    d <- ds |>
      filter(.data$dataset == "egresos", .data$event_year == y) |>
      select("event_year", "event_month", "sex", "age_years", "prov_res", "cant_res",
             "prov_occ", "cant_occ", "icd10") |>
      collect()
    d <- inectab::apply_geo_scheme(d, geo_scheme)   # provinces in current (2024) boundaries
    lk <- classify_codes(d$icd10, groups)
    d$group <- unname(lk[d$icd10]); d$group[is.na(d$group)] <- 0L
    d$group[d$group == 0L & inectab::icd_match(d$icd10, DELIVERY)] <- GROUP_DELIVERY
    d$age <- as.integer(pmin(floor(d$age_years), 100L))
    for (g in c("res", "occ")) {
      a <- count(tibble(year = d$event_year, province = d[[paste0("prov_", g)]], sex = d$sex,
                        age = d$age, group = d$group), year, province, sex, age, group, name = "n")
      annual[[length(annual) + 1]] <- mutate(a, geo = g, .before = 1)
    }
    monthly[[length(monthly) + 1]] <- count(tibble(year = d$event_year, month = d$event_month, age = d$age,
                                                   group = d$group), year, month, age, group, name = "n")
    h <- d[d$group %in% groups$group, ]
    icd[[length(icd) + 1]] <- count(tibble(year = h$event_year, province = h$prov_res, sex = h$sex, age = h$age,
                                           group = h$group, icd10 = h$icd10),
                                    year, province, sex, age, group, icd10, name = "n")
    lg <- bind_rows(
      tibble(year = y, issue = "3-character code whose category is only partly ACSC (counted as non-ACSC)",
             value = d$icd10[nchar(d$icd10) == 3 & d$icd10 %in% partial]),
      tibble(year = y, issue = "missing ICD-10 code", value = d$icd10[is.na(d$icd10)]),
      tibble(year = y, issue = "missing or unknown age", value = as.character(d$age[is.na(d$age)])),
      tibble(year = y, issue = "sex not Hombre/Mujer (excluded from rates)", value = d$sex[is.na(d$sex) | !d$sex %in% c("1", "2")]),
      tibble(year = y, issue = "province of residence outside the 24 provinces (excluded from rates)",
             value = d$prov_res[is.na(d$prov_res) | !d$prov_res %in% PROVS]),
      tibble(year = y, issue = "missing discharge month", value = as.character(d$event_month[is.na(d$event_month)])))
    lg$value[is.na(lg$value)] <- "NA"
    logs[[length(logs) + 1]] <- count(lg, year, issue, value, name = "n")
    message(y, ": ", nrow(d), " records, ", sum(d$group %in% groups$group), " ACSC")
  }
  list(annual = bind_rows(annual), monthly = bind_rows(monthly), icd = bind_rows(icd), log = bind_rows(logs))
}

# Provincial geometry; neighbours on the full-resolution polygons, simplified for display.
prepare_geometry <- function(shp, tolerance = 250) {
  sf::sf_use_s2(FALSE)
  s <- sf::st_read(shp, quiet = TRUE, options = "ENCODING=UTF-8")
  s$province <- substr(s$ADM1_PCODE, 3, 4)
  s <- sf::st_make_valid(s[, c("province", "ADM1_ES")]); s <- s[order(s$province), ]
  nb_queen <- spdep::poly2nb(s, queen = TRUE, snap = 1e-4)
  cent <- suppressWarnings(sf::st_coordinates(sf::st_point_on_surface(s)))
  nb_knn <- spdep::knn2nb(spdep::knearneigh(cent, k = 4, longlat = TRUE))
  sp <- sf::st_simplify(sf::st_transform(s, 32717), preserveTopology = TRUE, dTolerance = tolerance)
  list(geo = sf::st_transform(sp, 4326), nb = list(queen = nb_queen, knn4 = nb_knn), province = s$province)
}

cut_broad <- function(age) {
  factor(ifelse(is.na(age), "Unknown", ifelse(age < 15, "0-14", ifelse(age < 65, "15-64", "65+"))),
         levels = c("0-14", "15-64", "65+", "Unknown"))
}

# 2.2 Rates --------------------------------------------------------------------
# Direct age standardization with 18 five-year bands (85+ open), Ecuador 2022
# projected population as standard; gamma confidence intervals (Fay and Feuer 1997).

make_bands <- function(amin = 0, amax = 100, width = 5, top = 85) {
  lower <- sort(unique(c(amin, seq(ceiling(amin / width) * width, min(top, amax), width))))
  lower <- lower[lower <= amax]
  upper <- c(lower[-1] - 1, amax)
  lab <- ifelse(upper >= 100 & amax == 100 & seq_along(lower) == length(lower), paste0(lower, "+"),
                ifelse(lower == upper, as.character(lower), paste0(lower, "-", upper)))
  tibble(lower = lower, upper = upper, label = lab)
}

band_of <- function(age, bands) {
  out <- bands$label[pmax(findInterval(age, bands$lower), 1)]
  out[is.na(age)] <- NA
  factor(out, levels = bands$label)
}

std_weights <- function(standard, bands) {
  s <- STD[STD$standard == standard & STD$age >= min(bands$lower) & STD$age <= max(bands$upper), ]
  s$age_band <- as.character(band_of(s$age, bands))
  w <- tapply(s$weight, s$age_band, sum)
  tibble(age_band = names(w), weight = as.numeric(w) / sum(w))
}

pois_ci <- function(n, pop, per, level = 0.95) {
  a <- (1 - level) / 2
  lo <- ifelse(n == 0, 0, stats::qchisq(a, 2 * n) / 2)
  hi <- stats::qchisq(1 - a, 2 * (n + 1)) / 2
  list(lcl = lo / pop * per, ucl = hi / pop * per)
}

# Records usable for rates: known age and sex, residence in one of the 24 provinces.
rate_records <- function(cnt) cnt[!is.na(cnt$age) & cnt$sex %in% c("1", "2") & cnt$province %in% PROVS, ]

# cnt: counts with year, province, sex, age, n (+ columns named in `by`)
# by : output dimensions (year, period, province, region, sex, group, category, age_band)
compute_rates <- function(cnt, pop, by, levels, bands, standard = "ecuador_2022", per = 1e5) {
  cnt <- rate_records(cnt)
  cnt$age_band <- as.character(band_of(cnt$age, bands))
  pop$age_band <- as.character(band_of(pop$age, bands))
  popkeys <- intersect(by, c("year", "period", "province", "region", "sex"))
  ck <- unique(c(by, "age_band"))
  c2 <- summarise(group_by(cnt, across(all_of(ck))), n = sum(.data$n), .groups = "drop")
  p2 <- summarise(group_by(pop, across(all_of(unique(c(popkeys, "age_band"))))), pop = sum(.data$pop), .groups = "drop")
  lv <- levels[setdiff(by, "age_band")]; lv$age_band <- bands$label
  d <- left_join(as_tibble(expand.grid(lv, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)), c2, by = ck)
  d$n[is.na(d$n)] <- 0
  d <- left_join(d, p2, by = unique(c(popkeys, "age_band")))
  if ("age_band" %in% by) {
    ci <- pois_ci(d$n, d$pop, per)
    d$rate <- d$n / d$pop * per; d$rate_lcl <- ci$lcl; d$rate_ucl <- ci$ucl
    d$age_band <- factor(d$age_band, levels = bands$label)
    return(arrange(d, across(all_of(by))))
  }
  d <- left_join(d, std_weights(standard, bands), by = "age_band")
  out <- summarise(group_by(d, across(all_of(by))),
                   asr = sum(.data$weight * .data$n / .data$pop) * per,
                   var = sum((.data$weight / .data$pop)^2 * .data$n) * per^2,
                   wm = max(.data$weight / .data$pop) * per,
                   n_ = sum(.data$n), pop_ = sum(.data$pop), .groups = "drop")
  out$asr_lcl <- ifelse(out$var > 0, stats::qgamma(0.025, shape = out$asr^2 / out$var, scale = out$var / out$asr), 0)
  out$asr_ucl <- stats::qgamma(0.975, shape = (out$asr + out$wm)^2 / (out$var + out$wm^2),
                               scale = (out$var + out$wm^2) / (out$asr + out$wm))
  out$asr_se <- sqrt(out$var)
  out <- rename(out, n = "n_", pop = "pop_")
  ci <- pois_ci(out$n, out$pop, per)
  out$crude <- out$n / out$pop * per; out$crude_lcl <- ci$lcl; out$crude_ucl <- ci$ucl
  out <- select(out, all_of(by), "n", "pop", "crude", "crude_lcl", "crude_ucl", "asr", "asr_lcl", "asr_ucl", "asr_se")
  arrange(out, across(all_of(by)))
}

# Indirectly standardized expected counts by area (for the systematic component of variation).
expected_counts <- function(cnt, pop, area = "province", bands) {
  cnt <- rate_records(cnt)
  cnt$age_band <- as.character(band_of(cnt$age, bands)); pop$age_band <- as.character(band_of(pop$age, bands))
  ref <- left_join(summarise(group_by(cnt, .data$age_band, .data$sex), n = sum(.data$n), .groups = "drop"),
                   summarise(group_by(pop, .data$age_band, .data$sex), pop = sum(.data$pop), .groups = "drop"),
                   by = c("age_band", "sex"))
  ref$r <- ref$n / ref$pop
  pa <- summarise(group_by(pop, across(all_of(c(area, "age_band", "sex")))), pop = sum(.data$pop), .groups = "drop")
  pa <- left_join(pa, ref[, c("age_band", "sex", "r")], by = c("age_band", "sex")); pa$r[is.na(pa$r)] <- 0
  e <- summarise(group_by(pa, across(all_of(area))), expected = sum(.data$pop * .data$r), .groups = "drop")
  o <- summarise(group_by(cnt, across(all_of(area))), observed = sum(.data$n), .groups = "drop")
  out <- left_join(e, o, by = area); out$observed[is.na(out$observed)] <- 0
  out$sir <- out$observed / out$expected
  out
}

# 2.3 Joinpoint regression -------------------------------------------------------
# Weighted log-linear model, ln(rate) = b0 + b1 t + sum_k d_k (t - tau_k)+, weights
# 1/var(ln rate). Joinpoints at observed years, at least two observations from each end
# and between joinpoints; number of joinpoints by sequential permutation tests with
# Bonferroni correction (Kim et al., Stat Med 2000;19:335-51); maximum number of
# joinpoints as in NCI Joinpoint defaults. APC and AAPC with t-based 95% CI.

jp_design <- function(x, tau) { X <- cbind(1, x); for (t in tau) X <- cbind(X, pmax(x - t, 0)); X }

jp_candidates <- function(x, k, min_end = 2, min_between = 2) {
  n <- length(x)
  if (k == 0) return(list(integer()))
  idx <- (min_end + 1):(n - min_end)
  if (length(idx) < k) return(list())
  Filter(function(p) k == 1 || all(diff(p) >= min_between + 1), utils::combn(idx, k, simplify = FALSE))
}

wls_rss <- function(X, y, w) { sw <- sqrt(w); sum(.lm.fit(X * sw, y * sw)$residuals^2) }

jp_fit_k <- function(x, y, w, k, cands = jp_candidates(x, k)) {
  best <- NULL; best_rss <- Inf
  for (p in cands) { r <- wls_rss(jp_design(x, x[p]), y, w); if (r < best_rss) { best_rss <- r; best <- p } }
  list(k = k, pos = best, tau = x[best], rss = best_rss)
}

jp_full <- function(x, y, w, tau) {
  X <- jp_design(x, tau)
  fit <- stats::lm.wfit(X, y, w)
  df <- length(y) - ncol(X) - length(tau)        # estimated joinpoints count as parameters
  s2 <- sum(w * fit$residuals^2) / df
  list(coef = fit$coefficients, V = s2 * solve(crossprod(X * sqrt(w))), df = df,
       fitted = drop(X %*% fit$coefficients), rss = sum(w * fit$residuals^2))
}

# Minimum residual sum of squares over candidate designs, for every permuted series (columns of Yw).
jp_min_rss_batch <- function(Qlist, Yw) {
  best <- rep(Inf, ncol(Yw))
  for (Q in Qlist) best <- pmin(best, colSums((Yw - Q %*% crossprod(Q, Yw))^2))
  best
}

jp_select <- function(x, y, w, kmax = 4, nperm = 4499, alpha = 0.05, seed = SEED) {
  set.seed(seed)
  cands <- lapply(0:kmax, function(k) jp_candidates(x, k))
  fits <- lapply(0:kmax, function(k) jp_fit_k(x, y, w, k, cands[[k + 1]]))
  sw <- sqrt(w)
  Qs <- lapply(cands, function(cl) lapply(cl, function(p) qr.Q(qr(jp_design(x, x[p]) * sw))))
  a <- 0; b <- kmax; tests <- list()
  while (a < b) {
    f0 <- fits[[a + 1]]; f1 <- fits[[b + 1]]
    stat <- (f0$rss - f1$rss) / f1$rss
    full0 <- jp_full(x, y, w, f0$tau)
    r_std <- (y - full0$fitted) * sw
    Yw <- vapply(seq_len(nperm), function(i) full0$fitted * sw + sample(r_std), numeric(length(y)))
    r0 <- jp_min_rss_batch(Qs[[a + 1]], Yw); r1 <- jp_min_rss_batch(Qs[[b + 1]], Yw)
    p <- (sum((r0 - r1) / r1 >= stat) + 1) / (nperm + 1)
    lev <- alpha / kmax
    tests[[length(tests) + 1]] <- data.frame(H0 = a, H1 = b, stat = stat, p = p, level = lev,
                                             decision = if (p < lev) "reject H0" else "keep H0")
    if (p < lev) a <- a + 1 else b <- b - 1
  }
  list(k = a, fits = fits, tests = do.call(rbind, tests))
}

jp_kmax <- function(n) if (n < 7) 0L else if (n <= 11) 1L else if (n <= 16) 2L else if (n <= 21) 3L else if (n <= 26) 4L else 5L

jp_summary <- function(x, y, w, tau, label = "", periods = NULL) {
  f <- jp_full(x, y, w, tau); k <- length(tau); tq <- stats::qt(0.975, f$df)
  L <- matrix(0, k + 1, length(f$coef)); for (j in 0:k) L[j + 1, 2:(2 + j)] <- 1
  slope <- drop(L %*% f$coef); se <- sqrt(diag(L %*% f$V %*% t(L)))
  starts <- c(min(x), tau); ends <- c(tau, max(x))
  seg <- data.frame(series = label, segment = seq_len(k + 1), start = starts, end = ends,
                    apc = 100 * (exp(slope) - 1), apc_lcl = 100 * (exp(slope - tq * se) - 1),
                    apc_ucl = 100 * (exp(slope + tq * se) - 1), p = 2 * stats::pt(-abs(slope / se), f$df))
  aapc_of <- function(from, to) {
    wl <- pmax(0, pmin(ends, to) - pmax(starts, from)); l <- drop(wl %*% L) / sum(wl)
    b <- sum(l * f$coef); s <- sqrt(drop(t(l) %*% f$V %*% l))
    data.frame(series = label, period = paste0(from, "-", to), aapc = 100 * (exp(b) - 1),
               aapc_lcl = 100 * (exp(b - tq * s) - 1), aapc_ucl = 100 * (exp(b + tq * s) - 1),
               p = 2 * stats::pt(-abs(b / s), f$df))
  }
  if (is.null(periods)) periods <- list(c(min(x), max(x)))
  list(segments = seg, aapc = do.call(rbind, lapply(periods, function(pr) aapc_of(pr[1], pr[2]))),
       fitted = data.frame(series = label, year = x, observed = exp(y), fitted = exp(f$fitted)))
}

# d: year, asr, asr_se
run_joinpoint <- function(d, label = "", nperm = NPERM, kmax = NULL, periods = NULL) {
  d <- d[is.finite(d$asr) & d$asr > 0 & d$asr_se > 0, ]; d <- d[order(d$year), ]
  x <- d$year; y <- log(d$asr); w <- (d$asr / d$asr_se)^2
  if (is.null(kmax)) kmax <- jp_kmax(length(x))
  sel <- jp_select(x, y, w, kmax = kmax, nperm = nperm)
  out <- jp_summary(x, y, w, sel$fits[[sel$k + 1]]$tau, label, periods)
  out$tests <- cbind(series = label, sel$tests); out$k <- sel$k
  out
}

# 2.4 Segmented regression (interrupted time series) ---------------------------
# log E[y_t] = b0 + b1 t + b2 post_t + b3 time_since_t + Fourier(month, 2) + log(offset_t)
# Quasi-Poisson, Newey-West standard errors (lag 12). Level change = exp(b2);
# change in trend = exp(12 b3) per year.

monthly_population <- function(pop_year) {           # mid-year values interpolated to mid-month
  pop_year <- pop_year[order(pop_year$year), ]
  x <- pop_year$year + 0.5
  grid <- expand.grid(month = 1:12, year = min(pop_year$year):max(pop_year$year))
  tt <- grid$year + (grid$month - 0.5) / 12
  grid$pop <- stats::approx(x, pop_year$pop, xout = tt, rule = 2)$y
  lo <- tt < min(x); hi <- tt > max(x)
  s1 <- diff(pop_year$pop[1:2]); s2 <- diff(utils::tail(pop_year$pop, 2))
  grid$pop[lo] <- pop_year$pop[1] + s1 * (tt[lo] - x[1])
  grid$pop[hi] <- utils::tail(pop_year$pop, 1) + s2 * (tt[hi] - max(x))
  grid
}

its_frame <- function(series, pop_month, start, k_fourier = 2, exclude_months = 0) {
  d <- left_join(pop_month, series, by = c("year", "month")); d$n[is.na(d$n)] <- 0
  d <- d[order(d$year, d$month), ]
  d$t <- seq_len(nrow(d)); d$date <- as.Date(sprintf("%d-%02d-15", d$year, d$month))
  t0 <- which(d$year == start[1] & d$month == start[2])
  d$post <- as.integer(d$t >= t0); d$time_since <- pmax(0, d$t - t0)
  for (k in seq_len(k_fourier)) {
    d[[paste0("sin", k)]] <- sin(2 * pi * k * d$month / 12)
    d[[paste0("cos", k)]] <- cos(2 * pi * k * d$month / 12)
  }
  d$excluded <- d$t >= t0 & d$t < t0 + exclude_months
  attr(d, "t0") <- t0
  d
}

fit_its <- function(d, k_fourier = 2, nw_lag = 12, per = 1e5, label = "") {
  fr <- c("t", "post", "time_since", paste0(rep(c("sin", "cos"), k_fourier), rep(seq_len(k_fourier), each = 2)))
  m <- stats::glm(stats::as.formula(paste("n ~", paste(fr, collapse = " + "), "+ offset(log(pop))")),
                  family = stats::quasipoisson(), data = d[!d$excluded, ])
  V <- sandwich::NeweyWest(m, lag = nw_lag, prewhite = FALSE, adjust = TRUE)
  b <- stats::coef(m); se <- sqrt(diag(V)); z <- stats::qnorm(0.975)
  eff <- function(term, mult = 1)
    data.frame(estimate = 100 * (exp(mult * b[term]) - 1), lcl = 100 * (exp(mult * (b[term] - z * se[term])) - 1),
               ucl = 100 * (exp(mult * (b[term] + z * se[term])) - 1), p = 2 * stats::pnorm(-abs(b[term] / se[term])))
  res <- rbind(cbind(term = "Pre-interruption trend (% per year)", eff("t", 12)),
               cbind(term = "Immediate level change (%)", eff("post")),
               cbind(term = "Change in trend after interruption (% per year)", eff("time_since", 12)))
  l <- c(t = 12, time_since = 12)
  bb <- sum(b[names(l)] * l); ss <- sqrt(drop(t(l) %*% V[names(l), names(l)] %*% l))
  res <- rbind(res, data.frame(term = "Post-interruption trend (% per year)", estimate = 100 * (exp(bb) - 1),
                               lcl = 100 * (exp(bb - z * ss) - 1), ucl = 100 * (exp(bb + z * ss) - 1),
                               p = 2 * stats::pnorm(-abs(bb / ss))))
  res <- cbind(series = label, res); rownames(res) <- NULL
  d$fitted <- stats::predict(m, newdata = d, type = "response")
  cf <- d; cf$post <- 0; cf$time_since <- 0
  d$counterfactual <- stats::predict(m, newdata = cf, type = "response")
  ds <- d; for (v in grep("^(sin|cos)", names(ds), value = TRUE)) ds[[v]] <- 0
  d$trend <- stats::predict(m, newdata = ds, type = "response")
  ds$post <- 0; ds$time_since <- 0
  d$trend_cf <- stats::predict(m, newdata = ds, type = "response")
  for (v in c("n", "fitted", "counterfactual", "trend", "trend_cf")) d[[paste0(v, "_rate")]] <- d[[v]] / d$pop * per
  post <- d[d$post == 1, ]
  list(effects = res, data = d, label = label, dispersion = summary(m)$dispersion,
       ljung_box_p = stats::Box.test(stats::residuals(m, type = "pearson"), lag = 12, type = "Ljung-Box")$p.value,
       post_observed = sum(post$n), post_expected = sum(post$counterfactual))
}

# 2.5 Geographic variation and spatial autocorrelation --------------------------
# Coefficient of variation, extremal quotient and systematic component of variation
# (McPherson et al. 1982); global and local Moran's I (row-standardized weights).

variation_indicators <- function(asr, exp_tab) {
  k <- nrow(exp_tab)
  scv <- (sum(((exp_tab$observed - exp_tab$expected) / exp_tab$expected)^2) / k - sum(1 / exp_tab$expected) / k) * 100
  data.frame(indicator = c("Provinces", "Minimum ASR", "Maximum ASR", "Extremal quotient (max/min)",
                           "Coefficient of variation (%)", "Systematic component of variation (x100)"),
             value = c(length(asr), min(asr), max(asr), max(asr) / min(asr), 100 * stats::sd(asr) / mean(asr), scv))
}

moran_analysis <- function(values, geo_obj, nb_type = "queen", nsim = 999, seed = SEED) {
  nb <- geo_obj$nb[[nb_type]]
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  set.seed(seed)
  gm <- spdep::moran.mc(values, lw, nsim = nsim, zero.policy = TRUE)
  lm <- spdep::localmoran(values, lw, zero.policy = TRUE)
  q <- as.character(attr(lm, "quadr")$mean); p <- lm[, 5]
  cl <- ifelse(spdep::card(nb) == 0, "No neighbours", ifelse(p < 0.05, q, "Not significant"))
  list(global = data.frame(statistic = "Global Moran's I", value = unname(gm$statistic), p_value = gm$p.value,
                           nsim = nsim, neighbours = nb_type),
       local = data.frame(province = geo_obj$province, Ii = lm[, 1], p = p, cluster = cl))
}

# 2.6 Helpers for the results layer -------------------------------------------

SETTINGS <- list(geo = "res", years = YEARS, amin = 0, amax = 100, groups = GROUPS$group,
                 standard = "ecuador_2022", per = 1e5, nb_type = "queen")

prep_annual <- function(dat, s) {
  a <- dat$annual
  a <- a[a$geo == s$geo & a$year %in% s$years, ]
  a$hecsap <- a$group %in% s$groups                     # ACSC indicator
  a$category <- GROUPS$category[match(a$group, GROUPS$group)]
  a$region <- REGION_TAB$region[match(a$province, REGION_TAB$province)]
  a
}

level_list <- function(s) list(year = s$years, province = PROVS, region = REGIONS, sex = c("1", "2"),
                               group = intersect(GROUPS$group, s$groups), category = CATEGORIES)

province_rates <- function(a, pop, s, by_year = FALSE) {
  r <- compute_rates(a[a$hecsap, ], pop, by = c(if (by_year) "year", "province"), levels = level_list(s),
                     bands = make_bands(), standard = s$standard, per = s$per)
  if (!by_year) r$pop <- r$pop / length(s$years)        # mean annual population
  r$province_name <- REGION_TAB$province_name[match(r$province, REGION_TAB$province)]
  r$region <- REGION_TAB$region[match(r$province, REGION_TAB$province)]
  r
}

province_analysis <- function(a, pop, s, geo_obj) {
  pr <- province_rates(a, pop, s)
  ex <- expected_counts(a[a$hecsap, ], pop, "province", make_bands())
  ex <- ex[match(pr$province, ex$province), ]
  pr$expected <- ex$expected; pr$sir <- ex$sir
  var <- variation_indicators(pr$asr, ex)
  mor <- moran_analysis(pr$asr[match(geo_obj$province, pr$province)], geo_obj, s$nb_type)
  pr <- left_join(pr, mor$local, by = "province")
  pr$quintile <- cut(pr$asr, stats::quantile(pr$asr, 0:5 / 5), include.lowest = TRUE, labels = paste0("Q", 1:5))
  pr <- pr[order(-pr$asr), ]; pr$rank <- seq_len(nrow(pr))
  list(table = pr, variation = var, moran = mor$global)
}

age_sex_rates <- function(a, pop, s) {
  compute_rates(a[a$hecsap, ], pop, by = c("sex", "age_band"), levels = level_list(s),
                bands = make_bands(top = 80), standard = s$standard, per = s$per)
}

# Five most frequent diagnoses by broad age group and sex. Three-character roll-up,
# except categories only partly in the list (J15.9, N39.0, K92.2 ...).
top_diagnoses <- function(dat, pop, s, top = 5) {
  x <- dat$icd
  x <- x[x$year %in% s$years & x$group %in% s$groups & !is.na(x$age), ]
  x$age_broad <- as.character(cut_broad(x$age))
  part <- partial_categories()
  x$code <- ifelse(substr(x$icd10, 1, 3) %in% part, x$icd10, substr(x$icd10, 1, 3))
  t <- summarise(group_by(x, .data$age_broad, .data$sex, .data$code), n = sum(.data$n), .groups = "drop")
  t <- t[t$sex %in% c("1", "2"), ]
  tot <- summarise(group_by(t, .data$age_broad, .data$sex), total = sum(.data$n), .groups = "drop")
  pop$age_broad <- as.character(cut_broad(pop$age))
  pp <- summarise(group_by(pop[pop$province %in% PROVS, ], .data$age_broad, .data$sex), pop = sum(.data$pop), .groups = "drop")
  t <- left_join(left_join(t, tot, by = c("age_broad", "sex")), pp, by = c("age_broad", "sex"))
  t <- ungroup(mutate(slice_head(arrange(group_by(t, .data$age_broad, .data$sex), desc(.data$n)), n = top), rank = row_number()))
  t$label <- ICD_LABELS$label[match(t$code, ICD_LABELS$code)]
  t$pct <- 100 * t$n / t$total; t$rate <- t$n / t$pop * s$per; t$sex <- SEX_LAB[t$sex]
  t
}


# ---- 3. Optional: rebuild data/ from the INEC microdata ----------------------
# Needs the inectab package (harmonized INEC files in a parquet cache) and the shapefile.

if (REBUILD_FROM_MICRODATA) {
  agg <- build_aggregates(CACHE_DIR)
  readr::write_csv(agg$annual, gzfile(file.path(DATA, "acsc_annual.csv.gz")), na = "")
  readr::write_csv(agg$monthly, gzfile(file.path(DATA, "acsc_monthly.csv.gz")), na = "")
  readr::write_csv(agg$icd, gzfile(file.path(DATA, "acsc_icd.csv.gz")), na = "")
  readr::write_csv(agg$log, file.path(DATA, "record_log.csv"), na = "")
  p <- inectab::inec_population(); readr::write_csv(p[p$year %in% YEARS, ], gzfile(file.path(DATA, "population.csv.gz")))
  readr::write_csv(bind_rows(cbind(standard = "ecuador_2022", inectab::standard_population("ecuador_2022")),
                             cbind(standard = "who", inectab::standard_population("who"))),
                   file.path(DATA, "standard_population.csv"))
  ref <- inectab::icd10_reference(); codes <- unique(c(substr(agg$icd$icd10, 1, 3), agg$icd$icd10))
  readr::write_csv(ref[ref$code %in% codes, ], file.path(DATA, "icd10_labels.csv"))
  saveRDS(prepare_geometry(SHAPEFILE), file.path(DATA, "provinces_geo.rds"))
  rm(agg)
}


# ---- 4. Load the aggregated data ----------------------------------------------
# annual : geo (res = residence, occ = occurrence), year, province, sex, age, group, n
#          group 1-20 = ACSC group, 99 = delivery (O80-O84), 0 = other causes
# monthly: year, month, age, group, n (national)
# icd    : year, province (residence), sex, age, group, icd10, n (ACSC records only)

dat <- list(
  annual = readr::read_csv(file.path(DATA, "acsc_annual.csv.gz"), show_col_types = FALSE,
                           col_types = cols(geo = "c", year = "i", province = "c", sex = "c", age = "i", group = "i", n = "d")),
  monthly = readr::read_csv(file.path(DATA, "acsc_monthly.csv.gz"), show_col_types = FALSE,
                            col_types = cols(year = "i", month = "i", age = "i", group = "i", n = "d")),
  icd = readr::read_csv(file.path(DATA, "acsc_icd.csv.gz"), show_col_types = FALSE,
                        col_types = cols(year = "i", province = "c", sex = "c", age = "i", group = "i", icd10 = "c", n = "d")),
  log = readr::read_csv(file.path(DATA, "record_log.csv"), show_col_types = FALSE, na = character(), col_types = cols(.default = "c", year = "i", n = "d")))
pop <- readr::read_csv(file.path(DATA, "population.csv.gz"), show_col_types = FALSE,
                       col_types = cols(year = "i", prov = "c", sex = "c", age = "i", pop = "d")) |>
  rename(province = "prov") |> filter(year %in% YEARS)
pop$region <- REGION_TAB$region[match(pop$province, REGION_TAB$province)]
STD <- readr::read_csv(file.path(DATA, "standard_population.csv"), show_col_types = FALSE)
ICD_LABELS <- readr::read_csv(file.path(DATA, "icd10_labels.csv"), show_col_types = FALSE)
geo <- readRDS(file.path(DATA, "provinces_geo.rds"))

s <- SETTINGS
a <- prep_annual(dat, s)
PER <- c("2000\u20132007" = 2000, "2008\u20132014" = 2008, "2015\u20132019" = 2015, "2020\u20132024" = 2020)
period_of <- function(y) names(PER)[findInterval(y, PER)]
a$period <- period_of(a$year); pop$period <- period_of(pop$year)
LV <- c(level_list(s), list(period = names(PER)))
bands <- make_bands()
rate <- function(cnt, by) compute_rates(cnt, pop, by, LV, bands, s$standard)
h <- a[a$hecsap, ]                                   # ACSC records, residence-based
nonp <- function(x) x[x$group != GROUP_DELIVERY, ]


# ---- 5. Descriptive tables -------------------------------------------------------

t1 <- list()
add <- function(char, cat, sel) {
  for (p in c(names(PER), "Total")) {
    hh <- if (p == "Total") h else h[h$period == p, ]
    t1[[length(t1) + 1]] <<- data.frame(char = char, cat = cat, period = p, n = sum(sel(hh)$n), tot = sum(hh$n))
  }
}
add("Total", "All", identity)
add("Sex", "Men", function(x) x[x$sex == "1", ]); add("Sex", "Women", function(x) x[x$sex == "2", ])
add("Age", "0-4", function(x) x[!is.na(x$age) & x$age < 5, ])
add("Age", "5-14", function(x) x[!is.na(x$age) & x$age >= 5 & x$age < 15, ])
add("Age", "15-64", function(x) x[!is.na(x$age) & x$age >= 15 & x$age < 65, ])
add("Age", "65+", function(x) x[!is.na(x$age) & x$age >= 65, ])
for (cc in CATEGORIES) add("Category", cc, function(x) x[x$category %in% cc, ])
for (rg in REGIONS) add("Region", rg, function(x) x[x$region %in% rg, ])
add("Region", "Abroad/ZND", function(x) x[!x$province %in% PROVS, ])
t1 <- bind_rows(t1); t1$pct <- 100 * t1$n / t1$tot
write_res(t1, "table1_counts")

sumrow <- function(sel_years, label) {
  aa <- a[a$year %in% sel_years, ]; pp <- pop[pop$year %in% sel_years & pop$province %in% PROVS, ]
  hh <- aa[aa$hecsap, ]
  r1 <- compute_rates(hh, pp, character(), list(), bands, s$standard)
  data.frame(period = label, acsc = sum(hh$n), all = sum(aa$n), all_nodeliv = sum(nonp(aa)$n),
             share = 100 * sum(hh$n) / sum(aa$n), share_nodeliv = 100 * sum(hh$n) / sum(nonp(aa)$n),
             disch_per_1000 = 1000 * sum(aa$n[aa$province %in% PROVS]) / sum(pp$pop),
             crude = r1$crude, asr = r1$asr, asr_lcl = r1$asr_lcl, asr_ucl = r1$asr_ucl, years = length(sel_years))
}
write_res(bind_rows(lapply(names(PER), function(p) sumrow(s$years[period_of(s$years) == p], p)), sumrow(s$years, "Total")),
          "table1_summary")

ann <- a |> group_by(year) |>
  summarise(acsc = sum(n[hecsap]), all = sum(n), deliv = sum(n[group == GROUP_DELIVERY]), all_res = sum(n[province %in% PROVS])) |>
  left_join(pop |> filter(province %in% PROVS) |> group_by(year) |> summarise(pop = sum(pop)), by = "year") |>
  mutate(share = 100 * acsc / all, share_nodeliv = 100 * acsc / (all - deliv), disch_per_1000 = 1000 * all_res / pop)
write_res(ann, "annual_totals")

crude_p <- function(sel) {
  sapply(c(names(PER), "Total"), function(p) {
    hh <- if (p == "Total") h else h[h$period == p, ]; pp <- if (p == "Total") pop else pop[pop$period == p, ]
    hh <- sel(hh[hh$province %in% PROVS & hh$sex %in% c("1", "2") & !is.na(hh$age), ]); pp <- sel(pp[pp$province %in% PROVS, ])
    sum(hh$n) / sum(pp$pop) * 1e5 })
}
cr <- list(All = crude_p(identity),
           Men = crude_p(function(x) x[x$sex == "1", ]), Women = crude_p(function(x) x[x$sex == "2", ]),
           `0-4` = crude_p(function(x) x[x$age < 5, ]), `5-14` = crude_p(function(x) x[x$age >= 5 & x$age < 15, ]),
           `15-64` = crude_p(function(x) x[x$age >= 15 & x$age < 65, ]), `65+` = crude_p(function(x) x[x$age >= 65, ]))
write_res(bind_rows(lapply(names(cr), function(k) data.frame(stratum = k, t(cr[[k]]), check.names = FALSE))), "table1_crude")

gc <- h |> group_by(group) |> summarise(n = sum(n)) |> mutate(pct = 100 * n / sum(n))
write_res(gc, "group_counts")
write_res(GROUPS[, c("group", "name", "category", "icd10")], "tableS1_codes")
write_res(dat$log, "record_log")
write_res(dat$icd |> filter(group == 10) |> mutate(c3 = substr(icd10, 1, 3)) |> group_by(c3) |> summarise(n = sum(n)), "group10_codes")


# ---- 6. National series and joinpoint regression -------------------------------

nat  <- rate(h, "year");               write_res(nat, "nat_total")
nsex <- rate(h, c("year", "sex"));      write_res(nsex, "nat_sex")
ncat <- rate(h, c("year", "category")); write_res(ncat, "nat_category")
ngrp <- rate(h, c("year", "group"));    write_res(ngrp, "nat_group")
write_res(rate(h, c("year", "region")), "nat_region")
nage <- h |> filter(province %in% PROVS, sex %in% c("1", "2"), !is.na(age)) |> mutate(ab = as.character(cut_broad(age))) |>
  group_by(year, ab) |> summarise(n = sum(n), .groups = "drop") |>
  left_join(pop |> filter(province %in% PROVS) |> mutate(ab = as.character(cut_broad(age))) |>
              group_by(year, ab) |> summarise(pop = sum(pop), .groups = "drop"), by = c("year", "ab")) |>
  mutate(rate = n / pop * 1e5)
write_res(nage, "nat_agebroad")
write_res(rate(h, "period"), "period_total"); write_res(rate(h, c("period", "sex")), "period_sex")
write_res(rate(h, c("period", "category")), "period_category"); write_res(rate(h, c("period", "group")), "period_group")
write_res(rate(h, c("period", "region")), "period_region")

PERIODS_JP <- list(c(2000, 2024), c(2000, 2007), c(2008, 2019), c(2015, 2024))
jp_run <- function(d, label, excl = integer(), nperm = NPERM)
  run_joinpoint(d[!d$year %in% excl, ], label, nperm = nperm, periods = PERIODS_JP)
collect_jp <- function(excl, tag) {
  res <- list(jp_run(nat, "Total", excl))
  for (sx in c("1", "2")) res[[length(res) + 1]] <- jp_run(nsex[nsex$sex == sx, ], SEX_LAB[[sx]], excl)
  for (cc in CATEGORIES) res[[length(res) + 1]] <- jp_run(ncat[ncat$category == cc, ], cc, excl)
  for (g in GROUPS$group) res[[length(res) + 1]] <- jp_run(ngrp[ngrp$group == g, ], GROUPS$label[GROUPS$group == g], excl, nperm = 999)
  # ACSC share of hospitalizations excluding deliveries; binomial standard error
  sh <- ann |> transmute(year, asr = 100 * acsc / (all - deliv),
                         asr_se = 100 * sqrt((acsc / (all - deliv)) * (1 - acsc / (all - deliv)) / (all - deliv)))
  res[[length(res) + 1]] <- jp_run(sh, "Share of hospitalizations", excl)
  write_res(bind_rows(lapply(res, `[[`, "segments")), paste0("jp_segments", tag))
  write_res(bind_rows(lapply(res, `[[`, "aapc")), paste0("jp_aapc", tag))
  write_res(bind_rows(lapply(res, `[[`, "fitted")), paste0("jp_fitted", tag))
  write_res(bind_rows(lapply(res, `[[`, "tests")), paste0("jp_tests", tag))
}
message("Joinpoint regression")
collect_jp(integer(), "")
collect_jp(c(2020L, 2021L), "_excl2021")          # sensitivity: 2020 and 2021 left out


# ---- 7. Segmented regression: 2008 reform and COVID-19 --------------------------

monthly_all <- function(years, grp) {
  m <- dat$monthly; m <- m[m$year %in% years & m$group %in% grp, ]
  summarise(group_by(m, .data$year, .data$month), n = sum(.data$n), .groups = "drop")
}
# offset "pop": population (rates); offset "hosp": all hospitalizations excluding deliveries (share)
its_run <- function(years, start, grp, label, offset = c("pop", "hosp"), excl = 0L) {
  offset <- match.arg(offset)
  ser <- monthly_all(years, grp)
  pm <- if (offset == "pop") monthly_population(pop |> filter(year %in% years) |> group_by(year) |> summarise(pop = sum(pop)))
        else (monthly_all(years, c(0L, GROUPS$group)) |> rename(pop = n))[, c("month", "year", "pop")]
  f <- fit_its(its_frame(ser, pm, start = start, k_fourier = 2, exclude_months = excl), 2, 12, 1e5, label)
  f$effects$offset <- offset
  f
}
CATG <- lapply(setNames(CATEGORIES, CATEGORIES), function(cc) GROUPS$group[GROUPS$category == cc])
designs <- list(reform = list(years = 2000:2019, start = c(2008L, 1L)),
                covid  = list(years = 2012:2024, start = c(2020L, 3L)))
message("Segmented regression")
fits <- list(); eff <- list(); diag <- list()
keep_diag <- function(dn, f, offset) data.frame(design = dn, series = f$label, offset = offset, dispersion = f$dispersion,
                                                lb = f$ljung_box_p, obs = f$post_observed, exp = f$post_expected)
for (dn in names(designs)) {
  ds_ <- designs[[dn]]
  sers <- c(list(`All ACSC` = GROUPS$group), CATG, list(`Non-ACSC (excluding deliveries)` = 0L))
  for (sn in names(sers)) {
    f <- its_run(ds_$years, ds_$start, sers[[sn]], sn, "pop")
    fits[[paste(dn, sn, "pop")]] <- f
    eff[[length(eff) + 1]] <- cbind(design = dn, f$effects); diag[[length(diag) + 1]] <- keep_diag(dn, f, "pop")
  }
  f <- its_run(ds_$years, ds_$start, GROUPS$group, "ACSC share of hospitalizations", "hosp")
  fits[[paste(dn, "share")]] <- f
  eff[[length(eff) + 1]] <- cbind(design = dn, f$effects); diag[[length(diag) + 1]] <- keep_diag(dn, f, "hosp")
  for (g in GROUPS$group) {
    f <- tryCatch(its_run(ds_$years, ds_$start, g, GROUPS$label[GROUPS$group == g], "pop"), error = function(e) NULL)
    if (is.null(f)) next
    eff[[length(eff) + 1]] <- cbind(design = dn, f$effects); diag[[length(diag) + 1]] <- keep_diag(dn, f, "pop")
  }
}
write_res(bind_rows(eff), "its_effects"); write_res(bind_rows(diag), "its_diag")
sens <- list(
  covid_2015  = its_run(2015:2024, c(2020L, 3L), GROUPS$group, "All ACSC, 2015-2024 window"),
  covid_excl3 = its_run(2012:2024, c(2020L, 3L), GROUPS$group, "All ACSC, March-May 2020 left out", excl = 3L),
  reform_2014 = its_run(2000:2014, c(2008L, 1L), GROUPS$group, "All ACSC, 2000-2014 window"),
  reform_oct  = its_run(2000:2019, c(2008L, 10L), GROUPS$group, "All ACSC, interruption October 2008"))
write_res(bind_rows(lapply(names(sens), function(k) cbind(sensitivity = k, sens[[k]]$effects))), "its_sensitivity")


# ---- 8. Provinces ---------------------------------------------------------------

message("Provincial analysis")
prov_period <- function(p) {
  sp <- s; sp$years <- if (p == "Total") s$years else s$years[period_of(s$years) == p]
  pa <- province_analysis(a[a$year %in% sp$years, ], pop[pop$year %in% sp$years, ], sp, geo)
  mk <- moran_analysis(pa$table$asr[match(geo$province, pa$table$province)], geo, "knn4")$global   # k = 4 nearest neighbours
  list(table = cbind(period = p, pa$table), variation = cbind(period = p, pa$variation),
       moran = cbind(period = p, rbind(pa$moran, mk)))
}
pp_all <- lapply(c(names(PER), "Total"), prov_period)
write_res(bind_rows(lapply(pp_all, `[[`, "table")), "prov_period")
write_res(bind_rows(lapply(pp_all, `[[`, "variation")), "prov_variation")
write_res(bind_rows(lapply(pp_all, `[[`, "moran")), "prov_moran")
py <- province_rates(a, pop, s, by_year = TRUE); write_res(py, "province_year")
vy <- bind_rows(lapply(s$years, function(y) {
  t <- py[py$year == y, ]
  ex <- expected_counts(a[a$year == y & a$hecsap, ], pop[pop$year == y, ], "province", bands)
  scv <- (sum(((ex$observed - ex$expected) / ex$expected)^2) / nrow(ex) - sum(1 / ex$expected) / nrow(ex)) * 100
  mi <- moran_analysis(t$asr[match(geo$province, t$province)], geo, "queen", nsim = 499)$global
  data.frame(year = y, min = min(t$asr), max = max(t$asr), eq = max(t$asr) / min(t$asr),
             cv = 100 * sd(t$asr) / mean(t$asr), scv = scv, moran = mi$value, moran_p = mi$p_value)
}))
write_res(vy, "variation_by_year")
# provincial log-linear trends 2000-2011 and 2011-2019 (weighted least squares)
ptrend <- py |> group_by(province, province_name, region) |> group_modify(function(d, k) {
  f <- function(yrs) { dd <- d[d$year %in% yrs, ]; m <- lm(log(asr) ~ year, data = dd, weights = (asr / asr_se)^2)
    c(100 * (exp(coef(m)["year"]) - 1), 100 * (exp(confint(m)["year", ]) - 1)) }
  a1 <- f(2000:2011); a2 <- f(2011:2019)
  data.frame(apc_2000_2011 = a1[1], lcl1 = a1[2], ucl1 = a1[3], apc_2011_2019 = a2[1], lcl2 = a2[2], ucl2 = a2[3])
}) |> ungroup()
write_res(ptrend, "province_trend")
# occurrence-based rates, 2002-2024
so <- s; so$geo <- "occ"; so$years <- 2002:2024
write_res(province_analysis(prep_annual(dat, so), pop[pop$year %in% so$years, ], so, geo)$table, "prov_occurrence_2002_2024")


# ---- 9. Age, sex and leading diagnoses --------------------------------------------

age_sex_p <- function(p) {
  yy <- s$years[period_of(s$years) == p]; sp <- s; sp$years <- yy
  cbind(period = p, age_sex_rates(a[a$year %in% yy, ], pop[pop$year %in% yy, ], sp))
}
write_res(bind_rows(lapply(names(PER), age_sex_p)), "age_sex_period")
write_res(age_sex_rates(a, pop, s), "age_sex")
write_res(top_diagnoses(dat, pop, s), "table_top_diagnoses")
s20 <- s; s20$years <- 2020:2024; write_res(top_diagnoses(dat, pop[pop$year %in% 2020:2024, ], s20), "table_top_diagnoses_2020_2024")
s00 <- s; s00$years <- 2000:2007; write_res(top_diagnoses(dat, pop[pop$year %in% 2000:2007, ], s00), "table_top_diagnoses_2000_2007")


# ---- 10. Sensitivity analyses --------------------------------------------------------

write_res(compute_rates(h, pop, "year", LV, bands, "who"), "nat_total_who")      # WHO world standard
s10 <- s; s10$groups <- setdiff(s$groups, 10L)                                    # group 10 left out
a10 <- prep_annual(dat, s10); a10$period <- period_of(a10$year)
n10 <- compute_rates(a10[a10$hecsap, ], pop, "year", LV, bands, s$standard); write_res(n10, "nat_total_no_g10")
write_res(compute_rates(a10[a10$hecsap, ], pop, "sex", LV, bands, s$standard), "sex_no_g10")
write_res(compute_rates(h, pop, "sex", LV, bands, s$standard), "sex_total")
jp10 <- run_joinpoint(n10, "No group 10", periods = PERIODS_JP)
write_res(jp10$segments, "jp_no_g10_segments"); write_res(jp10$aapc, "jp_no_g10_aapc")
pa10 <- province_analysis(a10, pop, s10, geo); write_res(pa10$table, "prov_no_g10"); write_res(pa10$moran, "moran_no_g10")


# ---- 10b. Hospital bed supply ------------------------------------------------------------
# Available hospital beds by province of the establishment, 2004-2024 (INEC, Registro
# Estadistico de Camas Hospitalarias, series 1.5.2 and 1.5.3; national series 1.4.2 from 1997).
# Provincial series start in 2004, so the bed analyses cover 2004-2024.
# Bed counts are smoothed with a five-year running median, which removes one- and two-year
# registration spikes (e.g. Manabi 2008, Esmeraldas 2014-2015, Galapagos 2008 with zero beds);
# the reported counts are used in a sensitivity analysis.

message("Hospital bed supply")
BED_YEARS <- 2004:2024
beds <- readr::read_csv(file.path(DATA, "beds_province.csv"), show_col_types = FALSE,
                        col_types = cols(year = "i", province = "c", beds_available = "d", beds_normal = "d"))
clean_beds <- function(y, b, k = 5, tol = 1.25) {
  # Five-year running median (Tukey end rule): removes one- and two-year registration spikes and
  # keeps changes that last three years or more. `flag` marks values changed by more than 25%.
  o <- order(y); b <- b[o]
  bc <- as.numeric(stats::runmed(b, k = k, endrule = "median"))
  list(beds = bc, flag = b <= 0 | bc > tol * b | bc < b / tol)
}
beds <- beds |> group_by(province) |> arrange(year, .by_group = TRUE) |>
  mutate(cl = list(clean_beds(year, beds_available)), beds = cl[[1]]$beds, flagged = cl[[1]]$flag,
         beds_n = clean_beds(year, beds_normal)$beds) |>
  select(-cl) |> ungroup()
pop_py <- pop |> filter(province %in% PROVS) |> group_by(year, province) |> summarise(pop = sum(pop), .groups = "drop")
beds <- left_join(beds, pop_py, by = c("year", "province")) |>
  mutate(b10k = beds / pop * 1e4, b10k_raw = beds_available / pop * 1e4, b10k_normal = beds_n / pop * 1e4,
         province_name = REGION_TAB$province_name[match(province, REGION_TAB$province)],
         region = REGION_TAB$region[match(province, REGION_TAB$province)])
write_res(beds, "beds_province_year")
write_res(beds |> filter(flagged) |> select(year, province, province_name, reported = beds_available, smoothed = beds),
          "beds_flagged_values")
bn <- readr::read_csv(file.path(DATA, "beds_national.csv"), show_col_types = FALSE) |>
  mutate(b10k = beds_available / population * 1e4)
write_res(bn, "beds_national")

# Province-year counts and expected counts (indirect age-sex standardization to the
# 2004-2024 national rates of each outcome). Outcomes: ACSC, non-ACSC excluding deliveries.
province_year_panel <- function(geo_basis = "res") {
  sp <- s; sp$geo <- geo_basis; sp$years <- BED_YEARS
  aa <- rate_records(prep_annual(dat, sp))
  pp <- pop[pop$year %in% BED_YEARS, ]
  out <- list()
  for (o in c("acsc", "nonacsc")) {
    x <- if (o == "acsc") aa[aa$hecsap, ] else aa[aa$group == 0L, ]
    x$age_band <- as.character(band_of(x$age, bands)); pp$age_band <- as.character(band_of(pp$age, bands))
    ref <- left_join(summarise(group_by(x, age_band, sex), n = sum(n), .groups = "drop"),
                     summarise(group_by(pp, age_band, sex), pop = sum(pop), .groups = "drop"), by = c("age_band", "sex"))
    ref$r <- ref$n / ref$pop
    e <- pp |> group_by(year, province, age_band, sex) |> summarise(pop = sum(pop), .groups = "drop") |>
      left_join(ref[, c("age_band", "sex", "r")], by = c("age_band", "sex")) |>
      group_by(year, province) |> summarise(expected = sum(pop * r), .groups = "drop")
    ob <- x |> group_by(year, province) |> summarise(observed = sum(n), .groups = "drop")
    out[[o]] <- left_join(e, ob, by = c("year", "province")) |> mutate(observed = coalesce(observed, 0), outcome = o)
  }
  bind_rows(out) |> left_join(beds, by = c("year", "province"))
}
panel <- province_year_panel("res")
panel_occ <- province_year_panel("occ")
write_res(panel, "beds_panel")

# Fixed-effects quasi-Poisson models: log E[y_pt] = a_p + g_t + b * log(beds per 10,000_pt) + log(offset_pt).
# Offsets: expected counts (rates) or all non-delivery hospitalizations (ACSC share).
# Standard errors clustered by province. Effect reported per 10% more beds (elasticity).
fe_model <- function(d, xvar = "b10k", offset = "expected", label = "", lag = 0) {
  d <- d[order(d$province, d$year), ]
  if (lag > 0) d <- d |> group_by(province) |> mutate(x = dplyr::lag(.data[[xvar]], lag)) |> ungroup() else d$x <- d[[xvar]]
  d <- d[!is.na(d$x) & d$x > 0 & d[[offset]] > 0, ]
  m <- stats::glm(observed ~ log(x) + factor(province) + factor(year) + offset(log(d[[offset]])),
                  family = stats::quasipoisson(), data = d)
  V <- sandwich::vcovCL(m, cluster = d$province)
  b <- stats::coef(m)["log(x)"]; se <- sqrt(V["log(x)", "log(x)"]); k <- log(1.1); z <- stats::qnorm(0.975)
  data.frame(model = label, n_obs = nrow(d), estimate = 100 * (exp(k * b) - 1), lcl = 100 * (exp(k * (b - z * se)) - 1),
             ucl = 100 * (exp(k * (b + z * se)) - 1), p = 2 * stats::pnorm(-abs(b / se)), beta = unname(b), se = se)
}
acsc_p <- panel[panel$outcome == "acsc", ]; non_p <- panel[panel$outcome == "nonacsc", ]
share_p <- acsc_p |> left_join(non_p |> transmute(year, province, other = observed), by = c("year", "province")) |>
  mutate(all_nodeliv = observed + other)
occ_p <- panel_occ[panel_occ$outcome == "acsc", ]
bed_models <- bind_rows(
  fe_model(acsc_p, label = "ACSC hospitalization rate"),
  fe_model(non_p, label = "Non-ACSC hospitalization rate (excluding deliveries)"),
  fe_model(share_p, offset = "all_nodeliv", label = "ACSC share of hospitalizations"),
  fe_model(acsc_p, "b10k_raw", label = "ACSC rate, reported bed counts"),
  fe_model(acsc_p, "b10k_normal", label = "ACSC rate, beds of normal allocation"),
  fe_model(acsc_p, lag = 1, label = "ACSC rate, beds lagged one year"),
  fe_model(acsc_p[!acsc_p$year %in% 2020:2021, ], label = "ACSC rate, 2020-2021 left out"),
  fe_model(occ_p, label = "ACSC rate by province of the hospital"),
  fe_model(share_p, "b10k_raw", offset = "all_nodeliv", label = "ACSC share, reported bed counts"))
write_res(bed_models, "beds_models")

# Between-province association: provincial ASR and beds per 10,000 by period (Spearman),
# and the share of between-province variance of log ASR explained by log beds (2004-2024),
# with Moran's I before and after adjustment.
BED_PER <- list("2004–2007" = 2004:2007, "2008–2014" = 2008:2014, "2015–2019" = 2015:2019,
                "2020–2024" = 2020:2024, "2004–2024" = BED_YEARS)
bed_between <- bind_rows(lapply(names(BED_PER), function(p) {
  yy <- BED_PER[[p]]; sp <- s; sp$years <- yy
  r <- province_rates(a[a$year %in% yy, ], pop[pop$year %in% yy, ], sp)
  bb <- beds |> filter(year %in% yy) |> group_by(province) |> summarise(b10k = sum(beds) / sum(pop) * 1e4)
  d <- left_join(r[, c("province", "province_name", "asr")], bb, by = "province")
  ct <- suppressWarnings(stats::cor.test(d$asr, d$b10k, method = "spearman"))
  fit <- stats::lm(log(asr) ~ log(b10k), data = d)
  v <- d$asr[match(geo$province, d$province)]; res <- stats::residuals(fit)[match(geo$province, d$province)]
  data.frame(period = p, rho = unname(ct$estimate), rho_p = ct$p.value,
             slope_10pct = 100 * (1.1^stats::coef(fit)[2] - 1), r2 = summary(fit)$r.squared,
             moran_asr = moran_analysis(log(v), geo)$global$value, moran_asr_p = moran_analysis(log(v), geo)$global$p_value,
             moran_resid = moran_analysis(res, geo)$global$value, moran_resid_p = moran_analysis(res, geo)$global$p_value,
             b10k_min = min(d$b10k), b10k_max = max(d$b10k))
}))
write_res(bed_between, "beds_between_province")
write_res(bind_rows(lapply(names(BED_PER), function(p) {
  yy <- BED_PER[[p]]; sp <- s; sp$years <- yy
  r <- province_rates(a[a$year %in% yy, ], pop[pop$year %in% yy, ], sp)
  bb <- beds |> filter(year %in% yy) |> group_by(province) |> summarise(b10k = sum(beds) / sum(pop) * 1e4)
  cbind(period = p, left_join(r[, c("province", "province_name", "region", "asr", "asr_lcl", "asr_ucl")], bb, by = "province"))
})), "beds_asr_province_period")

# Contribution of the change in bed supply to the change in ACSC hospitalizations, 2011-2024:
# predicted national ACSC count in 2024 with provincial beds as observed versus held at 2011 levels.
bmain <- bed_models$beta[1]
ch <- acsc_p |> filter(year %in% c(2011, 2024)) |> select(year, province, b10k, observed, expected) |>
  tidyr::pivot_wider(names_from = year, values_from = c(b10k, observed, expected))
bed_contrib <- data.frame(
  beds_per_10k_2011 = sum(beds$beds[beds$year == 2011]) / sum(beds$pop[beds$year == 2011]) * 1e4,
  beds_per_10k_2024 = sum(beds$beds[beds$year == 2024]) / sum(beds$pop[beds$year == 2024]) * 1e4,
  sir_change_pct = 100 * ((sum(ch$observed_2024) / sum(ch$expected_2024)) / (sum(ch$observed_2011) / sum(ch$expected_2011)) - 1),
  attributable_pct = 100 * (sum(ch$observed_2024) / sum(ch$observed_2024 * (ch$b10k_2011 / ch$b10k_2024)^bmain) - 1))
bed_contrib$provinces_fewer_beds <- sum(ch$b10k_2024 < ch$b10k_2011)
write_res(bed_contrib, "beds_contribution_2011_2024")


# ---- 11. Figures ------------------------------------------------------------------------

message("Figures")
rd <- function(n) readr::read_csv(file.path(RES, paste0(n, ".csv")), show_col_types = FALSE)
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".tif")), p, width = w, height = h, dpi = 300, device = ragg::agg_tiff, compression = "lzw", bg = "white")
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 300, device = ragg::agg_png, bg = "white")
}
theme_set(theme_classic(base_size = 9) +
  theme(legend.position = "bottom", plot.tag = element_text(size = 11, face = "bold"), axis.line = element_line(linewidth = 0.3),
        strip.background = element_blank(), strip.text = element_text(size = 9, face = "bold")))
CAT_EN <- c("Chronic" = "Chronic conditions", "Infectious" = "Infectious conditions",
            "Vaccine-preventable" = "Vaccine-preventable conditions", "Other acute" = "Other acute conditions")
CAT_COL <- setNames(c("#1f5f99", "#d9822b", "#2e8b57", "#8064a2"), CAT_EN)
SEX_COL <- c("Total" = "#2b2b2b", "Men" = "#2c7fb8", "Women" = "#d95f0e")
REF <- list(geom_vline(xintercept = 2008, linetype = "dotted", colour = "grey45"),
            geom_vline(xintercept = 2020 + 2 / 12, linetype = "dashed", colour = "grey45"))
YRS <- scale_x_continuous(breaks = seq(2000, 2024, 4))
jpf <- rd("jp_fitted"); annf <- rd("annual_totals")

# Figure 1. A: age-standardized rates by sex with joinpoint models. B: ACSC share and hospitalization rate.
d2 <- bind_rows(mutate(rd("nat_total"), series = "Total"), mutate(rd("nat_sex"), series = ifelse(sex == 1, "Men", "Women")))
ft <- jpf |> filter(series %in% c("Total", "Hombre", "Mujer")) |> mutate(series = recode(series, Hombre = "Men", Mujer = "Women"))
d2$series <- factor(d2$series, names(SEX_COL)); ft$series <- factor(ft$series, names(SEX_COL))
pA <- ggplot(d2, aes(year, asr, colour = series)) + REF + geom_line(data = ft, aes(year, fitted), linewidth = 0.7) +
  geom_errorbar(aes(ymin = asr_lcl, ymax = asr_ucl), width = 0.25) + geom_point(size = 1.4) +
  scale_colour_manual(values = SEX_COL, name = NULL) + YRS + scale_y_continuous(labels = scales::label_comma(), limits = c(450, NA)) +
  labs(x = NULL, y = "Age-standardized rate per 100,000")
d1 <- bind_rows(annf |> transmute(year, value = share, what = "ACSC, % of all hospitalizations"),
                annf |> transmute(year, value = share_nodeliv, what = "ACSC, % of hospitalizations excluding deliveries"),
                annf |> transmute(year, value = disch_per_1000 / 4, what = "All hospitalizations per 1,000 population (right axis)"))
pB <- ggplot(d1, aes(year, value, colour = what, linetype = what)) + REF + geom_line(linewidth = 0.7) + geom_point(size = 1) +
  scale_colour_manual(values = c("#1f5f99", "#7fa7cc", "grey30"), name = NULL) +
  scale_linetype_manual(values = c("solid", "solid", "22"), name = NULL) +
  scale_y_continuous(name = "% of hospitalizations", limits = c(0, 25), sec.axis = sec_axis(~ . * 4, name = "Hospitalizations per 1,000")) +
  YRS + guides(colour = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1)) + labs(x = NULL)
save_fig(pA / pB + plot_annotation(tag_levels = "A"), "Figure1", 7, 8)

# Figure 2. Provincial age-standardized rates by period (common quintiles), High-High clusters outlined.
sf::sf_use_s2(FALSE)
g <- geo$geo; gm <- sf::st_geometry(g); i <- g$province == "20"; gm[i] <- gm[i] + c(7, 0)   # Galapagos moved 7 degrees east
sf::st_geometry(g) <- sf::st_set_crs(gm, sf::st_crs(geo$geo))
pp <- rd("prov_period") |> mutate(province = sprintf("%02d", as.integer(province)))
pp4 <- pp |> filter(period != "Total")
br <- quantile(pp4$asr, 0:5 / 5); br[1] <- br[1] - 1; br[6] <- br[6] + 1
labs_q <- paste0(formatC(round(br[-6]), format = "d", big.mark = ","), "\u2013", formatC(round(br[-1]), format = "d", big.mark = ","))
gm4 <- merge(g, pp4, by = "province")
gm4$q <- cut(gm4$asr, br, labels = labs_q, include.lowest = TRUE)
mtheme <- theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.title = element_blank(),
                legend.position = "bottom", legend.key.size = unit(3.5, "mm"))
f2 <- ggplot(gm4) + geom_sf(aes(fill = q), colour = "grey35", linewidth = 0.12) +
  geom_sf(data = gm4[gm4$cluster == "High-High", ], fill = NA, colour = "black", linewidth = 0.6) +
  facet_wrap(~ period, ncol = 2) + coord_sf(expand = FALSE) +
  scale_fill_brewer(palette = "YlOrRd", name = "ASR per 100,000") + guides(fill = guide_legend(nrow = 1)) + mtheme
save_fig(f2, "Figure2", 7.5, 7.5)

# Figure S1. Rates by condition category with joinpoint models.
fc <- jpf |> filter(series %in% names(CAT_EN)) |> mutate(cat = factor(CAT_EN[series], CAT_EN))
fs1 <- rd("nat_category") |> mutate(cat = factor(CAT_EN[category], CAT_EN)) |>
  ggplot(aes(year, asr, colour = cat)) + REF +
  geom_ribbon(aes(ymin = asr_lcl, ymax = asr_ucl, fill = cat), alpha = 0.2, colour = NA) +
  geom_point(size = 1.2) + geom_line(data = fc, aes(year, fitted), linewidth = 0.7) +
  facet_wrap(~ cat, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = CAT_COL, guide = "none") + scale_fill_manual(values = CAT_COL, guide = "none") +
  YRS + expand_limits(y = 0) + labs(x = NULL, y = "Age-standardized rate per 100,000")
save_fig(fs1, "FigureS1", 7.5, 5.5)

# Figure S2. Rates by ACSC group with joinpoint models.
ng <- rd("nat_group"); fg <- jpf |> filter(grepl("^\\d\\d ", series)) |> mutate(group = as.integer(substr(series, 1, 2)))
labf <- function(gg) factor(sprintf("%d. %s", gg, GROUPS_EN[gg]), sprintf("%d. %s", 1:20, GROUPS_EN))
ng$lab <- labf(ng$group); fg$lab <- labf(fg$group)
ng$cat <- CAT_EN[GROUPS$category[match(ng$group, GROUPS$group)]]; fg$cat <- CAT_EN[GROUPS$category[match(fg$group, GROUPS$group)]]
fs2 <- ggplot(ng, aes(year, asr, colour = cat)) + REF + geom_point(size = 0.6) + geom_line(data = fg, aes(year, fitted), linewidth = 0.5) +
  facet_wrap(~ lab, scales = "free_y", ncol = 4, labeller = label_wrap_gen(30)) +
  scale_colour_manual(values = CAT_COL, name = NULL) + scale_x_continuous(breaks = c(2000, 2012, 2024)) + expand_limits(y = 0) +
  labs(x = NULL, y = "Age-standardized rate per 100,000") +
  theme(strip.text = element_text(size = 7, face = "plain"), axis.text = element_text(size = 6.5))
save_fig(fs2, "FigureS2", 7.5, 8.75)

# Figure S3. Geographic variation over time.
vyf <- rd("variation_by_year")
d5 <- bind_rows(vyf |> transmute(year, value = cv, what = "Coefficient of variation (%)"),
                vyf |> transmute(year, value = scv, what = "Systematic component of variation (\u00d7100)"),
                vyf |> transmute(year, value = moran * 100, what = "Global Moran's I (\u00d7100)"))
fs3 <- ggplot(d5, aes(year, value)) + REF + geom_line(colour = "#1f5f99", linewidth = 0.7) + geom_point(colour = "#1f5f99", size = 1.2) +
  geom_point(data = vyf |> filter(moran_p < 0.05) |> transmute(year, value = moran * 100, what = "Global Moran's I (\u00d7100)"),
             shape = 21, size = 2.6, colour = "#c0392b", stroke = 0.7) +
  facet_wrap(~ what, ncol = 1, scales = "free_y") + YRS + expand_limits(y = 0) + labs(x = NULL, y = NULL)
save_fig(fs3, "FigureS3", 6.5, 6.5)

# Figure S4. Segmented regression around the 2008 reform and the onset of COVID-19.
its_panel <- function(f, title, ylab = "Hospitalizations per 100,000 per month", scale = 1) {
  d <- f$data; t0 <- d$date[attr(d, "t0")]; d$seg <- ifelse(d$post == 1, "post", "pre")
  ggplot(d, aes(date)) + geom_vline(xintercept = t0, colour = "grey40", linewidth = 0.4) +
    geom_point(aes(y = n_rate * scale), colour = "grey35", size = 0.6) +
    geom_line(data = d[d$post == 1, ], aes(y = trend_cf_rate * scale, linetype = "Projected without interruption"), colour = "#c0392b", linewidth = 0.6) +
    geom_line(aes(y = trend_rate * scale, group = seg, linetype = "Fitted trend (deseasonalized)"), colour = "#1f5f99", linewidth = 0.7) +
    scale_linetype_manual(values = c("Fitted trend (deseasonalized)" = "solid", "Projected without interruption" = "dashed"), name = NULL) +
    scale_x_date(date_breaks = "2 years", date_labels = "%Y") + expand_limits(y = 0) + labs(x = NULL, y = ylab, subtitle = title)
}
fs4 <- (its_panel(fits[["reform All ACSC pop"]], "2008 reform: ACSC hospitalizations") |
        its_panel(fits[["reform share"]], "2008 reform: ACSC share of hospitalizations", "ACSC, % of hospitalizations", 1e-3)) /
       (its_panel(fits[["covid All ACSC pop"]], "COVID-19: ACSC hospitalizations") |
        its_panel(fits[["covid Non-ACSC (excluding deliveries) pop"]], "COVID-19: non-ACSC hospitalizations")) +
  plot_annotation(tag_levels = "A") + plot_layout(guides = "collect") & theme(legend.position = "bottom")
save_fig(fs4, "FigureS4", 7.5, 7)

# Figure S5. Age-specific rates by sex and period.
asp <- rd("age_sex_period"); asp$age_band <- factor(asp$age_band, make_bands(top = 80)$label)
asp$sexl <- ifelse(asp$sex == 1, "Men", "Women")
PCOL <- setNames(c("#9ecae1", "#4292c6", "#08519c", "#d95f0e"), names(PER))
fs5 <- ggplot(asp, aes(age_band, rate, colour = period, group = period)) + geom_line(linewidth = 0.7) + geom_point(size = 1.1) +
  facet_wrap(~ sexl) + scale_colour_manual(values = PCOL, name = NULL) +
  scale_y_log10(labels = scales::label_comma(), breaks = c(200, 500, 1000, 2000, 5000)) +
  labs(x = "Age group (years)", y = "ACSC hospitalizations per 100,000 person-years (log scale)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7))
save_fig(fs5, "FigureS5", 7.5, 4.5)

# Figure S6. Annual provincial rates (heat map).
pyf <- rd("province_year"); tot <- pp |> filter(period == "Total")
pyf$province_name <- factor(pyf$province_name, tot$province_name[order(tot$asr)])
fs6 <- ggplot(pyf, aes(year, province_name, fill = asr)) + geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, name = "ASR per\n100,000", labels = scales::label_comma()) +
  scale_x_continuous(breaks = seq(2000, 2024, 2), expand = c(0, 0)) + labs(x = NULL, y = NULL) +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, vjust = 0.5), axis.line = element_blank())
save_fig(fs6, "FigureS6", 7.5, 6)

# Figure S7. Local Moran's I clusters by period.
LISA_COL <- c("High-High" = "#c0392b", "Low-Low" = "#2471a3", "High-Low" = "#f1948a", "Low-High" = "#85c1e9",
              "Not significant" = "#e5e5e5", "No neighbours" = "#ffffff")
gm4$cluster <- factor(gm4$cluster, names(LISA_COL))
fs7 <- ggplot(gm4) + geom_sf(aes(fill = cluster), colour = "grey35", linewidth = 0.12) + facet_wrap(~ period, ncol = 2) +
  coord_sf(expand = FALSE) + scale_fill_manual(values = LISA_COL, name = "Local Moran's I", drop = FALSE) + mtheme
save_fig(fs7, "FigureS7", 7.5, 7.5)

# Figure S8. Hospital beds per 10,000 population by province (2004 and 2024) and provincial
# ACSC rates against bed supply by period.
bp <- rd("beds_province_year") |> filter(year %in% c(2004, 2024)) |>
  mutate(province_name = factor(province_name, province_name[year == 2024][order(b10k[year == 2024])]))
pS8a <- ggplot(bp, aes(b10k, province_name)) +
  geom_line(aes(group = province_name), colour = "grey70", linewidth = 0.6) +
  geom_point(aes(colour = factor(year)), size = 1.8) +
  scale_colour_manual(values = c("2004" = "#9ecae1", "2024" = "#08519c"), name = NULL) +
  labs(x = "Beds per 10,000 population", y = NULL)
bap <- rd("beds_asr_province_period") |> filter(period != "2004\u20132024")
bap$period <- factor(bap$period, unique(bap$period))
pS8b <- ggplot(bap, aes(b10k, asr, colour = period)) + geom_point(size = 1.4) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5) +
  scale_x_log10() + scale_y_log10(labels = scales::label_comma()) +
  scale_colour_manual(values = setNames(c("#9ecae1", "#4292c6", "#08519c", "#d95f0e"), levels(bap$period)), name = NULL) +
  guides(colour = guide_legend(nrow = 2)) +
  labs(x = "Beds per 10,000 (log scale)", y = "ACSC age-standardized rate per 100,000 (log scale)")
save_fig(pS8a + pS8b + plot_annotation(tag_levels = "A") + plot_layout(widths = c(1, 1.1)), "FigureS8", 7.5, 4.8)

writeLines(capture.output(sessionInfo()), file.path(RES, "session_info.txt"))
message("Done. Results in ", RES, "/, figures in ", FIG, "/")
