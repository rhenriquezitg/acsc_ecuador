# =============================================================================
# Audit of "Hospitalizations for ambulatory care-sensitive conditions in Ecuador,
# 2000-2024" (Martinez-Abarca AL, Henriquez-Trujillo MC, Bonilla JE,
# Henriquez-Trujillo AR).
#
# Purpose: allow an independent reviewer to verify every number reported in the
# manuscript. The script
#   1. checks the integrity of the aggregated input data (data/) against totals
#      published by INEC;
#   2. re-implements each method from scratch, without the functions of
#      acsc_ecuador_analysis.R: direct age standardization and gamma confidence
#      intervals, joinpoint models, segmented (interrupted time series) regression,
#      variation indicators, Moran's I, and the fixed-effects models of hospital beds;
#   3. compares the re-implementation with the outputs of the main analysis
#      (results/), and
#   4. checks each number quoted in the manuscript (claims register).
#
# Requirements: run acsc_ecuador_analysis.R first (it writes results/). Then, from
# the repository root:   Rscript acsc_ecuador_audit.R
# Output: audit/audit_checks.csv, audit/audit_claims.csv and audit/audit_report.txt.
# On Linux or macOS use a UTF-8 locale (e.g. LANG=C.UTF-8).
#
# What is NOT re-derived here: the classification of 25 million records (done from
# the INEC microdata with the inectab package; see section 3 of the main script and
# the record log), and the joinpoint permutation tests (the selected models are
# refitted and the location of the joinpoints is checked by exhaustive search; the
# test statistics are reported from results/jp_tests.csv).
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(sandwich); library(spdep) })
DATA <- "data"; RES <- "results"; OUT <- "audit"
dir.create(OUT, showWarnings = FALSE)
PROVS <- sprintf("%02d", 1:24)
ACSC <- 1:20; DELIV <- 99L
CAT <- list(Chronic = c(2, 3, 5, 7, 8, 13, 17), Infectious = c(10, 12, 14, 15, 16),
            `Vaccine-preventable` = c(9, 18, 20), `Other acute` = c(1, 4, 6, 11, 19))
SEED <- 2026

# ---- registers ------------------------------------------------------------------
checks <- list(); claims <- list()
check <- function(area, what, value, target, tol, note = "") {
  ok <- isTRUE(all(abs(value - target) <= tol, na.rm = FALSE))
  checks[[length(checks) + 1]] <<- data.frame(area = area, check = what,
    max_abs_difference = signif(max(abs(value - target)), 3), tolerance = tol,
    status = if (ok) "PASS" else "FAIL", note = note)
}
# reported: number as printed in the manuscript; digits: decimals printed
claim <- function(section, text, reported, value, digits, source = "recomputed") {
  ok <- abs(round(value, digits) - reported) < 1e-9 || abs(value - reported) <= 0.5 * 10^-digits + 1e-9
  claims[[length(claims) + 1]] <<- data.frame(section = section, claim = text, reported = reported,
    recomputed = round(value, digits + 2), source = source, status = if (ok) "PASS" else "CHECK")
}
rd <- function(f) readr::read_csv(file.path(RES, paste0(f, ".csv")), show_col_types = FALSE)

# ---- data ---------------------------------------------------------------------------
annual <- read_csv(file.path(DATA, "acsc_annual.csv.gz"), show_col_types = FALSE,
                   col_types = cols(geo = "c", year = "i", province = "c", sex = "c", age = "i", group = "i", n = "d"))
monthly <- read_csv(file.path(DATA, "acsc_monthly.csv.gz"), show_col_types = FALSE,
                    col_types = cols(year = "i", month = "i", age = "i", group = "i", n = "d"))
icd <- read_csv(file.path(DATA, "acsc_icd.csv.gz"), show_col_types = FALSE,
                col_types = cols(year = "i", province = "c", sex = "c", age = "i", group = "i", icd10 = "c", n = "d"))
pop <- read_csv(file.path(DATA, "population.csv.gz"), show_col_types = FALSE,
                col_types = cols(year = "i", prov = "c", sex = "c", age = "i", pop = "d")) |> rename(province = prov)
std <- read_csv(file.path(DATA, "standard_population.csv"), show_col_types = FALSE)
bnat <- read_csv(file.path(DATA, "beds_national.csv"), show_col_types = FALSE)
bprov <- read_csv(file.path(DATA, "beds_province.csv"), show_col_types = FALSE, col_types = cols(province = "c"))
obst <- read_csv(file.path(DATA, "obstetric_discharges_2000_2024.csv"), show_col_types = FALSE)
covid <- read_csv(file.path(DATA, "covid_codes_2019_2024.csv"), show_col_types = FALSE)
rlog <- read_csv(file.path(DATA, "record_log.csv"), show_col_types = FALSE, na = character())
geo <- readRDS(file.path(DATA, "provinces_geo.rds"))
res <- annual[annual$geo == "res", ]
period_of <- function(y) cut(y, c(1999, 2007, 2014, 2019, 2024), labels = c("2000–2007", "2008–2014", "2015–2019", "2020–2024"))
PER <- levels(period_of(2000))

# =============================================================================
# 1. Data integrity
# =============================================================================
tot_res <- res |> group_by(year) |> summarise(n = sum(n))
tot_occ <- annual[annual$geo == "occ", ] |> group_by(year) |> summarise(n = sum(n))
tot_mon <- monthly |> group_by(year) |> summarise(n = sum(n))
check("Data", "Annual discharges equal INEC published totals (series 1.4.2)", tot_res$n, bnat$discharges[match(tot_res$year, bnat$year)], 1,
      "The 2014 microdata file has 1,192,750 records; INEC publishes 1,192,749")
check("Data", "Residence and occurrence files have the same totals", tot_res$n, tot_occ$n, 0)
check("Data", "Monthly file adds up to the annual file", tot_mon$n, tot_res$n, 0)
check("Data", "Obstetric file totals equal the annual file", obst$all_discharges, tot_res$n, 0)
acsc_res <- res |> filter(group %in% ACSC) |> group_by(year) |> summarise(n = sum(n))
acsc_icd <- icd |> group_by(year) |> summarise(n = sum(n))
check("Data", "ICD-10 detail file adds up to ACSC totals", acsc_icd$n, acsc_res$n, 0)
pnat <- pop |> group_by(year) |> summarise(pop = sum(pop))
check("Data", "Population equals INEC projections used in series 1.4.2", pnat$pop, bnat$population[match(pnat$year, bnat$year)], 1)
bsum <- bprov |> group_by(year) |> summarise(b = sum(beds_available))
ua <- read_csv(file.path(DATA, "beds_unassigned_2000_2003.csv"), show_col_types = FALSE)
bsum$b[bsum$year %in% ua$year] <- bsum$b[bsum$year %in% ua$year] + ua$beds_available[match(bsum$year[bsum$year %in% ua$year], ua$year)]
check("Data", "Provincial beds (plus unassigned) equal national available beds", bsum$b, bnat$beds_available[match(bsum$year, bnat$year)], 3,
      "Differences of 1-2 beds come from rounding the 2004 shares used for Santo Domingo and Santa Elena in 2000-2003")
check("Data", "Every ACSC group 1-20 present every year", nrow(distinct(res[res$group %in% ACSC, c("year", "group")])), 25 * 20, 0)

# =============================================================================
# 2. Independent re-implementation of the methods
# =============================================================================
# 2.1 Direct age standardization, 18 five-year bands (85+), Ecuador 2022 standard,
#     gamma confidence intervals (Fay and Feuer 1997)
band <- function(age) pmin(age %/% 5, 17)
w_std <- function(standard = "ecuador_2022") {
  s <- std[std$standard == standard, ]; w <- tapply(s$weight, band(s$age), sum); w / sum(w)
}
W <- w_std()
usable <- function(d) d[!is.na(d$age) & d$sex %in% c("1", "2") & d$province %in% PROVS, ]
asr_of <- function(cnt, pp, w = W) {
  n <- tapply(cnt$n, factor(band(cnt$age), 0:17), sum); n[is.na(n)] <- 0
  P <- tapply(pp$pop, factor(band(pp$age), 0:17), sum)
  r <- n / P; asr <- sum(w * r) * 1e5
  v <- sum((w / P)^2 * n) * 1e10; wm <- max(w / P) * 1e5
  c(asr = asr, lcl = if (v > 0) qgamma(0.025, asr^2 / v, scale = v / asr) else 0,
    ucl = qgamma(0.975, (asr + wm)^2 / (v + wm^2), scale = (v + wm^2) / (asr + wm)), se = sqrt(v), n = sum(n))
}
U <- usable(res)
series_asr <- function(sel, years = 2000:2024, extra = NULL) {
  t(sapply(years, function(y) asr_of(U[U$year == y & sel(U), ], if (is.null(extra)) pop[pop$year == y, ] else extra(pop[pop$year == y, ]))))
}
A_tot <- series_asr(function(d) d$group %in% ACSC)
nt <- rd("nat_total")
check("Rates", "National ASR 2000-2024 (independent)", A_tot[, "asr"], nt$asr, 1e-6)
check("Rates", "Gamma 95% CI of national ASR", c(A_tot[, "lcl"], A_tot[, "ucl"]), c(nt$asr_lcl, nt$asr_ucl), 1e-6)
A_men <- series_asr(function(d) d$group %in% ACSC & d$sex == "1", extra = function(p) p[p$sex == "1", ])
A_wom <- series_asr(function(d) d$group %in% ACSC & d$sex == "2", extra = function(p) p[p$sex == "2", ])
ns <- rd("nat_sex")
check("Rates", "ASR by sex", c(A_men[, "asr"], A_wom[, "asr"]), c(ns$asr[ns$sex == 1], ns$asr[ns$sex == 2]), 1e-6)
A_cat <- lapply(CAT, function(g) series_asr(function(d) d$group %in% g))
nc <- rd("nat_category")
check("Rates", "ASR by condition category", unlist(lapply(names(CAT), function(k) A_cat[[k]][, "asr"])),
      unlist(lapply(names(CAT), function(k) nc$asr[nc$category == k])), 1e-6)
A_grp <- lapply(ACSC, function(g) series_asr(function(d) d$group == g))
ng <- rd("nat_group")
check("Rates", "ASR by ACSC group (20 groups)", unlist(lapply(ACSC, function(g) A_grp[[g]][, "asr"])),
      unlist(lapply(ACSC, function(g) ng$asr[ng$group == g])), 1e-6)
W_who <- w_std("who")
A_who <- t(sapply(2000:2024, function(y) asr_of(U[U$year == y & U$group %in% ACSC, ], pop[pop$year == y, ], W_who)))
check("Rates", "ASR with WHO standard (sensitivity)", A_who[, "asr"], rd("nat_total_who")$asr, 1e-6)

# annual totals and shares (all discharges as denominator)
ann <- res |> group_by(year) |> summarise(acsc = sum(n[group %in% ACSC]), all = sum(n), all_res = sum(n[province %in% PROVS])) |>
  left_join(pop |> group_by(year) |> summarise(pop = sum(pop)), by = "year") |>
  mutate(share = 100 * acsc / all, per1000 = 1000 * all_res / pop)
at <- rd("annual_totals")
check("Rates", "Annual ACSC share of all discharges", ann$share, at$share, 1e-9)
check("Rates", "Annual discharges per 1,000 population", ann$per1000, at$disch_per_1000, 1e-9)

# 2.2 Joinpoint models: weighted log-linear regression, weights (rate/SE)^2; the
#     selected models are refitted with the joinpoints of results/jp_segments.csv;
#     df = n - parameters - joinpoints (NCI convention)
jp_refit <- function(year, rate, se, knots) {
  y <- log(rate); w <- (rate / se)^2
  X <- cbind(1, year); for (k in knots) X <- cbind(X, pmax(year - k, 0))
  f <- lm.wfit(X, y, w); b <- f$coefficients
  df <- length(y) - ncol(X) - length(knots); s2 <- sum(w * f$residuals^2) / df
  V <- s2 * solve(crossprod(X * sqrt(w))); tq <- qt(0.975, df)
  starts <- c(min(year), knots); ends <- c(knots, max(year)); out <- NULL
  for (j in seq_along(starts)) {
    l <- c(0, 1, rep(1, j - 1), rep(0, length(knots) - j + 1))
    sl <- sum(l * b); se_ <- sqrt(drop(t(l) %*% V %*% l))
    out <- rbind(out, data.frame(start = starts[j], end = ends[j], apc = 100 * (exp(sl) - 1),
                                 lcl = 100 * (exp(sl - tq * se_) - 1), ucl = 100 * (exp(sl + tq * se_) - 1), p = 2 * pt(-abs(sl / se_), df)))
  }
  wl <- ends - starts; L <- sapply(seq_along(starts), function(j) c(0, 1, rep(1, j - 1), rep(0, length(knots) - j + 1)))
  l <- drop(L %*% wl) / sum(wl); a <- sum(l * b); sa <- sqrt(drop(t(l) %*% V %*% l))
  attr(out, "aapc") <- c(aapc = 100 * (exp(a) - 1), lcl = 100 * (exp(a - tq * sa) - 1), ucl = 100 * (exp(a + tq * sa) - 1), p = 2 * pt(-abs(a / sa), df))
  attr(out, "rss") <- sum(w * f$residuals^2)
  out
}
seg <- rd("jp_segments"); aap <- rd("jp_aapc")
GLAB <- unique(seg$series[grepl("^\\d\\d ", seg$series)])
SER <- c(list(Total = A_tot, Hombre = A_men, Mujer = A_wom), A_cat, setNames(A_grp, GLAB[order(as.integer(substr(GLAB, 1, 2)))]))
sh_se <- 100 * sqrt((ann$acsc / ann$all) * (1 - ann$acsc / ann$all) / ann$all)
diffs <- c(); diffs_a <- c(); JP <- list()
for (k in c(names(SER), "Share of hospitalizations")) {
  sg <- seg[seg$series == k, ]; knots <- sg$start[-1]
  if (k == "Share of hospitalizations") { r <- ann$share; s_ <- sh_se } else { r <- SER[[k]][, "asr"]; s_ <- SER[[k]][, "se"] }
  ok <- r > 0 & s_ > 0
  f <- jp_refit((2000:2024)[ok], r[ok], s_[ok], knots); JP[[k]] <- f
  diffs <- c(diffs, as.matrix(f[, c("apc", "lcl", "ucl")]) - as.matrix(sg[, c("apc", "apc_lcl", "apc_ucl")]))
  a0 <- aap[aap$series == k & aap$period == "2000-2024", ]
  diffs_a <- c(diffs_a, attr(f, "aapc")[1:3] - unlist(a0[, c("aapc", "aapc_lcl", "aapc_ucl")]))
}
check("Joinpoint", "APC and 95% CI of every segment, all series (refit)", diffs, 0, 1e-6)
check("Joinpoint", "AAPC 2000-2024 and 95% CI, all series (refit)", diffs_a, 0, 1e-6)
# location of the joinpoint of the national series: exhaustive search, >= 2 points from each end
yy <- 2000:2024; rss1 <- sapply(2002:2022, function(k) attr(jp_refit(yy, A_tot[, "asr"], A_tot[, "se"], k), "rss"))
check("Joinpoint", "Best single joinpoint of the national ASR is 2011", (2002:2022)[which.min(rss1)], 2011, 0)
jt <- rd("jp_tests")
check("Joinpoint", "Permutation tests (results/jp_tests.csv) select 1 joinpoint for the national series",
      sum(jt$series == "Total" & jt$decision == "reject H0"), 1, 0, "Sequential tests at Bonferroni-corrected levels; not re-run here")

# 2.3 Segmented regression (ITS): quasi-Poisson, 2 Fourier pairs, Newey-West (lag 12)
mpop <- function(years) {
  py <- pnat[pnat$year %in% years, ]; x <- py$year + 0.5
  g <- expand.grid(month = 1:12, year = years); tt <- g$year + (g$month - 0.5) / 12
  p <- approx(x, py$pop, xout = tt, rule = 2)$y
  s1 <- diff(py$pop[1:2]); s2 <- diff(tail(py$pop, 2))
  p[tt < min(x)] <- py$pop[1] + s1 * (tt[tt < min(x)] - x[1]); p[tt > max(x)] <- tail(py$pop, 1) + s2 * (tt[tt > max(x)] - max(x))
  g$off <- p; g
}
mser <- function(years, grp) monthly |> filter(year %in% years, group %in% grp) |> group_by(year, month) |> summarise(n = sum(n), .groups = "drop")
its_fit <- function(years, start, grp, offset = "pop", excl = 0) {
  d <- if (offset == "pop") mpop(years) else rename(mser(years, c(0L, ACSC, DELIV)), off = n)
  d <- left_join(d, mser(years, grp), by = c("year", "month")); d$n[is.na(d$n)] <- 0
  d <- d[order(d$year, d$month), ]; d$t <- seq_len(nrow(d)); t0 <- which(d$year == start[1] & d$month == start[2])
  d$post <- as.integer(d$t >= t0); d$since <- pmax(0, d$t - t0)
  for (k in 1:2) { d[[paste0("s", k)]] <- sin(2 * pi * k * d$month / 12); d[[paste0("c", k)]] <- cos(2 * pi * k * d$month / 12) }
  fit_d <- d[!(d$t >= t0 & d$t < t0 + excl), ]
  m <- glm(n ~ t + post + since + s1 + c1 + s2 + c2 + offset(log(off)), family = quasipoisson(), data = fit_d)
  V <- NeweyWest(m, lag = 12, prewhite = FALSE, adjust = TRUE); b <- coef(m); z <- qnorm(0.975)
  eff <- function(l) { e <- sum(l * b[names(l)]); s <- sqrt(drop(t(l) %*% V[names(l), names(l)] %*% l)); 100 * (exp(c(e, e - z * s, e + z * s)) - 1) }
  cf <- d; cf$post <- 0; cf$since <- 0; post <- d$post == 1
  list(pre = eff(c(t = 12)), level = eff(c(post = 1)), change = eff(c(since = 12)), after = eff(c(t = 12, since = 12)),
       obs = sum(d$n[post]), exp = sum(predict(m, newdata = cf, type = "response")[post]))
}
ITS <- list(
  reform = list(years = 2000:2019, start = c(2008L, 1L)), covid = list(years = 2012:2024, start = c(2020L, 3L)))
SERIES <- c(list(`All ACSC` = ACSC), lapply(CAT, identity), list(`Non-ACSC` = c(0L, DELIV)))
e <- rd("its_effects"); dg <- rd("its_diag"); dd <- c(); oe <- c(); FITS <- list()
for (dn in names(ITS)) {
  for (sn in c(names(SERIES), "ACSC share of hospitalizations")) {
    f <- if (sn == "ACSC share of hospitalizations") its_fit(ITS[[dn]]$years, ITS[[dn]]$start, ACSC, "hosp") else its_fit(ITS[[dn]]$years, ITS[[dn]]$start, SERIES[[sn]])
    FITS[[paste(dn, sn)]] <- f
    r <- e[e$design == dn & e$series == sn, ]
    mine <- rbind(f$pre, f$level, f$change, f$after)
    theirs <- as.matrix(r[match(c("Pre-interruption trend (% per year)", "Immediate level change (%)",
                                  "Change in trend after interruption (% per year)", "Post-interruption trend (% per year)"), r$term), c("estimate", "lcl", "ucl")])
    dd <- c(dd, mine - theirs)
    g <- dg[dg$design == dn & dg$series == sn, ]; oe <- c(oe, f$obs - g$obs, f$exp - g$exp)
  }
}
check("Segmented regression", "Trend, level and slope changes with 95% CI, 2008 reform and COVID-19 (refit)", dd, 0, 1e-6)
check("Segmented regression", "Observed and counterfactual hospitalizations after each interruption", oe, 0, 1e-3)
sens <- rd("its_sensitivity")
S12 <- list(acsc = its_fit(2000:2019, c(2012L, 1L), ACSC), non = its_fit(2000:2019, c(2012L, 1L), c(0L, DELIV)),
            share = its_fit(2000:2019, c(2012L, 1L), ACSC, "hosp"), oct = its_fit(2000:2019, c(2008L, 10L), ACSC))
ct <- function(sn) unlist(sens[sens$sensitivity == sn & sens$term == "Change in trend after interruption (% per year)", c("estimate", "lcl", "ucl")])
check("Segmented regression", "Sensitivity: interruptions in January 2012 and October 2008 (refit)",
      c(S12$acsc$change - ct("reform_2012"), S12$non$change - ct("reform_2012_non"), S12$share$change - ct("reform_2012_share"), S12$oct$change - ct("reform_oct")), 0, 1e-6)

# 2.4 Provinces: ASR by period, CV, extremal quotient, SCV, Moran's I
prov_asr <- function(years) sapply(PROVS, function(p) asr_of(U[U$year %in% years & U$province == p & U$group %in% ACSC, ],
                                                            pop[pop$year %in% years & pop$province == p, ])["asr"])
scv_of <- function(years) {
  h <- U[U$year %in% years & U$group %in% ACSC, ]; pp <- pop[pop$year %in% years & pop$province %in% PROVS, ]
  h$b <- band(h$age); pp$b <- band(pp$age)
  r <- merge(aggregate(n ~ b + sex, h, sum), aggregate(pop ~ b + sex, pp, sum)); r$rate <- r$n / r$pop
  pe <- merge(aggregate(pop ~ province + b + sex, pp, sum), r[, c("b", "sex", "rate")])
  E <- tapply(pe$pop * pe$rate, pe$province, sum)[PROVS]; O <- tapply(h$n, h$province, sum)[PROVS]
  (mean(((O - E) / E)^2) - mean(1 / E)) * 100
}
moran_I <- function(x, nb) {
  # row-standardized weights; as in spdep (adjust.n = TRUE), n excludes units without neighbors
  n <- length(x); z <- x - mean(x); Wm <- matrix(0, n, n)
  for (i in seq_len(n)) if (length(nb[[i]]) && nb[[i]][1] != 0) Wm[i, nb[[i]]] <- 1 / length(nb[[i]])
  n_eff <- sum(rowSums(Wm) > 0)
  (n_eff / sum(Wm)) * drop(t(z) %*% Wm %*% z) / sum(z^2)
}
pv <- rd("prov_variation"); mo <- rd("prov_moran"); ppd <- rd("prov_period") |> mutate(province = sprintf("%02d", as.integer(province)))
PA <- list(); pdiff <- c(); vdiff <- c(); mdiff <- c()
for (p in c(PER, "Total")) {
  yrs <- if (p == "Total") 2000:2024 else (2000:2024)[period_of(2000:2024) == p]
  a <- prov_asr(yrs); names(a) <- PROVS; PA[[p]] <- a
  pdiff <- c(pdiff, a - ppd$asr[ppd$period == p][match(PROVS, ppd$province[ppd$period == p])])
  vi <- c(cv = 100 * sd(a) / mean(a), eq = max(a) / min(a), scv = scv_of(yrs))
  v0 <- pv[pv$period == p, ]
  vdiff <- c(vdiff, vi - v0$value[match(c("Coefficient of variation (%)", "Extremal quotient (max/min)", "Systematic component of variation (x100)"), v0$indicator)])
  x <- a[match(geo$province, PROVS)]
  mdiff <- c(mdiff, moran_I(x, geo$nb$queen) - mo$value[mo$period == p & mo$neighbours == "queen"],
             moran_I(x, geo$nb$knn4) - mo$value[mo$period == p & mo$neighbours == "knn4"])
}
check("Provinces", "Provincial ASR by period (independent)", pdiff, 0, 1e-6)
check("Provinces", "CV, extremal quotient and SCV by period (independent)", vdiff, 0, 1e-6)
check("Provinces", "Global Moran's I, queen and 4-nearest neighbors (computed by hand)", mdiff, 0, 1e-9)
# permutation p-values, empirical Bayes Moran's I and LISA reproduced with the same seed
lw <- nb2listw(geo$nb$queen, style = "W", zero.policy = TRUE); pdiff2 <- c(); lisa <- list()
for (p in c(PER, "Total")) {
  x <- PA[[p]][match(geo$province, PROVS)]
  set.seed(SEED); gm <- moran.mc(x, lw, nsim = 999, zero.policy = TRUE)
  set.seed(SEED); lp <- localmoran_perm(x, lw, nsim = 9999, zero.policy = TRUE)
  pp_ <- lp[, "Pr(z != E(Ii)) Sim"]; pp_[card(geo$nb$queen) == 0] <- NA
  lisa[[p]] <- data.frame(period = p, province = geo$province, quad = as.character(attr(lp, "quadr")$mean), p = pp_, fdr = p.adjust(pp_, "BH"))
  pdiff2 <- c(pdiff2, gm$p.value - mo$p_value[mo$period == p & mo$neighbours == "queen"],
              na.omit(pp_ - ppd$p[ppd$period == p][match(geo$province, ppd$province[ppd$period == p])]))
}
check("Provinces", "Moran's I permutation p-values and LISA p-values (same seed)", pdiff2, 0, 5e-4,
      "Provinces with one or two neighbors produce permutations that tie with the observed value; last-digit differences in the rates can move up to 4 of 9,999 ties")
lisa <- bind_rows(lisa)
# provincial log-linear trends (weighted least squares on annual ASR)
PY <- sapply(PROVS, function(p) sapply(2000:2024, function(y) asr_of(U[U$year == y & U$province == p & U$group %in% ACSC, ], pop[pop$year == y & pop$province == p, ])[c("asr", "se")]), simplify = FALSE)
trend <- function(p, yrs) {
  a <- PY[[p]]["asr", ]; s <- PY[[p]]["se", ]; i <- (2000:2024) %in% yrs
  m <- lm(log(a[i]) ~ yrs, weights = (a[i] / s[i])^2); ci <- confint(m)[2, ]; 100 * (exp(c(coef(m)[2], ci)) - 1)
}
TR <- t(sapply(PROVS, function(p) c(trend(p, 2000:2011), trend(p, 2011:2019))))
ptr <- rd("province_trend") |> mutate(province = sprintf("%02d", as.integer(province)))
check("Provinces", "Provincial trends 2000-2011 and 2011-2019 (independent)",
      as.vector(TR) - as.vector(as.matrix(ptr[match(PROVS, ptr$province), c("apc_2000_2011", "lcl1", "ucl1", "apc_2011_2019", "lcl2", "ucl2")])), 0, 1e-6)

# 2.5 Hospital beds: smoothing, fixed-effects quasi-Poisson models, cluster-robust SE
smooth <- function(b) { x <- b; repeat { nx <- as.numeric(runmed(x, 5, endrule = "median")); if (all(nx == x)) break; x <- nx }; x }
bp <- bprov |> arrange(province, year) |> group_by(province) |> mutate(bs = smooth(beds_available), bn = smooth(beds_normal)) |> ungroup()
bp <- bp |> left_join(pop |> group_by(year, province) |> summarise(pop = sum(pop), .groups = "drop"), by = c("year", "province")) |>
  mutate(b10k = bs / pop * 1e4, b10k_raw = beds_available / pop * 1e4, flag = beds_available <= 0 | bs > 1.25 * beds_available | bs < beds_available / 1.25)
bpy <- rd("beds_province_year") |> mutate(province = sprintf("%02d", as.integer(province)))
check("Beds", "Smoothed bed density by province and year (independent)", bp$b10k, bpy$b10k[match(paste(bp$year, bp$province), paste(bpy$year, bpy$province))], 1e-9)
panel_of <- function(grp, basis = "res") {
  x <- annual[annual$geo == basis, ]; x <- usable(x[x$group %in% grp, ]); x$b <- band(x$age)
  pp <- pop[pop$province %in% PROVS, ]; pp$b <- band(pp$age)
  r <- merge(aggregate(n ~ b + sex, x, sum), aggregate(pop ~ b + sex, pp, sum)); r$rate <- r$n / r$pop
  pe <- merge(aggregate(pop ~ year + province + b + sex, pp, sum), r[, c("b", "sex", "rate")])
  E <- aggregate(cbind(expected = pop * rate) ~ year + province, pe, sum)
  O <- aggregate(cbind(observed = n) ~ year + province, x, sum)
  merge(E, O, all.x = TRUE) |> mutate(observed = coalesce(observed, 0)) |> left_join(bp, by = c("year", "province"))
}
cr_glm <- function(d, xv = "b10k", offset = "expected", lag = 0) {
  d <- d[order(d$province, d$year), ]
  d$x <- if (lag) ave(d[[xv]], d$province, FUN = function(v) c(rep(NA, lag), head(v, -lag))) else d[[xv]]
  d <- d[!is.na(d$x) & d$x > 0 & d[[offset]] > 0, ]
  m <- glm(observed ~ log(x) + factor(province) + factor(year) + offset(log(d[[offset]])), family = quasipoisson(), data = d)
  X <- model.matrix(m); mu <- fitted(m); u <- X * (d$observed - mu)
  B <- solve(crossprod(X * sqrt(mu))); G <- length(unique(d$province)); N <- nrow(X); K <- ncol(X)
  meat <- crossprod(rowsum(u, d$province)); V <- (G / (G - 1)) * B %*% meat %*% B   # HC0 with cluster adjustment G/(G-1), as sandwich::vcovCL for glm
  b <- coef(m)["log(x)"]; se <- sqrt(V["log(x)", "log(x)"]); k <- log(1.1)
  c(est = 100 * (exp(k * b) - 1), lcl = 100 * (exp(k * (b - 1.96 * se)) - 1), ucl = 100 * (exp(k * (b + 1.96 * se)) - 1), beta = unname(b))
}
pa <- panel_of(ACSC); pn <- panel_of(c(0L, DELIV)); po <- panel_of(ACSC, "occ")
psh <- pa |> left_join(pn |> select(year, province, other = observed), by = c("year", "province")) |> mutate(all_hosp = observed + other)
BM <- rbind(cr_glm(pa), cr_glm(pn), cr_glm(psh, offset = "all_hosp"), cr_glm(pa, "b10k_raw"),
            cr_glm(pa |> mutate(b10k_normal = bn / pop * 1e4), "b10k_normal"), cr_glm(pa, lag = 1),
            cr_glm(pa[!pa$year %in% 2020:2021, ]), cr_glm(pa[pa$year >= 2004, ]),
            cr_glm(pa[!(pa$year <= 2003 & pa$province %in% c("09", "17", "23", "24")), ]), cr_glm(po), cr_glm(psh, "b10k_raw", "all_hosp"))
bm <- rd("beds_models")
check("Beds", "Fixed-effects models: % change per 10% more beds with 95% CI, 11 models (independent cluster-robust sandwich)",
      as.vector(BM[, 1:3]) - as.vector(as.matrix(bm[, c("estimate", "lcl", "ucl")])), 0, 0.005,
      "Cluster-robust variance computed by hand (Liang-Zeger, cluster adjustment G/(G-1)); tolerance 0.005 percentage points")
ch <- pa |> filter(year %in% c(2011, 2024)) |> select(year, province, b10k, observed, expected) |> pivot_wider(names_from = year, values_from = c(b10k, observed, expected))
contrib <- 100 * (sum(ch$observed_2024) / sum(ch$observed_2024 * (ch$b10k_2011 / ch$b10k_2024)^BM[1, "beta"]) - 1)
sir_ch <- 100 * ((sum(ch$observed_2024) / sum(ch$expected_2024)) / (sum(ch$observed_2011) / sum(ch$expected_2011)) - 1)
bc <- rd("beds_contribution_2011_2024")
check("Beds", "Contribution of bed changes 2011-2024 and change in standardized ratio", c(contrib, sir_ch), c(bc$attributable_pct, bc$sir_change_pct), 1e-6)

# =============================================================================
# 3. Claims register: numbers quoted in the manuscript
# =============================================================================
S <- "Abstract/Results"
claim("Methods", "Total discharges 2000-2024", 25074262, sum(ann$all), 0)
claim("Methods", "Records with partial three-character codes", 237, sum(as.numeric(rlog$n[grepl("partly", rlog$issue)])), 0)
claim("Methods", "Records with missing age", 7337, sum(as.numeric(rlog$n[rlog$issue == "missing or unknown age"])), 0)
claim("Methods", "Records with missing or unknown sex", 12, sum(as.numeric(rlog$n[grepl("^sex", rlog$issue)])), 0)
claim("Methods", "Residents abroad or in non-delimited zones", 39270, sum(as.numeric(rlog$n[grepl("^province", rlog$issue)])), 0)
claim("Methods", "Missing age, % of discharges", 0.03, 100 * 7337 / sum(ann$all), 2)
claim("Methods", "Residents abroad/ZND, % of discharges", 0.16, 100 * 39270 / sum(ann$all), 2)
claim("Methods", "O80-O84 as % of obstetric discharges, 2000", 71, 100 * obst$deliveries_O80_O84[1] / obst$obstetric_O00_O99[1], 0)
claim("Methods", "O80-O84 as % of obstetric discharges, 2024", 36, 100 * obst$deliveries_O80_O84[25] / obst$obstetric_O00_O99[25], 0)
claim("Methods", "Bed values changed >25% by smoothing (of 600)", 21, sum(bp$flag), 0)
claim("Introduction", "Discharges 2000 (millions)", 0.67, ann$all[1] / 1e6, 2)
claim("Introduction", "Discharges 2011 (millions)", 1.13, ann$all[12] / 1e6, 2)
claim("Introduction", "Population 2000 (millions)", 12.5, bnat$population[bnat$year == 2000] / 1e6, 1)
claim("Introduction", "Population 2024 (millions)", 18.0, bnat$population[bnat$year == 2024] / 1e6, 1)
claim(S, "ACSC hospitalizations", 3792560, sum(ann$acsc), 0)
claim(S, "ACSC share of all discharges, 2000-2024 (%)", 15.1, 100 * sum(ann$acsc) / sum(ann$all), 1)
claim(S, "Discharges per 1,000, 2000", 53.7, ann$per1000[1], 1)
claim(S, "Discharges per 1,000, 2011", 74.6, ann$per1000[12], 1)
claim(S, "Discharges per 1,000, 2024", 63.0, ann$per1000[25], 1)
claim(S, "ACSC share 2000 (%)", 17.0, ann$share[1], 1)
claim(S, "ACSC share 2011 (%)", 14.9, ann$share[12], 1)
claim(S, "ACSC share 2024 (%)", 14.4, ann$share[25], 1)
hc <- res[res$group %in% ACSC, ]; hc$cat <- names(CAT)[sapply(hc$group, function(g) which(sapply(CAT, function(v) g %in% v)))]
claim(S, "Infectious conditions, % of ACSC", 53.5, 100 * sum(hc$n[hc$cat == "Infectious"]) / sum(hc$n), 1)
claim(S, "Chronic conditions, % of ACSC", 30.5, 100 * sum(hc$n[hc$cat == "Chronic"]) / sum(hc$n), 1)
gp <- function(g) 100 * sum(hc$n[hc$group == g]) / sum(hc$n)
claim(S, "Infectious gastroenteritis (group 12), % of ACSC", 26.1, gp(12), 1)
claim(S, "Kidney and urinary tract infections (group 16), % of ACSC", 11.9, gp(16), 1)
claim(S, "Diabetes (group 5), % of ACSC", 9.4, gp(5), 1)
claim(S, "Prenatal care conditions (group 10), % of ACSC", 6.7, gp(10), 1)
hc$per <- period_of(hc$year)
ag <- function(per, lo, hi) { x <- hc[hc$per == per, ]; 100 * sum(x$n[!is.na(x$age) & x$age >= lo & x$age <= hi]) / sum(x$n) }
claim(S, "Children under five, % of ACSC, 2000-2007", 23.7, ag(PER[1], 0, 4), 1)
claim(S, "Children under five, % of ACSC, 2020-2024", 17.3, ag(PER[4], 0, 4), 1)
claim(S, "Aged 65+, % of ACSC, 2000-2007", 19.8, ag(PER[1], 65, 100), 1)
claim(S, "Aged 65+, % of ACSC, 2020-2024", 28.6, ag(PER[4], 65, 100), 1)
claim(S, "ASR 2000", 1016.6, A_tot[1, "asr"], 1); claim(S, "ASR 2010 (peak)", 1195.9, max(A_tot[, "asr"]), 1)
claim(S, "Year of peak ASR", 2010, (2000:2024)[which.max(A_tot[, "asr"])], 0)
claim(S, "ASR 2024", 897.2, A_tot[25, "asr"], 1); claim(S, "ASR 2024 below peak (%)", 25.0, 100 * (1 - A_tot[25, "asr"] / max(A_tot[, "asr"])), 1)
j <- JP$Total
claim(S, "APC 2000-2011", 2.0, j$apc[1], 1); claim(S, "APC 2000-2011 lower CI", 0.4, j$lcl[1], 1); claim(S, "APC 2000-2011 upper CI", 3.7, j$ucl[1], 1)
claim(S, "APC 2000-2011 p-value", 0.018, j$p[1], 3)
claim(S, "APC 2011-2024", -3.1, j$apc[2], 1); claim(S, "APC 2011-2024 lower CI", -4.3, j$lcl[2], 1); claim(S, "APC 2011-2024 upper CI", -1.9, j$ucl[2], 1)
claim(S, "AAPC 2000-2024", -0.8, attr(j, "aapc")[1], 1); claim(S, "AAPC lower CI", -1.4, attr(j, "aapc")[2], 1); claim(S, "AAPC upper CI", -0.1, attr(j, "aapc")[3], 1)
claim(S, "Joinpoint in men", 2010, JP$Hombre$start[2], 0)
js <- JP$`Share of hospitalizations`
claim(S, "Share: number of joinpoints", 0, nrow(js) - 1, 0)
claim(S, "Share APC 2000-2024", -1.1, js$apc[1], 1); claim(S, "Share APC lower CI", -1.4, js$lcl[1], 1); claim(S, "Share APC upper CI", -0.7, js$ucl[1], 1)
claim(S, "Vaccine-preventable APC until 2010", -11.7, JP$`Vaccine-preventable`$apc[1], 1)
claim(S, "Vaccine-preventable first joinpoint", 2010, JP$`Vaccine-preventable`$end[1], 0)
g12 <- JP[[GLAB[as.integer(substr(GLAB, 1, 2)) == 12]]]
claim(S, "Gastroenteritis: joinpoint year", 2007, g12$start[2], 0); claim(S, "Gastroenteritis APC after 2007", -5.9, g12$apc[2], 1)
claim(S, "Diabetes: decline starts", 2011, JP[[GLAB[as.integer(substr(GLAB, 1, 2)) == 5]]]$start[2], 0)
claim(S, "Urinary tract infections: decline starts", 2013, JP[[GLAB[as.integer(substr(GLAB, 1, 2)) == 16]]]$start[2], 0)
j159 <- icd |> filter(icd10 == "J159") |> group_by(year) |> summarise(n = sum(n))
claim(S, "J15.9 hospitalizations, 2000", 106, j159$n[j159$year == 2000], 0)
claim(S, "J15.9 hospitalizations, 2018 (between 9,500 and 10,000)", 1, as.numeric(j159$n[j159$year == 2018] >= 9500 & j159$n[j159$year == 2018] <= 10000), 0)
claim(S, "J15.9 hospitalizations, 2019 (between 9,500 and 10,000)", 1, as.numeric(j159$n[j159$year == 2019] >= 9500 & j159$n[j159$year == 2019] <= 10000), 0)
ra <- function(sx, lo) { x <- U[U$group %in% ACSC & U$sex == sx & U$age >= lo & U$age <= lo + 4, ]; sum(x$n) / sum(pop$pop[pop$sex == sx & pop$age >= lo & pop$age <= lo + 4]) }
rr <- sapply(c(15, 20, 25, 30), function(lo) ra("2", lo) / ra("1", lo))
claim(S, "Women/men rate ratio 15-34 years, minimum", 2.3, min(rr), 1); claim(S, "Women/men rate ratio 15-34 years, maximum", 3.1, max(rr), 1)
# provinces
P <- "Provinces"; a <- PA$Total
claim(P, "Lowest provincial ASR 2000-2024 (Orellana)", 800.3, a["22"], 1); claim(P, "Orellana is the lowest", 1, as.numeric(names(which.min(a)) == "22"), 0)
claim(P, "Highest provincial ASR 2000-2024 (Morona Santiago)", 2162.4, a["14"], 1)
claim(P, "Five highest: Morona Santiago, Pastaza, Napo, Zamora Chinchipe, Canar", 1, as.numeric(identical(names(sort(a, decreasing = TRUE))[1:5], c("14", "16", "15", "19", "03"))), 0)
cvs <- sapply(PER, function(p) 100 * sd(PA[[p]]) / mean(PA[[p]])); scvs <- sapply(PER, function(p) scv_of((2000:2024)[period_of(2000:2024) == p]))
claim(P, "CV 2000-2007", 38.3, cvs[1], 1); claim(P, "CV later periods, minimum", 29.8, min(cvs[-1]), 1); claim(P, "CV later periods, maximum", 30.9, max(cvs[-1]), 1)
claim(P, "SCV 2000-2007", 16.6, scvs[1], 1); claim(P, "SCV later periods, minimum", 10.8, min(scvs[-1]), 1); claim(P, "SCV later periods, maximum", 13.4, max(scvs[-1]), 1)
vy <- rd("variation_by_year")
cv_y <- sapply(2000:2024, function(y) { v <- sapply(PROVS, function(p) PY[[p]]["asr", y - 1999]); 100 * sd(v) / mean(v) })
check("Provinces", "Annual CV (independent) equals results/variation_by_year.csv", cv_y, vy$cv, 1e-6)
claim(P, "Highest annual CV (2003)", 42.5, max(cv_y), 1); claim(P, "Year of highest annual CV", 2003, (2000:2024)[which.max(cv_y)], 0)
claim(P, "Lowest annual CV (2011)", 24.8, min(cv_y), 1); claim(P, "Year of lowest annual CV", 2011, (2000:2024)[which.min(cv_y)], 0)
claim(P, "Provinces with significant rise 2000-2011", 15, sum(TR[, 2] > 0), 0)
claim(P, "Provinces with significant fall 2011-2019", 11, sum(TR[, 6] < 0), 0)
claim(P, "Only Zamora Chinchipe rose significantly 2011-2019", 1, as.numeric(identical(PROVS[TR[, 5] > 0], "19")), 0)
mq <- mo[mo$neighbours == "queen", ]; me <- mo[mo$neighbours == "queen, empirical Bayes", ]; mk <- mo[mo$neighbours == "knn4", ]
claim(P, "Moran's I 2000-2007", 0.08, moran_I(PA[[PER[1]]][match(geo$province, PROVS)], geo$nb$queen), 2)
claim(P, "Moran's I 2000-2007 p-value", 0.103, mq$p_value[mq$period == PER[1]], 3, "results (permutation, reproduced above)")
mI <- sapply(PER[-1], function(p) moran_I(PA[[p]][match(geo$province, PROVS)], geo$nb$queen))
claim(P, "Moran's I later periods, minimum", 0.24, min(mI), 2); claim(P, "Moran's I later periods, maximum", 0.30, max(mI), 2)
claim(P, "Moran's I later periods, largest p-value <= 0.013", 1, as.numeric(max(mq$p_value[mq$period %in% PER[-1]]) <= 0.013), 0, "results (permutation, reproduced above)")
claim(P, "Annual Moran's I significant every year from 2009", 1, as.numeric(all(vy$moran_p[vy$year >= 2009] < 0.05) && vy$moran_p[vy$year == 2008] >= 0.05), 0, "results/variation_by_year.csv (499 permutations)")
claim(P, "Empirical Bayes Moran's I, minimum 2008-2024 periods", 0.29, min(me$value[me$period %in% PER[-1]]), 2, "results (reproduced above)")
claim(P, "Empirical Bayes Moran's I, maximum 2008-2024 periods", 0.36, max(me$value[me$period %in% PER[-1]]), 2, "results (reproduced above)")
claim(P, "Moran's I with 4 nearest neighbors not significant in any period", 1, as.numeric(all(mk$p_value[mk$period %in% PER] >= 0.05)), 0, "results")
hh <- lisa[lisa$period %in% PER & !is.na(lisa$p) & lisa$p < 0.05 & lisa$quad == "High-High", ]
claim(P, "LISA high-high: Morona Santiago 2015-2019 and 2020-2024, Pastaza 2008-2014 only", 1,
      as.numeric(setequal(paste(hh$period, hh$province), c(paste(PER[3], "14"), paste(PER[4], "14"), paste(PER[2], "16")))), 0)
claim(P, "No province significant after FDR correction", 0, sum(lisa$fdr < 0.05, na.rm = TRUE), 0)
# beds
B <- "Beds"; b10n <- bnat$beds_available / bnat$population * 1e4
claim(B, "National beds per 10,000, 2000", 15.7, b10n[bnat$year == 2000], 1)
claim(B, "National beds per 10,000, 2011", 16.2, b10n[bnat$year == 2011], 1)
claim(B, "National beds per 10,000, 2024", 13.2, b10n[bnat$year == 2024], 1)
w11 <- bp[bp$year == 2011, ]; w24 <- bp[bp$year == 2024, ]
claim(B, "Provinces with fewer beds per inhabitant in 2024 than 2011", 13, sum(w24$b10k[match(w11$province, w24$province)] < w11$b10k), 0)
pb <- bp |> group_by(province) |> summarise(b = sum(bs) / sum(pop) * 1e4)
claim(B, "Lowest bed density 2000-2024 (Orellana)", 7.7, pb$b[pb$province == "22"], 1); claim(B, "Orellana lowest", 1, as.numeric(pb$province[which.min(pb$b)] == "22"), 0)
claim(B, "Highest bed density 2000-2024 (Tungurahua)", 18.9, pb$b[pb$province == "18"], 1); claim(B, "Tungurahua highest", 1, as.numeric(pb$province[which.max(pb$b)] == "18"), 0)
claim(B, "ACSC: % per 10% more beds", 4.5, BM[1, 1], 1); claim(B, "ACSC lower CI", 3.3, BM[1, 2], 1); claim(B, "ACSC upper CI", 5.8, BM[1, 3], 1)
claim(B, "Non-ACSC: % per 10% more beds", 4.8, BM[2, 1], 1); claim(B, "Non-ACSC lower CI", 3.6, BM[2, 2], 1); claim(B, "Non-ACSC upper CI", 6.1, BM[2, 3], 1)
claim(B, "Share: % per 10% more beds", -0.4, BM[3, 1], 1); claim(B, "Share lower CI", -1.4, BM[3, 2], 1); claim(B, "Share upper CI", 0.6, BM[3, 3], 1)
claim(B, "ACSC sensitivity analyses, minimum", 2.7, min(BM[c(4:10), 1]), 1); claim(B, "ACSC sensitivity analyses, maximum", 6.2, max(BM[c(4:10), 1]), 1)
claim(B, "Reduction attributable to beds 2011-2024 (%)", -5.7, contrib, 1); claim(B, "Observed change in standardized ratio 2011-2024 (%)", -22.1, sir_ch, 1)
claim(B, "Share of the decline attributable to beds (about a quarter)", 0.25, contrib / sir_ch, 1)
BP <- setNames(list(2000:2007, 2008:2014, 2015:2019, 2020:2024, 2000:2024), c(PER, "2000\u20132024"))
bt <- t(sapply(names(BP), function(p) {
  bb <- bp |> filter(year %in% BP[[p]]) |> group_by(province) |> summarise(b = sum(bs) / sum(pop) * 1e4)
  a <- if (p == "2000–2024") PA$Total else PA[[p]]
  ct_ <- suppressWarnings(cor.test(a[bb$province], bb$b, method = "spearman")); f <- lm(log(a[bb$province]) ~ log(bb$b))
  rz <- residuals(f)[match(geo$province, bb$province)]
  c(rho = unname(ct_$estimate), p = ct_$p.value, r2 = summary(f)$r.squared, mres = moran_I(rz, geo$nb$queen))
}))
claim(B, "Spearman rho by period, minimum", 0.25, min(bt[1:4, "rho"]), 2); claim(B, "Spearman rho by period, maximum", 0.43, max(bt[1:4, "rho"]), 2)
claim(B, "Spearman p-value 2000-2007", 0.035, bt[1, "p"], 3); claim(B, "Spearman not significant after 2007", 1, as.numeric(all(bt[2:4, "p"] >= 0.05)), 0)
claim(B, "Variance of provincial rates explained by beds, 2000-2024 (%)", 6, 100 * bt[5, "r2"], 0)
claim(B, "Moran's I of residuals 2008-2024 periods, minimum", 0.18, min(bt[2:4, "mres"]), 2)
claim(B, "Moran's I of residuals 2008-2024 periods, maximum", 0.22, max(bt[2:4, "mres"]), 2)
# segmented regression
R <- "Segmented regression"; F_ <- function(k) FITS[[k]]
claim(R, "Reform: ACSC pre-trend (% per year)", 2.8, F_("reform All ACSC")$pre[1], 1)
claim(R, "Reform: ACSC level change (%)", 7.1, F_("reform All ACSC")$level[1], 1)
claim(R, "Reform: ACSC level change lower CI", 1.4, F_("reform All ACSC")$level[2], 1); claim(R, "Reform: ACSC level change upper CI", 13.0, F_("reform All ACSC")$level[3], 1)
claim(R, "Reform: ACSC change in trend", -3.7, F_("reform All ACSC")$change[1], 1)
claim(R, "Reform: ACSC change in trend lower CI", -4.7, F_("reform All ACSC")$change[2], 1); claim(R, "Reform: ACSC change in trend upper CI", -2.7, F_("reform All ACSC")$change[3], 1)
claim(R, "Reform: ACSC post-interruption trend", -1.0, F_("reform All ACSC")$after[1], 1)
claim(R, "Reform: non-ACSC level change", 15.4, F_("reform Non-ACSC")$level[1], 1)
claim(R, "Reform: non-ACSC level change lower CI", 8.5, F_("reform Non-ACSC")$level[2], 1); claim(R, "Reform: non-ACSC level change upper CI", 22.7, F_("reform Non-ACSC")$level[3], 1)
claim(R, "Reform: non-ACSC change in trend", -3.2, F_("reform Non-ACSC")$change[1], 1)
claim(R, "Reform: non-ACSC change in trend lower CI", -4.0, F_("reform Non-ACSC")$change[2], 1); claim(R, "Reform: non-ACSC change in trend upper CI", -2.5, F_("reform Non-ACSC")$change[3], 1)
claim(R, "Reform: chronic level change", 19.2, F_("reform Chronic")$level[1], 1)
claim(R, "Reform: chronic level change lower CI", 10.3, F_("reform Chronic")$level[2], 1); claim(R, "Reform: chronic level change upper CI", 28.7, F_("reform Chronic")$level[3], 1)
fs <- F_("reform ACSC share of hospitalizations")
claim(R, "Reform: share pre-trend", -0.1, fs$pre[1], 1); claim(R, "Reform: share level change", -6.0, fs$level[1], 1)
claim(R, "Reform: share level change CI (3.6)", -3.6, fs$level[3], 1); claim(R, "Reform: share level change CI (8.4)", -8.4, fs$level[2], 1)
claim(R, "Reform: share post-interruption trend", -0.5, fs$after[1], 1)
claim(R, "Reform: ACSC below pre-reform share, 2008-2019 (%)", 8.5, 100 * (1 - fs$obs / fs$exp), 1)
claim(R, "Reform: share O/E (Table 3)", 0.92, fs$obs / fs$exp, 2)
claim(R, "January 2012: ACSC change in trend", -4.2, S12$acsc$change[1], 1)
claim(R, "January 2012: ACSC change in trend CI", -4.9, S12$acsc$change[2], 1); claim(R, "January 2012: ACSC change in trend CI ", -3.5, S12$acsc$change[3], 1)
claim(R, "January 2012: non-ACSC change in trend", -5.1, S12$non$change[1], 1)
claim(R, "January 2012: non-ACSC change in trend CI", -5.7, S12$non$change[2], 1); claim(R, "January 2012: non-ACSC change in trend CI ", -4.5, S12$non$change[3], 1)
claim(R, "January 2012: share change in trend", 0.8, S12$share$change[1], 1)
claim(R, "January 2012: share change in trend CI", 0.4, S12$share$change[2], 1); claim(R, "January 2012: share change in trend CI ", 1.2, S12$share$change[3], 1)
fc <- F_("covid All ACSC"); fn <- F_("covid Non-ACSC"); fcs <- F_("covid ACSC share of hospitalizations")
claim(R, "COVID: ACSC pre-trend", -1.7, fc$pre[1], 1)
claim(R, "COVID: ACSC level change", -44.0, fc$level[1], 1); claim(R, "COVID: ACSC level CI (30.6)", -30.6, fc$level[3], 1); claim(R, "COVID: ACSC level CI (54.9)", -54.9, fc$level[2], 1)
claim(R, "COVID: non-ACSC level change", -17.3, fn$level[1], 1); claim(R, "COVID: non-ACSC level CI (6.1)", -6.1, fn$level[3], 1); claim(R, "COVID: non-ACSC level CI (27.1)", -27.1, fn$level[2], 1)
claim(R, "COVID: share level change", -29.7, fcs$level[1], 1); claim(R, "COVID: share level CI (21.5)", -21.5, fcs$level[3], 1); claim(R, "COVID: share level CI (37.1)", -37.1, fcs$level[2], 1)
claim(R, "COVID: ACSC observed Mar 2020-Dec 2024", 671555, fc$obs, 0); claim(R, "COVID: ACSC expected", 818625, fc$exp, 0)
claim(R, "COVID: ACSC below expected (%)", 18.0, 100 * (1 - fc$obs / fc$exp), 1)
claim(R, "COVID: non-ACSC below expected (%)", 3.8, 100 * (1 - fn$obs / fn$exp), 1)
claim(R, "COVID: vaccine-preventable level change", -69.6, F_("covid Vaccine-preventable")$level[1], 1)
claim(R, "COVID: infectious level change", -48.7, F_("covid Infectious")$level[1], 1)
claim(R, "COVID: chronic level change", -33.8, F_("covid Chronic")$level[1], 1)
claim(R, "COVID: chronic below expected (%)", 4, 100 * (1 - F_("covid Chronic")$obs / F_("covid Chronic")$exp), 0)
claim(R, "COVID-coded discharges 2020", 55076, covid$u07_u10[covid$year == 2020], 0, "data (microdata count)")
claim(R, "COVID-coded discharges 2021", 64491, covid$u07_u10[covid$year == 2021], 0, "data (microdata count)")
claim(R, "ACSC pneumonia codes 2019", 10330, covid$acsc_pneumonia_codes[covid$year == 2019], 0, "data (microdata count)")
claim(R, "ACSC pneumonia codes 2020", 4800, covid$acsc_pneumonia_codes[covid$year == 2020], 0, "data (microdata count)")
claim(R, "ACSC pneumonia codes 2021", 2477, covid$acsc_pneumonia_codes[covid$year == 2021], 0, "data (microdata count)")
# discussion
D <- "Discussion"
claim(D, "Rise in discharges per 1,000, 2000-2011 (%)", 39, 100 * (ann$per1000[12] / ann$per1000[1] - 1), 0)
claim(D, "Growth in discharges 2000-2014 (%)", 78, 100 * (ann$all[15] / ann$all[1] - 1), 0)
g18 <- res |> filter(group == 18) |> group_by(year) |> summarise(n = sum(n))
claim(D, "Bacterial pneumonia, fold increase 2010-2018 (about fivefold)", 5, g18$n[g18$year == 2018] / g18$n[g18$year == 2010], 0)
claim(D, "COVID-coded discharges 2020-2021 (almost 120,000)", 120, sum(covid$u07_u10[covid$year %in% 2020:2021]) / 1000, 0)
claim(D, "Free care raised other hospitalizations (15%)", 15, F_("reform Non-ACSC")$level[1], 0)
claim(D, "and ACSC (7%)", 7, F_("reform All ACSC")$level[1], 0)

# =============================================================================
# 4. Report
# =============================================================================
checks <- bind_rows(checks); claims <- bind_rows(claims)
write_csv(checks, file.path(OUT, "audit_checks.csv")); write_csv(claims, file.path(OUT, "audit_claims.csv"))
rep <- c(sprintf("Audit run: %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
         sprintf("Method replication checks: %d of %d PASS", sum(checks$status == "PASS"), nrow(checks)),
         sprintf("Manuscript claims: %d of %d PASS", sum(claims$status == "PASS"), nrow(claims)), "",
         "Checks not passed:", capture.output(print(as.data.frame(checks[checks$status != "PASS", ]), row.names = FALSE)), "",
         "Claims to review:", capture.output(print(as.data.frame(claims[claims$status != "PASS", ]), row.names = FALSE)), "",
         capture.output(sessionInfo()))
writeLines(rep, file.path(OUT, "audit_report.txt"))
cat(rep[1:3], sep = "\n")
if (any(checks$status != "PASS") || any(claims$status != "PASS")) {
  cat("\nItems to review:\n"); print(as.data.frame(checks[checks$status != "PASS", c("area", "check", "max_abs_difference")]), row.names = FALSE)
  print(as.data.frame(claims[claims$status != "PASS", ]), row.names = FALSE)
}
