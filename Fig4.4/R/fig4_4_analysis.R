## =============================================================================
## Figure 4.4  First-progression pattern and post-progression cause of death
##             HCC first-line cohort (HAIC_Systemic / TACE_Systemic / Systemic)
## -----------------------------------------------------------------------------
## 4.4A  Absolute risk (no reference group needed)
##   A1  Adjusted 12-month CIF matrix (cause-specific Cox + G-computation)
##   A2  ARD vs population-standardised average
##   A3  12-month cause-of-death composition Q
##   A4  CIF-based proportional mortality ratio (PMR)
## 4.4B  Relative risk (reference = IH-NL)
##   B1  all-cause HR / csHR / sHR matrices + overall Wald test (df = 8)
##   B2  Lunn-McNeil cross-cause heterogeneity test + weighted-effect follow-up
## Validation (IH-TG reference) and sensitivity analyses
## -----------------------------------------------------------------------------
## Usage (from the Fig4.4 folder; the script file is UTF-8 encoded):
##   Rscript R/fig4_4_analysis.R
##   or in R: source("R/fig4_4_analysis.R", encoding = "UTF-8")
## Environment variables: N_BOOT (default 1000), N_CORES (default cores - 1;
##   forced to 1 on Windows), DATA_DIR (default "data"), OUT_DIR ("output")
## Inputs in DATA_DIR: raw.xlsx (2.1 response/PFS/OS), baseline.xlsx (1.3),
##   cod.xlsx (3 cause of mortality, sheet Cause_of_mortality_clean)
## Packages: survival, riskRegression, prodlim, coxphf, ggplot2, patchwork,
##   openxlsx, readxl, gridExtra, scales
## =============================================================================

suppressPackageStartupMessages({
  library(readxl); library(survival); library(riskRegression); library(prodlim)
  library(coxphf); library(ggplot2); library(patchwork); library(openxlsx)
  library(gridExtra); library(grid); library(parallel); library(scales)
})
options(stringsAsFactors = FALSE, warn = 1, width = 150)

DATA_DIR <- Sys.getenv("DATA_DIR", "data")
OUT_DIR  <- Sys.getenv("OUT_DIR", "output")
N_BOOT   <- as.integer(Sys.getenv("N_BOOT", "1000"))
N_CORES  <- as.integer(Sys.getenv("N_CORES", max(1, detectCores() - 1)))
if (.Platform$OS.type == "windows") N_CORES <- 1L
TAUS     <- c(6, 12, 24)          # months after first progression
TAU_MAIN <- 12
SEED     <- 20261010
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "png"), showWarnings = FALSE)

PAT_LEVELS <- c("IH-TG", "IH-NL", "IH-TG+NL", "EH-TG", "EH-NL", "EH-TG+NL",
                "Sys-TG", "Sys-NL", "Sys-TG+NL")
REF_MAIN   <- "IH-NL"
CAUSES     <- c(IH = 1, EH = 2, Sys = 3, NonHCC = 4, Undet = 5)
CAUSE_LAB  <- c(IH = "Intrahepatic HCC", EH = "Extrahepatic HCC", Sys = "Systemic HCC",
                NonHCC = "Non-HCC", Undet = "Undetermined", All = "All-cause")
COL_CAUSE  <- c(IH = "#2a78d6", EH = "#eb6834", Sys = "#1baf7a", NonHCC = "#eda100",
                Undet = "#b4b2a9")
INK <- "#0b0b0b"; INK2 <- "#52514e"; MUTED <- "#898781"; GRID_COL <- "#e1e0d9"

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

## -----------------------------------------------------------------------------
## 1. Load data and verify linkage
## -----------------------------------------------------------------------------
parse_date_mixed <- function(x) {          # Excel serials mixed with text dates
  x <- as.character(x)
  n <- suppressWarnings(as.numeric(x))
  out <- as.Date(n, origin = "1899-12-30")
  i <- is.na(n) & !is.na(x)
  out[i] <- as.Date(gsub("/", "-", substr(x[i], 1, 10)))
  out
}
rd <- function(f, sheet) suppressWarnings(read_excel(file.path(DATA_DIR, f), sheet = sheet, guess_max = 1e4))
raw <- rd("raw.xlsx", 1)
bl  <- rd("baseline.xlsx", 1)
cod <- rd("cod.xlsx", "Cause_of_mortality_clean")

dead <- raw$Death_status_0930 == 1
qc <- data.frame(
  check = c("rows (raw / baseline / cod)",
            "raw.序号 == baseline.no == cod.no",
            "Group identical in 3 tables",
            "baseline.patient_sn == cod.patient_sn",
            "raw.Baseline_date == baseline.首次目标治疗日期",
            "deaths (Death_status_0930 = 1)",
            "deaths with non-missing Cause_of_Mortality",
            "deaths with missing Cause_of_Mortality",
            "alive with a recorded Cause_of_Mortality",
            "PPOS_status = 1 but Death_status_0930 = 0 (progressors)"),
  result = c(paste(nrow(raw), nrow(bl), nrow(cod), sep = " / "),
             all(raw[["序号"]] == bl$no & bl$no == cod$no),
             all(raw$Group == bl$Group & bl$Group == cod$Group),
             all(bl$patient_sn == cod$patient_sn),
             sum(parse_date_mixed(raw$Baseline_date) != parse_date_mixed(bl[["治疗基线时间点_首次目标治疗日期"]]), na.rm = TRUE) == 0,
             sum(dead),
             sum(dead & !is.na(cod$Cause_of_Mortality)),
             sum(dead & is.na(cod$Cause_of_Mortality)),
             sum(!dead & !is.na(cod$Cause_of_Mortality)),
             sum(raw$PFS_status == 1 & raw$PFS_model != "Death" & raw$PPOS_status %in% 1 & !dead)))
stopifnot(qc$result[2:5] == "TRUE")
print(qc, right = FALSE)

## -----------------------------------------------------------------------------
## 2. Analysis cohort, progression pattern, outcome coding, covariates
## -----------------------------------------------------------------------------
## Pattern definition (classification document "解剖分类 x 肿瘤模式分类"):
##   components: IHG 肝内病灶增大, TTU 癌栓升级 (-> IHG), NIH 肝内新发病灶, nMVI 新发血管侵犯 (-> NIH),
##               EHG 肝外病灶增大, NEH 肝外新发转移, STG 系统性肿瘤增大
##   (1) IH-TG only IHG          (2) IH-NL only NIH          (3) IH-TG+NL IHG + NIH
##   (4) EH-TG only EHG          (5) EH-NL only NEH          (6) EH-TG+NL EHG + NEH
##   (7) Sys-TG only STG         (8) Sys-NL NIH + NEH        (9) Sys-TG+NL any other combination
COMP <- c("肝内病灶增大" = "IHG", "癌栓升级" = "IHG", "肝内新发病灶" = "NIH", "新发血管侵犯" = "NIH",
          "肝外病灶增大" = "EHG", "肝外新发转移" = "NEH", "系统性肿瘤增大" = "STG")
PAT_KEY <- c("IHG" = "IH-TG", "NIH" = "IH-NL", "IHG+NIH" = "IH-TG+NL",
             "EHG" = "EH-TG", "NEH" = "EH-NL", "EHG+NEH" = "EH-TG+NL",
             "STG" = "Sys-TG", "NEH+NIH" = "Sys-NL")
map_pattern <- function(s) vapply(strsplit(s, ";"), function(t) {
  t <- trimws(t)
  unk <- setdiff(t, names(COMP))
  if (length(unk)) stop("Unknown PFS_model component: ", paste(unk, collapse = ","))
  key <- paste(sort(unique(COMP[t])), collapse = "+")
  if (key %in% names(PAT_KEY)) PAT_KEY[[key]] else "Sys-TG+NL"
}, character(1))
## earlier site x form rule (sensitivity): any IH + EH involvement or STG = Sys; TG-only combinations
## such as STG + TTU or EHG + TTU were Sys-TG (this reproduces the original figure's Sys-TG n = 74)
map_pattern_v0 <- function(s) vapply(strsplit(s, ";"), function(t) {
  cmp <- unique(COMP[trimws(t)])
  ih <- any(cmp %in% c("IHG", "NIH")); eh <- any(cmp %in% c("EHG", "NEH")); sy <- "STG" %in% cmp
  tg <- any(cmp %in% c("IHG", "EHG", "STG")); nl <- any(cmp %in% c("NIH", "NEH"))
  site <- if (sy || (ih && eh)) "Sys" else if (ih) "IH" else "EH"
  form <- if (tg && nl) "TG+NL" else if (tg) "TG" else "NL"
  paste0(site, "-", form)
}, character(1))

cause_map <- c("Intrahepatic HCC-related Cause"  = 1,
               "Extrahepatic HCC-related Cause"  = 2,
               "Systemic HCC-related Cause"      = 3,
               "Liver-related Cause"             = 4,
               "Other Causes"                    = 4,
               "Treatment-related Cause"         = 4,
               "Unclassified HCC-related Cause"  = 5,
               "Unknown"                         = 5)
stopifnot(all(na.omit(unique(cod$Cause_of_Mortality)) %in% names(cause_map)))

all0 <- data.frame(
  id = raw[["序号"]], Group = raw$Group, PFS_status = raw$PFS_status, PFS_model = raw$PFS_model,
  TTP = raw$PFS, time = raw$PPOS, death = raw$Death_status_0930, PPOS_status = raw$PPOS_status,
  OS = raw$OS_month_0930, cause_raw = cod$Cause_of_Mortality, delete_cod = cod$Delete,
  delete_bl = bl[["是否需要删除"]],
  Age = bl$Age, sex = bl$sex, Etiology_HBV = bl$Etiology_HBV, ECOG_PS = bl$ECOG_PS,
  ALBI.Class = bl$ALBI.Class, AFP_400 = bl$AFP_400, PIVKA_1000 = bl$PIVKA_1000,
  Tumor_Number_3 = bl$Tumor_Number_3, Tumor_Diameter_class = bl$Tumor_Diameter_class,
  Vp_class = bl$Vp_class, PVTT = bl$PVTT, HVTT = bl$HVTT, EHS = bl$EHS, cirhosis = bl$cirhosis)

## cohort flow
f_prog   <- all0$PFS_status %in% 1 & all0$PFS_model != "Death"
f_ppos   <- f_prog & !is.na(all0$time)
del_flag <- (!is.na(all0$delete_cod) & !grepl("可纳入", all0$delete_cod)) | all0$delete_bl %in% 1
f_final  <- f_ppos & !del_flag
flow <- data.frame(
  step = c("Total patients", "PFS event (progression or death)", "  excluded: death without progression",
           "First radiological progression", "  excluded: PPOS missing (PFS > OS, data error)",
           "  excluded: flagged for deletion (non-HCC histology / other cancer / external treatment / long LTFU)",
           "Analysis cohort"),
  n = c(nrow(all0), sum(all0$PFS_status %in% 1), sum(all0$PFS_model %in% "Death"),
        sum(f_prog), sum(f_prog & is.na(all0$time)), sum(f_ppos & del_flag), sum(f_final)))
print(flow, right = FALSE)

prep <- function(x, rule = map_pattern) {
  x$pattern <- factor(rule(x$PFS_model), levels = PAT_LEVELS)
  x$site <- factor(sub("-.*", "", as.character(x$pattern)), levels = c("IH", "EH", "Sys"))
  x$form <- factor(sub("^[A-Za-z]+-", "", as.character(x$pattern)), levels = c("TG", "NL", "TG+NL"))
  x$status <- ifelse(x$death == 0, 0L,
                     ifelse(is.na(x$cause_raw), 5L, as.integer(cause_map[x$cause_raw])))
  ## 5 patients censored at PPOS = 0 -> half a day so that counting-process data are valid
  x$time <- pmax(x$time, 0.5 / 30.4375)
  x$Group <- factor(x$Group, levels = c("HAIC_Systemic", "TACE_Systemic", "Systemic"))
  x$Age10 <- x$Age / 10
  x$Male  <- as.integer(x$sex == "Male")
  x$HBV   <- as.integer(x$Etiology_HBV)
  x$ECOG  <- factor(x$ECOG_PS, levels = 0:2)
  x$ALBI  <- factor(x$ALBI.Class, levels = c("1/2a", "2b/3"))
  x$AFP400 <- as.integer(x$AFP_400); x$PIVKA1000 <- as.integer(x$PIVKA_1000)
  x$TumorNum3 <- as.integer(x$Tumor_Number_3)
  x$TumorSize <- factor(x$Tumor_Diameter_class, levels = c("(0, 5]", "(5, 10]", ">10", "Infiltrative"))
  ## Vp_class missing (n = 9): Absent when PVTT = 0, otherwise the modal Vp class among PVTT = 1
  vp <- x$Vp_class
  mode_pvtt <- names(which.max(table(vp[x$PVTT %in% 1])))
  vp[is.na(vp)] <- ifelse(x$PVTT[is.na(vp)] %in% 1, mode_pvtt, "Absent")
  x$Vp <- factor(vp, levels = c("Absent", "Vp1-Vp2", "Vp3-Vp4"))
  x$HVTT <- as.integer(x$HVTT); x$EHS <- as.integer(x$EHS); x$Cirrhosis <- as.integer(x$cirhosis)
  x
}
COVS <- c("TTP", "Age10", "Male", "HBV", "ECOG", "ALBI", "AFP400", "PIVKA1000", "TumorNum3",
          "TumorSize", "Vp", "HVTT", "EHS", "Cirrhosis")
d      <- prep(all0[f_final, ])
d_del  <- prep(all0[f_ppos, ])        # sensitivity: flagged patients kept
d_v0   <- prep(all0[f_final, ], map_pattern_v0)   # sensitivity: earlier pattern rule
pattern_check <- with(all0[f_final, ], {
  combo <- vapply(strsplit(PFS_model, ";"), function(t) paste(sort(unique(trimws(t))), collapse = ";"), "")
  z <- as.data.frame(table(PFS_model_components = combo, pattern = map_pattern(PFS_model),
                           earlier_rule = map_pattern_v0(PFS_model)), stringsAsFactors = FALSE)
  z <- z[z$Freq > 0, ]; z[order(factor(z$pattern, PAT_LEVELS), -z$Freq), ] })
log_msg("Pattern rule: ", sum(d$pattern != d_v0$pattern), " patients classified differently from the earlier rule")
stopifnot(!anyNA(d[, c("time", "status", "pattern", COVS, "Group")]))
log_msg("Analysis cohort n = ", nrow(d), "; deaths = ", sum(d$status > 0))

## events per pattern x cause (whole follow-up, and within each tau)
ev_tab <- function(x, tau = Inf) {
  st <- ifelse(x$time <= tau, x$status, 0L)
  t <- table(factor(x$pattern, PAT_LEVELS), factor(st, 0:5))
  data.frame(pattern = PAT_LEVELS, n = as.vector(table(factor(x$pattern, PAT_LEVELS))),
             IH = t[, "1"], EH = t[, "2"], Sys = t[, "3"], NonHCC = t[, "4"], Undet = t[, "5"],
             All = rowSums(t[, -1]), row.names = NULL)
}
events_all <- ev_tab(d)
events_tau <- lapply(setNames(TAUS, paste0("m", TAUS)), function(t) ev_tab(d, t))
print(events_all)

## follow-up and numbers at risk
fu_km <- survfit(Surv(time, death == 0) ~ 1, data = d)
fu_median <- unname(summary(fu_km)$table[c("median", "0.95LCL", "0.95UCL")])
risk_tab <- with(summary(survfit(Surv(time, status > 0) ~ pattern, data = d), times = TAUS, extend = TRUE),
                 data.frame(pattern = sub("pattern=", "", strata), month = time, n_risk = n.risk))
risk_tab <- reshape(risk_tab, idvar = "pattern", timevar = "month", direction = "wide")
log_msg(sprintf("Median post-progression follow-up (reverse KM): %.1f months (%.1f-%.1f)",
                fu_median[1], fu_median[2], fu_median[3]))

## -----------------------------------------------------------------------------
## 3. CIF engine: stratified cause-specific Cox + Breslow baseline + G-computation
##    F_j(t | x) = sum_{s <= t} S(s- | x) dLambda_j(s | x),  S = exp(-sum_j Lambda_j)
##    Causes 1-4: full model (pattern + TTP + covariates, strata = Group)
##    Cause 5 (undetermined, 16 deaths): Nelson-Aalen within Group (nuisance; keeps
##    sum of cause-specific CIFs equal to all-cause death)
## -----------------------------------------------------------------------------
cov_formula <- as.formula(paste("~", paste(COVS, collapse = " + ")))
make_X <- function(x, pat = NULL) {
  p <- if (is.null(pat)) x$pattern else factor(rep(pat, nrow(x)), levels = PAT_LEVELS)
  P <- sapply(setdiff(PAT_LEVELS, REF_MAIN), function(k) as.numeric(p == k))
  if (is.null(dim(P))) P <- matrix(P, nrow = 1, dimnames = list(NULL, names(P)))
  C <- model.matrix(cov_formula, x)[, -1, drop = FALSE]
  cbind(P, C)
}
breslow_inc <- function(time, ev, eta, grid) {     # Breslow increments at grid times
  o <- order(time); ts <- time[o]; rc <- rev(cumsum(rev(eta[o])))
  R <- rc[findInterval(grid, ts, left.open = TRUE) + 1]
  dN <- tabulate(match(time[ev], grid), nbins = length(grid))
  ifelse(dN > 0, dN / R, 0)
}
cif_engine <- function(x, taus = TAUS, full_causes = 1:4, na_causes = 5,
                       pats = PAT_LEVELS, ties = "breslow", return_fits = FALSE) {
  X <- make_X(x)
  betas <- lapply(full_causes, function(j) {
    f <- suppressWarnings(coxph(Surv(x$time, x$status == j) ~ X + strata(x$Group), ties = ties))
    b <- coef(f); b[is.na(b)] <- 0; names(b) <- colnames(X); b
  })
  cn <- c(names(CAUSES)[c(full_causes, na_causes)])
  out <- array(0, dim = c(length(pats), length(cn) + 1, length(taus)),
               dimnames = list(pats, c(cn, "All"), paste0("m", taus)))
  Xc <- X[, -(1:(length(PAT_LEVELS) - 1)), drop = FALSE]
  bcov <- lapply(betas, function(b) b[colnames(Xc)])
  for (g in levels(x$Group)) {
    ig <- which(x$Group == g); xg <- x[ig, ]
    grid <- sort(unique(xg$time[xg$status > 0 & xg$time <= max(taus)]))
    if (!length(grid)) next
    M <- length(grid); keep <- outer(grid, taus, "<=") * 1
    eta_obs <- lapply(seq_along(full_causes), function(m) exp(drop(X[ig, , drop = FALSE] %*% betas[[m]])))
    dL0 <- c(lapply(seq_along(full_causes), function(m)
                 breslow_inc(xg$time, xg$status == full_causes[m], eta_obs[[m]], grid)),
             lapply(na_causes, function(j) breslow_inc(xg$time, xg$status == j, rep(1, nrow(xg)), grid)))
    L0 <- lapply(dL0, cumsum)
    eta_cov <- lapply(seq_along(full_causes), function(m) drop(Xc[ig, , drop = FALSE] %*% bcov[[m]]))
    for (k in pats) {
      E <- lapply(seq_along(full_causes), function(m) {
        bp <- if (k == REF_MAIN) 0 else betas[[m]][k]
        exp(eta_cov[[m]] + bp)
      })
      E <- c(E, lapply(na_causes, function(j) rep(1, length(ig))))
      H <- Reduce(`+`, Map(function(e, l) outer(e, l), E, L0))
      Sm <- exp(-cbind(0, H[, -M, drop = FALSE]))
      for (m in seq_along(E)) {
        cif <- E[[m]] * (Sm %*% (dL0[[m]] * keep))       # n_g x length(taus)
        out[k, m, ] <- out[k, m, ] + colSums(cif)
      }
    }
  }
  out <- out / nrow(x)
  out[, "All", ] <- apply(out[, cn, , drop = FALSE], c(1, 3), sum)
  if (return_fits) attr(out, "betas") <- betas
  out
}
derive <- function(F, w) {
  cn <- setdiff(dimnames(F)[[2]], "All")
  Fbar <- apply(F, 2:3, function(v) sum(w * v))                     # (cause+All) x tau
  ARD  <- sweep(F, 2:3, Fbar)
  Q    <- sweep(F[, cn, , drop = FALSE], c(1, 3), F[, "All", ], "/")
  Qbar <- sweep(Fbar[cn, , drop = FALSE], 2, Fbar["All", ], "/")
  PMR  <- sweep(Q, 2:3, Qbar, "/")
  list(F = F, Fbar = Fbar, ARD = ARD, Q = Q, Qbar = Qbar, PMR = PMR)
}
pat_w <- function(x) as.vector(table(factor(x$pattern, PAT_LEVELS))) / nrow(x)

## ---- 3a. self-check of the engine against riskRegression::CSC + predictRisk ----
local({
  dv <- d; dv$status[dv$status == 5] <- 4L                 # 4 fully-modelled causes
  ## engine without G-computation: predict for the first 40 patients at their own pattern
  sub <- dv[1:40, ]
  csc <- CSC(as.formula(paste("Hist(time, status) ~ pattern +", paste(COVS, collapse = "+"), "+ strata(Group)")),
             data = dv, ties = "breslow")
  rr <- sapply(1:4, function(j) predictRisk(csc, newdata = sub, times = TAU_MAIN, cause = j, product.limit = FALSE))
  mine <- t(sapply(seq_len(nrow(sub)), function(i) {
    tmp <- dv; Fi <- NULL
    ## exploit the engine by standardising over a one-row population with the full data fit
    X <- make_X(dv); b <- lapply(1:4, function(j) coef(coxph(Surv(dv$time, dv$status == j) ~ X + strata(dv$Group), ties = "breslow")))
    g <- sub$Group[i]; ig <- which(dv$Group == g); xg <- dv[ig, ]
    grid <- sort(unique(xg$time[xg$status > 0 & xg$time <= TAU_MAIN]))
    dL0 <- lapply(1:4, function(j) breslow_inc(xg$time, xg$status == j, exp(drop(X[ig, ] %*% b[[j]])), grid))
    xi <- make_X(sub[i, , drop = FALSE])
    E <- sapply(1:4, function(j) exp(sum(xi * b[[j]])))
    H <- Reduce(`+`, lapply(1:4, function(j) E[j] * cumsum(dL0[[j]])))
    Sm <- exp(-c(0, H[-length(H)]))
    sapply(1:4, function(j) sum(Sm * E[j] * dL0[[j]]))
  }))
  diff <- max(abs(mine - rr))
  log_msg(sprintf("Engine vs riskRegression::predictRisk: max |diff| in 12-month CIF = %.2e", diff))
  if (diff > 1e-6) warning("CIF engine differs from riskRegression; check versions")
  assign("engine_check", diff, envir = .GlobalEnv)
})

## ---- 3b. point estimates + bootstrap ----
log_msg("Point estimates ...")
F_hat <- cif_engine(d)
est   <- derive(F_hat, pat_w(d))

log_msg("Bootstrap: B = ", N_BOOT, ", cores = ", N_CORES)
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED)
boot_fun <- function(b) {
  idx <- sample.int(nrow(d), nrow(d), replace = TRUE); db <- d[idx, ]
  tryCatch({ Fb <- cif_engine(db); derive(Fb, pat_w(db)) }, error = function(e) NULL)
}
t_boot <- system.time(
  boots <- if (N_CORES > 1) mclapply(seq_len(N_BOOT), boot_fun, mc.cores = N_CORES, mc.set.seed = TRUE)
           else lapply(seq_len(N_BOOT), boot_fun))
boots <- Filter(Negate(is.null), boots)
log_msg("Bootstrap done: ", length(boots), " successful replicates in ", round(t_boot[3] / 60, 1), " min")

boot_ci <- function(name) {
  A <- simplify2array(lapply(boots, `[[`, name))
  nd <- length(dim(A))
  list(lo = apply(A, seq_len(nd - 1), quantile, 0.025, na.rm = TRUE),
       hi = apply(A, seq_len(nd - 1), quantile, 0.975, na.rm = TRUE))
}
ci <- lapply(setNames(nm = c("F", "Fbar", "ARD", "Q", "Qbar", "PMR")), boot_ci)

## tidy long tables -------------------------------------------------------------
long_tab <- function(name, causes = c("IH", "EH", "Sys", "NonHCC", "Undet", "All")) {
  pe <- est[[name]]; lo <- ci[[name]]$lo; hi <- ci[[name]]$hi
  causes <- intersect(causes, dimnames(pe)[[2]])
  do.call(rbind, lapply(paste0("m", TAUS), function(tt) do.call(rbind, lapply(causes, function(cc)
    data.frame(measure = name, month = as.numeric(sub("m", "", tt)), pattern = PAT_LEVELS, cause = cc,
               est = pe[, cc, tt], lo = lo[, cc, tt], hi = hi[, cc, tt],
               events_by_tau = events_tau[[tt]][[cc]], row.names = NULL)))))
}
tab_F <- long_tab("F"); tab_ARD <- long_tab("ARD"); tab_Q <- long_tab("Q"); tab_PMR <- long_tab("PMR")
tab_Fbar <- do.call(rbind, lapply(paste0("m", TAUS), function(tt) data.frame(
  month = as.numeric(sub("m", "", tt)), cause = dimnames(est$Fbar)[[1]], Fbar = est$Fbar[, tt],
  lo = ci$Fbar$lo[, tt], hi = ci$Fbar$hi[, tt], row.names = NULL)))

## -----------------------------------------------------------------------------
## 4. Unadjusted Aalen-Johansen CIF
## -----------------------------------------------------------------------------
d$fstat <- factor(d$status, 0:5, c("censor", "IH", "EH", "Sys", "NonHCC", "Undet"))
aj <- survfit(Surv(time, fstat) ~ pattern, data = d)
saj <- summary(aj, times = TAUS, extend = TRUE)
km_all <- summary(survfit(Surv(time, status > 0) ~ pattern, data = d), times = TAUS, extend = TRUE)
st_cols <- c("IH", "EH", "Sys", "NonHCC", "Undet")
tab_AJ <- do.call(rbind, lapply(st_cols, function(cc) data.frame(
  month = saj$time, pattern = sub("pattern=", "", saj$strata), cause = cc,
  est = saj$pstate[, match(cc, saj$states)], lo = saj$lower[, match(cc, saj$states)],
  hi = saj$upper[, match(cc, saj$states)])))
tab_AJ <- rbind(tab_AJ, data.frame(month = km_all$time, pattern = sub("pattern=", "", km_all$strata),
                                   cause = "All", est = 1 - km_all$surv, lo = 1 - km_all$upper, hi = 1 - km_all$lower))
tab_AJ$events_by_tau <- mapply(function(p, cc, m) events_tau[[paste0("m", m)]][[cc]][PAT_LEVELS == p],
                               tab_AJ$pattern, tab_AJ$cause, tab_AJ$month)

## -----------------------------------------------------------------------------
## 5. 4.4B-B1  HR / csHR / sHR matrices
## -----------------------------------------------------------------------------
NE_MIN <- 3
wald_overall <- function(fit, idx) {
  b <- coef(fit)[idx]; V <- vcov(fit)[idx, idx]
  ok <- is.finite(b)
  chi <- as.numeric(t(b[ok]) %*% solve(V[ok, ok], b[ok]))
  c(chisq = chi, df = sum(ok), p = pchisq(chi, sum(ok), lower.tail = FALSE))
}
hr_matrix <- function(x, ref = REF_MAIN, adjust = TRUE, type = c("cs", "fg", "firth"),
                      outcomes = c("All", "IH", "EH", "Sys", "NonHCC")) {
  type <- match.arg(type)
  x$pattern <- relevel(factor(x$pattern, PAT_LEVELS), ref)
  covs <- if (adjust) COVS else character(0)
  res <- list(); ovl <- list(); zph <- list()
  for (o in outcomes) {
    if (type == "fg" && o == "All") next
    ev <- if (o == "All") as.integer(x$status > 0) else as.integer(x$status == CAUSES[[o]])
    nev <- tapply(ev, x$pattern, sum)
    if (type == "cs") {
      fit <- coxph(as.formula(paste("Surv(time, ev) ~", paste(c("pattern", covs, "strata(Group)"), collapse = " + "))),
                   data = cbind(x, ev = ev), ties = "efron")
      idx <- grep("^pattern", names(coef(fit)))
      b <- coef(fit)[idx]; se <- sqrt(diag(vcov(fit)))[idx]
      lo <- exp(b - 1.96 * se); hi <- exp(b + 1.96 * se); p <- 2 * pnorm(-abs(b / se))
      ovl[[o]] <- wald_overall(fit, idx)
      z <- cox.zph(fit)$table
      zph[[o]] <- data.frame(outcome = o, term = rownames(z), chisq = z[, "chisq"], df = z[, "df"], p = z[, "p"])
    } else if (type == "fg") {
      x$fstat <- factor(x$status, 0:5, c("censor", "IH", "EH", "Sys", "NonHCC", "Undet"))
      fg <- finegray(as.formula(paste("Surv(time, fstat) ~", paste(c("pattern", covs, "id", "Group", "strata(Group)"), collapse = " + "))),
                     data = x, etype = o)
      fit <- coxph(as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~", paste(c("pattern", covs, "strata(Group)"), collapse = " + "))),
                   data = fg, weights = fgwt, cluster = id)
      idx <- grep("^pattern", names(coef(fit)))
      b <- coef(fit)[idx]; se <- sqrt(diag(vcov(fit)))[idx]
      lo <- exp(b - 1.96 * se); hi <- exp(b + 1.96 * se); p <- 2 * pnorm(-abs(b / se))
      ovl[[o]] <- wald_overall(fit, idx)
    } else {                                   # Firth (coxphf has no strata: Group as covariate)
      fit <- coxphf(as.formula(paste("Surv(time, ev) ~", paste(c("pattern", covs, "Group"), collapse = " + "))),
                    data = cbind(x, ev = ev), pl = FALSE, maxit = 200)
      idx <- grep("^pattern", names(coef(fit)))
      b <- coef(fit)[idx]; lo <- fit$ci.lower[idx]; hi <- fit$ci.upper[idx]; p <- fit$prob[idx]
      V <- fit$var[idx, idx]; chi <- as.numeric(t(b) %*% solve(V, b))
      ovl[[o]] <- c(chisq = chi, df = length(b), p = pchisq(chi, length(b), lower.tail = FALSE))
    }
    pats <- sub("^pattern", "", names(b))
    r <- data.frame(model = type, outcome = o, pattern = c(ref, pats), events = as.vector(nev[c(ref, pats)]),
                    n = as.vector(table(x$pattern)[c(ref, pats)]),
                    HR = c(1, exp(b)), lo = c(NA, lo), hi = c(NA, hi), p = c(NA, p), row.names = NULL)
    r$NE <- r$pattern != ref & (r$events < NE_MIN | !is.finite(r$hi) | r$hi > 1e3)
    res[[o]] <- r
  }
  res <- do.call(rbind, res)
  res$q_BH <- NA
  ok <- !is.na(res$p) & !res$NE
  res$q_BH[ok] <- p.adjust(res$p[ok], "BH")
  res$pattern <- factor(res$pattern, PAT_LEVELS)
  list(tab = res[order(res$outcome, res$pattern), ],
       overall = data.frame(model = type, outcome = names(ovl), do.call(rbind, ovl), row.names = NULL),
       zph = if (length(zph)) do.call(rbind, zph) else NULL)
}
log_msg("B1: cause-specific Cox ...")
B_cs <- hr_matrix(d, type = "cs")
log_msg("B1: Fine-Gray ...")
B_fg <- hr_matrix(d, type = "fg", outcomes = c("IH", "EH", "Sys", "NonHCC"))
log_msg("Sensitivity: Firth ...")
B_firth <- hr_matrix(d, type = "firth")

## -----------------------------------------------------------------------------
## 6. 4.4B-B2  Lunn-McNeil heterogeneity across the 3 HCC-related causes
## -----------------------------------------------------------------------------
log_msg("B2: Lunn-McNeil ...")
lm_dat <- do.call(rbind, lapply(1:3, function(cc) {
  z <- d; z$cause <- cc; z$ev <- as.integer(d$status == cc); z }))
lm_dat$cause <- factor(lm_dat$cause, 1:3, c("IH", "EH", "Sys"))
Pm <- make_X(lm_dat)[, 1:(length(PAT_LEVELS) - 1)]
Cm <- model.matrix(cov_formula, lm_dat)[, -1]
Ic <- sapply(levels(lm_dat$cause), function(cc) as.numeric(lm_dat$cause == cc))
by_cause <- function(M, causes = colnames(Ic)) do.call(cbind, lapply(causes, function(cc) {
  z <- M * Ic[, cc]; colnames(z) <- paste0(colnames(M), ":", cc); z }))
Sd <- cbind(EH = as.numeric(lm_dat$site == "EH"), Sys = as.numeric(lm_dat$site == "Sys"))
Fd <- cbind(NL = as.numeric(lm_dat$form == "NL"), `TG+NL` = as.numeric(lm_dat$form == "TG+NL"))
X_m0   <- cbind(Pm, by_cause(Cm))
X_m1   <- cbind(by_cause(Pm), by_cause(Cm))
X_site <- cbind(X_m0, by_cause(Sd, c("EH", "Sys")))
X_form <- cbind(X_m0, by_cause(Fd, c("EH", "Sys")))
lmfit <- function(X) coxph(Surv(lm_dat$time, lm_dat$ev) ~ X + strata(lm_dat$cause, lm_dat$Group), ties = "efron")
m0 <- lmfit(X_m0); m1 <- lmfit(X_m1); ms <- lmfit(X_site); mf <- lmfit(X_form)
lrt <- function(a, b, label) {
  chi <- 2 * (b$loglik[2] - a$loglik[2]); df <- length(coef(b)) - length(coef(a))
  data.frame(test = label, chisq = chi, df = df, p = pchisq(chi, df, lower.tail = FALSE))
}
LM <- rbind(lrt(m0, m1, "pattern x cause (9 patterns; IH / EH / Systemic HCC death)"),
            lrt(m0, ms, "site x cause (IH / EH / Sys site), pattern main effect kept"),
            lrt(m0, mf, "form x cause (TG / NL / TG+NL), pattern main effect kept"))
## check: m1 equals 3 separate cause-specific fits
chk_lm <- max(abs(coef(m1)[paste0("X", colnames(Pm), ":Sys")] -
                  coef(coxph(Surv(time, status == 3) ~ make_X(d) + strata(Group), data = d))[1:8]))
log_msg(sprintf("Lunn-McNeil m1 vs separate Sys csCox: max |diff beta| = %.1e", chk_lm))
print(LM)

## follow-up: weighted effect coding, per pattern Wald test of equal effect across causes (df = 2)
w <- pat_w(d); b1 <- coef(m1); V1 <- vcov(m1); nm1 <- names(b1)
theta_L <- function(k, cc) {                   # theta_kc = beta_kc - sum_l w_l beta_lc
  L <- setNames(numeric(length(b1)), nm1)
  for (l in setdiff(PAT_LEVELS, REF_MAIN)) {
    a <- (l == k) - w[PAT_LEVELS == l]
    L[paste0("X", l, ":", cc)] <- a
  }
  L
}
FU <- do.call(rbind, lapply(PAT_LEVELS, function(k) {
  Ls <- lapply(c("IH", "EH", "Sys"), function(cc) theta_L(k, cc))
  th <- sapply(Ls, function(L) sum(L * b1)); se <- sapply(Ls, function(L) sqrt(drop(t(L) %*% V1 %*% L)))
  C <- rbind(Ls[[1]] - Ls[[2]], Ls[[1]] - Ls[[3]])
  cb <- drop(C %*% b1); chi <- drop(t(cb) %*% solve(C %*% V1 %*% t(C), cb))
  data.frame(pattern = k,
             eHR_IH = exp(th[1]), eHR_IH_lo = exp(th[1] - 1.96 * se[1]), eHR_IH_hi = exp(th[1] + 1.96 * se[1]),
             eHR_EH = exp(th[2]), eHR_EH_lo = exp(th[2] - 1.96 * se[2]), eHR_EH_hi = exp(th[2] + 1.96 * se[2]),
             eHR_Sys = exp(th[3]), eHR_Sys_lo = exp(th[3] - 1.96 * se[3]), eHR_Sys_hi = exp(th[3] + 1.96 * se[3]),
             chisq = chi, df = 2, p = pchisq(chi, 2, lower.tail = FALSE))
}))
FU$p_Holm <- p.adjust(FU$p, "holm")
print(FU[, c("pattern", "chisq", "p", "p_Holm")])

## -----------------------------------------------------------------------------
## 7. Validation against the original IH-TG-referenced matrix
## -----------------------------------------------------------------------------
log_msg("Validation ...")
target <- c(allcause_SysTG = 1.96, csHR_Sys_SysTG = 4.11, allcause_IHNL = 0.76)
val_one <- function(x, label, adjust, event = c("death", "ppos"), covs = COVS) {
  event <- match.arg(event)
  if (event == "ppos") { x$status[x$death == 0 & x$PPOS_status %in% 1] <- 5L }
  x$pattern <- relevel(factor(x$pattern, PAT_LEVELS), "IH-TG")
  rhs <- paste(c("pattern", if (adjust) covs, "strata(Group)"), collapse = " + ")
  fa <- coxph(as.formula(paste("Surv(time, status > 0) ~", rhs)), data = x)
  fs <- coxph(as.formula(paste("Surv(time, status == 3) ~", rhs)), data = x)
  g <- function(f, k) { b <- coef(f)[paste0("pattern", k)]; s <- sqrt(vcov(f)[paste0("pattern", k), paste0("pattern", k)])
    sprintf("%.2f (%.2f-%.2f)", exp(b), exp(b - 1.96 * s), exp(b + 1.96 * s)) }
  data.frame(variant = label, n = nrow(x), deaths = sum(x$status > 0),
             allcause_SysTG = g(fa, "Sys-TG"), csHR_Sys_SysTG = g(fs, "Sys-TG"), allcause_IHNL = g(fa, "IH-NL"))
}
VAL <- rbind(
  data.frame(variant = "ORIGINAL matrix (target)", n = NA, deaths = NA, allcause_SysTG = "1.96 (1.48-2.61)",
             csHR_Sys_SysTG = "4.11", allcause_IHNL = "0.76"),
  val_one(d, "Main cohort, adjusted, strata(Group)", TRUE),
  val_one(d_v0, "Earlier pattern rule (Sys-TG n=74), unadjusted", FALSE),
  val_one(d_v0, "Earlier pattern rule, baseline covariates without TTP", TRUE, covs = setdiff(COVS, "TTP")),
  val_one(d, "Main cohort, unadjusted, strata(Group)", FALSE),
  val_one(d_del, "Incl. flagged pts, adjusted", TRUE),
  val_one(d_del, "Incl. flagged pts, unadjusted", FALSE),
  val_one(d_del, "Incl. flagged, unadjusted, PPOS_status events", FALSE, "ppos"),
  val_one(d, "Main cohort, adjusted, PPOS_status events", TRUE, "ppos"),
  val_one(d, "Main cohort, baseline covariates without TTP", TRUE, covs = setdiff(COVS, "TTP")),
  val_one(d_del, "Incl. flagged pts, baseline covariates without TTP", TRUE, covs = setdiff(COVS, "TTP")))
print(VAL, right = FALSE)

## -----------------------------------------------------------------------------
## 8. Sensitivity analyses
## -----------------------------------------------------------------------------
log_msg("Sensitivity analyses ...")
## S5 exclude undetermined deaths
d_nound <- d[d$status != 5, ]
S_undet_hr  <- hr_matrix(d_nound, type = "cs")
S_undet_cif <- derive(cif_engine(d_nound, full_causes = 1:4, na_causes = integer(0)), pat_w(d_nound))
## include flagged patients
S_del_hr <- hr_matrix(d_del, type = "cs")
## earlier pattern rule
S_v0_hr <- hr_matrix(d_v0, type = "cs")
## treatment-group interaction for post-progression all-cause death
rhs0 <- paste(c("pattern", COVS, "strata(Group)"), collapse = " + ")
fi0 <- coxph(as.formula(paste("Surv(time, status > 0) ~", rhs0)), data = d)
fi_site <- coxph(as.formula(paste("Surv(time, status > 0) ~", rhs0, "+ site:Group")), data = d)
fi_pat  <- coxph(as.formula(paste("Surv(time, status > 0) ~", rhs0, "+ pattern:Group")), data = d)
INTER <- rbind(lrt(fi0, fi_site, "site x treatment group (all-cause death after progression)"),
               lrt(fi0, fi_pat,  "pattern x treatment group (all-cause death after progression)"))
INTER$df <- c(sum(!is.na(coef(fi_site))) - sum(!is.na(coef(fi0))), sum(!is.na(coef(fi_pat))) - sum(!is.na(coef(fi0))))
INTER$p <- pchisq(INTER$chisq, INTER$df, lower.tail = FALSE)
print(INTER)

## =============================================================================
## 9. Figures
## =============================================================================
fmt_ci  <- function(e, l, h, d = 1, s = 100) sprintf(paste0("%.", d, "f\n(%.", d, "f–%.", d, "f)"), s * e, s * l, s * h)
fmt_hr  <- function(e, l, h) sprintf("%.2f\n(%.2f–%.2f)", e, l, h)
wrap    <- function(s, w = 170) paste(strwrap(s, w), collapse = "\n")
stars   <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
fmt_p   <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
OUT_LAB <- c(All = "All-cause", IH = "Intrahepatic\nHCC death", EH = "Extrahepatic\nHCC death",
             Sys = "Systemic\nHCC death", NonHCC = "Non-HCC\ndeath", Undet = "Undetermined")
theme_mat <- function() theme_minimal(base_size = 10) + theme(
  panel.grid = element_blank(), axis.title = element_blank(),
  axis.text = element_text(colour = INK), axis.text.x = element_text(face = "bold", size = 9),
  axis.text.y = element_text(face = "bold", size = 9.5), legend.position = "bottom",
  plot.title = element_text(face = "bold", size = 12, colour = INK),
  plot.subtitle = element_text(size = 8.5, colour = INK2),
  plot.caption = element_text(size = 7.5, colour = INK2, hjust = 0),
  plot.background = element_rect(fill = "white", colour = NA))

## generic CIF matrix (rows = pattern, cols = outcome)
plot_cif_matrix <- function(tab, month, title, subtitle = NULL, row_order = NULL, caption = NULL, legend = TRUE) {
  x <- tab[tab$month == month & tab$cause %in% c("All", "IH", "EH", "Sys", "NonHCC"), ]
  if (is.null(row_order)) row_order <- PAT_LEVELS
  x$pattern <- factor(x$pattern, rev(row_order))
  x$cause <- factor(x$cause, c("All", "IH", "EH", "Sys", "NonHCC"), OUT_LAB[c("All", "IH", "EH", "Sys", "NonHCC")])
  x$zero <- x$events_by_tau == 0
  x$lab <- ifelse(x$zero, "0 events", fmt_ci(x$est, x$lo, x$hi))
  x$nlab <- ifelse(x$zero, "", paste0("n=", x$events_by_tau))
  x$fill <- ifelse(x$zero, NA, 100 * x$est)
  x$txtcol <- ifelse(!is.na(x$fill) & x$fill > 45, "white", INK)
  ggplot(x, aes(cause, pattern)) +
    geom_tile(aes(fill = fill), colour = "white", linewidth = 1.2) +
    geom_text(aes(label = lab, colour = txtcol), size = 2.85, lineheight = 0.9, nudge_y = 0.08) +
    geom_text(aes(label = nlab, colour = txtcol), size = 2.1, nudge_y = -0.34, alpha = 0.85) +
    scale_colour_identity() +
    scale_fill_gradientn(colours = c("#f3f7fd", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"),
                         values = rescale(sqrt(c(0, 2, 8, 20, 45, 90))), limits = c(0, 90),
                         trans = "identity", na.value = "#f0efec", name = sprintf("%d-month CIF (%%)", month),
                         guide = if (legend) guide_colourbar(barwidth = 10, barheight = 0.5, title.position = "top") else "none") +
    scale_x_discrete(position = "top") +
    labs(title = title, subtitle = if (!is.null(subtitle)) wrap(subtitle, 150),
         caption = if (!is.null(caption)) wrap(caption, 190)) + theme_mat()
}
## composition bars
plot_comp <- function(est_obj, month, row_order) {
  Q <- est_obj$Q[, , paste0("m", month)]
  x <- do.call(rbind, lapply(colnames(Q), function(cc) data.frame(pattern = rownames(Q), cause = cc, Q = Q[, cc])))
  x$pattern <- factor(x$pattern, rev(row_order))
  x$cause <- factor(x$cause, c("IH", "EH", "Sys", "NonHCC", "Undet"))
  x <- x[order(x$pattern, x$cause), ]
  x <- do.call(rbind, lapply(split(x, x$pattern), function(z) {
    z$xend <- cumsum(z$Q); z$xstart <- z$xend - z$Q; z$xmid <- z$xend - z$Q / 2; z }))
  x$y <- as.numeric(x$pattern)
  x$lab <- ifelse(x$Q >= 0.08, sprintf("%.0f", 100 * x$Q), "")
  ggplot(x) +
    geom_rect(aes(xmin = xstart, xmax = xend, ymin = y - 0.36, ymax = y + 0.36, fill = cause),
              colour = "white", linewidth = 0.4) +
    geom_text(aes(x = xmid, y = y, label = lab), size = 2.6, colour = "white", fontface = "bold") +
    scale_y_continuous(breaks = seq_along(levels(x$pattern)), labels = levels(x$pattern),
                       expand = expansion(add = 0.5)) +
    scale_fill_manual(values = COL_CAUSE, breaks = c("IH", "EH", "Sys", "NonHCC", "Undet"),
                      labels = CAUSE_LAB[c("IH", "EH", "Sys", "NonHCC", "Undet")], name = NULL) +
    scale_x_continuous(labels = percent, expand = c(0, 0), breaks = c(0, .5, 1), position = "top") +
    labs(title = sprintf("Cause composition at %d months", month), x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), axis.text.y = element_blank(), legend.position = "bottom",
          legend.text = element_text(size = 7.5), legend.key.size = unit(0.35, "cm"),
          plot.title = element_text(size = 10, face = "bold", colour = INK)) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
}
## HR matrix
plot_hr_matrix <- function(B, title, subtitle = NULL, caption = NULL, ref = REF_MAIN,
                           outcomes = c("All", "IH", "EH", "Sys", "NonHCC"), hr_lab = "HR") {
  x <- B$tab[B$tab$outcome %in% outcomes, ]
  x$pattern <- factor(as.character(x$pattern), rev(PAT_LEVELS))
  x$isref <- as.character(x$pattern) == ref
  x$lab <- ifelse(x$isref, "Reference", ifelse(x$NE, sprintf("NE\n(%d event%s)", x$events, ifelse(x$events == 1, "", "s")),
                  paste0(fmt_hr(x$HR, x$lo, x$hi), stars(x$p))))
  x$fill <- ifelse(x$isref | x$NE, NA, pmax(pmin(log2(x$HR), 2.5), -2.5))
  x$txtcol <- ifelse(!is.na(x$fill) & abs(x$fill) > 1.6, "white", INK)
  ov <- B$overall[B$overall$outcome %in% outcomes, ]
  ovl <- data.frame(outcome = ov$outcome, pattern = "Overall Wald P (df=8)",
                    lab = ifelse(ov$p < 0.001, "P < 0.001", paste0("P = ", fmt_p(ov$p))))
  lv_out <- outcomes
  x$outcome <- factor(x$outcome, lv_out, OUT_LAB[lv_out]); ovl$outcome <- factor(ovl$outcome, lv_out, OUT_LAB[lv_out])
  lv_rows <- c("Overall Wald P (df=8)", rev(PAT_LEVELS))
  x$pattern <- factor(as.character(x$pattern), lv_rows); ovl$pattern <- factor(ovl$pattern, lv_rows)
  ggplot(x, aes(outcome, pattern)) +
    geom_tile(aes(fill = fill), colour = "white", linewidth = 1.2) +
    geom_text(aes(label = lab, colour = txtcol), size = 2.85, lineheight = 0.9) +
    geom_text(data = ovl, aes(label = lab), size = 2.9, fontface = "bold", colour = INK) +
    geom_hline(yintercept = 1.5, colour = MUTED, linewidth = 0.4) +
    scale_colour_identity() +
    scale_fill_gradientn(colours = c("#104281", "#3987e5", "#9ec5f4", "#f0efec", "#f2a7a6", "#e34948", "#9e2322"),
                         values = rescale(c(-2.5, -1.5, -0.6, 0, 0.6, 1.5, 2.5)), limits = c(-2.5, 2.5),
                         na.value = "#f0efec", name = bquote(log[2](.(hr_lab))),
                         guide = guide_colourbar(barwidth = 10, barheight = 0.5, title.position = "top")) +
    scale_x_discrete(position = "top") + scale_y_discrete(limits = lv_rows) +
    labs(title = title, subtitle = if (!is.null(subtitle)) wrap(subtitle, 150),
         caption = if (!is.null(caption)) wrap(caption, 190)) + theme_mat()
}

row_order_A <- PAT_LEVELS[order(-est$F[, "All", paste0("m", TAU_MAIN)])]
lm_txt <- sprintf("Lunn–McNeil heterogeneity across IH/EH/Systemic HCC death: pattern × cause χ²(%d) = %.1f, P %s; site × cause χ²(%d) = %.1f, P %s; form × cause χ²(%d) = %.1f, P %s",
                  LM$df[1], LM$chisq[1], ifelse(LM$p[1] < 0.001, "< 0.001", paste("=", fmt_p(LM$p[1]))),
                  LM$df[2], LM$chisq[2], ifelse(LM$p[2] < 0.001, "< 0.001", paste("=", fmt_p(LM$p[2]))),
                  LM$df[3], LM$chisq[3], ifelse(LM$p[3] < 0.001, "< 0.001", paste("=", fmt_p(LM$p[3]))))
adj_note <- "Adjusted for TTP, age, sex, HBV, ECOG PS, ALBI grade, AFP ≥400, PIVKA-II ≥1000, tumour number >3, tumour size, Vp class, HVTT, EHS, cirrhosis; stratified by treatment group."

figA <- plot_cif_matrix(tab_F, TAU_MAIN,
          title = sprintf("A  Adjusted %d-month cumulative incidence of death after first progression", TAU_MAIN),
          subtitle = sprintf("Cause-specific Cox + G-computation (standardised to all %d progressors); %% (95%% bootstrap CI, B = %d); n = events within %d months",
                             nrow(d), length(boots), TAU_MAIN),
          row_order = row_order_A) +
  plot_comp(est, TAU_MAIN, row_order_A) +
  plot_layout(widths = c(5, 1.7)) +
  plot_annotation(caption = wrap(paste(adj_note, "Rows ordered by adjusted all-cause CIF. Cause-specific CIFs plus undetermined deaths sum to all-cause.",
                                  "Composition labels show % of deaths by 12 months (segments <8% unlabelled)."), 230),
                  theme = theme(plot.caption = element_text(size = 7.5, colour = INK2, hjust = 0)))
figB <- plot_hr_matrix(B_cs,
          title = "B  Post-progression death by first-progression pattern (reference: IH-NL)",
          subtitle = paste0("All-cause: Cox HR; cause-specific: csHR (cause-specific Cox). ", lm_txt),
          caption = paste(adj_note, "NE: <3 events in the cell. * P<0.05, ** P<0.01, *** P<0.001 (unadjusted; BH-FDR q values in Table S-HR).",
                          "Overall Wald test of the 8 pattern terms does not depend on the reference."))

## supplementary figures
figS_sHR <- plot_hr_matrix(B_fg, title = "Supplementary: subdistribution HR (Fine–Gray) by progression pattern (reference: IH-NL)",
                           subtitle = "Fine–Gray via survival::finegray with stratum-specific censoring weights; robust SE",
                           caption = adj_note, outcomes = c("IH", "EH", "Sys", "NonHCC"), hr_lab = "sHR")
figS_firth <- plot_hr_matrix(B_firth, title = "Supplementary: Firth-penalised HR / csHR (reference: IH-NL)",
                             subtitle = "coxphf (Firth-penalised partial likelihood, Wald CI from penalised information); treatment group entered as covariate (coxphf has no strata)",
                             caption = adj_note)
figS_undet <- plot_hr_matrix(S_undet_hr, title = "Supplementary: HR / csHR excluding undetermined deaths",
                             subtitle = sprintf("n = %d (%d deaths with unclassified/unknown/missing cause removed)", nrow(d_nound), sum(d$status == 5)), caption = adj_note)
figS_v0 <- plot_hr_matrix(S_v0_hr, title = "Supplementary: HR / csHR with the earlier pattern rule",
                          subtitle = sprintf("STG + TTU and EHG + TTU classified as Sys-TG instead of Sys-TG+NL (%d patients reclassified)",
                                             sum(d$pattern != d_v0$pattern)), caption = adj_note)
figS_del <- plot_hr_matrix(S_del_hr, title = "Supplementary: HR / csHR including patients flagged for deletion",
                           subtitle = sprintf("n = %d", nrow(d_del)), caption = adj_note)

plot_ard <- function(month) {
  x <- tab_ARD[tab_ARD$month == month & tab_ARD$cause %in% c("All", "IH", "EH", "Sys", "NonHCC"), ]
  x$pattern <- factor(x$pattern, rev(row_order_A))
  x$cause <- factor(x$cause, c("All", "IH", "EH", "Sys", "NonHCC"), OUT_LAB[c("All", "IH", "EH", "Sys", "NonHCC")])
  x$sig <- x$lo > 0 | x$hi < 0
  x$lab <- sprintf("%+.1f\n(%+.1f to %+.1f)", 100 * x$est, 100 * x$lo, 100 * x$hi)
  x$fill <- pmax(pmin(100 * x$est, 25), -25)
  ggplot(x, aes(cause, pattern)) +
    geom_tile(aes(fill = fill), colour = "white", linewidth = 1.2) +
    geom_tile(data = x[x$sig, ], fill = NA, colour = INK, linewidth = 0.5, width = 0.94, height = 0.9) +
    geom_text(aes(label = lab), size = 2.7, lineheight = 0.9, colour = INK) +
    scale_fill_gradientn(colours = c("#1c5cab", "#86b6ef", "#f0efec", "#f2a7a6", "#c62f2e"),
                         values = rescale(c(-25, -10, 0, 10, 25)), limits = c(-25, 25),
                         name = "ARD (percentage points)", guide = guide_colourbar(barwidth = 10, barheight = 0.5, title.position = "top")) +
    scale_x_discrete(position = "top") +
    labs(title = sprintf("Supplementary: absolute risk difference at %d months vs population-standardised average", month),
         subtitle = "ARD = standardised CIF of the pattern − weighted average over patterns (weights = observed pattern shares); boxed: 95% CI excludes 0",
         caption = wrap(adj_note, 190)) + theme_mat()
}
figS_ARD <- plot_ard(TAU_MAIN)

plot_pmr <- function(month) {
  x <- tab_PMR[tab_PMR$month == month & tab_PMR$cause %in% c("IH", "EH", "Sys", "NonHCC"), ]
  lab_map <- setNames(sprintf("%s  (P[Holm] = %s)", FU$pattern, fmt_p(FU$p_Holm)), FU$pattern)
  x$pattern_lab <- factor(lab_map[x$pattern], rev(lab_map[row_order_A]))
  x$cause <- factor(x$cause, c("IH", "EH", "Sys", "NonHCC"), CAUSE_LAB[c("IH", "EH", "Sys", "NonHCC")])
  x$zero <- x$events_by_tau == 0
  x$hi_c <- pmin(x$hi, 8); x$lo_c <- pmax(x$lo, 1 / 8); x$est_c <- pmin(pmax(x$est, 1 / 8), 8)
  ggplot(x[!x$zero, ], aes(est_c, pattern_lab, colour = cause)) +
    geom_vline(xintercept = 1, colour = MUTED, linewidth = 0.4) +
    geom_errorbarh(aes(xmin = lo_c, xmax = hi_c), height = 0, linewidth = 0.6) +
    geom_point(size = 2.2) +
    facet_wrap(~cause, nrow = 1) +
    scale_colour_manual(values = setNames(COL_CAUSE[c("IH", "EH", "Sys", "NonHCC")], CAUSE_LAB[c("IH", "EH", "Sys", "NonHCC")]), guide = "none") +
    scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4), labels = c("0.25", "0.5", "1", "2", "4")) +
    coord_cartesian(xlim = c(1 / 8, 8)) +
    labs(x = "CIF-based proportional mortality ratio (log scale)", y = NULL,
         title = sprintf("Supplementary: cause-of-death enrichment (PMR) at %d months", month),
         subtitle = "PMR = share of deaths from cause j within the pattern / share in all progressors. 95% bootstrap CI (clipped at 1/8 and 8).\nP[Holm]: weighted-effect-coded Wald test (df = 2) that the pattern's csHR is equal across IH/EH/Systemic HCC death (from Lunn–McNeil model).",
         caption = "PMR is a proportion metric, not a risk metric: interpret alongside the absolute CIF matrix; PMRs within a pattern are mutually constrained.") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 8.5, colour = INK2), plot.caption = element_text(size = 7.5, colour = INK2, hjust = 0),
          plot.background = element_rect(fill = "white", colour = NA))
}
figS_PMR <- plot_pmr(TAU_MAIN)
figS_AJ <- plot_cif_matrix(tab_AJ, TAU_MAIN,
             title = sprintf("Supplementary: unadjusted %d-month cumulative incidence (Aalen–Johansen)", TAU_MAIN),
             subtitle = "% (95% CI); all-cause = 1 − Kaplan–Meier; n = events within 12 months", row_order = row_order_A)
figS_F6  <- plot_cif_matrix(tab_F, 6, title = "Supplementary: adjusted 6-month CIF after first progression",
                            subtitle = "Same models as Figure 4.4A; n = events within 6 months", row_order = row_order_A)
figS_F24 <- plot_cif_matrix(tab_F, 24, title = "Supplementary: adjusted 24-month CIF after first progression",
                            subtitle = "Same models as Figure 4.4A; n = events within 24 months", row_order = row_order_A)

save_fig <- function(p, name, w, h) {
  ggsave(file.path(OUT_DIR, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(OUT_DIR, "png", paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
}
log_msg("Saving figures ...")
save_fig(figA, "Fig4.4A_adjusted_CIF_matrix", 13, 7.2)
save_fig(figB, "Fig4.4B_csHR_matrix_refIHNL", 11.5, 8)
supp <- list(figS_sHR = figS_sHR, figS_ARD = figS_ARD, figS_PMR = figS_PMR, figS_AJ = figS_AJ,
             figS_F6 = figS_F6, figS_F24 = figS_F24, figS_firth = figS_firth,
             figS_undet = figS_undet, figS_del = figS_del, figS_v0 = figS_v0)
supp_dims <- list(c(10.5, 8), c(11, 7.5), c(12, 6), c(11, 7.5), c(11, 7.5), c(11, 7.5), c(11.5, 8), c(11.5, 8), c(11.5, 8), c(11.5, 8))
cairo_pdf(file.path(OUT_DIR, "Fig4.4_Supplementary.pdf"), width = 12, height = 8.5, onefile = TRUE)
for (p in supp) print(p)
dev.off()
for (i in seq_along(supp)) ggsave(file.path(OUT_DIR, "png", paste0(names(supp)[i], ".png")), supp[[i]],
                                  width = supp_dims[[i]][1], height = supp_dims[[i]][2], dpi = 200, bg = "white")

## =============================================================================
## 10. Tables (xlsx) and a combined PDF report
## =============================================================================
r3 <- function(z) { z[] <- lapply(z, function(v) if (is.numeric(v)) signif(v, 4) else v); z }
hr_tab_fmt <- function(B) { z <- B$tab; z$HR_CI <- ifelse(z$NE, "NE", ifelse(is.na(z$lo), "Reference", sprintf("%.2f (%.2f-%.2f)", z$HR, z$lo, z$hi))); r3(z) }
sheets <- list(
  QC_linkage = qc, Cohort_flow = flow,
  Events_total = events_all, Events_by_12m = events_tau$m12,
  N_at_risk = risk_tab,
  CIF_adjusted = r3(tab_F), CIF_pop_average = r3(tab_Fbar), CIF_unadjusted_AJ = r3(tab_AJ),
  ARD = r3(tab_ARD), Composition_Q = r3(tab_Q), PMR = r3(tab_PMR),
  HR_csHR_refIHNL = hr_tab_fmt(B_cs), sHR_FineGray = hr_tab_fmt(B_fg), HR_Firth = hr_tab_fmt(B_firth),
  Overall_Wald = r3(rbind(B_cs$overall, B_fg$overall, B_firth$overall)),
  LunnMcNeil = r3(LM), LM_followup_effectcoded = r3(FU),
  PH_Schoenfeld = r3(B_cs$zph),
  Validation_refIHTG = VAL,
  Sens_excl_undet_HR = hr_tab_fmt(S_undet_hr),
  Sens_excl_undet_CIF12 = r3(data.frame(pattern = PAT_LEVELS, S_undet_cif$F[, , "m12"], check.names = FALSE)),
  Sens_incl_flagged_HR = hr_tab_fmt(S_del_hr),
  Sens_earlier_rule_HR = hr_tab_fmt(S_v0_hr),
  Pattern_classification = pattern_check,
  Interaction_treatment = r3(INTER),
  Run_info = data.frame(item = c("date", "R", "survival", "riskRegression", "coxphf", "bootstrap B (successful)",
                                 "seed", "engine vs predictRisk max diff", "LM m1 vs separate fit max diff",
                                 "median post-progression follow-up (months, reverse KM)"),
                        value = c(as.character(Sys.time()), R.version.string, as.character(packageVersion("survival")),
                                  as.character(packageVersion("riskRegression")), as.character(packageVersion("coxphf")),
                                  length(boots), SEED, signif(engine_check, 3), signif(chk_lm, 3),
                                  sprintf("%.1f (%.1f-%.1f)", fu_median[1], fu_median[2], fu_median[3]))))
wb <- createWorkbook()
for (s in names(sheets)) { addWorksheet(wb, s); writeData(wb, s, sheets[[s]]); freezePane(wb, s, firstRow = TRUE)
  setColWidths(wb, s, cols = seq_len(ncol(sheets[[s]])), widths = "auto") }
saveWorkbook(wb, file.path(OUT_DIR, "Fig4.4_Tables.xlsx"), overwrite = TRUE)
saveRDS(list(d = d, est = est, ci = ci, boots_n = length(boots), B_cs = B_cs, B_fg = B_fg, B_firth = B_firth,
             LM = LM, FU = FU, VAL = VAL, INTER = INTER), file.path(OUT_DIR, "Fig4.4_results.rds"))

## combined report PDF ---------------------------------------------------------
table_page <- function(df, title, note = NULL, base = 7) {
  base <- base * 1.35
  tg <- tableGrob(df, rows = NULL, theme = ttheme_minimal(base_size = base,
         core = list(fg_params = list(hjust = 0, x = 0.03)), colhead = list(fg_params = list(hjust = 0, x = 0.03, fontface = "bold"))))
  grid.newpage()
  grid.text(title, x = 0.03, y = 0.96, hjust = 0, gp = gpar(fontsize = 13, fontface = "bold"))
  if (!is.null(note)) grid.text(note, x = 0.03, y = 0.925, hjust = 0, vjust = 1, gp = gpar(fontsize = 8, col = INK2))
  tg$vp <- viewport(x = 0.03, y = if (is.null(note)) 0.91 else 0.87, just = c("left", "top"),
                    width = sum(tg$widths), height = sum(tg$heights))
  grid.draw(tg)
}
ev_show <- events_all; names(ev_show) <- c("Pattern", "N", "IH death", "EH death", "Sys death", "Non-HCC", "Undetermined", "All deaths")
risk_show <- risk_tab; names(risk_show) <- c("Pattern", paste("At risk", TAUS, "mo"))
ev_show <- merge(ev_show, risk_show, by.x = "Pattern", by.y = "Pattern", sort = FALSE)
ev_show <- ev_show[match(PAT_LEVELS, ev_show$Pattern), ]
lm_show <- transform(LM, chisq = sprintf("%.2f", chisq), p = fmt_p(p))
fu_show <- data.frame(Pattern = FU$pattern,
  `IH effect-coded HR` = sprintf("%.2f (%.2f-%.2f)", FU$eHR_IH, FU$eHR_IH_lo, FU$eHR_IH_hi),
  `EH effect-coded HR` = sprintf("%.2f (%.2f-%.2f)", FU$eHR_EH, FU$eHR_EH_lo, FU$eHR_EH_hi),
  `Sys effect-coded HR` = sprintf("%.2f (%.2f-%.2f)", FU$eHR_Sys, FU$eHR_Sys_lo, FU$eHR_Sys_hi),
  `Wald chi2 (df=2)` = sprintf("%.2f", FU$chisq), P = fmt_p(FU$p), `P Holm` = fmt_p(FU$p_Holm), check.names = FALSE)
zph_show <- B_cs$zph[B_cs$zph$term %in% c("pattern", "GLOBAL"), ]
zph_show <- transform(zph_show, chisq = sprintf("%.2f", chisq), p = fmt_p(p))
ov_show <- rbind(B_cs$overall, B_fg$overall, B_firth$overall); ov_show <- transform(ov_show, chisq = sprintf("%.2f", chisq), p = fmt_p(p))
inter_show <- transform(INTER, chisq = sprintf("%.2f", chisq), p = fmt_p(p))
qc_show <- qc; qc_show$check <- gsub("[^ -~]", "?", qc_show$check)   # ASCII-safe

cairo_pdf(file.path(OUT_DIR, "Fig4.4_Full_Report.pdf"), width = 13, height = 8.5, onefile = TRUE)
grid.newpage()
grid.text("Figure 4.4  Progression pattern and post-progression cause of death", 0.05, 0.85, hjust = 0, gp = gpar(fontsize = 18, fontface = "bold"))
grid.text(paste(c(
  sprintf("Analysis cohort: %d patients with first radiological progression (HAIC-Sys %d, TACE-Sys %d, Systemic %d); %d deaths after progression.",
          nrow(d), sum(d$Group == "HAIC_Systemic"), sum(d$Group == "TACE_Systemic"), sum(d$Group == "Systemic"), sum(d$status > 0)),
  sprintf("Deaths by cause: intrahepatic HCC %d, extrahepatic HCC %d, systemic HCC %d, non-HCC %d, undetermined %d.",
          sum(d$status == 1), sum(d$status == 2), sum(d$status == 3), sum(d$status == 4), sum(d$status == 5)),
  sprintf("Median post-progression follow-up (reverse KM): %.1f months. Time zero = date of first progression; main time point %d months.", fu_median[1], TAU_MAIN),
  sprintf("Bootstrap: %d successful patient-level resamples (all models refitted in each). Engine check vs riskRegression: max diff %.1e.", length(boots), engine_check),
  "", "Contents:", "  1. Figure 4.4A  adjusted 12-month CIF matrix + cause composition",
  "  2. Figure 4.4B  HR / csHR matrix (reference IH-NL) + overall Wald + Lunn-McNeil",
  "  3. Cohort flow, events and numbers at risk;  4. Lunn-McNeil and follow-up tests;  5. Overall Wald tests, PH tests",
  "  6. Validation vs original IH-TG matrix (incl. earlier pattern rule);  7. Treatment interaction;  8+. Supplementary figures (sHR, ARD, PMR, unadjusted AJ, 6/24 months, Firth, sensitivity)"),
  collapse = "\n"), 0.05, 0.72, hjust = 0, vjust = 1, gp = gpar(fontsize = 10.5, lineheight = 1.25))
print(figA); print(figB)
table_page(flow, "Cohort flow")
pc_show <- pattern_check
pc_show$PFS_model_components <- vapply(strsplit(pc_show$PFS_model_components, ";"), function(t)
  paste(c("肝内病灶增大" = "IHG", "癌栓升级" = "TTU", "肝内新发病灶" = "NIH", "新发血管侵犯" = "nMVI",
          "肝外病灶增大" = "EHG", "肝外新发转移" = "NEH", "系统性肿瘤增大" = "STG")[t], collapse = " + "), "")
pc_show$changed <- ifelse(pc_show$pattern != pc_show$earlier_rule, "<- reclassified", "")
names(pc_show) <- c("Components (TTU->IHG, nMVI->NIH)", "Pattern", "Earlier rule", "n", "")
for (ch in split(seq_len(nrow(pc_show)), ceiling(seq_len(nrow(pc_show)) / 24)))
table_page(pc_show[ch, ], "Progression-pattern classification by component combination (analysis cohort)", base = 6,
           note = "(1) IH-TG only IHG; (2) IH-NL only NIH; (3) IH-TG+NL IHG+NIH; (4) EH-TG only EHG; (5) EH-NL only NEH; (6) EH-TG+NL EHG+NEH; (7) Sys-TG only STG; (8) Sys-NL NIH+NEH; (9) Sys-TG+NL all other combinations")
table_page(ev_show, "Deaths after first progression by pattern and cause (whole follow-up) and numbers at risk",
           note = "IH/EH/Sys = intrahepatic / extrahepatic / systemic HCC-related death; Non-HCC = liver-related + other + treatment-related; Undetermined = unclassified HCC-related + unknown + missing cause", base = 8)
table_page(lm_show, "Lunn-McNeil cross-cause heterogeneity tests (likelihood ratio)", base = 9,
           note = "Expanded data (one row per patient per HCC-related cause); non-HCC and undetermined deaths censored; strata = cause x treatment group; covariates interact with cause")
table_page(fu_show, "Follow-up: pattern effect relative to weighted average pattern (weighted effect coding), per cause", base = 8,
           note = "Wald test (df = 2) that the effect-coded csHR is equal across IH/EH/Systemic HCC death; Holm-adjusted across 9 patterns. Descriptive counterpart: PMR (supplementary).")
table_page(ov_show, "Overall Wald tests of the 8 pattern terms (independent of reference group)", base = 9)
table_page(zph_show, "Proportional hazards (Schoenfeld residuals; cs-Cox / all-cause Cox)", base = 9,
           note = "Where the pattern term violates PH, the HR is a follow-up-averaged effect; the fixed-time CIF/ARD (4.4A) do not rely on PH for pattern and are the primary absolute-risk summaries.")
table_page(VAL, "Validation against the original IH-TG-referenced matrix", base = 8,
           note = "Target values from the original PDF: Sys-TG all-cause HR 1.96 (1.48-2.61), Sys-TG csHR systemic death 4.11, IH-NL all-cause HR 0.76")
table_page(inter_show, "Sensitivity: interaction with treatment group (justifies pooled analysis)", base = 9)
table_page(qc_show, "Data linkage checks", base = 9)
for (p in supp) print(p)
dev.off()
log_msg("All outputs written to ", normalizePath(OUT_DIR))
