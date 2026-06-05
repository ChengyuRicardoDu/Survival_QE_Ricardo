library(survival)
library(splines)
library(rstpm2)
library(ggplot2)
library(patchwork)

# Checks used for the report:
# functional form in the clinical Cox model,
# PH diagnostics for Cox models,
# AIC within each model family,
# Cox-Snell residuals for Cox/RP/AFT models,
# deviance residuals for the Cox models.

dat <- readRDS("output/dat_clean.rds")

dat$nyha_nominal <- factor(
  dat$nyha_class,
  levels = c("I", "II", "III", "IV"),
  ordered = FALSE
)

dat$nyha_II <- as.numeric(dat$nyha_nominal == "II")
dat$nyha_III <- as.numeric(dat$nyha_nominal == "III")
dat$nyha_IV <- as.numeric(dat$nyha_nominal == "IV")

dir.create("tables/goodness_of_fit", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/goodness_of_fit", recursive = TRUE, showWarnings = FALSE)

surv_obj <- Surv(dat$time, dat$status)

cox_clinical <- coxph(
  Surv(time, status) ~ nyha_nominal +
    age + gender + bmi + smoking + diabetes +
    ejection_fraction + ns(log(creatinine), df = 2) + sodium,
  data = dat
)

cox_gene <- coxph(
  Surv(time, status) ~ nyha_nominal + age + gender + bmi +
    smoking + diabetes + ejection_fraction + ns(log(creatinine), df = 2) + sodium +
    GENE_001 + GENE_002 + GENE_005 + GENE_024,
  data = dat
)

cox_clinical_summary <- summary(cox_clinical)
cox_coef <- as.data.frame(cox_clinical_summary$coefficients)
cox_ci <- as.data.frame(cox_clinical_summary$conf.int)

cox_results <- data.frame(
  term = rownames(cox_coef),
  coef = cox_coef$coef,
  hazard_ratio = cox_ci$`exp(coef)`,
  lower_95 = cox_ci$`lower .95`,
  upper_95 = cox_ci$`upper .95`,
  p_value = cox_coef$`Pr(>|z|)`
)

write.csv(
  cox_results,
  "tables/goodness_of_fit/clinical_cox_results.csv",
  row.names = FALSE
)

cox_gene_summary <- summary(cox_gene)
cox_gene_coef <- as.data.frame(cox_gene_summary$coefficients)
cox_gene_ci <- as.data.frame(cox_gene_summary$conf.int)

cox_gene_results <- data.frame(
  term = rownames(cox_gene_coef),
  coef = cox_gene_coef$coef,
  hazard_ratio = cox_gene_ci$`exp(coef)`,
  lower_95 = cox_gene_ci$`lower .95`,
  upper_95 = cox_gene_ci$`upper .95`,
  p_value = cox_gene_coef$`Pr(>|z|)`
)

write.csv(
  cox_gene_results,
  "tables/goodness_of_fit/selected_gene_cox_results.csv",
  row.names = FALSE
)

cox_linear <- coxph(
  Surv(time, status) ~ nyha_nominal +
    age + gender + bmi + smoking + diabetes +
    ejection_fraction + log(creatinine) + sodium,
  data = dat
)

cox_age_df2 <- coxph(
  Surv(time, status) ~ nyha_nominal +
    ns(age, df = 2) + gender + bmi + smoking + diabetes +
    ejection_fraction + log(creatinine) + sodium,
  data = dat
)

cox_bmi_df2 <- coxph(
  Surv(time, status) ~ nyha_nominal +
    age + gender + ns(bmi, df = 2) + smoking + diabetes +
    ejection_fraction + log(creatinine) + sodium,
  data = dat
)

cox_ef_df2 <- coxph(
  Surv(time, status) ~ nyha_nominal +
    age + gender + bmi + smoking + diabetes +
    ns(ejection_fraction, df = 2) + log(creatinine) + sodium,
  data = dat
)

cox_creatinine_df2 <- coxph(
  Surv(time, status) ~ nyha_nominal +
    age + gender + bmi + smoking + diabetes +
    ejection_fraction + ns(log(creatinine), df = 2) + sodium,
  data = dat
)

cox_sodium_df2 <- coxph(
  Surv(time, status) ~ nyha_nominal +
    age + gender + bmi + smoking + diabetes +
    ejection_fraction + log(creatinine) + ns(sodium, df = 2),
  data = dat
)

cox_creatinine_df3 <- coxph(
  Surv(time, status) ~ nyha_nominal +
    age + gender + bmi + smoking + diabetes +
    ejection_fraction + ns(log(creatinine), df = 3) + sodium,
  data = dat
)

functional_form <- data.frame(
  variable = c("age", "bmi", "ejection_fraction", "log(creatinine)", "sodium", "log(creatinine) df3 vs df2"),
  comparison = c(
    "linear vs df2 spline",
    "linear vs df2 spline",
    "linear vs df2 spline",
    "linear vs df2 spline",
    "linear vs df2 spline",
    "df2 spline vs df3 spline"
  ),
  AIC_simple = c(
    AIC(cox_linear),
    AIC(cox_linear),
    AIC(cox_linear),
    AIC(cox_linear),
    AIC(cox_linear),
    AIC(cox_creatinine_df2)
  ),
  AIC_flexible = c(
    AIC(cox_age_df2),
    AIC(cox_bmi_df2),
    AIC(cox_ef_df2),
    AIC(cox_creatinine_df2),
    AIC(cox_sodium_df2),
    AIC(cox_creatinine_df3)
  ),
  LR_p_value = c(
    anova(cox_linear, cox_age_df2, test = "LRT")[2, "Pr(>|Chi|)"],
    anova(cox_linear, cox_bmi_df2, test = "LRT")[2, "Pr(>|Chi|)"],
    anova(cox_linear, cox_ef_df2, test = "LRT")[2, "Pr(>|Chi|)"],
    anova(cox_linear, cox_creatinine_df2, test = "LRT")[2, "Pr(>|Chi|)"],
    anova(cox_linear, cox_sodium_df2, test = "LRT")[2, "Pr(>|Chi|)"],
    anova(cox_creatinine_df2, cox_creatinine_df3, test = "LRT")[2, "Pr(>|Chi|)"]
  )
)

write.csv(
  functional_form,
  "tables/goodness_of_fit/functional_form_checks.csv",
  row.names = FALSE
)

ph_clinical <- cox.zph(cox_clinical, terms = FALSE)
ph_gene <- cox.zph(cox_gene, terms = FALSE)

ph_tests <- rbind(
  data.frame(model = "Clinical Cox", variable = rownames(ph_clinical$table), ph_clinical$table),
  data.frame(model = "Selected-gene Cox", variable = rownames(ph_gene$table), ph_gene$table)
)

write.csv(
  ph_tests,
  "tables/goodness_of_fit/ph_tests.csv",
  row.names = FALSE
)

cox_tv_t <- coxph(
  Surv(time, status) ~ nyha_II + nyha_III + nyha_IV +
    age + gender + bmi + smoking + diabetes +
    ejection_fraction + ns(log(creatinine), df = 2) + sodium +
    tt(nyha_II) + tt(nyha_III) + tt(nyha_IV),
  data = dat,
  tt = function(x, t, ...) x * t
)

cox_tv_logt <- coxph(
  Surv(time, status) ~ nyha_II + nyha_III + nyha_IV +
    age + gender + bmi + smoking + diabetes +
    ejection_fraction + ns(log(creatinine), df = 2) + sodium +
    tt(nyha_II) + tt(nyha_III) + tt(nyha_IV),
  data = dat,
  tt = function(x, t, ...) x * log(t)
)

cox_tv_log1t <- coxph(
  Surv(time, status) ~ nyha_II + nyha_III + nyha_IV +
    age + gender + bmi + smoking + diabetes +
    ejection_fraction + ns(log(creatinine), df = 2) + sodium +
    tt(nyha_II) + tt(nyha_III) + tt(nyha_IV),
  data = dat,
  tt = function(x, t, ...) x * log1p(t)
)

tv_compare <- data.frame(
  model = c("NYHA x t", "NYHA x log(t)", "NYHA x log(1+t)"),
  logLik = c(
    as.numeric(logLik(cox_tv_t)),
    as.numeric(logLik(cox_tv_logt)),
    as.numeric(logLik(cox_tv_log1t))
  ),
  LR = c(
    2 * (as.numeric(logLik(cox_tv_t)) - as.numeric(logLik(cox_clinical))),
    2 * (as.numeric(logLik(cox_tv_logt)) - as.numeric(logLik(cox_clinical))),
    2 * (as.numeric(logLik(cox_tv_log1t)) - as.numeric(logLik(cox_clinical)))
  ),
  df = 3,
  AIC = c(
    AIC(cox_tv_t),
    AIC(cox_tv_logt),
    AIC(cox_tv_log1t)
  )
)

tv_compare$p_value <- pchisq(tv_compare$LR, df = tv_compare$df, lower.tail = FALSE)

write.csv(
  tv_compare,
  "tables/goodness_of_fit/time_varying_cox_comparison.csv",
  row.names = FALSE
)

rp_base <- stpm2(
  Surv(time, status) ~ nyha_II + nyha_III + nyha_IV +
    age + gender + bmi + smoking + diabetes +
    ejection_fraction + ns(log(creatinine), df = 2) + sodium,
  data = dat,
  df = 4
)

rp_tvc1 <- stpm2(
  Surv(time, status) ~ nyha_II + nyha_III + nyha_IV +
    age + gender + bmi + smoking + diabetes +
    ejection_fraction + ns(log(creatinine), df = 2) + sodium,
  data = dat,
  df = 4,
  tvc = list(nyha_II = 1, nyha_III = 1, nyha_IV = 1)
)

rp_tvc2 <- stpm2(
  Surv(time, status) ~ nyha_II + nyha_III + nyha_IV +
    age + gender + bmi + smoking + diabetes +
    ejection_fraction + ns(log(creatinine), df = 2) + sodium,
  data = dat,
  df = 4,
  tvc = list(nyha_II = 2, nyha_III = 2, nyha_IV = 2)
)

rp_tvc3 <- stpm2(
  Surv(time, status) ~ nyha_II + nyha_III + nyha_IV +
    age + gender + bmi + smoking + diabetes +
    ejection_fraction + ns(log(creatinine), df = 2) + sodium,
  data = dat,
  df = 4,
  tvc = list(nyha_II = 3, nyha_III = 3, nyha_IV = 3)
)

rp_tvc_gene <- stpm2(
  Surv(time, status) ~ nyha_II + nyha_III + nyha_IV +
    age + gender + bmi + smoking + diabetes +
    ejection_fraction + ns(log(creatinine), df = 2) + sodium +
    GENE_001 + GENE_002 + GENE_005 + GENE_024,
  data = dat,
  df = 4,
  tvc = list(nyha_II = 1, nyha_III = 1, nyha_IV = 1)
)

rp_compare <- data.frame(
  model = c("RP time-fixed", "RP tvc df1", "RP tvc df2", "RP tvc df3", "RP tvc df1 + selected genes"),
  AIC = c(
    AIC(rp_base),
    AIC(rp_tvc1),
    AIC(rp_tvc2),
    AIC(rp_tvc3),
    AIC(rp_tvc_gene)
  ),
  df = c(18, 21, 24, 27, 25)
)

write.csv(
  rp_compare,
  "tables/goodness_of_fit/rp_model_comparison.csv",
  row.names = FALSE
)

aft_weibull <- survreg(
  surv_obj ~ nyha_nominal +
    age + gender + bmi + smoking + diabetes +
    ejection_fraction + ns(log(creatinine), df = 2) + sodium,
  data = dat,
  dist = "weibull"
)

aft_exponential <- survreg(
  surv_obj ~ nyha_nominal +
    age + gender + bmi + smoking + diabetes +
    ejection_fraction + ns(log(creatinine), df = 2) + sodium,
  data = dat,
  dist = "exponential"
)

aft_lognormal <- survreg(
  surv_obj ~ nyha_nominal +
    age + gender + bmi + smoking + diabetes +
    ejection_fraction + ns(log(creatinine), df = 2) + sodium,
  data = dat,
  dist = "lognormal"
)

aft_loglogistic <- survreg(
  surv_obj ~ nyha_nominal +
    age + gender + bmi + smoking + diabetes +
    ejection_fraction + ns(log(creatinine), df = 2) + sodium,
  data = dat,
  dist = "loglogistic"
)

aft_coef <- as.data.frame(summary(aft_lognormal)$table)

aft_time_ratio <- data.frame(
  term = rownames(aft_coef),
  coef = aft_coef$Value,
  se = aft_coef$`Std. Error`,
  z = aft_coef$z,
  p_value = aft_coef$p,
  time_ratio = exp(aft_coef$Value),
  lower_95 = exp(aft_coef$Value - 1.96 * aft_coef$`Std. Error`),
  upper_95 = exp(aft_coef$Value + 1.96 * aft_coef$`Std. Error`)
)

write.csv(
  aft_time_ratio,
  "tables/goodness_of_fit/lognormal_aft_time_ratios.csv",
  row.names = FALSE
)

model_fit <- data.frame(
  family = c("Cox", "Cox", "Royston-Parmar", "Royston-Parmar", "AFT", "AFT", "AFT", "AFT"),
  model = c(
    "Clinical Cox",
    "Selected-gene Cox",
    "RP time-fixed",
    "RP time-varying NYHA",
    "Weibull",
    "Exponential",
    "Log-normal",
    "Log-logistic"
  ),
  df = c(
    attr(logLik(cox_clinical), "df"),
    attr(logLik(cox_gene), "df"),
    18,
    21,
    attr(logLik(aft_weibull), "df"),
    attr(logLik(aft_exponential), "df"),
    attr(logLik(aft_lognormal), "df"),
    attr(logLik(aft_loglogistic), "df")
  ),
  logLik = c(
    as.numeric(logLik(cox_clinical)),
    as.numeric(logLik(cox_gene)),
    NA,
    NA,
    as.numeric(logLik(aft_weibull)),
    as.numeric(logLik(aft_exponential)),
    as.numeric(logLik(aft_lognormal)),
    as.numeric(logLik(aft_loglogistic))
  ),
  AIC = c(
    AIC(cox_clinical),
    AIC(cox_gene),
    AIC(rp_base),
    AIC(rp_tvc1),
    AIC(aft_weibull),
    AIC(aft_exponential),
    AIC(aft_lognormal),
    AIC(aft_loglogistic)
  ),
  C_index = c(
    summary(cox_clinical)$concordance[1],
    summary(cox_gene)$concordance[1],
    NA, NA, NA, NA, NA, NA
  )
)

write.csv(
  model_fit,
  "tables/goodness_of_fit/model_fit_summary.csv",
  row.names = FALSE
)

cs_curve <- function(cs, status, model_name) {
  cs <- pmax(cs, 1e-8)
  cs_fit <- survfit(Surv(cs, status) ~ 1)
  data.frame(
    residual = cs_fit$time,
    cumulative_hazard = -log(pmax(cs_fit$surv, 1e-8)),
    model = model_name
  )
}

cs_clinical <- dat$status - residuals(cox_clinical, type = "martingale")
cs_gene <- dat$status - residuals(cox_gene, type = "martingale")

rp_surv <- predict(rp_tvc1, newdata = dat, type = "surv")
rp_surv <- as.numeric(rp_surv)
rp_surv <- pmin(pmax(rp_surv, 1e-8), 1 - 1e-8)
cs_rp <- -log(rp_surv)

lp_weibull <- predict(aft_weibull, type = "lp")
surv_weibull <- exp(-exp((log(dat$time) - lp_weibull) / aft_weibull$scale))

lp_exponential <- predict(aft_exponential, type = "lp")
surv_exponential <- exp(-exp(log(dat$time) - lp_exponential))

lp_lognormal <- predict(aft_lognormal, type = "lp")
surv_lognormal <- 1 - pnorm((log(dat$time) - lp_lognormal) / aft_lognormal$scale)

lp_loglogistic <- predict(aft_loglogistic, type = "lp")
surv_loglogistic <- 1 / (1 + exp((log(dat$time) - lp_loglogistic) / aft_loglogistic$scale))

cs_aft_weibull <- -log(pmin(pmax(surv_weibull, 1e-8), 1 - 1e-8))
cs_aft_exponential <- -log(pmin(pmax(surv_exponential, 1e-8), 1 - 1e-8))
cs_aft_lognormal <- -log(pmin(pmax(surv_lognormal, 1e-8), 1 - 1e-8))
cs_aft_loglogistic <- -log(pmin(pmax(surv_loglogistic, 1e-8), 1 - 1e-8))

cs_dat <- rbind(
  cs_curve(cs_clinical, dat$status, "Clinical Cox"),
  cs_curve(cs_gene, dat$status, "Selected-gene Cox"),
  cs_curve(cs_rp, dat$status, "RP time-varying NYHA"),
  cs_curve(cs_aft_weibull, dat$status, "AFT Weibull"),
  cs_curve(cs_aft_exponential, dat$status, "AFT Exponential"),
  cs_curve(cs_aft_lognormal, dat$status, "AFT Log-normal"),
  cs_curve(cs_aft_loglogistic, dat$status, "AFT Log-logistic")
)

cs_main <- cs_dat[
  cs_dat$residual <= quantile(cs_dat$residual, 0.98, na.rm = TRUE),
]

cs_summary <- aggregate(
  abs(cumulative_hazard - residual) ~ model,
  data = cs_main,
  FUN = function(x) c(mean = mean(x), max = max(x))
)

cs_summary <- data.frame(
  model = cs_summary$model,
  mean_abs_difference = cs_summary[, 2][, "mean"],
  max_abs_difference = cs_summary[, 2][, "max"]
)

write.csv(
  cs_summary,
  "tables/goodness_of_fit/cox_snell_summary.csv",
  row.names = FALSE
)

cs_plot <- ggplot(cs_dat, aes(x = residual, y = cumulative_hazard)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.35, colour = "grey45") +
  geom_step(linewidth = 0.55, colour = "black") +
  facet_wrap(~ model, scales = "free", ncol = 3) +
  coord_cartesian(xlim = c(0, quantile(cs_dat$residual, 0.98, na.rm = TRUE))) +
  labs(
    x = "Cox-Snell residual",
    y = "Nelson-Aalen cumulative hazard"
  ) +
  theme_classic(base_size = 8) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))

ggsave(
  "figures/goodness_of_fit/cox_snell_residuals.pdf",
  cs_plot,
  width = 8.5,
  height = 6.3
)

ggsave(
  "figures/goodness_of_fit/cox_snell_residuals.png",
  cs_plot,
  width = 8.5,
  height = 6.3,
  dpi = 300
)

dev_dat <- rbind(
  data.frame(
    linear_predictor = predict(cox_clinical, type = "lp"),
    residual = residuals(cox_clinical, type = "deviance"),
    model = "Clinical Cox"
  ),
  data.frame(
    linear_predictor = predict(cox_gene, type = "lp"),
    residual = residuals(cox_gene, type = "deviance"),
    model = "Selected-gene Cox"
  )
)

write.csv(
  dev_dat,
  "tables/goodness_of_fit/deviance_residuals.csv",
  row.names = FALSE
)

deviance_summary <- aggregate(
  residual ~ model,
  data = dev_dat,
  FUN = function(x) c(
    min = min(x),
    median = median(x),
    max = max(x),
    n_abs_gt_2 = sum(abs(x) > 2),
    n_abs_gt_3 = sum(abs(x) > 3)
  )
)

deviance_summary <- data.frame(
  model = deviance_summary$model,
  min = deviance_summary[, 2][, "min"],
  median = deviance_summary[, 2][, "median"],
  max = deviance_summary[, 2][, "max"],
  n_abs_gt_2 = deviance_summary[, 2][, "n_abs_gt_2"],
  n_abs_gt_3 = deviance_summary[, 2][, "n_abs_gt_3"]
)

write.csv(
  deviance_summary,
  "tables/goodness_of_fit/deviance_residual_summary.csv",
  row.names = FALSE
)

dev_plot <- ggplot(dev_dat, aes(x = linear_predictor, y = residual)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_hline(yintercept = c(-2, 2), linetype = "dotted", linewidth = 0.3, colour = "grey50") +
  geom_point(size = 0.8, alpha = 0.55) +
  facet_wrap(~ model, scales = "free_x") +
  labs(
    x = "Linear predictor",
    y = "Deviance residual"
  ) +
  theme_classic(base_size = 8) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))

ggsave(
  "figures/goodness_of_fit/deviance_residuals.pdf",
  dev_plot,
  width = 7.2,
  height = 3.3
)

ggsave(
  "figures/goodness_of_fit/deviance_residuals.png",
  dev_plot,
  width = 7.2,
  height = 3.3,
  dpi = 300
)

print(functional_form)
print(ph_tests)
print(model_fit)
print(cs_summary)
print(deviance_summary)

q(save = "no", status = 0, runLast = FALSE)
