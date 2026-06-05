library(survival)
library(splines)
library(ggplot2)
library(survminer)
library(patchwork)
library(rstpm2)

dat <- readRDS("output/dat_clean.rds")

dat$nyha_nominal <- factor(
  dat$nyha_class,
  levels = c("I", "II", "III", "IV"),
  ordered = FALSE
)

dat$nyha_II <- as.numeric(dat$nyha_nominal == "II")
dat$nyha_III <- as.numeric(dat$nyha_nominal == "III")
dat$nyha_IV <- as.numeric(dat$nyha_nominal == "IV")

dir.create("figures", showWarnings = FALSE)

cols_nyha <- c("#4C78A8", "#72A87E", "#E3A35C", "#C75B5B")
cols_contrast <- c("#72A87E", "#E3A35C", "#C95F5F")

surv_obj <- Surv(time = dat$time, event = dat$status)

km1 <- survfit(surv_obj ~ 1, data = dat)
km2 <- survfit(surv_obj ~ nyha_class, data = dat)

km_overall <- ggsurvplot(
  km1,
  data = dat,
  conf.int = FALSE,
  censor = TRUE,
  censor.size = 1.4,
  risk.table = TRUE,
  risk.table.fontsize = 3.0,
  risk.table.title = "Number at risk",
  risk.table.y.text = FALSE,
  xlab = "Time (months)",
  ylab = "Survival probability",
  xlim = c(0, 72),
  break.time.by = 12,
  palette = "black",
  ggtheme = theme_classic(base_size = 9)
)

km_nyha <- ggsurvplot(
  km2,
  data = dat,
  conf.int = FALSE,
  censor = TRUE,
  censor.size = 1.4,
  risk.table = TRUE,
  risk.table.fontsize = 3.0,
  risk.table.title = "Number at risk",
  risk.table.y.text.col = TRUE,
  pval = TRUE,
  pval.method = TRUE,
  pval.size = 3.2,
  pval.method.size = 3.2,
  pval.coord = c(36, 0.92),
  pval.method.coord = c(36, 0.86),
  palette = cols_nyha,
  legend.title = "NYHA class",
  legend.labs = c("I", "II", "III", "IV"),
  xlab = "Time (months)",
  ylab = "Survival probability",
  xlim = c(0, 72),
  break.time.by = 12,
  ggtheme = theme_classic(base_size = 9)
)

km_overall$plot <- km_overall$plot +
  ggtitle("A  Overall survival") +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")

km_nyha$plot <- km_nyha$plot +
  ggtitle("Survival by NYHA class") +
  theme(plot.title = element_text(face = "bold"), legend.position = "right")

km_overall$table <- km_overall$table +
  theme_classic(base_size = 7) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

km_nyha$table <- km_nyha$table +
  theme_classic(base_size = 7) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

left_km <- (km_overall$plot / km_overall$table) +
  plot_layout(heights = c(3, 0.95))

right_km <- (km_nyha$plot / km_nyha$table) +
  plot_layout(heights = c(3, 0.95))

p_km <- (left_km | right_km) +
  plot_layout(widths = c(0.9, 1.25))

cox_clinical <- coxph(
  Surv(time, status) ~ nyha_nominal +
    age + gender + bmi + smoking + diabetes +
    ejection_fraction + ns(log(creatinine), df = 2) + sodium,
  data = dat
)

ph_clinical <- cox.zph(cox_clinical, terms = FALSE)

ph_dat <- data.frame(
  time = rep(ph_clinical$time, 3),
  beta = c(
    ph_clinical$y[, "nyha_nominalII"],
    ph_clinical$y[, "nyha_nominalIII"],
    ph_clinical$y[, "nyha_nominalIV"]
  ),
  contrast = rep(
    c("II vs I", "III vs I", "IV vs I"),
    each = length(ph_clinical$time)
  )
)

p_ph <- ggplot(ph_dat, aes(x = time, y = beta)) +
  geom_point(size = 0.5, shape = 1, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.55, colour = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey45") +
  facet_wrap(~ contrast, nrow = 1) +
  scale_x_log10(
    breaks = c(1, 3, 6, 12, 24, 48),
    labels = c("1", "3", "6", "12", "24", "48")
  ) +
  labs(
    title = "B  Schoenfeld residual diagnostics",
    x = "Time (months)",
    y = "Beta(t)"
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", size = 9),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 8)
  )

rp_tvc1 <- stpm2(
  Surv(time, status) ~ nyha_II + nyha_III + nyha_IV +
    age + gender + bmi + smoking + diabetes +
    ejection_fraction + ns(log(creatinine), df = 2) + sodium,
  data = dat,
  df = 4,
  tvc = list(
    nyha_II = 1,
    nyha_III = 1,
    nyha_IV = 1
  )
)

time_grid <- seq(1, 60, length.out = 160)

pred_ref <- data.frame(
  time = time_grid,
  nyha_II = 0,
  nyha_III = 0,
  nyha_IV = 0,
  age = median(dat$age),
  gender = factor("Female", levels = levels(dat$gender)),
  bmi = median(dat$bmi),
  smoking = factor("Never", levels = levels(dat$smoking)),
  diabetes = factor("No", levels = levels(dat$diabetes)),
  ejection_fraction = median(dat$ejection_fraction),
  creatinine = median(dat$creatinine),
  sodium = median(dat$sodium)
)

hr_II <- predict(
  rp_tvc1,
  newdata = pred_ref,
  type = "hr",
  exposed = function(x) {
    x$nyha_II <- 1
    x
  },
  se.fit = TRUE,
  full = TRUE
)
hr_II$contrast <- "II vs I"

hr_III <- predict(
  rp_tvc1,
  newdata = pred_ref,
  type = "hr",
  exposed = function(x) {
    x$nyha_III <- 1
    x
  },
  se.fit = TRUE,
  full = TRUE
)
hr_III$contrast <- "III vs I"

hr_IV <- predict(
  rp_tvc1,
  newdata = pred_ref,
  type = "hr",
  exposed = function(x) {
    x$nyha_IV <- 1
    x
  },
  se.fit = TRUE,
  full = TRUE
)
hr_IV$contrast <- "IV vs I"

hr_dat <- rbind(hr_II, hr_III, hr_IV)

p_rp <- ggplot(hr_dat, aes(x = time, y = Estimate, colour = contrast, fill = contrast)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey45") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.13, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = cols_contrast) +
  scale_fill_manual(values = cols_contrast) +
  scale_y_log10(
    breaks = c(1, 2, 4, 8, 16, 32),
    labels = c("1", "2", "4", "8", "16", "32")
  ) +
  labs(
    title = "C  Flexible parametric time-varying NYHA effect",
    x = "Time (months)",
    y = "Adjusted hazard ratio",
    colour = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", size = 9),
    legend.position = "bottom"
  )

p_final <- p_km / (p_ph | p_rp) +
  plot_layout(heights = c(1.2, 0.8))

ggsave(
  "figures/figure_survival_modeling.pdf",
  p_final,
  width = 10.5,
  height = 7.8
)

ggsave(
  "figures/figure_survival_modeling.png",
  p_final,
  width = 10.5,
  height = 7.8,
  dpi = 300
)
