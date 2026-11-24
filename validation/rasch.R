# Run with Rscript --vanilla validation/rasch.R after installing glmmplanr.
library(glmmplanr)
library(lme4)

design <- design_crossed(list(), allocation = "crossed")
person_sd <- 0.8
item_sd <- 0.6
dgp <- calibrate_dgp(design, cell_probs = 0.65, chance = 0,
                     subj_sd = person_sd, item_sd = item_sd)
beta <- unname(dgp$beta["(Intercept)"])
control <- glmerControl(optimizer = "bobyqa")

rasch_replication <- function(schedule, seed) {
  set.seed(seed)
  data <- simulate_dgp_data(dgp, schedule)
  # Replay the independent draws on the IRT scale. The origin is E(theta) = 0.
  set.seed(seed)
  ability <- rnorm(schedule$n_subj, sd = person_sd)
  difficulty <- -beta - rnorm(schedule$n_item, sd = item_sd)
  probability <- plogis(ability[data$subj_id] - difficulty[data$item_id])
  stopifnot(identical(data$y, rbinom(nrow(data), 1L, probability)))

  fit <- glmer(y ~ 1 + (1 | subj) + (1 | item), data, binomial(), control = control)
  conv <- fit@optinfo$conv
  stopifnot(conv$opt == 0, is.null(conv$lme4$code))
  estimates <- ranef(fit)
  variances <- VarCorr(fit)
  theta_hat <- estimates$subj[["(Intercept)"]]
  b_hat <- -unname(fixef(fit)[1]) - estimates$item[["(Intercept)"]]
  c(person_rmse = sqrt(mean((theta_hat - ability)^2)),
    item_rmse = sqrt(mean((b_hat - difficulty)^2)),
    person_sd = sqrt(variances$subj[1, 1]), item_sd = sqrt(variances$item[1, 1]),
    singular = as.numeric(isSingular(fit)))
}

for (n_item in c(12L, 24L)) {
  schedule <- generate_design_schedule(design, n_subj = 1728L %/% n_item, n_item = n_item)
  seeds <- 461000L + 100L * n_item + seq_len(32L)
  runs <- vapply(seeds, function(seed) rasch_replication(schedule, seed), numeric(5))
  rmse <- runs[c("person_rmse", "item_rmse"), , drop = FALSE]
  print(data.frame(n_subj = schedule$n_subj, n_item = n_item, replications = ncol(runs),
                   singular_fits = sum(runs["singular", ]), parameter = rownames(rmse),
                   mean_rmse = rowMeans(rmse), mcse = apply(rmse, 1, sd) / sqrt(ncol(runs))),
        row.names = FALSE)
  sd_error <- runs[c("person_sd", "item_sd"), , drop = FALSE] - c(person_sd, item_sd)
  print(data.frame(parameter = rownames(sd_error), bias = rowMeans(sd_error),
                   rmse = sqrt(rowMeans(sd_error^2))), row.names = FALSE)
}

# Independently extend the schedule to four random intercepts. Class IDs are
# local to a school. All classrooms use the same bank of crossed items.
data <- generate_design_schedule(design, n_subj = 288L, n_item = 16L)$trials
data$school <- factor((data$subj_id - 1L) %/% 24L + 1L)
class_number <- (data$subj_id - 1L) %/% 8L %% 3L + 1L
data$class_id <- interaction(data$school, class_number, drop = TRUE)
sds <- c(school = 0.4, class_id = 0.5, subj = 0.8, item = 0.6)
stopifnot(nlevels(data$school) == 12L, nlevels(data$class_id) == 36L)
set.seed(472031L)
# The conditional intercept is fixed here. Extra variance changes the marginal mean.
eta <- rep(beta, nrow(data))
for (group in names(sds)) {
  id <- data[[group]]
  eta <- eta + rnorm(nlevels(id), sd = sds[[group]])[id]
}
data$y <- rbinom(nrow(data), 1L, plogis(eta))
formula <- reformulate(sprintf("(1 | %s)", names(sds)), response = "y")
fit <- glmer(formula, data, binomial(), control = control)
stopifnot(fit@optinfo$conv$opt == 0, is.null(fit@optinfo$conv$lme4$code),
          setequal(names(VarCorr(fit)), names(sds)))
print(data.frame(group = names(sds), groups = unname(ngrps(fit)[names(sds)]),
                 input_sd = unname(sds),
                 fitted_sd = vapply(VarCorr(fit)[names(sds)], function(x) sqrt(x[1, 1]), numeric(1))),
      row.names = FALSE)
