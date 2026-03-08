library(survival)
library(survminer)

df <- read.csv("surv_data_wrisk_clinicalvars.csv")

#calculate HR for risk score (binary groups)
df[["COHORT"]] <- factor(df[["COHORT"]], levels = c("LowRisk", "HighRisk"))

cox_formula <- as.formula(paste0("Surv(", "Duration.Months", ", ", "Event", ") ~ ", "COHORT"))
cox_model   <- coxph(cox_formula, data = df)

print(summary(cox_model))
hr_table <- as.data.frame(
  cbind(
    HR    = exp(coef(cox_model)),
    CI_lo = exp(confint(cox_model)[, 1]),
    CI_hi = exp(confint(cox_model)[, 2]),
    p_value = summary(cox_model)$coefficients[, "Pr(>|z|)"]
    )
  )
print(round(hr_table, 4))

#Calculate HR for grade
df[["Grade"]] <- factor(df[["Grade"]], levels = c(2, 3))

cox_formula <- as.formula(paste0("Surv(", "Duration.Months", ", ", "Event", ") ~ ", "Grade"))
cox_model   <- coxph(cox_formula, data = df)

print(summary(cox_model))
hr_table <- as.data.frame(
  cbind(
    HR    = exp(coef(cox_model)),
    CI_lo = exp(confint(cox_model)[, 1]),
    CI_hi = exp(confint(cox_model)[, 2]),
    p_value = summary(cox_model)$coefficients[, "Pr(>|z|)"]
  )
)
print(round(hr_table, 4))


#Calculate HR for subtypes (astrocytoma vs oligodendroglioma)
df[["Subtype"]] <- factor(df[["Subtype"]], levels = c("ODG", "AST"))

cox_formula <- as.formula(paste0("Surv(", "Duration.Months", ", ", "Event", ") ~ ", "Histotype"))
cox_model   <- coxph(cox_formula, data = df)

print(summary(cox_model))
hr_table <- as.data.frame(
  cbind(
    HR    = exp(coef(cox_model)),
    CI_lo = exp(confint(cox_model)[, 1]),
    CI_hi = exp(confint(cox_model)[, 2]),
    p_value = summary(cox_model)$coefficients[, "Pr(>|z|)"]
  )
)
print(round(hr_table, 4))
