# Engine functions for feature association analysis.

# Fit a Cox proportional hazards regression model for each protein.
exprCox <- function(expr, time, status) {
  pro.name <- rownames(expr)

  univ_formulas <- as.formula(paste('Surv(time, status)~', "expr"))
  univ_models <- lapply(pro.name, function(x) {
    clinic_clus <- data.frame(time = time, status = status,
                              expr = unlist(expr[x, ]))
    coxph(univ_formulas, data = clinic_clus)
  })

  univ_results <- lapply(univ_models, function(x){
    x <- summary(x)
    pval <- x$wald["pvalue"]
    wald.test <- x$wald["test"]
    beta <- x$coef[1]   #coefficient beta
    HR <- x$coef[2]   #exp(beta)
    HR.confint.lower <- x$conf.int[, "lower .95"]
    HR.confint.upper <- x$conf.int[, "upper .95"]

    res <- c(beta, HR, HR.confint.lower, HR.confint.upper, wald.test, pval)
    names(res) <- c("beta", "HR", "HR.lower", "HR.upper",
                    "wald.test", "pval")
    res
  })
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  res <- as.data.frame(res)
  res$padj <- p.adjust(res$pval, method = "BH")

  rownames(res) <- pro.name

  res
}

# Conduct Welch’s one-way test or ANOVA for categorical feature.
exprCatFeature <- function(expr, feature, varPCut = 0.05) {
  pro.name <- rownames(expr)

  library(snowfall)
  sfInit(parallel = TRUE, cpus = 4)
  sfExport("expr", "pro.name",
           "feature", "varPCut")

  mean_res <- sfLapply(pro.name, function(x) {
    if(sum(table(feature[which(!is.na(unlist(expr[x, ])))]) == 0) == 0 &
       sum(table(feature[which(!is.na(unlist(expr[x, ])))]) < 2) == 0 &
       length(table(feature[which(!is.na(unlist(expr[x, ])))])) > 1) {
      varTestP <- bartlett.test(unlist(expr[x, ]) ~ feature)$p.value

      if (varTestP >= varPCut) {
        varEqual <- T
        anova_res <- summary(aov(unlist(expr[x, ]) ~ feature))
        pval <- anova_res[[1]][["Pr(>F)"]][1]
      } else {
        varEqual <- F
        pval <- oneway.test(unlist(expr[x, ]) ~ feature, var.equal = FALSE)$p.val
      }
    } else {
      varEqual <- NA
      pval <- NA
    }

    res <- c(varEqual, pval)
    names(res) <- c("varEqual", "pval")
    res
  })
  sfStop()

  res <- t(as.data.frame(mean_res, check.names = FALSE))
  res <- as.data.frame(res)
  res$padj <- p.adjust(res$pval, method = "BH")

  rownames(res) <- pro.name

  res
}
