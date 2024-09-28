library(data.table)
library(rcompanion)

# Inputs
trait <- snakemake@input[["trait"]]
covariates <- snakemake@input[["covariates"]]
info <- snakemake@input[["info"]]
prs_score <- snakemake@input[["prs_score"]]

# Outputs
coefficients <- snakemake@output[["coefficients"]]
r2 <- snakemake@output[["r2"]]

# Wildcard
wildcard <- snakemake@wildcards[["BaseTarget"]]
print(wildcard)

# Load inputs
df_trait <- fread(trait)
df_cov <- fread(covariates)
df_info <- fread(info)
df_score <- fread(prs_score)


# Column names for df_trait and df_cov
colnames(df_trait) <- c("FID", "IID", "Trait")
colnames(df_cov)[c(1, 2)] <- c("FID", "IID")

# Merge phenotype and covariate matrix
df <- merge(df_trait, df_cov, by= c("FID", "IID"))


# Define if Target Trait is continuous or binary
variable <- as.character(df_info[df_info$BaseTarget == wildcard, "BinaryTargetOption"]) 
print(variable)

if (variable == "F") {
  print("Continuous trait: linear regression starting ...")
  # Calculate null model (without PRS)
  null.model <- lm(Trait~., data=df[,!c("FID","IID")], na.action = na.omit)
  null.r2 <- summary(null.model)$r.squared # R2 of the null model
  # Model with PRS: 
  # Merge with phenotype and covariates matrix
  prs <- merge(df, df_score[,c("FID","IID", "SCORESUM")], by=c("FID", "IID"))
  model <- lm(Trait~., data=prs[,!c("FID","IID")]) 
  model.r2 <- summary(model)$r.squared # model R2
  # Option 1: R2 of PRS is calculated as the model R2 minus the null R2
  prs.r2 <- model.r2-null.r2
  # Option 2: R2 is calculated in the same way as PRSice2:
  prsice.r2 <- 1 - ((1- model.r2) / (1 - null.r2))
  # Option 3: R2 partial 
  SSE.null<-sum(null.model$residuals**2)  
  SSE.full<-sum(model$residuals**2)  
  R2.partial<-(SSE.null-SSE.full)/SSE.null 
  # Coeficients and p-values of association 
  prs.coef <- data.frame(summary(model)$coeff)
  colnames(prs.coef) <- c("Estimate", "Std.Error", "t.value", "p.value")
  prs.coef$P.value.0.05 <- ifelse(prs.coef$p.value <0.05, "*", "")
  # Store r2 results
  prs.result <- data.frame("Null_Model_R2" = null.r2, "PRS_Model_R2" = model.r2, 
                           "PRS_R2" = prs.r2, "PRSice2_like_r2" = prsice.r2, "R2_partial" = R2.partial)
} else {
  print("Binary trait: glm starting ...")
  # Calculate null model (without PRS)
  null.model <- glm(Trait~., data=df[,!c("FID","IID")], family = binomial(link='logit'), na.action = na.omit)
  # Model with PRS: 
  # Merge with phenotype and covariates matrix
  prs <- merge(df, df_score[,c("FID","IID", "SCORESUM")], by=c("FID", "IID"))
  model <- glm(Trait~., data=prs[,!c("FID","IID")], family = binomial(link='logit'))
  # Pseudo R2 Nagelkerke
  Adj.R2<-nagelkerke(model, null=null.model)$Pseudo.R.squared.for.model.vs.null[3]
  # Store r2 result
  prs.result <- data.frame("Nagelkerke_Pseudo_R2" = Adj.R2)
  # Coeficients and p-values of association
  or <- exp(cbind(OR = coef(model),confint(model)))
  coef <- summary(model)$coefficients
  colnames(coef) <- c("Estimate", "Std.Error", "z.value", "p.value")
  prs.coef <- as.data.frame(cbind(coef, or))
  prs.coef$P.value.0.05 <- ifelse(prs.coef$p.value<0.05, "*", "")
}


# Save prs.results and prs.coef
fwrite(prs.coef, coefficients, quote = FALSE, col.names = TRUE, row.names = TRUE, sep = "\t")
fwrite(prs.result, r2, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")



