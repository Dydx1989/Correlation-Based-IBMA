# Load required libraries
library(survival)
library(plotly)
library(heatmaply)
library(BMA)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("iterativeBMAsurv")
library(iterativeBMAsurv)

# Load data
Datachopp = read.csv("chopp.csv")
dim(Datachopp)
head(Datachopp[1:4])

# Generate Heatmap
dat_X = Datachopp[, -c(1, 2)]
head(dat_X[1:4])
heatmaply(dat_X[1:30], k_row = 2, k_col = 2)
heatmaply(cor(dat_X[1:30]), k_row = 2, k_col = 2, colors = c("green", "blue"))

# Cumulative Hazard function
calcna = function(time, status) {
  na.fit = survfit(coxph(Surv(time, status) ~ 1), type = "aalen")
  jumps = c(0, na.fit$time, max(time))
  surv = c(1, na.fit$surv, na.fit$surv[length(na.fit$surv)])
  neglogsurv = -log(surv)
  naest = numeric(length(time))
  for (i in 2:length(jumps)) {
    naest[which(time >= jumps[i - 1] & time <= jumps[i])] = neglogsurv[i - 1]
  }
  return(naest) 
}

# Calculate Hazard function
Hazard = calcna(Datachopp$survtime, Datachopp$status)
y = Hazard
Datachopp2 = data.frame(y, dat_X)
head(Datachopp2[1:4])

# Step One Filtering
chopp.selection <- apply(Datachopp2[-1], 2, function(x) cor.test(x, Datachopp2[, 1])$p.value)
chopp.selection2 <- sort(chopp.selection[chopp.selection <= 0.001])
length(chopp.selection2)
print(length(chopp.selection2))

# Generate correlation details
selected_genes <- names(chopp.selection2)
results <- data.frame()
for (gene in selected_genes) {
  test <- cor.test(Datachopp2[, gene], Datachopp2[, 1])
  results <- rbind(results, data.frame(
    Gene = gene,
    Correlation = test$estimate,
    P_Value = test$p.value,
    Statistic = test$statistic,
    Conf_Lower = test$conf.int[1],
    Conf_Upper = test$conf.int[2],
    Method = test$method
  ))
}
print(results)

# Correlation Plot
library(ggplot2)
Genes_selected <- c(0.001, 0.01, 0.02, 0.03, 0.04, 0.1, 0.2, 0.3, 0.4, 0.5)
alpha_values <- c(27, 89, 158, 208, 247, 522, 907, 1287, 1647, 2032)
data <- data.frame(Alpha = alpha_values, Genes_selected)

ggplot(data, aes(y = Genes_selected, x = Alpha)) +
  geom_point(size = 3, color = "red", shape = 19, alpha = 0.8) +
  geom_line(color = "blue", size = 1) +
  geom_smooth(method = "lm", se = FALSE, color = "green") +
  xlab("Number of Genes Selected") +
  ylab("Alpha Values") + theme_minimal(15)

# Selected Genes with Hazard
Chopp_subset2 <- subset(Datachopp, select = c("survtime", "status", names(chopp.selection2)))
print(length(Chopp_subset2))
head(Chopp_subset2[1:4])
write.csv(Chopp_subset2, "Chopp_subset2.csv")
heatmaply(Chopp_subset2[-c(1, 2)], k_row = 2, k_col = 2)
heatmaply(cor(Chopp_subset2[-c(1, 2)]), k_row = 2, k_col = 2, colors = c("green", "blue"))

# Gene selection via iterativeBMAsurv
ret.list <- iterateBMAsurv.train.wrapper(
  x = Chopp_subset2[-c(1, 2)], 
  surv.time = Chopp_subset2[, 1], 
  cens.vec = Chopp_subset2[, 2], 
  nbest = 10, maxNvar = 15
)
print(ret.list)
ret.bic.surv <- ret.list$obj
gene.names <- ret.list$curr.names
top.gene.names <- gene.names[ret.bic.surv$probne0 > 0.0]
length(top.gene.names)
top.gene.names

# Variable Importance Plot
Datachop1 = read.csv("Chopp_subset2.csv")
Datachop3 = Datachop1[-c(1:3)]
Datachop_M <- apply(Datachop3, MARGIN = 2, FUN = median, na.rm = TRUE)
Datachop_Or <- order(Datachop_M, decreasing = FALSE)
library(RColorBrewer)
n <- ncol(Datachop3)
qual_col_palsD1 <- brewer.pal.info[brewer.pal.info$category == "qual", ]
col_vectorD1 <- unlist(mapply(brewer.pal, qual_col_palsD1$maxcolors, rownames(qual_col_palsD1)))
boxplot(Datachop3[, Datachop_Or], 
        col = sample(col_vectorD1, n), 
        las = 2, 
        cex.axis = 0.7, 
        ylab = "Ranger Normalized Permutation Importance", 
        boxwex = 0.6, 
        main = "DLBCL: Variable Importance")

# Cox Proportional Hazard Model
fitcox <- coxph(Surv(Hazard) ~ X1558999_x_at + X229839_at + X1553317_s_at +
                  X240777_at + X237797_at + X1569344_a_at + X244434_at + X242758_x_at + X243713_at +
                  X1557366_at + X237515_at + X205908_s_at +
                  X237493_at + X244346_at + X1563643_at, data = Datachopp, x = TRUE)
summary(fitcox)
AIC(fitcox)
BIC(fitcox)

# Other Survival Models
fitweibull <- survreg(Surv(Hazard) ~ X1558999_x_at + X229839_at + X1553317_s_at +
                        X240777_at + X237797_at + X1569344_a_at + X244434_at + X242758_x_at + X243713_at +
                        X1557366_at + X237515_at + X205908_s_at +
                        X237493_at + X244346_at + X1563643_at, data = Datachopp, dist = "weibull")
summary(fitweibull)
AIC(fitweibull)
BIC(fitweibull)

fitexponential <- survreg(Surv(Hazard) ~ X1558999_x_at + X229839_at + X1553317_s_at +
                            X240777_at + X237797_at + X1569344_a_at + X244434_at + X242758_x_at + X243713_at +
                            X1557366_at + X237515_at + X205908_s_at +
                            X237493_at + X244346_at + X1563643_at, data = Datachopp, dist = "exponential")
AIC(fitexponential)
BIC(fitexponential)

fitlognormal <- survreg(Surv(Hazard) ~ X1558999_x_at + X229839_at + X1553317_s_at +
                          X240777_at + X237797_at + X1569344_a_at + X244434_at + X242758_x_at + X243713_at +
                          X1557366_at + X237515_at + X205908_s_at +
                          X237493_at + X244346_at + X1563643_at, data = Datachopp, dist = "lognormal")
AIC(fitlognormal)
BIC(fitlognormal)

fitloglogistic <- survreg(Surv(Hazard) ~ X1558999_x_at + X229839_at + X1553317_s_at +
                            X240777_at + X237797_at + X1569344_a_at + X244434_at + X242758_x_at + X243713_at +
                            X1557366_at + X237515_at + X205908_s_at +
                            X237493_at + X244346_at + X1563643_at, data = Datachopp, dist = "loglogistic")
AIC(fitloglogistic)
BIC(fitloglogistic)
