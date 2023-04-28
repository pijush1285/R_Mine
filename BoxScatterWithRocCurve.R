###############################################################################
# This is a Box Scatter plot and ROC attached aside.
# The input data set is given below:

# Mutationtype	has
# Missense	    1.462177
# Missense	    3.068568
# Truncation	  11.492844
# Truncation	  129.645902

# Developed By Dr. Pijush Das

###############################################################################

library(readxl)
library(ROCR)

daataa27 <- read_excel("daataa27.xlsx")
View(daataa27)

#Adding the label.
daataa27$Label <- ifelse(daataa27$Mutationtype == "Missense", 1, 0)

#Converting the FPKM value into Log2 value.
Missense_values <- log2(daataa27$has[daataa27$Mutationtype == "Missense"])
Truncation_values <- log2(daataa27$has[daataa27$Mutationtype == "Truncation"])

# Calculating the p-value
Missense_values1 <- daataa27$has[daataa27$Mutationtype == "Missense"]
Truncation_values1 <- daataa27$has[daataa27$Mutationtype == "Truncation"]
Ans <- t.test(Missense_values1,Truncation_values1)
pval <- Ans$p.value


#Calculating the ROC curve value.
#Create a binary classification data set
labels <- factor(daataa27$Label)
scores <- c(daataa27$has)
# Calculate the prediction performance
pred <- prediction(scores, labels)
perf <- performance(pred, "tpr", "fpr")
#Calculate the AUC value
auc_value <- auc(labels, scores)

#Plotting Box Scatter plot
par(mfrow = c(1, 2))
mytitle<-'hsa-mir-196b'
boxplot(Missense_values,Truncation_values,col = c("lavender","lightpink1"),main = mytitle, xlab = paste0("p-value = ", round(pval, 4)),ylab = "Log2(FPKM)", names = c("Missense","Truncation"))
stripchart(Missense_values, vertical = TRUE,method = "jitter", add = TRUE, pch = 16, col = c('blue'))
stripchart(Truncation_values, vertical = TRUE, at = 2,method = "jitter", add = TRUE, pch = 16, col = c('red'))
#text(x = 1, y = 0, labels = paste0("p = ", round(pval, 4)), pos = 4)


#Plot the ROC curve
plot(perf@x.values[[1]], 1 - perf@y.values[[1]], main = "ROC Curve", 
     xlab = "1 - Specificity", ylab = "Sensitivity", type = "l", 
     col = "blue", lwd = 2)
text(x = 0.6, y = 0.9, labels = paste0("AUC = ", round(auc_value, 2)), pos = 4)
#Add a cutoff line
abline(0, 1, lty = 2, col = "red")

