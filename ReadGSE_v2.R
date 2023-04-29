
# The code is developed for Harpreet's Analysis. 
# The number of cell lines is given below.
# ID      SF268   SF295   SF539   SNB19   SNB75   U251


library(readxl)
library(GEOquery)
library(preprocessCore)
library(stats)
library(samr)


set.seed(84048) 
# Download GSE1000 data set
gse <- getGEO("GSE7505")
#or
#load("~/Desktop/Uma/3. Task by Hp 10.04.2023/RCode/Gse.RData")

expression_matrix <- exprs(gse[[1]])
dim(expression_matrix)

Cell_names <- c("786O","A498","A549","ACHN","BT549","CAKI1","CCRFCEM","Colo205","Du145","EKVX","HCC2998",
                "HCT116","HCT15","HL60","HOP62","HOP92","HS578T","HT29","IGROV1","K562","KM12","LoxIVMI",
                "M14","Malme3M","MCF7","MDAMB231","MDAMB435","MDAN","ML1","Molt4","NCIADR","NCIH226","NCIH23",
                "NCIH322M","NCIH460","NCIH522","NH32","OVCAR3","OVCAR4","OVCAR5","OVCAR8","PC3","RPMI8226","RXR393",
                "SF268","SF295","SF593","SKMel2","SKMel28","SKMel5","SKOV3","SN12C","SNB19","SNB75","SR",
                "SW620","T47D","TK10","TK6","U251","UACC257","UACC62","UO31")
colnames(expression_matrix) <- Cell_names
#There are few cell lines whose sf2 value is not present. Discarding them from the expression matrix.
not_present <- c("ACHN", "CAKI1", "RXR393", "SN12C", "TK10", "UO31")
# Remove the specified columns using the subsetting operator
expression_matrix1 <- expression_matrix[, -which(colnames(expression_matrix) %in% not_present)]
dim(expression_matrix1)

# Normalize the data using quantile normalization and log2 transformation
norm_data <- normalize.quantiles(expression_matrix1)
colnames(norm_data) <- colnames(expression_matrix1)

# GBM cell line 
List <- c("SF268",   "SF295",   "SF593",   "SNB19",   "SNB75",   "U251")
New_expression <- norm_data[ ,List]

#Read the data from excel.

Table_From_Ref1 <- read_excel("~/Desktop/Uma/3. Task by Hp 10.04.2023/GSE7505/Table From Ref2.xlsx")

# Now according to the column name of expression_matrix1 get the SF2 value from Table_From_Ref1
lookup_data <- data.frame(id = c(colnames(New_expression)))
merge_data <- data.frame(id = c(Table_From_Ref1$`Cell line`), SF2 = c(Table_From_Ref1$SF2))
merged_data <- merge(lookup_data, merge_data, by = "id")
X <- New_expression    #Or Quantile normalization is working fine.
Y <- merged_data$SF2
a<-SAM(X,Y,resp.type="Quantitative",nperms=50,fdr.output=.5)

print(a) #plot results 
plot(a)

Signature <- a$siggenes.table
Signature$genes.up
Signature$genes.lo



# Create a list to store the genes.up values for each iteration
genes_up_list <- vector(mode = "list", length = ncol(X))
genes_lo_list <- vector(mode = "list", length = ncol(X))

# Perform leave-one-out cross-validation
for (i in 1:ncol(X)) {
  # Remove the i-th observation from the dataset
  X_temp <- X[,-i]
  Y_temp <- Y[-i]
  
  # Run the SAM function on the modified dataset
  a <- SAM(X_temp, Y_temp, resp.type="Quantitative", nperms=50, fdr.output=.5)
  
  # Store the genes.up values for the iteration
  genes_up_list[[i]] <- a$siggenes.table$genes.up
  genes_lo_list[[i]] <- a$siggenes.table$genes.lo
}



# Now find the common gene names maximally found among the all iteration.
# Combine the gene lists into one vector
all_genes_up <- c()
for (i in 1:length(genes_up_list)) {
  all_genes_up <- c(all_genes_up, genes_up_list[[i]][,1])
}

all_genes_lo <- c()
for (i in 1:length(genes_up_list)) {
  all_genes_lo <- c(all_genes_lo, genes_lo_list[[i]][,1])
}

# Count the number of occurrences of each gene
gene_counts_up <- table(all_genes_up)
gene_counts_lo <- table(all_genes_lo)

# Sort the genes by their counts in descending order
sorted_genes_up <- sort(gene_counts_up, decreasing = TRUE)
sorted_genes_lo <- sort(gene_counts_lo, decreasing = TRUE)

# Print the top 10 most common genes
head(sorted_genes_up, n = 20)
head(sorted_genes_lo, n = 20)

#Name of top 10 genes.
Top10 <- names(sorted_genes_up[1:10])
Top10 <- names(sorted_genes_lo[1:10])

# Suppose we are using leave one out cross validation.
# Each time one sample will be removed and it will be treated as test sample.
# Let us the ithe sample is removed.

i = 2
X_temp <- X[,-i]
Y_temp <- Y[-i]

#Name of top 10 genes.
Top10 <- names(sorted_genes_up[1:10])
#Top10 <- names(sorted_genes_lo[1:10])
names_no_g <- gsub("g", "", Top10)
names_no_quotes <- as.numeric(names_no_g)


# Now the data matrix will be 
X_temp1 <-  X_temp[c(names_no_quotes), ]

# create a data frame containing the gene expression data and response variable
gene_data <- data.frame(gene1 = X_temp1[1,], 
                        gene2 = X_temp1[2,], 
                        gene3 = X_temp1[3,], 
                        gene4 = X_temp1[4,], 
                        gene5 = X_temp1[5,], 
                        gene6 = X_temp1[6,], 
                        gene7 = X_temp1[7,], 
                        gene8 = X_temp1[8,], 
                        gene9 = X_temp1[9,], 
                        gene10 = X_temp1[10,], 
                        response = Y_temp)

# fit a multiple linear regression model
model <- lm(response ~ gene1 + gene2 + gene3 + gene4 + gene5 + gene6 + gene7 + gene8 + gene9 + gene10, data = gene_data)
model <- lm(response ~ gene1 + gene2 + gene3 + gene4, data = gene_data)
model <- lm(response ~ gene1 + gene2 + gene8 , data = gene_data)

# print the model summary
summary(model)

# Now predict the test sample.
Test_Sample <- X[c(names_no_quotes),i]

# create a data frame containing new gene expression data for which to make a single prediction+
new_gene_data <- data.frame(gene1 = Test_Sample[1], 
                            gene2 = Test_Sample[2], 
                            gene3 = Test_Sample[3], 
                            gene4 = Test_Sample[4], 
                            gene5 = Test_Sample[5], 
                            gene6 = Test_Sample[6], 
                            gene7 = Test_Sample[7], 
                            gene8 = Test_Sample[8], 
                            gene9 = Test_Sample[9], 
                            gene10 = Test_Sample[10])

# make a single prediction using the model
predicted_response <- predict(model, newdata = new_gene_data)

# print the predicted response value
print(predicted_response)
#Original redio response value of the test sample. 
Y[i]



#Check for high correlation between the predictor variables.

#Multicollinearity occurs when the predictor variables are highly correlated with 
#each other, leading to unstable and imprecise coefficient estimates.


# How to find the best combination of genes for multiple linear regression model?
#   
# To find the best combination of genes for a multiple linear regression model, you can use 
# a technique called feature selection. Feature selection is the process of identifying a 
# subset of the most relevant features (i.e., genes) that contribute the most to the 
# prediction accuracy of the model while minimizing overfitting.
# 
# Here are some common approaches to feature selection:
#   
# Forward selection: Start with an empty model and iteratively add one gene at a time, 
# choosing the gene that provides the best improvement in the model's prediction accuracy. 
# This process continues until no further improvement is obtained.
# 
# Backward elimination: Start with a full model including all genes and iteratively remove 
# one gene at a time, choosing the gene whose removal has the least impact on the model's 
# prediction accuracy. This process continues until removing any additional gene results 
# in a significant drop in the model's prediction accuracy.
# 
# Recursive feature elimination: This is a more sophisticated approach that involves training 
# the model on all genes, ranking the genes by their importance, removing the least 
# important gene, and repeating the process until the desired number of genes is reached.
# 
# Regularization: This approach adds a penalty term to the model's objective function, 
# which encourages the model to select a sparse set of genes by setting the coefficients 
# of irrelevant genes to zero. Lasso and Ridge regression are examples of regularization techniques.
# 
# These techniques are not mutually exclusive, and you can use them in combination to 
# obtain the best subset of genes for your multiple linear regression model. However, 
# keep in mind that the choice of the feature selection technique can depend on the size 
# of the dataset, the number of genes, and the complexity of the biological system being studied.



#Forward selection using R

library(leaps)
df <- data.frame(y = c(1, 2, 3, 4, 5), 
                 x1 = c(2, 3, 4, 5, 6), 
                 x2 = c(3, 4, 5, 6, 7),
                 x3 = c(4, 5, 6, 7, 8),
                 x4 = c(5, 6, 7, 8, 9),
                 x5 = c(6, 7, 8, 9, 10))
fit <- regsubsets(y ~ ., data = df, method = "forward")
summary(fit)

plot(fit, scale = "r2")


df <- data.frame(gene1 = X_temp1[1,], 
                        gene2 = X_temp1[2,], 
                        gene3 = X_temp1[3,], 
                        gene4 = X_temp1[4,], 
                        gene5 = X_temp1[5,], 
                        gene6 = X_temp1[6,], 
                        gene7 = X_temp1[7,], 
                        gene8 = X_temp1[8,], 
                        gene9 = X_temp1[9,], 
                        gene10 = X_temp1[10,], 
                        y = Y_temp)


fit <- regsubsets(y ~ ., data = df, method = "forward")
summary(fit)

plot(fit, scale = "r2")


model <- lm(response ~ gene1 + gene2 + gene3 , data = gene_data)
model <- lm(response ~ gene1 + gene2 + gene8 , data = gene_data)

summary(model)


# Now predict the test sample.
Test_Sample <- X[c(names_no_quotes),i]

# create a data frame containing new gene expression data for which to make a single prediction
new_gene_data <- data.frame(gene1 = Test_Sample[1], 
                            gene2 = Test_Sample[2], 
                            gene3 = Test_Sample[3] 
                            )

# make a single prediction using the model
predicted_response <- predict(model, newdata = new_gene_data)

# print the predicted response value
print(predicted_response)
#Original redio response value of the test sample. 
Y[i]


################################################################################
#How to do the Backward selection using R 


#model <- lm(mpg ~ ., data = mtcars)
# Previously the mdel is created.

# Perform backward selection using the step function
selected_model <- step(model, direction = "backward")
