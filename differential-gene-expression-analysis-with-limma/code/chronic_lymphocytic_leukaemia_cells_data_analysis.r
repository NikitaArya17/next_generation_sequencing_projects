BiocManager::install(c("limma", "GO.db", "org.Hs.eg.db"))

library(Biobase)
library(limma)
library(GO.db)
library(org.Hs.eg.db)

cll_data <- readRDS('cll-eset.rds')

# Extract the phenotypic data, feature data and expression data

p_cll <- pData(cll_data)
f_cll <- fData(cll_data)
x_cll <- exprs(cll_data)

# We begin by visualising the expression levels of the first gene for 
#patients with stable and progressive leukaemia.
#This will be done with boxplots.

# Create a boxplot of the first gene in the expression matrix
boxplot(x_cll[1,] ~ p_cll[,"Disease"], main = f_cll[1,"symbol"]) 
#The `main` argument is to specify the title.

#The 3 different datasets containing the expression matrix, 
#feature data and phenotypic data can be combined into a single object 
#using the Bioconductor class ExpressionSet. 
#This allows for ease of subsetting.

# Create ExpressionSet object
eset_cll <- ExpressionSet(assayData = x_cll,
                      phenoData = AnnotatedDataFrame(p_cll),
                      featureData = AnnotatedDataFrame(f_cll))

# View the number of features (rows) and samples (columns)
dim(eset_cll)

# Subset to only include the first 10 samples (columns)
eset_cll_sub <- eset_cll[, 1:10]

# Check the dimensions of the subset
dim(eset_cll_sub)

# Create a boxplot of the 1000th gene in eset_sub
boxplot(exprs(eset_cll_sub)[1000, ] ~ pData(eset_cll_sub)[, "Disease"],
        main = fData(eset_cll_sub)[1000, "symbol"])

#We must now specify a linear model for R to use to compare the 2 groups 
#of differentially expressed leukaemia genes.
#We begin by constructing a design matrix with the intercept coefficient 
#and the coefficient indicating the disease status.

# Create design matrix for leukemia study
design_cll <- model.matrix(~Disease, data = pData(eset_cll))

# Count the number of samples modeled by each coefficient
colSums(design_cll)

#We can now test for differential expression between the two groups using limma


# Fit the model
fit_cll <- lmFit(eset_cll, design_cll)

# Calculate the t-statistics
fit_cll <- eBayes(fit_cll)

# Summarize results
results_cll <- decideTests(fit_cll[, "Diseasestable"])
summary(results_cll)

#Testing leukaemia data for differential expression 
#using the traditional treatment-contrasts parameterisation

# Create design matrix with no intercept
design_cll_no_intercept <- model.matrix(~0 + Disease, data = pData(eset_cll))

# Count the number of samples modeled by each coefficient
colSums(design_cll_no_intercept)

#This contrasts with the group means parameterisation 
#which we will now use.
#The coefficients in the linear model now represent the mean expression levels 
#for each of the two groups of samples.
#So we now need to specify a custom contrast 
#to test for differences between the two groups


# Create a contrasts matrix
cm_cll <- makeContrasts(status = Diseaseprogres. - Diseasestable,
                    levels = design_cll_no_intercept)

# View the contrasts matrix
cm_cll

#We have defined two required matrices, 
#the design matrix and the contrast matrix.
#We can now test for differential expression.

# Fit the model
fit_cll_no_intercept <- lmFit(eset_cll, design_cll_no_intercept)

# Fit the contrasts
fit2_cll <- contrasts.fit(fit_cll_no_intercept, contrasts = cm_cll)

# Calculate the t-statistics for the contrasts
fit2_cll <- eBayes(fit2_cll)

# Summarize results
results_cll_no_intercept <- decideTests(fit2_cll)
summary(results_cll_no_intercept)

#We have analysed and found results. Now we need to visualise them.
#The first method is to plot a histogram of p-values.

# Obtain the summary statistics for every gene
stats_cll <- topTable(fit2_cll, number = nrow(fit2_cll), sort.by = "none")

# Plot a histogram of the p-values
hist(stats_cll[, "P.Value"])

#A volcano plot visualises the extent of differential expression in the study. 
#It displays log odds on the y axis and log fold change on the x axis.

# Create a volcano plot. Highlight the top 5 genes
volcanoplot(fit2_cll, highlight = 5, names = fit2_cll$genes[, "symbol"])

#For further interpretation of differential expression results, 
#a common technique is to test for enrichment in known gene sets, 
#by exploring datasets of genes involved in the same biological pathway, 
#that were not included in the study. 
#The KEGG database is a good source of such datasets.

# Extract the entrez gene IDs
entrez_cll <- fit2_cll$genes[, "entrez"]

# Test for enriched KEGG Pathways
enrich_kegg_cll <- kegga(fit2_cll, geneid = entrez_cll, species = "Hs")

# View the top 20 enriched KEGG pathways
top_KEGG_pathways <- topKEGG(enrich_kegg_cll)
write.csv(top_KEGG_pathways, "top_KEGG_pathways.csv")

#We will now test for enrichment in gene sets 
#that are known to influence the same biological process or molecular function. 
#These are two of the three Gene Ontology categories, 
#the third being cellular component.

# Test for enriched GO categories
enrich_go <- goana(fit2_cll[1:500, ], geneid = entrez_cll[1:500], species = "Hs")

# View the top 20 enriched GO Biological Processes
top_GO_BP <- topGO(enrich_go, ontology = "BP")
write.csv(top_GO_BP, "top_GO_BP.csv")