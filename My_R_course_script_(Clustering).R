#---------------------------------------------------------------
# My_R_course_script_(Clustering).R

# The script is my practical part (Clustering) of our R course.

#---------------------------------------------------------------


rm(list = ls())  # cleaning all the variables in the Environment



### Libraries and methods
#install.packages("pvclust",dependencies = TRUE)
#install.packages("gplots",dependencies = TRUE)
library(pvclust)
library(gplots)


### Pre-processing
## Input
#folder_input <- '/home/vasilyromanov/Documents/R course/Clustering part/'   # for Linux
folder_input <- 'C:\\Users\\Vasily Romanov\\Desktop\\Work\\R course\\Clustering part\\'   # for Windows
matrix.org <- read.table(paste(folder_input, 'covert.withsymbols.log.Dec05.tab', sep = ""),
                         sep = "\t", header = TRUE, dec = ".", row.names = 1)


## Output folder
setwd(paste(folder_input, 'Plots', sep = ""))



## To exclude genes, which have a lower expression and which show only a slight variance of expression:
# SDs
variance <- apply(matrix.org, MARGIN = 1, FUN = sd)
variance.rank <- rank(1 / variance)
png(filename = "variance.png", width = 35, height = 20, units = "cm", res = 300)   # .png 1
plot(variance.rank, variance, xlab = "Genes", ylab = "Standard Deviation")
dev.off()

# Means 
mean.expression <- apply(matrix.org, MARGIN = 1, FUN = mean)
int.rank <- rank(1 / mean.expression)
png(filename = "means.png", width = 35, height = 20, units = "cm", res = 300)   # .png 2
plot(int.rank, mean.expression, xlab = "Genes", ylab = "Expression")
dev.off()


# indexes of variance and expression, which are above the threshold:
variance.sel <- which(variance > quantile(variance, 0.50))
expression.sel <- which(mean.expression > quantile(mean.expression, 0.50))

# Select those genes, which fulfill both criteria:
gene.sel <- intersect(variance.sel, expression.sel)
matrix.sel <- matrix.org[gene.sel, ]

# Plot the expression of the genes versus the standard deviation and highlight the selected genes:
png(filename = "Selection.png", width = 35, height = 25, units = "cm", res = 300)   # .png 3
plot(mean.expression, variance, xlab = "Expression", ylab = "Standard deviation", pch = ".")
points(mean.expression[gene.sel], variance[gene.sel], col = "red", pch = "+", cex = 0.7)
dev.off()




### Hierarchical clustering
method <- "average"
#method <- "single"
distance <- "euclidean"
#distance <- "manhattan"

# apply the hierarchical clustering method for the samples (microarrays) and display the associated dendrogram:
clustering <- hclust(dist(t(matrix.sel), method = distance), method = method)
png(filename = "dendrogram.png", width = 35, height = 25, units = "cm", res = 300)   # .png 4
plot(clustering)
dev.off()


# split the dendrogram into 4 clusters and obtain the members of each cluster:
groups <- cutree(clustering, 4)
samples <- colnames(matrix.sel)

for (i in c(1:4)) {
  print(samples[groups == i])
  print("--------------------")
}



# heatmap:
png(filename = "heatmap.png", width = 35, height = 25, units = "cm", res = 300)   # .png 5
heatmap.2(as.matrix(matrix.sel), Rowv = TRUE, Colv = as.dendrogram(clustering), key = TRUE, scale = "row", 
          col = redgreen(75), symkey = FALSE, density.info = "histogram", trace = "none", labRow = FALSE)

#heatmap.2(as.matrix(matrix.sel), Rowv = TRUE, Colv = as.dendrogram(clustering), key = TRUE, scale = "row", 
#          col = redgreen(75), symkey = FALSE, density.info = "histogram", trace = "none", 
#          labRow = rownames(as.matrix(matrix.sel[c(1:13),])))
dev.off()




### Cluster stability (sample clusters)
# Based on a bootstrapping approach, many datasets are "generated" by repeated drawing from the existing data. If 
# these many datasets generate the same clusters, then the clusters are "robust".
# p-values give an estimate how robust the clusters are supported by the data.
set.seed(5)
pv <- pvclust(matrix.sel, method.hclust = method, method.dist = distance, nboot = 10)
p.value <- 0.05
# Draw the clustering including the cluster stability and highlight the stable clusters with a red box:
png(filename = "Cluster_stability.png", width = 35, height = 25, units = "cm", res = 300)   # .png 6
plot(pv)
pvrect(pv, alpha= 1 - p.value)
dev.off()

# to remove outliers:
number_col_1 <- grep("ec_aer_appY_O_c", colnames(matrix.sel))
number_col_2 <- grep("ec_aer_fnr_nO_a", colnames(matrix.sel))
new.matrix <- matrix.sel[,-c(number_col_1, number_col_2)]
#new.matrix <- matrix.sel[,!(names(matrix.sel) %in% c("ec_aer_appY_O_c", "ec_aer_fnr_nO_a"))]

# hierarchical clustering without outliers:
new.clustering <- hclust(dist(t(new.matrix), method = distance), method = method)
png(filename = "dendrogram_no_outliers.png", width = 35, height = 25, units = "cm", res = 300)   # .png 7
plot(new.clustering)
dev.off()

# repeating the stability check after outlier removal:
set.seed(5)
pv <- pvclust(new.matrix, method.hclust = method, method.dist = distance, nboot = 10)
png(filename = "Cluster_stability_no_outliers.png", width = 35, height = 25, units = "cm", res = 300)   # .png 8
plot(pv)
pvrect(pv, alpha= 1 - p.value)
dev.off()


## Cluster stability (gene clusters)
set.seed(5)
pv <-  pvclust(t(matrix.sel), method.hclust = method, method.dist = distance, nboot = 10)
png(filename = "Cluster_genes.png", width = 60, height = 45, units = "cm", res = 300)   # .png 9
plot(pv)
dev.off()
# select stable clusters
pv.pp <- pvpick(pv, alpha= 1 - p.value)
n_cl  <-  length(pv.pp$clusters)
print(n_cl)

# Note: a low number of bootstraps (see nboot parameter of the 'pvclust()' function) doesn't return a reliable result.
# We can observe that the number of stable clusters is different when repeating the clustering of the genes a couple
# of times. You can avoid this, by increasing the nboot parameter to i.e. 100 or 1000 (but the clustering also takes
# much longer!):
#pv <- pvclust(t(matrix.sel), method.hclust = method, method.dist = distance, nboot = 1000)




### Interpretation (functional analysis)
# display the gene names of the chosen gene cluster:
#pv <- pvclust(t(matrix.sel), method.hclust = method, method.dist = distance, nboot = 100)
#pv.pp <- pvpick(pv, alpha= 1 - p.value)
cluster.sel  <- 1
pv.pp[[1]][[cluster.sel]]

# the same using 'hclust()' function:
clustering.genes <- hclust(dist(matrix.sel, method = distance), method = method)
png(filename = "Cluster_genes_2.png", width = 60, height = 45, units = "cm", res = 300)   # .png 10
plot(clustering.genes)
dev.off()

groups <- cutree(clustering.genes, k = 2)
gene <- rownames(matrix.sel)
genes_cluster <- gene[!groups == cluster.sel]   # Note: group numbers are different, compared to 'pvclust()'
print(genes_cluster)



### k-means clustering
# set the number of clusters to 2 and perform k-means clustering for all samples:
number_of_clusters <- 2
set.seed(7)
km <- kmeans(t(matrix.sel), centers = number_of_clusters)
table(km$cluster)

# mark the clusters of the k-means method with different colors
my.col <- palette()[2:3]
col <- my.col[as.vector(km$cluster)]
samp.col <- col[order(km$cluster)]

# heatmap:
png(filename = "heatmap_and_k_means_1.png", width = 35, height = 25, units = "cm", res = 300)   # .png 11
heatmap.2(as.matrix(matrix.sel), Rowv = FALSE, Colv = as.dendrogram(clustering), scale = "row", dendrogram = "col", 
          ColSideColors = samp.col, col = redgreen(75), density.info = "histogram", trace = "none", labRow = FALSE)
dev.off()


## re-setting the number of clusters and perform k-means clustering for samples (with outliers):
number_of_clusters <- 14
set.seed(7)
km <- kmeans(t(matrix.sel), centers = number_of_clusters, iter.max = 100, nstart = 25)
table(km$cluster)

# re-marking the clusters of the k-means method with different colors
my.col_25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
#pie(rep(1, 25), col = my.col_25)

col <- my.col_25[as.vector(km$cluster)]
samp.col <- col[order(km$cluster)]

# heatmap:
png(filename = "heatmap_and_k_means_2.png", width = 35, height = 25, units = "cm", res = 300)   # .png 12
heatmap.2(as.matrix(matrix.sel), Rowv = FALSE, Colv = as.dendrogram(clustering), scale = "row", dendrogram = "col", 
          ColSideColors = samp.col, col = redgreen(75), density.info = "histogram", trace = "none", labRow = FALSE)
dev.off()


## re-setting the number of clusters and perform k-means clustering for samples (without outliers)
# using samples "b" as initial cluster centers:
number_of_clusters <- 14
set.seed(7)
#km <- kmeans(t(new.matrix), centers=number_of_clusters, iter.max = 300000, nstart = 50)
km <- kmeans(t(new.matrix), centers = t(new.matrix[,c(grep("_b$", colnames(new.matrix), ignore.case = FALSE))]),
             iter.max = 100, nstart = 25)   #samples 'b' as initial cluster centers
table(km$cluster)

col <- my.col_25[as.vector(km$cluster)]
samp.col <- col[order(km$cluster)]

# heatmap:
png(filename = "heatmap_and_k_means_3.png", width = 35, height = 25, units = "cm", res = 300)   # .png 14
heatmap.2(as.matrix(new.matrix), Rowv = FALSE, Colv = as.dendrogram(new.clustering), scale = "row", dendrogram = "col", 
          ColSideColors = samp.col, col = redgreen(75), density.info = "histogram", trace = "none", labRow = FALSE)
dev.off()