library(clValid)
library(NbClust)
library(factoextra)

par(mar=c(1,1,1,1))

data = read.table('./output/norm_condtion_df.csv', sep=",", row.names=1,
                  header=TRUE)
#intern <- clValid(data, seq(from = 10, to = 20, length.out = 10), 
#                  clMethods=c("hierarchical","kmeans","pam"),
#                  validation=c("internal"))

q = NbClust(data, diss = NULL, distance = "euclidean",
        min.nc = 5, max.nc = 30, method = 'ward.D')

fviz_nbclust(data, hcut, method=c( "wss"), k.min=3, k.max=25)
