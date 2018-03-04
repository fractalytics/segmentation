library(FactoMineR)
library(tsne)

data=read.csv2("data-costumers.csv",row.names = 1,header = TRUE)

data=data[,1:8]



res.pca = PCA(data, scale.unit=TRUE, ncp=ncol(df), graph=F)
res.tsne=tsne(scale(data),k = 2,max_iter = 500,epoch = 100,perplexity = 65)

km.tsne=kmeans(res.tsne,3)

par(mfrow=c(1,2))
plot(res.tsne,col=km.tsne$cluster,main="tsne")
plot(res.pca$ind$coord[,1:2],main="pca")
