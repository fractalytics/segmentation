install.packages("FactoMineR")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("DT")
install.packages("factoextra")
install.packages("jsonlite")

library(FactoMineR)
library(factoextra)
library(ggplot2)
library(ggrepel)
library(DT)
library(jsonlite)



# data importation
df=read.csv("data/wellness.csv", sep=";",dec=",",row.names = 1, header= TRUE,check.names=FALSE)


################################################################################
#                                                                              #
#                      PCA, and factors map                                    #
#                                                                              #
################################################################################
res_acp = PCA(df, scale.unit=TRUE, ncp=ncol(df), graph=F)
df_acp=as.data.frame(res_acp$ind$coord)
fviz_pca(res_acp, col.var ="blue")



##############################################################################################
#                                                                                            #
#clustering, Hierarchical agglomerative clustering (HAC),on pca result                       #                                               
#                                                                                            #
##############################################################################################


CAH=hclust(dist(df_acp[, 1:4]),method="ward.D2")
plot(CAH)


################################################################################
#                                                                              #
#                            nb   optimal  cluster                             #
#                                                                              #
################################################################################


alpha= 0.20
opt_calc=fviz_nbclust(df_acp[,c(1,4)], hcut, method = "wss")

for ( i in 1:10){
   
  if(opt_calc$data$y[i+1]/opt_calc$data$y[1]<alpha){
    nb_optimal=i+1
    break
  }
}

fviz_nbclust(df_acp[,c(1,4)], hcut, method = "wss")+geom_vline(xintercept = nb_optimal, linetype = 2)
cluster=cutree(CAH, nb_optimal)
cluster=as.data.frame(cluster)
df_cluster=cbind(df_acp,cluster)


 


################################################################################
#                                                                              #
#      visualization of Individuals factor map  whith cluster assigned        #
#                                                                              #
################################################################################





set.seed(42)
ggplot(df_cluster) +
  geom_point(aes(Dim.1, Dim.2), size = 5, color = 'grey') +
  geom_label_repel(
    aes(Dim.1, Dim.2, fill = factor(cluster), label = rownames(df_cluster)),
    fontface = 'bold', color = 'white',
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines"),
    segment.color = 'grey50'
  ) +
  theme_classic(base_size = 16)


################################################################################
#                                                                              #
#                              cluster characterization                        #
#                                                                              #
################################################################################

z.test = function(k, g){
  
  n.k = length(k)
  n.g = length(g)
  Sk=((n.g-n.k)*var(g))/((n.g-1)*n.k)
  zeta = (mean(k) - mean(g)) / sqrt(Sk)
  
  return(zeta)
}

threshold=0.975
threshold.val=qnorm(threshold)

cluster_info <- data.frame(matrix(nrow=ncol(df)))
temp=data.frame()
df_c=cbind(df,cluster) 



for ( j in 1:nrow(unique(cluster))){ 
  cat("\n","cluster ", j, "\n", "\n", "\n" )
  for (i in 1:ncol(df)) {
    statistic=z.test(subset(df_c[,i], cluster==j),df_c[,i]) 
    
    temp=rbind(temp,statistic) 
    
    if (statistic >=threshold.val || statistic <= -threshold.val) {
      
      cat(colnames(df)[i],":",statistic[1],"\n")
      
    } 
    
  }
  colnames(temp)=paste("cluster",toString(j[1]) , sep=" ")  
  cluster_info=cbind(cluster_info,temp)
  temp=data.frame()
  
}






q=as.data.frame(apply(df_c,2,mean))
q=q[-nrow(q),]
mean=as.data.frame(q)
colnames(mean)=paste("average global")


clust=cbind(cluster_info,mean)

for ( i in 1:nrow(unique(cluster))){

  q=as.data.frame(apply(df_c[df_c$cluster==i,],2,mean))
  q=q[-nrow(q),]
  mean=as.data.frame(q)
  colnames(mean)=paste("moyenne cluster",toString(i[1]))
  clust=cbind(clust,mean)
}  




cluster_info=cluster_info[,1:nrow(unique(cluster))+1]

rownames(cluster_info)=colnames(df)




color_from_middle <- function (data, color1,color2) 
{
  max_val=max(abs(data))
  JS(sprintf("isNaN(parseFloat(value)) || value < 0 ? 'linear-gradient(90deg, transparent, transparent ' + (50 + value/%s * 50) + '%%, %s ' + (50 + value/%s * 50) + '%%,%s  50%%,transparent 50%%)': 'linear-gradient(90deg, transparent, transparent 50%%, %s 50%%, %s ' + (50 + value/%s * 50) + '%%, transparent ' + (50 + value/%s * 50) + '%%)'",
             max_val,color1,max_val,color1,color2,color2,max_val,max_val))
} 



for ( i in 1:nrow(unique(cluster))){
  dfc=subset(cluster_info[abs(cluster_info[,i])>threshold.val,], select = paste("cluster",toString(i[1]) , sep=" "))
  
  
  a=datatable(round(dfc,2)) %>%
    formatStyle(paste("cluster",toString(i[1]) , sep=" "),background=color_from_middle(dfc[,1],'tomato','deepskyblue'))
  print(a)
}
  


################################################################################
#                                                                              #
#                            individuals and distance                          #
#                                                                              #
################################################################################



for(i in 1:nb_optimal){
  
  
  k=df_cluster[df_cluster$cluster==i,c('Dim.1','Dim.2','cluster')]
  
  Distance=sqrt((k$Dim.1-mean(k$Dim.1))^2+(k$Dim.2-mean(k$Dim.2))^2)
  
  
  
  classe=as.data.frame(Distance,row.names = rownames(k))
  
  classe=classe[order(classe['Distance']),, drop = FALSE] 
  
  cat("\n\n\t","cluster:",i,"\n\n")
  print(classe)
  
}

      
################################################################################
#                                                                              #
#                             json generation for D3.js                        #
#                                                                              #
################################################################################


path_json="viz/outfile.json"
 
if (file.exists(path_json)) file.remove(path_json)

cat("{",'"features":',"[",file=path_json,append=TRUE)
for(i in 1:ncol(cluster_info)){
  
  cat('{\n"cluster":','"',i,'",','\n','"vars"',": [",'\n',file=path_json,append=TRUE)
  cluster_info=cluster_info[order(cluster_info[,i],decreasing = TRUE), , drop = FALSE]
  for(z in 1:nrow(cluster_info)){
    
    
    if(abs(cluster_info[z,i])>=threshold.val){
      cat('{','\n','"label" :','"',rownames(cluster_info)[z],'",','\n','"value":',round(cluster_info[z,i],digit=2),'\n }',file=path_json,append=TRUE)
      if(max(which(abs(cluster_info[,i])>=threshold.val))!=z) { cat(",",file=path_json,append=TRUE) }
      cat('\n',file=path_json,append=TRUE)
    }
  }
  
  cat("] \n } ",file=path_json,append=TRUE)
  if(ncol(cluster_info)!=i) { cat(",",file=path_json,append=TRUE) }
  cat("\n"  ,file=path_json,append=TRUE)
}
cat("],",'"observations": [',file=path_json,append=TRUE)
res1=df_cluster[, c(1,2,length(df_cluster))]
res1=cbind(ENTRY=rownames(df_cluster),df_cluster[, c(1,2,length(df_cluster))])
colnames(res1) <- c("name","Dim_1","Dim_2","cluster")
row.names(res1)<- NULL
x <- toJSON(unname(split(res1, 1:nrow(res1))))
x=gsub("[[]", "", x)
x=gsub("[]]", "", x)
cat(x,file=path_json,append=TRUE)
cat("\n ] \n }",file=path_json,append=TRUE)





