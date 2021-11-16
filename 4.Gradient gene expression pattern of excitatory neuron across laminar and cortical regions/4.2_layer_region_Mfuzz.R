layer <- c("L23it","L45it","L5it","L6it","L56itcar3","L56np","L6b","L6ct")
for (i in layer) {
  assign(i,readRDS(paste0("F:/工作/猴脑项目/2104/0401/layer_average_rds/",i,"average.rds")))
}
library("Mfuzz")
layerlist <- list(L23it,L45it,L5it,L6it,L56itcar3,L56np,L6b,L6ct)
for (i in 1:8) {
  p1<-layerlist[[i]]$RNA
  p1_pseudo<-p1+0.000001
  ave_ZDD=p1_pseudo
  ave_ZDD$max_exp=apply(ave_ZDD,1,max)
  sub_ave_ZDD=subset(ave_ZDD,ave_ZDD$max_exp>1)
  sub_ave_ZDD$Gene=rownames(sub_ave_ZDD)
  p1_pseudo <- sub_ave_ZDD[,1:3]
  #加一个极小值
  count <- data.matrix(p1_pseudo)
  eset <- new("ExpressionSet",exprs = count)
  eset <- filter.std(eset,min.std=0)
  eset <- standardise(eset)
  m <- mestimate(eset)
  cl <- mfuzz(eset, c =6, m = m)
  pdf(paste0("F:/工作/猴脑项目/2104/0408/layer_region_pattern/",layer[i],"mfuzz.pdf"))
  mfuzz.plot(eset,cl,mfrow=c(3,2),new.window= FALSE,min.mem=0.6,time.labels = c("PFC","M1","V1"))
  dev.off()
  membership=as.data.frame(cl$membership)
  membership$max=apply(membership,1,max)
  membership_sub=subset(membership,membership$max>0.5)
  n=ncol(membership_sub)+1
  a=ncol(membership_sub)
  for(j in 1:nrow(membership_sub)){
    membership_sub[j,n]=paste0("Cluster",colnames(membership_sub)[which(membership_sub[j,]==membership_sub[j,a])[1]])
  }
  colnames(membership_sub)[n]="Cluster"
  membership_sub$Gene=rownames(membership_sub)
  
  write.table(membership_sub,file = paste0("F:/工作/猴脑项目/2104/0408/layer_region_pattern/",layer[i],"_Gene_ClusterList_Stage_new.csv"),sep = ",")
}
