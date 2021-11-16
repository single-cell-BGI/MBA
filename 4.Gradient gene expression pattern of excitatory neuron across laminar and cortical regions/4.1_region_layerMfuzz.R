region <- c("PFC","M1","V1")
for (i in region) {
  assign(i,readRDS(paste0("F:/工作/猴脑项目/2104/0401/region_average_rds/",i,"average.rds")))
}

library("Mfuzz")
regionlist <- list(PFC,M1,V1)
for (i in 1:3) {
  p1<-regionlist[[i]]$RNA
  p1 <- p1[,1:4]
  p1_pseudo<-p1+0.000001
  #加一个极小值
  ave_ZDD=p1_pseudo
  ave_ZDD$max_exp=apply(ave_ZDD,1,max)
  sub_ave_ZDD=subset(ave_ZDD,ave_ZDD$max_exp > 1)
  sub_ave_ZDD$Gene=rownames(sub_ave_ZDD)
  p1_pseudo <- sub_ave_ZDD[,1:4]
  count <- data.matrix(p1_pseudo)
  eset <- new("ExpressionSet",exprs = count)
  eset <- filter.std(eset,min.std=0)
  eset <- standardise(eset)
  m <- mestimate(eset)
  cl <- mfuzz(eset, c = 5, m = m)
  pdf(paste0("F:/工作/猴脑项目/2104/0408/region_layer_pattern/",region[i],"mfuzz.pdf"),width = 9,height = 6)
  mfuzz.plot(eset,cl,mfrow=c(2,3),new.window= FALSE,min.mem=0.6,time.labels = c("L23it","L45it","L5it","L6it"))
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
  
  write.table(membership_sub,file = paste0("F:/工作/猴脑项目/2104/0408/region_layer_pattern/",region[i],"_Gene_ClusterList_Stage_new.csv"),sep = ",")
}

region <- c("PFC","M1","V1")
for (i in region) {
  assign(i,read.table(paste0("F:/工作/猴脑项目/2104/0408/region_layer_pattern/",i,"_Gene_ClusterList_Stage_new.csv"),sep = ",",header = T,row.names = 1))
}


library(VennDiagram)
pfc_m1 <- PFC[PFC$Gene %in% M1$Gene,]
pfc_v1 <- PFC[PFC$Gene %in% V1$Gene,]
m1_v1 <- M1[M1$Gene %in% V1$Gene,]
pfc_m1_v1 <- pfc_m1[pfc_m1$Gene %in% V1$Gene,]
library(RColorBrewer)
library(futile.logger)
library(VennDiagram)
pdf("F:/工作/猴脑项目/2104/0409/mfuzz_3re_vein.pdf")
venn.plot <- draw.triple.venn(
  area1 = nrow(PFC),
  area2 = nrow(M1),
  area3 = nrow(V1),
  n12 = nrow(pfc_m1),
  n23 = nrow(m1_v1),
  n13 = nrow(pfc_v1),
  n123 = nrow(pfc_m1_v1),
  category = c("PFC", "M1", "V1"),
  fill = c("indianred1", "orange1", "lightblue2"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("indianred1", "orange1", "lightblue2")
)
dev.off()

