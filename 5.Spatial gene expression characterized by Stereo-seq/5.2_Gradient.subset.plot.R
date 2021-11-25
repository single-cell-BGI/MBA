library(Seurat)
library(RColorBrewer)
library(data.table)
library(ggplot2)
ps<-brewer.pal(8,"Set1")
ps<-ps[-6]

##############################################PFC############################################################
rds<-readRDS("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/8.Macaca_brain/4.gradient/data/PFC_DZ24_SPOTlight.rds")
rds@meta.data$coor_x<-as.numeric(rds@meta.data$coor_x)
rds@meta.data$coor_y<-as.numeric(rds@meta.data$coor_y)

#tmp<-rds@meta.data
#tmp$coor_x<-as.numeric(tmp$coor_x)
#tmp$coor_y<-as.numeric(tmp$coor_y)

rds<-subset(rds,subset= coor_x >250)
rds<-subset(rds,subset= coor_y >75)
rds<-subset(rds,subset= coor_y <100)

#saveRDS(rds,file="PFC.cutlayer.rds")

EX<-c("L23it","L4it","L5.6.IT.Car3","L56np","L5et","L5it","L6b","L6ct","L6it")
rds
rds<-subset(rds,subset=CellType %in% EX)
rds

saveRDS(rds,file="PFC.cutlayer.rds")

tmp<-rds@meta.data
tmp$coor_x<-as.numeric(tmp$coor_x)
tmp$coor_y<-as.numeric(tmp$coor_y)
tmp$coor_x<-(-1*tmp$coor_x)
tmp<-as.data.table(tmp)

genelist<-readLines(con="/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/8.Macaca_brain/4.gradient/script/gene.list")
genelist<-genelist[genelist %in% rownames(rds)]
data.list<-list()
p<-list()
pdf("PFC.gene.sub.list_final.pdf")

for (i in c(1:length(genelist))){
tmp$Gene_exp<-as.numeric(rds@assays$SCT@data[paste(genelist[i]),])
data.list[[i]]<-tmp[, mean(Gene_exp), by = .(coor_x)]
p[[i]]<-ggplot(data.list[[i]],aes(x=as.numeric(coor_x),data.list[[i]]$V1))+geom_smooth(colour=ps[1],size=2)+theme_bw()+ggtitle(genelist[[i]])
print(p[[i]])
}
dev.off()


write.table(tmp,file="PFC.list.xls",quote=F,sep="\t")
pk<-ggplot(tmp,aes(x=as.numeric(coor_x),y=as.numeric(coor_y),color=CellType))+scale_color_manual(values=ps)+geom_point()+theme_bw()
pdf("PFC.layer.sub.list2.pdf")
pk
dev.off()

##########################################################################################################

##############################################M1############################################################
rds<-readRDS("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/8.Macaca_brain/4.gradient/data/M1_DZ22_SPOTlight.rds")
rds@meta.data$coor_x<-as.numeric(rds@meta.data$coor_x)
rds@meta.data$coor_y<-as.numeric(rds@meta.data$coor_y)

#tmp<-rds@meta.data
#tmp$coor_x<-as.numeric(tmp$coor_x)
#tmp$coor_y<-as.numeric(tmp$coor_y)

rds<-subset(rds,subset= coor_x <160)
rds<-subset(rds,subset= coor_y >275)
rds<-subset(rds,subset= coor_y <300)

#rds$coor_x<-(-1*rds$coor_x)

EX<-c("L23it","L4it","L5.6.IT.Car3","L56np","L5et","L5it","L6b","L6ct","L6it")
rds
rds<-subset(rds,subset=CellType %in% EX)
rds

saveRDS(rds,file="M1.cutlayer.rds")

tmp<-rds@meta.data
tmp$coor_x<-as.numeric(tmp$coor_x)
tmp$coor_y<-as.numeric(tmp$coor_y)
#tmp$coor_x<-(-1*tmp$coor_x)

tmp<-as.data.table(tmp)

genelist<-readLines(con="/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/8.Macaca_brain/4.gradient/script/gene.list")
genelist<-genelist[genelist %in% rownames(rds)]
data.list<-list()
p<-list()
pdf("M1.gene.sub.list2.pdf")

for (i in c(1:length(genelist))){
tmp$Gene_exp<-as.numeric(rds@assays$SCT@data[paste(genelist[i]),])
data.list[[i]]<-tmp[, mean(Gene_exp), by = .(coor_x)]
p[[i]]<-ggplot(data.list[[i]],aes(x=as.numeric(coor_x),data.list[[i]]$V1))+geom_smooth(colour=ps[1],size=2)+theme_bw()+ggtitle(genelist[[i]])
print(p[[i]])
}
dev.off()


pk<-ggplot(tmp,aes(x=as.numeric(coor_x),y=as.numeric(coor_y),color=CellType))+scale_color_manual(values=ps)+geom_point()+theme_bw()
pdf("M1.layer.sub.list2.pdf")
pk
dev.off()
##########################################################################################################

##############################################V1##########################################################
rds<-readRDS("/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/8.Macaca_brain/4.gradient/data/V1_DT21_SPOTlight.rds")
rds@meta.data$coor_x<-as.numeric(rds@meta.data$coor_x)
rds@meta.data$coor_y<-as.numeric(rds@meta.data$coor_y)

#tmp<-rds@meta.data
#tmp$coor_x<-as.numeric(tmp$coor_x)
#tmp$coor_y<-as.numeric(tmp$coor_y)

rds<-subset(rds,subset= coor_x >225)
rds<-subset(rds,subset= coor_y >250)
rds<-subset(rds,subset= coor_y <275)

EX<-c("L23it","L4it","L5.6.IT.Car3","L56np","L5et","L5it","L6b","L6ct","L6it")
rds
rds<-subset(rds,subset=CellType %in% EX)
rds

saveRDS(rds,file="V1.cutlayer.rds")

tmp<-rds@meta.data
tmp$coor_x<-as.numeric(tmp$coor_x)
tmp$coor_y<-as.numeric(tmp$coor_y)
tmp$coor_x<-(-1*tmp$coor_x)
tmp<-as.data.table(tmp)

genelist<-readLines(con="/jdfssz1/ST_SUPERCELLS/P20Z10200N0059/8.Macaca_brain/4.gradient/script/gene.list")
genelist<-genelist[genelist %in% rownames(rds)]
data.list<-list()
p<-list()
pdf("V1.gene.sub.list_final.pdf")

for (i in c(1:length(genelist))){
tmp$Gene_exp<-as.numeric(rds@assays$SCT@data[paste(genelist[i]),])
data.list[[i]]<-tmp[, mean(Gene_exp), by = .(coor_x)]
p[[i]]<-ggplot(data.list[[i]],aes(x=as.numeric(coor_x),data.list[[i]]$V1))+geom_smooth(colour=ps[1],size=2)+theme_bw()+ggtitle(genelist[[i]])
print(p[[i]])
}
dev.off()

write.table(tmp,file="V1.list.xls",quote=F,sep="\t")
pk<-ggplot(tmp,aes(x=as.numeric(coor_x),y=as.numeric(coor_y),color=CellType))+scale_color_manual(values=ps)+geom_point()+theme_bw()
pdf("V1.layer.sub.list2.pdf")
pk
dev.off()

##########################################################################################################

