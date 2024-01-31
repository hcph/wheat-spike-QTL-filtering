##############################integrate network########################################
##PPI
ppi <- read.csv("s1s2_compare_result.csv")
#add symbol
sy <- read.table("D:/小麦信息/motif2gene_mappingR3.txt",sep="\t")
colnames(sy) <- c("source_symbol","source")
ppi1 <- merge(ppi,sy,by="source",all=F)
colnames(sy) <- c("target_symbol","target")
ppi2 <- merge(ppi1,sy,by="target",all=F)
ppi2 <- ppi2[,c(2,1,9,10,6)]
ppi2$type <- gsub("s1_unique","S1",ppi2$type)
ppi2$type <- gsub("s2_unique","S2",ppi2$type)
ppi2$type <- gsub("s1s2_same","same",ppi2$type)
ppi2$group <- "PPI"
head(ppi2)
##TF network
ft <- read.csv("F印迹网络_S1_S2_compare.csv",row.names = 1)
ft$group <- "TFfootprint"
head(ft)
##eQTLs network
eq <- read.csv("落在OCR的eqtls.csv")
eq <- eq[,c(4,5)]
colnames(eq) <- c("source","target")
colnames(sy) <- c("source_symbol","source")
eq1 <- merge(eq,sy,by="source",all=F)
colnames(sy) <- c("target_symbol","target")
eq2 <- merge(eq1,sy,by="target",all=F)
eq2 <- eq2[,c(2,1,3,4)]
eq2$type <- "same"
eq2$group <- "eQTL"
###integrate the two network
tt <- rbind(ppi2,ft,eq2)
head(tt)
write.csv(tt,"total_network.csv",row.names=F,quote=F)
#分别导出网络
write.csv(unique(tt[which(tt$type == "S1" | tt$type == "same"),]),"SM_total_network.csv",row.names=F,quote=F)
write.csv(unique(tt[which(tt$type == "S2" | tt$type == "same"),]),"FM_total_network.csv",row.names=F,quote=F)
dim(unique(tt[which(tt$type == "S1" | tt$type == "same"),]))
dim(unique(tt[which(tt$type == "S2" | tt$type == "same"),]))

#################################different stages（PPI, TFGRN, eQTL）###########################################
net <- read.csv("total_network.csv")
head(net)
###S1
net1 <- net[which(net$type == "S1" | net$type == "same"),]
net1 <- net1[which(net1$group != "eQTL"),]
write.csv(net1,"S1_total_network.csv",row.names = F,quote=F)
a1 <- aggregate(net1$source, by=list(net1$group,net1$source), FUN=length)
a2 <- aggregate(net1$source, by=list(net1$group,net1$target), FUN=length)
head(a1)
head(a2)
colnames(a1) <- c("group","source","source_degree")
colnames(a2) <- c("group","target","target_degree")
a1$fz <- paste(a1$group,a1$source,sep="_")
a2$fz <- paste(a2$group,a2$target,sep="_")
a <- merge(a1,a2,by="fz",all=T)
head(a)
a[is.na(a)] = 0
b1 <- a[which(a$group.x == 0),]
b2 <- a[which(a$group.x != 0),]
head(b1)
head(b2)
b1$group.x <- b1$group.y
b <- rbind(b1,b2)
b1 <- b[which(b$source == 0),]
b2 <- b[which(b$source != 0),]
head(b1)
head(b2)
b1$source <- b1$target
b <- rbind(b1,b2)
head(b)
c1 <- b[which(b$group.x != "PPI"),]
c2 <- b[which(b$group.x == "PPI"),]
head(c1)
#calculate indegree
d1 <- aggregate(c1$source_degree,by=list(c1$source),FUN=sum)
head(d1)
d2 <- aggregate(c1$target_degree,by=list(c1$source),FUN=sum)
head(d2)
colnames(d1) <- c("gene","outdegree")
colnames(d2) <- c("gene","indegree")
head(c2)
d3 <- c2[,c(3,4,7)]
d3$nodirect <- d3$source_degree+d3$target_degree
head(d3)
d3 <- d3[,-2:-3]
colnames(d3) <- c("gene","nodirect")
e1 <- merge(d1,d2,by="gene",all=T)
head(e1)
e2 <- merge(e1,d3,by="gene",all=T)
head(e2)
e2[is.na(e2)] = 0
write.csv(e2,"s1_gene_degree.csv")
###S2
net2 <- net[which(net$type == "S2" | net$type == "same"),]
net2 <- net2[which(net2$group != "eQTL"),]
write.csv(net2,"S2_total_network.csv",row.names = F,quote=F)
a1 <- aggregate(net2$source, by=list(net2$group,net2$source), FUN=length)
a2 <- aggregate(net2$source, by=list(net2$group,net2$target), FUN=length)
head(a1)
head(a2)
colnames(a1) <- c("group","source","source_degree")
colnames(a2) <- c("group","target","target_degree")
a1$fz <- paste(a1$group,a1$source,sep="_")
a2$fz <- paste(a2$group,a2$target,sep="_")
a <- merge(a1,a2,by="fz",all=T)
head(a)
a[is.na(a)] = 0
b1 <- a[which(a$group.x == 0),]
b2 <- a[which(a$group.x != 0),]
head(b1)
head(b2)
b1$group.x <- b1$group.y
b <- rbind(b1,b2)
b1 <- b[which(b$source == 0),]
b2 <- b[which(b$source != 0),]
head(b1)
head(b2)
b1$source <- b1$target
b <- rbind(b1,b2)
head(b)
c1 <- b[which(b$group.x != "PPI"),]
c2 <- b[which(b$group.x == "PPI"),]
head(c1)
#calculate indegree
d1 <- aggregate(c1$source_degree,by=list(c1$source),FUN=sum)
head(d1)
d2 <- aggregate(c1$target_degree,by=list(c1$source),FUN=sum)
head(d2)
colnames(d1) <- c("gene","outdegree")
colnames(d2) <- c("gene","indegree")
head(c2)
d3 <- c2[,c(3,4,7)]
d3$nodirect <- d3$source_degree+d3$target_degree
head(d3)
d3 <- d3[,-2:-3]
colnames(d3) <- c("gene","nodirect")
e1 <- merge(d1,d2,by="gene",all=T)
head(e1)
e2 <- merge(e1,d3,by="gene",all=T)
head(e2)
e2[is.na(e2)] = 0
write.csv(e2,"s2_gene_degree.csv")

#################################################produce node file##################################################
##node
net <- read.csv("total_network.csv")
x1 <- net[,c(1,3)]
x2 <- net[,c(2,4)]
colnames(x2) <- colnames(x1)
node <- rbind(x1,x2)
head(node)
##add known functional genes
library(openxlsx)
know <- read.xlsx("最终穗部性状基因R3.xlsx")
know <- know[,c(1,3)]
colnames(know) <- c("source","Gene_name")
head(know)
nodee <- unique(merge(node,know,by="source",all.x=T))
head(nodee)
dim(nodee)
nodee[is.na(nodee)] <- "no"
x1 <- nodee[which(nodee$Gene_name == "no"),]
x1$known <- "no"
x2 <- nodee[which(nodee$Gene_name != "no"),]
x2$known <- "yes"
nodeee <- rbind(x1,x2)
dim(nodeee)
##add degree
#merge
e1 <- read.csv("s1_gene_degree.csv",row.names=1)
e2 <- read.csv("s2_gene_degree.csv",row.names = 1)
colnames(e1) <- c("source","SM_outdegree","SM_indegree","SM_nodirect")
colnames(e2) <- c("source","FM_outdegree","FM_indegree","FM_nodirect")
node1 <- merge(nodeee,e1,by="source",all.x=T)
node2 <- merge(node1,e2,by="source",all.x=T)
node2[is.na(node2)] <- 0
head(node2)
##add putative 
putative <- read.csv("正向区段及对应的表达基因.csv")
putative$putative <- "yes"
putative <- putative[,c(2,22)]
colnames(putative) <- c("source","putative")
node2 <- merge(node2,putative,by="source",all.x=TRUE)
node2[is.na(node2)] <- "no"
head(node2)
##add putative hapotype information
hy <- read.csv("total_genev1.1.csv",row.names=1)
hy$haplotype <- "significate"
colnames(hy) <- c("source","triat","haplotype")
head(hy)
nodeeee <- merge(node2,hy,by="source",all.x=TRUE)
head(nodeeee) 
nodeeee[is.na(nodeeee)] <- "no"
write.csv(nodeeee,"node_total_network.csv",row.names = F)
head(nodeeee)
##add degree by stages
node1 <- nodeeee[,c(1,2,3,4,5,6,7,11,12,13)]
node2 <- nodeeee[,c(1,2,3,4,8,9,10,11,12,13)]
head(node1)
head(node2)
##add hub
i=200
hub1 <- node1[which(node1$SM_outdegree > i | node1$SM_indegree > i),]
hub2 <- node2[which(node2$FM_outdegree > i | node2$FM_indegree > i),]
head(hub1)
head(hub2)
dim(hub1)
dim(hub2)
#save result
write.csv(hub1,"s1_hub_node.csv",row.names = F)
write.csv(hub2,"s2_hub_node.csv",row.names = F)
#add hub into node
hub1$hub <- "yes"
hub2$hub <- "yes"
hub1 <- hub1[,c(1,11)]
hub2 <- hub2[,c(1,11)]
node1 <- merge(node1,hub1,by="source",all.x=TRUE)
node2 <- merge(node2,hub2,by="source",all.x=TRUE)
node1[is.na(node1)] <- "no"
node2[is.na(node2)] <- "no"
head(node1)
head(node2)
write.csv(node1,"S1_node_total_network.csv",row.names = F)
write.csv(node2,"S2_node_total_network.csv",row.names = F)

###########################heatmap - degree higer than 200#####################
i=200
de1 <- read.csv("S1_node_total_network.csv",check.names = F)
de2 <- read.csv("S2_node_total_network.csv",check.names = F)
head(de1)
head(de2)
de12t <- merge(de1,de2,by="source",all=T)
de12t[is.na(de12t)] <- 0
head(de12t)
a1 <- de12t[which(de12t$Gene_name.x == "0" | de12t$Gene_name.x == "no"),]
a1$Gene_name.x <-  a1$Gene_name.y
a2 <- de12t[which(de12t$Gene_name.x != "0" & de12t$Gene_name.x != "no"),]
de12t <- rbind(a1,a2)
head(de12t)
de12 <- de12t[,c(1,5,6,15,16)]
head(de12)
de12$max1 <- rowMaxs(as.matrix(de12[,2:3]))
de12$max2 <- rowMaxs(as.matrix(de12[,4:5]))
head(de12)
de12 <- de12[which(de12$max1 > i | de12$max2 > i),1:5]
head(de12)
dim(de12)
row.names(de12) <- de12[,1]
de12 <- de12[,-1]
head(de12)
de12 <- log10(de12+1)
colnames(de12) <- c("SM_outdegree","SM_indegree","FM_outdegree","FM_indegree")
heat <- Heatmap(as.matrix(de12),
                #km=2,
                col = colorRampPalette(colors = c("white","red"))(20), #定义热图由低值到高值的渐变颜色
                heatmap_legend_param = list(grid_height = unit(10,'mm')),  #图例高度设置
                show_row_names = F,
                cluster_columns = F,
                cluster_rows = T)
#dev.off()
name <- de12t[,c(1,3)]
colnames(name) <- c("gene","symbol")
head(name)
row.names(name) <- name[,1]
name <- unique(name[which(name$symbol != 0 & name$symbol != "no"),])
de12name <- merge(de12,name,by="row.names",all=F)
head(de12name)
de12name <- de12name[,c(6,7)]
head(de12name)
#name <- as.data.frame(x[1:2,1])
#colnames(name)<-c("v1")
colnames(de12name)<-c("gene","symbol")
pdf("hub_gene_degree.pdf",family="ArialMT")
heat + rowAnnotation(link = anno_mark(at = which(rownames(de12) %in% de12name$gene), 
                                      #labels = paste(name$gene,name$symbol,sep="-"),
                                      labels = de12name$symbol, 
                                      labels_gp = gpar(fontsize = 10)))
dev.off()


####################################number of 的putative gene in the integrated network###############################
#s1
net <- read.csv("total_network.csv")
head(net)
###S1
net1 <- net[which(net$type == "S1" | net$type == "same"),]
net1 <- net1[which(net1$group != "eQTL"),]
head(net1)
ppi1 <- unique(c(net1[which(net1$group == "PPI"),1],net1[which(net1$group == "PPI"),2]))
tfgrn1 <- unique(c(net1[which(net1$group == "TFfootprint"),1],net1[which(net1$group == "TFfootprint"),2]))
net2 <- net[which(net$type == "S2" | net$type == "same"),]
net2 <- net2[which(net2$group != "eQTL"),]
head(net2)
ppi2 <- unique(c(net2[which(net2$group == "PPI"),1],net2[which(net2$group == "PPI"),2]))
tfgrn2 <- unique(c(net2[which(net2$group == "TFfootprint"),1],net2[which(net2$group == "TFfootprint"),2]))
head(ppi1)
head(tfgrn1)
putative <- read.csv("正向区段及对应的表达基因.csv")
head(putative)
library(VennDiagram)
venn_list <- list(SM_TFfootprint = tfgrn1,SM_PPI=ppi1,Putative = unique(putative$V7),FM_TFfootprint = tfgrn2,FM_PPI =ppi2)
venn.plot<-venn.diagram(venn_list, filename = NULL, 
                        fill = c('#7FC97F', "#B2DF8A","#ab549c",'#FF7F00', "#FDBF6F"), alpha = 0.50, 
                        cat.col = c('#7FC97F', "#B2DF8A","#ab549c",'#FF7F00', "#FDBF6F"), cat.cex = 1.5, cat.fontfamily = 'serif',
                        col = c('#7FC97F', "#B2DF8A","#ab549c",'#FF7F00', "#FDBF6F"), cex = 1.5, fontfamily = 'serif')
pdf("SM_FM_putative_overlapR2.pdf",family="ArialMT")
grid.draw(venn.plot)
dev.off()
#SM
putative <- putative[,c(1,2)]
colnames(putative) <- c("x","source")
ne1_p <- merge(net1,putative,by="source",all=F)
colnames(putative) <- c("x","target")
ne1_q <- merge(net1,putative,by="target",all=F)
head(ne1_p)
head(ne1_q)
ne1_p <- ne1_p[,c(1,6)]
ne1_q <- ne1_q[,c(1,6)]
colnames(ne1_p) <- colnames(ne1_q)
ne1 <- unique(rbind(ne1_p,ne1_q))
head(ne1)
ne1_number <- aggregate(ne1$target,list(ne1$group),length)
head(ne1_number) #PPI 1269  TFGRN 2765
#FM
colnames(putative) <- c("x","source")
ne2_p <- merge(net2,putative,by="source",all=F)
colnames(putative) <- c("x","target")
ne2_q <- merge(net2,putative,by="target",all=F)
head(ne2_p)
head(ne2_q)
ne2_p <- ne2_p[,c(1,6)]
ne2_q <- ne2_q[,c(1,6)]
colnames(ne2_p) <- colnames(ne2_q)
ne2 <- unique(rbind(ne2_p,ne2_q))
head(ne2)
ne2_number <- aggregate(ne2$target,list(ne2$group),length)
head(ne2_number) #PPI 1285 TFGRN 2572

############################################## generate the nwtwork associated with putative genes #####################################################################
net <- read.csv("total_network.csv")
#related to putative genes
putative <- read.csv("正向区段及对应的表达基因.csv")
putative <- data.frame(putative[,c(2)])
colnames(putative) <- "source"
pn1 <- merge(net,putative,by="source",all=F)
colnames(putative) <- "target"
pn2 <- merge(net,putative,by="target",all=F)
pn2 <- pn2[,c(2,1,3,4,5,6)]
pnt <- unique(rbind(pn1,pn2))
head(pnt)
write.csv(pnt,"putative_related_network.csv",row.names = F)
##count
#PPI 
length(unique(pnt[which((pnt$type == "same" | pnt$type == "S1" ) & (pnt$group == "PPI")),1])) #1843
length(unique(pnt[which((pnt$type == "same" | pnt$type == "S1" ) & (pnt$group == "PPI")),2])) #1555
nrow(unique(pnt[which((pnt$type == "same" | pnt$type == "S1" ) & (pnt$group == "PPI")),]))  #3841

length(unique(pnt[which((pnt$type == "same" | pnt$type == "S2" ) & (pnt$group == "PPI")),1])) #1890
length(unique(pnt[which((pnt$type == "same" | pnt$type == "S2" ) & (pnt$group == "PPI")),2])) #1670
nrow(unique(pnt[which((pnt$type == "same" | pnt$type == "S2" ) & (pnt$group == "PPI")),]))  #4165

length(unique(pnt[which((pnt$type == "same" | pnt$type == "S1" ) & (pnt$group == "TFfootprint")),1])) #762
length(unique(pnt[which((pnt$type == "same" | pnt$type == "S1" ) & (pnt$group == "TFfootprint")),2])) #12655
nrow(unique(pnt[which((pnt$type == "same" | pnt$type == "S1" ) & (pnt$group == "TFfootprint")),]))  #290908

length(unique(pnt[which((pnt$type == "same" | pnt$type == "S2" ) & (pnt$group == "TFfootprint")),1])) #738
length(unique(pnt[which((pnt$type == "same" | pnt$type == "S2" ) & (pnt$group == "TFfootprint")),2])) #11572
nrow(unique(pnt[which((pnt$type == "same" | pnt$type == "S2" ) & (pnt$group == "TFfootprint")),]))  #241903


#################################### filtering the above network using haplotype signicant################
putative <- read.csv("正向区段及对应的表达基因.csv")
putative <- as.data.frame(putative$V7)
pnt <- read.csv("putative_related_network.csv")
s1_putative <- pnt[which((pnt$type == "same" | pnt$type == "S1" )),]
s2_putative <- pnt[which((pnt$type == "same" | pnt$type == "S2" )),]
#SM
s1_source <- data.frame(s1_putative$source)
s1_source$type <- "source"
s1_target <- data.frame(s1_putative$target)
s1_target$type <- "target"
colnames(s1_source) <- colnames(s1_target) <- c("source","type")
s1_node <- rbind(s1_source,s1_target)
s1_node$sample <- "SM"
#FM
s2_source <- data.frame(s2_putative$source)
s2_source$type <- "source"
s2_target <- data.frame(s2_putative$target)
s2_target$type <- "target"
colnames(s2_source) <- colnames(s2_target) <- c("source","type")
s2_node <- rbind(s2_source,s2_target)
s2_node$sample <- "FM"
#merge
s12_node <- unique(rbind(s1_node,s2_node))
head(s12_node)
#add whether significant
hy <- read.csv("../GWAS/total_genev1.1.csv",row.names=1)
hy$haplotype <- "significate"
colnames(hy) <- c("source","triat","haplotype")
hy <- hy[,c(1,3)]
head(hy)
s12_node_hy <- merge(s12_node,hy,by="source",all.x=TRUE)
s12_node_hy[is.na(s12_node_hy)] <- "no"
s12_node_hy <- unique(s12_node_hy[,-2])
head(s12_node_hy)
#add whether located in QTL regions
colnames(putative) <- "source"
s12_node_hy_p <- unique(merge(s12_node_hy,putative,by="source",all=F))
head(s12_node_hy_p)
#count
s12_node_hy_num <- aggregate(s12_node_hy_p$source,list(s12_node_hy_p$sample,s12_node_hy_p$haplotype),length)
head(s12_node_hy_num)
pdf("网络筛选后的putative基因的单倍型显著性情况.pdf",family="ArialMT")
ggplot(s12_node_hy_num,mapping = aes(x=Group.1,y=x,fill=Group.2))+
  geom_bar(stat='identity',position='stack',show.legend = TRUE) +
  #geom_bar(stat='identity',position = position_dodge(),show.legend = TRUE) +
  labs(y = 'Number of gene') +
  scale_fill_manual(values=brewer.pal(10,"Paired")) +
  #scale_fill_manual(values=c("#CAB2D6","#6A3D9A")) +
  theme(panel.grid=element_blank(), panel.background=element_rect(fill='transparent', color='black'),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5)) +
  #guides(fill = guide_legend(title = 'Type')) +
  geom_text(aes(label = x), size = 3, hjust = 0.5, vjust = 1, position = "stack")
dev.off()
###count the number of the signicant genes
aa <- data.frame(unique(s12_node_hy_p[which(s12_node_hy_p$haplotype == "significate"),1]))
colnames(aa) <- "source"
node <- read.csv("node_total_networkR1.csv")
aa <-merge(aa,node,by="source",all=F)
head(aa)
bb <- aa[which(aa$hub == "yes" & aa$putative == "yes"),]
dim(unique(aa[,c(15,16,17)]))
write.table(aa,"TableS31.txt",row.names = F, quote = F,sep="\t")

#################################Generation of sub-networks of haplotype-significant and known hub gene (8) associations known to be located within segments###########################
##Extraction of hub genes in positive segments with significant haplotypes
node <- read.csv("node_total_networkR1.csv")
hub_node <- node[which(node$hub == "yes" & node$putative == "yes" & node$Gene_name != "no"),]
hub_node <- as.data.frame(hub_node$source)
###提取hub基因关联的网络
##SM
net <- read.csv("putative_related_network.csv")
colnames(hub_node) <- "source"
net1_a <- merge(net[which(net$type != "S2"),],hub_node,by="source",all=F)
colnames(hub_node) <- "target"
net1_b <- merge(net[which(net$type != "S2"),],hub_node,by="target",all=F)
net1_b <- net1_b[,c(2,1,3,4,5,6)]
net1 <- rbind(net1_a,net1_b)
dim(net1)
write.csv(net1,"s1_new_networkR2.csv",row.names = F,quote = F)
##FM
colnames(hub_node) <- "source"
net2_a <- merge(net[which(net$type != "S1"),],hub_node,by="source",all=F)
colnames(hub_node) <- "target"
net2_b <- merge(net[which(net$type != "S1"),],hub_node,by="target",all=F)
net2_b <- net1_b[,c(2,1,3,4,5,6)]
net2 <- rbind(net2_a,net2_b)
dim(net2)
write.csv(net2,"s2_new_networkR2.csv",row.names = F,quote = F)




################## 8个hub基因关联的网络中，有多少基因是hub（区间内单倍型显著169个）################
setwd("D:/AG_幼穗/第一轮审稿/总网络")
node <- read.csv("node_total_networkR1.csv")
hub_sig <- node[which(node$hub == "yes" & node$putative == "yes" & node$haplotype == "significate"),]
write.csv(hub_sig,"hub_putative_haplotype.csv",row.names=F)
hub_net1 <- read.csv("s1_new_networkR2.csv")
hub_net2 <- read.csv("s2_new_networkR2.csv")
head(hub_net1)
#####绘制韦恩图,
library(VennDiagram)
venn_list <- list(SM = unique(c(hub_net1$source,hub_net1$target)),
                  #Significate = tp, 
                  Hub_significate=hub_sig$source,
                  FM = unique(c(hub_net2$source,hub_net2$target)))
venn.plot<-venn.diagram(venn_list, filename = NULL, 
                        fill = c('#64cdcc', "#1a9850",'#9f8f12'), alpha = 0.50, 
                        cat.col = c('#64cdcc', "#1a9850",'#9f8f12'), cat.cex = 1.5, cat.fontfamily = 'serif',
                        col = c('#64cdcc', "#1a9850",'#9f8f12'), cex = 1.5, fontfamily = 'serif')
pdf("putative_hub_8hub_overlap.pdf",family="ArialMT")
grid.draw(venn.plot)
dev.off()


#####合并的总网络与郭的overlap
setwd("D:/AG_幼穗/第一轮审稿/总网络")
net <- read.csv("total_network.csv")
gll <- read.csv("../../综合的网络/与郭的网络进行比较/郭伟龙穗发育调控网络.csv")
head(net)
head(gll)
#郭伟龙网络中有多少在两个时期表达的基因
tpm <- read.csv("D:/AG_幼穗/RNA-seq/result/S1vsS2_sleuth_tpm_norm_gene_level.csv")
tpm <- tpm[,c(1,8,9)]
colnames(tpm) <- c("source_id","SM_tpm","FM_tpm")
head(tpm)
sgll <- gll[,c(1,3)]
colnames(sgll) <- c("source_id","v2")
tgll <- gll[,c(2,6)]
colnames(tgll) <- c("source_id","v2")
sgll <- merge(sgll,tpm,by="source_id",all=F)
tgll <- merge(tgll,tpm,by="source_id",all=F)
head(sgll)
head(tgll)
#绘制柱形图
sgll$type <- "source"
tgll$type <- "target"
ggll <- rbind(sgll,tgll)
ggll <- unique(ggll)
ggll <- ggll %>% 
  rownames_to_column(var = 'sample') %>% 
  pivot_longer( cols =  c("SM_tpm":"FM_tpm"),
                names_to = 'stage',
                values_to = 'expr')
head(ggll)
ggll$log2TPM <- log2(ggll$expr+1)
ggboxplot(ggll, x = "stage", y ="log2TPM", fill = "type",
          #add = "mean_sd", error.plot = "crossbar") +
          #add = "boxplot") +
          #outlier.shape = NA,
          add = "none") +
  scale_fill_manual(values=c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a","#d62728","#ff9896","#6A3D9A","#CAB2D6")) +
  ggtitle("S1") + 
  theme(plot.title = element_text(hjust = 0.5,size=18),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab(NULL) + ylab("log2 (TPM+1)")
###表达量在两个时期都小于0.5的基因数量，大于0.5的基因数目
head(ggll)
ggll1 <- ggll[which(ggll$expr <= 0.5),]
ggll1$group <- "<=0.5"
ggll2 <- ggll[which(ggll$expr > 0.5),]
ggll2$group <- ">0.5"
ggll <- rbind(ggll1,ggll2)
head(ggll)
ggll_number <- aggregate(ggll$source_id,list(ggll$type,ggll$stage,ggll$group),length)
head(ggll_number)
#柱形图
#pdf("表达的基因数目_在AI里换一下数值.pdf",family="ArialMT")
ggbarplot(ggll_number, x = "Group.1", y = "x", fill = "Group.3",
          color = "Group.3",
          palette = c("#aec7e8","#ffbb78","#98df8a","#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#CAB2D6","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928"),
          add = c("none")) +
  #scale_fill_manual(values=c("#aec7e8","#ffbb78","#98df8a","#ff9896","#c5b0d5","#c49c94","#f7b6d2","lightgrey")) +
  #scale_fill_manual(values=c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#CAB2D6","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928")) +
  ggtitle("Number of genes") + 
  theme(plot.title = element_text(hjust = 0.5,size=18),
        axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1)) +
  geom_text(aes(label = x), size = 3, vjust =-0.2, position = position_dodge(0.9)) +
  xlab(NULL) + ylab("Number")  +
  facet_grid(~Group.2)
#dev.off()


###网络中总覆盖的基因数目
library(VennDiagram)
head(net)
gglee <- ggll[which(ggll$expr > 0.5),]
head(gglee)
venn_list <- list(AG = unique(c(net$source,net$target)), GWL = unique(gglee$source_id))
venn.plot<-venn.diagram(venn_list, filename = NULL, 
                        fill = c('#7FC97F', '#FF7F00'), alpha = 0.50, 
                        cat.col = c('#7FC97F','#FF7F00'), cat.cex = 1.5, cat.fontfamily = 'serif',
                        col = c('#7FC97F', '#FF7F00'), cex = 1.5, fontfamily = 'serif')
pdf("与郭网络中总基因的overlap.pdf",family = "ArialMT")
grid.draw(venn.plot)
dev.off()

###两个网络与已知基因的overlap
gene <- read.xlsx("../../最终穗部性状基因R3.xlsx")
head(gene)
venn_list <- list(AG = unique(c(net$source,net$target)), GWL = unique(gglee$source_id), known = unique(gene$GeneID))
venn.plot<-venn.diagram(venn_list, filename = NULL, 
                        fill = c("#18499E","#6AB82D","#EFEA3C"), alpha = 0.50, 
                        cat.col = c("#18499E","#6AB82D","#EFEA3C"), cat.cex = 1.5, cat.fontfamily = 'serif',
                        col = c("#18499E","#6AB82D","#EFEA3C"), cex = 1.5, fontfamily = 'serif')
pdf("已知基因与郭网络中总基因的overlap.pdf",family = "ArialMT")
grid.draw(venn.plot)
dev.off()
