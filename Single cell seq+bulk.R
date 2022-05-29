library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(monocle)

logFCfilter=1               #logFC?Ĺ???????
adjPvalFilter=0.05          #????????pvalue?Ĺ???????
inputFile="treat.txt"       #?????ļ?
#setwd("###")        #???ù???Ŀ¼

#??ȡ?ļ????????????ļ?????????
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#??????ת??ΪSeurat???󣬲??????ݽ??й???
pbmc=CreateSeuratObject(counts = data,project = "seurat", min.cells=3, min.features=50, names.delim = "_")
#ʹ??PercentageFeatureSet???????????????????İٷֱ?
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
#???ƻ?????????С????ͼ
pdf(file="01.featureViolin.pdf", width=10, height=6)
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
pbmc=subset(x = pbmc, subset = nFeature_RNA > 50 & percent.mt < 5)    #?????ݽ??й???

#???????ȵ???????ͼ
pdf(file="01.featureCor.pdf",width=10,height=6)
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",,pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

#?????ݽ??б?׼??
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#??ȡϸ????????ϵ???ϴ??Ļ???
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)
#????????????ͼ
top10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf(file="01.featureVar.pdf",width=10,height=6)
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()





###################################02.PCA???ɷַ???###################################
##PCA????
pbmc=ScaleData(pbmc)          #PCA??ά֮ǰ?ı?׼Ԥ????????
pbmc=RunPCA(object= pbmc,npcs = 20,pc.genes=VariableFeatures(object = pbmc))     #PCA????

#????ÿ??PCA?ɷֵ?????????
pdf(file="02.pcaGene.pdf",width=10,height=8)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca",nfeatures = 20)
dev.off()

#???????ɷַ???ͼ??
pdf(file="02.PCA.pdf",width=6.5,height=6)
DimPlot(object = pbmc, reduction = "pca")
dev.off()

#???ɷַ?????ͼ
pdf(file="02.pcaHeatmap.pdf",width=10,height=8)
DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE,nfeatures = 30,ncol=2)
dev.off()

#?õ?ÿ??PC??pֵ?ֲ?
pbmc <- JackStraw(object = pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:15)
pdf(file="02.pcaJackStraw.pdf",width=8,height=6)
JackStrawPlot(object = pbmc, dims = 1:15)
dev.off()





###################################03.TSNE??????????marker????###################################
##TSNE????????
pcSelect=14
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)       #?????ڽӾ???
pbmc <- FindClusters(object = pbmc, resolution = 0.5)         #??ϸ??????,??ϸ????׼ģ?黯
pbmc <- RunTSNE(object = pbmc, dims = 1:pcSelect)             #TSNE????
pdf(file="03.TSNE.pdf",width=6.5,height=6)
TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)    #TSNE???ӻ?
dev.off()
write.table(pbmc$seurat_clusters,file="03.tsneCluster.txt",quote=F,sep="\t",col.names=F)

##????ÿ???????Ĳ???????
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="03.clusterMarkers.txt",sep="\t",row.names=F,quote=F)

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#????marker?ڸ???cluster????ͼ
pdf(file="03.tsneHeatmap.pdf",width=12,height=9)
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
dev.off()

#????marker??С????ͼ
pdf(file="03.markerViolin.pdf",width=10,height=6)
VlnPlot(object = pbmc, features = row.names(sig.markers)[1:2])
dev.off()

#??Ҫչʾ?Ļ??򣬿????޸?
showGenes=c("OLFM4","PDIA2","CPS1","LGALS2","TMED3","CLDN2","MSMB","AQP5","SPINK4","PIGR") 

#????marker?ڸ???cluster??ɢ??ͼ
pdf(file="03.markerScatter.pdf",width=10,height=6)
FeaturePlot(object = pbmc, features = showGenes, cols = c("green", "red"))
dev.off()

#????marker?ڸ???cluster??????ͼ
pdf(file="03.markerBubble.pdf",width=12,height=6)
cluster10Marker=showGenes
DotPlot(object = pbmc, features = cluster10Marker)
dev.off()




###################################04.SingleR R??ע??ϸ??????###################################
counts<-pbmc@assays$RNA@counts
clusters<-pbmc@meta.data$seurat_clusters
ann=pbmc@meta.data$orig.ident
#ref=get(load("ref_Human_all.RData"))
ref=celldex::HumanPrimaryCellAtlasData()
singler=SingleR(test=counts, ref =ref,
                labels=ref$label.main, clusters = clusters)
clusterAnn=as.data.frame(singler)
clusterAnn=cbind(id=row.names(clusterAnn), clusterAnn)
clusterAnn=clusterAnn[,c("id", "labels")]
write.table(clusterAnn,file="04.clusterAnn.txt",quote=F,sep="\t", row.names=F)
singler2=SingleR(test=counts, ref =ref, 
                labels=ref$label.main)
cellAnn=as.data.frame(singler2)
cellAnn=cbind(id=row.names(cellAnn), cellAnn)
cellAnn=cellAnn[,c("id", "labels")]
write.table(cellAnn, file="04.cellAnn.txt", quote=F, sep="\t", row.names=F)

#clusterע?ͺ??Ŀ??ӻ?
newLabels=singler$labels
names(newLabels)=levels(pbmc)
pbmc=RenameIdents(pbmc, newLabels)
pdf(file="04.TSNE.pdf",width=6.5,height=6)
TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)    #TSNE???ӻ?
dev.off()

##clusterע?ͺ??Ĳ???????
pbmc.markers=FindAllMarkers(object = pbmc,
                            only.pos = FALSE,
                            min.pct = 0.25,
                            logfc.threshold = logFCfilter)
sig.cellMarkers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.cellMarkers,file="04.cellMarkers.txt",sep="\t",row.names=F,quote=F)



###################################05.monocle R??ϸ???켣????###################################
#׼??ϸ???켣??????Ҫ???ļ?
monocle.matrix=as.matrix(pbmc@assays$RNA@data)
monocle.sample=pbmc@meta.data
monocle.geneAnn=data.frame(gene_short_name = row.names(monocle.matrix), row.names = row.names(monocle.matrix))
monocle.clusterAnn=clusterAnn
monocle.markers=sig.markers

#??Seurat????ת??Ϊmonocle??Ҫ??ϸ????????ϸ??ע?ͱ????ͻ???ע?ͱ???
data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"
pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])

#????ϸ??????????
clusterAnn=as.character(monocle.clusterAnn[,2])
names(clusterAnn)=paste0("cluster",monocle.clusterAnn[,1])
pData(cds)$cell_type2 <- plyr::revalue(as.character(pData(cds)$Cluster),clusterAnn)

#ϸ???켣????????
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- setOrderingFilter(cds, as.vector(sig.markers$gene))
#plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree')
cds <- orderCells(cds)
#??????֦??ϸ???켣ͼ
pdf(file="05.trajectory.State.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "State")
dev.off()
#????ʱ????ϸ???켣ͼ
pdf(file="05.trajectory.Pseudotime.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "Pseudotime")
dev.off()
#????ϸ?????Ƶ?ϸ???켣ͼ
pdf(file="05.trajectory.cellType.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "cell_type2")
dev.off()
#??????????ϸ???켣ͼ
pdf(file="05.trajectory.cluster.pdf",width=6.5,height=6)
plot_cell_trajectory(cds, color_by = "Cluster")
dev.off()

#ϸ???켣????????
groups=subset(pData(cds),select='State')
pbmc=AddMetaData(object=pbmc, metadata=groups, col.name="group")
geneList=list()
for(i in levels(factor(groups$State))){
	pbmc.markers=FindMarkers(pbmc, ident.1 = i, group.by = 'group')
	sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
	sig.markers=cbind(Gene=row.names(sig.markers), sig.markers)
	write.table(sig.markers,file=paste0("05.monocleDiff.", i, ".txt"),sep="\t",row.names=F,quote=F)
	geneList[[i]]=row.names(sig.markers)
}
#???潻??????
unionGenes=Reduce(union,geneList)
write.table(file="05.monocleDiff.union.txt",unionGenes,sep="\t",quote=F,col.names=F,row.names=F)


library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(monocle)

#读取对照组的数据,并对数据进行整理
rt=read.table("control.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
conData=avereps(data)
colnames(conData)=paste0("C.", colnames(conData))

#读取实验组的数据,并对数据进行整理
rt=read.table("treat.txt", header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
treatData=avereps(data)
colnames(treatData)=paste0("T.", colnames(treatData))

#数据合并
sameGene=intersect(row.names(conData), row.names(treatData))
data=cbind(conData[sameGene,], treatData[sameGene,])

#将矩阵转换为Seurat对象，并对数据进行过滤
pbmc <- CreateSeuratObject(counts = data,project = "seurat", min.cells = 3, min.features = 50, names.delim = "_")
#使用PercentageFeatureSet函数计算线粒体基因的百分比
pbmc[["percent.mt"]]=PercentageFeatureSet(object = pbmc, pattern = "^MT-")
pbmcCon=subset(x = pbmc, subset = nFeature_RNA > 50 & percent.mt < 5)    #对数据进行过滤

#对数据进行标准化
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#提取那些在细胞间变异系数较大的基因
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)

#组间差异分析
logFCfilter=1
adjPvalFilter=0.05
groups=gsub("(.*?)\\..*", "\\1", colnames(pbmc))
names(groups)=colnames(pbmc)
pbmc=AddMetaData(object=pbmc, metadata=groups, col.name="group")
pbmc.markers=FindMarkers(pbmc, ident.1 = "T", ident.2 = "C", group.by = 'group')
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
sig.markers=cbind(Gene=row.names(sig.markers), sig.markers)
write.table(sig.markers,file="diffGene.txt",sep="\t",row.names=F,quote=F)

library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05       #p值过滤条件
qvalueFilter=0.05       #矫正后的p值过滤条件
inputFile="05.monocleDiff.1.txt"      #输入文件名称

#读取输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F)

#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

#基因名字转换为基因id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#GO富集分析
kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
#保存富集结果
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

#定义显示GO的数目
showNum=10
if(nrow(GO)<30){
	showNum=nrow(GO)
}

#柱状图
pdf(file="barplot.pdf", width=9, height=7)
bar=barplot(kk, drop=TRUE, showCategory=showNum, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()
		
#气泡图
pdf(file="bubble.pdf", width=9, height=7)
bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bub)
dev.off()

library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05         #p值过滤条件
qvalueFilter=0.05         #矫正后的p值过滤条件
inputFile="05.monocleDiff.1.txt"      #输入文件名称

#读取输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F)

#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

#基因名字转换为基因id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=data.frame(genes, entrezID=entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]        #去除基因id为NA的基因
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

#kegg富集分析
kk=enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$genes[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
#保存富集结果
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#定义显示通路的数目
showNum=30
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

#柱状图
pdf(file="barplot.pdf", width=9, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, color=colorSel)
dev.off()

#气泡图
pdf(file="bubble.pdf", width = 9, height = 7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", color=colorSel)
dev.off()


library(limma)
library(ConsensusClusterPlus)
expFile="geoMatrix.txt"                 #表达数据文件
geneFile="05.monocleDiff.union.txt"     #基因列表文件

#读取输入文件，并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(data)

#如果数据没有取log2,会自动对数据取log2
qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
    rt[rt<0]=0
    rt=log2(rt+1)}

#对数据进行矫正
data=normalizeBetweenArrays(rt)
outData=rbind(ID=colnames(data), data)
write.table(outData, file="geoNormalize.txt", sep="\t", quote=F, col.names=F)

#获取基因的表达量
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(data))
data=data[sameGene,]

#聚类
maxK=9
results=ConsensusClusterPlus(data,
              maxK=maxK,
              reps=50,
              pItem=0.8,
              pFeature=1,
              title=workDir,
              clusterAlg="pam",
              distance="euclidean",
              seed=123456,
              plot="png")


#输出分型结果
clusterNum=3        #分几类，根据判断标准判断
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
cluster$cluster=paste0("C", cluster$cluster)
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="cluster.txt", sep="\t", quote=F, col.names=F)

library(survival)
library(survminer)
clusterFile="cluster.txt"       #分型结果文件
cliFile="time.txt"              #生存数据文件

#读取输入文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
#rownames(cluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(cluster))
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
colnames(cli)=c("futime", "fustat")
cli$futime=cli$futime/365

#数据合并
sameSample=intersect(row.names(cluster), row.names(cli))
rt=cbind(cli[sameSample,,drop=F], cluster[sameSample,,drop=F])

#生存差异统计
length=length(levels(factor(rt$cluster)))
diff=survdiff(Surv(futime, fustat) ~ cluster, data = rt)
pValue=1-pchisq(diff$chisq, df=length-1)
if(pValue<0.001){
	pValue="p<0.001"
}else{
	pValue=paste0("p=",sprintf("%.03f",pValue))
}
fit <- survfit(Surv(futime, fustat) ~ cluster, data = rt)
#print(surv_median(fit))

#绘制生存曲线
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length]
surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=F,
		           pval=pValue,
		           pval.size=6,
		           legend.title="Cluster",
		           legend.labs=levels(factor(rt[,"cluster"])),
		           legend = c(0.8, 0.8),
		           font.legend=10,
		           xlab="Time(years)",
		           break.time.by = 1,
		           palette = bioCol,
		           surv.median.line = "hv",
		           risk.table=T,
		           cumevents=F,
		           risk.table.height=.25)
pdf(file="survival.pdf",onefile = FALSE,width=7,height=5.5)
print(surPlot)
dev.off()

library(ggplot2)


#读取输入文件
clinical=read.table("clinical.txt", header=T, sep="\t", check.names=F, row.names=1)
cluster=read.table("cluster.txt", header=T, sep="\t", check.names=F, row.names=1)

#合并数据
sameSample=intersect(row.names(clinical), row.names(cluster))
clinical=clinical[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
data=cbind(clinical, cluster)

#绘制柱状图
for(i in colnames(data)[1:(ncol(data)-1)]){
    rt=data[,c("cluster",i)]
    colnames(rt)=c("Cluster","type")
    #统计检验
    tableStat=table(rt)
	stat=chisq.test(tableStat)
	pvalue=stat$p.value
	if(pvalue<0.001){
		pvalue="p<0.001"
	}else{
		pvalue=paste0("p=",sprintf("%.03f",pvalue))
	}
	#绘制柱状图
	p=ggplot(rt, aes(Cluster)) + 
	    geom_bar(aes(fill=type), position="fill")+
	    labs(x = 'Cluster',y = '',title=paste0(i," (",pvalue,")"))+
	    guides(fill = guide_legend(title =i))+
	    theme_bw()+
	    theme(plot.title = element_text(hjust = 0.5))
    pdf(file=paste0("cliCor.",i,".pdf"),width=4.5,height=6)	
	print(p)
	dev.off()
}


library(limma)
library(reshape2)
library(ggpubr)
expFile="geoNormalize.txt"            #表达数据文件
geneFile="05.monocleDiff.2.txt"       #基因列表文件
cluFile="cluster.txt"                 #分型结果文件

#读取输入文件，并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#获取基因的表达量
geneRT=read.table(geneFile, header=T, sep="\t", check.names=F)
up=geneRT[geneRT[,3]>0,]
down=geneRT[geneRT[,3]<0,]
sameGene=intersect(as.vector(geneRT[,1]), rownames(data))
data=t(data[sameGene,])

#读取分型的结果文件
cluster=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(cluster))
rt1=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])

#绘制箱线图
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
data=melt(rt1, id.vars=c("cluster"))
colnames(data)=c("Cluster","Gene","Expression")
data$Cluster=factor(data$Cluster, levels=levels(factor(data$Cluster)))
yIndex=gsub(".+\\.(.+?)\\.txt", "\\1", geneFile)

#绘制上调基因的箱线图
data1=data[which(data[,2] %in% as.vector(up[,1])),]
p=ggboxplot(data1, x="Cluster", y="Expression", fill = "Cluster",
	     ylab=paste0("Up-regulated gene in Branch ",yIndex), add = "none", xlab="Cluster", palette =bioCol[1:length(levels(factor(data$Cluster)))])
#输出图片文件
pdf(file=paste0(yIndex,".up.pdf"), width=4.5, height=6)
print(p)
dev.off()

#绘制下调基因的箱线图
data2=data[which(data[,2] %in% as.vector(down[,1])),]
p=ggboxplot(data2, x="Cluster", y="Expression", fill = "Cluster",
	     ylab=paste0("Down-regulated gene in Branch ",yIndex), add = "none", xlab="Cluster", palette =bioCol[1:length(levels(factor(data$Cluster)))])
#输出图片文件
pdf(file=paste0(yIndex,".down.pdf"), width=4.5, height=6)
print(p)
dev.off()

library(limma)
library(estimate)
expFile="geoNormalize.txt"      #表达数据文件

#读取输入文件，并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#输出整理后的矩阵文件
out=rbind(ID=colnames(data), data)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)

#运行estimate包
filterCommonGenes(input.f="uniq.symbol.txt", 
                  output.f="commonGenes.gct", 
                  id="GeneSymbol")

estimateScore(input.ds = "commonGenes.gct",
              output.ds="estimateScore.gct")

#输出每个样品的打分
scores=read.table("estimateScore.gct", skip=2, header=T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
rownames(scores)=gsub("\\.", "\\-", rownames(scores))
out=rbind(ID=colnames(scores), scores)
write.table(out, file="scores.txt", sep="\t", quote=F, col.names=F)

library(ggpubr)      #引用包
cluFile="cluster.txt"       #分型结果文件
socreFile="scores.txt"      #肿瘤微环境打分文件

#读取分型结果数据
cluster=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)

#读取肿瘤微环境打分文件
score=read.table("scores.txt", header=T, sep="\t", check.names=F, row.names=1)

#数据合并
sameSample=intersect(row.names(cluster), row.names(score))
data=cbind(score[sameSample,,drop=F], cluster[sameSample,,drop=F])

#设置比较组
group=levels(factor(data$cluster))
data$cluster=factor(data$cluster, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#绘制小提琴图
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
for(i in colnames(data)[1:(ncol(data)-1)]){
	violin=ggviolin(data, x="cluster", y=i, fill = "cluster",
	         xlab="Cluster", ylab=i,
	         legend.title="Cluster",
	         palette=bioCol, 
	         add="boxplot", add.params = list(fill="white"))+ 
		stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
	pdf(file=paste0("vioplot.", i, ".pdf"), width=6, height=5.5)
	print(violin)
	dev.off()
}

library("limma")         #引用包
expFile="geoNormalize.txt"       #表达数据文件

#读取输入文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#输出整理后的矩阵文件
out=rbind(ID=colnames(data), data)
write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)

#运行CIBERSORT，得到免疫细胞浸润的结果
source("sCell22.CIBERSORT.R")
results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=1000, QN=TRUE)

library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
immFile="CIBERSORT-Results.txt"     #免疫细胞浸润结果文件
cluFile="cluster.txt"               #分型结果文件
pFilter=0.05            #免疫细胞浸润结果的过滤条件
setwd("C:\\Users\\lexb4\\Desktop\\sCell\\23.immunePlot")     #设置工作目录

#读取免疫细胞结果文件，并对数据进行整理
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<pFilter,]
data=as.matrix(immune[,1:(ncol(immune)-3)])

#读取分型文件
cluster=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并
sameSample=intersect(row.names(data), row.names(cluster))
data=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])
data=data[order(data$cluster),]
gaps=c(1, as.vector(cumsum(table(data$cluster))))
xlabels=levels(factor(data$cluster))


##################绘制柱状图##################
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
data1=t(as.matrix(data[,-ncol(data)]))
pdf("barplot.pdf",height=10,width=18)
col=rainbow(nrow(data1),s=0.7,v=0.7)
par(las=1,mar=c(8,5,4,16),mgp=c(3,0.1,0),cex.axis=1.5)
a1=barplot(data1,col=col,yaxt="n",ylab="Relative Percent",xaxt="n",cex.lab=1.8)
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
par(srt=0,xpd=T)
for(i in 1:length(gaps)){
	j=i+1
	rect(xleft=a1[gaps[i]], ybottom = -0.01, xright = a1[gaps[j]], ytop= -0.06, col=bioCol[i])
	text((a1[gaps[i]]+a1[gaps[j]])/2,-0.035,xlabels[i],cex=2)
}
ytick2 = cumsum(data1[,ncol(data1)])
ytick1 = c(0,ytick2[-length(ytick2)])
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data1),col=col,pch=15,bty="n",cex=1.3)
dev.off()

##################绘制箱线图##################
#把数据转换成ggplot2输入文件
data=melt(data,id.vars=c("cluster"))
colnames(data)=c("cluster", "Immune", "Expression")

#绘制箱线图
group=levels(factor(data$cluster))
data$cluster=factor(data$cluster, levels=group)
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="cluster",
				  xlab="",
				  ylab="Fraction",
				  legend.title="Cluster",
				  width=0.8,
				  palette=bioCol)+
				  rotate_x_text(50)+
	stat_compare_means(aes(group=cluster),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")

#输出图片
pdf(file="immune.diff.pdf", width=7, height=6)
print(boxplot)
dev.off()

library(limma)
library(survival)
library(survminer)
immuneFile="CIBERSORT-Results.txt"     #免疫浸润结果文件
surFile="time.txt"                     #生存数据文件
pFilter=0.05                           #CIBERSORT结果过滤条件
setwd("C:\\Users\\lexb4\\Desktop\\sCell\\24.immuneSur")      #设置工作目录

#读取免疫浸润结果文件，并对数据进行整理
immune=read.table(immuneFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<pFilter,]
data=as.matrix(immune[,1:(ncol(immune)-3)])

#读取生存数据
surTime=read.table(surFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并
sameSample=intersect(row.names(data),row.names(surTime))
data=data[sameSample,]
surTime=surTime[sameSample,]
rt=cbind(surTime, data)
rt$futime=rt$futime/365

outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
	if(sd(rt[,i])<0.01){next}
	#KM分析
	group=ifelse(rt[,i]>median(rt[,i]), "high", "low")
	diff=survdiff(Surv(futime, fustat) ~group, data = rt)
	pValue=1-pchisq(diff$chisq, df=1)

	#对p<0.05的免疫细胞绘制生存曲线
	if(pValue<pFilter){
		outVector=cbind(i, pValue)
		outTab=rbind(outTab, outVector)
		#绘制生存曲线	    
		if(pValue<0.001){
			pValue="p<0.001"
		}else{
			pValue=paste0("p=",sprintf("%.03f",pValue))
		}
		fit <- survfit(Surv(futime, fustat) ~ group, data = rt)
		surPlot=ggsurvplot(fit, 
				           data=rt,
				           pval=pValue,
				           pval.size=6,
				           conf.int=T,
				           legend.title=i,
				           legend.labs=c("High", "Low"),
				           palette=c("red", "blue"),
				           xlab="Time(years)",
				           break.time.by=1,
				           risk.table=F,
				           risk.table.title="",
				           risk.table.height=.35)
		pdf(file=paste0("sur.",i,".pdf"), onefile = FALSE, width=6, height=5)
		print(surPlot)
		dev.off()
	}
}
#输出免疫细胞和p值表格文件
write.table(outTab,file="immuneSur.result.txt",sep="\t",row.names=F,quote=F)


library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
expFile="geoNormalize.txt"      #表达数据文件
cluFile="cluster.txt"           #分型结果文件
geneFile="gene.txt"             #免疫检查点的基因文件


#读取基因表达文件,并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
	
#读取基因文件
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data),as.vector(gene[,1]))
data=t(data[sameGene,])

#合并分型数据
cluster=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(cluster))
rt1=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])

#把数据转换成ggplot2输入文件
data=melt(rt1,id.vars=c("cluster"))
colnames(data)=c("cluster", "Gene", "Expression")

#绘制箱线图
group=levels(factor(data$cluster))
data$cluster=factor(data$cluster, levels=group)
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Gene", y="Expression", fill="cluster",
				  orientation="horizontal",
				  xlab="",
				  ylab="Gene expression",
				  legend.title="Cluster",
				  width=1,
				  palette=bioCol)+
				  rotate_x_text(50)+
	stat_compare_means(aes(group=cluster),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")

#输出图片
pdf(file="checkpoint.diff.pdf", width=7, height=8)
print(boxplot)
dev.off()


library(limma)
library(survival)
library(survminer)
expFile="geoNormalize.txt"      #表达数据文件
geneFile="gene.txt"             #免疫检查点的基因文件
surFile="time.txt"              #生存数据文件

#读取基因表达文件,并对数据进行处理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
	
#读取基因文件
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(row.names(data),as.vector(gene[,1]))
data=t(data[sameGene,])

#读取生存数据
surTime=read.table(surFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并
sameSample=intersect(row.names(data),row.names(surTime))
data=data[sameSample,]
surTime=surTime[sameSample,]
rt=cbind(surTime, data)
rt$futime=rt$futime/365

outTab=data.frame()
for(i in colnames(rt[,3:ncol(rt)])){
	if(sd(rt[,i])<0.1){next}
	#KM分析
	group=ifelse(rt[,i]>median(rt[,i]), "high", "low")
	diff=survdiff(Surv(futime, fustat) ~group, data = rt)
	pValue=1-pchisq(diff$chisq, df=1)

	#对p<0.05的免疫细胞绘制生存曲线
	if(pValue<0.05){
		outVector=cbind(i, pValue)
		outTab=rbind(outTab, outVector)
		#绘制生存曲线	    
		if(pValue<0.001){
			pValue="p<0.001"
		}else{
			pValue=paste0("p=",sprintf("%.03f",pValue))
		}
		fit <- survfit(Surv(futime, fustat) ~ group, data = rt)
		surPlot=ggsurvplot(fit, 
				           data=rt,
				           pval=pValue,
				           pval.size=6,
				           conf.int=T,
				           legend.title=i,
				           legend.labs=c("High", "Low"),
				           palette=c("red", "blue"),
				           xlab="Time(years)",
				           break.time.by=1,
				           risk.table=F,
				           risk.table.title="",
				           risk.table.height=.35)
		pdf(file=paste0("sur.",i,".pdf"), onefile = FALSE, width=6, height=5)
		print(surPlot)
		dev.off()
	}
}
#输出基因和p值表格文件
write.table(outTab,file="geneSur.result.txt",sep="\t",row.names=F,quote=F)

library(limma)
library(sva)
tcgaExpFile="symbol.txt"                 #TCGA表达数据文件
geoExpFile="geoNormalize.txt"            #GEO表达数据文件
geneFile="05.monocleDiff.union.txt"      #基因列表文件

#读取TCGA基因表达文件,并对数据进行处理
rt=read.table(tcgaExpFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
tcga=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
tcga=avereps(tcga)

#FPKM转换为TPM
fpkmToTpm=function(fpkm){
    exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
tpm=apply(tcga, 2, fpkmToTpm)
tcga=log2(tpm+1)

#读取geo基因表达文件,并对数据进行处理
rt = read.table(geoExpFile,header=T,sep="\t",check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
geo=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
geo=avereps(geo)

#对基因取交集,分别得到交集基因在TCGA矩阵和GEO矩阵的表达量
sameGene=intersect(row.names(tcga),row.names(geo))
tcgaOut=tcga[sameGene,]
geoOut=geo[sameGene,]

#批次矫正
all=cbind(tcgaOut,geoOut)
batchType=c(rep(1,ncol(tcgaOut)),rep(2,ncol(geoOut)))
outTab=ComBat(all, batchType, par.prior=TRUE)
tcgaOut=outTab[,colnames(tcgaOut)]
tcgaOut[tcgaOut<0]=0
geoOut=outTab[,colnames(geoOut)]
geoOut[geoOut<0]=0

#获取基因的表达量
gene=read.table(geneFile, header=F, sep="\t", check.names=F)
sameGene=intersect(as.vector(gene[,1]), rownames(tcgaOut))
tcgaShareExp=tcgaOut[sameGene,]
geoShareExp=geoOut[sameGene,]

#输出细胞轨迹基因的表达量
tcgaShareExp=rbind(ID=colnames(tcgaShareExp),tcgaShareExp)
write.table(tcgaShareExp,file="tcga.share.txt",sep="\t",quote=F,col.names=F)
geoShareExp=rbind(ID=colnames(geoShareExp),geoShareExp)
write.table(geoShareExp,file="geo.share.txt",sep="\t",quote=F,col.names=F)

#引用包
library(limma)
library(WGCNA)
expFile="tcga.share.txt"      #表达数据文件
cliFile="clinical.txt"        #临床数据文件
setwd("C:\\Users\\lexb4\\Desktop\\sCell\\28.WGCNA")      #设置工作目录

#读取输入文件，并对输入文件进行整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[apply(data,1,sd)>0,]     #删除波动小的基因

#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
data=data[,group==0,drop=F]
data=t(data)
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
datExpr0=avereps(data)

###检查缺失值
gsg = goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK){
	# Optionally, print the gene and sample names that were removed:
	if (sum(!gsg$goodGenes)>0)
	    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
	if (sum(!gsg$goodSamples)>0)
	    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
	# Remove the offending genes and samples from the data:
	datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

###样品聚类
sampleTree = hclust(dist(datExpr0), method = "average")
pdf(file = "1_sample_cluster.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
###剪切线
abline(h = 10000, col = "red")
dev.off()

###删除剪切线以下的样品
clust = cutreeStatic(sampleTree, cutHeight = 10000, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr0 = datExpr0[keepSamples, ]

###power值散点图
enableWGCNAThreads()   #多线程工作
powers = c(1:20)       #幂指数范围1:20
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
pdf(file="2_scale_independence.pdf",width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9
###拟合指数与power值散点图
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red") #可以修改
###平均连通性与power值散点图
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


###邻接矩阵转换
sft #查看最佳power值
softPower =sft$powerEstimate   #最佳power值
adjacency = adjacency(datExpr0, power = softPower)
softPower

###TOM矩阵
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

###基因聚类
geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file="3_gene_clustering.pdf",width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()


###动态剪切模块识别
minModuleSize=30      #模块基因数目,每个模块最少包含多少个基因
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
pdf(file="4_Dynamic_Tree.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()


###查找相似模块聚类
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")
pdf(file="5_Clustering_module.pdf",width=7,height=7)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.3  #剪切高度可修改
abline(h=MEDissThres, col = "red")
dev.off()


###相似模块合并
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
pdf(file="6_merged_dynamic.pdf", width = 9, height = 6)
plotDendroAndColors(geneTree, mergedColors,"Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()
moduleColors = mergedColors
table(moduleColors)
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs


###模块与临床性状的相关型热图
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(cli), row.names(MEs))
MEs=MEs[sameSample,]
datTraits=cli[sameSample,]
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
pdf(file="7_Module_trait.pdf", width=6.5, height=5.2)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

###基因所在的模块
probes = colnames(datExpr0)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "allModules.txt",sep="\t",row.names=F,quote=F)


###输出每个模块的基因
for (mod in 1:nrow(table(moduleColors))){  
	modules = names(table(moduleColors))[mod]
	probes = colnames(datExpr0)
	inModule = (moduleColors == modules)
	modGenes = probes[inModule]
	write.table(modGenes, file =paste0("module_",modules,".txt"),sep="\t",row.names=F,col.names=F,quote=F)
}

#引用包
library(limma)
library(ggplot2)
library(pheatmap)
expFile="tcga.share.txt"       #表达数据文件
logFCfilter=1                  #logFC临界值
fdrFilter=0.05                 #fdr临界值
setwd("C:\\Users\\lexb4\\Desktop\\sCell\\29.diff")   #设置工作目录

#读取输入文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#读取模块中的基因
geneVec=c()
geneFiles=dir()                          #获取目录下所有文件
geneFiles=grep("^module_", geneFiles, value=T)     #提取module_开头的文件
for(geneFile in geneFiles){
	gene=read.table(geneFile, header=F, sep="\t", check.names=F)
	geneVec=c(geneVec, as.vector(gene[,1]))
}
sameGene=intersect(geneVec, row.names(data))
data=data[sameGene,]

#正常和肿瘤数目
group=sapply(strsplit(colnames(data),"\\-"), "[", 4)
group=sapply(strsplit(group,""), "[", 1)
group=gsub("2", "1", group)
conNum=length(group[group==1])       #正常组样品数目
treatNum=length(group[group==0])     #肿瘤组样品数目
Type=c(rep(1,conNum), rep(2,treatNum))

#差异分析
outTab=data.frame()
for(i in row.names(data)){
	rt=data.frame(expression=data[i,], Type=Type)
	wilcoxTest=wilcox.test(expression ~ Type, data=rt)
	pvalue=wilcoxTest$p.value
	conGeneMeans=mean(data[i,1:conNum])
	treatGeneMeans=mean(data[i,(conNum+1):ncol(data)])
	logFC=treatGeneMeans-conGeneMeans
	conMed=median(data[i,1:conNum])
	treatMed=median(data[i,(conNum+1):ncol(data)])
	diffMed=treatMed-conMed
	if( ((logFC>0) & (diffMed>0)) | ((logFC<0) & (diffMed<0)) ){  
		outTab=rbind(outTab,cbind(gene=i,conMean=conGeneMeans,treatMean=treatGeneMeans,logFC=logFC,pValue=pvalue))
	}
}
pValue=outTab[,"pValue"]
fdr=p.adjust(as.numeric(as.vector(pValue)), method="fdr")
outTab=cbind(outTab, fdr=fdr)

#输出所有基因的差异情况
write.table(outTab,file="tcga.all.txt",sep="\t",row.names=F,quote=F)

#输出差异表格
outDiff=outTab[( abs(as.numeric(as.vector(outTab$logFC)))>logFCfilter & as.numeric(as.vector(outTab$fdr))<fdrFilter),]
write.table(outDiff,file="tcga.diff.txt",sep="\t",row.names=F,quote=F)

#输出差异基因的表达文件
diffExp=rbind(ID=colnames(data[as.vector(outDiff[,1]),]),data[as.vector(outDiff[,1]),])
write.table(diffExp,file="tcga.diffExp.txt",sep="\t",col.names=F,quote=F)

#绘制差异基因热图
hmExp=data[as.vector(outDiff[,1]),]
Type=c(rep("N",conNum),rep("T",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="heatmap.pdf",width=10,height=8)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c(rep("blue",3), "white", rep("red",3)))(50),
         cluster_cols =F,
         scale="row",
         show_colnames = F,
         show_rownames = T,
         fontsize = 8,
         fontsize_row=3,
         fontsize_col=8)
dev.off()

#定义显著性
outTab$fdr=as.numeric(outTab$fdr)
outTab$logFC=as.numeric(outTab$logFC)
Significant=ifelse((outTab$fdr<fdrFilter & abs(outTab$logFC)>logFCfilter), ifelse(outTab$logFC>logFCfilter,"Up","Down"), "Not")
#绘制火山图
p = ggplot(outTab, aes(logFC, -log10(fdr)))+
    geom_point(aes(col=Significant))+
    scale_color_manual(values=c("green", "black", "red"))+
    labs(title = " ")+
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
p=p+theme_bw()
#保存为图片
pdf("vol.pdf", width=6.2, height=5.5)
print(p)
dev.off()


#引用包
library(limma)
library(survival)
expFile="tcga.diffExp.txt"     #差异基因表达文件
cliFile="time.txt"             #生存数据文件
setwd("C:\\Users\\lexb4\\Desktop\\sCell\\30.uniCox")     #设置工作目录

#读取表达文件，并对输入文件整理
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]

#删掉正常样品
group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2","1",group)
data=t(data[,group==0])
rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data=avereps(data)

#读取生存数据
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)     #读取临床文件
cli$futime=cli$futime/365

#数据合并
sameSample=intersect(row.names(data), row.names(cli))
data=data[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
rt=cbind(cli, data)

#对基因进行循环，找出预后相关的基因
outTab=data.frame()
sigGenes=c("futime","fustat")
for(i in colnames(rt[,3:ncol(rt)])){
	#cox分析
	cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
	coxSummary = summary(cox)
	coxP=coxSummary$coefficients[,"Pr(>|z|)"]
	if(coxP<0.05){
		sigGenes=c(sigGenes,i)
		outTab=rbind(outTab,
				         cbind(id=i,
				         HR=coxSummary$conf.int[,"exp(coef)"],
				         HR.95L=coxSummary$conf.int[,"lower .95"],
				         HR.95H=coxSummary$conf.int[,"upper .95"],
				         pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
				        )
	}
}

#输出单因素的结果
write.table(outTab,file="uniCox.txt",sep="\t",row.names=F,quote=F)

#保存单因素显著基因的表达量
sigGeneExp=rt[,sigGenes]
sigGeneExp=rbind(id=colnames(sigGeneExp), sigGeneExp)
write.table(sigGeneExp, file="tcga.uniSigExp.txt", sep="\t", quote=F, col.names=F)

############绘制森林图函数############
bioForest=function(coxFile=null, forestFile=null, forestCol=null){
	#读取输入文件
	rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
	gene <- rownames(rt)
	hr <- sprintf("%.3f",rt$"HR")
	hrLow  <- sprintf("%.3f",rt$"HR.95L")
	hrHigh <- sprintf("%.3f",rt$"HR.95H")
	Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
	pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
		
	#输出图形
	pdf(file=forestFile, width=6.5, height=5)
	n <- nrow(rt)
	nRow <- n+1
	ylim <- c(1,nRow)
	layout(matrix(c(1,2),nc=2),width=c(3,2.5))
		
	#绘制森林图左边的基因信息
	xlim = c(0,3)
	par(mar=c(4,2.5,2,1))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
	text.cex=0.8
	text(0,n:1,gene,adj=0,cex=text.cex)
	text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
	text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)
	
	#绘制森林图
	par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
	xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
	arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
	abline(v=1,col="black",lty=2,lwd=2)
	boxcolor = ifelse(as.numeric(hr) > 1, forestCol[1], forestCol[2])
	points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.6)
	axis(1)
	dev.off()
}

bioForest(coxFile="uniCox.txt", forestFile="forest.pdf", forestCol=c("red","green"))




library(limma)                #引用包
expFile="geo.share.txt"       #表达数据文件
cliFile="time.txt"            #临床数据
setwd("C:\\Users\\lexb4\\Desktop\\sCell\\31.geoMergeTime")    #设置工作目录

#读取表达文件，并对输入文件整理
rt=read.table(expFile,header=T,sep="\t",check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data=t(data)

#读取生存数据
cli=read.table(cliFile,header=T,sep="\t",check.names=F,row.names=1)

#数据合并并输出结果
sameSample=intersect(row.names(data),row.names(cli))
data=data[sameSample,]
cli=cli[sameSample,]
out=cbind(cli,data)
out=cbind(id=row.names(out),out)
write.table(out,file="geo.expTime.txt",sep="\t",row.names=F,quote=F)

#引用包
library(glmnet)
library(survival)
trainFile="tcga.uniSigExp.txt"      #train组输入文件
testFile="geo.expTime.txt"          #test组输入文件
setwd("C:\\Users\\lexb4\\Desktop\\sCell\\32.model")                        #设置工作目录
rt=read.table(trainFile, header=T, sep="\t", row.names=1,check.names=F)    #读取train组输入文件

#COX模型构建
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

#输出模型相关信息
outMultiTab=data.frame()
outMultiTab=cbind(
		          coef=multiCoxSum$coefficients[,"coef"],
		          HR=multiCoxSum$conf.int[,"exp(coef)"],
		          HR.95L=multiCoxSum$conf.int[,"lower .95"],
		          HR.95H=multiCoxSum$conf.int[,"upper .95"],
		          pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
write.table(outMultiTab,file="multiCox.txt",sep="\t",row.names=F,quote=F)

#输出train组风险值
trainScore=predict(multiCox,type="risk",newdata=rt)
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="trainRisk.txt",sep="\t",quote=F,row.names=F)

#输出test组风险值
rt=read.table(testFile, header=T, sep="\t", row.names=1, check.names=F)
rt$futime=rt$futime/365
testFinalGeneExp=rt[,coxGene]
testScore=predict(multiCox,type="risk",newdata=rt)
outCol=c("futime","fustat",coxGene)
risk=as.vector(ifelse(testScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="testRisk.txt",sep="\t",quote=F,row.names=F)

library(survival)
library(survminer)
setwd("C:\\Users\\lexb4\\Desktop\\sCell\\33.survival")      #设置工作目录

bioSurvival=function(inputFile=null, outFile=null){
	#读取输入文件
	rt=read.table(inputFile, header=T, sep="\t", check.names=F)
	#比较高低风险组生存差异，得到显著性p值
	diff=survdiff(Surv(futime, fustat) ~risk,data = rt)
	pValue=1-pchisq(diff$chisq,df=1)
	if(pValue<0.001){
		pValue="p<0.001"
	}else{
		pValue=paste0("p=",sprintf("%.03f",pValue))
	}
	fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
		
	#绘制生存曲线
	surPlot=ggsurvplot(fit, 
		           data=rt,
		           conf.int=TRUE,
		           pval=pValue,
		           pval.size=5,
		           legend.title="Risk",
		           legend.labs=c("High risk", "Low risk"),
		           xlab="Time(years)",
		           break.time.by = 1,
		           palette=c("red", "blue"),
		           risk.table=F,
		           risk.table.title="",
		           risk.table.height=.25)
	pdf(file=outFile,onefile = FALSE,width = 5,height =4.5)
	print(surPlot)
	dev.off()
}

#调用函数，绘制生存曲线
bioSurvival(inputFile="trainRisk.txt", outFile="train.survival.pdf")
bioSurvival(inputFile="testRisk.txt", outFile="test.survival.pdf")

#引用包
library(survival)
library(survminer)
library(timeROC)
setwd("C:\\Users\\lexb4\\Desktop\\sCell\\34.ROC")      #设置工作目录

#定义绘制ROC曲线函数
bioROC=function(inputFile=null, rocFile=null){
	#读取输入文件
	rt=read.table(inputFile, header=T, sep="\t", check.names=F)
	#ROC曲线
	ROC_rt=timeROC(T=rt$futime,delta=rt$fustat,
	               marker=rt$riskScore,cause=1,
	               weighting='aalen',
	               times=c(1,3,5),ROC=TRUE)
	pdf(file=rocFile,width=5,height=5)
	plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
	plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
	plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
	legend('bottomright',
	        c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
	          paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
	          paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
	        col=c("green",'blue','red'),lwd=2,bty = 'n')
	dev.off()
}

#调用函数,绘制ROC曲线
bioROC(inputFile="trainRisk.txt", rocFile="train.ROC.pdf")
bioROC(inputFile="testRisk.txt", rocFile="test.ROC.pdf")

library(survival)       #引用包
setwd("C:\\Users\\lexb4\\Desktop\\sCell\\35.indep")     #设置工作目录

############绘制森林图函数############
bioForest=function(coxFile=null, forestFile=null, forestCol=null){
	#读取输入文件
	rt <- read.table(coxFile, header=T, sep="\t", check.names=F, row.names=1)
	gene <- rownames(rt)
	hr <- sprintf("%.3f",rt$"HR")
	hrLow  <- sprintf("%.3f",rt$"HR.95L")
	hrHigh <- sprintf("%.3f",rt$"HR.95H")
	Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
	pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
		
	#输出图形
	pdf(file=forestFile, width=6.5, height=4.5)
	n <- nrow(rt)
	nRow <- n+1
	ylim <- c(1,nRow)
	layout(matrix(c(1,2),nc=2),width=c(3,2.5))
		
	#绘制森林图左边的临床信息
	xlim = c(0,3)
	par(mar=c(4,2.5,2,1))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
	text.cex=0.8
	text(0,n:1,gene,adj=0,cex=text.cex)
	text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
	text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
		
	#绘制右边的森林图
	par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
	xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
	arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=3)
	abline(v=1, col="black", lty=2, lwd=2)
	boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
	points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=2)
	axis(1)
	dev.off()
}
############绘制森林图函数############

#定义独立预后分析函数
indep=function(riskFile=null,cliFile=null,uniOutFile=null,multiOutFile=null,uniForest=null,multiForest=null){
	risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)    #读取风险文件
	cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)      #读取临床文件
	
	#数据合并
	sameSample=intersect(row.names(cli),row.names(risk))
	risk=risk[sameSample,]
	cli=cli[sameSample,]
	rt=cbind(futime=risk[,1], fustat=risk[,2], cli, riskScore=risk[,(ncol(risk)-1)])
	
	#单因素独立预后分析
	uniTab=data.frame()
	for(i in colnames(rt[,3:ncol(rt)])){
		 cox <- coxph(Surv(futime, fustat) ~ rt[,i], data = rt)
		 coxSummary = summary(cox)
		 uniTab=rbind(uniTab,
		              cbind(id=i,
		              HR=coxSummary$conf.int[,"exp(coef)"],
		              HR.95L=coxSummary$conf.int[,"lower .95"],
		              HR.95H=coxSummary$conf.int[,"upper .95"],
		              pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
		              )
	}
	write.table(uniTab,file=uniOutFile,sep="\t",row.names=F,quote=F)
	bioForest(coxFile=uniOutFile, forestFile=uniForest, forestCol="green")

	#多因素独立预后分析
	uniTab=uniTab[as.numeric(uniTab[,"pvalue"])<1,]
	rt1=rt[,c("futime", "fustat", as.vector(uniTab[,"id"]))]
	multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
	multiCoxSum=summary(multiCox)
	multiTab=data.frame()
	multiTab=cbind(
	             HR=multiCoxSum$conf.int[,"exp(coef)"],
	             HR.95L=multiCoxSum$conf.int[,"lower .95"],
	             HR.95H=multiCoxSum$conf.int[,"upper .95"],
	             pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
	multiTab=cbind(id=row.names(multiTab),multiTab)
	write.table(multiTab,file=multiOutFile,sep="\t",row.names=F,quote=F)
	bioForest(coxFile=multiOutFile, forestFile=multiForest, forestCol="red")
}

#调用函数，进行独立预后分析
indep(riskFile="trainRisk.txt",
      cliFile="clinical.txt",
      uniOutFile="uniCox.txt",
      multiOutFile="multiCox.txt",
      uniForest="uniForest.pdf",
      multiForest="multiForest.pdf")


library(rms)                   #引用包
riskFile="trainRisk.txt"       #风险输入文件
cliFile="clinical.txt"         #临床数据文件
setwd("C:\\Users\\lexb4\\Desktop\\sCell\\36.Nomo")     #修改工作目录

#读取风险输入文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
risk=risk[,c("futime", "fustat", "riskScore")]

#读取临床数据文件
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
samSample=intersect(row.names(risk), row.names(cli))
risk1=risk[samSample,,drop=F]
cli=cli[samSample,,drop=F]
rt=cbind(risk1, cli)
paste(colnames(rt)[3:ncol(rt)],collapse="+")

#数据打包
dd <- datadist(rt)
options(datadist="dd")
#生成函数
f <- cph(Surv(futime, fustat) ~ riskScore+Age+Gender+Grade+Stage+T+M+N, x=T, y=T, surv=T, data=rt, time.inc=1)
surv <- Survival(f)
#建立nomogram
nom <- nomogram(f, fun=list(function(x) surv(1, x), function(x) surv(3, x), function(x) surv(5, x)), 
    lp=F, funlabel=c("1-year survival", "3-year survival", "5-year survival"), 
    maxscale=100, 
    fun.at=c(0.99, 0.9, 0.8, 0.7, 0.5, 0.3,0.1,0.01))  

#列线图可视化
pdf(file="Nomogram.pdf", width=9.5, height=7.5)
plot(nom)
dev.off()

#校准曲线
time=3    #预测年限
f <- cph(Surv(futime, fustat) ~ riskScore+Age+Gender+Grade+Stage+T+M+N, x=T, y=T, surv=T, data=rt, time.inc=time)
cal <- calibrate(f, cmethod="KM", method="boot", u=time, m=80, B=1000)
pdf(file="calibration.pdf", width=9, height=8.5)
plot(cal,
	 xlim=c(0,1),
	 ylim=c(0,1),
	 xlab=paste0("Nomogram-Predicted Probability of ", time, "-Year OS"),
	 ylab=paste0("Actual ", time, "-Year OS(proportion)"), lwd=1.5,
	 col="red", sub=T)
dev.off()




