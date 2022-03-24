library(limma)
setwd("D:\\summer.xu\\GSE13052x\\GSE13052data")
target <- read.table(file = "mRNA_exp13052.txt", sep = "\t", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
group1 <-  read.table(file = "group13052.txt", sep = "\t", header = TRUE,  stringsAsFactors = FALSE)
GSMgroup <- read.table(file = "GSM-group.txt", sep = "\t", header = TRUE,  stringsAsFactors = FALSE)
GSMNU <- read.table(file = "GSM-NU.txt", sep = "\t", header = TRUE,  stringsAsFactors = FALSE)

aa <- merge(GSMgroup,GSMNU,by = "GSM")

target <- target[,c(9,1,3,14:15,17:19,23:25,2,4:8,10:12,16,20:22,26:30,13)]
target <- t(apply(target, 1, function(x){log2(x+1)}))
target <- log2(target)
group <- factor(c(rep("Healthy",11), rep("Dengue",19)))
group2 <- group1[,2]
design <- model.matrix(~0+factor(group))  
colnames(design) <- levels(factor(group))
rownames(design) <- group1[,1]
contrast.matrix<-makeContrasts("Healthy-Dengue",levels=design)
fit <- lmFit(target,design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2) 
tempOutput = topTable(fit2, coef=1, n=Inf)#
nrDEG = na.omit(tempOutput) 
head(nrDEG)
write.csv(nrDEG,"final_limma_GSE13052.csv",quote = F)
