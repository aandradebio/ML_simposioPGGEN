#### Minicurso de Introdução a Aprendizagem de máquina para Bioinformática
### Msc. Amanda Araújo Serrão de Andrade 


####dendograma colorido por classe. Usa o info para colorir
tree <- read.tree(file="tree")
info <- read.csv(file="length.csv", header=T)
p <- ggtree(tree,branch.length="none") + layout_dendrogram() ##cria um dendograma simples
p2 <- p %<+% info + geom_tippoint(aes(color=classes)) + scale_color_manual(values=c("Arboviruses"="#66C2A5", "Mosquito-specific viruses"="#FC8D62", "Other viruses"="#8DA0CB"))
p3 <- p2 %<+% info + geom_text(aes(x=branch, label=species),size=3,hjust=0,color="black",check_overlap = T,family='Arial') ##acrescentar os nomes das espécies


library(data.table)
library(seqinr)
library(kmer)
library(dplyr)
library(e1071)
library(dplyr)
library(FactoMineR)
library(factoextra)
library(treeio)
library("FactoMineR")
library("factoextra")
library("wordspace")
library(stats)
library(pheatmap)
library(ape)
library(NbClust)
library(ggtree)


###funciona
fa <- seqinr::read.fasta("/home/amanda/Downloads/dataset_ML/arbo_mos/new_dataset/arbo/aino.fasta", set.attributes = T) ###fasta de entrada
kd <- kdistance(fa, k=5,method="euclidean") ###matriz de distancia diretamente dos kmers
tree <- as.dendrogram(hclust(kd, "complete")) ###clustering hierarquico feito a partir da matriz de distancia
plot(tree)

###posso usar as demais features. só dar o csv das freq entrada. sem as classes
kd <- read.csv(file="freq.csv", header=T)
dis <- distance(kd, method = "euclidean", use.row.names = TRUE, as.dist.obj = TRUE) ###matriz quadrada
hc <- hclust(kd, method = "complete")

plot(as.phylo(tree)) ###tree
plot(as.phylo(tree), type = "fan",cex = 0.6)
plot(tree, cex = 0.6, labels=F)
heatmap(as.matrix(tree)) 

###escreve a arvore
tree <- as.phylo(tree)
write.tree(tree, file="tree", tree.names = FALSE)


###organizar os dados pro k-means
k2 <- read.csv(file="k2_vImp.csv")
k2 <- subset(k2, select = -c(V1,classes))
k2 <- t(k2)
colnames(k2) <- k2[17,]
k2 <- k2[-17,] 
k2 <- as.data.frame(k2)

km <- rownames(k2)
k2 <- subset(k2, select = c(non.mosquito, mosquito))
k2 <- lapply(k2, function(x) as.numeric(as.character(x)))
k2 <- as.data.frame(k2)
row <- k2[,18]
rownames(k2) <- row
row <- make.names(row, unique = TRUE)

r <- as.data.frame(r)
subseqk2 <- k2[unlist(lapply(r[,1], function(x) grep(x, k2$V1, fixed = TRUE))),]
write.csv(subseqk2, file="subseq.csv", col.names=T, row.names=T)
subseqk2 <- read.csv(file="subseq.csv", stringsAsFactors = T)
subseqk2 <- subseqk2[!duplicated(subseqk2), ]
subseqk2 <- subseqk2[,-1]

###K-Means

dists=dist(t(k2))
mds=cmdscale(dists)
plot(mds,pch=19,col=rainbow(5))

set.seed(101)
kclu=kmeans(t(k2),centers=5)  
table(kclu$cluster)

getAnywhere("cluster_mat")
edit(pheatmap)
pheatmap(test)

###PCA
prcomp(k2, scale = FALSE) ###PCA
princomp(k2, cor = FALSE, scores = TRUE)
res.pca <- prcomp(k2, scale = TRUE)
fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)

###Escolher o melhor K do K-Means
fviz_nbclust(k2, kmeans, method = "wss") +  geom_vline(xintercept = 4, linetype = 2)+  labs(subtitle = "Elbow method")

k2 <- (t(k2))
kclu=kmeans(t(k2),centers=4)  
table(kclu$cluster)
type2kmclu = data.frame(Class=substr(colnames(k2),1,3), cluster=kclu$cluster)
table(type2kmclu)
kmclu=cluster::pam(t(k2),k=4)
dists=dist(t(k2))

# calculate MDS
mds=cmdscale(dists)

# plot the patients in the 2D space
plot(mds,pch=19,col=rainbow(4)[kclu$cluster])




# Qualquer dúvida ou problema na instalação, pode me mandar um email aandradebio@gmail.com.
# Dúvidas específicas devem ser acompanhadas do código, print da mensagem de erro e se possível de um conjunto de dados para teste.
