library(dplyr)
library(data.table)
library(caret)
#install.packages(factoextra)
#library(factoextra)

load("./data/handling_data_Aula1.RData")

### Explorar os k-mers 

mosquito_vir <- read.csv(file="./data/mosquito_vir_k2.csv", header=T, stringsAsFactors = T)
names(mosquito_vir)[1] <- "ID"
mosquito_vir <- as.data.frame(mosquito_vir)
plant_vir <- read.csv(file="./data/plant_vir_k2.csv", header=T, stringsAsFactors = T)
names(plant_vir)[1] <- "ID"
plant_vir <- as.data.frame(plant_vir)

#### Atribuição das classes (Classificação supervisionada)

classes <- replicate(nrow(mosquito_vir), "mosquito.vir")
mosquito_vir <- cbind(mosquito_vir, classes)
classes <- replicate(nrow(plant_vir), "plant.vir")
plant_vir <- cbind(plant_vir, classes)

# O rotulo das classes não deve ter espaços ou underscore. Deve ser separado apenas por ponto. 

### Construção da matriz de características 

featureMatrix <- rbind(mosquito_vir, plant_vir)# Mesma matriz dos slides
featureMatrix <- as.data.frame(featureMatrix)

### Explorando as matrizes

head(featureMatrix)
prop.table(table(featureMatrix$classes)) # Porcentagem de vírus por classe
barplot(prop.table(table(featureMatrix$classes)))

heatmap(as.matrix(select_if(featureMatrix, is.numeric))) # Variação das contagens
heatmap(as.matrix(select_if(featureMatrix[featureMatrix$classes == 'mosquito.vir',], is.numeric)))
heatmap(as.matrix(select_if(featureMatrix[featureMatrix$classes == 'plant.vir',], is.numeric)))

barplot((as.matrix(select_if(featureMatrix, is.numeric))), main="Contagens de k-mers", xlab="k-mers", ylab = "Contagens")
barplot((as.matrix(select_if(featureMatrix[featureMatrix$classes == 'mosquito.vir',], is.numeric))), main="Contagens de k-mers para os vírus de mosquito", xlab="k-mers", ylab = "Contagens")

plot(at ~ ac, data=as.matrix(select_if(featureMatrix, is.numeric)), main="Correlação entre as contagens do k-mer AC e do AT", xlab="AC", 
     ylab="AT", col=c("blue", "red")) ### Scatter plot entre dois k-mers

### Filtragem e pré-processamento da matriz 

sapply(featureMatrix, function(x) sum(is.na(featureMatrix))) # Quantidade de NAs por coluna
featureMatrix <- featureMatrix[complete.cases(featureMatrix), ] ## remove NAs
featureMatrix_dup <- featureMatrix[duplicated(featureMatrix), ] ## remove instâncias duplicadas
featureMatrix_subset <- subset(featureMatrix, select = c(ID, aa, ac, ag, at, classes)) # seleciona apenas algumas colunas
featureMatrix_filt <- featureMatrix[,colSums(featureMatrix[,2:17]) > 2] # Filtra k-mers que apresentam baixa contagens (pouco informativos)

### Imputation of NA values (substituir NA por outros valores em vez de os remover)
# Utilizando a média
featureMatrix_imp <- featureMatrix %>% mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))

### Remoção de características altamente correlacionadas

num <- subset(featureMatrix, select = -c(ID,classes))
num <- lapply(num, function(x) as.numeric(as.character(x)))
num <- as.data.frame(num)
descrCor <- cor(num)
summary(descrCor[upper.tri(descrCor)])
highlyCorrelated <- findCorrelation(descrCor, cutoff=1)
highlyCorCol <- colnames(num)[highlyCorrelated]
featureMatrix_corr <- featureMatrix[, -which(colnames(featureMatrix) %in% highlyCorCol)]
dim(featureMatrix_corr)

### Características com pouca variação entre as classes (baixo poder preditivo)

preproc = preProcess(featureMatrix, method = "nzv", freqCut = 2, uniqueCut = 20)
featureMatrix_remove = predict(preproc, featureMatrix)
dim(featureMatrix_remove)

### Identificação das características mais importantes

control <- trainControl(method="repeatedcv", number=10, repeats=3)
model <- train(classes~., data=featureMatrix[,-1], method="lvq", preProcess="scale", trControl=control)
importance <- varImp(model, scale=FALSE) # estimar a importância das características
print(importance)
plot(importance)

### Seleçao de característica tendo como base métricas de classificação

featureMatrix <- featureMatrix[,-1]
train_index <- createDataPartition(featureMatrix$classes, p=0.8, list=FALSE)
train <- featureMatrix[train_index, ]
test <- featureMatrix[-train_index, ]
subsets <- c(16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1) ## Subconjuntos das características
subsets <- c(2)

rctrl1 <- rfeControl(method = "cv", number = 10, returnResamp = "all", functions = caretFuncs, saveDetails = TRUE)
lmProfile <- rfe(classes ~ ., data = train, sizes = subsets, method = "rf", trControl = trainControl(method = "cv",classProbs = TRUE), rfeControl = rctrl1)
lmProfile$fit
lmProfile$results
lmProfile
lmProfile$optVariables

print(lmProfile)

ggplot(data = lmProfile, metric = "Accuracy") + theme_bw()

varimp_data <- data.frame(feature = row.names(varImp(lmProfile))[1:10],importance = varImp(lmProfile)[1:10, 1])
p <- ggplot(data = varimp_data, aes(x = reorder(feature, importance), y = importance, fill = feature)) +
  geom_bar(stat="identity") + labs(x = "Features", y = "Variable Importance") + 
  geom_text(aes(label = round(importance, 2)), vjust=1.6, color="black", size=4) + 
  theme_classic() + theme(legend.position = "none") + coord_flip()
p

### Salvar a matriz de características final
write.csv(featureMatrix, file="featureMatrix.csv", row.names=F)

save.image("./data/handling_data_Aula1.RData")
