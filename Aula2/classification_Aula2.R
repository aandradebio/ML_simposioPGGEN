library(data.table)
library(caret)
library(caretEnsemble)
library(dplyr)
library(ranger)
library(e1071)
library(adabag)
library(ISLR)
library(caretEnsemble)

## Remove a notação científica do R
options(scipen=999)

## Carrega as variáveis
load("classification_Aula2.RData")

## Carrega a matriz de características
featureMatrix <- read.csv(file="featureMatrix.csv", header=T, stringsAsFactors = T)

### Classificação supervisionada 
# Escolha do método de validação cruzada
crtl <- trainControl(method="repeatedcv", number=10, repeats=3,classProbs=T, savePredictions = 'final')

# Divisão das instâncias em treino e teste
train_index <- createDataPartition(featureMatrix$classes, p=0.8, list=FALSE)
train <- featureMatrix[train_index, ]
test <- featureMatrix[-train_index, ]

# Classificação com Random Forest
model_rf <- train(classes ~ ., data=train, method="rf", preProc=c("center", "scale"), trControl=crtl)

# Predição de dados desconhecidos
pred_rf <- predict(model_rf, test)

# Obtenção das métricas da predição
metrics_rf <- confusionMatrix(pred_rf, as.factor(test$classes))
metrics_rf$overall
metrics_rf$byClass

# Classificação com K-NN 
model_knn <- train(classes ~ ., data=train, method="knn", preProc=c("center", "scale"), trControl=crtl)
pred <- predict(model_knn, test)
metrics_knn <- confusionMatrix(pred, as.factor(test$classes))
metrics_knn$overall
metrics_knn$byClass

### Método Ensemble

model_list <- caretList(classes~., data=train,trControl=crtl,preProcess="center",methodList=c("ranger", "knn"))
model_list
stackControl <- trainControl(method="repeatedcv", number=10, repeats=3, savePredictions='final', classProbs=TRUE)
model_ensemble <- caretStack(model_list, method="rf", metric="Accuracy", trControl=stackControl)
model_ensemble
pred_ensemble <- predict(model_ensemble, test)
metrics_ensemble <- confusionMatrix(pred_ensemble, as.factor(test$classes))
metrics_ensemble$byClass
metrics_ensemble$overall
results <- resamples(model_list)
summary(results)
modelCor(results)
dotplot(results)
splom(results)

#save.image("classification_Aula2.RData")
