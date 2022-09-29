library(data.table)
library(caret)
library(caretEnsemble)
library(dplyr)
library(ranger)
library(ISLR)
library(caretEnsemble)
library(MLeval)
library(seqinr)
library(kmer)

set.seed(7)

## Remove a notação científica do R
options(scipen=999)

## Carregar as variáveis da aula
load("performance_Aula3.RData")

## Pacote MLeval de avaliação das métricas do Machine Learning
# É preciso dar de entrada o modelo gerado 
# ctrl <- trainControl(method="cv", summaryFunction=twoClassSummary, classProbs=T, savePredictions = T)
# model_rf <- train(classes ~ ., data=train, method="rf", preProc=c("center", "scale"), trControl=ctrl)
# pred_rf <- predict(model_rf, test)
# metrics_rf <- confusionMatrix(pred_rf, as.factor(test$classes))
# model_knn <- train(classes ~ ., data=train, method="knn", preProc=c("center", "scale"), trControl=ctrl)
# pred <- predict(model_knn, test)
# metrics_knn <- confusionMatrix(pred, as.factor(test$classes))

# Comparação dos resultados
results <- resamples(list(rf=model_rf, knn=model_knn))
summary(results)
dotplot(results)
scales <- list(x=list(relation="free"), y=list(relation="free"))
densityplot(results, scales=scales, pch = "|")

# Curvas ROCs e valores de AUC
res <- evalm(list(model_rf, model_knn))
res$roc ###plot curva roc e AUC
res$proc ###plot precision/recall
res$prg ###plot precision/recall gain
res$cc ###plot calibration curve
res$probs
res$optres
res$stdres

### Teste com outras sequências (sequência aleatória, sequência de outro vírus e sequência de homo sapiens)
#BioSequences to numbers (igual ao script da Aula 1)

new_test <- read.fasta(file="./Aula3/test.fasta", set.attributes = T)
AccessionNum_new_test <- getAnnot(new_test)
AccessionNum_new_test <- unname(unlist(AccessionNum_new_test))
AccessionNum_new_test <- gsub('>','',AccessionNum_new_test)
new_test <- kcount(new_test, k = 2, residues = NULL, gap = "-", named = TRUE,compress = TRUE, encode = FALSE) 

pred_rf <- predict(model_rf, new_test)
pred_knn <- predict(model_knn, new_test)

# Organizar a saída
out <- cbind(AccessionNum_new_test, pred_rf)
out <- cbind(out, pred_knn)
out <- gsub('1','mosquito.vir',out)
out <- gsub('2','plant.vir',out)
out

#save.image("performance_Aula3.RData")
