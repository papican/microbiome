# Comparison with other classifiers
### Remove rare OTUs
bac.clean.nolog.100reads <- prune_taxa(taxa_sums(bac.clean.ss) > 500, bac.clean.nolog)

# Random Forest Classification
## classification
#### All machine learning
# (1) Random forest
library(randomForest)

#Control vs RA
table1 <- otu_table(bac.clean.nolog.100reads)
table1 <- t(table1)
rownames(table1)

result1 <- c(rep(0,8),rep(1, 8)) 

table1 <- cbind(table1, result1)
dataset_diagnosis <- as.data.frame(table1)

## get rid of X in front
#colnames(dataset) <- gsub("^X(.{32})", "\\1",colnames(dataset))
sort(colnames(dataset_diagnosis))

#write.xlsx(dataset,'rf_dataset.xlsx')

# Encoding the target feature as factor
dataset_diagnosis$result1 = factor(dataset_diagnosis$result1, levels = c(0, 1))
dataset_diagnosis$result1

# Splitting the dataset into the Training set and Test set
# install.packages('caTools')
library(caTools)

set.seed(123)
split = sample.split(dataset_diagnosis$result, SplitRatio = 0.66)
training_set_diagnosis = subset(dataset_diagnosis, split == TRUE)
test_set_diagnosis = subset(dataset_diagnosis, split == FALSE)
dim(training_set_diagnosis)   # when split ratio 0.66 -> 10
dim(test_set_diagnosis)  # when split ratio 0.33 -> 6


set.seed(123)
classifier_diagnosis = randomForest(x = training_set_diagnosis[-ncol(dataset_diagnosis)],
                                    y = training_set_diagnosis$result, ntree=1000)
rf_pred_diagnosis = predict(classifier_diagnosis, newdata = test_set_diagnosis[-ncol(dataset_diagnosis)])
cm = table(test_set_diagnosis[, ncol(dataset_diagnosis)], rf_pred_diagnosis)
print(cm)

##Accuracy calculation based on AUC
library(ROCR)
predictions.rf=as.vector(as.numeric(rf_pred_diagnosis))
pred.rf=prediction(predictions.rf, test_set_diagnosis[,ncol(dataset_diagnosis)])
AUC.rf=performance(pred.rf,"auc") #Calculate the AUC value
(AUC.rf=AUC.rf@y.values[[1]]) #0.5


##K-fold cross validation
## (1) random forest cross validation
library(caret)
folds = createFolds(training_set_diagnosis$result, k = 10)
cv = lapply(folds, function(x) {
  training_fold = training_set_diagnosis[-x, ]
  test_fold = training_set_diagnosis[x, ]
  classifier = randomForest(x = training_fold[-ncol(dataset_diagnosis)],
                            y = training_fold$result, ntree=2000)
  
  # classifier = svm(formula = result ~ .,
  #                  data = training_fold,
  #                  type = 'C-classification',
  #                  kernel = 'radial')
  y_pred = predict(classifier, newdata = test_fold[-ncol(dataset_diagnosis)])
  cm = table(test_fold[, ncol(dataset_diagnosis)], y_pred)
  accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})

accuracy = mean(as.numeric(cv)) 
accuracy   ## 0.7638889
(sd <- sd(as.numeric(cv))) ## 0.01464017



##Error rate
colnames(classifier_diagnosis$err.rate) <- c('OOB','Cynanchumwilfordii','Cynanchumauriculatum')
plot(classifier_diagnosis, main = 'Out-of-bag Error estimate')
legend("topright", colnames(classifier_diagnosis$err.rate),col=1:3,cex=0.8,fill=1:3)
oob <- classifier_diagnosis$err.rate[,1]
oob.oob <- oob[length(oob)]
legend("bottomright", paste0('OOB : ',round(as.numeric(oob.oob), digits = 4)*100,'%'))



### Fungal model
sample.order<-c("F1","F2","F3","F4","F5","F6","F7","F8","F9","F10","F11","F12",
                "F13","F14","F15","F16","F17","F18","F19","F20","F21","F22","F23",
                "F24","F25","F26","F27","F28","F29","F30","F36","F37","F38","F39",
                "F40","F41","F42","F43", "F44","F45","F46","F47","F48","F49","F50",
                "F51","F52","F53","F54","F55","F56","F57","F58","F59","F60","F61",
                "F62","F63","F64","F65","F71","F72","F73","F74","F75","F76","F77",
                "F78","F79","F81","F82","F83","F84","F87","F88","F89","F90","F91",
                "F93","F94","F96","F97","F98","F99","F100","F101","F102","F103","F104",
                "F105","F106","F107","F108","F109","F110","F111","F112","F113",
                "F114","F115","F116","F117","F118","F119","F120","F121","F122","F123",
                "F124","F125","F126","F128","F129","F130","F131","F132","F133","F134",
                "F135","F136","F137","F138","F139","F140","F141","F142","F143",
                "F144","F145")

table2 <- otu_table(fun.its1.clean.ss)[,sample.order]
table2 <- t(table2)
rownames(table2)
result2 <- c(rep(0,30),rep(1, 99)) 
# result2 <- c(rep(0,2),rep(1, 10), rep(0,1), rep(1,10), rep(0,1),rep(1,9),rep(0,1),rep(1,10),rep(0,1), 
#              rep(1,6),rep(0,18), rep(1,4),rep(0,1),rep(1,10), rep(0,1), rep(1,10), rep(0,1), rep(1,6), rep(0,1),rep(1,9),rep(0,1),rep(1,7), rep(0,1), rep(1,8)) 

table2 <- cbind(table2, result2)
dataset_diagnosis.f <- as.data.frame(table2)
dataset_diagnosis.f$result2
sort(colnames(dataset_diagnosis.f))

## get rid of X in front
#colnames(dataset) <- gsub("^X(.{32})", "\\1",colnames(dataset))
sort(colnames(dataset_diagnosis.f))

#write.xlsx(dataset,'rf_dataset.xlsx')

# Encoding the target feature as factor
dataset_diagnosis.f$result2 = factor(dataset_diagnosis.f$result2, levels = c(0, 1))
dataset_diagnosis.f$result2

# Splitting the dataset into the Training set and Test set
# install.packages('caTools')
library(caTools)

set.seed(123)
split.f = sample.split(dataset_diagnosis.f$result2, SplitRatio = 0.66)
training_set_diagnosis.f = subset(dataset_diagnosis.f, split.f == TRUE)
test_set_diagnosis.f = subset(dataset_diagnosis.f, split.f == FALSE)
dim(training_set_diagnosis.f) 
dim(test_set_diagnosis.f)  


set.seed(123)
classifier_diagnosis.f = randomForest(x = training_set_diagnosis.f[-ncol(dataset_diagnosis.f)],
                                      y = training_set_diagnosis.f$result2, ntree=2000)
rf_pred_diagnosis.f = predict(classifier_diagnosis.f, newdata = test_set_diagnosis.f[-ncol(dataset_diagnosis.f)])
cm.f = table(test_set_diagnosis.f[, ncol(dataset_diagnosis.f)], rf_pred_diagnosis.f)
print(cm.f)

##Accuracy calculation based on AUC
library(ROCR)
predictions.rf.f=as.vector(as.numeric(rf_pred_diagnosis.f))
pred.rf.f=prediction(predictions.rf.f, test_set_diagnosis.f[,ncol(dataset_diagnosis.f)])
AUC.rf.f=performance(pred.rf.f,"auc") #Calculate the AUC value
(AUC.rf.f=AUC.rf.f@y.values[[1]]) #0.5


##K-fold cross validation
## (1) random forest cross validation
library(caret)
folds.f = createFolds(training_set_diagnosis.f$result2, k = 10)
cv.f = lapply(folds.f, function(x) {
  training_fold = training_set_diagnosis.f[-x, ]
  test_fold = training_set_diagnosis.f[x, ]
  classifier = randomForest(x = training_fold[-ncol(dataset_diagnosis.f)],
                            y = training_fold$result2, ntree=2000)
  
  # classifier = svm(formula = result ~ .,
  #                  data = training_fold,
  #                  type = 'C-classification',
  #                  kernel = 'radial')
  y_pred = predict(classifier, newdata = test_fold[-ncol(dataset_diagnosis.f)])
  cm = table(test_fold[, ncol(dataset_diagnosis.f)], y_pred)
  accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})

accuracy.f = mean(as.numeric(cv.f)) 
accuracy.f  ## 0.7638889 (ITS2) / 0.7861111 (ITS1)
(sd <- sd(as.numeric(cv.f))) ## 0.01464017 (ITS2) / 0.05563266 (ITS1)



##Error rate
colnames(classifier_diagnosis.f$err.rate) <- c('OOB','Control','RA')
plot(classifier_diagnosis.f, main = 'Out-of-bag Error estimate')
legend("topright", colnames(classifier_diagnosis.f$err.rate),col=1:3,cex=0.8,fill=1:3)
oob <- classifier_diagnosis.f$err.rate[,1]
oob.oob <- oob[length(oob)]
legend("bottomright", paste0('OOB : ',round(as.numeric(oob.oob), digits = 4)*100,'%'))

dev.off()


### Mininum number of OTUs for explaining the differences between conventional farming and no fertilization
### bacteria
cv.bac <- rfcv(training_set_diagnosis[-ncol(dataset_diagnosis)],training_set_diagnosis$result1, cv.fold=10, step=0.9)
cv.bac$error.cv 
#332 299 269 242 218 196 176 159 143 129 116 104  94  84  76  68  62  55  50  45  40  36  33  29  26 
#0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 
#24  21  19  17  16  14  13  11  10   9   8   7   6   5   4   3   2   1 
#0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.0 0.1 0.1 0.2 0.2 0.3

mean(cv.bac$error.cv) #0.1069767

### fungi
cv2.fun <- rfcv(training_set_diagnosis.f[-ncol(dataset_diagnosis.f)],training_set_diagnosis.f$result2, cv.fold=10, step=0.9)
cv2.fun$error.cv
mean(cv2.fun$error.cv) #0.2334686

### plot

par(mfrow=c(1,1), mar=c(5,5,5,2)+0.1)

##Split plot
#Bacteria
plot(cv.bac$n.var,cv.bac$error.cv,type='o', lty=2, col=2, xlab='Number of OTUs included in RF model',
     ylab='Cross-validation error', ylim=c(0,0.30), xaxt='n')
axis(side = 1, at=seq(0,18292,2000))
abline(h=0.1249567, col="firebrick2")
abline(v=1000, col="limegreen")

#Fungi
plot(cv2.fun$n.var,cv2.fun$error.cv,type='o', lty=2, col=2, xlab='Number of OTUs included in RF model',
     ylab='Cross-validation error', ylim=c(0,0.30), xaxt='n')
axis(side = 1, at=seq(0,1410,70))
abline(h=mean(cv2.fun$error.cv), col="firebrick2")
abline(v=40, col="limegreen")


df.cv1 <- data.frame(cv.bac$n.var,cv.bac$error.cv)
df.cv2 <- data.frame(cv2.fun$n.var,cv2.fun$error.cv)
df.cv3 <- data.frame(cv2.fun$n.var,cv2.fun$error.cv)

write.xlsx(df.cv1, 'bacteria-otu-cverror.xlsx')
write.xlsx(df.cv2, 'fungi-otu-cverror.xlsx')
write.xlsx(df.cv3,'fungi-otu-cverror-its1.xlsx')





#### practice


## let's upgrade one time
# deepmodel <- h2o.deeplearning(y='result',
#                               training_frame = as.h2o(training_set),
#                               validation_frame = as.h2o(training_set),         
#                               #standardize = T,        
#                               model_id = "deep_model",        
#                               activation = "Rectifier",       
#                               epochs = 100,       
#                               seed = 1,       
#                               hidden = 5,         
#                               variable_importances = T,
#                               train_samples_per_iteration = -2)
# h2o.performance(deepmodel,valid = T)
# # Predicting the Test set results
# ann_pred = h2o.predict(deepmodel, newdata = as.h2o(test_set[-ncol(dataset)]))
# ann_pred
# 
# # Making the Confusion Matrix
# cm = table(test_set[, ncol(dataset)], as.vector(ann_pred$predict))   ## to make table, we need vector
# cm
# 
# #### second, parameter grid search
# 
# #set parameter space
# activation_opt <- c("Rectifier","RectifierWithDropout", "Maxout","MaxoutWithDropout")
# hidden_opt <- list(c(10,10),c(20,15),c(50,50,50))
# l1_opt <- c(0,1e-3,1e-5)
# l2_opt <- c(0,1e-3,1e-5)
# 
# hyper_params <- list( activation=activation_opt,
#                       hidden=hidden_opt,
#                       l1=l1_opt,
#                       l2=l2_opt )
# 
# #set search criteria
# search_criteria <- list(strategy = "RandomDiscrete", max_models=10)
# 
# #train model
# dl_grid <- h2o.grid("deeplearning"
#                     ,grid_id = "deep_learn"
#                     ,hyper_params = hyper_params
#                     ,search_criteria = search_criteria
#                     ,training_frame = as.h2o(training_set)
#                     ,y='result'
#                     ,nfolds = 5
#                     ,epochs = 100)
# 
# #get best model
# d_grid <- h2o.getGrid("deep_learn",sort_by = "accuracy")
# best_dl_model <- h2o.getModel(d_grid@model_ids[[1]])
# h2o.performance (best_dl_model,xval = T)
# 
# 
# ann_pred = h2o.predict(best_dl_model, newdata = as.h2o(test_set[-ncol(dataset)]))
# ann_pred
# 
# # Making the Confusion Matrix
# cm = table(test_set[, ncol(dataset)], as.vector(ann_pred$predict))   ## to make table, we need vector
# cm
