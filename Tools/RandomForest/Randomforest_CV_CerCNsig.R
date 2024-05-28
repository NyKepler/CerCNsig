## (1) Define the packages that will be needed
packages <- c('fs', 'dplyr', 'ggplot2', 'caret', 'ROCR', 'randomForest')

## (2) Install them if not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

## (3) Load the packages into R session
invisible(lapply(packages, library, character.only = TRUE))


## 1. Import data
folder.name <- "Randomforest_CNsig" 
setwd(paste(path_home_r(), "CerCNsig", folder.name, sep = "/"))

#' read RDS files and combine benign and RRSO group into one group called 'Benign'
Benign_VS_good.SxCMat <- readRDS("Benign_VS_good.SxCMat.rds") %>% as.data.frame()
Benign_VS_good.SxCMat$Class <- "Benign"
HGSC_VS_good.SxCMat <- readRDS("HGSC_VS_good.SxCMat.rds") %>% as.data.frame()
HGSC_VS_good.SxCMat$Class <- "HGSC"


VS <-rbind(Benign_VS_good.SxCMat, HGSC_VS_good.SxCMat) #' total 163 samples
#' VS$sample <- paste(VS$Class, row.names(VS), sep = "_")
#' VS$sample <- NULL
saveRDS(VS, "HGSC_Benign_combined_SxCMat.rds")

## 2. k-fold partition for cross validation 
#' Set seed for reproducibility
set.seed(1220)

#' Define repeated cross validation with 10 folds and 100 repeats
repeat_cv <- trainControl(method='repeatedcv', number=10, repeats=100, search = "grid")

#' Split the data so that we use 70% of it for training
train_index <- createDataPartition(y=VS$Class, p=0.7, list=FALSE)

#' Subset the data
training_set <- VS[train_index, ] #' 115
testing_set <- VS[-train_index, ] #' 48

## 3. Create randomforest

#' Train a random forest model
forest <- train(
  
  #' Source of data; 
  data=training_set, 
  
  #' Formula. We are using all variables to predict Species
  Class~., 
  
  #' Number of trees to grow. 
  #' Ensure that every input row gets predicted at least a few times.
  #' Here we set 10 times.
  #' ntree = sqrt (number of row * number of columns)/numberofcpu
  
  #' `rf` method for random forest
  method='rf', 

  #' Add repeated cross validation as trControl
  trControl=repeat_cv,
  
  #' Accuracy to measure the performance of the model
  metric='Accuracy')

#' Print out the details about the model
sink("RF_Model.txt")
forest 
print("Final Model")
forest$finalModel
sink()
saveRDS(forest, "RF_Model_CerCNsig.rds")

## Get variable importance, and turn into a data frame
var_imp <- varImp(forest, scale=FALSE)$importance
var_imp <- data.frame(variables=row.names(var_imp), importance=var_imp$Overall)

## Create a plot of variable importance
pdf("Importance.pdf")
var_imp %>%
  
  ## Sort the data by importance
  arrange(importance) %>%
  
  ## Create a ggplot object for aesthetic
  ggplot(aes(x=reorder(variables, importance), y=importance, fill=importance <1)) + 
  
  ## Plot the bar graph
  geom_bar(stat='identity') + 
  scale_fill_manual(values = c("#FC2D00","#008EFC"),
                    labels=c("TRUE"= "Importance < 1",
                             "FALSE"= "Importance > 1")) +
  
  ## Flip the graph to make a horizontal bar plot
  coord_flip() + 
  
  ## Add x-axis label
  xlab('Component') +
  ylab('Importance') +
  labs(fill = "") +
  
  ## Add a title
  labs(title='Prediction Importance') + 
  
  ## Some layout for the plot
  theme_minimal() + 
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 15), 
        plot.title = element_text(size = 20), 
  )
dev.off()

## Generate predictions
RF_prediction <- predict(
  
  ## Random forest object
  object=forest, 
  
  ## Data to use for predictions; remove the Class
  newdata=testing_set[, -40])

## Print the accuracy (AUC)
sink("RF_AUC.txt")
accuracy <- mean(RF_prediction == testing_set$Class)*100
cat('Accuracy on testing data: ', round(accuracy, 2), '%',  sep='')

##  AUC for the testing data

pred.test=predict(forest,
              type = "prob",
              newdata = testing_set[,-40])

ROCR.pred.test = prediction(pred.test[,2], testing_set$Class)
# 1. Area under curve
auc.test = performance(ROCR.pred.test, "auc")
print("AUC_Test")
auc.test@y.values 

#' 2. True Positive and Negative Rate
ROCR.perf.test = performance(ROCR.pred.test, "tpr","fpr")
saveRDS(ROCR.perf.test, "RF_Performance_Test.rds")


##  AUC for the training set
pred.train=predict(forest,
              type = "prob",
              newdata = training_set[, -40])

ROCR.pred.train = prediction(pred.train[,2], training_set$Class)
#' 1. Area under curve
auc.train = performance(ROCR.pred.train, "auc")
print("AUC_Train")
auc.train@y.values

#' 2. True Positive and Negative Rate
ROCR.perf.train = performance(ROCR.pred.train, "tpr","fpr")
sink()
saveRDS(ROCR.perf.train, "RF_Performance_Train.rds")

## calculate cutoff
sink("RF_Cutoff.txt")
cost_perf = performance(ROCR.pred.test, "cost")
print("1-Specificity")
ROCR.pred.test@cutoffs[[1]][which.min(cost_perf@y.values[[1]])] 
FPR<-ROCR.pred.test@cutoffs[[1]][which.min(cost_perf@y.values[[1]])] 

print("Sensitivity")
ROCR.pred.test@cutoffs[[1]][which(cost_perf@y.values[[1]]>FPR)] %>% min() #'TPR
sink()

threshold<-cbind(Actualvalue=testing_set$Class,Predictedvalue=pred.test)  
write.table(threshold, "RF_Prediction_Value_Test.csv", quote = F, sep = ",", 
            col.names = T, row.names = T)


## plot ROC 
pdf("ROC_Train_Test.pdf", onefile = T)
plot(ROCR.perf.train, main="ROC Curve for Random Forest",col=4,lwd=3, xlab = "1-Specificity", ylab= "Sensitivity", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(ROCR.perf.test, add= T, main="ROC Curve for Random Forest",col=3,lwd=3)
abline(a=0,b=1,lwd=2,lty=2,col="black")
dev.off()

pdf("ROC_Models_Comparison.pdf", onefile = T)
ROCR.perf.test.filt <- readRDS("/home/researcher/CerCNsig/Randomforest_CNsig_filt/RF_Performance_Test.rds")
plot(ROCR.perf.test, main="ROC Curve for Random Forest",col=3,lwd=3, xlab = "1-Specificity", ylab= "Sensitivity", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(ROCR.perf.test.filt, add= T, main="ROC Curve for Random Forest",col=2,lwd=3) # filtered
abline(a=0,b=1,lwd=2,lty=2,col="black")
dev.off()

#' Plot cut off and prediction 
pdf("ROC_with_cutoff_gradient.pdf", onefile = T)
plot(ROCR.perf.test,colorize=TRUE,print.cutoffs.at=seq(0.1,by=0.1))
plot(ROCR.perf.test,colorize=TRUE)
dev.off()

#threshold<-caret::confusionMatrix(testing_set,Actualvalue=testing_set$Class,Predictedvalue=pred.test>FPR) 

sessionInfo() %>% capture.output(file="session_info.txt")

