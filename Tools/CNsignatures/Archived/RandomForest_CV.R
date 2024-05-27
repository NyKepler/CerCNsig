## (1) Define the packages that will be needed
packages <- c('dplyr', 'ggplot2', 'caret', 'ROCR')

## (2) Install them if not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

## (3) Load the packages into R session
invisible(lapply(packages, library, character.only = TRUE))

#' K-fold partitioning
Benign_VS_good.SxCMat <- readRDS("~/CerCNsig/Randomforest/Benign_VS_good.SxCMat.rds") %>% as.data.frame()
Benign_VS_good.SxCMat$Class <- "Benign"
HGSC_VS_good.SxCMat <- readRDS("~/CerCNsig/Randomforest/HGSC_VS_good.SxCMat.rds") %>% as.data.frame()
HGSC_VS_good.SxCMat$Class <- "HGSC"
RRSO_VS.SxCMat <- readRDS("~/CerCNsig/Randomforest/RRSO_VS.SxCMat.rds") %>% as.data.frame()
RRSO_VS.SxCMat$Class <- "Benign"

VS <-rbind(Benign_VS_good.SxCMat, HGSC_VS_good.SxCMat) 
VS$sample <- paste(VS$Class, row.names(VS), sep = "_")
VS <- VS[, -41]
saveRDS(VS, "HGSC_Benign_combined_SxCMat.rds")

#' k-fold
## Set seed for reproducibility
set.seed(1220)

## Define repeated cross validation with 5 folds and three repeats
repeat_cv <- trainControl(method='repeatedcv', number=10, repeats=100, search = "grid")

## Set seed for reproducibility
set.seed(1220)

## Split the data so that we use 70% of it for training
train_index <- createDataPartition(y=VS$Class, p=0.7, list=FALSE)

## Subset the data
training_set <- VS[train_index, ]
testing_set <- VS[-train_index, ]

## Create randomforest
## Set seed for reproducibility
set.seed(1220)
library(randomForest)

## Train a random forest model
forest <- train(
  
  # Formula. We are using all variables to predict Species
  Class~., 
  
  # Source of data; remove the Species variable
  data=training_set, 
  
  # `rf` method for random forest
  method='rf', 
  
  # Add repeated cross validation as trControl
  trControl=repeat_cv,
  
  # Accuracy to measure the performance of the model
  metric='Accuracy')

## Print out the details about the model
forest 

forest$finalModel

## Get variable importance, and turn into a data frame
var_imp <- varImp(forest, scale=FALSE)$importance
var_imp <- data.frame(variables=row.names(var_imp), importance=var_imp$Overall)

## Create a plot of variable importance
var_imp %>%
  
  ## Sort the data by importance
  arrange(importance) %>%
  
  ## Create a ggplot object for aesthetic
  ggplot(aes(x=reorder(variables, importance), y=importance, fill=importance <1)) + 
  
  ## Plot the bar graph
  geom_bar(stat='identity') + 
  scale_fill_manual(values = c("#FC2D00","#008EFC"),
                    labels=c("TRUE"= "importance < 1",
                             "FALSE"= "importance > 1")) +
  
  ## Flip the graph to make a horizontal bar plot
  coord_flip() + 
  
  ## Add x-axis label
  xlab('Component') +
  labs(fill = "") +
  
  ## Add a title
  labs(title='CNcomponents importance') + 
  
  ## Some layout for the plot
  theme_minimal() + 
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 15), 
        plot.title = element_text(size = 20), 
  )

## Generate predictions
y_hats <- predict(
  
  ## Random forest object
  object=forest, 
  
  ## Data to use for predictions; remove the Class
  newdata=testing_set[, -40])

## Print the accuracy (AUC)
accuracy <- mean(y_hats == testing_set$Class)*100
cat('Accuracy on testing data: ', round(accuracy, 2), '%',  sep='')
## Accuracy on testing data: 83.33%

Actualvalue=test$PoorCare,Predictedvalue=pred_test>0.5

##  AUC for the testing data
pred.test=predict(forest,
              type = "prob",
              newdata = testing_set[,-40])

ROCR.perd.test = prediction(pred.test[,2], testing_set$Class)
# 1. Area under curve
auc.test = performance(ROCR.perd.test, "auc")
auc.test@y.values
# 2. True Positive and Negative Rate
ROCR.pref.test = performance(ROCR.perd.test, "tpr","fpr")

##  AUC for the training set
pred.train=predict(forest,
              type = "prob",
              newdata = training_set[, -40])

ROCR.perd.train = prediction(pred.train[,2], training_set$Class)
# 1. Area under curve
auc.train = performance(ROCR.perd.train, "auc")
auc.train@y.values
# 2. True Positive and Negative Rate
ROCR.pref.train = performance(ROCR.perd.train, "tpr","fpr")

## seperate then based on diagnostic time point on 
Randomization_240423_ny <- read_csv("CerCNsig/Randomization_240423_ny.csv")
diagnostic <-Randomization_240423_ny$Library[Randomization_240423_ny$Subgroup %in% "diagnostic"] %>% c()
prediagnostic.6 <-Randomization_240423_ny$Library[Randomization_240423_ny$Subgroup %in% "prediagnostic_>6"]
prediagnostic.0 <-Randomization_240423_ny$Library[Randomization_240423_ny$Subgroup %in% "prediagnostic_0-6"]

## prepare different validation/testing set
RRSO_VS.SxCMat
diag.HGSC.SxCMat <- HGSC_VS_good.SxCMat[row.names(HGSC_VS_good.SxCMat) %in% diagnostic,]
diag.6.HGSC.SxCMat <- HGSC_VS_good.SxCMat[row.names(HGSC_VS_good.SxCMat) %in% prediagnostic.6,]
diag.0.HGSC.SxCMat <- HGSC_VS_good.SxCMat[row.names(HGSC_VS_good.SxCMat) %in% prediagnostic.0,]

data <- rbind(RRSO_VS.SxCMat, diag.HGSC.SxCMat)


  pred=predict(forest,
                     type = "prob",
                     newdata = data[, -40])
  
  ROCR.pred = prediction(pred[,2], data$Class)
  # 1. Area under curve
  auc = performance(perf, "auc")
  auc@y.values
  # 2. True Positive and Negative Rate
  pref = performance(perf, "tpr","fpr")


## plot both 
plot(ROCR.pref.train, main="ROC Curve for Random Forest",col=4,lwd=5, xlab = "1-Specificity", ylab= "Sensitivity", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(ROCR.pref.test, add= T, main="ROC Curve for Random Forest",col=2,lwd=5)
plot(pref, add= T, main="ROC Curve for Random Forest",col=2,lwd=2)

abline(a=1,b=-1,lwd=2,lty=2,col="black")

## calculate cutoff
cost_perf = performance(ROCR.perd.test, "cost")
ROCR.perd.test@cutoffs[[1]][which.min(cost_perf@y.values[[1]])] # 0.546
ROCR.perd.test@cutoffs[[1]][which(cost_perf@y.values[[1]]>0.4)] # 0.546
#' 0.766


threshold<-cbind(Actualvalue=testing_set$Class,Predictedvalue=pred.test) # assuming thershold to be 0.29

#' use pROC
#' 
#' 
mroc<-roc_(testing_set, pred.train)
plot(ROCR.pref.test,colorize=TRUE,print.cutoffs.at=seq(0.1,by=0.05))
plot(ROCR.pref.test,colorize=TRUE,print.predictions.at=seq(0.1,by=0.1))

threshold<-caret::confusionMatrix(Actualvalue=testing_set$Class,Predictedvalue=pred.test>0.546) # assuming thershold to be 0.29
