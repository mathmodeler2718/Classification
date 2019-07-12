datafile <- read.table("GSE13425-Bt-T.txt",check.names = F)
df<- as.data.frame(t(datafile))

#pkgs <- c("rpart", "rpart.plot", "party", "randomForest", "e1071", "class", "data.table", "tibble")
#install.packages(pkgs, depend=TRUE)
library(data.table)
library(tibble)

df = setDT(df, keep.rownames = TRUE)[]

df = add_column(df, class = 1:190, .after = 1)

i=1
for(i in 1:190){
  df[i,2] = ifelse(df[i,1] == "T", df[i,2] <- 1, ifelse(df[i,1] == "Bt",df[i,2] <- 1, df[i,2] <- 0))
  i = i+1
}


#we no longer need the rn column, b/c we are only concerned if the observations
#are Bt or T, and therefore 1 or 0
df<- df[,-1]


#set two classes: Bt or T, Other
df$class <- factor(df$class, levels=c(1,0),
                   labels=c("Bt-T", "Other"))

#randomly divide into 70% for training, 30% for validation
set.seed(1234)
train <- sample(nrow(df), 0.7*nrow(df))
#store each fragment of data 
df.train <- df[train,]
df.validate <- df[-train,]

#visualize summaries of benign/malignant for each data subset
table(df.train$class)
table(df.validate$class)

#################### Perform classification ########################
#Decision (classification) trees
library(rpart)

#grow the tree
set.seed(1234)
dtree <- rpart(class ~ ., data=df.train, method="class", parms=list(split="information"))
dtree$cptable

#visualize
print(dtree)
summary(dtree)

#plots the cross-validated error
plotcp(dtree)

library(rpart.plot)
#plotting the decision tree
#no need to prune after looking at cptable
prp(dtree, type = 2, extra = 104, fallen.leaves = TRUE, main="Decision Tree")

#use predictor for the validation subset
dtree.pred <- predict(dtree, df.validate, type="class")
#create and display the "confusion matrix"
dtree.perf <- table(df.validate$class, dtree.pred, dnn=c("Bt-T", "Other"))
dtree.perf

#Exercise: choose a one row sample and perform prediction
df.exercise.tree <- df[3,]
df.exercise.tree
#what is your prediction? What if you use the original tree?

#let's compute
ex.pred.tree <- predict(dtree, df.exercise, type="class")
ex.perf.tree <- table(df.exercise$class, ex.pred, dnn=c("Bt-T", "Other"))
ex.perf.tree

############################################################################
#KNN Classification
library(class)
knn.pred <- knn(df.train[,-1], 
                df.validate[,-1], df.train$class, k = 1)
knn.pred

#checking the accuracy
mean(knn.pred == df.validate$class)

#constructing a confusion matrix for the knn model
knn.perf <- table(knn.pred, df.validate$class)
knn.perf

#testing a one row sample
df.exercise.knn <- df[3,]
df.exercise.knn
ex.pred.knn <- knn(df.train[,-1], df.exercise[,-1], df.train$class, k = 1)
ex.perf.knn <- table(ex.pred, df.exercise$class)
ex.perf.knn


#################################################################################
#2. SVM classification with svm()
library(e1071)
set.seed(1234)

#hTrain= replace(df.train, is.na(df.train),FALSE)
#hValidate= replace(df.validate, is.na(df.validate),FALSE)
fit.svm <- svm(class~., data=df.train)
fit.svm
#perform classification of validation data
svm.pred <- predict(fit.svm, na.omit(df.validate))

#evaluate the predictive accuracy
#a cross-tabulation of actual status and predicted status 
#(called a confusion matrix)
svm.perf <- table(na.omit(df.validate)$class,
                  svm.pred, dnn=c("Actual", "Predicted"))
svm.perf


#Let's perform tuning
set.seed(1234)
tuned <- tune.svm(class~., data=df.train,
                  gamma=10^(-6:1),
                  cost=10^(-10:10))
tuned


fit.svm <- svm(class~., data=df.train, gamma=.01, cost=1)
svm.pred <- predict(fit.svm, na.omit(df.validate))
svm.perf <- table(na.omit(df.validate)$class,
                  svm.pred, dnn=c("Actual", "Predicted"))
svm.perf

df.exercise.svm <- df[3,]
names(df.exercise.svm) <- c("class","STMN1|200783_s_at" ,"CD9|201005_at" ,"HEXB|201944_at" ,"S100A13|202598_at" ,
                        "PLCG1|202789_at" ,"RYK|202853_s_at" ,"CTSS|202902_s_at" ,"OFD1|203569_s_at" ,
                        "RBMS1|203748_x_at" ,"LRMP|204674_at" ,"BTG3|205548_s_at" ,"LY86|205859_at" ,
                        "LILRB1|207104_x_at" ,"MAPK6|207121_s_at" ,"RBMS1|207266_x_at" ,"FLNB|208613_s_at" ,
                        "PTP4A2|208616_s_at" ,"PDLIM1|208690_s_at" ,"CAST|208908_s_at" ,"PRKD2|209282_at" ,
                        "POP7|209482_at" ,"UBE2E3|210024_s_at" ,"LILRB1|211336_x_at" ,"CAST|212586_at" ,
                        "MAPKBP1|213394_at" ,"N4BP3|214775_at" ,"AIF1|215051_x_at" ,"CD44|217523_at" ,
                        "NGRN|217722_s_at" ,"NA|218017_s_at" ,"SCPEP1|218217_at" ,"GLT25D1|218473_s_at" ,
                        "ATP13A2|218608_at" ,"NOLA1|219110_at" ,"VPREB1|221349_at" ,"LAT2|221581_s_at" ,
                        "LRMP|35974_at" ,"PRKD2|38269_at" ,"PRR14|45687_at" ,"LOC90379|91952_at" )

df.exercise.svm

ex.pred.svm <- predict(fit.svm, df.exercise.svm, type="response")
ex.resp.svm <- table(df.exercise.svm$class, ex.pred.svm, dnn=c("Actual", "Predicted"))
ex.resp.svm

#################################################################################
#Logistic Regression

fit.lr <- glm(class~., data=df.train, family=binomial())
summary(fit.lr)

prob <- predict(fit.lr, df.validate, type="response")
lr.pred <- factor(prob > .5, levels=c(FALSE, TRUE), labels=c("Bt-T", "other"))

lr.resp <- table(df.validate$class, lr.pred, dnn=c("Actual", "Predicted"))
lr.resp

df.exercise.lr <- df[3,]
df.exercise.lr

ex.prob.lr <- predict(fit.lr, df.exercise.lr, type="response")
ex.pred.lr <- factor(ex.prob.lr > .5, levels=c(FALSE, TRUE), labels=c("Bt-T", "other"))
ex.resp.lr <- table(df.exercise.lr$class, ex.pred.lr, dnn=c("Actual", "Predicted"))
ex.resp.lr
