require(randomForest)
data(iris)
iris.rf <- randomForest(iris[,1:4], iris$Species)
print(iris.rf)

print(iris.rf$oob.times[1:10]) # oob = out of bag

sous_echant=c(25,75,135)
iris$Species[sous_echant]

iris.rf$oob.times[sous_echant]
m=margin(iris.rf)
print(m[sous_echant])     
print(iris.rf$predicted[1:10])
table(iris$Species, iris.rf$predicted)
iris.rf$confusion
error_rate=1-sum(diag(iris.rf$confusion))/sum(iris.rf$confusion)
print(error_rate)

library(readr)
Species_to_predict <- read.table(here::here("data","Species_to_predict.csv"),
                   sep=";", header=T, dec=".")

predicted=predict(iris.rf,newdata=Species_to_predict)
print(predicted[1:5])

print(iris.rf$importance)

iris.rf=randomForest(iris[,1:4], iris$Species, ntree=5000) 
plot(iris.rf$err.rate[,1], type="l")
