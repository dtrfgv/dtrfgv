data(iris)
names(iris)
table(iris$Species)
data<-iris[,c(5,1:4)]
names(data)
names(data)[1]<-"Y"
data$Y<-c(rep(1,75),rep(0,75))
group<-c(NA,1,1,2,2)
names(data)[-1]<-paste(names(data)[-1],"_G",group[-1],sep="")
names(data)

ind0 <- which(data$Y == "0")
ind1 <- which(data$Y == "1")
ind0_val<-sample(ind0, floor(length(ind0) * 1/3), F)
ind1_val<-sample(ind1, floor(length(ind1) * 1/3), F)
validation<-data[c(ind0_val,ind1_val),]
train<-data[setdiff(c(ind0,ind1),c(ind0_val,ind1_val)),]
table(validation$Y)
dim(validation)
table(train$Y)
dim(train)

#maximal_tree <-Tree_PLDA(train,group=group,grp.importance=TRUE)
#print(maximal_tree$tree)