require(TDA)

wise_train <- read.table('WISE_train.txt', head=TRUE, fill=TRUE)[,4:7]
wise_test1 <- read.table('WISE_test1.txt', head=TRUE, fill=TRUE)[,4:7]
wise_test2 <- read.table('WISE_test2.txt', head=TRUE, fill=TRUE)[,4:7]
wise_test3 <- read.table('WISE_test3.txt', head=TRUE, fill=TRUE)[,4:7]
wise_test4 <- read.table('WISE_test4.txt', head=TRUE, fill=TRUE)[,4:7]


index <- sample(1:1000, 150,replace=F)


Diagtrain <- ripsDiag(X = wise_train[index,], maxdimension = 2, maxscale = 5)
Diag1 <- ripsDiag(X = wise_test1[index,], maxdimension = 2, maxscale = 5)
Diag2 <- ripsDiag(X = wise_test2[index,], maxdimension = 2, maxscale = 5)
Diag3 <- ripsDiag(X = wise_test3[index,], maxdimension = 2, maxscale = 5)
Diag4 <- ripsDiag(X = wise_test4[index,], maxdimension = 2, maxscale = 5)



TD1 <- print(bottleneck(Diag1 = Diagtrain[["diagram"]], Diag2 = Diag1[["diagram"]],
                 dimension = 1))
TD2 <- print(bottleneck(Diag1 = Diagtrain[["diagram"]], Diag2 = Diag2[["diagram"]],
                        dimension = 1))
TD3 <- print(bottleneck(Diag1 = Diagtrain[["diagram"]], Diag2 = Diag3[["diagram"]],
                        dimension = 1))
TD4 <- print(bottleneck(Diag1 = Diagtrain[["diagram"]], Diag2 = Diag4[["diagram"]],
                        dimension = 1))




pdf("TDA.pdf",height=20,width = 15)
par(mfrow=c(5,4))

attach(wise_train)
plot(Diagtrain[["diagram"]],main="Model")
plot(w1mpro-w2mpro, w2mpro-w3mpro, pch=20, cex=0.5, main='Training set', col='#11111120') # infrared color-color diagram
plot(w2mpro-w3mpro, w3mpro-w4mpro, pch=20, cex=0.5, col='#11111120')  # infrared  color-color diagram
plot(w2mpro, w2mpro-w3mpro, pch=20, cex=0.2, col='#11111120')  # infrared color-magnitude diagram (with truncation)

attach(wise_test1)
plot(Diag1[["diagram"]],main=paste("Bottleneck distance to model = ",round(TD1,3)))
plot(w1mpro-w2mpro, w2mpro-w3mpro, pch=20, cex=0.5, main='Test set 1', col='#11111120') # infrared color-color diagram
plot(w2mpro-w3mpro, w3mpro-w4mpro, pch=20, cex=0.5, col='#11111120')  # infrared  color-color diagram
plot(w2mpro, w2mpro-w3mpro, pch=20, cex=0.5, col='#11111120')  # infrared color-magnitude diagram (with truncation)

attach(wise_test2)
plot(Diag2[["diagram"]],main=paste("Bottleneck distance to model = ",round(TD2,3)))
plot(w1mpro-w2mpro, w2mpro-w3mpro, pch=20, cex=0.5, main='Test set 2', col='#11111120') # infrared color-color diagram
plot(w2mpro-w3mpro, w3mpro-w4mpro, pch=20, cex=0.5, col='#11111120')  # infrared  color-color diagram
plot(w2mpro, w2mpro-w3mpro, pch=20, cex=0.5, col='#11111120')  # infrared color-magnitude diagram (with truncation)

attach(wise_test3)
plot(Diag3[["diagram"]],main=paste("Bottleneck distance to model = ",round(TD3,3)))
plot(w1mpro-w2mpro, w2mpro-w3mpro, pch=20, cex=0.5, main='Test set 3', col='#11111120') # infrared color-color diagram
plot(w2mpro-w3mpro, w3mpro-w4mpro, pch=20, cex=0.52, col='#11111120')  # infrared  color-color diagram
plot(w2mpro, w2mpro-w3mpro, pch=20, cex=0.5, col='#11111120')  # infrared color-magnitude diagram (with truncation)


attach(wise_test4)
plot(Diag4[["diagram"]],main=paste("Bottleneck distance to model = ",round(TD4,3)))

plot(w1mpro-w2mpro, w2mpro-w3mpro, pch=20, cex=0.5, main='Test set 4', col='#11111120') # infrared color-color diagram
plot(w2mpro-w3mpro, w3mpro-w4mpro, pch=20, cex=0.5, col='#11111120')  # infrared  color-color diagram
plot(w2mpro, w2mpro-w3mpro, pch=20, cex=0.5, col='#11111120')  # infrared color-magnitude diagram (with truncation)

dev.off()









 