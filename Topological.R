require(TDA)

wise_train <- read.table('WISE_train.txt', head=TRUE, fill=TRUE)[,4:7]
wise_test1 <- read.table('WISE_test1.txt', head=TRUE, fill=TRUE)[,4:7]
wise_test4 <- read.table('WISE_test4.txt', head=TRUE, fill=TRUE)[,4:7]


index <- sample(1:1000, 100,replace=F)


Diag1 <- ripsDiag(X = wise_train[index,], maxdimension = 2, maxscale = 5)
plot(Diag1[["diagram"]])

Diag2 <- ripsDiag(X = wise_test1[index,], maxdimension = 2, maxscale = 5)
plot(Diag2[["diagram"]])

Diag4 <- ripsDiag(X = wise_test4[index,], maxdimension = 2, maxscale = 5)
plot(Diag4[["diagram"]])

print(bottleneck(Diag1 = Diag1[["diagram"]], Diag2 = Diag4[["diagram"]],
                 dimension = 1))


print(wasserstein(Diag1 = Diag1[["diagram"]], Diag2 = Diag4[["diagram"]],
                  p = 2, dimension = 1))



Diag <- matrix(c(0, 0, 10, 1, 0, 3, 1, 3, 8), ncol = 3, byrow = TRUE)
DiagLim <- 10
colnames(Diag) <- c("dimension", "Birth", "Death")

#persistence landscape
tseq <- seq(0,DiagLim, length = 1000)
Land <- landscape(Diag, dimension = 1, KK = 1, tseq)

par(mfrow = c(1,2))
plot.diagram(Diag)
plot(tseq, Land, type = "l", xlab = "t", ylab = "landscape", asp = 1)


 