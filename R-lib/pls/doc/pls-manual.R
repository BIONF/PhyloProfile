### R code from vignette source 'pls-manual.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: pls-manual.Rnw:67-69
###################################################
pdf.options(pointsize=10)
options(digits = 4)


###################################################
### code chunk number 2: pls-manual.Rnw:304-305
###################################################
library(pls)


###################################################
### code chunk number 3: pls-manual.Rnw:332-335
###################################################
data(yarn)
data(oliveoil)
data(gasoline)


###################################################
### code chunk number 4: pls-manual.Rnw:343-349
###################################################
par(mar = c(2, 4, 0, 1) + 0.1)
matplot(t(gasoline$NIR), type = "l", lty = 1, ylab = "log(1/R)", xaxt = "n")
ind <- pretty(seq(from = 900, to = 1700, by = 2))
ind <- ind[ind >= 900 & ind <= 1700]
ind <- (ind - 898) / 2
axis(1, ind, colnames(gasoline$NIR)[ind])


###################################################
### code chunk number 5: pls-manual.Rnw:358-360
###################################################
gasTrain <- gasoline[1:50,]
gasTest <- gasoline[51:60,]


###################################################
### code chunk number 6: pls-manual.Rnw:363-364
###################################################
gas1 <- plsr(octane ~ NIR, ncomp = 10, data = gasTrain, validation = "LOO")


###################################################
### code chunk number 7: pls-manual.Rnw:369-370
###################################################
summary(gas1)


###################################################
### code chunk number 8: pls-manual.Rnw:379-380 (eval = FALSE)
###################################################
## plot(RMSEP(gas1), legendpos = "topright")


###################################################
### code chunk number 9: pls-manual.Rnw:385-387
###################################################
par(mar = c(4, 4, 2.5, 1) + 0.1)
plot(RMSEP(gas1), legendpos = "topright")


###################################################
### code chunk number 10: pls-manual.Rnw:406-407 (eval = FALSE)
###################################################
## plot(gas1, ncomp = 2, asp = 1, line = TRUE)


###################################################
### code chunk number 11: pls-manual.Rnw:413-415
###################################################
par(mar = c(4, 4, 2.5, 1) + 0.1)
plot(gas1, ncomp = 2, asp = 1, line = TRUE)


###################################################
### code chunk number 12: pls-manual.Rnw:429-430 (eval = FALSE)
###################################################
## plot(gas1, plottype = "scores", comps = 1:3)


###################################################
### code chunk number 13: pls-manual.Rnw:435-436
###################################################
plot(gas1, plottype = "scores", comps = 1:3)


###################################################
### code chunk number 14: pls-manual.Rnw:444-448
###################################################
par(mar = c(4, 4, 0.3, 1) + 0.1)
plot(gas1, "loadings", comps = 1:2, legendpos = "topleft",
     labels = "numbers", xlab = "nm")
abline(h = 0)


###################################################
### code chunk number 15: pls-manual.Rnw:462-463
###################################################
explvar(gas1)


###################################################
### code chunk number 16: pls-manual.Rnw:469-472 (eval = FALSE)
###################################################
## plot(gas1, "loadings", comps = 1:2, legendpos = "topleft",
##      labels = "numbers", xlab = "nm")
## abline(h = 0)


###################################################
### code chunk number 17: pls-manual.Rnw:481-482
###################################################
predict(gas1, ncomp = 2, newdata = gasTest)


###################################################
### code chunk number 18: pls-manual.Rnw:486-487
###################################################
RMSEP(gas1, newdata = gasTest)


###################################################
### code chunk number 19: pls-manual.Rnw:587-588
###################################################
dens1 <- plsr(density ~ NIR, ncomp = 5, data = yarn)


###################################################
### code chunk number 20: pls-manual.Rnw:592-594
###################################################
dim(oliveoil$sensory)
plsr(sensory ~ chemical, data = oliveoil)


###################################################
### code chunk number 21: pls-manual.Rnw:615-617
###################################################
trainind <- which(yarn$train == TRUE)
dens2 <- update(dens1, subset = trainind)


###################################################
### code chunk number 22: pls-manual.Rnw:621-622
###################################################
dens3 <- update(dens1, ncomp = 10)


###################################################
### code chunk number 23: pls-manual.Rnw:652-653
###################################################
olive1 <- plsr(sensory ~ chemical, scale = TRUE, data = oliveoil)


###################################################
### code chunk number 24: pls-manual.Rnw:658-659
###################################################
gas2 <- plsr(octane ~ msc(NIR), ncomp = 10, data = gasTrain)


###################################################
### code chunk number 25: pls-manual.Rnw:664-665 (eval = FALSE)
###################################################
## predict(gas2, ncomp = 3, newdata = gasTest)


###################################################
### code chunk number 26: pls-manual.Rnw:717-721 (eval = FALSE)
###################################################
## ncomp.onesigma <- selectNcomp(gas2, method = "onesigma", plot = TRUE,
##                               ylim = c(.18, .6))
## ncomp.permut <- selectNcomp(gas2, method = "randomization", plot = TRUE,
##                             ylim = c(.18, .6))


###################################################
### code chunk number 27: pls-manual.Rnw:735-740
###################################################
par(mfrow = c(1,2))
ncomp.onesigma <- selectNcomp(gas1, "onesigma", plot = TRUE,
                              ylim = c(.18, .6))
ncomp.permut <- selectNcomp(gas1, "randomization", plot = TRUE,
                            ylim = c(.18, .6))


###################################################
### code chunk number 28: pls-manual.Rnw:760-763
###################################################
gas2.cv <- crossval(gas2, segments = 10)
plot(MSEP(gas2.cv), legendpos="topright")
summary(gas2.cv, what = "validation")


###################################################
### code chunk number 29: pls-manual.Rnw:818-820 (eval = FALSE)
###################################################
## plot(gas1, plottype = "coef", ncomp=1:3, legendpos = "bottomleft",
##      labels = "numbers", xlab = "nm")


###################################################
### code chunk number 30: pls-manual.Rnw:824-827
###################################################
par(mar = c(4, 4, 2.5, 1) + 0.1)
plot(gas1, plottype = "coef", ncomp=1:3, legendpos = "bottomleft",
     labels = "numbers", xlab = "nm")


###################################################
### code chunk number 31: pls-manual.Rnw:859-861
###################################################
par(mar = c(4, 4, 0, 1) + 0.1)
plot(gas1, plottype = "correlation")


###################################################
### code chunk number 32: pls-manual.Rnw:930-931
###################################################
predict(gas1, ncomp = 2:3, newdata = gasTest[1:5,])


###################################################
### code chunk number 33: pls-manual.Rnw:943-944
###################################################
predict(gas1, comps = 2, newdata = gasTest[1:5,])


###################################################
### code chunk number 34: pls-manual.Rnw:961-962
###################################################
drop(predict(gas1, ncomp = 2:3, newdata = gasTest[1:5,]))


###################################################
### code chunk number 35: pls-manual.Rnw:1000-1001 (eval = FALSE)
###################################################
## predplot(gas1, ncomp = 2, newdata = gasTest, asp = 1, line = TRUE)


###################################################
### code chunk number 36: pls-manual.Rnw:1007-1009
###################################################
par(mar = c(4, 4, 2.5, 1))
predplot(gas1, ncomp = 2, newdata = gasTest, asp = 1, line = TRUE)


###################################################
### code chunk number 37: pls-manual.Rnw:1049-1050
###################################################
pls.options()


###################################################
### code chunk number 38: pls-manual.Rnw:1056-1057
###################################################
pls.options(plsralg = "oscorespls")


###################################################
### code chunk number 39: pls-manual.Rnw:1210-1219
###################################################
X <- gasTrain$NIR
Y <- gasTrain$octane
ncomp <- 5
cvPreds <- matrix(nrow = nrow(X), ncol = ncomp)
for (i in 1:nrow(X)) {
    fit <- simpls.fit(X[-i,], Y[-i], ncomp = ncomp, stripped = TRUE)
    cvPreds[i,] <- (X[i,] - fit$Xmeans) %*% drop(fit$coefficients) +
        fit$Ymeans
}


###################################################
### code chunk number 40: pls-manual.Rnw:1222-1223
###################################################
sqrt(colMeans((cvPreds - Y)^2))


