library(gllvm)

setwd("/home/bertv/github/GLLVM-workshop/data")
data(Skabbholmen, package  = "gllvm")
Y <- Skabbholmen$Y
X <- Skabbholmen$X

# We should standardise covariates in advance of model fitting
# This is always good practice when fitting models with numerical optimisation
X$Elevation.s <- c(scale(X$Elevation))
X$Year.s <- X$Year-min(X$Year)
# Omit species with few observations
# Technically, only needed in trade with complexity of the model (number of parameters per species)
# But we will do it here in advance
Y <- Y[,colSums(Y>0)>4]

# Place a frequent species in front
# This is not strictly necessary
# But can make model fitting a bit more painless
idx <- order(colSums(Y>0))[rep(mean(1:ncol(Y)),2)+c(0,1)]
Y <- Y[,c(idx,(1:ncol(Y))[-idx])]

TMB::openmp(parallel::detectCores()-1,autopar=TRUE, DLL = "gllvm")

# We fit a basic JSDM
model<-gllvm(Y,X,formula=~Elevation.s,num.lv=2,family="ordinal", seed = 123)
# alternative link functions for binomial models
# model<-gllvm(Y,X,formula=~ELEV.s,num.lv=2,family=binomial(link="probit"), seed = 2)
# model<-gllvm(Y,X,formula=~ELEV.s,num.lv=2,family=binomial(link="logit"), seed = 2)

coefplot(model)

SR <- predictSR(model)
plot(SR)

# partial prediction plot
# Draw grid for prediction
newELEV.s <- seq(min(X$Elevation.s), max(X$Elevation.s), length.out = 100)
newx <- newELEV.s * sd(X$Elevation) + mean(X$Elevation)
newLV <- matrix(0, ncol = 2, nrow = 100)

# Predict SR
SRpred <- predictSR(model, newX = data.frame(Elevation.s = newELEV.s),
                    newLV = newLV)
# Create plot
plot(y=rowSums(Y>0),X$Elevation, xlab = "Elevation", ylab = "Prediction", col = "gray", pch =  20)
# Plot prediction and intervals
lines(y = SRpred$expected$fit, x = newx, col = "red")
lines(y = SRpred$expected$lower, x = newx, col = "red", lty = "dashed")
lines(y = SRpred$expected$upper, x = newx, col = "red", lty = "dashed")

# Here, I have set the seed after running the model with seed = 1, n.init  = 10 to get a stable fit
model2<-gllvm(Y,X,formula=~Elevation.s+I(Elevation.s^2),num.lv=2,family="ordinal", seed = 8004,trace=TRUE)

# Predict SR
SRpred2 <- predictSR(model2, newX = data.frame(Elevation.s = newELEV.s),
                    newLV = newLV)
# Create plot
plot(y=rowSums(Y>0),X$Elevation, xlab = "Elevation", ylab = "Prediction", col = "gray", pch =  20)

# Plot prediction and intervals
lines(y = SRpred2$expected$fit, x = newx, col = "red")
lines(y = SRpred2$expected$lower, x = newx, col = "red", lty = "dashed")
lines(y = SRpred2$expected$upper, x = newx, col = "red", lty = "dashed")

# we could now fit the following model
# model3a <- gllvm(Y,X,formula=~(Elevation.s+I(Elevation.s^2))*Year.s,num.lv=2,family="ordinal", seed = 8004,trace=TRUE)
# but this has 8 parameters per species, plus more, and we seleted species by 4 observations or more
# So we are hitting our limit 
# here we do a constrained ordination instead
# it reduces the parameters relative to the JSDM, so we can easily use a lot more predictors
# we go from 1 parameter per covariate per species to `num.RR` parameters per covariate and per species

model3 <- gllvm(Y,X,lv.formula=~(Elevation.s+I(Elevation.s^2))*Year.s,num.RR=2,family="ordinal", seed = 1,trace=TRUE)

# Predict SR
SRpred3a <- predictSR(model3, newX = data.frame(Elevation.s = newELEV.s, Year.s = 0),# 1978
                     newLV = newLV)
SRpred3b <- predictSR(model3, newX = data.frame(Elevation.s = newELEV.s, Year.s = 6),# 1984
                     newLV = newLV)

# AUC
goodnessOfFit(model2, measure = "AUC");goodnessOfFit(model3,measure  ="AUC") # model 2 seems better
# but beware! We need to check without the unconstrained LVs!
goodnessOfFit(model2,pred=predict(model2,newLV=matrix(0,ncol=2,nrow=nrow(Y)), ordinal.cat = 1, type = "response"),measure="AUC")
goodnessOfFit(model3,measure="AUC")

# so that we conclude: constrained ordination captures the SR-covariate gradient slightly better.

# Compare with a naive  richness model
SRnaive <- rowSums(Y>0)
SRnaiveMod <- glm(SRnaive~(Elevation.s+I(Elevation.s^2))*Year.s,data=X,family="poisson")

predSRnaive1 <- predict(SRnaiveMod, newdata = data.frame(Elevation.s = newELEV.s, Year.s = 0), se.fit = TRUE)# 1978
predSRnaive2 <- predict(SRnaiveMod, newdata = data.frame(Elevation.s = newELEV.s, Year.s = 6), se.fit = TRUE)# 1984

par(mfrow = c(1, 2))

# Panel 1: year-specific richness predictions
plot(y = rowSums(Y > 0), X$Elevation, xlab = "Elevation", ylab = "Prediction",
     col = 3 + X$Year.s, pch = 20 + X$Year.s - 1, main = "Richness  (gllvm)")

polygon(c(newx, rev(newx)), c(SRpred3a$expected$lower, rev(SRpred3a$expected$upper)),
        col = adjustcolor(3, alpha.f = 0.25), border = NA)
polygon(c(newx, rev(newx)), c(SRpred3b$expected$lower, rev(SRpred3b$expected$upper)),
        col = adjustcolor(9, alpha.f = 0.25), border = NA)

lines(y = SRpred3a$expected$fit, x = newx, col = 3, lwd = 2)
lines(y = SRpred3b$expected$fit, x = newx, col = 9, lwd = 2)

points(y = rowSums(Y > 0), X$Elevation, col = 3 + X$Year.s, pch = 20 + X$Year.s - 1)

# Panel 2: naive richness, by year
plot(y = rowSums(Y > 0), X$Elevation, xlab = "Elevation", ylab = "Prediction",
     col = 3 + X$Year.s, pch = 20 + X$Year.s - 1, main = "Richness (naive)")

polygon(c(newx, rev(newx)),
        c(exp(predSRnaive1$fit - 1.96*predSRnaive1$se.fit), rev(exp(predSRnaive1$fit + 1.96*predSRnaive1$se.fit))),
        col = adjustcolor(3, alpha.f = 0.25), border = NA)
polygon(c(newx, rev(newx)),
        c(exp(predSRnaive2$fit - 1.96*predSRnaive2$se.fit), rev(exp(predSRnaive2$fit + 1.96*predSRnaive2$se.fit))),
        col = adjustcolor(9, alpha.f = 0.25), border = NA)

lines(y = exp(predSRnaive1$fit), x = newx, col = 3, lwd = 2)
lines(y = exp(predSRnaive2$fit), x = newx, col = 9, lwd = 2)

points(y = rowSums(Y > 0), X$Elevation, col = 3 + X$Year.s, pch = 20 + X$Year.s - 1)
