#R code from vignette

### An RRPP data frame
library(RRPP)
y <- matrix(rnorm(30), 10, 3)
x <- rnorm(10)
df <- data.frame(x = x, y = y)
df
rdf <- rrpp.data.frame(x = x, y = y)
rdf # looks more like a list

is.list(df)
is.list(rdf)

d <- dist(y) # distance matrix as data
rdf <- rrpp.data.frame(rdf, d = d) # works

#Example 1: Univariate ANOVA and Mixed-Models

data("PupfishHeads")
str(PupfishHeads)
PupfishHeads$logHeadSize <- log(PupfishHeads$headSize)
fit <- lm.rrpp(logHeadSize ~ sex + locality/year, 
               SS.type = "I", data = PupfishHeads, 
               print.progress = FALSE)
coef(fit)
summary(fit)
anova(fit, effect.type = "F") 

#Adjusting EMS for terms in model
anova(fit, effect.type = "F", 
  error = c("Residuals", "locality:year", "Residuals"))

#Model coefficients with bootstrap CI
coef(fit, test = TRUE)

#Visualize model predictions
sizeDF <- data.frame(sex = c("Female", "Male"))
rownames(sizeDF) <- c("Female", "Male")
sizePreds <- predict(fit, sizeDF)
plot(sizePreds)

plot(sizePreds, pch = 21, cex = 3, bg = c(2,4), lwd = 2)

#Toggling SS type I, II, III
fit2 <- lm.rrpp(logHeadSize ~ sex + locality/year, 
                SS.type = "II", data = PupfishHeads, print.progress = FALSE)
fit3 <- lm.rrpp(logHeadSize ~ sex + locality/year, 
                SS.type = "III", data = PupfishHeads, print.progress = FALSE)

anova(fit)
anova(fit2)
anova(fit3)

#Example 2: RRPP with High-Dimensional Data

data(Pupfish)
Pupfish$logSize <- log(Pupfish$CS) 
fit <- lm.rrpp(coords ~ logSize + Sex*Pop, SS.type = "I", 
               data = Pupfish, print.progress = FALSE) 
summary(fit, formula = FALSE)
anova(fit) 
coef(fit, test = TRUE) 

#Compare to parametric methods
fit$LM$data$coords <- Pupfish$coords
fit.par <- lm(fit$call$f1, data = fit$LM$data)
identical(fit$LM$coefficients, fit.par$coefficients)
summary(manova(fit.par))  


#Predicted values and ordinations
shapeDF <- expand.grid(Sex = levels(Pupfish$Sex), Pop = levels(Pupfish$Pop))
rownames(shapeDF) <- paste(shapeDF$Sex, shapeDF$Pop, sep = ".")
shapePreds <- predict(fit, shapeDF, confidence = 0.95)
plot(shapePreds, PC = TRUE, ellipse = TRUE) # generic 
plot(shapePreds, PC = TRUE, ellipse = TRUE, 
     pch = 19, col = 1:NROW(shapeDF)) # with added par arguments

groups <- interaction(Pupfish$Sex, Pupfish$Pop)
plot(fit, type = "PC") # generic
plot(fit, type = "PC", pch = 19, col = groups) # with added par arguments

#Regression plot of H-D data
plot(fit, type = "regression", reg.type = "PredLine", 
    predictor = Pupfish$logSize, pch=19,
    col = as.numeric(groups))

#Using a distance matrix as Y
D <- dist(Pupfish$coords) # inter-observation Euclidean distances
Pupfish$D <- D

fitD <- lm.rrpp(D ~ logSize + Sex*Pop, SS.type = "I", 
                data = Pupfish, print.progress = FALSE) 
anova(fitD)
anova(fit)

#Comparison to Vegan
fit.uni <- lm.rrpp(logSize ~ Sex*Pop, data=Pupfish, print.progress = FALSE)
anova(fit.uni)

library(vegan)
adonis(dist(Pupfish$logSize) ~ Pupfish$Sex*Pupfish$Pop)

fit.uni <- lm.rrpp(logSize ~ Sex*Pop, data=Pupfish,RRPP = FALSE, print.progress = FALSE)
anova(fit.uni)

#Example 3: Pairwise comparisons of groups, slopes, and dispersion
PWT <- pairwise(fit, groups = interaction(Pupfish$Sex, Pupfish$Pop))
summary(PWT, confidence = 0.95)

#Within-group dispersion analysis 
summary(PWT, confidence = 0.95, test.type = "var")

#Pairwise comparison of slopes
fit2 <- lm.rrpp(coords ~ logSize * Sex * Pop, SS.type = "I", 
                data = Pupfish, print.progress = FALSE, iter = 999) 
anova(fit, fit2, print.progress = FALSE)

PW2 <- pairwise(fit2, fit.null = fit, groups = groups, 
                covariate = Pupfish$logSize, print.progress = FALSE) 
PW2
summary(PW2, confidence = 0.95, 
        test.type = "dist") # distances between slope vector lengths
summary(PW2, confidence = 0.95, 
        test.type = "VC",
        angle.type = "deg") # correlation between slope vectors (and angles)

#Example 4: GLS Estimation
data(PlethMorph)
library(ape)
tree <- read.tree('plethtree.tre')
plot(tree)

fitOLS <- lm.rrpp(TailLength ~ SVL, 
                  data = PlethMorph,
                  print.progress = FALSE)
fitGLS <- lm.rrpp(TailLength ~ SVL, 
                  data = PlethMorph, 
                  Cov = PlethMorph$PhyCov,
                  print.progress = FALSE)

anova(fitOLS)
anova(fitGLS)

coef(fitOLS, test = TRUE)
coef(fitGLS, test = TRUE)

#Multivariate GLS
Y <- as.matrix(cbind(PlethMorph$TailLength,
PlethMorph$HeadLength,
PlethMorph$TailLength,
PlethMorph$Snout.eye,
PlethMorph$BodyWidth,
PlethMorph$Forelimb,
PlethMorph$Hindlimb))
PlethMorph <- rrpp.data.frame(PlethMorph, Y=Y)

fitOLSm <- lm.rrpp(Y ~ SVL, data = PlethMorph,
                   print.progress = FALSE)
fitGLSm <- lm.rrpp(Y ~ SVL, data = PlethMorph, 
                   Cov = PlethMorph$PhyCov,
                   print.progress = FALSE)

anova(fitOLSm)
anova(fitGLSm)

sizeDF <- data.frame(SVL = sort(PlethMorph$SVL))

plot(predict(fitOLSm, sizeDF), PC= TRUE) # Correlated error
plot(predict(fitGLSm, sizeDF), PC= TRUE) # Independent error

