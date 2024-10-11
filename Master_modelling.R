###########################
########## TO DO ##########
###########################

#variogram
#GAM section expand

###########################
##### Data exploration ####
###########################

#!!! for data exploration, must cite:
#Mixed effects models and extensions in ecology with R. (2009).
#Zuur, AF, Ieno, EN, Walker, N, Saveliev, AA, and Smith, GM. Springer.

#cont.var <- c("var1","var2")
#Outliers
Mydotplot <- function(DataSelected){
  P <- dotplot(as.matrix(as.matrix(DataSelected)),
               groups=FALSE,
               strip = strip.custom(bg = 'white',
                                    par.strip.text = list(cex = 1.2)),
               scales = list(x = list(relation = "free", draw = TRUE),
                             y = list(relation = "free", draw = FALSE)),
               col=1, cex  = 0.5, pch = 16,
               xlab = list(label = "Value of the variable", cex = 1.5),
               ylab = list(label = "Order of the data from text file", cex = 1.5))
  print(P)  
}
#Mydotplot(dat[,cont.var])

#Collinearity
#Continuous variables
Mypairs <- function(Z) {
  MyVarx <- colnames(Z)
  pairs(Z, labels = MyVarx,
        cex.labels =  2,
        lower.panel = function(x, y, digits=2, prefix="", cex.cor = 7) {
          panel.cor(x, y, digits, prefix, cex.cor)}, 
        upper.panel =  function(x, y) points(x, y, 
                                             pch = 16, cex = 0.8, 
                                             col = gray(0.1)))
  #print(P)
}

panel.cor <- function(x, y, digits=1, prefix="", cex.cor = 6)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r1=cor(x,y,use="pairwise.complete.obs")
  r <- abs(cor(x, y,use="pairwise.complete.obs"))
  txt <- format(c(r1, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) { cex <- 0.9/strwidth(txt) } else {
    cex = cex.cor}
  text(0.5, 0.5, txt, cex = cex * r)
}

#Mypairs(dat[, cont.var])

#Multipanel boxplots
#Continuous v categorical variables
Mybwplot <- function(Z, MyVar, TargetVar){
  AllY <- as.vector(as.matrix(Z[,MyVar]))
  AllX <- rep(Z[,TargetVar], length(MyVar))
  ID <- rep(MyVar, each = nrow(Z))
  P <- bwplot(AllY ~ factor(AllX) | ID, horizontal = FALSE,
              ylab = "", xlab = "",
              scales = list(alternating = TRUE,cex.lab = 1.5,
                            x = list(relation = "same",rot =90, abbreviate = TRUE, cex = 1.5),
                            y = list(relation = "free", draw = FALSE)),
              strip = strip.custom(bg = 'white',
                                   par.strip.text = list(cex = 1.2)),
              cex = .5,
              par.settings = list(
                box.rectangle = list(col = 1),
                box.umbrella  = list(col = 1),
                plot.symbol   = list(cex = .5, col = 1)))
  print(P)
}
#Mybwplot(dat, cont.var, "Sex")

#Multipanel scatterplots
Myxyplot <- function(Z, MyV, NameY1, MyXlab = "", MyYlab="") {
  AllX  <- as.vector(as.matrix(Z[,MyV]))
  AllY  <- rep(Z[,NameY1] , length(MyV))
  AllID <- rep(MyV, each = nrow(Z))
  library(mgcv)
  library(lattice)
  P <- xyplot(AllY ~ AllX|factor(AllID), col = 1,
              xlab = list(MyXlab, cex = 1.5),
              #ylab = list("Response variable", cex = 1.5),
              #ylab = list("Pearson residuals", cex = 1.5),
              ylab = list(MyYlab, cex = 1.5),
              #layout = c(2,2),   #Modify
              strip = function(bg='white', ...)
                strip.default(bg='white', ...),
              scales = list(alternating = TRUE,
                            x = list(relation = "free"),
                            y = list(relation = "same")),
              panel=function(x, y){
                panel.grid(h=-1, v= 2)
                panel.points(x, y, col = 1)
                panel.loess(x, y, span = 0.8,col = 1, lwd = 2)
              })
  print(P)
}
#Myxyplot(dat, var, "Wingcrd")

#Multipanel scatterplots
MyMultipanel.ggp2 <- function(Z, varx, vary, 
                              ylab = "Response variable",
                              addSmoother = FALSE,
                              addRegressionLine = FALSE,
                              addHorizontalLine = FALSE) {
  K <- length(varx)
  MyData <- data.frame(Y = rep(as.vector(as.matrix(Z[,vary])), K),
                       X = as.vector(as.matrix(Z[, varx])),
                       Var = rep(varx, each = nrow(Z))) 
  library(ggplot2)
  p <- ggplot(MyData, aes(y = Y, x = X))
  p <- p + geom_point(size = 0.5) + ylab(ylab) + xlab("Covariates")
  p <- p + theme(text = element_text(size=15))
  if (addSmoother == TRUE) {
    p <- p + geom_smooth(se = TRUE, col = "red", lwd = 1)
  }
  if (addRegressionLine == TRUE) {
    p <- p + geom_smooth(se = TRUE, col = "red", lwd = 1, method = "lm")
  }
  if (addRegressionLine == TRUE) {
    p <- p + geom_smooth(se = TRUE, col = "red", lwd = 1, method = "lm")
  }
  if (addHorizontalLine == TRUE) {
    p <- p + geom_hline(yintercept = 0)
  }
  p <- p + facet_wrap(~ Var, scales = "free_x")
  suppressMessages(print(p)) 	
}
#MyMultipanel.ggp2(dat, cont.var, vary = "vary", ylab = "label")

###########################
#### Model validation #####
###########################

#LM
#homogeneity: residuals v fitted
#independence: residuals v covariates in/out of model
#normality: of residuals
#influential observations: cooks distance

#GLM
#homogeneity: residuals v fitted (except for binomial)
#independence: residuals v covariates in/out of model
#influential observations: cooks distance
#overdispersion: >1 = overdispersed (except gamma, NB, GP, CMP; these built in correct for dispersion)
#spatial/temporal autocorrelation: variogram

#GAM
#homogeneity: residuals v fitted (except for binomial)
#independence: residuals v covariates in/out of model
#normality: of residuals
#??? overdispersion: unsure
#spatial/temporal autocorrelation: variogram

#homogeneity (residuals v fitted) function
resid_fit_plot <- function(mod,dat) {
  E2 <- resid(mod, type = "pearson")
  F2 <- fitted(mod)
  plot(x = F2, 
       y = E2, 
       xlab = "Fitted values",
       ylab = "Residuals")
  abline(h = 0) 
}
#resid_fit_plot(mod,dat)

#independence function continuous
indep_cont_plot <- function(mod,dat,var) {
  E2 <- resid(mod, type = "pearson")
  plot(x = var, 
       y = E2, 
       xlab = "var",
       ylab = "Pearson residuals", 
       pch = 16)
  abline(h = 0) 
}
#indep_cont_plot(mod,dat,df$var)

#independence function categorical
indep_cat_plot <- function(mod,dat,var) {
  E2 <- resid(mod, type = "pearson")
  boxplot(E2 ~ var, 
          data = dat, 
          cex.lab = 1.5,
          xlab = "var",
          ylab = "Pearson residuals")
  abline(h = 0) 
}
#indep_cat_plot(mod,dat,df$var)

#cooks distance function
cooks_dist_plot <- function(mod) {
  plot(cooks.distance(mod), 
       type = "h",
       ylim = c(0,1),
       xlab = "Observations",
       ylab = "Cook distance values")
  abline(h = 0, lwd = 2, lty = 2) 
}
#cooks_dist_plot(mod)

#normality function
resid_hist<- function(mod) {
  E2 <- resid(mod, type = "pearson")
  hist(E2, main = "", 
       xlab = "Residuals")
}
#resid_hist(mod)

#qqplots
resid_qqplot <- function(mod) {
  E2 <- rstandard(mod)
  qqnorm(E2, main = "")
  qqline(E2)
}
#resid_qqplot(mod)

#overdispersion function (poisson; var = mean)
disp_check <- function(mod,dat) {
  N <- nrow(dat)
  p <- length(coef(mod))
  E1 <- resid(mod, type = "pearson")
  Dispersion <- sum(E1^2)/ (N-p)
  return(Dispersion)
}
#disp_check(mod,dat)

#???variogram

###########################
#### DHARMa validation ####
###########################

#set.seed(556)
#testDispersion(mod) #test for over/under dispersion
#testZeroInflation(mod) #test zero inflation

#simulationOutput <- simulateResiduals(mod, plot = F) # calculate residuals (randomized quantile residuals)
#par(mfrow = c(2, 2)) # set panel arrangement
#plotResiduals(simulationOutput) #homogeneity (residuals v fitted)
#plotResiduals(simulationOutput, form = dat$var1) #independence continuous
#plotResiduals(simulationOutput, form = dat$var2) #independence categorical
#plotQQunif(simulationOutput) #qqplot

###########################
##### Model selection #####
###########################

#step(mod) #automated backward selection
#drop1(mod) #test whether vars matter

###########################
#### Other useful code ####
###########################

#GLM generalised R^2: 
#(null deviance - residual deviance)/ null deviance

#AIC(M1, M2) #compare multiple models
#anova(M1, M2, test = "F")

#100 * sum(dat$var== 0) / nrow(dat) #count zeros in datase

Spat.cor <- function(mod,dat,dist) {
  coords <- cbind(dat$longitude, dat$latitude)
  matrix.dist = as.matrix(dist(cbind(dat$longitude, dat$latitude)))
  matrix.dist[1:10, 1:10]
  matrix.dist.inv <- 1/matrix.dist
  matrix.dist.inv[1:10, 1:10]
  diag(matrix.dist.inv) <- 0
  matrix.dist.inv[1:10, 1:10]
  myDist = dist
  # calculate residuals autocovariate (RAC)
  rac <- autocov_dist(resid(mod), coords, nbs = myDist, type = "inverse", zero.policy = TRUE, style = "W", longlat=T)
  return(rac)
}
#Spat.cor(mod,dat,2000)

Spat.cor.rep <- function(mod,dat,dist) {
  coords <- cbind(dat$longitude, dat$latitude)+matrix(runif(2*nrow(dat), 0, 0.00001), nrow = nrow(dat), ncol = 2)
  matrix.dist = as.matrix(dist(cbind(dat$longitude, dat$latitude)))
  matrix.dist[1:10, 1:10]
  matrix.dist.inv <- 1/matrix.dist
  matrix.dist.inv[1:10, 1:10]
  diag(matrix.dist.inv) <- 0
  matrix.dist.inv[1:10, 1:10]
  myDist = dist
  # calculate residuals autocovariate (RAC)
  rac <- autocov_dist(resid(mod), coords, nbs = myDist, type = "inverse", zero.policy = TRUE, style = "W", longlat=T)
  return(rac)
}
#Spat.cor.rep(mod,dat,2000)