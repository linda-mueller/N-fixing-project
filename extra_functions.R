#additional functions for N-fixing project

###########################
##### Data exploration ####
###########################

#!!! for data exploration, must cite:
#Mixed effects models and extensions in ecology with R. (2009).
#Zuur, AF, Ieno, EN, Walker, N, Saveliev, AA, and Smith, GM. Springer.

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



########## Spatial correlation #########

correlogram <- function(mod, dat, filename){  

  #calc residuals
  resids.calc <- resid(mod)
  sp.dat <- dat %>%
    mutate(resid = resids.calc)
  
  #generate sp object
  sp <- ncf::spline.correlog(x = as.numeric(sp.dat$latitude),
                             y = as.numeric(sp.dat$longitude),
                             z = as.numeric(sp.dat$resid),
                             xmax = 4000, resamp = 100, latlon=TRUE)
  #save figure
  png(filename, width = 10, height = 10, units = 'in', res = 300)
  plot(sp)
  dev.off()
  
}

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


###########################
#### Model validation #####
###########################

#GLM
#homogeneity: residuals v fitted (except for binomial)
#independence: residuals v covariates in/out of model
#influential observations: cooks distance
#overdispersion: >1 = overdispersed (except gamma, NB, GP, CMP; these built in correct for dispersion)

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


###########################
#### DHARMa validation ####
###########################

#testDispersion(mod) #test for over/under dispersion
#testZeroInflation(mod) #test zero inflation

#simulationOutput <- simulateResiduals(mod, plot = F) # calculate residuals (randomized quantile residuals)
#par(mfrow = c(2, 2)) # set panel arrangement
#plotResiduals(simulationOutput) #homogeneity (residuals v fitted)
#plotResiduals(simulationOutput, form = dat$var1) #independence continuous
#plotResiduals(simulationOutput, form = dat$var2) #independence categorical
#plotQQunif(simulationOutput) #qqplot

#R-square function

r2_glm <- function(response, formula, dat) {
  
  # Fit the full model with all variables
  full_model <- glm(formula = formula, data = dat, family = binomial(link ="logit"))
  
  # Fit the null model (intercept only)
  null_model <- glm(response ~ 1, data = dat, family = binomial(link ="logit"))
  
  # Get deviances
  deviance_full <- deviance(full_model)
  deviance_null <- deviance(null_model)
  
  # Calculate R square
  R2<- (deviance_null - deviance_full) / deviance_null
  
  return(R2)
}


#Partial R-square function

partial_r2_glm <- function(response, formula, dat, variable) {
  
  # Ensure the variable exists in the formula
  if (!variable %in% all.vars(formula)) {
    stop(paste("Variable", variable, "not found in the model formula."))
  }
  
  # Fit the full model with all variables
  full_model <- glm(formula = formula, data = dat, family = binomial(link ="logit"))
  
  # Create the reduced formula by removing the specified variable
  reduced_formula <- update(formula, paste(". ~ . -", variable))
  
  # Fit the reduced model without the specified variable
  reduced_model <- glm(formula = reduced_formula, data = dat, family = binomial(link ="logit"))
  
  # Get deviances
  deviance_full <- deviance(full_model)
  deviance_reduced <- deviance(reduced_model)

  # Calculate partial R2
  partial_R2 <- (deviance_reduced - deviance_full) / deviance_reduced

  return(partial_R2)
}


partial_r2_glm_list <- function(response, formula, dat, variable_list) {
  
  # Initialize an empty list to store the partial R^2 results
  partial_r2_list <- list()
  
  for (variable in variable_list) {
    # Ensure the variable exists in the formula
    if (!grepl(variable, paste(deparse(formula), collapse = " "))) {
      warning(paste("Variable or interaction term", variable, "not found in the model formula. Skipping."))
      next
    }
  
  # Fit the full model with all variables
  full_model <- glm(formula = formula, data = dat, family = binomial(link ="logit"))
  
  # Create the reduced formula by removing the specified variable
  reduced_formula <- update(formula, paste(". ~ . -", variable))
  
  # Fit the reduced model without the specified variable
  reduced_model <- glm(formula = reduced_formula, data = dat, family = binomial(link ="logit"))
  
  # Get deviances
  deviance_full <- deviance(full_model)
  deviance_reduced <- deviance(reduced_model)
  
  # Calculate partial R2
  partial_R2 <- (deviance_reduced - deviance_full) / deviance_reduced
  
  # Store the partial R2 in the list, using the variable name as the key
  partial_r2_list[[variable]] <- partial_R2
  }
  
  return(partial_r2_list)
}
