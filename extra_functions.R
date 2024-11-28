#additional functions for analysis of N-fixing project

###########################
##### Data exploration ####
###########################

#for data exploration, must cite:
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


