#N-fixing project: Script version November 2024

####Libraries and initial dataset####
# Load the necessary libraries
library(dplyr);library(segmented);library(nlme);library(mgcv);library(gridExtra);library(betareg);
library(MASS);library(lme4);library(lmerTest);library(lsmeans);library(ggeffects);library(spdep);
library(ggplot2);library(effects);library(ncf);library(ape);library(sjPlot);library(MuMIn);
library(tidyverse);library(car); library(V8); library(DHARMa); library(broom);library(grateful)

source("extra_functions.R")    #load script with additional functions
options(na.action = "na.fail") #Change na. action

#Load dataset from data cleaning script
gdat.ref <- readRDS("data/fulldata_for_analysis_2024.rds") %>%
  dplyr::select(c("entity_ID","entity_class","nfix","nfixno", "latitude","abs.lat","longitude",                          # Select columns
                  "landtype","status","presence","area","elev_range","precipitation","dist", "temperature","urbanland")) 

#### M1 Broad Presence (across islands and mainlands including native and naturalized floras)####

#dataset for broad presence model including mainlands
gdat.ml.pres <- gdat.ref %>%
  dplyr::select(c("entity_ID","nfix","nfixno", "latitude","abs.lat","longitude",                               # Select columns
                  "landtype","status","presence","area","elev_range","precipitation","temperature"))%>%
  filter(!landtype=="other_island")%>%                                                                         #remove islands with undetermined type
  mutate(area = as.vector(log10((area)+.01)))%>%                                                               #log transform area
  drop_na()%>%
  mutate_at(c("abs.lat","area","elev_range","precipitation", "temperature"), scale)                            #scale all explanatory variables

#check correlations
dat <- gdat.ml.pres
names(gdat.ml.pres)
cont.var <- c("abs.lat", "elev_range","area","precipitation","temperature")
Mypairs(dat[,cont.var]) # area and elevation range 0.56, abs.lat and temperature -0.88


#### M1 Broad Presence (across islands and mainlands including native and naturalized floras)
broad.pres.model.full<- glm(presence~landtype*status+abs.lat+area+elev_range+precipitation+temperature, data =gdat.ml.pres, family= binomial(link ="logit"))
summary(broad.pres.model.full)

#model including spatial correlation variable
rac <- Spat.cor.rep(broad.pres.model.full,gdat.ml.pres,2000)
broad.pres.model.rac.full <- glm(presence~landtype*status+abs.lat+area+elev_range+precipitation+temperature+rac, data = gdat.ml.pres, family = binomial(link ="logit"))
summary(broad.pres.model.rac.full)

#M1 with only selected variables after stepwise regression
broad.pres.model<- glm(presence~landtype*status+area+elev_range+precipitation+temperature, data =gdat.ml.pres, family= binomial(link ="logit"))
summary(broad.pres.model)

#Correlogram to test distance of spatial autocorrelation
correlogram(broad.pres.model, gdat.ml.pres, "figures/M1_broadpres_correlogram.jpg")

#model including spatial correlation variable
rac <- Spat.cor.rep(broad.pres.model,gdat.ml.pres,2000)
broad.pres.model.rac <- glm(presence~landtype*status+area+elev_range+precipitation+temperature+rac, data = gdat.ml.pres, family = binomial(link ="logit"))
summary(broad.pres.model.rac)

#add rac variable to dataframe to include in lsmeans
gdat.ml.pres$rac <- rac

ref.land.stat.rac<-lsmeans(broad.pres.model.rac,pairwise~landtype*status, data = gdat.ml.pres, type="response")
ref.table.land.stat.rac<-as.data.frame(ref.land.stat.rac$lsmeans) 

#check variance inflation factor (should be below 5 for all variables)
vif(broad.pres.model.rac)

#check model assumptions
simulationOutput <- simulateResiduals(fittedModel = broad.pres.model.rac, plot = F)
#access residuals:
residuals(simulationOutput)
#qqplot
plotQQunif(simulationOutput)
#residuals
plotResiduals(simulationOutput)
#overdispersion
testDispersion(simulationOutput)
#zero inflation
testZeroInflation(simulationOutput)

#contrasts
means <- emmeans(broad.pres.model.rac, ~landtype*status)
#look at means order to determine how to write contrasts:
means

#write contrasts:
contrasts <- list("mainland v all islands native" = c(-2,1,1,0,0,0),
                  "mainland v nonoceanic islands native" = c(-1,1,0,0,0,0),
                  "mainland v oceanic islands native" = c(-1,0,1,0,0,0),
                  "nonoceanic v oceanic islands native" = c(0,-1,1,0,0,0),
                  "mainland v all islands naturalized" = c(0,0,0,-2,1,1),
                  "mainland v oceanic islands naturalized" = c(0,0,0,-1,0,1),
                  "mainland v nonoceanic islands naturalized" = c(0,0,0,-1,1,0),
                  "nonoceanic v oceanic islands naturalized" = c(0,0,0,0,-1,1),
                  "mainland v all islands * status" = c(-2,1,1,2,-1,-1),
                  "mainland v oceanic islands * status" = c(-1,0,1,1,0,-1),
                  "mainland v nonoceanic islands * status" =  c(-1,1,0,1,-1,0),
                  "nonoceanic v oceanic island * status" = c(0,-1,1,0,1,-1),
                  "native v naturalized only on mainlands" =  c(-1,0,0,1,0,0),
                  "native v naturalized only on oceanic islands" =  c(0,0,-1,0,0,1),
                  "native v naturalized only on nonoceanic islands" =  c(0,-1,0,0,1,0),
                  "native v naturalized on both island types" =  c(0,-1,-1,0,1,1)
)

#extract results to df:
results <- lsmeans::contrast(means,contrasts)
results.df <- as.data.frame(results)
results.df

#M1 plot Presence naturalized flora
gdat.ml.ntz<- gdat.ml.pres%>%filter(status=="naturalized")                             #subset data to only include naturalized for plot
ref.table.land.stat.rac.ntz <- ref.table.land.stat.rac%>%filter(status=="naturalized") #subset lsmeans to only include naturalized for plot 

#Create a custom color scale
colScale <- scale_colour_manual(values=c("coral"))
fillScale <- scale_fill_manual(values=c("coral"))

#plot the model results for M1 presence of naturalized flora across islands and mainlands
broad.pres.ntz<- 
  ggplot(ref.table.land.stat.rac.ntz, aes(x = landtype, y = lsmean, fill = status, color = status))  + 
  geom_point(position=position_dodge(1), size = 2) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.4,size=1, position=position_dodge(1)) + 
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  geom_bar(stat="identity",position=position_dodge(1), alpha=0.4)+
  geom_point(data = gdat.ml.ntz,aes(x=landtype,y=presence,color=status),position=position_jitterdodge(dodge.width = 1),size=2,alpha=0.3)+
  coord_cartesian(ylim=c(0,1))+
  colScale+
  fillScale+
  labs(y="Probability of N-fixing Plants")+
  xlab(" ")+
  theme_classic(base_size = 30) +
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  theme(legend.position = "top")+
  theme(legend.key = element_rect(colour =c("coral") ))+
  theme(legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))+
  guides(color="none", fill = guide_legend(title="Floral status:",override.aes = list(shape = NA)))+
  theme(axis.text.y=element_text(size=50))+
  theme(axis.text.x=element_text(size=50))+
  ylim(0,1)

broad.pres.ntz

#Saving the plot as a png
png("figures/M1_broad_presence_naturalized.jpg", width=10, height= 12, units='in', res=300)
broad.pres.ntz
dev.off()

#M1 plot Presence native flora
gdat.ml.ntv<- gdat.ml.pres%>%filter(status=="native")                             #subset data to only include native species for plot
ref.table.land.stat.rac.ntv <- ref.table.land.stat.rac%>%filter(status=="native") #subset lsmeans to only include native for plot

#Create a custom color scale
colScale <- scale_colour_manual(values=c("darkcyan"))
fillScale <- scale_fill_manual(values=c("darkcyan"))

#plot the model results for M1 presence of native flora across islands and mainlands
broad.pres.ntv<- 
  ggplot(ref.table.land.stat.rac.ntv, aes(x = landtype, y = lsmean, fill = status, color = status))  + 
  geom_point(position=position_dodge(1), size = 2) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.4,size=1, position=position_dodge(1)) + 
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  geom_bar(stat="identity",position=position_dodge(1), alpha=0.4)+
  geom_point(data = gdat.ml.ntv,aes(x=landtype,y=presence,color=status),position=position_jitterdodge(dodge.width = 1),size=2,alpha=0.3)+
  coord_cartesian(ylim=c(0,1))+
  colScale+
  fillScale+
  labs(y="Probability of N-fixing Plants")+
  xlab(" ")+
  theme_classic(base_size = 30) +
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  theme(legend.position = "top")+
  theme(legend.key = element_rect(colour =c("darkcyan") ))+
  theme(legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))+
  guides(color="none", fill = guide_legend(title="Floral status:",override.aes = list(shape = NA)))+
  theme(axis.text.y=element_text(size=50))+
  theme(axis.text.x=element_text(size=50))+
  theme(axis.title.y=element_text(size=50))

broad.pres.ntv

#Saving the plot as a png
png("figures/M1_broad_presence_native.jpg", width=10, height= 12, units='in', res=300)
broad.pres.ntv
dev.off()


####M2 Broad Proportion (across islands and mainlands including native and naturalized floras)####

#dataset for M2 broad proportion model including mainlands
gdat.ml.prop <- gdat.ref %>%
  dplyr::select(c("entity_ID","nfix","nfixno", "latitude","abs.lat","longitude",                          # Select columns
                  "landtype","status","presence","area","elev_range","precipitation", "temperature"))%>%
  filter(!landtype=="other_island")%>%                                                                    #remove islands with undetermined type
  mutate(area = as.vector(log10((area)+.01)))%>%                                                          #log transform area
  mutate(species = nfix + nfixno) %>%                                                                     #find out rows that have no species counts
  filter(species > 0) %>%                                                                                 #filter our rows that have no species counts
  drop_na() %>%
  mutate_at(c("abs.lat","area","elev_range","precipitation","temperature"), scale) %>%                    #scale all explanatory variables
  filter(precipitation < 3.5)                                                                             #remove outliers

#check correlations
dat <- gdat.ml.prop
names(gdat.ml.prop)
cont.var <- c("abs.lat", "elev_range","area","precipitation","temperature")
Mypairs(dat[,cont.var]) #area and elevation range 0.54, abs.lat and temperature -0.88

#### M2 Model Broad Proportion across islands and mainlands including native and naturalized flora 
broad.prop.model<- glm(cbind(nfix,nfixno)~landtype*status+abs.lat+area+elev_range+precipitation+temperature, data =gdat.ml.prop, family= binomial(link ="logit"))
summary(broad.prop.model)

#model including spatial correlation variable (rac)
rac <- Spat.cor.rep(broad.prop.model,gdat.ml.prop,2000)
broad.prop.model.rac <- glm(cbind(nfix,nfixno)~landtype*status+abs.lat+area+elev_range+precipitation+temperature+rac, data = gdat.ml.prop, family = binomial(link ="logit"))
summary(broad.prop.model.rac)

#M2 with only selected variables after stepwise regression because vif value of abs.lat is too high
broad.prop.model<- glm(cbind(nfix,nfixno)~landtype*status+area+elev_range+precipitation+temperature, data =gdat.ml.prop, family= binomial(link ="logit"))
summary(broad.prop.model)

#Correlogram to test distance of spatial autocorrelation
correlogram(broad.prop.model, gdat.ml.prop, "figures/M2_broadprop_correlogram.jpg")

#M2 including spatial correlation variable (rac)
rac <- Spat.cor.rep(broad.prop.model,gdat.ml.prop,2000)
broad.prop.model.rac <- glm(cbind(nfix,nfixno)~landtype*status+area+elev_range+precipitation+temperature+rac, data = gdat.ml.prop, family = binomial(link ="logit"))
summary(broad.prop.model.rac)

#add rac to df to include in lsmeans.
gdat.ml.prop$rac <- rac

ref.land.stat.rac<-lsmeans(broad.prop.model.rac,pairwise~landtype*status, data = gdat.ml.prop, type="response")
ref.table.land.stat.rac<-as.data.frame(ref.land.stat.rac$lsmeans) 

#check variance inflation factor (should be below 5 for all variables)
vif(broad.prop.model.rac)

#check model assumptions
simulationOutput <- simulateResiduals(fittedModel = broad.prop.model.rac, plot = F)
#access residuals:
residuals(simulationOutput)
#qqplot
plotQQunif(simulationOutput)
#residuals
plotResiduals(simulationOutput)
#overdispersion
testDispersion(simulationOutput)
#zero inflation
testZeroInflation(simulationOutput)

#check residuals for each variable separately
E2 <- resid(broad.prop.model.rac, type = "pearson")
#list of variables landtype, status, area, elev_range, precipitation, temperature
plot(E2 ~ precipitation, data = gdat.ml.prop)
plot(E2 ~ area, data = gdat.ml.prop)
plot(E2 ~ temperature, data = gdat.ml.prop)
plot(E2 ~ elev_range, data = gdat.ml.prop)

#contrasts
means <- emmeans(broad.prop.model.rac, ~landtype*status)
#look at means order to determine how to write contrasts:
means

#write contrasts:
contrasts <- list("mainland v all islands native" = c(-2,1,1,0,0,0),
                  "mainland v nonoceanic islands native" = c(-1,1,0,0,0,0),
                  "mainland v oceanic islands native" = c(-1,0,1,0,0,0),
                  "nonoceanic v oceanic islands native" = c(0,-1,1,0,0,0),
                  "mainland v all islands naturalized" = c(0,0,0,-2,1,1),
                  "mainland v oceanic islands naturalized" = c(0,0,0,-1,0,1),
                  "mainland v nonoceanic islands naturalized" = c(0,0,0,-1,1,0),
                  "nonoceanic v oceanic islands naturalized" = c(0,0,0,0,-1,1),
                  "mainland v all islands * status" = c(-2,1,1,2,-1,-1),
                  "mainland v oceanic islands * status" = c(-1,0,1,1,0,-1),
                  "mainland v nonoceanic islands * status" =  c(-1,1,0,1,-1,0),
                  "nonoceanic v oceanic island * status" = c(0,-1,1,0,1,-1),
                  "native v naturalized only on mainlands" =  c(-1,0,0,1,0,0),
                  "native v naturalized only on oceanic islands" =  c(0,0,-1,0,0,1),
                  "native v naturalized only on nonoceanic islands" =  c(0,-1,0,0,1,0),
                  "native v naturalized on both island types" =  c(0,-1,-1,0,1,1))

#extract results to df:
results <- lsmeans::contrast(means,contrasts)
results.df <- as.data.frame(results)
results.df


#M2 plot Proportion naturalized flora
gdat.prop.ntz<- gdat.ml.prop %>% mutate(propnfix=nfix/(nfix+nfixno))%>%filter(status=="naturalized")  #subset data to only include naturalized for plot
ref.table.land.stat.rac.ntz <- ref.table.land.stat.rac%>%filter(status=="naturalized")                #subset lsmeans to only include naturalized for plot 

#Create a custom color scale
colScale <- scale_colour_manual(values=c("coral","coral","coral","coral"))
fillScale <- scale_fill_manual(values=c("coral","coral","coral","coral"))

broadpropntz<- 
  ggplot(ref.table.land.stat.rac.ntz, aes(x = landtype, y = lsmean, fill = status, color = status))  + 
  geom_point(position=position_dodge(1), size = 2) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.4,size=1, position=position_dodge(1)) + 
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  geom_bar(stat="identity",position=position_dodge(1), alpha=0.4)+
  geom_point(data = gdat.prop.ntz,aes(x=landtype,y=propnfix,color=landtype),position=position_jitterdodge(dodge.width = 1),size=2,alpha=0.4)+
  coord_cartesian(ylim=c(0,0.02))+
  colScale+
  fillScale+
  labs(y="Proportion N-fixing Plant Species")+
  xlab(" ")+
  theme_classic(base_size = 30) +
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  theme(legend.position = "top")+
  theme(legend.key = element_rect(colour ="coral" ))+
  theme(legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))+
  guides(color="none", fill = guide_legend(title="Floral status:",override.aes = list(shape = NA,size = 10)))+
  theme(plot.title=element_text(size=20))

broadpropntz

#Saving the plot as a png
png("figures/M2_broad_proportion_naturalized.jpg", width=10, height= 10, units='in', res=300)
broadpropntz
dev.off()


#Plot M2 Proportion native flora 
gdat.prop.ntv<- gdat.ml.prop %>% mutate(propnfix=nfix/(nfix+nfixno))%>%filter(status=="native") #subset data to only include native for plot
ref.table.land.stat.rac.ntv <- ref.table.land.stat.rac%>%filter(status=="native")               #subset lsmeans to only include native for plot 

#Create a custom color scale
colScale <- scale_colour_manual(values=c("darkcyan","darkcyan"))
fillScale <- scale_fill_manual(values=c("darkcyan","darkcyan"))

broadpropntv<- 
  ggplot(ref.table.land.stat.rac.ntv, aes(x = landtype, y = lsmean, fill = status, color = status))  + 
  geom_point(position=position_dodge(1), size = 2) + 
  geom_errorbar(aes(ymin=lsmean-SE, ymax=lsmean+SE), width=0.4,size=1, position=position_dodge(1)) + 
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  geom_bar(stat="identity",position=position_dodge(1), alpha=0.4)+
  geom_point(data = gdat.prop.ntv,aes(x=landtype,y=propnfix,color=status),position=position_jitterdodge(dodge.width = 1),size=2,alpha=0.4)+
  coord_cartesian(ylim=c(0,0.16))+
  colScale+
  fillScale+
  labs(y="Proportion N-fixing Plant Species")+
  xlab(" ")+
  theme_classic(base_size = 30) +
  theme(legend.justification=c(1,1), legend.position=c(1,1))+
  theme(legend.position = "top")+
  theme(legend.key = element_rect(colour =c("darkcyan","darkcyan") ))+
  theme(legend.background = element_blank(),legend.box.background = element_rect(colour = "black"))+
  guides(color="none", fill = guide_legend(title="Floral status:",override.aes = list(shape = NA)))+
  theme(axis.text.y=element_text(size=30))+
  theme(axis.text.x=element_text(size=30))

broadpropntv

#Saving the plot as a png
png("figures/M2_broad_proportion_native.jpg", width=10, height= 10, units='in', res=300)
broadpropntv
dev.off()


####Data only including NATURALIZED species on OCEANIC ISLANDS (M3,M4)####

#Load dataset and subset to only naturalized species and oceanic islands
gdat.isl.ntz <- gdat.ref %>%
  dplyr::select(c("entity_ID","entity_class","nfix","nfixno", "latitude","abs.lat","longitude",                              # Select columns
                  "landtype","status","presence","area","dist","elev_range","precipitation", "temperature","urbanland"))%>%
  filter(status=="naturalized")%>%                                                                                           #filter data to only include naturalized species
  filter(landtype=="oceanic")%>%                                                                                             #filter data to only include oceanic islands
  mutate(area = as.vector(log10((area)+.01)))%>%                                                                             #log10 transformation of area
  drop_na()%>%                                                                                                               #take out NAs
  mutate_at(c("abs.lat","area","dist","elev_range","precipitation", "temperature","urbanland"), scale)                       #scale all explanatory variables

#check correlations
dat <- gdat.isl.ntz
names(gdat.isl.ntz)
cont.var <- c("abs.lat", "area","elev_range","precipitation","temperature","dist", "urbanland")
Mypairs(dat[,cont.var]) # area and elevation range 0.72, abs.lat and temperature 0.88


####M3 Presence of naturalized N-fixing species on oceanic islands####
model.oc.pres.full<- glm(presence~abs.lat +elev_range	+precipitation+temperature	+area*dist +urbanland*dist, data=gdat.isl.ntz, family=binomial(link ="logit"))
summary(model.oc.pres.full)

#model including spatial correlation variable (rac)
rac <- Spat.cor.rep(model.oc.pres.full,gdat.isl.ntz,2000)
model.oc.pres.rac.full <- glm(presence~abs.lat +elev_range+precipitation+ temperature+area*dist +urbanland*dist+rac, data = gdat.isl.ntz, family = binomial(link ="logit"))
summary(model.oc.pres.rac.full)

#model with selected variables after stepwise regression
model.oc.pres<- glm(presence~abs.lat+area+dist+area:dist, data=gdat.isl.ntz, family=binomial(link ="logit"))
summary(model.oc.pres)

#Correlogram to test distance of spatial autocorrelation
correlogram(model.oc.pres, gdat.isl.ntz, "figures/M3_ntzocpres_correlogram.jpg")

#model including spatial correlation variable (rac)
rac <- Spat.cor.rep(model.oc.pres,gdat.isl.ntz,2000)
model.oc.pres.rac<- glm(presence~ abs.lat+area+dist+area:dist+rac, data = gdat.isl.ntz, family = binomial(link ="logit"))
summary(model.oc.pres.rac)

#tidy and convert the model output into a df for estimates plot
coef_dataM3 <- tidy(model.oc.pres.rac, conf.int = TRUE)

#check variance inflation factor 
vif(model.oc.pres.rac)

#check model assumptions
simulationOutput <- simulateResiduals(fittedModel = model.oc.pres.rac, plot = F)
#access residuals:
residuals(simulationOutput)
#qqplot
plotQQunif(simulationOutput)
#residuals
plotResiduals(simulationOutput)
#overdispersion
testDispersion(simulationOutput)
#zero inflation
testZeroInflation(simulationOutput)


#check residuals for each variable separately
E2 <- resid(model.oc.pres.rac, type = "pearson")
#list of variables: abs.lat, area, dist
plot(E2 ~ abs.lat, data = gdat.isl.ntz)

var.plot<- ggeffect(model.oc.pres.rac, terms=c("abs.lat"), type="re")
effect.plot<- plot(var.plot, colors="coral")
effect.plot

#test area distance effect
var.plot<- ggeffect(model.oc.prop.rac, terms=c("dist","area"), type="re")
plot(var.plot)

#Figure area:distance interaction for presence of naturalized N-fixing species on oceanic islands
#create a color scale
colScale<-scale_color_manual(values = c("small" = "darkred", "medium" = "brown2", "large" = "coral"))
fillScale<-scale_fill_manual(values = c("small" = "darkred", "medium" = "brown2", "large" = "coral"))

sd(gdat.isl.ntz$area) #get standard deviation to create size groups (small, medium, large)

area.minmax <- data.frame(area = c(mean(gdat.isl.ntz$area)-sd(gdat.isl.ntz$area),mean(gdat.isl.ntz$area),mean(gdat.isl.ntz$area)+sd(gdat.isl.ntz$area)), size =c ("small", "medium","large")) #create size groups 

dist.range <- with(gdat.isl.ntz, expand.grid(dist = seq(min(dist), max(dist), length = 1000))) 

ex.grid <- expand.grid(area = area.minmax$area, dist = dist.range$dist)

pred.dat <- ex.grid %>% mutate(rac = mean(rac), abs.lat= mean(gdat.isl.ntz$abs.lat))
pred <- predict.glm(model.oc.pres.rac,type="response",newdata = pred.dat,se=TRUE)

pdat <- cbind(ex.grid,pred) %>%
  left_join(area.minmax)

#plot the area distance interaction position=position_jitterdodge(dodge.width = 1)
areadist_oc_pres_ntz_plot<- ggplot(pdat, aes(x = dist, y = fit, fill=size, color=size))+
  geom_line() +
  geom_ribbon(aes(ymin=fit-se.fit, ymax=fit+se.fit), alpha=0.4)+ 
  geom_point(data = gdat.isl.ntz,aes(x=dist,y=presence),color= "brown2", fill="brown2",alpha=0.2,show.legend = FALSE)+
  theme_minimal()+
  labs(y="Probability of N-fixing plant species")+
  xlab("Distance [km]")+
  theme(axis.text.y=element_text(size=30))+
  theme(axis.text.x=element_text(size=30))+
  theme(axis.title.y=element_text(size=30))+
  theme(axis.title.x=element_text(size=30,margin = margin(t = 10)))+
  colScale+
  fillScale+
  ylim(0,1)+
  guides(
    color = guide_legend(title = ""),  # Combine color legend
    fill = guide_legend(title = ""))+  # Ensure fill uses the same legend title
  theme(legend.key.size = unit(1, 'cm'),      #change legend key size
        legend.key.height = unit(1, 'cm'),    #change legend key height
        legend.key.width = unit(1, 'cm'),     #change legend key width
        legend.title = element_text(size=34), #change legend title font size
        legend.text = element_text(size=30),  #change legend text font size
        legend.position =  c(0.2, 0.85))       #change legend position

areadist_oc_pres_ntz_plot

#Saving the plot as a png
png("figures/M3_areadist_interaction.jpg", width=10, height= 10, units='in', res=300)
areadist_oc_pres_ntz_plot
dev.off()


####M4 Proportion of naturalized N-fixing species on oceanic islands####

#load data and filter out all locations without species counts
oceanic.prop<-gdat.isl.ntz%>%
  mutate(species = nfix + nfixno)%>%
  filter(species >0)%>%
  filter(precipitation< 4) #take out outlier (only do this for scaled version for modelling, remove when creating figures)

#check correlations
names(oceanic.prop)
cont.var <- c("abs.lat", "elev_range","area","dist","precipitation", "temperature", "urbanland")
Mypairs(oceanic.prop[,cont.var]) # area and elevation range 0.77, abs.lat temperature -0.87


#####M4  naturalized proportion on oceanic islands
model.oc.prop.full<- glm(cbind(nfix,nfixno)~abs.lat +area +dist +elev_range	+precipitation+temperature +urbanland	+area:dist +urbanland:dist, data=oceanic.prop, family=binomial(link ="logit"))
summary(model.oc.prop.full)

#model including spatial correlation variable
rac <- Spat.cor.rep(model.oc.prop.full,oceanic.prop,2000)
model.oc.prop.rac.full <- glm(cbind(nfix,nfixno)~abs.lat +area +dist +elev_range+precipitation+temperature +urbanland	+area:dist +urbanland:dist+rac, data = oceanic.prop, family = binomial(link ="logit"))
summary(model.oc.prop.rac.full)

#M4 with only selected variables after stepwise regression
model.oc.prop<- glm(cbind(nfix,nfixno)~abs.lat +area +dist, data=oceanic.prop, family=binomial(link ="logit"))
summary(model.oc.prop)

#Correlogram to test distance of spatial autocorrelation
correlogram(model.oc.prop, oceanic.prop, "figures/M4_ntzocprop_correlogram.jpg")

#model including spatial correlation variable (rac)
rac <- Spat.cor.rep(model.oc.prop,oceanic.prop,2000)
model.oc.prop.rac <- glm(cbind(nfix,nfixno)~abs.lat +area +dist+rac, data = oceanic.prop, family = binomial(link ="logit"))
summary(model.oc.prop.rac)

# tidy and convert the model output into a df for estimates plot
coef_dataM4 <- tidy(model.oc.prop.rac, conf.int = TRUE)

#check variance inflation factor (should be below 5 for all variables)
vif(model.oc.prop.rac)

#check model assumptions
simulationOutput <- simulateResiduals(fittedModel = model.oc.prop.rac, plot = F)
#access residuals:
residuals(simulationOutput)
#qqplot
plotQQunif(simulationOutput)
#residuals
plotResiduals(simulationOutput)
#overdispersion
testDispersion(simulationOutput)
#zero inflation
testZeroInflation(simulationOutput)

#access residuals:
residuals(model.oc.prop.rac)

#check residuals for each variable separately
E2 <- resid(model.oc.prop.rac, type = "pearson")
#list of variables: abs.lat, area, dist
plot(E2 ~ dist, data = oceanic.prop)


####Data only including NATIVE species on OCEANIC ISLANDS (M5,M6)####

#Load dataset and subset to only native species and oceanic islands
gdat.isl.ntv <- gdat.ref %>%
  dplyr::select(c("entity_ID","entity_class","nfix","nfixno", "latitude","abs.lat","longitude",                 # Select columns
                  "landtype","status","presence","area","dist","elev_range","precipitation", "temperature"))%>%
  filter(status=="native")%>%                                                                                   #filter data to only include naturalized species
  filter(landtype=="oceanic")%>%                                                                                #filter data to only include oceanic islands
  mutate(area = as.vector(log10((area)+.01)))%>%                                                                #log10 transformation of area for models only
  drop_na()%>%                                                                                                  #take out NAs
  mutate_at(c("abs.lat","area","dist","elev_range","precipitation", "temperature"), scale)                      #scale all explanatory variables


####M5 presence of native N-fixing species on oceanic islands####

#create data subset for presence analysis
oceanic.pres.ntv <- gdat.isl.ntv%>%
  filter(!entity_ID == "594")# one data point removed due to extreme residual outlier status (see below)

#check correlations
dat <- oceanic.pres.ntv
names(oceanic.pres.ntv)
cont.var <- c("abs.lat", "area","elev_range","precipitation","temperature","dist")
Mypairs(dat[,cont.var]) # area and elevation range 0.73, abs.lat and temperature -0.89

#M5 native presence on oceanic islands
model.oc.pres.full <- glm(presence~abs.lat + area + dist +elev_range +precipitation+temperature +area:dist , data=oceanic.pres.ntv, family=binomial(link ="logit"))
summary(model.oc.pres.full)

#model including spatial correlation variable (rac)
rac <- Spat.cor.rep(model.oc.pres.full,oceanic.pres.ntv,2000)
model.oc.pres.rac.full <- glm(presence~abs.lat + area + dist +elev_range +precipitation+temperature +area:dist +rac, data = oceanic.pres.ntv, family = binomial(link ="logit")) 
summary(model.oc.pres.rac.full)

#M5 with only selected variables after stepwise regression
model.oc.pres <- glm(presence~ area + dist +precipitation+temperature , data=oceanic.pres.ntv, family=binomial(link ="logit"))
summary(model.oc.pres)

#Correlogram to test distance of spatial autocorrelation
correlogram(model.oc.pres, oceanic.pres.ntv, "figures/M5_ntvocpres_correlogram.jpg")

#model including spatial correlation variable (rac)
rac <- Spat.cor.rep(model.oc.pres,oceanic.pres.ntv,2000)
model.oc.pres.rac <- glm(presence~ area + dist +precipitation+temperature +rac, data = oceanic.pres.ntv, family = binomial(link ="logit")) 
summary(model.oc.pres.rac)

# tidy and convert the model output into a df for estimates plot
coef_dataM5 <- tidy(model.oc.pres.rac, conf.int = TRUE)

#check variance inflation factor (should be below 5 for all variables)
vif(model.oc.pres.rac)

#check model assumptions
simulationOutput <- simulateResiduals(fittedModel = model.oc.pres.rac, plot = F)
#access residuals:
residuals(simulationOutput)
#qqplot
plotQQunif(simulationOutput)
#residuals
plotResiduals(simulationOutput)
#overdispersion
testDispersion(simulationOutput)
#zero inflation
testZeroInflation(simulationOutput)

#check residuals for each variable separately
E2 <- resid(model.oc.pres.rac, type = "pearson")
#list of variables: abs.lat, area, dist, elev_range, precipitation
plot(E2 ~ precipitation, data = oceanic.pres.ntv)
#outlier in precipitation

# pinpoint the outlier: entity_ID = 594
outlier <- oceanic.pres.ntv %>% filter(abs.lat < -1) %>% filter(dist < 1.1 & dist > 0.9) %>% filter(area < 1.1 & area > 0.9)

var.plot<- ggeffect(model.oc.pres.rac, terms=c("precipitation"), type="re")
precip.ntv.plot<- plot(var.plot, colors="darkcyan")
precip.ntv.plot

#distance plot
new.dat.oc <- with(oceanic.pres.ntv, expand.grid(dist = seq(min(dist), max(dist), length = 1000))) %>%
  mutate(rac = mean(rac), abs.lat= mean(oceanic.pres.ntv$abs.lat), area=mean(oceanic.pres.ntv$area), temperature= mean(oceanic.pres.ntv$temperature), precipitation=mean(oceanic.pres.ntv$precipitation)) 

pred.ml <- predict(model.oc.pres.rac, newdata = new.dat.oc, type = "response", se = TRUE) %>%
  as.data.frame() %>%
  mutate(dist = new.dat.oc$dist, landtype = "oceanic")


oc_pres_native_dist_plot<-ggplot(pred.ml, aes(x = dist, y = fit))+
  geom_line(color="darkcyan") +
  geom_ribbon(aes(ymin=fit-se.fit, ymax=fit+se.fit), alpha=0.3,color="darkcyan", fill="darkcyan")+ 
  geom_point(data = gdat.isl.ntv,aes(x=dist,y=presence), color="darkcyan", fill="darkcyan",alpha=0.4)+
  theme_minimal()+
  labs(y="Probability of N-fixing plants")+
  xlab("Distance [km]")+  
  theme(axis.text.y=element_text(size=30))+
  theme(axis.text.x=element_text(size=30))+
  theme(axis.title.y=element_text(size=30))+
  theme(axis.title.x=element_text(size=30,margin = margin(t = 10)))

oc_pres_native_dist_plot

#Saving the plot as a png
png("figures/M5_dist_oc_pres_ntv.jpg", width=10, height= 10, units='in', res=300)
oc_pres_native_dist_plot
dev.off()

####M6 proportion of native N-fixing species on oceanic islands ####

#subset data for proportion analysis
oceanic.prop.ntv <- gdat.isl.ntv%>%
  #remove outliers
  filter(!entity_ID == 675) %>%
  filter(!entity_ID == 921) %>%
  filter(!entity_ID == 11474)%>%
  mutate(species = nfix + nfixno)%>% #find out values that have no species counts
  filter(species > 0)%>%             #filter out values that have no species counts
  filter(precipitation<4)

#check correlations
dat <- oceanic.prop.ntv
names(oceanic.prop.ntv)
cont.var <- c("abs.lat", "elev_range","area","dist","precipitation","temperature")
Mypairs(dat[,cont.var]) # area and elevation range 0.73, abs.lat and temperature -0.89


#####M6 Model proportion of native N-fixing species on oceanic islands
model.oc.prop.full<- glm(cbind(nfix,nfixno)~abs.lat+area +dist +elev_range+precipitation +temperature+ area:dist, data=oceanic.prop.ntv, family=binomial(link ="logit"))
summary(model.oc.prop.full)

#model including spatial correlation variable (rac)
rac <- Spat.cor.rep(model.oc.prop.full,oceanic.prop.ntv,2000)
model.oc.prop.rac.full <- glm(cbind(nfix,nfixno)~abs.lat +area +dist+elev_range +precipitation+temperature + area:dist+rac, data = oceanic.prop.ntv, family = binomial(link ="logit"))
summary(model.oc.prop.rac.full)

#selected variables
model.oc.prop<- glm(cbind(nfix,nfixno)~abs.lat+area +dist+precipitation +temperature+ area:dist, data=oceanic.prop.ntv, family=binomial(link ="logit"))
summary(model.oc.prop)

#Correlogram to test distance of spatial autocorrelation
correlogram(model.oc.prop, oceanic.prop.ntv, "figures/M6_ntvocprop_correlogram.jpg")

#model including spatial correlation variable (rac)
rac <- Spat.cor.rep(model.oc.prop,oceanic.prop.ntv,2000)
model.oc.prop.rac <- glm(cbind(nfix,nfixno)~abs.lat+area +dist+precipitation+temperature + area:dist+rac, data = oceanic.prop.ntv, family = binomial(link ="logit"))
summary(model.oc.prop.rac)


#tidy and convert the model output into a df for estimates plot
coef_dataM6 <- tidy(model.oc.prop.rac, conf.int = TRUE)

#check variance inflation factor (should be below 5 for all variables)
vif(model.oc.prop.rac)

#check model assumptions
simulationOutput <- simulateResiduals(fittedModel = model.oc.prop.rac, plot = F)
#access residuals:
residuals(simulationOutput)
#qqplot
plotQQunif(simulationOutput)
#residuals
plotResiduals(simulationOutput)
#overdispersion
testDispersion(simulationOutput)
#zero inflation
testZeroInflation(simulationOutput)


#check residuals for each variable separately
E2 <- resid(model.oc.prop.rac, type = "pearson")
#list of variables: abs.lat, area, dist, elev_range, precipitation
plot(E2 ~ abs.lat, data = oceanic.prop.ntv)
#outlier in precipitation

var.plot<- ggeffect(model.oc.prop.rac, terms=c("precipitation"), type="re")
plot(var.plot, colors="darkcyan")


#####FIGURES: Plot of Estimates #####

###Presence Models M3 and M5

#add distinction between models
coef_dataM3$model <- "Model 3"
coef_dataM5$model <- "Model 5"


#combine estimates of both models into one df
combined_data <- rbind(coef_dataM3, coef_dataM5)
combined_data$term <- factor(combined_data$term, levels = c("rac","temperature", "precipitation","abs.lat","area:dist","dist","area", "(Intercept)")) #Specify the order of variables
combined_data<- combined_data%>% filter(!term=="(Intercept)") 

#Create the combined plot
oc_pres_est_plot <- ggplot(combined_data, aes(x = term, y = estimate, color = model)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +  # Dodge points horizontally
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                position = position_dodge(width = 0.5), 
                width = 0.2, size = 1.2) +  # Dodge error bars
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
  scale_color_manual(values = c("Model 3" = "coral", "Model 5" = "darkcyan")) +
  labs(title = "Coefficient Estimates for Variables in Models 3 and 5") +
  theme_minimal() +
  coord_flip() +  # Flip coordinates for horizontal terms
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.position = "top"  # Adjust legend position if needed
  )

oc_pres_est_plot

#Saving the plot as a png
png("figures/M3M5_estimates_oc_pres.jpg", width=10, height= 9, units='in', res=300)
oc_pres_est_plot
dev.off()


###Proportion Models M4 and M6

#add distinction between models
coef_dataM4$model <- "Model 4"
coef_dataM6$model <- "Model 6"

#combine estimates of both models into one df
combined_data <- rbind(coef_dataM4, coef_dataM6)
combined_data$term <- factor(combined_data$term, levels = c("rac","temperature", "precipitation","abs.lat","area:dist","dist","area", "(Intercept)")) #Specify the order of variables
combined_data<- combined_data%>% filter(!term=="(Intercept)")                                                                                         #filter out intercept because variables are scaled
  
#Create the combined plot
oc_prop_est_plot <- ggplot(combined_data, aes(x = term, y = estimate, color = model)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +  # Dodge points horizontally
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                position = position_dodge(width = 0.5), 
                width = 0.2, size = 1.2) +  # Dodge error bars
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
  scale_color_manual(values = c("Model 4" = "coral", "Model 6" = "darkcyan")) +
  labs(title = "Coefficient Estimates for Variables in Models 4 and 6") +
  theme_minimal() +
  coord_flip() +  # Flip coordinates for horizontal terms
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    legend.position = "top"  # Adjust legend position if needed
  )

oc_prop_est_plot

#Saving the plot as a png
png("figures/M4M6_estimates_oc_prop.jpg", width=10, height= 9, units='in', res=300)
oc_prop_est_plot
dev.off()


####Package versions and citations####

#get packages metadata

pkg_info <-get_pkgs_info(pkgs = c("dplyr","segmented","nlme","mgcv","gridExtra","betareg","MASS","Matrix","lme4","lmerTest","lsmeans","spdep","ggeffects","ggplot2","effects",
                       "ncf","ape","sjPlot","MuMIn","tidyverse","car","V8","DHARMa","broom"), out.dir = getwd())

# Extract package names from the data
used_pkgs <- pkg_info$pkg

# Generate citations only for used packages
cite_packages(pkgs = used_pkgs, output = "file", out.dir = getwd())

