#13.10.2024
#Linda Mueller


###structure; run datacleaning for native species, save resulting dataset, clear workspace, then repeat for naturalized species and combine the two

#load libraries
library(tidyverse)
library(tidyr)

####Native Species#################################################################

#Read datasets
#--------------------
geo.ref <- readRDS("data/GIFT_EXT_geo_Aug_2023.rds") %>%
  mutate(geology= ifelse(entity_class=="Island"&is.na(geology), "other_island", geology))%>%
  mutate(geology= ifelse(entity_class=="Mainland", "other_mainland", geology)) %>%
  filter(!entity_class== "undetermined")

geo.ml<-geo.ref%>% filter(entity_class=="Mainland")

geo <- geo.ref %>%
  mutate(entity_class2 = case_when(geology == "dev" ~ "oceanic",             
                                   geology == "nondev" ~ "nonoceanic",
                                   entity_class == "Mainland" ~ "mainland",
                                   geology== "other_island" ~ "other_island"))
      
ndat <- read.csv("data/Werner_NFix.csv", header = TRUE) %>%                                   
  dplyr::select(c("species", "family", "data_fixing"))                                            

species <- read.csv("data/GIFT_EXT_species_Aug_2023.csv") %>%
  dplyr::select("entity_ID",'species','family') %>%
  filter(entity_ID %in% geo$entity_ID)

geo <- geo %>%
  filter(entity_ID %in% species$entity_ID)

#Merge datasets
#--------------------
dat <- species %>% left_join(ndat, by = c("species", "family")) %>%                       # Merge species and ndat, keep all unassigned species
  left_join(geo, by = "entity_ID") %>%                                                    # Merge species2 with geo
  group_by(entity_ID) %>% distinct(species, .keep_all = TRUE) %>%                         # Remove duplicates of species within each region
  ungroup() # Ungroup

#Species counts
#---------------------------
#Count N-fixing per entity
Nfix.yes <- dat %>%                               # Take dat
  group_by(entity_ID, data_fixing) %>%            # Group by location and N-fixing
  summarise(n = n()) %>%                          # Count sum of data fixing by loc
  filter(data_fixing == "Yes") %>%                # Subset only N-fixing-yes
  dplyr::select(-data_fixing) %>%                 # Remove data_fixing cat
  rename(Nfixnum = n) %>%                         # Rename cols
  ungroup()                                       # Ungroup

#Count non-N-fixing per entity
Nfix.no <- dat %>%                                # Take dat
  group_by(entity_ID, data_fixing) %>%            # Group by location and data fixing
  summarise(n = n()) %>%                          # Count sum of data fixing by loc
  filter(data_fixing == "No") %>%                 # Subset only no
  dplyr:: select(-data_fixing) %>%                # Remove data_fixing cat
  rename(noNfixnum = n) %>%                       # Rename cols
  ungroup()                                       # Ungroup

#Join together
Nfix.merge <- Nfix.yes %>% 
  full_join(Nfix.no, by = "entity_ID") %>%                    # Merge yes and no counts
  mutate(Nfixnum  = ifelse(is.na(Nfixnum), 0, Nfixnum)) %>%   # Replace NA with 0
  mutate(noNfixnum = ifelse(is.na(noNfixnum), 0, noNfixnum))  # Replace NA with 0

#Family counts
#-------------------------------
#Count number of N-fixing species per family and proportion of N-fixing within family
Nfix.fam <- dat %>%
  distinct(species, .keep_all=TRUE) %>%             # Keep only distinct species across dataset
  drop_na(data_fixing) %>%                          # Remove NAs before counting
  group_by(family, data_fixing) %>%                 # Group by family and data_fixing
  summarise(n = n()) %>%                            # Count
  mutate(p = n / sum(n)) %>%                        # Get proportion
  filter(data_fixing == "Yes") %>%                  # Only take yes counts
  dplyr::select(-data_fixing) %>%                   # Remove data fixing col
  rename(Nfixnumfam = n, Nfixpropfam = p) %>%       # Rename cols
  ungroup()                                         # Ungroup

#Count number of non-N-fixing species per family and proportion of non-N-fixing within family
noNfix.fam <- dat %>%
  distinct(species, .keep_all=TRUE) %>%             # Keep only distinct species across dataset
  drop_na(data_fixing) %>%                          # Remove NAs before counting
  group_by(family, data_fixing) %>%                 # Group by family and data_fixing
  summarise(n = n()) %>%                            # Count
  mutate(p = n / sum(n)) %>%                        # Get proportion
  filter(data_fixing == "No") %>%                   # Only take no counts
  dplyr::select(-data_fixing) %>%                   # Remove data fixing col
  rename(noNfixnumfam = n, noNfixpropfam = p) %>%   # Rename cols
  ungroup()                                         # Ungroup

#Merge proportions of N-fixing and non-Nfixing per family
Nfix.fam.prop <- Nfix.fam %>%
  full_join(noNfix.fam, by = "family") %>%                                        # Merge Nfix and Non Nfix proportions
  dplyr::select(-c(Nfixnumfam, noNfixnumfam)) %>%                                 # Remove sums
  mutate(Nfixpropfam = ifelse(is.na(Nfixpropfam), 0, Nfixpropfam)) %>%            # Replace NA with 0
  mutate(noNfixpropfam = ifelse(is.na(noNfixpropfam), 0, noNfixpropfam)) %>%      # Replace NA with 0
  mutate(Nfixpropfam = ifelse(noNfixpropfam == 1, 0, Nfixpropfam)) %>%            # Fill all Nfixprop with zero if noNfixprop is 1
  mutate(noNfixpropfam = ifelse(Nfixpropfam == 1, 0, noNfixpropfam))              # Fill all noNfixprop with zero if Nfixprop is 1

#Save data
#---------------------------
saveRDS(Nfix.fam.prop, "data/Nfix.fam.prop_2024.rds")

#Find number of N-fixing families per entity
FamilypropNfix <- dat %>% 
  left_join(Nfix.fam.prop, by = "family") %>%             # Merge dat and family nprop (apply nprop to each species occurence)
  filter(is.na(data_fixing)) %>%                          # Find where species id resulted in NA data_fixing to apply fam prop
  group_by(entity_ID) %>%                                 # Group by entity_ID
  summarise(sum = sum(Nfixpropfam, na.rm = T)) %>%        # Calculate sum of Nfix prop
  mutate(sum = ceiling(sum)) %>%
  rename(NfixnumF = sum) %>%
  ungroup()

#Find number of non-N-fixing families per entity
FamilypropnoNfix <- dat %>% 
  left_join(Nfix.fam.prop, by = "family") %>%             # Merge dat and family nprop (apply nprop to each species occurence)
  filter(is.na(data_fixing)) %>%                          # Find where species id resulted in NA data_fixing to apply fam prop
  group_by(entity_ID) %>%                                 # Group by entity_ID
  summarise(sum = sum(noNfixpropfam, na.rm = T)) %>%      # Calculate sum of noNfix prop
  mutate(sum = ceiling(sum)) %>%
  rename(noNfixnumF = sum) %>%
  ungroup()

#Merge numbers of N-fixing and non-Nfixing families per entity
Nfix.fam.merge <- FamilypropNfix %>%
  full_join(FamilypropnoNfix, by = "entity_ID")

#Number of N-fix and non-N-fix per entity (species assigned and family proportion assigned)
combined.Nfix.sp.fam <- Nfix.fam.merge %>%
  left_join(Nfix.merge, by = "entity_ID") %>%                          # Combine these counts (sum) with species level dat:
  filter(!is.na(entity_ID)) %>%                                        # Remove na entity_ID
  replace(is.na(.), 0) %>%                                             # Replace na with 0
  mutate(combinedNfix = Nfixnum + NfixnumF) %>%                        # Sum Ncounts(fam) and Nfixnum (sp)
  mutate(combinednoNfix = noNfixnum + noNfixnumF) %>%                  # Sum NoNcounts(fam) and noNfixnum(sp)
  dplyr::select(c(entity_ID, combinedNfix, combinednoNfix)) %>%        # Select columns of interest
  rename(Nfixnum = combinedNfix, noNfixnum = combinednoNfix)           # Rename columns

#No ID species (NA)
#------------------
species.fam.na <- dat %>% 
  left_join(Nfix.fam.prop, by = "family") %>%               # Merge dat and family nprop (apply nprop to each species occurence)
  filter(is.na(noNfixpropfam)) %>%                          # Omit where known to fam prop
  filter(is.na(data_fixing)) %>%                            # Omit where known to species
  mutate(data_fixing = "No")                                # Replace NA with No; this is making all unknowns "No"s: 52004 species

Nfix.no.na <- species.fam.na %>%                            # Take species.fam.na
  group_by(entity_ID, data_fixing) %>%                      # Group by location and data fixing
  summarise(n = n()) %>%                                    # Count sum of data fixing by loc
  filter(data_fixing == "No")%>%                            # Subset only nos
  dplyr::select(-data_fixing) %>%                           # Remove data_fixing cat
  rename(noNfixnum.na = n)                                  # Rename cols

#Add NA to non-Nfixing
combined.Nfix <- combined.Nfix.sp.fam %>% 
  left_join(Nfix.no.na, by = "entity_ID") %>%                               # Join sp.fam with no.na fam
  mutate(noNfixnum.na = ifelse(is.na(noNfixnum.na), 0, noNfixnum.na)) %>%   # make na counts zero if NA
  mutate(noNfixnum_withna = noNfixnum + noNfixnum.na) %>%                   # Sum the sp.fam no with the no.na
  dplyr::select(-c(noNfixnum, noNfixnum.na)) %>%                            # Remove select columns
  replace(is.na(.), 0) %>%                                                  # Replace na with 0
  rename(noNfixnum = noNfixnum_withna)                                      # Rename columns
  
#Join Nfix counts with geo data
#-------------------------
drivers <- combined.Nfix %>% left_join(geo, by = "entity_ID")

#Save data
#---------------------------
saveRDS(drivers, "data/Nfixdrivers_speciesnative_2024.rds")

####CLEAR R#################################################################

rm(list = ls(all.names = TRUE))

####Naturalized species#################################################################

#Read datasets
#-----------------

geo.ref <- readRDS("data/GIFT_EXT_geo_Aug_2023.rds")%>%
  mutate(geology= ifelse(entity_class=="Island" & is.na(geology), "other_island", geology))%>%
  mutate(geology= ifelse(entity_class=="Mainland", "other_mainland", geology))%>%
  filter(!entity_class== "undetermined")

geo <- geo.ref %>%
  mutate(entity_class2 = case_when(geology == "dev" ~ "oceanic",             
                                   geology == "nondev" ~ "nonoceanic",
                                   entity_class == "Mainland" ~ "mainland",
                                   geology== "other_island" ~ "other_island"))

ndat <- read.csv("data/Werner_NFix.csv", header = TRUE) %>%                                   
  dplyr::select(c("species", "family", "data_fixing"))                                            

species <- read.csv("data/glonaf_newcompleteness.csv", header = TRUE) %>%               #Load in species table (GLONAF)
  dplyr::select(c("entity_ID", "species_name", "family_tpl", "completeness_id")) %>%    #Subset to necessary columns 
  rename(family = family_tpl, species = species_name) %>%                               #Rename columns
  filter(completeness_id == 2 | completeness_id == 3)  %>%                              #Select only completeness 2 or 3
  dplyr::select(-c("completeness_id")) %>%   
  filter(entity_ID %in% geo$entity_ID)

geo <- geo %>%
  filter(entity_ID %in% species$entity_ID)

#Merge datasets
#--------------------
dat <- species %>% left_join(ndat, by = c("species", "family")) %>%                       # Merge species and ndat, keep all unassigned species
  left_join(geo, by = "entity_ID") %>%                                                    # Merge species2 with geo
  group_by(entity_ID) %>% distinct(species, .keep_all = TRUE) %>%                         # Remove duplicates of species within each region
  ungroup() # Ungroup
                                                                            

#Species counts
#---------------------------
#Count N-fixing per entity
Nfix.yes <- dat %>%                               # Take dat
  group_by(entity_ID, data_fixing) %>%            # Group by location and N-fixing
  summarise(n = n()) %>%                          # Count sum of data fixing by loc
  filter(data_fixing == "Yes") %>%                # Subset only N-fixing-yes
  dplyr::select(-data_fixing) %>%                 # Remove data_fixing cat
  rename(Nfixnum = n) %>%                         # Rename cols
  ungroup()                                       # Ungroup

#Count non-N-fixing per entity
Nfix.no <- dat %>%                                # Take dat
  group_by(entity_ID, data_fixing) %>%            # Group by location and data fixing
  summarise(n = n()) %>%                          # Count sum of data fixing by loc
  filter(data_fixing == "No") %>%                 # Subset only no
  dplyr:: select(-data_fixing) %>%                # Remove data_fixing cat
  rename(noNfixnum = n) %>%                       # Rename cols
  ungroup()                                       # Ungroup

#Join together
Nfix.merge <- Nfix.yes %>% 
  full_join(Nfix.no, by = "entity_ID") %>%                    # Merge yes and no counts
  mutate(Nfixnum  = ifelse(is.na(Nfixnum), 0, Nfixnum)) %>%   # Replace NA with 0
  mutate(noNfixnum = ifelse(is.na(noNfixnum), 0, noNfixnum))  # Replace NA with 0

#Family counts
#Read native version:
Nfix.fam.prop <- readRDS("data/Nfix.fam.prop_2024.rds")

#Find number of N-fixing families per entity
FamilypropNfix <- dat %>% 
  left_join(Nfix.fam.prop, by = "family") %>%             # Merge dat and family nprop (apply nprop to each species occurence)
  filter(is.na(data_fixing)) %>%                          # Find where species id resulted in NA data_fixing to apply fam prop
  group_by(entity_ID) %>%                                 # Group by entity_ID
  summarise(sum = sum(Nfixpropfam, na.rm = T)) %>%        # Calculate sum of Nfix prop
  mutate(sum = ceiling(sum)) %>%
  rename(NfixnumF = sum) %>%
  ungroup()

#Find number of non-N-fixing families per entity
FamilypropnoNfix <- dat %>% 
  left_join(Nfix.fam.prop, by = "family") %>%             # Merge dat and family nprop (apply nprop to each species occurence)
  filter(is.na(data_fixing)) %>%                          # Find where species id resulted in NA data_fixing to apply fam prop
  group_by(entity_ID) %>%                                 # Group by entity_ID
  summarise(sum = sum(noNfixpropfam, na.rm = T)) %>%      # Calculate sum of noNfix prop
  mutate(sum = ceiling(sum)) %>%
  rename(noNfixnumF = sum) %>%
  ungroup()

#Merge numbers of N-fixing and non-Nfixing families per entity
Nfix.fam.merge <- FamilypropNfix %>%
  full_join(FamilypropnoNfix, by = "entity_ID")

#Number of N-fix and non-N-fix per entity (species assigned and family proportion assigned)
combined.Nfix.sp.fam <- Nfix.fam.merge %>%
  left_join(Nfix.merge, by = "entity_ID") %>%                          # Combine these counts (sum) with species level dat:
  filter(!is.na(entity_ID)) %>%                                        # Remove na entity_ID
  replace(is.na(.), 0) %>%                                             # Replace na with 0
  mutate(combinedNfix = Nfixnum + NfixnumF) %>%                        # Sum Ncounts(fam) and Nfixnum (sp)
  mutate(combinednoNfix = noNfixnum + noNfixnumF) %>%                  # Sum NoNcounts(fam) and noNfixnum(sp)
  dplyr::select(c(entity_ID, combinedNfix, combinednoNfix)) %>%        # Select columns of interest
  rename(Nfixnum = combinedNfix, noNfixnum = combinednoNfix)           # Rename columns

#No ID species (NA)
#------------------
species.fam.na <- dat %>% 
  left_join(Nfix.fam.prop, by = "family") %>%               # Merge dat and family nprop (apply nprop to each species occurence)
  filter(is.na(noNfixpropfam)) %>%                          # Omit where known to fam prop
  filter(is.na(data_fixing)) %>%                            # Omit where known to species
  mutate(data_fixing = "No")                                # Replace NA with No; this is making all unknowns "No"s: 52004 species

Nfix.no.na <- species.fam.na %>%                            # Take species.fam.na
  group_by(entity_ID, data_fixing) %>%                      # Group by location and data fixing
  summarise(n = n()) %>%                                    # Count sum of data fixing by loc
  filter(data_fixing == "No")%>%                            # Subset only nos
  dplyr::select(-data_fixing) %>%                           # Remove data_fixing cat
  rename(noNfixnum.na = n)                                  # Rename cols

#Add NA to non-Nfixing
combined.Nfix <- combined.Nfix.sp.fam %>% 
  left_join(Nfix.no.na, by = "entity_ID") %>%                               # Join sp.fam with no.na fam
  mutate(noNfixnum.na = ifelse(is.na(noNfixnum.na), 0, noNfixnum.na)) %>%   # make na counts zero if NA
  mutate(noNfixnum_withna = noNfixnum + noNfixnum.na) %>%                   # Sum the sp.fam no with the no.na
  dplyr::select(-c(noNfixnum, noNfixnum.na)) %>%                            # Remove select columns
  replace(is.na(.), 0) %>%                                                  # Replace na with 0
  rename(noNfixnum = noNfixnum_withna)                                      # Rename columns

#Join Nfix counts with geo data
#-------------------------
drivers <- combined.Nfix %>% left_join(geo, by = "entity_ID")

#Save data
#---------------------------
saveRDS(drivers, "data/Nfixdrivers_speciesnaturalized_2024.rds")

####CLEAR R#################################################################

rm(list = ls(all.names = TRUE))

####JOIN NATIVE AND NATURALIZED #################################################################

native.drivers <- readRDS("data/Nfixdrivers_speciesnative_2024.rds") %>%
  rename(nfix = Nfixnum, nfixno = noNfixnum)
naturalized.drivers <- readRDS("data/Nfixdrivers_speciesnaturalized_2024.rds") %>%
  dplyr::select("entity_ID","Nfixnum","noNfixnum") %>% rename(nfix.inv = Nfixnum, nfixno.inv = noNfixnum)
  
humans <- read.csv("data/geoentities_ice_LC_new.csv") %>%
  dplyr::select(c("entity_ID", "mean_consensus_full_class_7", "mean_consensus_full_class_9")) %>%
  rename(c(landuse = mean_consensus_full_class_7, urban = mean_consensus_full_class_9)) %>%
  mutate(urbanland = urban + landuse)%>%
  dplyr::select(-"urban")%>% dplyr::select(-"landuse")

#combine native and naturalized in a table like drivers, join with human land use
#change naturalized NAs to zeros
drivers.combined <- full_join(native.drivers, naturalized.drivers, by = "entity_ID") %>% 
  rename(temperature = CHELSA_annual_mean_Temp, precipitation = CHELSA_annual_Prec) %>%
  left_join(humans, by = "entity_ID") %>%
  mutate(nfix.inv = ifelse(is.na(nfix.inv), 0, nfix.inv))%>%
  mutate(nfixno.inv = ifelse(is.na(nfixno.inv), 0, nfixno.inv))%>%
  mutate(nfix = ifelse(is.na(nfix), 0, nfix)) %>%
  mutate(nfixno = ifelse(is.na(nfixno), 0, nfixno))

saveRDS(drivers.combined, "data/combined_drivers_2024.rds")

####CLEAR R#################################################################

rm(list = ls(all.names = TRUE))

#create reference database for modelling
gdat.reference <- readRDS("data/combined_drivers_2024.rds") %>%
  dplyr::select(c("entity_ID","entity_class","entity_class2","nfix.inv","nfixno.inv","nfix","nfixno", "latitude","longitude","area", # Select columns
                  "elev_range","temperature","precipitation","dist","urbanland","Popdensity")) %>%
  mutate(landtype = entity_class2) %>%
  mutate(elev_range = ifelse(elev_range== 0,1,elev_range)) %>%                        # Make 0 elevations 1                  
  mutate(elev_range = ifelse(is.na(elev_range),1,elev_range)) %>%                    # Make unknown elevations 1 
  mutate(abs.lat = abs(latitude))%>%                                                #create absolute latitude
  relocate(abs.lat, .after = latitude) %>%                                          #move column after latitude column
  dplyr::select(-c("entity_class2"))

#assign status native/naturalized
gdat.ref.ntv <- gdat.reference%>% dplyr::select(-c("nfix.inv", "nfixno.inv")) %>% mutate(status="native") %>% rename(nfix= nfix, nfixno=nfixno)     #create new column for status and save native species
gdat.ref.inv <- gdat.reference%>% dplyr::select(-c("nfix", "nfixno")) %>% mutate(status="naturalized") %>% rename(nfix= nfix.inv, nfixno=nfixno.inv) # create new column for status and save naturalized species

gdat.status <- rbind(gdat.ref.ntv, gdat.ref.inv) #combine native and naturalized subsets to get column for status

#create presence column
gdat.presence <- gdat.status %>%
  mutate(presence = case_when(nfix==0 ~ 0, nfix> 0 ~1)) 

saveRDS(gdat.presence, "data/fulldata_for_analysis_2024.rds")
