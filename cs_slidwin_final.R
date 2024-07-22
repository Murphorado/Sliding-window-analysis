# Set working directory and open relevant packages
setwd("C:/Users/murph/OneDrive - University of Edinburgh/Masters/Dissertation/R")
library(tidyverse) # Data manipulation and figures
library(lme4) # Mixed models
library(ggeffects) # Predictions from models
library(AICcmodavg) # AICc values
library(ggplot2)
library(dplyr) # data manipulation

# Import the blue tit phenology and clutch size data
birds <- read.csv("datasets/birds.csv")
View(birds)

# Drop any rows with zero values for species or clutch size. I did this to only include species which were confirmed and to eliminate unoccupied nestboxes.
birds_clean_cs <- birds %>% 
  drop_na(species, cs)

birds_clean_cs <- birds_clean_cs %>%
  filter(fki != "")

birds_clean_cs <- birds_clean_cs %>% 
  drop_na(fki)

# Filter the dataset to remove rows where cs < 2. Dropped any rows with less than 2 for clutch size to make sure this was actually a breeding attempt
birds_clean_cs <- birds_clean_cs %>%
  filter(cs >= 2)

# Only include rows with blue tit species.
birds_clean_cs <- birds_clean_cs %>%
  filter(species == "bluti")

# Create a "site_year" column within the "birds_clean_fed" dataset. This is so I can later combine with the temperature dataset in the sliding window analysis based on the unique values for each site and year. 
birds_clean_cs$site_year <- paste(birds_clean_cs$site, birds_clean_cs$year, sep = " ")

# Calculate the mean first egg date for each site_year
birds_clean_cs <- birds_clean_cs %>%
  group_by(site_year) %>%
  summarise(cs = mean(cs, na.rm = TRUE))

# Separate the site_year column into site and year columns
birds_clean_cs <- birds_clean_cs %>%
  separate(site_year, into = c("site", "year"), sep = " ", remove = FALSE)

# Remove the site and year columns
birds_clean_fed.2 <- birds_clean_fed %>%
  select(-site, -year)

birds_clean_cs <- birds_clean_cs %>% 
  left_join(birds_clean_fed.2, by = "site_year")

## Adding unique IDs for female birds to bird phenology data from adults dataset.

# Creates a "site_year_box" column within the "birds_clean_fed" dataset. This will allow me to later combine with the "adults" dataset.
birds_clean_cs$site_year_box <- paste(birds_clean_cs$site_year, birds_clean_cs$box, sep = " ")

adults <- read.csv("C:/Users/murph/Dropbox/master_data/blue tits/Adults.csv")
adults$site_year_box <- paste(adults$site, adults$year, adults$box, sep=" ")

# Combine datasets based on site_year_box and sex columns with "F" in them.

# Filter the adults dataset for female birds
adults_female <- adults %>% filter(sex == "F")

# Select only the relevant columns
adults_female <- adults_female %>% select(site_year_box, ring)

# Merge the datasets based on the site_year_box column
birds_clean_cs <- birds_clean_cs %>% 
  left_join(adults_female, by = "site_year_box")

# Assign unique IDs to NA values in the 'ring' column. I used ID_row number. So these individuals can still be included as unique replicates for the random effect of individual.
birds_clean_cs <- birds_clean_cs %>%
  mutate(ring = if_else(is.na(ring), paste0("ID_", row_number()), ring))

cs_start_col <- seq(339,2019,168)
# list of values from 290 to 2808 in increments of 168 - weekly start date intervals (I started from temperature windows where NA values stop)

cs_duration <- seq(168,2019,168) 
# list of values from 168 to 2808 in increments of 168 - duration of 1 week or more increasing by 1 week

cs_windows <- data.frame(start_col=rep(cs_start_col,1,each=length(cs_duration)),
                         duration=rep(cs_duration,length(cs_start_col)))
# this repeats every start date for the number of durations there are and vice versa to pair all options

cs_windows$end_col <- cs_windows$start_col+cs_windows$duration-1 
# working out the end column, -1 is included because the start date is included in the window

cs_windows <- cs_windows[-which(cs_windows$end_col>2019),]
# removing any windows that extend past the available data

# Give the windows an ID so it's clear which window they test
cs_windows$window_ID <- paste0(colnames(temperature)[cs_windows$start_col],
                               "_",cs_windows$duration,"hours") 
# Here we've taken the column name for the start date of the window and combined it with the duration of the window 
# The ID now says which ordinal date the window will start on and how long it is in days

# create and empty plot with x axis for the number of days of temp data and y for each window
plot(NA, xlim=c(0,2808), ylim=c(1,nrow(cs_windows)), xlab="Column number", ylab="Different windows") 
# Use a loop to plot each window
for(i in 1:nrow(cs_windows)){ 
  points(cs_windows[i,c("start_col","end_col")], c(i,i), type="l") 
}

cs_base_mod <- lmer(cs ~ 1 + (1|site) + (1|year) + (1|ring), birds_clean_cs, REML=F)
summary(cs_base_mod)

# Make an empty data frame that the results will go into
cs_slidwin <- data.frame(matrix(NA, ncol=8, nrow=nrow(cs_windows)))

# Name the columns
colnames(cs_slidwin) <- c("window_ID", "start_date", "end_date", "deltaAICc", "temp_coef", "temp_SE", "deviation_coef", "deviation_SE")

for(i in 1:nrow(cs_windows)){
  
  #Extract relevant temperature data
  temp_dat <- data.frame(site_year = temperature$site_year,  
                         window_temp = rowMeans(temperature[,cs_windows$start_col[i]
                                                            :cs_windows$end_col[i]]))
  
  # Join temperature and caterpillar data
  cs_windtemp <- left_join(birds_clean_cs, temp_dat, by="site_year")
  
  site_avg_temps_cs <- cs_windtemp %>%
    group_by(site) %>%
    summarize(site_avg_temps_cs = mean(window_temp, na.rm = TRUE))
  
  cs_windtemp <- left_join(cs_windtemp, site_avg_temps_cs, by = "site")
  
  cs_windtemp$deviations_cs<-cs_windtemp$window_temp - cs_windtemp$site_avg_temps_cs
  
  # Run the model 
  mod <- lmer(cs ~ 1 + site_avg_temps_cs + deviations_cs + (1|site) + (1|year) + (1|ring), cs_windtemp, REML=F)
  
  # Store the relevant information
  cs_slidwin$window_ID[i] <- cs_windows$window_ID[i] 
  cs_slidwin$start_date[i] <- cs_windows$start_col[i]
  cs_slidwin$end_date[i] <- cs_windows$end_col[i] 
  cs_slidwin$deltaAICc[i] <- AICc(mod)-AICc(cs_base_mod)  
  cs_slidwin$temp_coef[i] <- summary(mod)$coefficients["site_avg_temps_cs", "Estimate"]  
  cs_slidwin$temp_SE[i] <- summary(mod)$coefficients["site_avg_temps_cs","Std. Error"]
  cs_slidwin$deviation_coef[i] <- summary(mod)$coefficients["deviations_cs", "Estimate"]
  cs_slidwin$deviation_SE[i] <- summary(mod)$coefficients["deviations_cs", "Std. Error"]
  
  # remove elements that were specific to this run of the sliding window
  rm(temp_dat, cs_windtemp, mod)
  
}
View(cs_slidwin)

# blank plot with axis from min to max values for dates and AICc
plot(NA, xlim=c(min(cs_slidwin$start_date),max(cs_slidwin$end_date)),
     ylim=c(min(cs_slidwin$deltaAICc),max(cs_slidwin$deltaAICc)), 
     xlab="Ordinal Date (1 = 1st Jan)", ylab="deltaAICc") 

# use loop to draw the lines
for(i in 1:nrow(cs_slidwin)){
  points(c(cs_slidwin$start_date[i],cs_slidwin$end_date[i]),
         c(cs_slidwin$deltaAICc[i],cs_slidwin$deltaAICc[i]),type="l") 
} 

# line at 2 AICc above the lowest
abline(h=(min(cs_slidwin$deltaAICc)+2), col="red", lty="dashed")

# Row number for the window with the lowest AIC
cs_wind_row <- which(cs_slidwin$deltaAICc==min(cs_slidwin$deltaAICc))

# ID for the best window
cs_wind_ID <- cs_slidwin$window_ID[cs_wind_row] 

#The row number is the same in the cater_windows and cater_slidwin dataframes
cs_wind_row

which(cs_windows$window_ID==cs_wind_ID)

# Mean temperature during the identified window
cs_best_temp_wind <- data.frame(site_year = temperature$site_year,  
                                best_temp = rowMeans(temperature[,cs_windows$start_col[cs_wind_row]:cs_windows$end_col[cs_wind_row]])) 

# Join temperature and bird phenology data
cs_best_temp <- left_join(birds_clean_cs, cs_best_temp_wind, by="site_year")

cs_best_temp <- cs_best_temp %>% 
  drop_na(best_temp)

best_site_avg_temps_cs <- cs_best_temp %>%
  group_by(site) %>%
  summarize(best_site_avg_temps_cs = mean(best_temp, na.rm = TRUE))

cs_best_temp <- left_join(cs_best_temp, best_site_avg_temps_cs, by = "site")

cs_best_temp$deviations<-cs_best_temp$best_temp - cs_best_temp$best_site_avg_temps_cs

mean_fed <- birds_clean_cs %>%
  group_by(site_year) %>%
  summarize(mean_fed = mean(fed, na.rm = TRUE))

cs_best_temp <- left_join(cs_best_temp, mean_fed, by = "site_year")

cs_best_temp$fed_deviation <- cs_best_temp$fed - cs_best_temp$mean_fed

# Run the same model as before but with REML=TRUE
cs_mod <- lmer(cs~best_site_avg_temps_cs+deviations+fed_deviation+(1|year)+(1|site)+(1|ring), cs_best_temp, REML=FALSE) 
summary(cs_mod)

caterpillars <- read.csv("datasets/caterpillars.csv") # Import dataset
# Extract relevant columns
caterpillars <- data.frame(year=caterpillars$year, site=caterpillars$site, tree=caterpillars$tree, 
                           tree.species=caterpillars$tree.species, date=caterpillars$date, caterpillars=caterpillars$caterpillars, 
                           caterpillar.mass=caterpillars$caterpillar.mass,beats=caterpillars$beats)

View(caterpillars)

# Create site_year column
caterpillars$site_year <- paste(caterpillars$site, caterpillars$year, sep=" ")

# Calculate the total number of caterpillars and the number of unique trees for each site_year. This is so I can calculate a standardised "caterpillars per tree" metric for each "site_year".
caterpillars_std <- caterpillars %>%
  group_by(site_year) %>%
  summarize(
    total_caterpillars = sum(caterpillars, na.rm = TRUE),
    beats_per_siteyear = sum(beats)
  ) %>%
  mutate(caterpillars_per_tree = total_caterpillars / beats_per_siteyear)

# Add the caterpillars_per_tree to the best temperature window data to add it as a fixed effect in the model

cs_best_window.cat <- left_join(cs_best_temp, caterpillars_std, by="site_year")
cs_mod.cat <- lmer(cs~best_site_avg_temps_cs+deviations+fed_deviation+log1p(caterpillars_per_tree)+(1|year)+(1|site)+(1|ring), cs_best_window.cat, REML=TRUE)
summary(cs_mod.cat)

trees <- read.csv("datasets/trees.csv") # Import dataset

# Remove any missing values for tree fbb
View(trees)
trees$fbb <- as.numeric(trees$fbb)
trees <- trees %>% drop_na(fbb)

# Create a "site_year" column in the trees dataset to combine it later with the mean fbb.
trees$site_year <- paste(trees$site, trees$year, sep=" ")
trees <- data.frame(treeID=trees$treeID, 
                    speciesgroup=trees$speciesgroup, fbb=trees$fbb, cbb=trees$cbb, 
                    flf=trees$flf, clf=trees$clf, site_year=trees$site_year, site=trees$site, year=trees$year)

# Calculate mean first bud burst for each "site_year"
mean_fbb <- trees %>%
  group_by(site_year) %>%
  summarize(mean_fbb = mean(fbb, na.rm = TRUE))

# Join with trees dataset
trees <- left_join(trees, mean_fbb, by="site_year")

site_avg_fbb <- trees %>%
  group_by(site) %>%
  summarize(site_avg_fbb = mean(fbb, na.rm = TRUE))

trees <- left_join(trees, site_avg_fbb, by = "site")

trees$fbb_deviations<-trees$mean_fbb - trees$site_avg_fbb

# Assuming your dataset is called 'trees' and the first bud burst column is named 'fbb'
# Filter for birch species and calculate the mean fbb for each site_year
mean_fbb_birch <- trees %>%
  filter(speciesgroup == "Birch") %>%
  group_by(site_year) %>%
  summarize(mean_fbb_birch = mean(fbb, na.rm = TRUE))

trees <- left_join(trees, mean_fbb_birch, by="site_year")

fbb_deviations <- trees %>% group_by(site_year) %>% summarize(fbb_deviations = mean(fbb_deviations, na.rm=TRUE))
cs_best_window.cat.tree <- left_join(cs_best_window.cat, fbb_deviations, by="site_year")

mean_fbb <- trees %>% group_by(site_year) %>% summarize(mean_fbb = mean(mean_fbb, na.rm=TRUE))
cs_best_window.cat.tree <- left_join(cs_best_window.cat.tree, mean_fbb, by="site_year")

site_avg_fbb <- trees %>% group_by(site_year) %>% summarize(site_avg_fbb = mean(site_avg_fbb, na.rm=TRUE))
cs_best_window.cat.tree <- left_join(cs_best_window.cat.tree, site_avg_fbb, by="site_year")

mean_fbb_birch <- trees %>% group_by(site_year) %>% summarize(mean_fbb_birch = mean(mean_fbb_birch, na.rm=TRUE))
cs_best_window.cat.tree <- left_join(cs_best_window.cat.tree, mean_fbb_birch, by="site_year")

# Calculate the total number of trees in each site_year
total_trees_per_site_year <- trees %>%
  group_by(site_year) %>%
  summarize(total_trees = n(), .groups = 'drop')

# Calculate the number of Oak trees in each site_year
oak_trees_per_site_year <- trees %>%
  filter(speciesgroup == "Oak") %>%
  group_by(site_year) %>%
  summarize(oak_trees = n(), .groups = 'drop')

# Perform a full join to include all site_year values
proportion_oak_trees <- total_trees_per_site_year %>%
  full_join(oak_trees_per_site_year, by = "site_year") %>%
  # Replace NA values with 0
  mutate(oak_trees = ifelse(is.na(oak_trees), 0, oak_trees)) %>%
  # Calculate the proportion of Oak trees
  mutate(proportion_oak = oak_trees / total_trees)

# View the resulting data
print(proportion_oak_trees)

proportion_oak_trees <- data.frame(proportion_oak_trees)

oak_proportion <- proportion_oak_trees %>% select(site_year, proportion_oak)

trees <- left_join(trees, oak_proportion, by="site_year")

trees <- trees %>% select(speciesgroup, site_year, mean_fbb, mean_fbb_birch, site_avg_fbb, fbb_deviations, proportion_oak)

cs_best_window.cat.tree <- left_join(cs_best_window.cat, trees, by="site_year")

View(cs_best_window.cat.tree)

# Run model with both caterpillars per tree and average fbb as fixed effects
cs_mod.cat.tree <- lmer(cs~best_site_avg_temps_cs+deviations+fed_deviation+log1p(caterpillars_per_tree)+site_avg_fbb+fbb_deviations+(1|year)+(1|site)+(1|ring), cs_best_window.cat.tree, REML=TRUE)
summary(cs_mod.cat.tree)

testDispersion(cs_mod.cat.tree)
simulationOutput_cs <- simulateResiduals(fittedModel = cs_mod.cat.tree, plot = F)
residuals(simulationOutput_cs)
plot(simulationOutput_cs)
testOutliers(simulationOutput_cs, type="bootstrap")

cs_mod.cat.tree.birch <- lmer(cs~best_site_avg_temps_cs+deviations+fed_deviation+log1p(caterpillars_per_tree)+mean_fbb+mean_fbb_birch+(1|year)+(1|site)+(1|ring), cs_best_window.cat.tree, REML=TRUE)
summary(cs_mod.cat.tree.birch)

plot_model(cs_mod.cat.tree, type='diag')
sum(cs_best_window.cat.tree$caterpillars_per_tree > 9, na.rm=TRUE)
cs_best_window.cat.tree <- cs_best_window.cat.tree %>% filter(caterpillars_per_tree < 3)
summary(cs_mod.cat.tree.diag)
cs_mod.cat.tree.diag <- lmer(log1p(cs)~best_site_avg_temps_cs+deviations+fed_deviation+log1p(caterpillars_per_tree)+site_avg_fbb+fbb_deviations+mean_fbb_birch+(1|year)+(1|site)+(1|ring), cs_best_window.cat.tree.diag, REML=TRUE)
plot_model(cs_mod.cat.tree.diag, type='diag')

cs_best_window.cat.tree.diag <- cs_best_window.cat.tree.diag[cs_best_window.cat.tree.diag$caterpillars_per_tree <= 3, ]


cs_mod.cat.tree.2 <- lmer(cs~best_site_avg_temps_cs+best_temp+fed_deviation+log1p(caterpillars_per_tree)+site_avg_fbb+fbb_deviations+(1|year)+(1|site)+(1|ring), cs_best_window.cat.tree, REML=TRUE)
summary(cs_mod.cat.tree.2)

cs_mod.cat.tree.3 <- lmer(cs~best_site_avg_temps_cs+deviations+fed_deviation+log1p(caterpillars_per_tree)+site_avg_fbb+mean_fbb+(1|year)+(1|site)+(1|ring), cs_best_window.cat.tree, REML=TRUE)
summary(cs_mod.cat.tree.3)

# Store the slope and confidence intervals 
cs_cater_coef <- summary(cs_mod.cat.tree)$coefficients["log1p(caterpillars_per_tree)","Estimate"] 
cs_cater_confint <- confint(cs_mod.cat.tree)["log1p(caterpillars_per_tree)",]

cs_site_avg_fbb_coef <- summary(cs_mod.cat.tree.diag)$coefficients["site_avg_fbb","Estimate"]
cs_site_avg_fbb_confint <- confint(cs_mod.cat.tree.diag)["site_avg_fbb",]

cs_fbb_deviations_coef <- summary(cs_mod.cat.tree.diag)$coefficients["fbb_deviations","Estimate"]
cs_fbb_deviations_confint <- confint(cs_mod.cat.tree.diag)["fbb_deviations"]

# This shows the mean and 95% confidence intervals for the slope in units of days per°C
cs_cater_coef 
cs_cater_confint

cs_mean_fbb_coef
cs_mean_fbb_confint

# Use ggpredict to get estimates of clutch size across the range of temperatures included in the data
pred_cs_cater <- ggpredict(cs_mod.cat.tree, "caterpillars_per_tree", type="fixed")
pred_cs_fbb <- ggpredict(cs_mod.cat.tree.diag, "mean_fbb", type="fixed")

pred_cs_fbb_site <- ggpredict(cs_mod.cat.tree.diag, "site_avg_fbb")
pred_cs_fbb_dev <- ggpredict(cs_mod.cat.tree.diag, "fbb_deviations")

#Plot the mean prediction and CIs with the data
ggplot(pred_cs_cater, aes(x,predicted))+ 
  geom_line(lwd=1.2, colour="#006D5B")+ 
  geom_point(data=cs_best_window.cat.tree, aes(caterpillars_per_tree, cs), colour="#77DD77")+ 
  geom_ribbon(data=pred_cs_cater, aes(x=x, ymin=conf.low, ymax=conf.high), alpha=0.25, colour="#006D5B")+
  xlab("Caterpillars per beat")+ 
  ylab("Clutch size")+ 
  theme_bw()

#Plot the mean prediction and CIs with the data
cs_site_fbb_plot <- ggplot(pred_cs_fbb_site, aes(x,predicted))+ 
  geom_line(lwd=1.2, colour="darkblue")+ 
  geom_point(data=cs_best_window.cat.tree.diag, aes(site_avg_fbb, cs), colour="blue")+ 
  geom_ribbon(data=pred_cs_fbb_site, aes(x=x, ymin=conf.low, ymax=conf.high), alpha=0.25, colour="darkblue")+
  scale_fill_paletteer_d("nationalparkcolors::Acadia")+
  xlab("Site mean first bud burst date")+ 
  ylab("Clutch size")+ 
  theme_bw()

#Plot the mean prediction and CIs with the data
cs_dev_fbb_plot <- ggplot(pred_cs_fbb_dev, aes(x,predicted))+ 
  geom_line(lwd=1.2, colour="blue")+ 
  geom_point(data=cs_best_window.cat.tree.diag, aes(fbb_deviations, cs), colour="lightblue")+ 
  geom_ribbon(data=pred_cs_fbb_dev, aes(x=x, ymin=conf.low, ymax=conf.high), alpha=0.25, colour="blue")+
  scale_fill_paletteer_d("nationalparkcolors::Acadia")+
  xlab("Annual deviations from site mean first bud burst date")+ 
  ylab("Clutch size")+ 
  theme_bw()

cs_fbb_plots <- cs_site_fbb_plot + cs_dev_fbb_plot
print(cs_fbb_plots)

n_rand <- 200 # Number of randomisations

cs_rand_deltaAICc <- c() 
# Empty vector space to store the lowest AICc from each run of the full sliding window

for(j in 1:n_rand){
  
  cs_rand <- c() 
  # Empty vector space to store the AICc from each model in one sliding window
  
  
  temperature$rand_site_year <- sample(temperature$site_year, replace=FALSE) 
  # Randomise the site-years: we want each site-year to appear once so don't replace them after sampling
  
  #The sliding window loop is the same as before except now using the rand_site_year and only storing deltaAICc
  for(i in 1:nrow(cs_windows)){
    rand_temp_dat <- data.frame(site_year = temperature$rand_site_year,
                                window_temp = rowMeans(temperature[,cs_windows$start_col[i]
                                                                   :cs_windows$end_col[i]]))
    
    rand_cs_windtemp <- left_join(birds_clean_cs, rand_temp_dat, by="site_year")
    
    rand_site_avg_temps_cs <- rand_cs_windtemp %>%
      group_by(site) %>%
      summarize(rand_site_avg_temps_cs = mean(window_temp, na.rm = TRUE))
    
    rand_cs_windtemp <- left_join(rand_cs_windtemp, rand_site_avg_temps_cs, by = "site")
    
    rand_cs_windtemp$rand_deviations_cs<-rand_cs_windtemp$window_temp - rand_cs_windtemp$rand_site_avg_temps_cs
    
    mod <- lmer(cs ~ 1 + rand_site_avg_temps_cs + rand_deviations_cs + (1|site) + (1|year) + (1|ring), rand_cs_windtemp, REML=F)
    
    cs_rand[i] <- AICc(mod)-AICc(cs_base_mod)
    
    rm(rand_temp_dat, rand_cs_windtemp, mod)
  }
  
  # Store the minimum deltaAICc from each sliding window with randomised data
  cs_rand_deltaAICc[j] <- min(cs_rand)
  
  # Clear the used randomised data each time
  temperature$rand_site_year <- NULL 
  rm(cs_rand)
  
  #print(paste(j,"of",n_rand)) 
  # If you remove the hashtag at the start of the line above the loop will print how many randomisations have been run 
}

ggplot()+ 
  geom_histogram(aes(cs_rand_deltaAICc), binwidth=1,colour="black",fill="white")+
  xlab("Permutation test minimum deltaAICc values")+
  theme_bw()

ggplot()+ 
  geom_histogram(aes(cs_rand_deltaAICc), binwidth=1,colour="black",fill="white")+
  geom_vline(aes(xintercept=min(cs_slidwin$deltaAICc)), color="red", linetype="dashed")+ 
  xlab("Permutation test minimum deltaAICc values")+
  theme_bw()

probfunc <- ecdf(cs_rand_deltaAICc) 
# This computes a cumulative distribution function for the permutation AICc results

#You can see the cumulative distribution function by plotting it
plot(probfunc, xlab="deltaAICc", ylab="Percentile")

probfunc(min(cs_slidwin$deltaAICc)) # This gives the percentile of the permutation distribution that the true AICc falls in. In this case it also represents the probability of observing the result by chance.