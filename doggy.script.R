library(dplyr)
library(fitdistrplus)
library(sads)
library(poweRlaw)
library(ggplot2)
library(scales)
library(broom)
library(tidyr)

rm(list = ls())

setwd("/Users/yufeihui/Desktop/mini.project")
dog <- read.csv("dogbreedreg.csv")

#make the breed name format consistent
#change all breed names to lower case
dog$breed <- tolower(dog$breed)
#get rid of spaces
dog$breed <- gsub(" ", "", dog$breed)

#To do a Wright-Fisher model without mutation
#First remove the 0s in 1997
breeds_1997 <- dog %>% 
  filter(year == '1997' & n != 0)

#make a vector of breed list in 1997
breeds_1997 <- breeds_1997$breed

#remove new breeds that don't exist in 1997 in the following years
dog_no_mutation <- dog %>%
  filter(breed %in% breeds_1997)

str(dog_no_mutation)


#calculate total dog population and the frequency of each breed in each year
dog_no_mutation <- dog_no_mutation %>% 
  group_by(year) %>% 
  mutate(year_total = sum(n),
         frequency = n/sum(n))

#dynamics of breed frequency over the years
ggplot(dog_no_mutation, 
       aes(x = year, y = frequency, 
           fill = as.factor(breed)))+
       geom_area(linewidth = 0.1, colour = "black")+
       guides(colour = "none", fill = "none")+
       scale_fill_viridis_d() +
       theme_classic()+
       theme(aspect.ratio = 1)

#-----------------------------------------------
#log-series model of observed data in 2022
dog2022 <- dog_no_mutation %>% 
  filter(year == '2022')
dog2022counts <- dog2022$n
#remove 0
dog2022counts <- dog2022counts[dog2022counts >0]

#fit logseries distribution using SADs
#Trunc = 0, we have 0 
sads.logser <- fitsad(dog2022counts,'ls', trunc = 0)#logseries
sads.logser


# generate random individuals from the theoretical distributions, given the parameter estimates. 
#log series
logser.ran <- rls(n = 100000, #if more than this, take forever to run
                  N = coef(sads.logser)[1], 
                  alpha = coef(sads.logser)[2])
#write
write.csv(logser.ran, "dog2022_logserfit.csv", row.names=FALSE)

#read
c <- read.csv("dog2022_logserfit.csv", header=TRUE)
logser.ran <- c$x

# examine via a log2(hist) Preston octave plot.
hist(log2(logser.ran))

#subset to counts>=1
logser.ran <- logser.ran[logser.ran >= 1]
#make dataframes
dog2022counts.y<-log10(1-ecdf(dog2022counts)(sort(dog2022counts)))
dog2022counts.x<-log10(sort(dog2022counts))
dog2022counts.dist<-as.data.frame(cbind(dog2022counts.x, dog2022counts.y))

logser.y<-log10(1-ecdf(logser.ran)(sort(logser.ran)))
logser.x<-log10(sort(logser.ran))
logser.dist<-as.data.frame(cbind(logser.x, logser.y))

ggplot()+
  geom_point(data=dog2022counts.dist, 
             aes(x=dog2022counts.x, y=dog2022counts.y),
             shape = 16, size=2, alpha = 0.5,
             color = "grey60")+
  geom_line(data=logser.dist, 
            aes(x=logser.x, y=logser.y))+
  #geom_errorbar(data=rev.simple.summ, aes(x=rev.length, ymin=mean-CI,ymax=mean+CI), size=1, alpha=1, width=0.4)+
  #scale_colour_gradient(low = "grey",high = "red")+  #values=(c(1,3,5,7,9)-1)/(9-1)
  # here the theme
  ylab("Probability of a breed having an abundance > x")+
  xlab("x (breed abundance)")+
  #coord_cartesian(ylim = c(0, 1))+
  scale_x_continuous(breaks=pretty_breaks(n=10))+
  scale_y_continuous(breaks=pretty_breaks(n=10))+
  theme_classic(base_size = 12, base_family = "")+
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))
#------------------------------------

#distribution of the slopes in observed data
#extract all slopes
slopes_ob <- dog_no_mutation %>%
  group_by(breed) %>%
  do(tidy(lm(frequency ~ year, data = .))) %>%
  filter(term == "year") %>%
  dplyr::select(breed, slope = estimate)  #did not work without::

#absolute value of slopes
slopes_ob <-slopes_ob %>% 
  mutate(abs = abs(slope),
         case = "observation")

#summary statistics
summary(slopes_ob$abs)
sd(slopes_ob$abs)

med_ob <- median(slopes_ob$abs)

#histogram of slopes absolute value
ggplot(slopes_ob,
       aes(x = abs))+
  geom_histogram(binwidth = 0.00001)

#too concentrated, change the scale of y to make it look better.
y_trans <- function() {
  trans_new(name = 'ytrans',
            transform = function(x) ifelse(x <= 800, x, 800 + (x - 800)/100),
            inverse = function(x) ifelse(x <= 800, x, 800 + (x - 800) * 100))
}
#compress x axis
x_trans <- function() {
  trans_new(name = 'xtrans',
            transform = function(x) ifelse(x <= 0.0005, x, 0.0005 + (x - 0.0005)/10),
            inverse = function(x) ifelse(x <= 0.0005, x, 0.0005 + (x - 0.0005)*10))
}

#use median because mean is easily influenced by extreme value
ggplot(slopes_ob,
       aes(x = abs))+
  geom_density()+
  scale_y_continuous(trans = y_trans())+
  geom_vline(aes(
    xintercept = med_ob
  ))

  

#one simulation------------------------------------

#get 1997
dog1997 <- filter(dog_no_mutation, year == 1997)

#repeat names according to the count
dog1997 <- rep(dog1997$breed,
                times = dog1997$n
)

#year and population count
yearcnt <- dog_no_mutation %>% distinct(year, .keep_all = TRUE)

#create a list to store results
all_years <- list()

#simulate 1997-2022 population初始化第一年
thisyear <- dog1997
#store the first year in list
all_years[[1]] <- dog1997

#markov process for following years
for(row in 2:nrow(yearcnt)){
  
  size <- yearcnt[row, 'year_total']   #sample size
  
  nextyear <- sample(thisyear,
                     size = size[[1]], 
                     replace = TRUE)
  
  all_years[[row]] <- nextyear
  
  thisyear <- nextyear
}

#fill in NA in every sample—
#to match size of max size
maxrow <- max(yearcnt$year_total)

all_years <- lapply(all_years, 
                    function(x) {
                      if (length(x) < maxrow) {
                        c(x, rep(NA, 
                                 maxrow - length(x)))
                      } else {
                        x
                      }
                    })

#combine by column
#(because if I combine by row, it seems that
#I have to pivot and then pivot longer?)
allyears <- do.call(cbind, all_years)
allyears <- as.data.frame(allyears)

#rename column as years
colnames(allyears) <- 1997:2022

#pivot to long dataframe
allyears <- pivot_longer(allyears,
                         names_to = "year",
                         values_to = "breed",
                         cols = "1997":"2022")

#get rid of NA
allyears <- filter (allyears, 
                    !is.na(breed) )

#collapse into year and count
simulation <- allyears %>% 
  count(year, breed)

simulation <- simulation %>% 
  group_by(year) %>% 
  mutate(total = sum(n),
         prob = n/total)


head(simulation)
#convert year to numeric
simulation$year <- as.numeric(simulation$year)

test <-simulation %>% count(year)
#in 2022 breeds declines to 163 due to random drift

#--------------------------------------
#distribution of the slopes in simulation data
#extract all slopes
slopes_si <- simulation %>%
  group_by(breed) %>%
  do(tidy(lm(prob ~ year, data = .))) %>%
  filter(term == "year") %>%
  dplyr::select(breed, slope = estimate)

#absolute value of slopes
slopes_si <-slopes_si %>% 
  mutate(abs = abs(slope),
         case = "simulation")

#summary statistics 
summary(slopes_si$abs, na.rm=T)
sd(slopes_si$abs, na.rm = T)

med_si <- median(slopes_si$abs, na.rm = T)


#histogram of slopes absolute value in simulation
ggplot(slopes_si,
       aes(x = abs))+
  geom_histogram(binwidth = 0.00001)


#combine with observed slopes
slopes <- rbind(slopes_ob, slopes_si)


#draw on a graph of two cases
#the highest slope is not included in the graph to make it easier to visualize
ggplot(slopes,
       aes(x = abs, color = case))+
  geom_density()+
  scale_y_continuous(trans = y_trans(),
                     expand = c(0, 0),
                     labels = label_scientific())+
  scale_x_continuous(trans = x_trans(),
                     limits = c(0, 0.0035),
                     expand = c(0, 0),
                     labels = label_scientific())+
  xlab("Absolute values of slopes")+
  ylab("Density")+
  scale_color_manual(values=c(observation = "grey20",
                               simulation = "grey60"))+
  geom_vline(aes(xintercept = med_si), 
             linetype="dashed", color = "grey60")+
  geom_vline(aes(xintercept = med_ob),
             linetype="dashed", color = "grey20"    
  )+
  theme_classic(base_size = 12, base_family = "")+
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        legend.position = "none")

#compare simulation and observation slopes
compare1 <- slopes_ob %>% 
  dplyr::select(abs, case)
compare2 <- slopes_si %>% 
  dplyr::select(abs, case)
compare <- rbind(compare1, compare2)

t.test(abs ~ case, data = compare)


 
#now do it 100 times--------------------------------
#THIS TAKES 4 MINUTES
medians <- numeric(100)

for (i in 1:100){
  #create a list to store results
  all_years <- list()
  
  #simulate 1997-2022 population初始化第一年
  thisyear <- dog1997
  #store the first year in list
  all_years[[1]] <- dog1997
  
  #markov process for following years
  for(row in 2:nrow(yearcnt)){
    
    size <- yearcnt[row, 'year_total']   #sample size
    
    nextyear <- sample(thisyear,
                       size = size[[1]], 
                       replace = TRUE)
    
    all_years[[row]] <- nextyear
    
    thisyear <- nextyear
  }
  
  #fill in NA in every sample—
  #to match size of max size
  maxrow <- max(yearcnt$year_total)
  
  all_years <- lapply(all_years, 
                      function(x) {
                        if (length(x) < maxrow) {
                          c(x, rep(NA, 
                                   maxrow - length(x)))
                        } else {
                          x
                        }
                      })
  
  #combine by column
  #(because if I combine by row, it seems that
  #I have to pivot and then pivot longer?)
  allyears <- do.call(cbind, all_years)
  allyears <- as.data.frame(allyears)
  
  #rename column as years
  colnames(allyears) <- 1997:2022
  
  #pivot to long dataframe
  allyears <- pivot_longer(allyears,
                           names_to = "year",
                           values_to = "breed",
                           cols = "1997":"2022")
  
  #get rid of NA
  allyears <- filter (allyears, 
                      !is.na(breed) )
  
  #collapse into year and count
  simulation <- allyears %>% 
    count(year, breed)
  
  simulation <- simulation %>% 
    group_by(year) %>% 
    mutate(total = sum(n),
           prob = n/total)
  
  
  head(simulation)
  #convert year to numeric
  simulation$year <- as.numeric(simulation$year)
  
  #extract all slopes
  slopes_si <- simulation %>%
    group_by(breed) %>%
    do(tidy(lm(prob ~ year, data = .))) %>%
    filter(term == "year") %>%
    dplyr::select(breed, slope = estimate)
  
  #absolute value of slopes
  slopes_si <-slopes_si %>% 
    mutate(abs = abs(slope),
           case = "simulation")
  
  medians[i] <- median(slopes_si$abs, na.rm = T)
  
  
}

medians_df <- data.frame(medians)
#add case
medians_df <- mutate(medians_df, case = "all_si")

#match names to slopes_ob
colnames(medians_df) <- c("abs","case")

#combine and plot
test2 <- slopes_ob %>% 
  ungroup() %>% 
  dplyr::select(c(abs, case))

test2 <- rbind(test2, medians_df)

#density plots of 100medians and the observation
ggplot(test2,
       aes(x = abs, color = case))+
  geom_density()+
  scale_y_continuous(trans = y_trans(),
                     expand = c(0, 0),
                     labels = label_scientific())+
  scale_x_continuous(trans = x_trans(),
                     limits = c(0, 0.0035),
                     expand = c(0, 0),
                     labels = label_scientific())+
  xlab("Absolute values of slopes")+
  ylab("Density")+
  scale_color_manual(values=c(observation = "grey20",
                              all_si = "grey60"),
                     name = NULL,
                     labels = c("Simulation(s)", "Observation"))+
  theme_classic(base_size = 12, base_family = "")+
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        legend.position = "right")


#French bulldog and sausage dog
hotdog_ob <- dog_no_mutation %>% 
  filter(breed == "frenchbulldog" | breed == "dachshundminiaturesmoothhaired") %>% 
  mutate(case="ob")
hotdog_si <- simulation %>% 
  filter(breed == "frenchbulldog" | breed == "dachshundminiaturesmoothhaired") %>% 
  mutate(case="si")

#match names to combine
colnames(hotdog_si) <- c("year","breed","n","year_total","frequency", "case")
hotdog <- rbind(hotdog_ob,hotdog_si)

#plot in presentation but not essay
ggplot(hotdog,aes(x = year, y = frequency, color =interaction(breed, case)))+
  geom_point()+
  geom_smooth(size = 0.6)+
  scale_color_manual(values=c(
    dachshundminiaturesmoothhaired.ob = "dodgerblue4",
    dachshundminiaturesmoothhaired.si = "dodgerblue",
    frenchbulldog.ob = "orange4",
    frenchbulldog.si = "orange"),
                     name = NULL,
                     labels = c("Sausage dogs in observation",
                                "French bulldogs in observation",
                                "Sausage dogs in simulation",
                                "French bulldogs in simulation"))+
  ylab("Popularity")+
  xlab("Year")+
  theme_classic(base_size = 12, base_family = "")+
  theme(axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
        legend.key = element_blank())
