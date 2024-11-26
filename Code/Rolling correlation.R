library(dplyr)
library(tidyr)
library(TTR)
library(EpiEstim)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(data.table)
library(devtools)

### source pika
source_url("https://raw.githubusercontent.com/mrc-ide/pika/master/R/analysis_functions.R")
source_url("https://raw.githubusercontent.com/mrc-ide/pika/master/R/plot_functions.R")


### rolling correlation between Rt and mobility index
typelist <- c('smoothed_within-city movement',	'smoothed_inter-city inflow',	'smoothed_inter-city outflow')
rolling_days <- 14   ## rolling days


data <-read_xlsx("Result/Dataset_wave1.xlsx")
data$date <- as.Date(data$date,format="%Y-%m-%d")

# check city outbreak duration and exclude less than xx days
day_n <- table(data$city)
sort(day_n)
day_table <- data.table(sort(day_n))
temp_set <- subset(day_table,N < 7+21+rolling_days) # last less than 7+7 = 14 days; 7+14 = 21 days; 7+21 = 28 days;
list_exclude <- c(temp_set$V1) 

for (i in 1:length(list_exclude)){
  exclude_city = list_exclude[i]
  data = filter(data, data$city !=exclude_city)
  
}


## do the rolling correlation
for (i in 1:length(typelist)){ 
  mobility_type <- typelist[i]
  data_corr_real <- rolling_corr(dat = data,
                                 date_var = "date",
                                 grp_var = "city",
                                 x_var = "Mean.R.",
                                 y_var = typelist[i],
                                 n = rolling_days)
  
  
  # output <- rbind(data_corr_real)
  # write.table(output,paste0("Result/rolling correlation/wave1",typelist[i],"_",rolling_days,".csv"),sep=",",row.names=FALSE)
}






