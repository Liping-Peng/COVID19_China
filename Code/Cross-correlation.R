library(tseries)
library(boot)
library(ggplot2)
library(dplyr)
library(boot)
library(tidyr)
library(data.table)
library(readxl)



### function: temporal correlation using bootstrap data
temp_corr.fun <- function(tsb,factor1,factor2,corr_out) {
  temp_corr_out <- data.frame()
  data_roll <- data.frame(tsb)
  
  city_corr <- cor.test(data_roll[,1],data_roll[,2], alternative = "two.side", method = "pearson", conf.level = 0.95)
  print (city_corr)
  corr_out <- city_corr$estimate
  
  return (corr_out)

}



### Cross-correlation between Rt and mobility index
list_wave <- c("wave1","wave2")
list_mobility <- c('smoothed_within-city movement',	'smoothed_inter-city inflow',	'smoothed_inter-city outflow')


for (k in 1:length(list_wave)){
  wave_dataset <- list_wave[k]
  
  for (j in 1:length(list_mobility)){
    mobility_type <- list_mobility[j]
    
    ## load dataset
    path_in <- paste0("Result/Dataset_",wave_dataset,".xlsx")
    data <- read_xlsx(path_in)
    data$date <- as.Date(data$date,origin="1970-01-01") #format='%m/%d/%Y'
    
    
    output <- data.frame()
    citylist <- unique(data$city)
    for (i in 1:length(citylist)){ 
      city_name <-  citylist[i]
      
      ### prepare data
      city_data <- subset(data,city==city_name)
      # city_data$date <- as.Date(city_data$date,)
      city_data <- na.omit(city_data[c('Mean(R)',mobility_type,'date','city')])
      N_sample <- nrow(city_data)
      ts_city <- ts(city_data)
      
      ## bootstrap
      res_city <- tsboot(ts_city,temp_corr.fun,R=200,l=length(ts_city)^(1/3),sim="geom")
      
      
      ## manage output data 
      city_output <- data.frame(res_city$t)
      city_output <- c(city_name,N_sample,t(city_output))
      output <- rbind(output,city_output)
      
    }
    
    colnames(output)[1] <- "city"
    colnames(output)[2] <- "N_sample"
    colnames(output)[3:ncol(output)] <- 1:(ncol(output)-2)
    # path_out <- paste0("Result/",wave_dataset,"_",mobility_type,".csv")
    # write.table(output,path_out,sep=",",row.names = FALSE)
    
  }
  
}





####  check optimal lag 
data <- read_xlsx(paste0("Result/Dataset_wave1.xlsx"))
data$date <- as.Date(data$date,format='%m/%d/%Y',origin="1970-01-01")
df <- na.omit(data[c('date','city', 'Mean(R)', 'smoothed_within-city movement',	'smoothed_inter-city inflow',	'smoothed_inter-city outflow')])

par(mfrow = c(3, 1),mar = c(2, 3, 4, 2) + 0.1)
ccf1 <- ccf(df$`Mean(R)`, df$`smoothed_within-city movement`, lag.max = 20, main="Within-city movement")
ccf2 <- ccf(df$`Mean(R)`, df$`smoothed_inter-city inflow`, lag.max = 20, main="Inter-city inflow")
ccf3 <- ccf(df$`Mean(R)`, df$`smoothed_inter-city outflow`, lag.max = 20, main="Inter-city outflow")






### Pooled cross-correlation for each mobility index and each wave
wave_dataset <- "wave1"  #"wave1", 'wave2'
list_mobility <- c('within-city movement',	'inter-city inflow',	'inter-city outflow')


for (j in 1:length(list_mobility)){
  mobility_type <- list_mobility[j]
  
  path_in <- paste0("Result/cross correlation/",wave_dataset,"_",mobility_type,".csv")
  data <- read.csv(path_in)
  city_df <- data[,c("city","N_sample")]
  
  
  ## compute city correlation se
  df_se <- data.frame()
  for (h in 1:nrow(data)){
    city_name_se <- data[h,1]
    test02 <-t(data[h,3:ncol(data)])
    city_se <- sd(test02)
    
    vect_se <- c(city_name_se,city_se)
    df_se <- rbind(df_se,vect_se)
  }
  
  
  
  ## data fisher transformation
  data_z <- data.table(FisherZ(data[,3:ncol(data)]))
  
  ## compute combine average rolling corr, then inverse fisher transformation
  data_z[sapply(data_z, is.infinite)] <- NA  ## for "Heilongjiang_Jiamusi", fisher'z r = inf
  aver_mean_z <- data.frame(rowMeans(data_z,na.rm = TRUE))
  aver_mean <- FisherZInv(aver_mean_z)
  
  
  ## merge data
  city_data <- cbind(city_df,aver_mean_z,aver_mean)
  colnames(city_data) <- c("city","N_sample","mean_z","mean")
  
  
  output_semi <- data.frame()
  ## compute temporal correlation for each city (can use CorCI() or r.con() instead)
  
  for (k in 1:nrow(city_data)){
    ## compute combine average se, then compute CI, then inverse fisher transformation CI
    city_name <- city_data[k,1]
    N_sample <- city_data[k,2]
    city_mean_z <- city_data[k,3]
    city_mean <- city_data[k,4]
    
    city_se_z <- 1/sqrt((N_sample-3))
    city_lb_z <- city_mean_z - 1.96*city_se_z
    city_ub_z <- city_mean_z + 1.96*city_se_z
    city_lb <- FisherZInv(city_lb_z)
    city_ub <- FisherZInv(city_ub_z)
    
    vect <- c(city_name,city_mean_z,city_se_z,city_mean,city_lb,city_ub) 
    output_semi <- rbind(output_semi,vect)
    
  }
  
  output <-merge(output_semi,df_se)
  colnames(output)<- c("city","z_mean_corr","z_se","mean_corr","lb_corr","ub_corr","se")
  
  ## compute pooled correlation and se
  N_city <- nrow(output)
  pool_sum_mean_z <- sum(as.numeric(output$z_mean_corr))
  pool_mean_z <- pool_sum_mean_z/N_city
  pool_sum_se_z <- sqrt(sum(as.numeric(output$z_se)))
  pool_se_z <- pool_sum_se_z/N_city
  
  pool_lb_z <- pool_mean_z - 1.96*pool_se_z
  pool_ub_z <- pool_mean_z + 1.96*pool_se_z
  
  pool_mean <- FisherZInv(pool_mean_z)
  pool_lb <- FisherZInv(pool_lb_z)
  pool_ub <- FisherZInv(pool_ub_z)
  
  ## combine all data for output_final
  output_final <- output[,c("city","mean_corr","lb_corr","ub_corr","se")]
  pool_vect <- c("Pooled",pool_mean,pool_lb,pool_ub,NA)
  output_final <- rbind(output_final,pool_vect)
  
  ## add province column
  output_final$province <- strsplit(output_final$city, '_')[1]
  temp_test <- strsplit(output_final$city, '_')
  tenmp_df <- unlist(temp_test)[1:length(temp_test)*2-1]
  output_final$province <- unlist(temp_test)[1:length(temp_test)*2-1]
  
  # path_out <- paste0("Result/cross correlation/pooled_",wave_dataset,"_",mobility_type,".csv")
  # write.table(output_final,path_out,sep=",",row.names = FALSE)
  
}


