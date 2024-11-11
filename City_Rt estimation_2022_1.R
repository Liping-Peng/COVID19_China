library("rstan")
library(EpiNow2)
library(rstudioapi)
library(data.table)
library(TeachingDemos)
Sys.setlocale("LC_TIME", "English")



setwd("C:/Users/Liz Pang/Dropbox (Personal)/LIPING/China COVID-19/3 waves/Github_data/")


outbreak_date <- read.csv("Data/outbreak date_2022_1.csv",encoding = "UTF-8",header=T,sep=",")
outbreak_date$start_date <- as.Date(outbreak_date$start_date ,format='%m/%d/%Y',origin="1970-01-01")
outbreak_date$end_date <- as.Date(outbreak_date$end_date ,format='%m/%d/%Y',origin="1970-01-01")


### Rt estimation for each city with outbreak in Omicron wave1
citylist <- unique(outbreak_date$city)
for (i in 1:length(citylist)){
  # select city
  city_name <- citylist[i]
  
  # load case data
  data_path = paste0("Data/city cases_for Rt/",city_name,".csv")
  data_input <-read.csv(data_path,encoding = "UTF-8",header=T,sep=",")
  data_input$date <- as.Date(data_input$date,format='%m/%d/%Y',origin="1970-01-01")   ## '%Y-%m-%d'
  
  data_cases <- data_input[c('date','cases')]
  
  
  start_date <- subset(outbreak_date,city==city_name)$start_date
  end_date <- subset(outbreak_date,city==city_name)$end_date

  
  # prepare reported cases, for Rt estimation next step
  temp <- data.table(data_input)
  date_start <- temp[temp$date==start_date,which=TRUE]
  date_end <- temp[temp$date==end_date,which=TRUE]
  print (data.table(data_cases)[date_start])
  print (data.table(data_cases)[date_end])

    
  reported_cases <- data.table(data_cases)[date_start:(date_end)]   # nrow(data_cases)   
  nrow(reported_cases)
  reported_cases[,2] <-as.numeric((unlist(reported_cases[,2])))
  # reported_cases[,1] <-as.Date(unlist(reported_cases[,1]),format='%Y-%m-%d',origin="1970-01-01")
  head(reported_cases)
  
  
  d1 <- reported_cases
  sum(d1$cases)
  
  
  #### step 1
  ## use incidental to estimate infection time series from reported time series
  ## output 1000 estimated infection time series
  library(incidental)
  # incubation dist
  meanvec <- c(3.2,2.2)
  incpara <- meanvec
  incpara[1] <- log(meanvec[1])-0.5*log((meanvec[2]/meanvec[1])^2+1)
  incpara[2] <- sqrt(log((meanvec[2]/meanvec[1])^2+1))
  inc_vec <- (plnorm(1:15,incpara[1],incpara[2])-plnorm(0:14,incpara[1],incpara[2]))/plnorm(14,incpara[1],incpara[2])
  # delay dist
  infvec <- c(0.5,0.5)
  # the delay from inf to report
  w_dis1 <- rep(0,15)
  for (i in 1:15){
    for (j in 1:2){
      if (i+j>0&(i+j<=15)){  
        w_dis1[i+j] <- w_dis1[i+j] + inc_vec[i]*infvec[j]  
      }
    }  
  }
  w_dis1[w_dis1==0] <- 0.001
  w_dis1 <- w_dis1/sum(w_dis1)
  
  input <- as.matrix(d1[,2])
  d2 <- fit_incidence(input,w_dis1)
  
  
  ### after obtain the 1000 estimated time series in step 1 
  ### fit Rt for each estimated time series by EpiEstim
  library(EpiEstim)
  # incubation dist
  mu <- 3.3
  sd <- 2.4
  GI <- c(0,(pgamma(1:15,shape=(mu/sd)^2,scale=sd^2/mu)-pgamma(0:14,shape=(mu/sd)^2,scale=sd^2/mu))/pgamma(14,shape=(mu/sd)^2,scale=sd^2/mu))
  GI <- GI/sum(GI)
  
  Rt_record <- matrix(NA,1000,ncol(d2$Isamps)-7)
  Rt_record_sd <- matrix(NA,1000,ncol(d2$Isamps)-7)
  
  for (i in 1:1000){
    input2 <- data.frame(d1)
    input2[,2] <- d2$Isamps[i,]
    names(input2) <- c("dates","I")
    
    out <- estimate_R(
      input2,
      method = c("non_parametric_si"),
      si_data = NULL,
      si_sample = NULL,
      config = make_config(list(si_distr = GI))
    )
    
    Rt_record[i,] <- out$R$`Mean(R)`
    Rt_record_sd[i,] <- out$R$`Std(R)`
    
  }
  
  ## obtain final estimates
  library(matrixStats)
  Rt_out <- matrix(NA,ncol(d2$Isamps)-7,4)
  Rt_out[,1] <- colMeans(Rt_record)
  Rt_out[,2] <- sqrt(colMeans(Rt_record_sd)^2 + colSds(Rt_record)^2)
  Rt_out[,3] <- Rt_out[,1]-1.96*Rt_out[,2]
  Rt_out[,4] <- Rt_out[,1]+1.96*Rt_out[,2]
  
  
  # #########
  # keep <- Rt_out
  # sigma.test(diff(log(Rt_out[,1])))
  # sigma.test(diff(log(keep[,1])))
  
  
  # ### create new folder and output result
  # city_folder <- paste0('Result/Rt of city/2022_1/',city_name)
  # dir.create(city_folder)
  # save.image(file =paste0('Result/Rt of city/2022_1/',city_name,'/',city_name,'.RData'))
  # write.table(out$R,paste0('Result/Rt of city/2022_1/',city_name,'/',city_name,'.csv'),sep=",",row.names=FALSE)
  # write.table(d1,paste0('Result/Rt of city/2022_1/',city_name,'/',city_name,'_input.csv'),sep=",",row.names=FALSE)

  
}




