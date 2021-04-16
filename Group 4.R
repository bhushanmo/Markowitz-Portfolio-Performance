


library(dplyr)
library(tidyverse)
library(PortfolioAnalytics)
library(IntroCompFinR) # install via: install.packages("IntroCompFinR", repos="http://R-Forge.R-project.org")
library(openair)
library(timetk)
library(plotly)
library(matrixStats)

rf = 0.005/252
r.free <- 0.005/252#252 trading days in 2018
asset.names = c("AAPL", "AHT", "AMZN","GRUB","ATVI","BA",
                "GNE","BCX","CSCO","DB","DIS","F",
                "FB","GM","GMZ","GOOG","AIZ","JNJ",
                "KO","LMT","MA","MCD","MERC","MMM",
                "MRK","MS","MSFT","NB","NFLX","XOM")

#compute the inverse-variance portfolio
getIVP <- function(cov) {
  ivp <- 1/diag(as.matrix(cov))
  ivp <- ivp/sum(ivp)
  return(ivp)
}

#compute variance per cluster
getClusterVar <- function(cov, cItems) {
  #matrix slice
  cov_ <- cov[cItems, cItems]
  w_ <- getIVP(cov_)
  cVar <- t(w_) %*% as.matrix(cov_) %*% w_
  return(cVar)
}

#compute HRP alloc
getRecBipart <- function(covMat, sortIx) {
  return(getRecBipart_help(rep(1,ncol(covMat)), covMat, sortIx))
}

#helper-function for Rekursion => bisection and parse in pairs
getRecBipart_help <- function(w, covMat, sortIx) {
  
  #index vector [1,..., half the lenght of sortIx truncated toward 0]
  subIdx <- 1:trunc(length(sortIx)/2)
  
  #cluster 1
  cItems0 <- sortIx[subIdx]
  #cluster 2
  cItems1 <- sortIx[-subIdx]
  
  #computing cluster variance of sub-clusters
  cVar0 <- getClusterVar(covMat, cItems0)
  cVar1 <- getClusterVar(covMat, cItems1)
  
  alpha <- 1 - cVar0/(cVar0 + cVar1)
  
  #weight 1
  w[cItems0] <- w[cItems0] * alpha
  #weight 2
  w[cItems1] <- w[cItems1] * (1-alpha)
  
  if(length(cItems0) > 1) {
    w <- getRecBipart_help(w, covMat, cItems0)
  }
  
  if(length(cItems1) > 1) {
    w <- getRecBipart_help(w, covMat, cItems1)
  }
  return(w)
}

#creating hierarchical risk parity portfolios:
HRP <- function(stock_data){
  #output W=(w1,...,w30) weights
  
  #covariance matrix
  covMat <- data.frame(t(cov(stock_data[,2:31])))
  #correlation matrix
  corMat <- data.frame(t(cor(stock_data[,2:31])))
  
  #create cluster and get the order
  clustOrder <- hclust(dist(corMat), method = 'single')$order
  
  #return the weights:
  return(getRecBipart(covMat, clustOrder))
}

select_data <- function(data, year){
  
 df <- selectByDate(data.frame(data[1]), year=year)
 df <- select(df, "date", "ret")
 colnames(df) <- c("date", "ret_1")
 
 for (i in 2:30){
   temp<-selectByDate(data.frame(data[i]),year=year)
   temp<-select(temp, "date", "ret")
   ret_i<-paste("ret_", as.character(i), sep="")
   colnames(temp) <- c("date", ret_i)
   df<-merge(x = df, y=temp,by = "date")
 }
 return(df)
}

expected_ret<-function(data){
 exp_ret<-rep(0,30)
 for(i in 1:30){
    exp_ret[i]<-mean(data[,i+1])
 }
 return(exp_ret)
}

cov_mat<-function(data,names){
  covmat <- cov(data[,2:31])
  dimnames(covmat) = list(names, names)
  return(covmat)
} 

random_portf<-function(seed){
  pspec <- portfolio.spec(assets=30, weight_seq=generatesequence())
  pspec <- add.constraint(portfolio=pspec,
                          type="weight_sum",
                          min_sum=0,
                          max_sum=1)
  set.seed(seed)
  return(random_portfolios(portfolio=pspec, permutations = 1001, rp_method = "sample",
                           eliminate = TRUE))
}

MVP<-function(er, covmat){
  #Compute global minimum variance portfolio
  gmin.port = globalMin.portfolio(er, covmat)
  attributes(gmin.port)
  summary(gmin.port, risk.free=r.free)
  plot(gmin.port, col="blue")
  return(gmin.port[["weights"]])
}

plot_Portf<-function(rp,er,year, mvp, hrp,df){
  #Multiply returns of each stock by weighted average of each portfolio
  portfolios<- rp*er[col(rp)]
  #print(er[col(rp)])
  #print(portfolios)
  #Sum the total return of each stock
  returns <- data.frame(rowSums(portfolios))
  
  #Calculate the Standard deviation of each portfolio
  #Create xts object
  df_xts = tk_xts(df,date_col = date)
  for ( i in 1:1000){
    returns$std_dev[i]<-StdDev(df_xts, weights = rp[i,]) 
  }
  #print(returns$std_dev)
  # Rename column to expected return
  returns=returns%>% rename(expected_returns= rowSums.portfolios.)
  
  #calculate sharpe ratio column
  returns$sharpe_ratio = (returns$expected_returns- rf)/returns$std_dev
  
  a <- list(
    x = mvp$std_dev,
    y = mvp$expected_returns,
    text = "MVP",
    xref = "x",
    yref = "y",
    showarrow = TRUE,
    arrowhead = 7,
    ax = 20,
    ay = -40 )
  
  b <- list(
    x = hrp$std_dev,
    y = hrp$expected_returns,
    text = "HRP",
    xref = "x",
    yref = "y",
    showarrow = TRUE,
    arrowhead = 7,
    ax = 20,
    ay = -40 )
  
  x <- list(title = "standard deviation")
  y <- list(title = "expected returns")
  fig <- plot_ly(data = returns,x =returns$std_dev , y = returns$expected_returns,
                 color= returns$sharpe_ratio, type = 'scatter')
  fig <- fig %>% add_markers()
  fig <- fig %>% layout(title=as.character(year), xaxis=x,yaxis=y,annotations = a)
  fig <- fig %>% layout(title=as.character(year), xaxis=x,yaxis=y,annotations = b)
  show(fig)
}


#Obtain data from WRDS to R
res <- dbSendQuery(wrds, "select  b.permno, a.comnam, a.ticker, b.date, b.openprc, b.prc, b.ret
                   from crsp.stocknames a join crsp.dsf b
                   on a.permno = b.permno
                   where b.date between'2017-01-01'
                   and '2018-12-31'
                   and b.permno in (90319,84788,10107,13407,
                   79678,14593,89393,12369,57665,25785,91130,
                   75578,19561,21178,22592,26403,69032,59408,
                   14287,89199,12593,11850,14980,21936,22111,
                   11308,19502,43449,76076,92611,91233,22752,
                   14567,89803,90038,13019)")
data <- dbFetch(res, n = -1)
dbClearResult(res)

#Drop NA Values
data <- drop_na(data)

#Split data into dataframes by its ticker
split_data<-split(data,data$ticker)

#List the data in the global environment
list2env(split_data, envir = .GlobalEnv)

#Cleaning of the data
rm('V','XON','WAG','VWO','SHLX','PFE','NIKE','GOOGL','LK','MERCS','MWD','NCB','NCNB', 'NKE','WBA','AEZNS', 'BAC', 'J', 'ASIAS')
AAPL <- AAPL[AAPL$comnam == 'APPLE INC', ]
ATVI <- ATVI[ATVI$comnam == 'ACTIVISION BLIZZARD INC', ]
BA <- unique(BA)
DIS <- DIS[DIS$comnam == 'DISNEY WALT CO', ]
F <- unique(F)
JNJ <- unique(JNJ)
KO <- unique(KO)
MCD <- unique(MCD)
MS <- MS[MS$comnam == 'MORGAN STANLEY GROUP INC', ]
MMM <- MMM[MMM$comnam == '3M CO', ]
MRK <- MRK[MRK$comnam == 'MERCK & CO INC NEW', ]

df17 <- select_data(list(AAPL, AHT, AMZN,GRUB,ATVI,BA,
                       GNE,BCX,CSCO,DB,DIS,F,
                       FB,GM,GMZ,GOOG,AIZ,JNJ,
                       KO,LMT,MA,MCD,MERC,MMM,
                       MRK,MS,MSFT,NB,NFLX,XOM), "2017")

df18 <- select_data(list(AAPL, AHT, AMZN,GRUB,ATVI,BA,
                         GNE,BCX,CSCO,DB,DIS,F,
                         FB,GM,GMZ,GOOG,AIZ,JNJ,
                         KO,LMT,MA,MCD,MERC,MMM,
                         MRK,MS,MSFT,NB,NFLX,XOM), "2018")

er17<-expected_ret(df17)
names(er17) = asset.names
covmat17<-cov_mat(df17, asset.names)

# Expected returns vector
er18 <- expected_ret(df18)
names(er18) = asset.names
covmat18<-cov_mat(df18, asset.names)

#Calculate the Standard deviation of each portfolio
#Create xts object
df2017_xts = tk_xts(df17,date_col = date)
df2018_xts = tk_xts(df18,date_col = date)

#1000 random Portfolios
#2017
rp17<-random_portf(1)
portfolios17<- rp17*er17[col(rp17)]
#2018
rp18<-random_portf(1)
portfolios18<- rp18*er18[col(rp18)]

#MVP
print(" 2017 MVP portfolio weight: ")
print( MVP(er17, covmat17))
MVP_sd_17 <- StdDev(df2017_xts, weights = MVP(er17, covmat17))
MVP_sd_18 <- StdDev(df2018_xts, weights = MVP(er17, covmat17))
MVP_ret_17<- Return.portfolio(df2017_xts, weights = MVP(er17, covmat17))
MVP_ret_18<- Return.portfolio(df2018_xts, weights = MVP(er17, covmat17))
MVP_sharpe_ratio_17 = (mean(MVP_ret_17$portfolio.returns) - rf)/MVP_sd_17
MVP_sharpe_ratio_18 = (mean(MVP_ret_18$portfolio.returns) - rf)/MVP_sd_18

#HRP
print(" 2017 HRP portfolio weight: " )
print(HRP(df17))
HRP_sd_17 <- StdDev(df2017_xts, weights = HRP(df17))
HRP_sd_18 <- StdDev(df2018_xts, weights = HRP(df17))
HRP_ret_17<- Return.portfolio(df2017_xts, weights = HRP(df17))
HRP_ret_18<- Return.portfolio(df2018_xts, weights = HRP(df17))
HRP_sharpe_ratio_17 = (mean(HRP_ret_17$portfolio.returns) - rf)/HRP_sd_17
HRP_sharpe_ratio_18 = (mean(HRP_ret_18$portfolio.returns) - rf)/HRP_sd_18

# MVP 2017 results

print(paste(" 2017 MVP portfolio expected return in 2017:",mean(MVP_ret_17$portfolio.returns)))
print(paste(" 2017 MVP portfolio standard deviation in 2017 :",MVP_sd_17))
print(paste(" 2017 MVP portfolio sharpe ratio in 2017 :",MVP_sharpe_ratio_17))

#out of sample(2018) results for 2017 MVP portfolio

print(paste(" 2017 MVP portfolio expected return in 2018:",mean(MVP_ret_18$portfolio.returns)))
print(paste(" 2017 MVP portfolio standard deviation in 2018 :",MVP_sd_18))
print(paste(" 2017 MVP portfolio sharpe ratio in 2018 :",MVP_sharpe_ratio_18))

# HRP 2017 results

print(paste(" 2017 HRP portfolio expected return in 2017:",mean(HRP_ret_17$portfolio.returns)))
print(paste(" 2017 HRP portfolio standard deviation in 2017 :",HRP_sd_17))
print(paste(" 2017 HRP portfolio sharpe ratio in 2017 :",HRP_sharpe_ratio_17))

#out of sample(2018) results for 2017 HRP portfolio

print(paste(" 2017 HRP portfolio expected return in 2018:",mean(HRP_ret_18$portfolio.returns)))
print(paste(" 2017 HRP portfolio standard deviation in 2018 :",HRP_sd_18))
print(paste(" 2017 HRP portfolio sharpe ratio in 2018 :",HRP_sharpe_ratio_18))

#Plot 2017
plot_Portf(rp17 ,er17,2017, 
           data.frame("expected_returns" = mean(MVP_ret_17$portfolio.returns) , "std_dev" =MVP_sd_17, "sharpe_ratio" =MVP_sharpe_ratio_17), 
           data.frame("expected_returns" = mean(HRP_ret_17$portfolio.returns) , "std_dev" =HRP_sd_17,  "sharpe_ratio" =HRP_sharpe_ratio_17),df17)

#Plot 2018
plot_Portf(rp18,er18,2018, 
           data.frame("expected_returns" = mean(MVP_ret_18$portfolio.returns) , "std_dev" =MVP_sd_18, "sharpe_ratio" =MVP_sharpe_ratio_18), 
           data.frame("expected_returns" = mean(HRP_ret_18$portfolio.returns) , "std_dev" =HRP_sd_18,  "sharpe_ratio" =HRP_sharpe_ratio_18),df18)

