#Written by Vishwesha Guttal, Centre for Ecological Sciences, Indian Institute of Science, Bangalore, 560012, India.
#Edited by Nikunj Goel, Undergraduate Department, Indian Institute of Science, Bangalore, 560012, India.
# This code generates figure 3 equivalent of the main text and a pvalue table as an output in the Console

# Please read "Instructions" file.

# Note: Please make sure the following packages are installed before running the code
library('moments')
library('Kendall')
library('KernSmooth')

#Compile the function file.
source("functions_stock_analysis.R")

# Load data (DJI, SP500 and NASDAQ)
#d2=read.csv("dow jones index per day data.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)
#d2=read.csv("sp500 per day data.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)
d2=read.csv("nasdaq per day data.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)

d2$Date = as.Date(d2$Date,"%d/%m/%Y")

# Input a date prior to which indicator trends have to be calculated.
# Please note that the input date is not one of the days when the stock maket was closed.
# Stock market data for DJI, SP500 and NASDAQ is only available after 1896, 1928 and 1971 respectively. 
# stdcrashdate = c("03/09/1929");
# stdcrashdate = c("25/08/1987");
 stdcrashdate = c("14/01/2000");
# stdcrashdate = c("02/5/2008");
stdcrashdate = as.Date(stdcrashdate,"%d/%m/%Y");
crashyear = as.numeric(format(stdcrashdate,"%Y"))

# Parameter values to calculate indicator trends
# Note: choose parameter values such that, N > l_rw + l_ktau
l_rw=500 # Length of rolling window. Default: 500 (2 years)
bw=25 # Bandwith. Default: 25
N=1000; #length of stock-index to be analyzed. Default: 1000 (4 years)
l_ktau=250; #length over which kendall is estimated. Default: 250 (1 year)
N_realization=1000; # To generate pseudo time series to calculate pvalues. Default: 1000 


# Range of parameters to evaluate kendall-tau histograms.
# Same as in mauscript. Refer to sensitivity analysis in methods
# Note: choose parameter values such that, max(N_sensitivity) > max(rwrange) + max(l_kw) + max(k_end) 
rwrange=c( seq (375,625,by=25));# Length of rolling window. Default: from 375 to 625 in steps of 25
bwrange=c(seq (2.5,100,by=2.5));# Bandwith. Default: from 2.5 to 100 in steps of 2.5
N_sensitivity=1500; #length of stock-index to be analyzed Default: 1500 (6 years)
l_kw= c(seq(175,325,by = 5)) #length over which kendall is estimated. Default: from 175 to 325 in steps of 5
k_end = c(seq(0,200,by = 5))#Stand/end point of the kendall window. Default: from 0 to 200 in steps of 5

#ETA: 1 to 1.5 hrs depending on the machine, parameter values, input date and length of stock index time series
begtime = Sys.time(); #To keep track of how much it takes.

#A function to obtain obtain residuals and smoothened data for dow-precrash data
stock_precrash = precrashdata(stdcrashdate,d2,bw,N);

#A function to calculate indicator trends.
ews_trends = ews(stock_precrash$residuals, l_rw)

#A function to to calculate Kendall-tau coefficients.
kendalls = kendall_coefficient(ews_trends,l_ktau,0,N)

# Sensitivity Analysis (Histograms)
kendall_histogram = sensitivity_histograms(stdcrashdate,d2,rwrange,bwrange,N_sensitivity,l_kw,k_end)

#function to calculate pvalues
pvalues=significance(stdcrashdate,d2, stock_precrash$residuals, l_rw, l_ktau, kendalls, N_realization,bw)

#function to print pvalues and kendall values in a table format
pvalue_table(kendalls,pvalues,N)

# Plotting trends of indicators prior to a crash
plotting_graphs(stock_precrash,ews_trends,kendall_histogram)

# Estimating time taken to run the code
print("Time to run all null models and sensitivity analysis is given below")
print(Sys.time()-begtime)
