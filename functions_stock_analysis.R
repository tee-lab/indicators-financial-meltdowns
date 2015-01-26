#Contains all functions to analyze time series of stock markets. 

precrashdata <- function(stdcrashdate,d2,bw,N)
{
  
  stdindex = which(d2$Date == stdcrashdate)
  indices_precrash = seq(stdindex-N+1, stdindex);
  precrash = d2$value[indices_precrash]
  smoothdata=ksmooth(1:N,precrash,kernel=c("box","normal"),bandwidth=bw) # To calculate deterministic trends in time series of stock market
  smooth=smoothdata$y
  residuals=precrash-smooth # Residual data
  
  return(list(dates=d2$Date[indices_precrash],precrash=precrash,smooth=smooth,residuals=residuals))
  
}

ews <- function(residuals, l_rw) 
{
  
  #Declare arrays   
  N=length(residuals)
  var_residuals=numeric(N)
  acf_residuals=numeric(N)
  avgspec_residuals=numeric(N)
  
  #Apply rolling window to residuals.
  #Residual data analysis.
  for (i in 1:(N-l_rw+1))
  {
    rolldata = residuals[i:(i+l_rw-1)];
    
    var_residuals[i+l_rw-1] = var(rolldata);
    acf_residuals[i+l_rw-1] = acf(rolldata,plot=FALSE)$acf[2];
    
    spec_residuals = spectrum(rolldata,plot=FALSE)$spec;
    avgspec_residuals[i+l_rw-1] = mean(spec_residuals[2:floor(l_rw/8)]);
  }
  acf_residuals[1:l_rw]= NA
  var_residuals[1:l_rw]= NA
  avgspec_residuals[1:l_rw]= NA
  
  #Output results
  ews_trends = list(acf_residuals=acf_residuals, var_residuals=var_residuals, spec_residuals=avgspec_residuals)
  return(ews_trends)
  
}

kendall_coefficient <- function(ews_trends,l_kw,k_end,N)
{
  
  init_index= N - l_kw; 
  kendaltau_acf = Kendall(1:l_kw, ews_trends$acf_residuals[(init_index-k_end):(N-1-k_end)])$tau[1];
  kendaltau_var = Kendall(1:l_kw, ews_trends$var_residuals[(init_index-k_end):(N-1-k_end)])$tau[1];
  kendaltau_avgspec = Kendall(1:l_kw, ews_trends$spec_residuals[(init_index-k_end):(N-1-k_end)])$tau[1];
  kendalls = list(acf=kendaltau_acf,  var=kendaltau_var, spec=kendaltau_avgspec)
  return(kendalls)
  
  # All the kendall-tau values equal to one have been made to be 0.999 because of the error
  # it produces while generating histograms
  
}

sensitivity_histograms <- function(stdcrashdate,d2,rwrange,bwrange,N_sensitivity,l_kw,k_end)
{
  z = 0;
  allkendalls=array(0,dim=c(length(bwrange)*length(rwrange)*length(k_end)*length(l_kw),7))
  kendall_acf = array(0, dim=c(length(rwrange),length(bwrange)))
  kendall_var = array(0, dim=c(length(rwrange),length(bwrange)))
  kendall_avgspec = array(0, dim=c(length(rwrange),length(bwrange)))

  for (rw in rwrange)
  {
    
    for (bw in bwrange)
    {
      
      stock_precrash = precrashdata(stdcrashdate,d2,bw,N_sensitivity)
      ews_trends = ews(stock_precrash$residuals, l_rw)
      
      for (kw in l_kw)
      {
        
        for ( end  in k_end)
        { 
          
          z = z+1;
          kendalls = kendall_coefficient(ews_trends,kw,end,N_sensitivity)
          allkendalls[z,] = c(rw, bw,kw,end, kendalls$acf, kendalls$var, kendalls$spec);
          
        } 
      } 
    } 
  }
  
  return(list(acf =  allkendalls[,5],var =  allkendalls[,6],spec =  allkendalls[,7]))
  
}

plotting_graphs <- function(stock_precrash,ews_trends,kendall_histogram)
  
{ 
  
  # Creating a 5x1 pannel. Equivalent of figure 3 in main text.
  par(mfrow=c(4,2),mai=c(.8,.9,0.35,0.05))
  plot(stock_precrash$dates,stock_precrash$precrash,type='l',lwd = 5,xlab = 'Date',ylab = 'Stock Value')
  points(stock_precrash$dates,stock_precrash$smooth,type='l',lwd = 3,col='red')
  plot(stock_precrash$dates,stock_precrash$residuals,type='l',lwd = 5,xlab = 'Date',ylab = 'Residuals',
       col='blue',)
  plot(stock_precrash$dates,ews_trends$acf_residuals,type='l',lwd = 5,xlab = 'Date',ylab = 'Autocorrelation \n at Lag 1',
       col='green')
  h = hist(kendall_histogram$acf,breaks = seq(-1.05,1.05, by = 0.05),plot=F)
  h$counts = h$counts/sum(h$counts)
  plot(h, col="grey",xlab = 'Kendall-Tau',ylab = "Normalized \n Frequency",xlim = c(-1,1),
       main="Kendall-Tau histogram \n for Autocorrelation at lag one")
  plot(stock_precrash$dates,ews_trends$var_residuals,type='l',lwd = 5,xlab= 'Date',ylab = 'Variance',
       col='chartreuse4')
  h = hist(kendall_histogram$var,breaks = seq(-1.05,1.05, by = 0.05),plot=F);
  h$counts = h$counts/sum(h$counts)
  plot(h, col="grey",xlab = 'Kendall-Tau',ylab = "Normalized \n Frequency",xlim = c(-1,1),
       main="Kendall-Tau histogram \n for Variance")
  plot(stock_precrash$dates,ews_trends$spec_residuals,type='l',lwd = 5,xlab = 'Date', ylab = 'Mean Power Spectrum \n at Low Frequencies',
       col='deeppink')
  h = hist(kendall_histogram$spec,breaks = seq(-1.05,1.05, by = 0.05),plot=F);
  h$counts = h$counts/sum(h$counts)
  plot(h, col="grey",xlab = 'Kendall-Tau',ylab = "Normalized \n Frequency",xlim = c(-1,1),
       main="Kendall-Tau histogram \n for Mean Power Spectrum")

}

hyst <- function (stdcrashdate,d2, l_rw, l_ktau, kendalls,N,bw)
{
  years = as.numeric(format(d2$Date,"%Y"))
  crashyear = as.numeric(format(stdcrashdate,"%Y"))
  
  startindex = 1;
  endindex = max(which(years< floor(crashyear-N/250)))
  indices_inc = 5;

  Nts = floor((endindex - startindex - N + 1)/indices_inc) + 1
  istart = startindex - indices_inc  
  
  n = Nts
  allkendall = array(0,dim=c(n,3)) #store n number of 3 kendall's. 
  
  for (i in 1:n)
  {
    istart = istart + indices_inc
    indices = seq(istart, istart+N-1);
    stocks = d2$value[indices];
  
    #Detrending the data
    smoothdata=ksmooth(indices,stocks,kernel=c("box","normal"),bandwidth=bw)
    smooth=smoothdata$y
  
    #Find residual data
    residuals=stocks-smooth
    
    #Find historical Kendalls
    pseudo_kendall = kendall_coefficient(ews(residuals, l_rw),l_ktau,0,N)
    allkendall[i,] = c(pseudo_kendall$acf,pseudo_kendall$var,pseudo_kendall$spec)
    
  }  
  
  pacf = length(which(allkendall[,1] >= kendalls$acf))/n
  pvar = length(which(allkendall[,2] >= kendalls$var))/n
  pspec = length(which(allkendall[,3] >= kendalls$spec))/n
  
  pvalues = list(pacf=pacf, pvar=pvar, pspec=pspec)
  
}

spec <- function(residuals, l_rw, l_ktau, kendalls, N_realization)
{
  N=length(residuals)
  n=N_realization #number of realizations
  
  allkendall = array(0,dim=c(n,3)) #store n number of 3 kendall's. 
  
  for (i in 1:n) {
    
    resample_residuals = numeric(N) 
    
    fftresiduals = complex(N)
    fftresiduals_shift = complex(N)
    randphase = numeric(N)
    randmult = complex(N)
    
    fftresiduals = fft(residuals)
    
    #Create random phases.
    randphase[1]=0
    randphase[N/2+1]=0
    
    for (m in 2:(N/2))
    {
      randphase[m] = runif(1)*2*pi;
      randphase[N-m+2] = -randphase[m]; 
      #Above indexing done after careful checking of indices and whether final answer is real.
    }
    
    randmult = complex(real = cos(randphase), imaginary = sin(randphase))
    
    #Shifting fft by a random phase.
    fftresiduals_shift = fftresiduals*randmult;
    resample_residuals = fft(fftresiduals_shift,inverse=TRUE)/length(fftresiduals);
    
    pseudo_kendall = kendall_coefficient(ews(Re(resample_residuals), l_rw),l_ktau,0,N)
    allkendall[i,] = c(pseudo_kendall$acf,pseudo_kendall$var,pseudo_kendall$spec)

  }
  
  pacf = length(which(allkendall[,1] >= kendalls$acf))/n
  pvar = length(which(allkendall[,2] >= kendalls$var))/n
  pspec = length(which(allkendall[,3] >= kendalls$spec))/n
  pvalues = list(pacf=pacf, pvar=pvar, pspec=pspec)
  
  return(pvalues)
}

ar1 <- function(residuals, l_rw, l_ktau, kendalls, N_realization)
{
  N=length(residuals)
  n=N_realization #number of realizations
  
  allkendall = array(0,dim=c(n,3)) #store n number of 3 kendall's. 
  
  for (i in 1:n) {
    
    resample_residuals = numeric(N) 
    
    mu=mean(residuals);
    nu=var(residuals);
    acf1=acf(residuals,plot=FALSE)$acf[2]; #acf at lag1.
    
    a1=acf1;
    a0=mu*(1-acf1);
    sigma=sqrt(nu*(1-acf1^2));
    resample_residuals[1]=residuals[1];
    for (j in 2:N)
    {
      resample_residuals[j] = a0 + a1*resample_residuals[j-1] + sigma*rnorm(1);
    }
    
    pseudo_kendall = kendall_coefficient(ews(resample_residuals, l_rw),l_ktau,0,N)
    allkendall[i,] = c(pseudo_kendall$acf,pseudo_kendall$var,pseudo_kendall$spec)
  }
  
  pacf = length(which(allkendall[,1] >= kendalls$acf))/n
  pvar = length(which(allkendall[,2] >= kendalls$var))/n
  pspec = length(which(allkendall[,3] >= kendalls$spec))/n
  
  pvalues = list(pacf=pacf, pvar=pvar, pspec=pspec)
  
  return(pvalues)
}

bootstrap <- function(residuals, l_rw, l_ktau, kendalls, N_realization)
{
  N=length(residuals)
  n=N_realization #number of realizations
  
  allkendall = array(0,dim=c(n,3)) #store n number of 3 kendall's. 
  
  for (i in 1:n) {
    
    resample_residuals = numeric(N) 
    
    for (j in 1:N)
    {
      randindex = floor(runif(1)*N)+1;
      resample_residuals[j] = residuals[randindex];  	
    }

    pseudo_kendall = kendall_coefficient(ews(resample_residuals, l_rw),l_ktau,0,N)
    allkendall[i,] = c(pseudo_kendall$acf,pseudo_kendall$var,pseudo_kendall$spec)
  }
  
  pacf = length(which(allkendall[,1] >= kendalls$acf))/n
  pvar = length(which(allkendall[,2] >= kendalls$var))/n
  pspec = length(which(allkendall[,3] >= kendalls$spec))/n
  
  pvalues = list(pacf=pacf, pvar=pvar, pspec=pspec)
  
  return(pvalues)
}

significance <- function (stdcrashdate,d2, residuals, l_rw, l_ktau, kendalls, N_realization,bw)
{

  phist = hyst(stdcrashdate,d2, l_rw, l_ktau, kendalls,length(residuals),bw)
  pspec = spec(residuals, l_rw, l_ktau, kendalls, N_realization)
  par1 = ar1(residuals, l_rw, l_ktau, kendalls, N_realization)
  pboot = bootstrap(residuals, l_rw, l_ktau, kendalls , N_realization)
  
  return(list(phist=phist, pboot=pboot, par1=par1, pspec=pspec))
}

pvalue_table <- function(kendalls,pvalues,N)
{
  Table = array(0,dim=c(5,3));
  Table[1,1] = kendalls$acf;
  Table[1,2] = kendalls$var;
  Table[1,3] = kendalls$spec;

  Table[2,1] = pvalues$phist$pacf;
  Table[2,2] = pvalues$phist$pvar;
  Table[2,3] = pvalues$phist$pspec;
  
  Table[3,1] = pvalues$pspec$pacf;
  Table[3,2] = pvalues$pspec$pvar;
  Table[4,3] = pvalues$pspec$pspec;
  
  Table[4,1] = pvalues$par1$pacf;
  Table[4,2] = pvalues$par1$pvar;
  Table[4,3] = pvalues$par1$pspec;
  
  Table[5,1] = pvalues$pboot$pacf;
  Table[5,2] = pvalues$pboot$pvar;
  Table[5,3] = pvalues$pboot$pspec;
  
  # We cannot estimate pvalues smaller than 1/N because of the finite number of realizations
  Table[which(Table == 0)] = (1/N)
  
  Table = as.data.frame(Table)
  colnames(Table) <- c("Autocorrelation","      Variance"," Power Spectrum")
  rownames(Table) <- c("Kendall-Tau","      phist","      pspec","       par1","      pboot")
  print("Kendall-Tau coefficents and pvalues for")
  print(stdcrashdate)
  print("are given in the table below")
  print(Table)
}

