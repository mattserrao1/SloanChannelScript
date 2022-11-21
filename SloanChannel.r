#Package install and call
install.packages("dataRetrieval", 
                 repos=c("http://owi.usgs.gov/R",
                         getOption("repos")))
						 
library(zoo)
#library(raster) not needed, I just had added for something else
#library(rgdal) not needed, I just had added for something else
library(hydroTSM)
library(Kendall)
library(dataRetrieval)
library(measurements)
##Let's Pull Discharge, "Q" (00060);And Rainfall, "P" (00045)
site.code = "09419665" #sitecode for sloan channel via USGS
readNWISsite(site.code)  
what.data = whatNWISdata(siteNumber = site.code)
start.date = ""  # Blanks get all of the data, but if a specific range is needed enter here 
end.date = ""
x2 <- readNWISdv(site.code,c('00045','00060'),
                                       start.date, end.date, statCd=c('00006','00003'))
names(x2)[c(4,5,6,7)] = c("P.in","PA","Q.ft.s","QA")##rename
breaks <- seq(from=as.Date("1988-10-01"), to=as.Date("2022-09-30"), by="year") ##setting Oct-Sept year for region
years.breaks = as.numeric(format(breaks,"%Y"))
labels.wy = years.breaks[2:length(breaks)] ##creates label.wy 
x2$WaterYear <- cut(x2$Date, breaks,labels=labels.wy) ##adds water year column
fannual.mean.cfs = aggregate(x2$Q.ft.s,by=list(x2$WaterYear),FUN=mean) ##Creates annual mean cfs
fannual.mean.pin = aggregate(x2$P.in,by=list(x2$WaterYear),FUN=sum) ##Creates annual sum rainfall 
names(fannual.mean.cfs) = c("WYear","MeanQcfs")
names(fannual.mean.pin) = c("WYear","MeanPin")
fannual.mean.cfs$WYear = as.numeric(as.character(fannual.mean.cfs$WYear))
fannual.mean.cfs$MeanQcfs = as.numeric(as.character(fannual.mean.cfs$MeanQcfs))
fannual.mean.pin$WYear = as.numeric(as.character(fannual.mean.pin$WYear))
fannual.mean.pin$MeanPin = as.numeric(as.character(fannual.mean.pin$MeanPin))
class(fannual.mean.cfs$WYear) ##both of these are to check the class of the variable is numeric
class(fannual.mean.cfs$MeanQcfs)
#Convert QF/S to mm/yr, and inches to mm/yr
site.data = readNWISsite(site.code) ##setting variable site.data as the site data from the readNWISSite function, for the variable site.code that we set earlier
drain.area.mi2 = site.data$drain_area_va ##setting the variable drain.area.mi2 as the data pulled from variable site.data sub drain area
drain.area.mm2= conv_unit(drain.area.mi2, "mi2", "mm2")
print(drain.area.mm2)
ft3tommyr = 8.932128e+14/drain.area.mm2 #8.93212 is ft3/s to mm3/yr
fannual.mean.cfs$mmyr = fannual.mean.cfs$MeanQcfs/ft3tommyr
fannual.mean.pin$mmyr = fannual.mean.pin$MeanPin*25.4
####Now we have P in mm/yr and Q in mm/yr, so we can do trend tests but we need Q/P and to clean up the data
##First, let's find our runoff coefficient for each year, or Q/P
QPSloan = fannual.mean.cfs$mmyr/fannual.mean.pin$mmyr
fannual.mean.cfs$QP = QPSloan
fannual.mean.cfs$pmm = fannual.mean.pin$mmyr
completedtable = na.omit(fannual.mean.cfs)
completedtable1 = completedtable[-c(16,17), ]
##LM & Kendall for Q/P, Q, and P
#Q Linear Test, Mann Kendall, and Residual Plot
lm.sloan.annualq = lm(completedtable1$mmyr ~ completedtable1$WYear) 
km.sloan.annualq = MannKendall(completedtable1$mmyr)
fannualqstd = rstandard(lm.sloan.annualq)
summary(km.sloan.annualq)
summary(lm.sloan.annualq)
qqnorm(fannualqstd,ylab="Standardized Residuals",xlab="Normal Scores",main="Annual Q Residuals")
qqline(fannualqstd)
#P Linear Test, Mann Kendall, and Residual Plot
lm.sloan.annualp = lm(completedtable1$pmm ~ completedtable1$WYear) 
km.sloan.annualp = MannKendall(completedtable1$pmm)
fannualpstd = rstandard(lm.sloan.annualp)
summary(km.sloan.annualp)
summary(lm.sloan.annualp)
qqnorm(fannualpstd,ylab="Standardized Residuals",xlab="Normal Scores",main="Annual P Residuals")
qqline(fannualpstd)
#Q/P Linear Test
lm.sloan.annualqp = lm(completedtable1$QP ~ completedtable1$WYear) 
km.sloan.annualqp = MannKendall(completedtable1$QP)
fannualqpstd = rstandard(lm.sloan.annualqp)
summary(km.sloan.annualqp)
summary(lm.sloan.annualqp)
qqnorm(fannualqpstd,ylab="Standardized Residuals",xlab="Normal Scores",main="Annual Q/P Residuals")
qqline(fannualqpstd)
##3 panel residual
par(mfrow=c(2,2))
qqnorm(fannualqstd,ylab="Standardized Residuals",xlab="Normal Scores",main="Annual Q Residuals")
qqline(fannualqstd)
qqnorm(fannualpstd,ylab="Standardized Residuals",xlab="Normal Scores",main="Annual P Residuals")
qqline(fannualpstd)
qqnorm(fannualqpstd,ylab="Standardized Residuals",xlab="Normal Scores",main="Annual Q/P Residuals")
qqline(fannualqpstd)
##Plotting Annual P,Q,and Q/P
#3 Panel Plots Showing Annual P, Annual Q, and Q/P for the study period
par(mfrow=c(3,1))
par(mar=c(2,2,2,1),oma=c(3,3,1,1))#creates multipanel plot, c(3,1) makes the plot 3 row 1 column #mar is the margin between the plots, oma is outer margin
#Panel 1
plot(completedtable1$WYear,completedtable1$pmm,type="l")
axis(side=2,labels=FALSE)
mtext("Annual Precip,mm",side=2, line=2)
abline(lm.sloan.annualp,col="red")
#Panel 2
plot(completedtable1$WYear, completedtable1$mmyr,type="l")
axis(side=2,labels=FALSE)
abline(lm.sloan.annualq,col="red")
mtext("Annual Runoff,mm",side=2,line=2)
#Panel 3
plot(completedtable1$WYear,completedtable1$QP,xlab="Year",type="l")
abline(lm.sloan.annualqp,col="red")
mtext("Q/P ,mm",side=2,line=2)
##Title & add a label for the bottom 
mtext("Sloan Channel Watershed Precip & Runoff, 1989-2019", side=3,line=-2,outer=TRUE)
mtext("Water Year", side=1,line=0,outer=TRUE)
mtext("By:M.Serrao", side=1,line=1,at=0,cex=0.5,outer=TRUE,font=10)
legend("bottomright",c("Regression Line"), col=c("red"), pch=95)
#write.csv(completedtable1,file="sloandata1.csv") #outputs csv for tables
##Discharge over 1
#Separate out years
rainfall.2006wy = x2[x2$Date %in% seq(as.Date("2005-10-01"),as.Date("2006-09-30"),by="day"),] #Pulls all the data we need
rainfall.1989wy = x2[x2$Date %in% seq(as.Date("1988-10-01"),as.Date("1989-09-30"),by="day"),] #Pulls all the data we need
rainfall.2017wy = x2[x2$Date %in% seq(as.Date("2016-10-01"),as.Date("2017-09-30"),by="day"),] #Pulls all the data we need
##Flood Flow for 2017, regardless of rainfall due to flash-flood processes just in case:), do for 1989 and for 2006
#2017
discharge.over.1 = rainfall.2017wy[rainfall.2017wy$Q.ft.s> 1,] ## 0.8255 is the mean discharge for 2017 WY
floodindex = c(0,cumsum(diff(discharge.over.1$Date)>1))
#2006
discharge2006.over.1 = rainfall.2006wy[rainfall.2006wy$Q.ft.s>1,] ##0.7324 is the mean discharge of 2006 WY
floodindex2006 = c(0,cumsum(diff(discharge2006.over.1$Date)>1)) #not needed unless doing stormflow (P v Q, mm for x events)
#1989
discharge1989.over.1 = rainfall.1989wy[rainfall.1989wy$Q.ft.s>1,] ##0.0742 is mean discharge of 1989 WY
floodindex1989 = c(0,cumsum(diff(discharge1989.over.1$Date)>1))
#Triple Plot
par(mfrow=c(3,1))
par(mar=c(2,2,2,1),oma=c(3,3,1,1))
plot(discharge1989.over.1$Date,discharge1989.over.1$Q.ft.s,xlab="Date",ylab="Discharge, Qf/s",pch=19,col=460)
mtext("Sloan Channel Discharge > 1", side=3,line=-2,outer=TRUE)
mtext("1989",side=2,line=2,cex=0.85)
plot(discharge2006.over.1$Date,discharge2006.over.1$Q.ft.s,xlab="Date",ylab="Discharge, Qf/s",pch=19,col=460)
mtext("Discharge, Qft/s",side=2,line=3.5)
mtext("2006",side=2,line=2,cex=0.85)
plot(discharge.over.1$Date,discharge.over.1$Q.ft.s,xlab="Date",ylab="Discharge, Qf/s",pch=19,col=460)
mtext("2017",side=2,line=2,cex=0.85)
#FDC Curves and stuff
Q.1989.1999 = x2[(x2$Date >= as.Date("1988-10-01")) & (x2$Date <= as.Date("1999-09-30")),]
Q.2000.2010 = x2[(x2$Date >= as.Date("1999-10-01")) & (x2$Date <= as.Date("2010-09-30")),]
Q.2011.2019 = x2[(x2$Date >= as.Date("2010-10-01")) & (x2$Date <= as.Date("2019-09-30")),]
Q.2011.2015 = x2[(x2$Date >= as.Date("2010-10-01")) & (x2$Date <= as.Date("2015-04-03")),]
Q.2015.2019 = x2[(x2$Date >= as.Date("2015-04-04")) & (x2$Date <= as.Date("2019-10-01")),]
#Plot Daily Discharge, Colored by Break
begin.end = range(x2$Date)  # Get the date range for the whole dataset
plot(as.Date(Q.1989.1999$Date),Q.1989.1999$Q.ft.s,xlim=as.Date(c("1988-10-01","2019-09-30")),ylab="Q ft3/s",type="l",xlab="")
lines(Q.2000.2010$Date,Q.2000.2010$Q.ft.s,col="red")
lines(Q.2011.2019$Date,Q.2011.2019$Q.ft.s,col="blue")
lines(Q.2011.2015$Date,Q.2011.2015$Q.ft.s,col=460)
mtext("Sloan Channel Daily Discharge 1989-2019 Water Years", side=3,line=-2,outer=TRUE)
legend("topleft",c("1989-1999","2000-2010","2011-2015","2015-2019"),col=c("black","red",460,"blue"),lty=c(1,1,1))
#Flow Duration Curves
fdc.1989.1999 = fdc(Q.1989.1999$Q.ft.s,new=TRUE,thr.shw=FALSE,xlim=c(0,1),ylab="Q ft3/s")
fdc.2000.2010 = fdc(Q.2000.2010$Q.ft.s,new=FALSE,thr.shw=FALSE,col="red")
fdc.2011.2019 = fdc(Q.2011.2019$Q.ft.s,new=FALSE,thr.shw=FALSE,col="blue")
fdc.2015.2019 = fdc(Q.2015.2019$Q.ft.s,new=FALSE,thr.shw=FALSE,col=460)
legend("topright",c("1989-1999","2000-2010","2011-2019", "2015-2019"),col=c("black","red","blue",460),lty=c(1,1,1))
plot(completedtable$pmm,completedtable$Qmm, xlab="Rainfall, mm",ylab="Discharge, mm",pch=19)
mtext("Sloan Channel Stormflow Analysis 1989-2019 Water Years", side=3,line=-2,outer=TRUE)
#For finding mean
Q.2000.2010fixed = na.omit(Q.2000.2010)
