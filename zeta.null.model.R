##############################################################
#
#  Tests zeta decay by removing occurrences, but with a bias parameter to bias towards
# rare/common species
#
##############################################################



library(zetadiv)
library(vegan)
library(picante)



###Maximum number of transects to consider
max_t<-40

###Number of sample to average over
samp_s<-100


gombe<-read.csv("Gombe_138_71.csv", sep=",", header=T)

gombe<-read.csv("MahaleN_56_71.csv", sep=",", header=T)

	k<-read.csv("Kalilani_56_71.csv", sep=",", header=T)
	k2<-ifelse(k>0, 1, 0)
	k2<-as.data.frame(k2)

g2<-ifelse(gombe>0, 1, 0)
g2<-as.data.frame(g2)

g2<-g2[, colSums(g2)>0]
sp_occ<-colSums(g2)


#Get the number of species
nc<-ncol(g2)

#Get the number of surveys
nr<-nrow(g2)

#sp_loss sets the number of species occurrences to remove
sp_loss<-507

#bias determines the bias towards rare (bias<0) or common (bias>0) species
bias<-  -1.2

fname<-file("zeta.lossM.txt", open="a")

#cat(file=fname, "bias", "loss", "exp.aic", "exp.int", "exp.c", "pl.aic", "pl.int", "pl.c", sep="\t", fill=T)
#close(fname)

#The number of Monte Carlo simulations to run for each parameter set
iter=100

sum<-0

g3<-g2
total<-sum(g3)

#Get the occurrence frequency for each species
sp_freq<-(colSums(g3))

#Now weight the frequency by the bias parameter
sp_freq<-sp_freq^bias
sp_freq<-sp_freq

sp_freq<-as.vector(sp_freq)
sp_int<-cumsum(sp_freq)



##We will loop through a series of bias parameters

for(bias in seq(-2.0,-1, 0.1))
{
sum<-0
for(m in 1 : iter)
	{

	#Start each run from the pristine empirical dataset
	g3<-g2
	total<-sum(g3)

	sp_freq<-(colSums(g3))
	sp_freq<-sp_freq^bias
	sp_freq<-sp_freq

	sp_freq<-as.vector(sp_freq)
	sp_int<-cumsum(sp_freq)

	n<-0

	#Now choose a species at random, taking into account any bias towards rare/common species
	while(n<sp_loss)
		{

		rn<-runif(1, 0, sp_int[nc])

		for(j in 1:nc)
			if(sp_int[j]>rn)
				{
				ec=j
				break
				}

		flag<-0

	#Having chosen a species, now choose a survey occurrence to remove. All occurrences are equally likely to be chosen
		while(flag<1)
			{		
			
			er<-round(runif(1, 1,nr))
			if(sum(g3[,ec])<1)
				{
				flag=1
				}
			else if(g3[er,ec]==1)
				{
				g3[er,ec]=0
				n=n+1
				flag=1
				}
			}
		}	

#Having removed the necessary number of occurrences we now run the stnadard zeta analyses
#We record the AIC values for the power law and exponential model fits, and also the parameter values for the respective models

	gz<-Zeta.decline(g3, orders=1:max_t, sam=samp_s)
	gz$aic

	#Record how many times the exponential model is selected
	if(gz$aic[1,2]<gz$aic[2,2]){
	sum=sum+1
	}
#Now write the results to file
	cat(bias, sp_loss, gz$aic[1,2], as.vector(gz$zeta.exp$coefficients[1]), as.vector(gz$zeta.exp$coefficients[2]), gz$aic[2,2],  as.vector(gz$zeta.pl$coefficients[1]), as.vector(gz$zeta.pl$coefficients[2]), file=fname, sep="\t", fill=T, append=T)
	#close(fname)
		
	
#Print proportion of simulations that lead to exponential distribution selection
	cat(sum/m, "\n")
}
cat(bias, sp_loss, sum/m, file="zeta.loss.aicM.txt", sep="\t", fill=T, append=T)
}


sum/iter