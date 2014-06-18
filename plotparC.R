library(ggplot2); library(gridExtra);library(plyr);library(GGally)
params = c(paste("theta[", 1:8,"]", sep=""), "theta[1]/theta[2]", "theta[5]/theta[6]")
paramsfn = c(paste("theta", 1:8, sep=""), "theta1bytheta2", "theta5bytheta6")

NP = 10
NUM = 2

if (Sys.info()["sysname"] == "Darwin") setwd("~/Dropbox/GPC/working document/sim/prog")

sim_id <- "CPAR100000lin49x20T49.0TC"
thetas_all20 <- read.csv(paste(sim_id,"/thetas",sim_id,".csv",sep=""),header=FALSE)
thetas_all20 = cbind(thetas_all20, V12=thetas_all20[,1]/ thetas_all20[,2], V56=thetas_all20[,5]/ thetas_all20[,6])
sim_id <- "CPAR100000lin49x50T49.0TC"
thetas_all50 <- read.csv(paste(sim_id,"/thetas",sim_id,".csv",sep=""),header=FALSE)
thetas_all50 = cbind(thetas_all50,V12= thetas_all50[,1]/ thetas_all50[,2], V56= thetas_all50[,5]/ thetas_all50[,6])
thetas_all <- rbind(thetas_all20,thetas_all50)
colnames(thetas_all) <- params
d_all <- stack(thetas_all)
head(d_all)
N <- nrow(thetas_all20)
thetas_all$m <- rep(c('20','50'),each=N)
colnames(thetas_all) <- c(params, 'm')
d_all$m <- rep(rep(c('20','50'),each=N),NP)
names(d_all) <- c('value', 'parameter','m')


d_all$iterate <- rep(1:N,NUM*NP)
d_all$yintercept <- rep(c(0.1, 0.7, 0.35, 0.2, 0.1, 0.9, 0.3, 0.1, 0.1/0.7, 0.1/0.9), each=NUM*N)

textlarge <- theme(text = element_text(size=18))

iterates_first_fig <- ggplot(subset(d_all,iterate<=500), aes(x=iterate,y=value)) + geom_path() + facet_grid(parameter~m,scales='free',labeller= label_parsed)  + geom_hline(aes(yintercept=yintercept),color='red') + ylab("") + textlarge +theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))


d_all2 <- subset(d_all,iterate>500) # burnin removed
d_all2thin <- subset(d_all2,(iterate%%50)==0)
iterates_thin_fig <- ggplot(d_all2thin, aes(x=iterate,y=value)) + geom_path() + facet_grid(parameter~m,scales='free',labeller= label_parsed)  + geom_hline(aes(yintercept=yintercept),color='red') + ylab("") + textlarge+ textlarge +theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))

den_fig = list(NP)
# density plots
i = 1
for (i in 1:NP)
{
den_fig[[i]] <- ggplot(subset(d_all2,parameter==params[i]), aes(x=value,colour=m,linetype=m)) + geom_line(stat='density',size=1.5) + facet_grid(.~parameter,scales='free',labeller= label_parsed)   +xlab("") + textlarge + scale_colour_discrete(guide=FALSE) +scale_linetype_discrete(guide=FALSE) + facet_grid(.~parameter,scales='free_x',labeller= label_parsed)    
}
#denc2_fig <- ggplot(subset(d_all2,parameter=='c2'), aes(x=value,colour=m,linetype=m)) + 
  #geom_line(stat='density',size=1.5)+ facet_grid(.~parameter,scales='free_x',labeller= label_parsed)   +xlab("")+ textlarge+ scale_colour_discrete(guide=FALSE) +scale_linetype_discrete(guide=FALSE) 


#denc3_fig <- ggplot(subset(d_all2,parameter=='c3'), aes(x=value,colour=m,linetype=m)) + geom_line(stat='density',size=1.52) + facet_grid(.~parameter,scales='free_x',labeller= label_parsed)   +xlab("")+ textlarge




#### acf plots
lagmax <- 30

ind <- sort(unique(d_all2thin$iterate))
cs = c(1,3,7)
c1.acf <- c(
               acf(thetas_all20[ind,cs[1]],lagmax,plot=F)$acf,
               acf(thetas_all50[ind,cs[1]],lagmax,plot=F)$acf)
c2.acf <- c(
               acf(thetas_all20[ind,cs[2]],lagmax,plot=F)$acf,
               acf(thetas_all50[ind,cs[2]],lagmax,plot=F)$acf)
c3.acf <- c(
               acf(thetas_all20[ind,cs[3]],lagmax,plot=F)$acf,
               acf(thetas_all50[ind,cs[3]],lagmax,plot=F)$acf)
acf.df <- data.frame(lag=rep(0:lagmax,3*length(c)),
                     acf=c(c1.acf,c2.acf,c3.acf),
                     parameter=rep(params[cs],each=NUM*(lagmax+1  )))
acf.df$m <- rep(c(rep('20',lagmax+1),rep('50',lagmax+1)),3)
head(acf.df)
acf_fig <- ggplot(data=acf.df, aes(x=lag, y=acf)) + 
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0))+ facet_grid(parameter~m,labeller= label_parsed)+ textlarge



pdf('Cpairs_fig_par.pdf')
par(mfrow=c(1,NUM))
#pairs plots
prs= c(5,6)

t20 <- subset(thetas_all,m=='20')
colnames(t20) <- paramsfn
show(ggpairs(t20[seq(100, 100000,10),prs],title='m=20',diag=list(continuous='density', combo = "box"), axisLabels="show",params=list(cex=0.8)))
 
t50 <- subset(thetas_all,m=='50')
colnames(t50) <- paramsfn
show(ggpairs(t50[seq(100, 100000,10),prs],title='m=50',diag=list(continuous='density', combo = "box"), axisLabels="show",params=list(cex=0.8)))


dev.off();


ddply(d_all2, parameter~m, summarise, 
      mean = round(mean(value),4), 
      sd = round(sd(value),4) )


pdf('Citerates_first_fig_par.pdf')
show(iterates_first_fig)
dev.off()

pdf('Citerates_thin_fig_par.pdf')
show(iterates_thin_fig)
dev.off()

for (i in 1:NP)
{
  pdf(paste('Cden', paramsfn[i], 'fig_par.pdf',sep=''))
  show(den_fig[[i]])
  dev.off()
}

 

pdf('Cacf_fig_par.pdf')
show(acf_fig)
dev.off()


