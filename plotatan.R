library(ggplot2); library(gridExtra);library(plyr);library(GGally)


#if (Sys.info()["sysname"] == "Darwin") setwd("~/Dropbox/GPC/working document/sim/prog")
truth_id <- "Inno.Atan100T30.0"
TC = "TC"

sim_id <- paste("Inno.Atan10000lin100x10T30.0", TC, sep="")
thetas_all10 <- read.csv(paste(sim_id,"/thetas",sim_id,".csv",sep=""),header=FALSE)
sim_id <- paste("Inno.Atan10000lin100x1000T30.0", TC, sep="")
thetas_all100 <- read.csv(paste(sim_id,"/thetas",sim_id,".csv",sep=""),header=FALSE)
sim_id <- paste("Inno.Atan10000lin100x100T30.0", TC, sep="")
thetas_all1000 <- read.csv(paste(sim_id,"/thetas",sim_id,".csv",sep=""),header=FALSE)
thetas_all <- rbind(thetas_all10, thetas_all100,thetas_all1000)
colnames(thetas_all) <- c('alpha', 'beta', 'sigma')
d_all <- stack(thetas_all)
head(d_all)
N <- nrow(thetas_all10)
thetas_all$m <- rep(c('10','100','1000'),each=N)
colnames(thetas_all) <- c('alpha', 'beta', 'sigma', 'm')
d_all$m <- rep(rep(c('10','100','1000'),each=N),3)
names(d_all) <- c('value', 'parameter','m')


d_all$iterate <- rep(1:N,9)
d_all$yintercept <- rep(c(-2.0,0.,0.75), each=3*N)

textlarge <- theme(text = element_text(size=18))

iterates_fig <- ggplot(subset(d_all,iterate<=500), aes(x=iterate,y=value)) + geom_path() + facet_grid(parameter~m,scales='free',labeller= label_parsed)  + geom_hline(aes(yintercept=yintercept),color='red') + ylab("") + textlarge +theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))


d_all2 <- subset(d_all,iterate>500) # burnin removed

iterates_all_fig <- ggplot(subset(d_all2,(iterate%%5)==0), aes(x=iterate,y=value)) + geom_path() + facet_grid(parameter~m,scales='free',labeller= label_parsed)  + geom_hline(aes(yintercept=yintercept),color='red') + ylab("") + textlarge+ textlarge +theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))


# density plots
denalpha_fig <- ggplot(subset(d_all2,parameter=='alpha'), aes(x=value,colour=m,linetype=m)) + geom_line(stat='density',size=1.5) + facet_grid(.~parameter,scales='free',labeller= label_parsed)   +xlab("") + textlarge + scale_colour_discrete(guide=FALSE) +scale_linetype_discrete(guide=FALSE) 


denbeta_fig <- ggplot(subset(d_all2,parameter=='beta'), aes(x=value,colour=m,linetype=m)) + 
  geom_line(stat='density',size=1.5)+ facet_grid(.~parameter,scales='free_x',labeller= label_parsed)   +xlab("")+ textlarge+ scale_colour_discrete(guide=FALSE) +scale_linetype_discrete(guide=FALSE) 


densigma_fig <- ggplot(subset(d_all2,parameter=='sigma'), aes(x=value,colour=m,linetype=m)) + geom_line(stat='density',size=1.52) + facet_grid(.~parameter,scales='free_x',labeller= label_parsed)   +xlab("")+ textlarge



# runmean <- function(x,Nburnin) 
# { 
#   N <- length(x)
#   c(rep(NA,Nburnin),cumsum(x[seq(Nburnin+1,N)])/seq(1,N-Nburnin))
# }  
# 
# Nburnin <- 100
# 
# d_all$runmean <- c(runmean(thetas_all10[,1],Nburnin),runmean(thetas_all100[,1],Nburnin),runmean(thetas_all1000[,1],Nburnin),runmean(thetas_all10[,2],Nburnin),runmean(thetas_all100[,2],Nburnin),runmean(thetas_all1000[,2],Nburnin),runmean(thetas_all10[,3],Nburnin),runmean(thetas_all100[,3],Nburnin),runmean(thetas_all1000[,3],Nburnin))
# 
# 
# p1 <- ggplot(subset(d_all,iterate>Nburnin), aes(x=iterate,y=runmean)) + geom_path() + facet_grid(parameter~m,scales='free',labeller= label_parsed)  + geom_hline(aes(yintercept=yintercept),color='red') + ylab("") + ggtitle('Running means from data-augmentation algorithm')
# 
# pdf(paste(sim_id,"/",sim_id,".pdf",sep=""),width=5,height=5)
# show(p1)
# dev.off()


#### acf plots
lagmax <- 50

ind <- (1:500)
alpha.acf <- c(acf(thetas_all10[-ind,1],lagmax,plot=F)$acf,
               acf(thetas_all100[-ind,1],lagmax,plot=F)$acf,
               acf(thetas_all1000[-ind,1],lagmax,plot=F)$acf)
beta.acf <- c(acf(thetas_all10[-ind,2],lagmax,plot=F)$acf,
               acf(thetas_all100[-ind,2],lagmax,plot=F)$acf,
               acf(thetas_all1000[-ind,2],lagmax,plot=F)$acf)
sigma.acf <- c(acf(thetas_all10[-ind,3],lagmax,plot=F)$acf,
               acf(thetas_all100[-ind,3],lagmax,plot=F)$acf,
               acf(thetas_all1000[-ind,3],lagmax,plot=F)$acf)
acf.df <- data.frame(lag=rep(0:lagmax,9),
                     acf=c(alpha.acf,beta.acf,sigma.acf),
                     parameter=rep(c('alpha','beta','sigma'),each=3*(lagmax+1  )))
acf.df$m <- rep(c(rep('10',lagmax+1),rep('100',lagmax+1),rep('1000',lagmax+1)),3)
head(acf.df)
acf_fig <- ggplot(data=acf.df, aes(x=lag, y=acf)) + 
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0))+ facet_grid(parameter~m,labeller= label_parsed)+ textlarge




# 
# 
# xtrue <- read.csv(paste("xInno.Atan100T30.0", truth_id, ".csv",sep=""), header=F)
# colnames(xtrue) <- c('t','x_t')
# TT = xtrue$t[length(xtrue$t)]
# xobs <- read.csv(paste("xobs", truth_id, ".csv",sep=""), header=F)
# colnames(xobs) <- c('t','x_t')
# 
# p4 <- ggplot(xtrue, aes(x=t,y=x_t)) + geom_line() + geom_point(aes(x=t,y=x_t), data=xobs, color='red', size=1, shape=21,fill='white')
# show(p4)
# 
# pdf(paste(truth_id,"diffpath.pdf",sep=""),width=6,height=2.5)
# show(p4)
# dev.off()


# pairs plots
# t10 <- subset(thetas_all,m=='10')
# ggpairs(t10[, 1:3],title='m=10',diag=list(continuous='density', combo = "box"), axisLabels="show",params=list(cex=0.8))
# 
# t100 <- subset(thetas_all,m=='100')
# ggpairs(t100[, 1:3],title='m=100',diag=list(continuous='density', combo = "box"), axisLabels="show")
# 
# t1000 <- subset(thetas_all,m=='1000')
# ggpairs(t1000[, 1:3],title='m=1000',diag=list(continuous='density', combo = "box"), axisLabels="show")



# numerical summary of the results
#data.frame(QV_true=0.75^2,QV_cont=sum(diff(xtrue$x_t)^2)/TT,QV_obs=sum(diff(xobs$x_t)^2)/TT)


ddply(d_all2, parameter~m, summarise, 
      mean = round(mean(value),4), 
      sd = round(sd(value),4) )


pdf('iterates_fig.pdf')
iterates_fig
dev.off()

pdf('iterates_all_fig.pdf')
iterates_all_fig
dev.off()

pdf('denalpha_fig.pdf')
denalpha_fig
dev.off()

pdf('denbeta_fig.pdf')
denbeta_fig
dev.off()

pdf('densigma_fig.pdf')
densigma_fig
dev.off()

pdf('acf_fig.pdf')
acf_fig
dev.off()


