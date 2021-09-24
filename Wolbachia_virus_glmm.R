
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(MCMCglmm)
library(boot)
library(reshape2)
library(ggplot2)



priorb = list(R = list(V = 1, fix = 1),
              G = list(G1 = list(V = diag(2), nu = 0.002),
              G2 = list(V = 1, nu = 0.002)))

#model4 <- MCMCglmm(fixed=inf~wolbachia,
#                   random=~idh(wolbachia):virus+fly,
#                   family="categorical",
#                   data=inf_status2, 
#                   nitt=106000,thin=100,burnin=6000,
#                   pr=T,
#                   prior=priorb)

#save(model4, file = "mixed_model_fit.rda")
load("mixed_model_fit.rda")


########################extract parameter estimates of interest from the model#############
#posteriors of flies +/- wolbachia
#add intercept
wolb_posterior=model4$Sol[,3:24]+model4$Sol[,1]
#add wolb fixed effect
wolb_posterior[,12:22]=wolb_posterior[,12:22]+model4$Sol[,2]

#get posterior probabilities for an effect of Wolbachia on each virus
effect=wolb_posterior[,1:11]-wolb_posterior[,12:22]
HPDinterval(effect)
Pmcmc=function(x){length(which(x<0))/length(x)}
Pval=apply(effect,2,Pmcmc)

#estimate Risk Ratio 
wolb_posterior_proportion=inv.logit(wolb_posterior)
risk_ratio=wolb_posterior_proportion[,1:11]/wolb_posterior_proportion[,12:22]
rr=colMeans(risk_ratio)
x=HPDinterval(risk_ratio)
virus=sub("wolbachiaw-.virus.","",row.names(x))
rr2=cbind(virus,rr,x)
colnames(rr2)=c("virus","rr","lowCI","upCI")

#estimate prevalence
prevalence=colMeans(wolb_posterior_proportion)
cis=HPDinterval(wolb_posterior_proportion)
prevalence2=cbind(prevalence,cis)

##########################plot risk ratios#################################################


rr3=data.frame(rr2[order(rr2[,2]),])
rr3$rr=as.numeric(as.character(rr3$rr))
rr3$lowCI=as.numeric(as.character(rr3$lowCI))
rr3$upCI=as.numeric(as.character(rr3$upCI))
rr3$virus=as.factor(rr3$virus)
rr3$virus <- factor(rr3$virus, levels = rr2[order(rr2[,2]),1])

p=ggplot(rr3, aes(x=virus, y=rr))+
  geom_point(stat='identity',fill="grey50")+
  #  geom_bar(stat='identity',width=0.8,fill="grey50")+
  geom_errorbar(data=rr3, aes(ymin=lowCI, ymax=upCI), width=.4)+ 
  scale_x_discrete(labels=c("c.hill" = 'Craigie\'s Hill',"chaq" = "Galbut (Chaq)","dav" = "DAV","galbut"="Galbut",
                            "motts_mill"="Motts Mill",'kallithea'='Kallithea','la.jolla'='La Jolla','mottsmill'='Motts Mill',
                            'new.0370'='Vera','vera'='Vera',
                            'nora'='Nora','sigma'='Sigma','thika'='Thika'))+
  geom_hline(yintercept=1,colour="red")+
  
  
  ylab(expression('Risk Ratio'))+
  xlab('Virus')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,  color='black', size=10,hjust=1),
        axis.text.y = element_text(size=10),
        axis.title.y = element_text(colour = "black",  size=10),
        axis.title.x = element_text(colour = "black",  size=10),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))




pdf(file="risk_ratio_glmm_estimates.pdf",height=4,width=3.8)
p
dev.off()




##########################plot prevalence#################################################


virus=substr(rownames(prevalence2), 19, nchar(rownames(prevalence2)))
wolbachia=rep("Wolbachia",length(virus))
wolbachia[grep("w-",rownames(prevalence2))]="No Wolbachia"
prevalence3=cbind(prevalence2,virus,wolbachia)
prevalence4=data.frame(prevalence3)
prevalence4[,1]=as.numeric(as.character(prevalence4[,1]))
prevalence4[,2]=as.numeric(as.character(prevalence4[,2]))
prevalence4[,3]=as.numeric(as.character(prevalence4[,3]))
prevalence4=prevalence4[c(order(rr2[,2]),order(rr2[,2])+length(rr2[,2])),]
prevalence4$virus=factor(prevalence4$virus,levels=prevalence4$virus[1:(nrow(prevalence4)/2)])




p2=ggplot(data=prevalence4, aes(x=virus, y=prevalence, fill=wolbachia)) +
  
  geom_bar(stat="identity", position=position_dodge(), width=0.75)+
  scale_x_discrete(labels=c("c.hill" = 'Craigie\'s Hill',"chaq" = "Galbut (Chaq)","dav" = "DAV","galbut"="Galbut",
                            "grom"="Motts Mill",'kallithea'='Kallithea','la.jolla'='La Jolla','motts_mill'='Motts Mill',
                            'vera'='Vera','partitiv.'='Vera',
                            'nora'='Nora','sigma'='Sigma','thika'='Thika'))+
  scale_fill_manual(values=c('#2C5F2D','#97BC62FF'))+
  ylab(expression('Prevalence'))+
  xlab('Virus')+
  geom_text(x=8, y=0.8, label="No Wolbachia",colour = "#2C5F2D",hjust=0)+
  geom_text(x=8, y=0.73, label="Wolbachia",colour = "#97BC62FF",hjust=0)+
  annotate("text", x = 1:length(Pval), y = prevalence4$prevalence[1:(nrow(prevalence4)/2)]+0.05, 
           label = paste("p=",format(Pval[order(rr2[,2])]),sep=""),vjust=0,size=2.7)+
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,  color='black', size=10,hjust=1),
        axis.text.y = element_text(size=10),
        axis.title.y = element_text(colour = "black",  size=10),
        axis.title.x = element_text(colour = "black",  size=10),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))


p2

pdf(file="prevalence_glmm_estimates.pdf",height=4,width=5.8)
p2
dev.off()

raw_prev=read.csv("raw_virus_prevalences.csv")[,1:4]
prevalence5=merge(prevalence4,raw_prev)
write.csv(prevalence5,file="model_estimates_prevalence.csv")

write.csv(rr3,file="model_estimates_risk_ratios.csv")



