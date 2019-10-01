library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(lme4)
source('summarySE2.R')

##data files:
uneq = 'uneq.csv'
eqfo = 'eqfo.csv'

#set colors for plotting
red = brewer.pal(9,"Reds")
blue = brewer.pal(9,"Blues")
purple = brewer.pal(9,"Purples")
green = brewer.pal(9,"Greens")

###Accuracy (F1Ai)


#Run Behavioral Model 1
BM1 =lmer(RT~certN + certD +pref + half +  RL + #
              (pref*certN) + (half * certN)  + (half*pref) + 
              (pref*certD) + (half * certD) +
              # (acc*certN) + (acc*certD) +
              (1+certN+certD+half+pref+RL|subj) , data=uneq) # + # + certN+half+acc+pref+RL
summary(BM1)

###summarise parameters for plotting
sumsP0 = as.data.frame(coef(summary(BM1)))
sumsP0$Parameter = rownames(sumsP0)
sumsP0$se = sumsP0$`Std. Error`
sumsP0 = sumsP0[2:nrow(sumsP0),]
sumsP0$Parameter = c('Cert Sum', 'Cert Diff', 'Pref', 'Cue Remap', 'Right-side bias', 
                     'Cert Sum:Pref', 'Cert Sum:Cue Remap', 'Pref:Cue Remap',
                     'Cert Diff:Pref', 'Cert Diff:Cue Remap')
sumsP0$Parameter = ordered(sumsP0$Parameter, 
                           levels = rev(sumsP0$Parameter))
sumsP0$col = ifelse(sumsP0$Pr<0.001, red[7], ifelse(sumsP0$Pr<0.01, red[5], ifelse(sumsP0$Pr<0.05, red[3], '#FFFFFF')))

##plot in a new window
dev.new()
ggplot(data = sumsP0, aes(x=Parameter, y = Estimate, color=col, fill=col)) +
  coord_flip() +                         # set axis
  ggtitle('Predictors of RT \n(Unequal Choices)') +
  labs(x="", y = "") + 
  theme_classic(base_size = 28)+                            #set theme
  theme(legend.position = 'none'
        ,plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(colour = 'black', size = 30),
        axis.text.y = element_text(colour = 'black')) +
  geom_bar(stat="identity" , position = position_stack(reverse=T),      #[7,5,3,,,2,1]
           alpha = .35, size=1.1, )+        #bar   fill = 'white'
  geom_point( size = 3) +
  geom_errorbar(aes(ymin=Estimate-1*se, ymax=Estimate+1*se), width=.15,    #error
                 size = 1.05) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black', size = 1) +
  scale_color_manual(values=c(red[6], red[3], 'grey50')) + 
  scale_fill_manual(values=c(red[6], red[3], 'white')) 




#load equal and forced conditions
eqfo = read.csv('') ##eqfo.csv path

#run model 2:
library(lmerTest)
BM2 = lmer(RT~certN + choice + accF + halfF +  RL +  #
       (choice*accF*certN) + (choice*halfF * certN)  + (choice*halfF*accF) +
       (1+certN+choice+halfF+accF+RL|subj), data=eqfo) # + (0+half|subj) 

summary(BM2)

#create a summary table for plotting
sumsP1 = as.data.frame(coef(summary(BM2)))
sumsP1$Parameter = rownames(sumsP1)
sumsP1$se = sumsP1$`Std. Error`
sumsP1 = sumsP1[c(2:12, 14),]
sumsP1$Parameter = c('Cert', 'Choice', 'Pref', 'Cue Remap', 'Right-side bias', 
                     'Choice:Pref', 'Cert:Choice', 'Cert:Pref',
                     'Choice:Cue Remap','Cert:Cue Remap','Pref:Cue Remap',
                     'Cert:Choice:Cue Remap')
sumsP1$Parameter = ordered(sumsP1$Parameter, 
                           levels = rev(sumsP1$Parameter))
sumsP1$col = ifelse(sumsP1$Pr<0.001, red[7], ifelse(sumsP1$Pr<0.01, red[5], ifelse(sumsP1$Pr<0.05, red[3], '#FFFFFF')))
#plot model 2 parameters

dev.new()

ggplot(data = sumsP1, aes(x=Parameter, y = Estimate, col=col, fill=col)) +
  coord_flip() +    
  ggtitle('Predictors of RT \n(Equal+Forced)') +
  labs(x="", y = "") + # set axis
  theme_classic(base_size = 28)+                            #set theme
  theme(legend.position = 'none',
        axis.text.x = element_text(colour = 'black', size = 30, angle = 45, hjust = 1),
        axis.text.y = element_text(colour = 'black')
  ) +
  geom_bar(stat="identity" , position = position_stack(reverse=T),      #[7,5,3,,,2,1]
           alpha = .35,  size=1.1)+        #bar   color = 'grey50', 
  geom_point( size = 2.5) +
  geom_errorbar(aes(ymin=Estimate-1*se, ymax=Estimate+1*se), width=.15,    #error
                 size = 1.05) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'black', size = 1) +
  scale_fill_manual(values=c(red[6], 'black')) +
  scale_color_manual(values=c(red[6], red[3], 'grey50')) +
  scale_fill_manual(values=c(red[6], red[3], 'white')) 


