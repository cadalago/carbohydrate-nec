## calculate summary statistics for carbohydrate study
source('pig_utils.R')
library(nlme)
library(multcomp)
library(survival)

## get data from pigs used in this study; specify to include newborn pigs or not based on analysis being done
pigs <- getMJNPigs(incl.new = F)

## descriptive summaries
with(pigs, aggregate(weight.b ~ group + clin.nec, FUN=mean))
with(pigs, aggregate(weight.b ~ group + clin.nec, FUN=se))
# repeat for other morphometrics (i.e. weight.b, weight.gain, vh, cd)


# incidence
addmargins(xtabs(~ litter + group, data=pigs))
with(pigs, addmargins(table(litter, clin.nec, useNA='ifany')))
addmargins(table(pigs$group, pigs$clin.nec, useNA='ifany'))
addmargins(prop.table(xtabs(~ group + clin.nec, data=pigs, drop.unused.levels=T), 1), 2)
addmargins(prop.table(xtabs(~ litter + clin.nec, data=pigs, drop.unused.levels=T), 1), 2)
xtabs(~ litter + interaction(group, clin.nec), data=pigs, sparse=T)[, 4:6] / xtabs(~ litter + group, exclude='NEW', sparse=T, data=pigs)  # NEC rate x litter:group

# severity
with(pigs, aggregate(clin.tot ~ group, FUN=mean))
with(pigs, aggregate(clin.tot ~ group, FUN=range))
with(pigs, aggregate(clin.tot ~ group + clin.nec, FUN=mean))
with(pigs, aggregate(clin.tot ~ group + clin.nec, FUN=range))
with(pigs, aggregate(clin.tot ~ litter, FUN=mean))
with(pigs, aggregate(clin.tot ~ litter, FUN=range))
with(pigs, aggregate(clin.tot ~ litter + group, FUN=mean))
with(pigs, aggregate(clin.tot ~ litter + group, FUN=range))
with(pigs, aggregate(clin.tot ~ litter + clin.nec, FUN=mean))
with(pigs, aggregate(clin.tot ~ litter + clin.nec, FUN=range))

## simple stats
with(pigs, fisher.test(clin.nec.c, group))
with(pigs, summary(aov(clin.tot ~ group + weight.b + farm)))
with(pigs[which(pigs$clin.nec == 'Healthy'),], summary(aov(clin.tot ~ group)))
with(pigs[which(pigs$clin.nec == 'NEC'),], summary(aov(clin.tot ~ group)))

## more stats
# incidence
logreg <- glm(clin.nec.c ~ group + weight.b + farm, family=binomial, data=pigs)
summary(logreg)
summary(glht(logreg, linfct = mcp(group='Tukey')))  #  Tukey (all-vs-all) or Dunnett (all-vs-CTRL)

# all pigs severity
sevreg <- glm(hist.co ~ group + weight.b + farm, family=gaussian, data=pigs)
summary(sevreg)
summary(glht(sevreg, linfct = mcp(group='Tukey')))  #  Tukey (all-vs-all) or Dunnett (all-vs-CTRL)
# NEC pigs severity
sevreg <- glm(clin.st ~ group + weight.b + farm, family=gaussian, data=pigs[which(pigs$clin.nec == 'NEC'),])
summary(sevreg)
summary(glht(sevreg, linfct = mcp(group='Tukey')))  #  Tukey (all-vs-all) or Dunnett (all-vs-CTRL)
# Healthy pigs severity
sevreg <- glm(clin.je ~ group + weight.b + farm, family=gaussian, data=pigs[which(pigs$clin.nec == 'Healthy'),])
summary(sevreg)
summary(glht(sevreg, linfct = mcp(group='Tukey')))  #  Tukey (all-vs-all) or Dunnett (all-vs-CTRL)

rm(logreg, sevreg)
##-------------------------------------------------------------end descriptive summaries


## plots

colors=c('#4daf4a', '#377eb8', '#e41a1c')    # 'palegreen2', 'goldenrod1', 'red3')
hn_colors <- c('turquoise3', 'orangered1')

# survival curves
# count any NEC pig sacrificed before 9am (beginning of collection day) (aka 117 = 3 hrs less than 120) as an event
surv_model <- with(pigs, Surv(hrs.after.feed, as.numeric(hrs.after.feed<119 & clin.nec=='NEC')))

sf  <-  survfit(surv_model ~ group, conf.type="log", 
              conf.int=0.95, type="kaplan-meier", error="greenwood", data=pigs)

op <- par(font=2, cex.axis=1.5, cex.lab=1.8, mar = c(5, 7, 4, 2) + 0.1)
plot(sf, lty=1, lwd=4, font=2, font.lab=2, mark.time=TRUE, xlab='Hours of enteral feeding', frame.plot=F, xaxt='n', yaxt='n', col=colors)
  rect(0, 0, 124, 1, col = 'white', lty=0)
  axis(1, lwd=2, font=2, font.lab=2, at=c(0, 140), lwd.ticks=0, pos=c(0, 0))
  axis(1, at=seq(0, 130, by=24), lwd=0, lwd.ticks=1, font=2, pos=c(0, 0))
  axis(2, lwd=2, las=2, font=2, font.lab=2, at=seq(0, 1, 0.2), labels=paste(seq(0, 100, 20), "%", sep=""))
  lines(sf, lty=1, lwd=4, mark.time=F, col=colors)
  title(ylab='Non-NEC Survival', font.lab=2, cex.lab=1.8, line=4.5)
  legend(text.width=14, x=22, y=0.35, legend=c('LAC', 'MIX', 'CSS'), col=colors, lty=1, lwd=4, bty="y", bg='white', )
par(op)

# Log-rank test for differences in survival curves
survdiff(surv_model~group, data=pigs, rho=0)  ## overall
survdiff(surv_model~group, data=pigs, subset=group %in% c('LAC', 'CSS'))  ## pairwise
survdiff(surv_model~group, data=pigs, subset=group%in%c('LAC', 'MIX'))  ## pairwise
survdiff(surv_model~group, data=pigs, subset=group%in%c('CSS', 'MIX'))  ## pairwise

rm(surv_model, sf, op)

## plot severity scores across GIT segments
library(reshape2)
library(ggplot2)
library(Rmisc)
theme_nogrid <- function (base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(
      panel.grid = element_blank(),
      plot.background = element_blank()
    )   
}
theme_set(theme_nogrid())

idvars <- c('pigID', 'litter', 'farm', 'group', 'sex', 'hist.nec')
givars <- c('hist.pj', 'hist.di', 'hist.co') # or 'clin.st','clin.je','clin.il','clin.co,'clin.tot'
mpigs <- melt(pigs[, c(idvars, givars)], id.vars=idvars, variable.name='site')   # NEC or Healthy
tmp <- summarySE(mpigs, measurevar = 'value', groupvars = c('group', 'site'), na.rm=TRUE) # summarySE function below
with(tmp, ggplot(data=tmp, aes(x=site, fill=group, y=value)) + 
		 	 geom_errorbar(aes(ymin = value - se, ymax = value + se), width = 0.2, size = 1.5, position = position_dodge(0.9)) + 
		   geom_bar(stat='identity', position=position_dodge()) + 
		 	 
       # ylim(0, 6) + 
       scale_fill_manual(values=colors) + 
       labs(y='mean severity score'))

# mpigs <- melt(pigs[, c('pigID', 'group', givars)])
# with(mpigs, stripchart(value ~ variable, method='jitter', pch=16, col=rainbow(length(levels(variable))), vertical=T))

rm(idvars, givars, mpigs, tmp)

### plot NEC incidence
library(Rmisc)
library(scales)

t<-with(pigs,table(group,clin.nec))
pigs$clin.nec.logical <- pigs$clin.nec == 'NEC'
clin_nec_incidence_summary<-summarySE(pigs,measurevar='clin.nec.logical',groupvars='group')
ggplot(clin_nec_incidence_summary,aes(x=group,y=clin.nec.logical,fill=group)) + 
  scale_fill_manual(values=colors) + geom_bar(stat='identity',width=0.8) +
  annotate('text',x=1:3,y=0.05,label=paste(t[,2],clin_nec_incidence_summary$N,sep='/'),fontface='bold', color='white', cex=6) +
  labs(list(title='NEC Incidence',x='',y='')) + scale_y_continuous(element_blank()) +
  #   theme(legend.position='none',axis.title.x=element_text(vjust=-0.5),text=element_text(face='bold',size=18),axis.text=element_text(color='black')) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour='black'),axis.title.x=element_text(vjust=-0.5),axis.title.y=element_text(vjust=1.5),text=element_text(face='bold',size=18),axis.text=element_text(color='black'),legend.title=element_blank(),legend.position='none')
rm(t, clin_nec_incidence_summary)


### plot and test DI cytokines
pigs <- getMJNPigs() # need all pigs including NEW
exp = read.csv('../PCR/DI_samples_204-215/DI_cytokines_melt.csv', strip.white=T)
exp <- merge(pigs[, c('pigID', 'group', 'clin.nec')], exp)

brightness <- function(rgbcol, v) {
  conv <- as.list(as.data.frame(t(rgb2hsv(col2rgb(rgbcol)))))
  conv[[3]] <- v
  do.call(hsv, conv)
}  
new_colors <- c('#ffc0cb', brightness(colors[1], seq(0.7, 0.4, length = 2)), brightness(colors[2], seq(0.7, 0.4, length = 2)), brightness(colors[3], seq(0.9, 0.7, length = 2)))

# x=gene, dodge by groups and/or NEC
exp$group <- paste(exp$group, exp$clin.nec, sep='.')
exp$group <- factor(exp$group, levels = c('NEW.NA', 'LAC.Healthy', 'LAC.NEC', 'MIX.Healthy', 'MIX.NEC', 'CSS.Healthy', 'CSS.NEC'), ordered = T)
sum_c <- summarySE(exp[exp$clin.nec %in% c('NEC', 'NA'), ], measurevar = 'rel.exp', groupvars = c('group', 'clin.nec', 'gene'), na.rm=T)  # find mean and sd of given cytokine
ggplot(sum_c, aes(x=gene, y=rel.exp, fill = group)) +
  scale_fill_manual(values = colors) + scale_colour_manual(values = colors) + # use new_colors if dodge by NEC
  geom_errorbar(aes(ymin = rel.exp - se, ymax = rel.exp + se), width = 0.2, position = position_dodge(width = 0.9)) +
  geom_bar(stat = 'identity', aes(colour = group), position = position_dodge(width = 0.9)) +
  #scale_y_log10(breaks=c(c(0.1, 0.3, 1, 3), c(10, 30, 100))) +
  #scale_y_sqrt(trans_breaks('sqrt', function(x) x^2)(c(1, 100)), limits = c(0, 100)) +
  #scale_y_log10(trans_breaks('log10', function(x) 10^x)(c(1, 200000)), labels = trans_format('log10', math_format(10^.x))) +
  labs(list(title='DI Cytokines Expression (L204-215)', x=element_blank(), y='Relative Expression')) +
  theme(panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour='black'),axis.title.x=element_text(vjust=-0.5),axis.title.y=element_text(vjust=1.5),text=element_text(face='bold',size=18),axis.text=element_text(color='black'),legend.title=element_blank(),legend.position = 'right') +
  guides(colour=guide_legend(nrow=1,byrow=T))




## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
											conf.interval=.95, .drop=TRUE) {
	library(plyr)
	
	# New version of length which can handle NA's: if na.rm==T, don't count them
	length2 <- function (x, na.rm=FALSE) {
		if (na.rm) sum(!is.na(x))
		else       length(x)
	}
	
	# This does the summary. For each group's data frame, return a vector with
	# N, mean, and sd
	datac <- ddply(data, groupvars, .drop=.drop,
								 .fun = function(xx, col) {
								 	c(N    = length2(xx[[col]], na.rm=na.rm),
								 		mean = mean   (xx[[col]], na.rm=na.rm),
								 		sd   = sd     (xx[[col]], na.rm=na.rm)
								 	)
								 },
								 measurevar
	)
	
	# Rename the "mean" column    
	datac <- rename(datac, c("mean" = measurevar))
	
	datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
	
	# Confidence interval multiplier for standard error
	# Calculate t-statistic for confidence interval: 
	# e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
	ciMult <- qt(conf.interval/2 + .5, datac$N-1)
	datac$ci <- datac$se * ciMult
	
	return(datac)
}


summary(glht(glm(weight.b ~ group, family=gaussian, data=pigs), linfct = mcp(group='Tukey')))
summary(glht(glm(weight.d ~ clin.nec, family=gaussian, data=pigs), linfct = mcp(clin.nec='Tukey')))
summary(glht(glm(weight.gain ~ clin.nec, family=gaussian, data=pigs), linfct = mcp(clin.nec='Tukey')))
summary(glht(glm(pj.vh ~ clin.nec, family=gaussian, data=pigs), linfct = mcp(clin.nec='Tukey')))
summary(glht(glm(di.vh ~ clin.nec, family=gaussian, data=pigs), linfct = mcp(clin.nec='Tukey')))
summary(glht(glm(pj.cd ~ clin.nec, family=gaussian, data=pigs), linfct = mcp(clin.nec='Tukey')))
summary(glht(glm(di.cd ~ clin.nec, family=gaussian, data=pigs), linfct = mcp(clin.nec='Tukey')))
summary(glht(glm(co.cd ~ clin.nec, family=gaussian, data=pigs), linfct = mcp(clin.nec='Tukey')))

pigs.h <- pigs[which(pigs[, 'clin.nec'] == 'Healthy'), ]
l.228 <- which(pigs.h[, 'litter'] == 228)
summary(glht(glm(weight.b ~ group, family=gaussian, data=pigs.h[-l.228, ]), linfct = mcp(group='Tukey')))
summary(glht(glm(pj.vh ~ group, family=gaussian, data=pigs.h[-l.228, ]), linfct = mcp(group='Tukey')))
summary(glht(glm(di.vh ~ group, family=gaussian, data=pigs.h[-l.228, ]), linfct = mcp(group='Tukey')))
summary(glht(glm(pj.cd ~ group, family=gaussian, data=pigs.h[-l.228, ]), linfct = mcp(group='Tukey')))
summary(glht(glm(di.cd ~ group, family=gaussian, data=pigs.h[-l.228, ]), linfct = mcp(group='Tukey')))
summary(glht(glm(co.cd ~ group, family=gaussian, data=pigs.h[-l.228, ]), linfct = mcp(group='Tukey')))

with(pigs[which(pigs$group == 'CSS'), ], t.test(weight.b ~ clin.nec))
with(pigs[which(pigs$group == 'CSS'), ], t.test(weight.d ~ clin.nec))
with(pigs[which(pigs$group == 'CSS'), ], t.test(weight.gain ~ clin.nec))
with(pigs[which(pigs$group == 'CSS'), ], t.test(pj.vh ~ clin.nec))
with(pigs[which(pigs$group == 'CSS'), ], t.test(pj.cd ~ clin.nec))
with(pigs[which(pigs$group == 'CSS'), ], t.test(di.vh ~ clin.nec))
with(pigs[which(pigs$group == 'CSS'), ], t.test(di.cd ~ clin.nec))
with(pigs[which(pigs$group == 'CSS'), ], t.test(co.cd ~ clin.nec))
