## includes utility functions to be used in analysis of piglet NEC carbohydrate study

getMJNPigs <- function(incl.new = TRUE) {
		# read in data
		pigs <- read.csv('../data/morphometry_incidence_severity/pigs_anthropometrics_morphometry_severity.csv', strip.white=T)
		
		# parse 'pigID' into separate columns for litter, letter
		pigs <- cbind(factor(unlist(lapply(pigs[, 'pigID'], function(x){substr(x, 4, 4)}))), pigs)
		pigs <- cbind(factor(unlist(lapply(pigs[, 'pigID'], function(x){as.numeric(substr(x, 1, 3))}))), pigs)
		names(pigs)[1:2] <- c('litter', 'letter')
		
		# add in calculated columns to pigs
		pigs$clin.tot <- pigs$clin.st+pigs$clin.je+pigs$clin.il+pigs$clin.co
		pigs$hist.tot <- pigs$hist.pj+pigs$hist.di+pigs$hist.co
		pigs$clin.max <- apply(pigs[, c('clin.st', 'clin.je', 'clin.il', 'clin.co')], 1, FUN=max)
		pigs$hist.max <- apply(pigs[, c('hist.pj', 'hist.di', 'hist.co')], 1, FUN=max)
		pigs$clin.nec.c <- apply(pigs[, c('clin.st', 'clin.je', 'clin.il', 'clin.co')], 1, function(x){max(x) > 2})
		pigs$clin.nec.c <- factor(pigs$clin.nec.c, labels=c('Healthy', 'NEC'))
		
		pigs$dob <- as.Date(pigs$dob)
		pigs$dod <- as.Date(pigs$dod)
		
		# convert time of death (tod) to hours in decimal form
		pigs$tod <- sapply(strsplit(as.character(pigs[, 'tod']), ':'), function(x){
		  x <- as.numeric(x)
		  x[1]+x[2]/60
		})
		
		pigs$hrs.after.feed <- as.numeric((pigs$dod-pigs$dob-2)*24 + pigs$tod - 12) # 24 hrs for each day of feeding (first 2 days of life are TPN only) plus the time of death on final day minus 12 hrs b/c feeding begins at 12:00
		
		# weight.gain expressed as:  g / (kg * day)  from birthweight to deathweight
		pigs$weight.gain <- (pigs$weight.d - pigs$weight.b) / ((pigs$weight.b / 1000) * (pigs$hrs.after.feed / 24))
		
		# according to preference passed in function call, include or remove all newborn
		if(incl.new) {
			pigs$group <- factor(pigs$group, levels=c('NEW', 'LAC', 'MIX', 'CSS'), ordered=T)  # force use of this order L -> R	
		}
		else {
			pigs <- pigs[which(pigs[, 'group'] != 'NEW'), ]  
			pigs$group <- droplevels(pigs$group)
			pigs$group <- factor(pigs$group, levels=c('LAC', 'MIX', 'CSS'), ordered=T)  # force use of this order L -> R
		}
		
		# remove pigs from litters 234 and 248-249
		pigs <- pigs[!c(pigs$litter %in% c(234, 248, 249)), ]
		
		# remove pigs with perforated stomachs
		pigs <- pigs[!c(pigs$pigID %in% c('215E','215I','215L')), ]
		
		# remove pigs which received lactose or mix diet from litters 212, 213 (bad batch of formula)
		pigs <- pigs[!c(pigs$litter %in% c(212, 213) & pigs$group %in% c('LAC', 'MIX')), ]
		pigs$litter <- droplevels(pigs$litter)
		
		return(pigs)
}

se <- function(x) {
  
  sd(x) / sqrt(length(x))
}

