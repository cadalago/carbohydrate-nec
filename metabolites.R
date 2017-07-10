# preprocessing ====================================================================
source('pig_utils.R')
path <- '../data/metabolites/'

#######################  check cadaverine, 1,6-anydroglucose....something doesn't seem right!! ########

# load data
c.cmpds <- read.csv(paste(path, 'cecal/cecal_compounds.csv', sep = ''), strip.white = T)
cecal <- read.csv(paste(path, 'cecal/cecal.csv', sep = ''), strip.white = T)  # "COMP ID" uniquely identifies compounds
p.cmpds <- read.csv(paste(path, 'plasma/plasma_compounds.csv', sep = ''), strip.white = T)
plasma <- read.csv(paste(path, 'plasma/plasma.csv', sep = ''), strip.white = T)
pigs <- getMJNPigs(incl.new = T)
cecal <- merge(pigs, cecal)
plasma <- merge(pigs, plasma)

# get column names for compounds
c.compIDs <- names(cecal)[!(names(cecal) %in% c(names(pigs), 'mass.dry.g', 'bradford.protein'))]
p.compIDs <- names(plasma)[!(names(plasma) %in% names(pigs))]

# remove compounds that have no values (i.e. were only found in 'NEW' pigs)
#  cecal
c.compIDs.toFilter <- names(which(colSums(cecal[, c.compIDs], na.rm = T) == 0))
cecal <- cecal[!names(cecal) %in% c.compIDs.toFilter]
c.compIDs <- setdiff(c.compIDs, c.compIDs.toFilter)
#  plasma
p.compIDs.toFilter <- names(which(colSums(plasma[, p.compIDs], na.rm = T) == 0))
plasma <- plasma[!names(plasma) %in% p.compIDs.toFilter]
p.compIDs <- setdiff(p.compIDs, p.compIDs.toFilter)

# normalize concentration measurements for each cecal sample - divide by dry mass of sample
cecal[, c.compIDs] <- sweep(cecal[, c.compIDs], 1, cecal[, 'mass.dry.g'], '/')

# replace 0s and NA (missing values) w/ 1/2 column minimum
#  cecal
c.colmins <- apply(cecal[, c.compIDs], 2, function(x){ min(x, na.rm = T) })
c.nas <- which(is.na(cecal), arr.ind = T)
c.nas <- c.nas[which(c.nas[, 2] > ncol(cecal) - length(c.compIDs)), ]
cecal[c.nas] <- c.colmins[c.nas[, 2] - (ncol(cecal) - length(c.compIDs))] / 2
#  plasma
p.colmins <- apply(plasma[, p.compIDs], 2, function(x){ min(x, na.rm = T) })
p.nas <- which(is.na(plasma), arr.ind = T)
p.nas <- p.nas[which(p.nas[, 2] > ncol(plasma) - length(p.compIDs)), ]
plasma[p.nas] <- p.colmins[p.nas[, 2] - (ncol(plasma) - length(p.compIDs))] / 2

# # remove 10% of the features, based on highest relative standard deviation
# #  cecal
# c.rsd <- apply(cecal[, c.compIDs], 2, function(x) { sd(x) / mean(x) })
# c.compIDs.toFilter <- names(sort(c.rsd, decreasing = T)[1: round(length(c.compIDs) * 0.1)])
# cecal <- cecal[!names(cecal) %in% c.compIDs.toFilter]
# c.compIDs <- setdiff(c.compIDs, c.compIDs.toFilter)
# #  plasma
# p.rsd <- apply(plasma[, p.compIDs], 2, function(x) { sd(x) / mean(x) })
# p.compIDs.toFilter <- names(sort(p.rsd, decreasing = T)[1: round(length(p.compIDs) * 0.1)])
# plasma <- plasma[!names(plasma) %in% p.compIDs.toFilter]
# p.compIDs <- setdiff(p.compIDs, p.compIDs.toFilter)

# save absolute values to use for calculating fold-change
c.abs <- cecal
p.abs <- plasma

# transform values by generalized log
cecal[, c.compIDs] <- apply(cecal[, c.compIDs], 2, function(x) { log2((x + sqrt(x^2 + 1)) / 2) })
plasma[, p.compIDs] <- apply(plasma[, p.compIDs], 2, function(x) { log2((x + sqrt(x^2 + 1)) / 2) })

# scale values by autoscaling (force each variable to have mean = 0 and sd = 1)
cecal[, c.compIDs] <- apply(cecal[, c.compIDs], 2, function(x) { (x - mean(x)) / sd(x) })
plasma[, p.compIDs] <- apply(plasma[, p.compIDs], 2, function(x) { (x - mean(x)) / sd(x) })

# change column names from 'c.#####' or 'p.#####' format to biochemical name
#  cecal
c.compIDs.nums <- unlist(lapply(c.compIDs, function(x) { substr(x, 3, nchar(x)) }))
c.compIDs.names <- c.cmpds[c.cmpds$comp.id %in% c.compIDs.nums, 'biochemical']
colnames(cecal)[which(colnames(cecal) %in% c.compIDs)] <- as.character(c.compIDs.names)
colnames(c.abs)[which(colnames(c.abs) %in% c.compIDs)] <- as.character(c.compIDs.names)
#  plasma
p.compIDs.nums <- unlist(lapply(p.compIDs, function(x) { substr(x, 3, nchar(x)) }))
p.compIDs.names <- p.cmpds[p.cmpds$comp.id %in% p.compIDs.nums, 'biochemical']
colnames(plasma)[which(colnames(plasma) %in% p.compIDs)] <- as.character(p.compIDs.names)
colnames(p.abs)[which(colnames(p.abs) %in% p.compIDs)] <- as.character(p.compIDs.names)



# test for differences and print out results to csv file
library(multcomp)

# test for differences between health and NEC (fed pigs only)
c.fed <- droplevels(cecal[cecal$group %in% c('LAC', 'MIX', 'CSS'), ])
p.fed <- droplevels(plasma[plasma$group %in% c('LAC', 'MIX', 'CSS'), ])

# #   use welch's t-test
# #    cecal
# c.tt <- apply(c.fed[, as.character(c.compIDs.names)], 2, 
#             function(x) tryCatch(t.test(x ~ c.fed$hist.nec)$p.value, error = function(x) NA))
# c.tt <- p.adjust(c.tt[-which(is.na(c.tt))], method = 'fdr')
# median.h <- apply(c.abs[c.abs$hist.nec == 'Healthy', names(c.tt)], 2, median, na.rm = T)
# median.n <- apply(c.abs[c.abs$hist.nec == 'NEC', names(c.tt)], 2, median, na.rm = T)
# out.df <- cbind('median H/N' = median.h / median.n, 'FDR-adj p' = c.tt, 
#                 'comp.id' = c.cmpds[c.cmpds$biochemical %in% names(c.tt), 'comp.id'])
# out.df <- subset(merge(out.df, c.cmpds), select = -c(comp.id, pathway.sortorder, chemical.id, ri))
# write.csv(out.df[order(out.df$'FDR-adj p'), ],
#           '../output/metabolites/allfed_Healthy_vs_NEC_cecal_t-test.csv', row.names = F)
# #    plasma
# p.tt <- apply(p.fed[, as.character(p.compIDs.names)], 2, 
#             function(x) tryCatch(t.test(x ~ p.fed$hist.nec)$p.value, error = function(x) NA))
# p.tt <- p.adjust(p.tt[-which(is.na(p.tt))], method = 'fdr')
# median.h <- apply(p.abs[p.abs$hist.nec == 'Healthy', names(p.tt)], 2, median, na.rm = T)
# median.n <- apply(p.abs[p.abs$hist.nec == 'NEC', names(p.tt)], 2, median, na.rm = T)
# out.df <- cbind('median H/N' = median.h / median.n, 'FDR-adj p' = p.tt, 
#                 'comp.id' = p.cmpds[p.cmpds$biochemical %in% names(p.tt), 'comp.id'])
# out.df <- subset(merge(out.df, p.cmpds), select = -c(comp.id, pathway.sortorder, chemical.id, ri))
# write.csv(out.df[order(out.df$'FDR-adj p'), ],
#           '../output/metabolites/allfed_Healthy_vs_NEC_plasma_t-test.csv', row.names = F)
# 
# #   use linear model w/ diet as covariate 
# #    cecal
# c.reg.tt <- apply(c.fed[, as.character(c.compIDs.names)], 2, 
#             function(x) t.test(residuals(lm(x ~ c.fed$group)) ~ c.fed$hist.nec)$p.value)
# c.reg.tt <- p.adjust(c.reg.tt, method = 'fdr')
# median.h <- apply(c.abs[c.abs$hist.nec == 'Healthy', names(c.reg.tt)], 2, median, na.rm = T)
# median.n <- apply(c.abs[c.abs$hist.nec == 'NEC', names(c.reg.tt)], 2, median, na.rm = T)
# out.df <- cbind('median H/N' = median.h / median.n, 'FDR-adj p' = c.reg.tt, 
#                 'comp.id' = c.cmpds[c.cmpds$biochemical %in% names(c.reg.tt), 'comp.id'])
# out.df <- subset(merge(out.df, c.cmpds), select = -c(comp.id, pathway.sortorder, chemical.id, ri))
# write.csv(out.df[order(out.df$'FDR-adj p'), ],
#           '../output/metabolites/allfed_Healthy_vs_NEC_cecal_reg.csv', row.names = F)
# #    plasma
# p.reg.tt <- apply(p.fed[, as.character(p.compIDs.names)], 2, 
#             function(x) t.test(residuals(lm(x ~ p.fed$group)) ~ p.fed$hist.nec)$p.value)
# p.reg.tt <- p.adjust(p.reg.tt, method = 'fdr')
# median.h <- apply(p.abs[p.abs$hist.nec == 'Healthy', names(p.reg.tt)], 2, median, na.rm = T)
# median.n <- apply(p.abs[p.abs$hist.nec == 'NEC', names(p.reg.tt)], 2, median, na.rm = T)
# out.df <- cbind('median H/N' = median.h / median.n, 'FDR-adj p' = p.reg.tt, 
#                 'comp.id' = p.cmpds[p.cmpds$biochemical %in% names(p.reg.tt), 'comp.id'])
# out.df <- subset(merge(out.df, p.cmpds), select = -c(comp.id, pathway.sortorder, chemical.id, ri))
# write.csv(out.df[order(out.df$'FDR-adj p'), ],
#           '../output/metabolites/allfed_Healthy_vs_NEC_plasma_reg.csv', row.names = F)



#   use two-way ANOVA (type II, extract p-value for hist.nec effect after group main effect, as well
#                       as p-value for interaction)
#    cecal
c.anova <- as.data.frame(apply(c.fed[, as.character(c.compIDs.names)], 2,
                               function(x) anova(lm(x ~ c.fed$group * c.fed$hist.nec))$'Pr(>F)'[2:3]))
c.anova[1, ] <- p.adjust(c.anova[1, ], method = 'fdr')
c.anova[2, ] <- p.adjust(c.anova[2, ], method = 'fdr')
median.h <- apply(c.abs[c.abs$hist.nec == 'Healthy', names(c.anova)], 2, median, na.rm = T)
median.n <- apply(c.abs[c.abs$hist.nec == 'NEC', names(c.anova)], 2, median, na.rm = T)
out.df <- cbind('median H/N' = median.h / median.n, 'NEC, FDR-adj' = as.numeric(c.anova[1, ]),
                'group:NEC, FDR-adj' = as.numeric(c.anova[2, ]),
                'comp.id' = c.cmpds[c.cmpds$biochemical %in% names(c.anova), 'comp.id'])
out.df <- subset(merge(out.df, c.cmpds), select = -c(comp.id, pathway.sortorder, chemical.id, ri))
write.csv(out.df[order(out.df$'NEC, FDR-adj'), ],
          '../output/metabolites/allfed_Healthy_vs_NEC_cecal_anova.csv', row.names = F)
#    plasma
p.anova <- as.data.frame(apply(p.fed[, as.character(p.compIDs.names)], 2,
                               function(x) anova(lm(x ~ p.fed$group * p.fed$hist.nec))$'Pr(>F)'[2:3]))
p.anova[1, ] <- p.adjust(p.anova[1, ], method = 'fdr')
p.anova[2, ] <- p.adjust(p.anova[2, ], method = 'fdr')
median.h <- apply(p.abs[p.abs$hist.nec == 'Healthy', names(p.anova)], 2, median, na.rm = T)
median.n <- apply(p.abs[p.abs$hist.nec == 'NEC', names(p.anova)], 2, median, na.rm = T)
out.df <- cbind('median H/N' = median.h / median.n, 'NEC, FDR-adj' = as.numeric(p.anova[1, ]),
                'group:NEC, FDR-adj' = as.numeric(p.anova[2, ]), 
                'comp.id' = p.cmpds[p.cmpds$biochemical %in% names(p.anova), 'comp.id'])
out.df <- subset(merge(out.df, p.cmpds), select = -c(comp.id, pathway.sortorder, chemical.id, ri))
write.csv(out.df[order(out.df$'NEC, FDR-adj'), ],
          '../output/metabolites/allfed_Healthy_vs_NEC_plasma_anova.csv', row.names = F)

# test for difference between heatlhy LAC and CSS




# plot interesting metabolites
library(ggplot2)
library(gridExtra)
theme_set(theme_classic())
theme_update(plot.title = element_text(hjust = 0.5))
group_colors <- c('#000000', '#4daf4a', '#377eb8', '#e41a1c')    # aka 'black' (previously 'pink'(ffc0cb)), 'palegreen2', 'goldenrod1', 'red3'
group_NEC_colors <- c('#000000', '#4daf4a', '#316f2f', '#377eb8', '#1e4666', '#e41a1c', '#b01414')
diet_colors <- c('#4daf4a', '#377eb8', '#e41a1c')
hn_colors <- c('turquoise3', 'orangered1')

# add factor to allow plotting by group-NEC
cecal <- cbind(cecal, 'group-NEC' = paste(cecal[, 'group'], cecal[, 'hist.nec'], sep = '-'))
levels(cecal$'group-NEC')[levels(cecal$'group-NEC') == 'NEW-NA'] <- 'NEW'
cecal$'group-NEC' <- factor(cecal$'group-NEC', levels = c('NEW', 'LAC-Healthy', 'LAC-NEC', 'MIX-Healthy', 'MIX-NEC', 'CSS-Healthy', 'CSS-NEC'))
plasma <- cbind(plasma, 'group-NEC' = paste(plasma[, 'group'], plasma[, 'hist.nec'], sep = '-'))
levels(plasma$'group-NEC')[levels(plasma$'group-NEC') == 'NEW-NA'] <- 'NEW'
plasma$'group-NEC' <- factor(plasma$'group-NEC', levels = c('NEW', 'LAC-Healthy', 'LAC-NEC', 'MIX-Healthy', 'MIX-NEC', 'CSS-Healthy', 'CSS-NEC'))

# add factor to allow plotting by group-NEC
#   do for absolute values also
c.abs <- cbind(c.abs, 'group-NEC' = paste(c.abs[, 'group'], c.abs[, 'hist.nec'], sep = '-'))
levels(c.abs$'group-NEC')[levels(c.abs$'group-NEC') == 'NEW-NA'] <- 'NEW'
c.abs$'group-NEC' <- factor(c.abs$'group-NEC', levels = c('NEW', 'LAC-Healthy', 'LAC-NEC', 'MIX-Healthy', 'MIX-NEC', 'CSS-Healthy', 'CSS-NEC'))
p.abs <- cbind(p.abs, 'group-NEC' = paste(p.abs[, 'group'], p.abs[, 'hist.nec'], sep = '-'))
levels(p.abs$'group-NEC')[levels(p.abs$'group-NEC') == 'NEW-NA'] <- 'NEW'
p.abs$'group-NEC' <- factor(p.abs$'group-NEC', levels = c('NEW', 'LAC-Healthy', 'LAC-NEC', 'MIX-Healthy', 'MIX-NEC', 'CSS-Healthy', 'CSS-NEC'))


# plot all 7 groups

plot.df <- cecal  # cecal or plasma
plot.cmpd <- 'palmitoyl ethanolamide' # names(sort(c.anova[1, ]))[5]  # 'beta-muricholate'
ggplot(plot.df, aes(x = plot.df$'group-NEC', y = plot.df[, plot.cmpd], color = plot.df$'group-NEC')) +
  geom_boxplot(lwd = 1, position = position_dodge(width = 0.85), width = 0.7, outlier.shape = NA) +
  geom_point(shape = 19, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.85)) +
  scale_colour_manual(values = group_NEC_colors) + theme(legend.position = 'none') +
  labs(list(title = plot.cmpd, x = element_blank(), y = 'Relative Concentration')) +
  theme(axis.text.x = element_text(angle = -30, hjust = 0))

# 
# # plot Healthy/NEC and Diet group plots side-by-side
# plot.df <- p.fed  # c.fed or p.fed
# top <- names(sort(p.reg.tt))[1:10]  # c.reg.tt or p.reg.tt
# plot.cmpd <- 'lactate'# top[10]
# gg1 <- ggplot(plot.df, aes(x = plot.df$hist.nec, y = plot.df[, plot.cmpd], color = plot.df$hist.nec)) +
#   geom_boxplot(lwd = 1, position = position_dodge(width = 0.85), width = 0.7, outlier.shape = NA) +
#   geom_point(shape = 19, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.85)) +
#   scale_colour_manual(values = hn_colors) + theme(legend.position = 'none') +
#   labs(list(title = element_blank(), x = element_blank(), y = 'Relative Concentration'))
# gg2 <- ggplot(plot.df, aes(x = plot.df$group, y = plot.df[, plot.cmpd], color = plot.df$group)) +
#   geom_boxplot(lwd = 1, position = position_dodge(width = 0.85), width = 0.7, outlier.shape = NA) +
#   geom_point(shape = 19, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.85)) +
#   scale_colour_manual(values = diet_colors) + theme(legend.position = 'none') +
#   labs(list(title = element_blank(), x = element_blank(), y = 'Relative Concentration'))
# grid.arrange(gg1, gg2, top = plot.cmpd, ncol = 2)  



# plot heatmaps
library(RColorBrewer)
library(pheatmap)

c.fed <- droplevels(cecal[cecal$group %in% c('LAC', 'MIX', 'CSS'), ])
p.fed <- droplevels(plasma[plasma$group %in% c('LAC', 'MIX', 'CSS'), ])

# cecal
# get top metabolites to plot, according to type II SS for hist.nec
c.anova <- as.data.frame(apply(c.fed[, as.character(c.compIDs.names)], 2,
                               function(x) anova(lm(x ~ c.fed$group * c.fed$hist.nec))$'Pr(>F)'[2:3]))
c.anova[1, ] <- p.adjust(c.anova[1, ], method = 'fdr')
c.metabs.plot <- names(c.anova[which(c.anova[1, ] < 0.01)])[1:50]

# plot heatmap
c.mat <- t(as.matrix.data.frame(cecal[order(cecal$hist.nec, cecal$group), c.metabs.plot]))
colnames(c.mat) <- cecal[order(cecal$hist.nec, cecal$group), 'pigID']
annot.df <- cecal[, c('hist.nec', 'group')]
annot.df[, 'hist.nec'] <- as.character(annot.df[, 'hist.nec'])
annot.df[, 'group'] <- as.character(annot.df[, 'group'])
rownames(annot.df) <- cecal$pigID
ann.colors <- list(hist.nec = c('Healthy' = hn_colors[1], 'NEC' = hn_colors[2]),
                   group = c('LAC' = diet_colors[1], 'MIX' = diet_colors[2], 'CSS' = diet_colors[3]))
pheatmap(c.mat,
         cluster_rows = T, 
         cluster_cols = F,
         color = rev(colorRampPalette(brewer.pal(8, "RdBu"))(256)),
         annotation = annot.df,
         annotation_colors = ann.colors,
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         border_color = 'gray',
         scale = 'none',
         gaps_col = length(which(annot.df[, 'hist.nec'] == 'Healthy')),
         show_colnames = F,
         main = 'Cecal')


# plasma
# get top metabolites to plot, according to type II SS for hist.nec
p.anova <- as.data.frame(apply(p.fed[, as.character(p.compIDs.names)], 2,
                               function(x) anova(lm(x ~ p.fed$group * p.fed$hist.nec))$'Pr(>F)'[2:3]))
p.anova[1, ] <- p.adjust(p.anova[1, ], method = 'fdr')
p.metabs.plot <- names(p.anova[which(p.anova[1, ] < 0.01)])[1:50]

# plot heatmap
p.mat <- t(as.matrix.data.frame(plasma[order(plasma$hist.nec, plasma$group), p.metabs.plot]))
colnames(p.mat) <- plasma[order(plasma$hist.nec, plasma$group), 'pigID']
annot.df <- plasma[, c('hist.nec', 'group')]
annot.df[, 'hist.nec'] <- as.character(annot.df[, 'hist.nec'])
annot.df[, 'group'] <- as.character(annot.df[, 'group'])
rownames(annot.df) <- plasma$pigID
ann.colors <- list(hist.nec = c('Healthy' = hn_colors[1], 'NEC' = hn_colors[2]),
                   group = c('LAC' = diet_colors[1], 'MIX' = diet_colors[2], 'CSS' = diet_colors[3]))
pheatmap(p.mat,
         cluster_rows = T, 
         cluster_cols = F,
         color = rev(colorRampPalette(brewer.pal(8, "RdBu"))(256)),
         annotation = annot.df,
         annotation_colors = ann.colors,
         clustering_distance_rows = 'euclidean',
         clustering_method = 'ward.D',
         border_color = 'gray',
         scale = 'none',
         gaps_col = length(which(annot.df[, 'hist.nec'] == 'Healthy')),
         show_colnames = F,
         main = 'Plasma')

# 
# 
# # analysis =========================================================================
# # cecal<-cecal[,-c(193,218)]  # these 2 compounds have only 1 sample with detectable amounts
# 
# # to just compare healthy LAC v CSS
# cecal<-cecal[which(cecal[,'hist.nec']=='Healthy'),]
# cecal<-cecal[!cecal$group=='MIX',]
# 
# # remember to change name of output files accordingly
# # cecal_split<-split.data.frame(cecal,cecal$group)
# # cecal_CSS<-cecal_split$CSS
# # cecal_LAC<-cecal_split$LAC
# # rm(cecal_split)
# # 
# # cecal<-cecal_CSS
# # 
# bile<-cmpds[grep('bile',cmpds[,'sub.pathway'],ignore.case=T),]
# cecal<-cbind(cecal[1:8],cecal[,unlist(lapply(colnames(cecal),function(x){substr(x,2,nchar(x))})) %in% bile[,'comp.id']])
# # edit above code block to select different groups and compounds
# 
# n.healthy<-apply(cecal[,9:ncol(cecal)],2,function(x){length(which(!is.na(x[cecal$hist.nec=='Healthy'])))})
# n.nec<-apply(cecal[,9:ncol(cecal)],2,function(x){length(which(!is.na(x[cecal$hist.nec=='NEC'])))})
# 
# # impute missing values with minimum
# colmins<-apply(cecal[,9:ncol(cecal)],2,function(x){min(x,na.rm=T)})
# for(r in 1:nrow(cecal)) {
#   for (c in 9:ncol(cecal)) {
#     if(is.na(cecal[r,c])) 
#       cecal[r,c]<-colmins[c-8]
#   }
# }
# rm(c,r,colmins)
# 
# tt<-apply(cecal[,9:ncol(cecal)],2,function(x){tryCatch(t.test(x~cecal$hist.nec,var.equal=F)$p.value, error=function(x) NA)}) # Welch's t-test; could also use Mann-Whitney U (wilcox.test)
# tt<-p.adjust(tt,method='fdr')
# mean.healthy<-apply(cecal[,9:ncol(cecal)],2,function(x){mean(x[cecal$hist.nec=='Healthy'],na.rm=T)})
# mean.nec<-apply(cecal[,9:ncol(cecal)],2,function(x){mean(x[cecal$hist.nec=='NEC'],na.rm=T)})
# stdev.healthy<-apply(cecal[,9:ncol(cecal)],2,function(x){sd(x[cecal$hist.nec=='Healthy'],na.rm=T)})
# stdev.nec<-apply(cecal[,9:ncol(cecal)],2,function(x){sd(x[cecal$hist.nec=='NEC'],na.rm=T)})
# 
# tt<-cbind(tt,mean.healthy,stdev.healthy,n.healthy,mean.nec,stdev.nec,n.nec)
# 
# tt<-cbind(unlist(lapply(row.names(tt),function(x){substr(x,2,nchar(x))})),tt)
# colnames(tt)[1]<-'comp.id'
# tt<-merge(tt,cmpds,by='comp.id') # add compound info
# tt<-tt[,c(2:8,10:12,1,9,13:ncol(tt))] # reorder columns
# tt[,1]<-as.numeric(as.character(tt[,1])) # needed in order to sort properly
# 
# write.table(matrix(c('# q-value','','','N non-imp','','','N non-imp'),nrow=1),paste(path,'cecal_NECvHealthy_raw.csv',sep=''),row.names=F,col.names=F,sep=',')
# write.table(tt[order(tt[,1]),,drop=F],paste(path,'cecal_NECvHealthy_raw.csv',sep=''),row.names=F,col.names=T,quote=F,sep=',',append=T)
# 
# rm(list=ls())
# 
# # ==================================================================================