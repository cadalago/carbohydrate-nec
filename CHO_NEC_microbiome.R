## Microbiome stats and figs for CHO/NEC study
library('phyloseq')
library('ggplot2')
library('multcomp')
library('reshape2')
library('gridExtra')

# set plotting theme
theme_set(theme_classic())
group_colors=c('#000000', '#4daf4a', '#377eb8', '#e41a1c')    # aka 'black' (previously 'pink'(ffc0cb)), 'palegreen2', 'goldenrod1', 'red3'
hn_colors<-c('turquoise3', 'orangered1')

# define file locations
biom_file <- '../data/microbes/latest/str_MJN_w_newborn_samples/str_MJN_w_newborn_samples_removed___.biom' # ...samples_contents.biom or ...samples_tissue.biom
tree_file <- '../data/microbes/latest/silva123_V4.tre'
mapping_file <- '../data/microbes/latest/mapping_files/MJN_w_newborn_samples_mappingfile_allmd_new_as_healthy_phyloseq.csv'

# create phyloseq object from input files
ps <- import_biom(biom_file, tree_file)
sample_names(ps) <- make.names(sample_names(ps))
sd <- sample_data(read.csv(mapping_file, row.names = 1))
row.names(sd) <- make.names(row.names(sd))
ps <- merge_phyloseq(ps, sd)
colnames(tax_table(ps)) <- c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
sample_data(ps)$group <- factor(sample_data(ps)$group, levels = c('NEW', 'LAC', 'MIX', 'CSS'))
sample_data(ps)$site <- factor(sample_data(ps)$site, levels = c('Stomach', 'Ileum', 'Colon'))
sample_data(ps)$NEC.group <- factor(sample_data(ps)$NEC.group, levels = c('NEW', 'Healthy-LAC', 'Healthy-MIX', 'Healthy-CSS', 'NEC-LAC', 'NEC-MIX', 'NEC-CSS'))


# remove spurious taxa
ps <- subset_taxa(ps, Domain == 'Bacteria')

# Make a data frame with a column for the read counts of each sample
windows()
sample_sum_df <- data.frame(sum = sample_sums(ps)) # try: prune_samples(sample_sums(ps) >= 1000, ps)))

# Histogram of sample read counts
#   by sample type (tissue v contents)
ggplot(merge(sample_sum_df, sample_data(ps), by = 'row.names'), aes(x = sum)) +  facet_grid( ~ type) +
  geom_histogram(binwidth = 500) +
  ggtitle('Distribution of sample sequencing depth (red line = 1000)') + 
  xlab('Read counts (bin=500)') +
  xlim(-500, NA) +  # -500 to ensure all data are displayed on chart; NA to set limit to data max
  geom_vline(xintercept = 1000, linetype = 'dashed', color = 'red') +
  theme(axis.title.y = element_blank())
#   by all categories (site, type, group, NEC)
ggplot(merge(sample_sum_df, sample_data(ps), by = 'row.names'), aes(x = sum)) +  facet_grid(type + group ~ site + NEC) +
  geom_histogram(binwidth = 500) + theme_bw() +
  ggtitle('Distribution of sample sequencing depth (red line = 1000)') + 
  xlab('Read counts (bin=500)') +
  xlim(-500, NA) +  # -500 to ensure all data are displayed on chart; NA to set limit to data max
  geom_vline(xintercept = 1000, linetype = 'dashed', color = 'red') +
  theme(axis.title.y = element_blank())

# remove samples with low read counts
ps <- prune_samples(sample_sums(ps) >= 1000, ps)

# subset into tissue and contents
ps.t <- subset_samples(ps, type == 'Tissue')
ps.c <- subset_samples(ps, type == 'Contents')

# subset into Healthy and NEC
ps.t.h <- subset_samples(ps.t, NEC == 'Healthy')
ps.t.n <- subset_samples(ps.t, NEC == 'NEC')
ps.c.h <- subset_samples(ps.c, NEC == 'Healthy')
ps.c.n <- subset_samples(ps.c, NEC == 'NEC')

# remove OTUs which are present in 0 samples
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
ps.t <- prune_taxa(taxa_sums(ps.t) > 0, ps.t)
ps.t.h <- prune_taxa(taxa_sums(ps.t.h) > 0, ps.t.h)
ps.t.n <- prune_taxa(taxa_sums(ps.t.n) > 0, ps.t.n)
ps.c <- prune_taxa(taxa_sums(ps.c) > 0, ps.c)
ps.c.h <- prune_taxa(taxa_sums(ps.c.h) > 0, ps.c.h)
ps.c.n <- prune_taxa(taxa_sums(ps.c.n) > 0, ps.c.n)

# plot richness and alpha diversity
# rich.c.n <- plot_richness(ps.c.h, x='site', measures = c('Observed', 'InvSimpson'), color = 'group', shape = NA) +
#   geom_boxplot(lwd = 1, position = position_dodge(width = 0.85), width = 0.7, outlier.shape = NA) + #facet_grid(. ~ NEC) +
#   geom_point(shape = 19, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.85)) +
#   scale_colour_manual(values = group_colors) +
#   labs(list(title = 'Healthy, Contents', x = element_blank(), y = element_blank())) +
#   ylim(0, NA) +
#   guides(color=guide_legend())

rich.c <- merge(estimate_richness(ps.c, measures = c('Observed', 'InvSimpson')), sample_data(ps.c), by = 'row.names')
rich.t <- merge(estimate_richness(ps.t, measures = c('Observed', 'InvSimpson')), sample_data(ps.t), by = 'row.names')
rich.c.h <- merge(estimate_richness(ps.c.h, measures = c('Observed', 'InvSimpson')), sample_data(ps.c.h), by = 'row.names')
rich.c.n <- merge(estimate_richness(ps.c.n, measures = c('Observed', 'InvSimpson')), sample_data(ps.c.n), by = 'row.names')
rich.t.h <- merge(estimate_richness(ps.t.h, measures = c('Observed', 'InvSimpson')), sample_data(ps.t.h), by = 'row.names')
rich.t.n <- merge(estimate_richness(ps.t.n, measures = c('Observed', 'InvSimpson')), sample_data(ps.t.n), by = 'row.names')

gg1 <- ggplot(rich.c.h, aes(x = site, y = Observed, color = group)) +
  geom_boxplot(lwd = 1, position = position_dodge(width = 0.85), width = 0.7, outlier.shape = NA) +
  geom_point(shape = 19, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.85)) +
  scale_colour_manual(values = group_colors) +
  labs(list(title = 'Obs, Healthy, Contents', x = element_blank(), y = element_blank())) +
  ylim(0, 80) + guides(color = FALSE)  # shape=19, then copy as metafile, paste special "windows metafile", ungroup
gg2 <- ggplot(rich.c.n, aes(x = site, y = Observed, color = group)) +
  geom_boxplot(lwd = 1, position = position_dodge(width = 0.85), width = 0.7, outlier.shape = NA) +
  geom_point(shape = 19, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.85)) +
  scale_colour_manual(values = group_colors[2:4]) + scale_y_continuous(breaks = NULL, limits = c(0, 80)) +
  labs(list(title = 'Obs, NEC, Contents', x = element_blank(), y = element_blank()))
grid.arrange(gg1, gg2, ncol = 2)  
#plot_grid(gg.rich.c.h, gg.rich.c.n)

gg1 <- ggplot(rich.c.h, aes(x = site, y = InvSimpson, color = group)) +
  geom_boxplot(lwd = 1, position = position_dodge(width = 0.85), width = 0.7, outlier.shape = NA) +
  geom_point(shape = 19, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.85)) +
  scale_colour_manual(values = group_colors) +
  labs(list(title = 'InvS, Healthy, Contents', x = element_blank(), y = element_blank())) +
  ylim(0, 15) + guides(color = FALSE)
gg2 <- ggplot(rich.c.n, aes(x = site, y = InvSimpson, color = group)) +
  geom_boxplot(lwd = 1, position = position_dodge(width = 0.85), width = 0.7, outlier.shape = NA) +
  geom_point(shape = 19, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.85)) +
  scale_colour_manual(values = group_colors[2:4]) + scale_y_continuous(breaks = NULL, limits = c(0, 15)) +
  labs(list(title = 'InvS, NEC, Contents', x = element_blank(), y = element_blank()))
grid.arrange(gg1, gg2, ncol = 2)

gg1 <- ggplot(rich.t.h, aes(x = site, y = Observed, color = group)) +
  geom_boxplot(lwd = 1, position = position_dodge(width = 0.85), width = 0.7, outlier.shape = NA) +
  geom_point(shape = 19, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.85)) +
  scale_colour_manual(values = group_colors) +
  labs(list(title = 'Obs, Healthy, Tissue', x = element_blank(), y = element_blank())) +
  ylim(0, 80) + guides(color = FALSE)
gg2 <- ggplot(rich.t.n, aes(x = site, y = Observed, color = group)) +
  geom_boxplot(lwd = 1, position = position_dodge(width = 0.85), width = 0.7, outlier.shape = NA) +
  geom_point(shape = 19, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.85)) +
  scale_colour_manual(values = group_colors[2:4]) + scale_y_continuous(breaks = NULL, limits = c(0, 80)) +
  labs(list(title = 'Obs, NEC, Tissue', x = element_blank(), y = element_blank()))
grid.arrange(gg1, gg2, ncol = 2)  

gg1 <- ggplot(rich.t.h, aes(x = site, y = InvSimpson, color = group)) +
  geom_boxplot(lwd = 1, position = position_dodge(width = 0.85), width = 0.7, outlier.shape = NA) +
  geom_point(shape = 19, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.85)) +
  scale_colour_manual(values = group_colors) +
  labs(list(title = 'InvS, Healthy, Tissue', x = element_blank(), y = element_blank())) +
  ylim(0, 15) + guides(color = FALSE)
gg2 <- ggplot(rich.t.n, aes(x = site, y = InvSimpson, color = group)) +
  geom_boxplot(lwd = 1, position = position_dodge(width = 0.85), width = 0.7, outlier.shape = NA) +
  geom_point(shape = 19, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.85)) +
  scale_colour_manual(values = group_colors[2:4]) + scale_y_continuous(breaks = NULL, limits = c(0, 15)) +
  labs(list(title = 'InvS, NEC, Tissue', x = element_blank(), y = element_blank()))
grid.arrange(gg1, gg2, ncol = 2)

#gg.rich.c <- ggplot(rich.c, aes(x = site, y = Observed, color = group)) + facet_grid(NEC ~ ., scale = 'free_x') +
#  geom_boxplot(lwd = 1, position = position_dodge(width = 0.85), width = 0.7, outlier.shape = NA) +
#  geom_point(shape = 19, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.85)) +
#  scale_colour_manual(values = group_colors) +
#  labs(list(title = 'Observed, Contents', x = element_blank(), y = element_blank())) +
#  ylim(0, NA)
  
# test for significant differences in richness and alpha diversity
summary(glht(glm(Observed ~ group, family = gaussian, data = rich.c.h), linfct = mcp(group = 'Tukey')))
summary(glht(glm(Observed ~ group, family = gaussian, data = rich.c.n), linfct = mcp(group = 'Tukey')))
summary(glht(glm(InvSimpson ~ group, family = gaussian, data = rich.c.h), linfct = mcp(group = 'Tukey')))
summary(glht(glm(InvSimpson ~ group, family = gaussian, data = rich.c.n), linfct = mcp(group = 'Tukey')))
summary(glht(glm(Observed ~ group, family = gaussian, data = rich.t.h), linfct = mcp(group = 'Tukey')))
summary(glht(glm(Observed ~ group, family = gaussian, data = rich.t.n), linfct = mcp(group = 'Tukey')))
summary(glht(glm(InvSimpson ~ group, family = gaussian, data = rich.t.h), linfct = mcp(group = 'Tukey')))
summary(glht(glm(InvSimpson ~ group, family = gaussian, data = rich.t.n), linfct = mcp(group = 'Tukey')))

# # prevalence exploration
# # compute prevalence of each feature, store as data.frame
# prevdf = apply(X = otu_table(ps),
#                MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
#                FUN = function(x){sum(x > 0)})
# # Add taxonomy and total read counts to this data.frame
# prevdf = data.frame(Prevalence = prevdf,
#                     TotalAbundance = taxa_sums(ps),
#                     tax_table(ps))
# 
# # Filter taxa based on above computed prevalence
# df2 <- plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
# filter_phyla <- df2[df2$'1' <= 1, 'Phylum']
# ps1 <- subset_taxa(ps, !Phylum %in% filter_phyla)
# 
# prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps1, "Phylum"))
# ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
#   # Include a guess for parameter
#   geom_hline(yintercept = 0.1, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
#   scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
#   facet_wrap(~Phylum) + theme(legend.position="none")
# 

# Group comparisons, across diets (only healthy, incl. NEW)

# subset by GI tract site
ps.t.h.il <- subset_samples(ps.t.h, site == 'Ileum')
ps.t.h.co <- subset_samples(ps.t.h, site == 'Colon')
ps.t.h.ilco <- subset_samples(ps.t.h, site %in% c('Ileum', 'Colon'))
ps.c.h.st <- subset_samples(ps.c.h, site == 'Stomach')
ps.c.h.il <- subset_samples(ps.c.h, site == 'Ileum')
ps.c.h.co <- subset_samples(ps.c.h, site == 'Colon')
ps.c.h.ilco <- subset_samples(ps.c.h, site %in% c('Ileum', 'Colon'))

# transform abundances from absolute to relative (proportions)
#ps.c.h.co <- subset_samples(ps.c.h, site == 'Colon')
ps.t.h.il.rel <- transform_sample_counts(ps.t.h.il, function(x) x / sum(x))
ps.t.h.co.rel <- transform_sample_counts(ps.t.h.co, function(x) x / sum(x))
ps.t.h.ilco.rel <- transform_sample_counts(ps.t.h.ilco, function(x) x / sum(x))
ps.c.h.st.rel <- transform_sample_counts(ps.c.h.st, function(x) x / sum(x))
ps.c.h.il.rel <- transform_sample_counts(ps.c.h.il, function(x) x / sum(x))
ps.c.h.co.rel <- transform_sample_counts(ps.c.h.co, function(x) x / sum(x))
ps.c.h.ilco.rel <- transform_sample_counts(ps.c.h.ilco, function(x) x / sum(x))

# remove OTUs whose relative abundance is not at least 0.01 in at least one sample
#    (i.e., remove OTUs whose max relative abundance is less than 0.01)
ps.t.h.il.rel.fil <- prune_taxa(apply(otu_table(ps.t.h.il.rel), 1, max) > 0.01, ps.t.h.il.rel)
ps.t.h.co.rel.fil <- prune_taxa(apply(otu_table(ps.t.h.co.rel), 1, max) > 0.01, ps.t.h.co.rel)
ps.t.h.ilco.rel.fil <- prune_taxa(apply(otu_table(ps.t.h.ilco.rel), 1, max) > 0.01, ps.t.h.ilco.rel)
ps.c.h.st.rel.fil <- prune_taxa(apply(otu_table(ps.c.h.st.rel), 1, max) > 0.01, ps.c.h.st.rel)
ps.c.h.il.rel.fil <- prune_taxa(apply(otu_table(ps.c.h.il.rel), 1, max) > 0.01, ps.c.h.il.rel)
ps.c.h.co.rel.fil <- prune_taxa(apply(otu_table(ps.c.h.co.rel), 1, max) > 0.01, ps.c.h.co.rel)
ps.c.h.ilco.rel.fil <- prune_taxa(apply(otu_table(ps.c.h.ilco.rel), 1, max) > 0.01, ps.c.h.ilco.rel)

# agglomerate OTUs by genus
names.b.e = c('NA', '', ' ', '\t', 'uncultured', 'g', 'Incertae_Sedis')
ps.t.h.il.rel.fil.genus <- tax_glom(ps.t.h.il.rel.fil, 'Genus', bad_empty = names.b.e)
ps.t.h.co.rel.fil.genus <- tax_glom(ps.t.h.co.rel.fil, 'Genus', bad_empty = names.b.e)
ps.t.h.ilco.rel.fil.genus <- tax_glom(ps.t.h.ilco.rel.fil, 'Genus', bad_empty = names.b.e)
ps.c.h.st.rel.fil.genus <- tax_glom(ps.c.h.st.rel.fil, 'Genus', bad_empty = names.b.e)
ps.c.h.il.rel.fil.genus <- tax_glom(ps.c.h.il.rel.fil, 'Genus', bad_empty = names.b.e)
ps.c.h.co.rel.fil.genus <- tax_glom(ps.c.h.co.rel.fil, 'Genus', bad_empty = names.b.e)
ps.c.h.ilco.rel.fil.genus <- tax_glom(ps.c.h.ilco.rel.fil, 'Genus', bad_empty = names.b.e)


# plot heatmaps

sort.df <- get_variable(ps.t.h.il.rel.fil.genus, c('pigID', 'NEC.group'))
p <- plot_heatmap(ps.t.h.il.rel.fil.genus, 'PCoA', 'jaccard', sample.label = 'NEC.group', taxa.label = 'Genus',
             sample.order = row.names(sort.df[order(sort.df$NEC.group), ]), low = '#C0FFFF', high = '#000099', na.value = 'white', 
             title = 'Healthy, Tissue, Ileum') + theme(text = element_text(size = 10),
                                                       axis.text.x = element_text(size = 10), 
                                                       axis.text.y = element_text(size = 10))
ggsave('D:/Desktop/rstudio_out/healthy_tissue_ileum.pdf', p, units = 'in',
       width = 0.25 * dim(otu_table(ps.t.h.il.rel.fil.genus))[2] + 0.1 * max(nchar(tax_table(ps.t.h.il.rel.fil.genus)[, 'Genus'])),
       height = 0.2 * dim(otu_table(ps.t.h.il.rel.fil.genus))[1])

sort.df <- get_variable(ps.t.h.co.rel.fil.genus, c('pigID', 'NEC.group'))
p <- plot_heatmap(ps.t.h.co.rel.fil.genus, 'PCoA', 'jaccard', sample.label = 'NEC.group', taxa.label = 'Genus',
             sample.order = row.names(sort.df[order(sort.df$NEC.group), ]), low = '#C0FFFF', high = '#000099', na.value = 'white', 
             title = 'Healthy, Tissue, Colon') +  theme(text = element_text(size = 10),
                                                        axis.text.x = element_text(size = 10), 
                                                        axis.text.y = element_text(size = 10))
ggsave('D:/Desktop/rstudio_out/healthy_tissue_colon.pdf', p, units = 'in',
       width = 0.25 * dim(otu_table(ps.t.h.co.rel.fil.genus))[2] + 0.1 * max(nchar(tax_table(ps.t.h.co.rel.fil.genus)[, 'Genus'])),
       height = 0.2 * dim(otu_table(ps.t.h.co.rel.fil.genus))[1])

sort.df <- get_variable(ps.t.h.ilco.rel.fil.genus, c('pigID', 'NEC.group'))
p <- plot_heatmap(ps.t.h.ilco.rel.fil.genus, 'PCoA', 'jaccard', sample.label = 'NEC.group', taxa.label = 'Genus',
             sample.order = row.names(sort.df[order(sort.df$NEC.group), ]), low = '#C0FFFF', high = '#000099', na.value = 'white', 
             title = 'Healthy, Tissue, Ileum + Colon') +  theme(text = element_text(size = 10),
                                                                axis.text.x = element_text(size = 10), 
                                                                axis.text.y = element_text(size = 10))
ggsave('D:/Desktop/rstudio_out/healthy_tissue_ileum_colon.pdf', p, units = 'in',
       width = 0.25 * dim(otu_table(ps.t.h.ilco.rel.fil.genus))[2] + 0.1 * max(nchar(tax_table(ps.t.h.ilco.rel.fil.genus)[, 'Genus'])),
       height = 0.2 * dim(otu_table(ps.t.h.ilco.rel.fil.genus))[1])

sort.df <- get_variable(ps.c.h.st.rel.fil.genus, c('pigID', 'NEC.group'))
p <- plot_heatmap(ps.c.h.st.rel.fil.genus, 'PCoA', 'jaccard', sample.label = 'NEC.group', taxa.label = 'Genus',
             sample.order = row.names(sort.df[order(sort.df$NEC.group), ]), low = '#C0FFFF', high = '#000099', na.value = 'white', 
             title = 'Healthy, Contents, Stomach') +  theme(text = element_text(size = 10),
                                                            axis.text.x = element_text(size = 10), 
                                                            axis.text.y = element_text(size = 10))
ggsave('D:/Desktop/rstudio_out/healthy_contents_stom.pdf', p, units = 'in',
       width = 0.25 * dim(otu_table(ps.c.h.st.rel.fil.genus))[2] + 0.1 * max(nchar(tax_table(ps.c.h.st.rel.fil.genus)[, 'Genus'])),
       height = 0.2 * dim(otu_table(ps.c.h.st.rel.fil.genus))[1])

sort.df <- get_variable(ps.c.h.il.rel.fil.genus, c('pigID', 'NEC.group'))
p <- plot_heatmap(ps.c.h.il.rel.fil.genus, 'PCoA', 'jaccard', sample.label = 'NEC.group', taxa.label = 'Genus',
             sample.order = row.names(sort.df[order(sort.df$NEC.group), ]), low = '#C0FFFF', high = '#000099', na.value = 'white', 
             title = 'Healthy, Contents, Ileum') +  theme(text = element_text(size = 10),
                                                          axis.text.x = element_text(size = 10), 
                                                          axis.text.y = element_text(size = 10))
ggsave('D:/Desktop/rstudio_out/healthy_contents_ileum.pdf', p, units = 'in',
       width = 0.25 * dim(otu_table(ps.c.h.il.rel.fil.genus))[2] + 0.1 * max(nchar(tax_table(ps.c.h.il.rel.fil.genus)[, 'Genus'])),
       height = 0.2 * dim(otu_table(ps.c.h.il.rel.fil.genus))[1])

sort.df <- get_variable(ps.c.h.co.rel.fil.genus, c('pigID', 'NEC.group'))
p <- plot_heatmap(ps.c.h.co.rel.fil.genus, 'PCoA', 'jaccard', sample.label = 'NEC.group', taxa.label = 'Genus',
             sample.order = row.names(sort.df[order(sort.df$NEC.group), ]), low = '#C0FFFF', high = '#000099', na.value = 'white', 
             title = 'Healthy, Contents, Colon') +  theme(text = element_text(size = 10),
                                                          axis.text.x = element_text(size = 10), 
                                                          axis.text.y = element_text(size = 10))
ggsave('D:/Desktop/rstudio_out/healthy_contents_colon.pdf', p, units = 'in',
       width = 0.25 * dim(otu_table(ps.c.h.co.rel.fil.genus))[2] + 0.1 * max(nchar(tax_table(ps.c.h.co.rel.fil.genus)[, 'Genus'])),
       height = 0.2 * dim(otu_table(ps.c.h.co.rel.fil.genus))[1])

sort.df <- get_variable(ps.c.h.ilco.rel.fil.genus, c('pigID', 'NEC.group'))
p <- plot_heatmap(ps.c.h.ilco.rel.fil.genus, 'PCoA', 'jaccard', sample.label = 'NEC.group', taxa.label = 'Genus',
             sample.order = row.names(sort.df[order(sort.df$NEC.group), ]), low = '#C0FFFF', high = '#000099', na.value = 'white', 
             title = 'Healthy, Contents, Ileum + Colon') +  theme(text = element_text(size = 10),
                                                                  axis.text.x = element_text(size = 10), 
                                                                  axis.text.y = element_text(size = 10))
ggsave('D:/Desktop/rstudio_out/healthy_contents_ileum_colon.pdf', p, units = 'in',
       width = 0.25 * dim(otu_table(ps.c.h.ilco.rel.fil.genus))[2] + 0.1 * max(nchar(tax_table(ps.c.h.ilco.rel.fil.genus)[, 'Genus'])),
       height = 0.2 * dim(otu_table(ps.c.h.ilco.rel.fil.genus))[1])


# phenotype comparisons, Healthy v NEC

# subset by GI tract site
ps.t.il <- subset_samples(ps.t, site == 'Ileum')
ps.t.co <- subset_samples(ps.t, site == 'Colon')
ps.t.ilco <- subset_samples(ps.t, site %in% c('Ileum', 'Colon'))
ps.c.st <- subset_samples(ps.c, site == 'Stomach')
ps.c.il <- subset_samples(ps.c, site == 'Ileum')
ps.c.co <- subset_samples(ps.c, site == 'Colon')
ps.c.ilco <- subset_samples(ps.c, site %in% c('Ileum', 'Colon'))

# transform abundances from absolute to relative (proportions)
ps.t.il.rel <- transform_sample_counts(ps.t.il, function(x) x / sum(x))
ps.t.co.rel <- transform_sample_counts(ps.t.co, function(x) x / sum(x))
ps.t.ilco.rel <- transform_sample_counts(ps.t.ilco, function(x) x / sum(x))
ps.c.st.rel <- transform_sample_counts(ps.c.st, function(x) x / sum(x))
ps.c.il.rel <- transform_sample_counts(ps.c.il, function(x) x / sum(x))
ps.c.co.rel <- transform_sample_counts(ps.c.co, function(x) x / sum(x))
ps.c.ilco.rel <- transform_sample_counts(ps.c.ilco, function(x) x / sum(x))

# remove OTUs whose relative abundance is not at least 0.01 in at least one sample
#    (i.e., remove OTUs whose max relative abundance is less than 0.01)
ps.t.il.rel.fil <- prune_taxa(apply(otu_table(ps.t.il.rel), 1, max) > 0.01, ps.t.il.rel)
ps.t.co.rel.fil <- prune_taxa(apply(otu_table(ps.t.co.rel), 1, max) > 0.01, ps.t.co.rel)
ps.t.ilco.rel.fil <- prune_taxa(apply(otu_table(ps.t.ilco.rel), 1, max) > 0.01, ps.t.ilco.rel)
ps.c.st.rel.fil <- prune_taxa(apply(otu_table(ps.c.st.rel), 1, max) > 0.01, ps.c.st.rel)
ps.c.il.rel.fil <- prune_taxa(apply(otu_table(ps.c.il.rel), 1, max) > 0.01, ps.c.il.rel)
ps.c.co.rel.fil <- prune_taxa(apply(otu_table(ps.c.co.rel), 1, max) > 0.01, ps.c.co.rel)
ps.c.ilco.rel.fil <- prune_taxa(apply(otu_table(ps.c.ilco.rel), 1, max) > 0.01, ps.c.ilco.rel)

# agglomerate OTUs by genus
names.b.e = c('NA', '', ' ', '\t', 'uncultured', 'g', 'Incertae_Sedis')
ps.t.il.rel.fil.genus <- tax_glom(ps.t.il.rel.fil, 'Genus', bad_empty = names.b.e)
ps.t.co.rel.fil.genus <- tax_glom(ps.t.co.rel.fil, 'Genus', bad_empty = names.b.e)
ps.t.ilco.rel.fil.genus <- tax_glom(ps.t.ilco.rel.fil, 'Genus', bad_empty = names.b.e)
ps.c.st.rel.fil.genus <- tax_glom(ps.c.st.rel.fil, 'Genus', bad_empty = names.b.e)
ps.c.il.rel.fil.genus <- tax_glom(ps.c.il.rel.fil, 'Genus', bad_empty = names.b.e)
ps.c.co.rel.fil.genus <- tax_glom(ps.c.co.rel.fil, 'Genus', bad_empty = names.b.e)
ps.c.ilco.rel.fil.genus <- tax_glom(ps.c.ilco.rel.fil, 'Genus', bad_empty = names.b.e)

# plot heatmaps
sort.df <- get_variable(ps.t.il.rel.fil.genus, c('pigID', 'NEC.group'))
p <- plot_heatmap(ps.t.il.rel.fil.genus, 'PCoA', 'jaccard', sample.label = 'NEC.group', taxa.label = 'Genus',
                  sample.order = row.names(sort.df[order(sort.df$NEC.group), ]), low = '#C0FFFF', high = '#000099', na.value = 'white', 
                  title = 'All, Tissue, Ileum') + theme(text = element_text(size = 10),
                                                            axis.text.x = element_text(size = 10), 
                                                            axis.text.y = element_text(size = 10))
ggsave('D:/Desktop/rstudio_out/All_tissue_ileum.pdf', p, units = 'in',
       width = 0.25 * dim(otu_table(ps.t.il.rel.fil.genus))[2] + 0.1 * max(nchar(tax_table(ps.t.il.rel.fil.genus)[, 'Genus'])),
       height = 0.2 * dim(otu_table(ps.t.il.rel.fil.genus))[1])

sort.df <- get_variable(ps.t.co.rel.fil.genus, c('pigID', 'NEC.group'))
p <- plot_heatmap(ps.t.co.rel.fil.genus, 'PCoA', 'jaccard', sample.label = 'NEC.group', taxa.label = 'Genus',
                  sample.order = row.names(sort.df[order(sort.df$NEC.group), ]), low = '#C0FFFF', high = '#000099', na.value = 'white', 
                  title = 'All, Tissue, Colon') +  theme(text = element_text(size = 10),
                                                             axis.text.x = element_text(size = 10), 
                                                             axis.text.y = element_text(size = 10))
ggsave('D:/Desktop/rstudio_out/All_tissue_colon.pdf', p, units = 'in',
       width = 0.25 * dim(otu_table(ps.t.co.rel.fil.genus))[2] + 0.1 * max(nchar(tax_table(ps.t.co.rel.fil.genus)[, 'Genus'])),
       height = 0.2 * dim(otu_table(ps.t.co.rel.fil.genus))[1])

sort.df <- get_variable(ps.t.ilco.rel.fil.genus, c('pigID', 'NEC.group'))
p <- plot_heatmap(ps.t.ilco.rel.fil.genus, 'PCoA', 'jaccard', sample.label = 'NEC.group', taxa.label = 'Genus',
                  sample.order = row.names(sort.df[order(sort.df$NEC.group), ]), low = '#C0FFFF', high = '#000099', na.value = 'white', 
                  title = 'All, Tissue, Ileum + Colon') +  theme(text = element_text(size = 10),
                                                                     axis.text.x = element_text(size = 10), 
                                                                     axis.text.y = element_text(size = 10))
ggsave('D:/Desktop/rstudio_out/All_tissue_ileum_colon.pdf', p, units = 'in',
       width = 0.25 * dim(otu_table(ps.t.ilco.rel.fil.genus))[2] + 0.1 * max(nchar(tax_table(ps.t.ilco.rel.fil.genus)[, 'Genus'])),
       height = 0.2 * dim(otu_table(ps.t.ilco.rel.fil.genus))[1])

sort.df <- get_variable(ps.c.st.rel.fil.genus, c('pigID', 'NEC.group'))
p <- plot_heatmap(ps.c.st.rel.fil.genus, 'PCoA', 'jaccard', sample.label = 'NEC.group', taxa.label = 'Genus',
                  sample.order = row.names(sort.df[order(sort.df$NEC.group), ]), low = '#C0FFFF', high = '#000099', na.value = 'white', 
                  title = 'All, Contents, Stomach') +  theme(text = element_text(size = 10),
                                                                 axis.text.x = element_text(size = 10), 
                                                                 axis.text.y = element_text(size = 10))
ggsave('D:/Desktop/rstudio_out/All_contents_stom.pdf', p, units = 'in',
       width = 0.25 * dim(otu_table(ps.c.st.rel.fil.genus))[2] + 0.1 * max(nchar(tax_table(ps.c.st.rel.fil.genus)[, 'Genus'])),
       height = 0.2 * dim(otu_table(ps.c.st.rel.fil.genus))[1])

sort.df <- get_variable(ps.c.il.rel.fil.genus, c('pigID', 'NEC.group'))
p <- plot_heatmap(ps.c.il.rel.fil.genus, 'PCoA', 'jaccard', sample.label = 'NEC.group', taxa.label = 'Genus',
                  sample.order = row.names(sort.df[order(sort.df$NEC.group), ]), low = '#C0FFFF', high = '#000099', na.value = 'white', 
                  title = 'All, Contents, Ileum') +  theme(text = element_text(size = 10),
                                                               axis.text.x = element_text(size = 10), 
                                                               axis.text.y = element_text(size = 10))
ggsave('D:/Desktop/rstudio_out/All_contents_ileum.pdf', p, units = 'in',
       width = 0.25 * dim(otu_table(ps.c.il.rel.fil.genus))[2] + 0.1 * max(nchar(tax_table(ps.c.il.rel.fil.genus)[, 'Genus'])),
       height = 0.2 * dim(otu_table(ps.c.il.rel.fil.genus))[1])

sort.df <- get_variable(ps.c.co.rel.fil.genus, c('pigID', 'NEC.group'))
p <- plot_heatmap(ps.c.co.rel.fil.genus, 'PCoA', 'jaccard', sample.label = 'NEC.group', taxa.label = 'Genus',
                  sample.order = row.names(sort.df[order(sort.df$NEC.group), ]), low = '#C0FFFF', high = '#000099', na.value = 'white', 
                  title = 'All, Contents, Colon') +  theme(text = element_text(size = 10),
                                                               axis.text.x = element_text(size = 10), 
                                                               axis.text.y = element_text(size = 10))
ggsave('D:/Desktop/rstudio_out/All_contents_colon.pdf', p, units = 'in',
       width = 0.25 * dim(otu_table(ps.c.co.rel.fil.genus))[2] + 0.1 * max(nchar(tax_table(ps.c.co.rel.fil.genus)[, 'Genus'])),
       height = 0.2 * dim(otu_table(ps.c.co.rel.fil.genus))[1])

sort.df <- get_variable(ps.c.ilco.rel.fil.genus, c('pigID', 'NEC.group'))
p <- plot_heatmap(ps.c.ilco.rel.fil.genus, 'PCoA', 'jaccard', sample.label = 'NEC.group', taxa.label = 'Genus',
                  sample.order = row.names(sort.df[order(sort.df$NEC.group), ]), low = '#C0FFFF', high = '#000099', na.value = 'white', 
                  title = 'All, Contents, Ileum + Colon') +  theme(text = element_text(size = 10),
                                                                       axis.text.x = element_text(size = 10), 
                                                                       axis.text.y = element_text(size = 10))
ggsave('D:/Desktop/rstudio_out/All_contents_ileum_colon.pdf', p, units = 'in',
       width = 0.25 * dim(otu_table(ps.c.ilco.rel.fil.genus))[2] + 0.1 * max(nchar(tax_table(ps.c.ilco.rel.fil.genus)[, 'Genus'])),
       height = 0.2 * dim(otu_table(ps.c.ilco.rel.fil.genus))[1])

# 
# # or use pheatmap:
# dfgroup <- data.frame(row.names = sample_names(ps.c.h.rel.fil.genus),
#                       group = sample_data(ps.c.h.rel.fil.genus)$group,
#                       site = sample_data(ps.c.h.rel.fil.genus)$site)
# mdist <- distance(ps.c.h.rel.fil.genus, method = 'jaccard')
# pheatmap(otu_table(ps.c.h.rel.fil.genus), border_color = 'NA', annotation_col = dfgroup,
#          annotation_colors = list(group = c(NEW = '#ffc0cb', LAC = '#4daf4a', MIX = '#377eb8', CSS = '#e41a1c')),
#          show_colnames = F, treeheight_row = 0, treeheight_col = 0, scale = 'none', color = colorRampPalette(c('#000033', '#FF3300'))(100),
#          labels_row = tax_table(ps.c.h.rel.fil.genus)[, 'Genus'], 
#          clustering_distance_cols = mdist, main = 'Healthy, Contents')
# 

# try making sure cell size is equal across heatmaps, as in 
#    http://stackoverflow.com/questions/25937355/draw-two-heatmap-with-the-same-cell-size
#    or http://stackoverflow.com/questions/20638294/geom-tile-and-facet-grid-facet-wrap-for-same-height-of-tiles/20639481#20639481
#    or https://github.com/baptiste/gridextra/wiki/arranging-ggplot#fix-the-panel-size

# try adding color coded sample labels, as in http://sebastianraschka.com/Articles/heatmaps_in_r.html

# try reducing LOC by using loops, as in https://joey711.github.io/phyloseq-demo/unifrac.html
#    or http://www.reed.edu/data-at-reed/resources/R/loops_with_ggplot2.html



# Boxplots for top genera
# Healthy by diet
# CSS by NEC


