# set work directory
setwd('F:/Code')


# load packages
library(rstatix)
library(ggplot2)
library(agricolae)


# import data
taxonomy = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
treatment <- c('Ctrl', '-N', '-P', '-K')
stage = c('D0', 'D1', 'D4', 'D7', 'D14', 'D28', 'D42', 'D60', 'D72')

qmp_taxonomy_data <- read.csv('qmp_taxonomy_data.csv', check.names = F)

metadata <- read.csv('metadata.csv', row.names = 1)
metadata$Treatment <- factor(metadata$Treatment, levels = treatment)
metadata$Stage <- factor(metadata$Stage, levels = stage)


# calculating data
absolute_df <- data.frame(Value = colSums(qmp_taxonomy_data[!names(qmp_taxonomy_data) %in% taxonomy]))


# BS
group_bs <- metadata[metadata$Compartment == 'BS',]
tmp_df <- merge(group_bs, absolute_df, by = 'row.names', all.x = T)

dunn_res <- dunn_test(tmp_df, Value ~ Treatment, p.adjust.method = 'fdr')
dunn_res_df <- data.frame(dunn_res[-1])

for (n in 1:nrow(dunn_res_df)) {
    if (dunn_res_df$p.adj[n] <= 0.001) 
        dunn_res_df$p.adj.signif[n] <- '***'
    else if (dunn_res_df$p.adj[n] <= 0.01) 
        dunn_res_df$p.adj.signif[n] <- '**'
    else if (dunn_res_df$p.adj[n] <= 0.05) 
        dunn_res_df$p.adj.signif[n] <- '*'
    else 
        dunn_res_df$p.adj.signif[n] <- ''
}

Label <- c('', dunn_res_df$p.adj.signif[1:3])

tmp_df$Value <- log10(tmp_df$Value)
p <- ggplot(tmp_df, aes(x = Treatment, y = Value, color = Treatment)) + 
    geom_boxplot(width = 0.68, outlier.shape = NA) +
    geom_point(size = 0.5) + 
    labs(
        x = '', 
        y = 'Absolute abundance (log10copies g-1)'
    ) + 
    facet_grid( ~ Compartment, scales = 'free_x', switch = 'x') + 
    theme_bw() + 
    theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
          axis.title = element_text(size = 7, color = 'black'), 
          axis.text = element_text(size = 6, color = 'black'), 
          legend.title = element_text(size = 7, color = 'black'), 
          legend.text = element_text(size = 6,  color = 'black'), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          legend.position = 'none')

width = 2.2; height = 6
ggsave('BS_absolute_boxplot.pdf', p, width = width, height = height, units = 'cm')


# R&RS
for (niche in c('R', 'RS')) {
    group_df <- metadata[metadata$Compartment == niche,]
    tmp_df <- merge(group_df, absolute_df, by = 'row.names', all.x = T)
    
    dunn_res <- dunn_test(tmp_df, Value ~ Stage, p.adjust.method = 'fdr')
    dunn_res_df <- data.frame(dunn_res[-1])
    
    mean_df <- aggregate(tmp_df['Value'], by = list(Group = tmp_df$Stage), FUN = mean)
    n <- nrow(mean_df)
    pvalue_df <- matrix(1, ncol = n, nrow = n)
    k <- 0
    for(i in 1:(n - 1)) { 
        for(z in (i + 1):n){ 
            k <- k + 1
            pvalue_df[i,z] <- dunn_res_df$p.adj[k]
            pvalue_df[z,i] <- dunn_res_df$p.adj[k]
        }
    }
    
    sig_df <- orderPvalue(mean_df$Group, mean_df$Value, 0.05, pvalue_df, console = TRUE)
    sig_df <- sig_df[stage[-1],]
    print(sig_df)
    
    tmp_df$Value <- log10(tmp_df$Value)
    p <- ggplot(tmp_df, aes(x = Treatment, y = Value, color = Treatment)) + 
        geom_boxplot(width = 0.68, outlier.shape = NA) +
        geom_point(size = 0.5) + 
        labs(
            x = '', 
            y = 'Absolute abundance (log10copies g-1)'
        ) + 
        facet_grid( ~ Stage, scales = 'free_x', switch = 'x') + 
        theme_bw() + 
        theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
              plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
              axis.title = element_text(size = 7, color = 'black'), 
              axis.text = element_text(size = 6, color = 'black'), 
              legend.title = element_text(size = 7, color = 'black'), 
              legend.text = element_text(size = 6,  color = 'black'), 
              axis.text.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              strip.text.x = element_text(size = 6, margin = margin(0.05, 0, 0.05, 0, 'cm')),
              panel.spacing.x = unit(0.1, 'lines'),
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              legend.key.size = unit(0.25, 'cm'),
              legend.position = 'none')
    
    width = 6.5; height = 6
    ggsave(paste0(niche, '_absolute_boxplot.pdf'), p, width = width, height = height, units = 'cm')
}
