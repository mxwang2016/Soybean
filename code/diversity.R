# set work directory
setwd('F:/Code')


# load packages
library(vegan)
library(Rmisc)
library(rstatix)
library(ggplot2)


# import data
treatment <- c('Ctrl', '-N', '-P', '-K')
stage = c('D0', 'D1', 'D4', 'D7', 'D14', 'D28', 'D42', 'D60', 'D72')

metadata <- read.csv('metadata.csv', row.names = 1)
metadata$Treatment <- factor(metadata$Treatment, levels = treatment)
metadata$Stage <- factor(metadata$Stage, levels = stage)

qmp_data <- read.csv('qmp_data.csv', row.names = 1, check.names = F)


# 1 alpha diversity
shannon <- data.frame(Shannon = diversity(as.data.frame(t(qmp_data)), index = 'shannon', base = exp(1)))
    

## 1.1 BS
group_bs <- metadata[metadata$Compartment == 'BS',]

shannon_df <- merge(group_bs, shannon, by = 'row.names', all.x = TRUE)
shannon_df_summary <- summarySE(shannon_df, measurevar = 'Shannon', groupvars = 'Treatment')

for (i in 1:4) {
    shannon_df_summary$max <- max(shannon_df[shannon_df$Treatment == shannon_df_summary$Treatment[i], 'Shannon'])
}

dunn_res <- dunn_test(shannon_df, Shannon ~ Treatment, p.adjust.method = 'fdr')
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

p <- ggplot(shannon_df, aes(x = Treatment, y = Shannon, color = Treatment)) + 
    geom_boxplot(width = 0.68, outlier.shape = NA) +
    geom_point(size = 0.5) + 
    geom_text(data = shannon_df_summary, aes(y = max * 1.05, label = Label), position = position_dodge(0.9), size = 2.5) +
    labs(x = '', y = 'Shannon index') + 
    scale_y_continuous(limits = c(6, 7.5), labels = c(6.0, 6.5, 7.0, 7.5)) +
    theme_bw() + 
    theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
          axis.title = element_text(size = 7, color = 'black'), 
          axis.text = element_text(size = 6, color = 'black'), 
          legend.title = element_text(size = 7, color = 'black'), 
          legend.text = element_text(size = 6,  color = 'black'), 
          legend.position = 'none')

width = 3.4; height = 3.7
ggsave('BS_alpha_diversity.pdf', p, width = width, height = height, units = 'cm')


## 1.2 R&RS
for (niche in c('R', 'RS')) {
    group_df <- metadata[metadata$Compartment == niche,]
    
    df <- merge(group_df, shannon, by = 'row.names', all.x = TRUE)
    df_summary <- summarySE(df, measurevar = 'Shannon', groupvars = c('Day', 'Treatment'))
    
    p <- ggplot(df_summary, aes(x = Day, y = Shannon, fill = Treatment, color = Treatment, group = Treatment)) +
        geom_line(size = 0.8) +
        labs(
            x = 'Days post germination (d)', 
            y = 'Shannon index',
        ) + 
        geom_ribbon(aes(ymin = Shannon - sd, ymax = Shannon + sd, fill = Treatment), alpha = 0.1, colour = NA) +
        scale_x_continuous(breaks = c(1, 4, 7, 14, 28, 42, 60, 72)) +
        theme_bw() + 
        theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
              plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
              axis.title = element_text(size = 7, color = 'black'), 
              axis.text = element_text(size = 6, color = 'black'), 
              axis.text.x = element_text(angle = 0), 
              legend.title = element_text(size = 7, color = 'black'), 
              legend.text = element_text(size = 6,  color = 'black'), 
              panel.grid.minor = element_blank(), 
              legend.position = 'none')
    
    width = 4.7; height = 3.7
    ggsave(paste0(niche, '_alpha_diversity.pdf'), p, width = width, height = height, units = 'cm')
}

# 2 beta diversity
bray_dist <- vegdist(as.data.frame(t(qmp_data)), method = 'bray')
bray_dist <- as.data.frame(as.matrix(bray_dist))

## 2.1 BS
tmp_bd <- bray_dist[row.names(group_bs), row.names(group_bs)]
tmp_dist <- as.dist(tmp_bd, diag = FALSE, upper = FALSE)

tmp_pcoa <- cmdscale(tmp_dist, k = 3, eig = TRUE)
tmp_eig <- round(tmp_pcoa$eig / sum(tmp_pcoa$eig) * 100, 2)

tmp_points <- as.data.frame(tmp_pcoa$points)
names(tmp_points) <- paste0('PCoA', 1:3)
points <- cbind(tmp_points, group_bs)

sig_res <- adonis2(tmp_dist ~ Treatment, data = group_bs)

p <- ggplot(points, aes(x = PCoA1, y = PCoA2, shape = Treatment, color = Treatment)) + 
    geom_point(size = 1, alpha = 0.75)+
    labs(
        subtitle = paste0('Treatment: R2 = ', round(sig_res$R2[1],3)), 
        x = paste0('PCoA1 (', tmp_eig[1], '%)'), 
        y = paste0('PCoA2 (', tmp_eig[2], '%)')
    ) + 
    scale_shape_manual(values = c(1, 15, 16, 17)) +
    theme_bw() + 
    stat_ellipse(level = 0.8) +
    theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
          axis.title = element_text(size = 7, color = 'black'), 
          axis.text = element_text(size = 6, color = 'black'), 
          legend.title = element_text(size = 7, color = 'black'), 
          legend.text = element_text(size = 6,  color = 'black'), 
          legend.position = 'none')

width = 3.5; height = 4.5
ggsave('BS_beta_diversity.pdf', p, width = width, height = height, units = 'cm')


## 2.2 R&RS
for (niche in c('RS', 'R')) {
    group_df <- metadata[metadata$Compartment == niche,]
    
    tmp_bd <- bray_dist[row.names(group_df), row.names(group_df)]
    tmp_dist <- as.dist(tmp_bd, diag = FALSE, upper = FALSE)
    
    tmp_pcoa <- cmdscale(tmp_dist, k = 3, eig = TRUE)
    tmp_eig <- round(tmp_pcoa$eig / sum(tmp_pcoa$eig) * 100, 2)
    
    tmp_points <- as.data.frame(tmp_pcoa$points)
    names(tmp_points) <- paste0('PCoA', 1:3)
    points <- cbind(tmp_points, group_df)
    
    sig_res <- adonis2(tmp_dist ~ Stage * Treatment, data = group_df)
    
    p <- ggplot(points, aes(x = PCoA1, y = PCoA2, shape = Treatment, color = Stage)) + 
        geom_point(size = 1, alpha = 0.75)+
        labs(

            subtitle = paste0('Stage: R2 = ', round(sig_res$R2[1],3), '   Treatment: R2 = ', round(sig_res$R2[2],3)), 
            x = paste0('PCoA1 (', tmp_eig[1], '%)'), 
            y = paste0('PCoA2 (', tmp_eig[2], '%)')
        ) + 
        scale_shape_manual(values = c(1, 15, 16, 17)) +
        scale_color_manual(values = c('#00512c', '#00a04b', '#84d775', '#c2ebba', '#ffc680', '#ff8c00', '#f34a00', '#e60700')) +
        theme_bw() + 
        theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
              plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
              axis.title = element_text(size = 7, color = 'black'), 
              axis.text = element_text(size = 6, color = 'black'), 
              legend.title = element_text(size = 7, color = 'black'), 
              legend.text = element_text(size = 6,  color = 'black'), 
              strip.text.x = element_text(size = 6, margin = margin(0.05, 0, 0.05, 0, 'cm')),
              panel.spacing.x = unit(0.1, 'lines'),
              legend.key.size = unit(0.25, 'cm'), 
              legend.position = 'none')
    
    width = 5; height = 4.5
    ggsave(paste0(niche, '_beta_diversity.pdf'), p, width = width, height = height, units = 'cm')
}

