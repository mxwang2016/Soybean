# set work directory
setwd('F:/Code')


# load packages
library(reshape2)
library(ggplot2)


# load function
data_processing <- function(relative_df, absolute_df, group_df) {
    relative_tmp_df <- relative_df
    relative_tmp_df <- relative_tmp_df[rowSums(relative_tmp_df) > 0,]
    relative_tmp_df$Sum <- apply(relative_tmp_df, 1, sum)
    relative_tmp_df <- relative_tmp_df[order(-relative_tmp_df$Sum),]
    relative_tmp_df <- relative_tmp_df[-ncol(relative_tmp_df)]
    
    if ('unidentified' %in% row.names(relative_tmp_df)) {
        unidentified_df <- relative_tmp_df['unidentified',]
        relative_tmp_df <- relative_tmp_df[row.names(relative_tmp_df) != 'unidentified',]
        relative_tmp_df <- rbind(relative_tmp_df, unidentified_df)
    }
    
    if (nrow(relative_tmp_df) <= 10) {
        relative_top_df <- relative_tmp_df
        absolute_top_df <- absolute_df[row.names(relative_top_df),]
        if ('unidentified' %in% row.names(relative_tmp_df)) {
            row.names(relative_top_df)[nrow(relative_top_df)] <- 'Others'
            row.names(absolute_top_df)[nrow(relative_top_df)] <- 'Others'
            color_manual <- color[1:nrow(relative_top_df)]
        } else {
            color_manual <- color[1:(nrow(relative_top_df) + 1)]
        }
        
    } else {
        relative_top_df <- relative_tmp_df[1:9,]
        relative_top_df[10,] <- apply(relative_tmp_df[10:nrow(relative_tmp_df),], 2, sum)
        row.names(relative_top_df)[10] <- 'Others'
        absolute_top_df <- absolute_df[row.names(relative_tmp_df)[1:9],]
        absolute_top_df[10,] <- apply(absolute_df[row.names(relative_tmp_df)[10:nrow(relative_tmp_df)],], 2, sum)
        row.names(absolute_top_df)[10] <- 'Others'
        color_manual <- color
    }
    
    relative_top_df$Taxonomy <- row.names(relative_top_df)
    relative_top_df_melt <- melt(relative_top_df, id.vars = 'Taxonomy', variable.name = 'Sample', value.name = 'Value')
    relative_fin_df <- merge(relative_top_df_melt, group_df, by.x = 'Sample', by.y = 'row.names', all.x = TRUE)
    relative_fin_df$Taxonomy <- factor(relative_fin_df$Taxonomy, levels = relative_top_df$Taxonomy)
    
    absolute_top_df$Taxonomy <- row.names(absolute_top_df)
    absolute_top_df_melt <- melt(absolute_top_df, id.vars = 'Taxonomy', variable.name = 'Sample', value.name = 'Value')
    absolute_fin_df <- merge(absolute_top_df_melt, group_df, by.x = 'Sample', by.y = 'row.names', all.x = TRUE)
    absolute_fin_df$Taxonomy <- factor(absolute_fin_df$Taxonomy, levels = absolute_top_df$Taxonomy)
    
    color_manual <- color_manual[1:nrow(relative_top_df)]
    top_color <- data.frame(Taxonomy = relative_top_df$Taxonomy, Color = color_manual)
    
    data_df <- list(relative_fin_df, absolute_fin_df, color_manual, top_color)
    return(data_df)
}



# import data
taxonomy = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
color = c('#e53525', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#2bffd5', '#f781bf', '#a65628', '#999999')
treatment <- c('Ctrl', '-N', '-P', '-K')
stage = c('D0', 'D1', 'D4', 'D7', 'D14', 'D28', 'D42', 'D60', 'D72')

rmp_taxonomy_data <- read.csv('rmp_taxonomy_data.csv', check.names = F)
qmp_taxonomy_data <- read.csv('qmp_taxonomy_data.csv', check.names = F)

metadata <- read.csv('metadata.csv', row.names = 1)
metadata$Treatment <- factor(metadata$Treatment, levels = treatment)
metadata$Stage <- factor(metadata$Stage, levels = stage)
metadata <- metadata[order(metadata$Stage, metadata$Treatment),]


# calculating data
relative_df <- aggregate(rmp_taxonomy_data[!(names(rmp_taxonomy_data) %in% taxonomy)], by = list(Taxonomy = rmp_taxonomy_data[,'Phylum']), sum)
row.names(relative_df) <- relative_df$Taxonomy
relative_df <- relative_df[-1]
relative_df <- data.frame(apply(relative_df, 2, function(x) {x / sum(x)}), check.names = F)

absolute_df <- aggregate(qmp_taxonomy_data[!(names(qmp_taxonomy_data) %in% taxonomy)], by = list(Taxonomy = qmp_taxonomy_data[,'Phylum']), sum)
row.names(absolute_df) <- absolute_df$Taxonomy
absolute_df<- absolute_df[-1]

data_df <- data_processing(relative_df, absolute_df, metadata)
relative_fin_df <- data_df[[1]]
relative_fin_df$Sample <- factor(relative_fin_df$Sample, levels = row.names(metadata))

absolute_fin_df <- data_df[[2]]
absolute_fin_df$Sample <- factor(absolute_fin_df$Sample, levels = row.names(metadata))

color_manual <- data_df[[3]]
top_color <- data_df[[4]]
names(top_color)[1] <- 'Phylum'

# 1 RMP
## 1.1 BS
tmp_df <- relative_fin_df[relative_fin_df$Compartment == 'BS',]

p <- ggplot(tmp_df, aes(x = Sample, y = Value, fill = Taxonomy)) + 
    geom_bar(stat = 'identity', position = position_fill(reverse = TRUE), width = 1) +
    labs(x = '', 
         y = 'Relative abundance') + 
    scale_y_continuous(labels = scales::percent) + 
    facet_grid( ~ Compartment, scales = 'free_x', switch = 'x') + 
    scale_fill_manual(values = color_manual) + 
    theme_bw() + 
    theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
          axis.title = element_text(size = 7, color = 'black'), 
          axis.line = element_blank(),
          axis.text = element_text(size = 6, color = 'black'), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          strip.text.x = element_text(size = 6, margin = margin(0.05, 0, 0.05, 0, 'cm')),
          strip.text.y = element_text(size = 6, margin = margin(0, 0.05, 0, 0.05, 'cm')),
          panel.spacing.x = unit(0.1, 'lines'),
          legend.title = element_text(size = 7, color = 'black'), 
          legend.text = element_text(size = 6,  color = 'black'), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          legend.position = 'none')

width = 2.2; height = 8
ggsave('BS_relative_barplot.pdf', p, width = width, height = height, units = 'cm')


# 1.2 R&RS
for (niche in c('R', 'RS')) {
    tmp_df <- relative_fin_df[relative_fin_df$Compartment == niche,]
    p <- ggplot(tmp_df, aes(x = Sample, y = Value, fill = Taxonomy)) + 
        geom_bar(stat = 'identity', position = position_fill(reverse = TRUE), width = 1) +
        labs(x = '', 
             y = 'Relative abundance') + 
        scale_y_continuous(labels = scales::percent) + 
        facet_grid( ~ Stage, scales = 'free_x', switch = 'x') + 
        scale_fill_manual(values = color_manual) + 
        theme_bw() + 
        theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
              plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
              axis.title = element_text(size = 7, color = 'black'), 
              axis.line = element_blank(),
              axis.text = element_text(size = 6, color = 'black'), 
              axis.text.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              strip.text.x = element_text(size = 6, margin = margin(0.05, 0, 0.05, 0, 'cm')),
              strip.text.y = element_text(size = 6, margin = margin(0, 0.05, 0, 0.05, 'cm')),
              panel.spacing.x = unit(0.1, 'lines'),
              legend.title = element_text(size = 7, color = 'black'), 
              legend.text = element_text(size = 6,  color = 'black'), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              legend.key.size = unit(0.25, 'cm'), 
              legend.position = 'none')
    
    width = 6.7; height = 8
    ggsave(paste0(niche, '_relative_barplot.pdf'), p, width = width, height = height, units = 'cm')
}


# 2 QMP
## 2.1 BS
tmp_df <- absolute_fin_df[absolute_fin_df$Compartment == 'BS',]

p <- ggplot(tmp_df, aes(x = Sample, y = Value, fill = Taxonomy)) + 
    geom_bar(stat = 'identity', position = position_stack(), width = 1) +
    labs(x = '', 
         y = 'Absolute abundance (log10copies g-1)') + 
    facet_grid( ~ Compartment, scales = 'free_x', switch = 'x') + 
    scale_fill_manual(values = color_manual) + 
    theme_bw() + 
    theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
          plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
          axis.title = element_text(size = 7, color = 'black'), 
          axis.line = element_blank(),
          axis.text = element_text(size = 6, color = 'black'), 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          strip.text.x = element_text(size = 6, margin = margin(0.05, 0, 0.05, 0, 'cm')),
          strip.text.y = element_text(size = 6, margin = margin(0, 0.05, 0, 0.05, 'cm')),
          panel.spacing.x = unit(0.1, 'lines'),
          legend.title = element_text(size = 7, color = 'black'), 
          legend.text = element_text(size = 6,  color = 'black'), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(), 
          legend.position = 'none')

width = 3.5; height = 6
ggsave('BS_absolute_barplot.pdf', p, width = width, height = height, units = 'cm')


# 2.2 R&RS
for (niche in c('R', 'RS')) {
    tmp_df <- absolute_fin_df[absolute_fin_df$Compartment == niche,]
    p <- ggplot(tmp_df, aes(x = Sample, y = Value, fill = Taxonomy)) + 
        geom_bar(stat = 'identity', position = position_stack(), width = 1) +
        labs(x = '', 
             y = 'Absolute abundance (log10copies g-1)') + 
        facet_grid( ~ Stage, scales = 'free_x', switch = 'x') + 
        scale_fill_manual(values = color_manual) + 
        theme_bw() + 
        theme(plot.title = element_text(size = 7, color = 'black', hjust = 0.5), 
              plot.subtitle = element_text(size = 6, color = 'black', hjust = 0.5), 
              axis.title = element_text(size = 7, color = 'black'), 
              axis.line = element_blank(),
              axis.text = element_text(size = 6, color = 'black'), 
              axis.text.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              strip.text.x = element_text(size = 6, margin = margin(0.05, 0, 0.05, 0, 'cm')),
              strip.text.y = element_text(size = 6, margin = margin(0, 0.05, 0, 0.05, 'cm')),
              panel.spacing.x = unit(0.1, 'lines'),
              legend.title = element_text(size = 7, color = 'black'), 
              legend.text = element_text(size = 6,  color = 'black'), 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              legend.key.size = unit(0.25, 'cm'), 
              legend.position = 'none')
    
    width = 14; height = 8
    ggsave(paste0(niche, '_absolute_barplot.pdf'), p, width = width, height = height, units = 'cm')
}

