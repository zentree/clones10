# Analysis of clonal trials - Ring ranking plot
# Luis A. Apiolaza
# School of Forestry
#

require(ggplot2)
setwd('~/Documents/Research/2012/clones10')

load('GeneticValues15Clones.Rdata')
plot.d <- data.frame(clone = clone.thr$clone[1:15],
                     gvgold = clone.thr$effect[1:15],
                     gvwait = clone.thr$effect[16:30])

plot.d$gvmean <- with(plot.d, (gvgold + gvwait)/2)
plot.d$gvmean <- plot.d$gvmean + 8

plot.d$class <- factor(cut(plot.d$gvmean, breaks = c(0, 7.80, 8.15, 9), 
                           labels = c('top', 'middle', 'bottom')))

# Color-blind-friendly palette
# http://wiki.stdout.org/rcookbook/Graphs/Colors%20(ggplot2)/
cbfPalette <- scale_colour_manual(values=c("#56B4E9", "#000000", "#D55E00"))

ring.plot <- ggplot(plot.d, aes(x = gvmean, y = 1, colour = class, label = clone)) 
ring.plot <- ring.plot + geom_text(size = 3, angle = 90) + 
             scale_x_continuous('First ring to reach 7 GPa threshold') +
             scale_y_continuous('') + cbfPalette +
             opts(axis.title.x = theme_text(size = 12),
             axis.text.x = theme_text(size = 10, colour = 'black'),
             axis.text.y = theme_text(colour = 'white'),
             legend.position = 'none')

pdf('ring-ranking-plot.pdf', width = 7, height = 2)
ring.plot
dev.off()