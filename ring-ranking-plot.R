# Analysis of clonal trials - Ring ranking plot
# Luis A. Apiolaza
# School of Forestry
#

require(ggplot2)
setwd('~/Documents/Research/2012/clones10')


# Load Harewood ranking
load('~/Documents/Research/2012/harewood/newmoe-genetic-values.Rdata')

load('GeneticValues15Clones.Rdata')

clone.thr <- merge(clone.thr, moe.gv[, c('TENclone', 'class')], by.x = 'clone', by.y = 'TENclone')
clone.thr <- clone.thr[order(clone.thr$site, clone.thr$clone),]

plot.d <- data.frame(clone = clone.thr$clone[1:15],
                     gvgold = clone.thr$effect[1:15],
                     gvwait = clone.thr$effect[16:30],
                     class = clone.thr$class[1:15])

plot.d$gvmean <- with(plot.d, (gvgold + gvwait)/2)
plot.d$gvmean <- plot.d$gvmean + 8

# Color-blind-friendly palette
# http://wiki.stdout.org/rcookbook/Graphs/Colors%20(ggplot2)/
cbfPalette <- scale_colour_manual(values=c("#D55E00", "#000000", "#56B4E9"))

ring.plot <- ggplot(plot.d, aes(x = 1, y = gvmean, colour = class, label = clone)) 
ring.plot <- ring.plot + geom_text(size = 3) + 
             scale_y_reverse('First ring to reach 7 GPa threshold', angle = 90) +
             scale_x_continuous('') + cbfPalette +
             opts(axis.title.x = theme_text(size = 12),
             axis.text.y = theme_text(size = 10, colour = 'black'),
             axis.text.x = theme_text(colour = 'white'),
             legend.position = 'none')

ring.plot

pdf('ring-ranking-plot.pdf', width = 2, height = 7)
ring.plot
dev.off()