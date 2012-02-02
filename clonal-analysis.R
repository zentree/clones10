# Analysis of clonal trials
# Luis A. Apiolaza
# School of Forestry
#

# options
options(stringsAsFactors = FALSE)

# libraries
require(lme4)
require(lattice)
require(Hmisc) #for upData, describe, fancy summary
require(plyr)  # for ddply

# deafult working directory and options
setwd('~/Documents/Research/2012/clones10')
critical <- 6 # Quality threshold in GPa

# Reading and sorting out data
silv <- read.csv('Master.Data.csv', skip = 1, row.names = 1, header= FALSE)
names(silv) <- c('site', 'comp', 'rep', 'clone', 'core', 
                 'age', 'rwidth', 'dens', 'year', 'mfa', 'moe')
                
head(silv)
silv$site <- factor(silv$site)
silv$comp <- factor(silv$comp)
silv$rep <- factor(silv$rep)
silv$clone <- factor(silv$clone)

# New age counting from the bark side
silv$new.age <- silv$year - 2000

silv <- upData(silv,
               labels = list(comp = 'compartment', rep = 'replicate', rwidth = 'ring width',
                             dens = 'density', mfa = 'microfibril angle', moe = 'modulus of elasticity'),
               units = list(rwidth = 'mm', dens = 'kg m-3', mfa = 'degrees', moe = 'GPa', new.age = 'years'))

describe(silv)

# One compartment in Golden Downs (63), Two in Waitarere (47, 80)
xtabs(~ comp + clone + site, data = silv)
# Clones 253, 259, 272 and 286 are only in Waitarere
xtabs(~ clone + site, data = silv)

# Now we need to sort out groups of records that constitute a
# uniquely identifiable core
# This table shows a range of rings going from 4 to 11, with most
# trees having 8 ~ 10. There is a problem with clone 166, rep 1 in
# Golden Downs, which has 20 values for core a and 20 for core x
xtabs(~ core + clone + rep + site, data = silv)

# In a totally arbitrary decision I am changing the third and fourth
# core to other types. We can change this later. First figure out rows
# to change
silv[silv$site == 'Golden Downs' & silv$rep == 1 & silv$clone ==166,]
silv$core[115:124] <- 'b'
silv$core[125:134] <- 'v'

# And now I'm dropping clone 272 in rep4 of Waitarere, which has only four rings
silv <- subset(silv, !(clone == 272 & rep == 4 & site == 'Waitarere Beach'))

# and know make core into a factor
silv$core <- factor(silv$core)
silv$tree.core <- with(silv, factor(paste(site, rep, clone, core, sep = ':')))

# Quick look at all trees per site and rep
xyplot(mfa ~ new.age|site*rep, group = tree.core, type = 'l', data = silv)
xyplot(moe ~ new.age|site*rep, group = tree.core, type = 'l', data = silv)
xyplot(dens ~ new.age|site*rep, group = tree.core, type = 'l', data = silv)

# Looking at clonal level
xyplot(mfa ~ new.age|clone, group = tree.core, type = 'l', 
       data = silv, subset = site == 'Golden Downs', main = 'Golden Downs')
xyplot(mfa ~ new.age|clone, group = tree.core, type = 'l', 
       data = silv, subset = site == 'Waitarere Beach', main = 'Waitarere Beach')

xyplot(dens ~ new.age|clone, group = tree.core, type = 'l', 
       data = silv, subset = site == 'Golden Downs', main = 'Golden Downs')
xyplot(dens ~ new.age|clone, group = tree.core, type = 'l', 
       data = silv, subset = site == 'Waitarere Beach', main = 'Waitarere Beach')

xyplot(moe ~ new.age|clone, group = tree.core, type = 'l', 
       data = silv, subset = site == 'Golden Downs', main = 'Golden Downs')
xyplot(moe ~ new.age|clone, group = tree.core, type = 'l', 
       data = silv, subset = site == 'Waitarere Beach', main = 'Waitarere Beach')


#### And here we start with the analyses ####


# Threshold function (how early can we reach 6 GPa?) 
moe.thres <- function(df) {
  ring <- min(df$new.age[df$moe >= critical])
}

threshold <- ddply(silv, .(site, rep, clone, core), moe.thres)
names(threshold)[5] <- 'ring'
threshold$ring <- with(threshold, ifelse(ring > 11, 12, ring))

# Description of thresholds. There doesn't seem to be a huge difference
# for thresholds (too rough a measurement?)
histogram(~ring | site, data = threshold)
xtabs(~ clone + site, data = threshold)
summary(ring ~ site + clone, data = threshold)

thres.lmer1 <- lmer(ring ~ site + (1|rep %in% site) + (1|clone), data = threshold)
summary(thres.lmer1)

thres.lmer2 <- lmer(ring ~ site + (1|rep %in% site) + (1|clone) + (1|clone:site), data = threshold)
summary(thres.lmer2)

anova(thres.lmer1, thres.lmer2)