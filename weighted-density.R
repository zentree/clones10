# Analysis of clonal trials
# Luis A. Apiolaza
# School of Forestry
#


# libraries
require(lattice)
require(Hmisc)    # for upData, describe, fancy summary
require(plyr)     # for ddply
require(asreml)   # for decent REML
require(ggplot2)


#### default working directory and local options ####
setwd('~/Documents/Research/2012/clones10')
critical <- 7          # Quality threshold in GPa (I'm now using RPBC's 7)
drop.comp.80 <- TRUE   # Should we drop compartment 80?
options(stringsAsFactors = FALSE)


#### Reading and sorting out data ####
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

# Compartment 80 in Waitarere contains four clones, which are
# not available anywhere else
xtabs(~  comp + rep + site, data = silv)

# Now we need to sort out groups of records that constitute a
# uniquely identifiable core
# This table shows a range of rings going from 4 to 11, with most
# trees having 8 ~ 10. There is a problem with clone 166, rep 1 in
# Golden Downs, which has 20 values for core a and 20 for core x
xtabs(~ core + clone + rep + site, data = silv)

# In a totally arbitrary decision I am changing the third and fourth
# core to other types. We can change this later. First figure out rows
# to change. Update 2012-02-07: The updated codes reflect Paul's
# notation
silv[silv$site == 'Golden Downs' & silv$rep == 1 & silv$clone ==166,]
silv$core[115:124] <- 'v'
silv$core[125:134] <- 'c'

# And now I'm dropping clone 272 in rep4 of Waitarere, which has only four rings
silv <- subset(silv, !(clone == 272 & rep == 4 & site == 'Waitarere Beach'))

# Should we drop compartment 80, which has only four clones and drop levels?
if(drop.comp.80 == TRUE) {
  silv <- subset(silv, !(comp == 80))
  silv$clone <- silv$clone[drop = TRUE]
}


# and know make core into a factor
silv$core <- factor(silv$core)
silv$tree.core <- with(silv, factor(paste(site, rep, clone, core, sep = ':')))

# Quick look at all trees per site and rep
xyplot(dens ~ age|site*rep, group = tree.core, type = 'l', data = silv)


# Weighted average function
xtabs(~  core + clone + rep + site, data = silv)

weighted.dens <- function(df){
  n <- length(df$rwidth)
  cum.radius <- cumsum(df$rwidth)
  prev.radius <- c(0, cum.radius[1:(n - 1)])
  weight <- pi * (cum.radius^2 - prev.radius^2)
  return(sum(weight*df$dens)/sum(weight))
}

weighdens <- ddply(silv, .(site, rep, clone, core), weighted.dens, .drop = TRUE)
weighdens <- ddply(silv, .(tree.core), weighted.dens)