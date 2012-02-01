# Analysis of clonal trials
# Luis A. Apiolaza
# School of Forestry
#

require(lattice)

setwd('~/Documents/Research/2012/clones10')


# Reading and sorting out data
silv <- read.csv('Master.Data.csv', skip = 1, row.names = 1, header= FALSE)
names(silv) <- c('site', 'comp', 'rep', 'clone', 'core', 
                 'age', 'rwidth', 'bdens', 'year', 'mfa', 'moe')
head(silv)
silv$comp <- factor(silv$comp)
silv$rep <- factor(silv$rep)
silv$clone <- factor(silv$clone)

# One compartment in Golden Downs (63), Two in Waitarere (47, 80)
xtabs(~ comp + clone + site, data = silv)
# Clones 253, 259, 272 and 286 are only in Waitarere
xtabs(~ clone + site, data = silv)

# Now we need to sort out groups of records that constitute a
# uniquely identifiable core
xtabs(~ core + clone + rep + site, data = silv)

# This table shows a range of rings going from 4 to 11, with most
# trees having 8 ~ 10. There is a problem with clone 166, rep 1 in
# Golden Downs, which has 20 values for core a and 20 for core x
silv$tree <- with(silv, site:rep:clone:core)

# Quick look at all trees per site and rep
xyplot(mfa ~ age|site*rep, group = tree, type = 'l', data = silv)
xyplot(moe ~ age|site*rep, group = tree, type = 'l', data = silv)
xyplot(bdens ~ age|site*rep, group = tree, type = 'l', data = silv)

# Looking at clonal level
xyplot(mfa ~ age|clone, group = tree, type = 'l', 
       data = silv, subset = site == 'Golden Downs', main = 'Golden Downs')
xyplot(mfa ~ age|clone, group = tree, type = 'l', 
       data = silv, subset = site == 'Waitarere Beach', main = 'Waitarere Beach')

xyplot(bdens ~ age|clone, group = tree, type = 'l', 
       data = silv, subset = site == 'Golden Downs', main = 'Golden Downs')
xyplot(bdens ~ age|clone, group = tree, type = 'l', 
       data = silv, subset = site == 'Waitarere Beach', main = 'Waitarere Beach')

xyplot(moe ~ age|clone, group = tree, type = 'l', 
       data = silv, subset = site == 'Golden Downs', main = 'Golden Downs')
xyplot(moe ~ age|clone, group = tree, type = 'l', 
       data = silv, subset = site == 'Waitarere Beach', main = 'Waitarere Beach')


