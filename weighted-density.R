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

write.csv(weighdens, '~/Dropbox/research/2012/forgen-density-per-core.csv', 
          quote = FALSE, row.names = FALSE)

names(weighdens)[5] <- 'wdens'


# Simple tree means
tree.den.means <- ddply(weighdens, .(site, rep, clone), function(df) mean(df$wdens), .drop = TRUE)
tree.den.means
names(tree.den.means)[4] <- 'wdens'




# ASReml analysis
# First model: univariate with common rep, clone, core and 
# residual variances across sites
wden.as1 <- asreml(wdens ~ site, random = ~ rep %in% site + clone + core %in% clone, 
                   data = weighdens)
summary(wden.as1)
# $call
# asreml(fixed = wdens ~ site, random = ~rep %in% site + clone + 
#   core %in% clone, data = weighdens)
# 
# $loglik
# [1] -971.7942
# 
# $nedf
# [1] 272
# 
# $sigma
# [1] 19.18688
# 
# $varcomp
# gamma   component std.error     z.ratio constraint
# rep:site!rep.var     0.0985530461  36.2809612  25.34827  1.43129958   Positive
# clone!clone.var      1.2707362576 467.8042406 185.02157  2.52837671   Positive
# clone:core!clone.var 0.0005332254   0.1962996  19.10362  0.01027552   Positive
# R!variance           1.0000000000 368.1363759  34.89964 10.54843051   Positive

# Second model: equivalent to two univariate analyses 
# at site level, allowing for different rep, clone, core and residual
# variances.
wden.as2 <- asreml(wdens ~ site, random = ~ diag(site):rep + diag(site):clone + 
                   diag(site):core %in% clone, rcov = ~ at(site):units, 
                   data = weighdens)
summary(wden.as2)
# $call
# asreml(fixed = wdens ~ site, random = ~diag(site):rep + diag(site):clone + 
#   diag(site):core %in% clone, rcov = ~at(site):units, data = weighdens)
# 
# $loglik
# [1] -976.6759
# 
# $nedf
# [1] 272
# 
# $sigma
# [1] 1
# 
# $varcomp
# gamma    component std.error   z.ratio constraint
# site:rep!site.Golden Downs.var           3.977233e+01 3.977233e+01  39.83252 0.9984890   Positive
# site:rep!site.Waitarere Beach.var        4.971817e+01 4.971817e+01  43.42615 1.1448903   Positive
# site:clone!site.Golden Downs.var         5.646406e+02 5.646406e+02 236.75349 2.3849304   Positive
# site:clone!site.Waitarere Beach.var      3.732754e+02 3.732754e+02 152.53599 2.4471299   Positive
# site:clone:core!site.Golden Downs.var    2.122704e+01 2.122704e+01  55.12762 0.3850528   Positive
# site:clone:core!site.Waitarere Beach.var 4.350783e-05 4.350783e-05        NA        NA   Boundary
# site_Golden Downs!variance               4.566106e+02 4.566106e+02  64.67653 7.0599119   Positive
# site_Waitarere Beach!variance            2.392417e+02 2.392417e+02  32.42754 7.3777319   Positive


# Third-model: full multivariate analysis
# There is no data to estimate rep or residual correlations, because we
# have two separate sites. However, we can estimate the correlation between
# clones
wden.as3 <- asreml(wdens ~ site, random = ~ diag(site):rep + corgh(site):clone + 
                  diag(site):core %in% clone, rcov = ~ at(site):units, 
                  data = weighdens)
summary(wden.as3)                  
# $call
# asreml(fixed = wdens ~ site, random = ~diag(site):rep + corgh(site):clone + 
#   diag(site):core %in% clone, rcov = ~at(site):units, data = weighdens)
# 
# $loglik
# [1] -962.7753
# 
# $nedf
# [1] 272
# 
# $sigma
# [1] 1
# 
# $varcomp
# gamma    component std.error   z.ratio constraint
# site:rep!site.Golden Downs.var                         3.106230e+01 3.106230e+01  33.68226 0.9222154   Positive
# site:rep!site.Waitarere Beach.var                      4.125628e+01 4.125628e+01  36.72439 1.1234026   Positive
# site:clone!site.Waitarere Beach:!site.Golden Downs.cor 9.990000e-01 9.990000e-01        NA        NA   Boundary
# site:clone!site.Golden Downs                           5.820398e+02 5.820398e+02 239.27583 2.4325058   Positive
# site:clone!site.Waitarere Beach                        3.609671e+02 3.609671e+02 147.80540 2.4421777   Positive
# site:clone:core!site.Golden Downs.var                  5.252370e-05 5.252370e-05        NA        NA   Boundary
# site:clone:core!site.Waitarere Beach.var               4.350783e-05 4.350783e-05        NA        NA   Boundary
# site_Golden Downs!variance                             4.678842e+02 4.678842e+02  57.83802 8.0895607   Positive
# site_Waitarere Beach!variance                          2.398936e+02 2.398936e+02  32.00942 7.4944703   Positive

# H2 Golden Downs: > 582/(582 + 31 + 468)
# [1] 0.5383904
# 
# H2 Waitarere Beach: > 361/(361 + 41 + 240)
# [1] 0.5623053
