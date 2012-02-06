# Analysis of clonal trials
# Luis A. Apiolaza
# School of Forestry
#

# options
options(stringsAsFactors = FALSE)

# libraries
require(lme4)
require(lattice)
require(Hmisc)  # for upData, describe, fancy summary
require(plyr)   # for ddply
require(asreml) # for decent REML

# default working directory and options
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
# to change. Update 2012-02-07: The updated codes reflect Paul's
# notation
silv[silv$site == 'Golden Downs' & silv$rep == 1 & silv$clone ==166,]
silv$core[115:124] <- 'v'
silv$core[125:134] <- 'c'

# And now I'm dropping clone 272 in rep4 of Waitarere, which has only four rings
# as well as compartment 80, which has only four clones and drop unused levels
silv <- subset(silv, !(clone == 272 & rep == 4 & site == 'Waitarere Beach'))
silv <- subset(silv, !(comp == 80))
silv$clone <- silv$clone[drop = TRUE]

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
summary(ring ~ site + clone, data = threshold, na.rm = TRUE)
xyplot(jitter(ring) ~ clone | site*rep, group = core, data = threshold)



#### Using asreml-R for data analysis ####

# First model: univariate with common rep, clone, core and 
# residual variances across sites
moe.as1 <- asreml(ring ~ site, random = ~ rep %in% site + clone + core %in% clone, 
                  data = threshold)
summary(moe.as1)

#$loglik
#[1] -159.8502

#$varcomp
#gamma    component    std.error   z.ratio constraint
#rep:site!rep.var     1.983905e-01 1.982330e-01 1.228581e-01  1.613513   Positive
#clone!clone.var      2.421327e-01 2.419405e-01 1.133691e-01  2.134096   Positive
#clone:core!clone.var 1.011929e-07 1.011126e-07 9.047404e-09 11.175868   Boundary
#R!variance           1.000000e+00 9.992065e-01 8.940750e-02 11.175868   Positive

# Second model: equivalent to two univariate analyses 
# at site level, allowing for different rep, clone, core and residual
# variances.
moe.as2 <- asreml(ring ~ site, random = ~ diag(site):rep + diag(site):clone + 
                  diag(site):core %in% clone, rcov = ~ at(site):units, data = threshold)
summary(moe.as2)

#$loglik
#[1] -151.603

#$varcomp
#gamma    component  std.error   z.ratio
#site:rep!site.Golden Downs.var           2.371610e-02 2.371610e-02 0.03297339 0.7192496
#site:rep!site.Waitarere Beach.var        5.064329e-01 5.064329e-01 0.40641894 1.2460858
#site:clone!site.Golden Downs.var         2.647989e-01 2.647989e-01 0.13714140 1.9308459
#site:clone!site.Waitarere Beach.var      2.945326e-01 2.945326e-01 0.17046140 1.7278552
#site:clone:core!site.Golden Downs.var    3.266999e-02 3.266999e-02 0.06350567 0.5144421
#site:clone:core!site.Waitarere Beach.var 7.931322e-08 7.931322e-08         NA        NA
#site_Golden Downs!variance               6.431478e-01 6.431478e-01 0.08745054 7.3544178
#site_Waitarere Beach!variance            1.283638e+00 1.283638e+00 0.17370290 7.3898496
#
#constraint
#site:rep!site.Golden Downs.var             Positive
#site:rep!site.Waitarere Beach.var          Positive
#site:clone!site.Golden Downs.var           Positive
#site:clone!site.Waitarere Beach.var        Positive
#site:clone:core!site.Golden Downs.var      Positive
#site:clone:core!site.Waitarere Beach.var   Boundary
#site_Golden Downs!variance                 Positive
#site_Waitarere Beach!variance              Positive


# Third-model: full multivariate analysis
# There is no data to estimate rep or residual correlations, because we
# have two separate sites. However, we can estimate the correlation between
# clones
moe.as3 <- asreml(ring ~ site, random = ~ diag(site):rep + corgh(site):clone + 
                  diag(site):core %in% clone, rcov = ~ at(site):units, 
                  data = threshold)
summary(moe.as3)                  

#$loglik
#[1] -149.6627

#$varcomp
#gamma    component  std.error
#site:rep!site.Golden Downs.var                         2.178308e-02 2.178308e-02 0.03153581
#site:rep!site.Waitarere Beach.var                      4.360756e-01 4.360756e-01 0.35665131
#site:clone!site.Waitarere Beach:!site.Golden Downs.cor 7.471796e-01 7.471796e-01 0.27120679
#site:clone!site.Golden Downs                           2.697226e-01 2.697226e-01 0.13955495
#site:clone!site.Waitarere Beach                        2.726850e-01 2.726850e-01 0.16278103
#site:clone:core!site.Golden Downs.var                  3.578696e-02 3.578696e-02 0.06470654
#site:clone:core!site.Waitarere Beach.var               7.931322e-08 7.931322e-08         NA
#site_Golden Downs!variance                             6.417365e-01 6.417365e-01 0.08720998
#site_Waitarere Beach!variance                          1.300994e+00 1.300994e+00 0.17715740
#
#z.ratio    constraint
#site:rep!site.Golden Downs.var                         0.6907413      Positive
#site:rep!site.Waitarere Beach.var                      1.2226946      Positive
#site:clone!site.Waitarere Beach:!site.Golden Downs.cor 2.7550179 Unconstrained
#site:clone!site.Golden Downs                           1.9327338      Positive
#site:clone!site.Waitarere Beach                        1.6751642      Positive
#site:clone:core!site.Golden Downs.var                  0.5530656      Positive
#site:clone:core!site.Waitarere Beach.var                      NA      Boundary
#site_Golden Downs!variance                             7.3585216      Positive
#site_Waitarere Beach!variance                          7.3437159      Positive

# Total genetic correlation between sites for threshold: 0.75
# Broad sense heritabilities for thresholds:
# H2 = sigma_g^2/(sigma_g^2 + sigma_rep^2 + sigma_core^2 + sigma_e^2)
#  Golden Downs: 2.697226e-01/(2.697226e-01 + 2.178308e-02 + 3.578696e-02 + 6.417365e-01) = 0.28
#  Waitarere: 4.360756e-01/(4.360756e-01 + 4.360756e-01 + 7.931322e-08 + 1.300994e+00) = 0.20

# We treat clones as random so we can estimate their variances and regress 
# their values towards the mean by considering the inheritance of the trait.
# That is, if H2 is 1 there is no regression (equivalent to fixed effect) and
# if H2 is 0 the clonal values are totally shrink to zero


# Clonal genetic values. First extract all random effects
# and then keep only the ones referring to clones (these
# also keeps the interaction, that well drop).

# It is always tricky to sort out the results from
# a multivariate evaluation. I should write some functions
# for that.
clone.thr <- data.frame(coef(moe.as3, pattern = 'clone'))
clone.thr$site <- apply(data.frame(rownames(clone.thr)), 1, 
                        FUN = function(x) unlist(strsplit(x, ':'))[1])
clone.thr$clone <- apply(data.frame(rownames(clone.thr)), 1, 
                         FUN = function(x) unlist(strsplit(x, ':'))[2])

clone.thr$site <- substr(clone.thr$site, 6, 20)
clone.thr$clone <- substr(clone.thr$clone, 7, 20)
clone.thr <- clone.thr[1:30, ]


# Just in case, these results should add up to zero.
# Yes, they do.
sum(clone.thr$effect)
#[1] -2.592239e-12

xyplot(effect[1:15] ~ effect[16:30], data = clone.thr)
