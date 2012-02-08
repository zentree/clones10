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
drop.comp.80 <- FALSE   # Should we drop compartment 80?
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


# Threshold function (how early can we reach the critical GPa value?) 
moe.thres <- function(df) {
  ring <- min(df$new.age[df$moe >= critical])
}

threshold <- ddply(silv, .(site, rep, clone, core), moe.thres)
names(threshold)[5] <- 'ring'
threshold$ring <- with(threshold, ifelse(ring > 11, 12, ring))

# Description of thresholds. There is a site difference
histogram(~ring | site, data = threshold)
xtabs(~ clone + site, data = threshold)
summary(ring ~ site + clone, data = threshold, na.rm = TRUE)
xyplot(jitter(ring) ~ clone | site*rep, group = core, data = threshold)



#### Using asreml-R for data analysis ####

# We treat clones as random so we can estimate their variances and regress 
# their values towards the mean by considering the inheritance of the trait.
# That is, if H2 is 1 there is no regression (equivalent to fixed effect) and
# if H2 is 0 the clonal values are totally shrink to zero


if(drop.comp.80 == TRUE) {

# First model: univariate with common rep, clone, core and 
# residual variances across sites
moe.as1 <- asreml(ring ~ site, random = ~ rep %in% site + clone + core %in% clone, 
                  data = threshold)
summary(moe.as1)

#$loglik
#[1] -150.6036

#$varcomp
#gamma    component    std.error   z.ratio constraint
#rep:site!rep.var     2.824173e-01 2.576539e-01 1.509078e-01  1.707359   Positive
#clone!clone.var      3.362507e-01 3.067670e-01 1.357569e-01  2.259678   Positive
#clone:core!clone.var 1.011929e-07 9.231991e-08 8.261133e-09 11.175212   Boundary
#R!variance           1.000000e+00 9.123162e-01 8.163748e-02 11.175212   Positive


# Second model: equivalent to two univariate analyses 
# at site level, allowing for different rep, clone, core and residual
# variances.
moe.as2 <- asreml(ring ~ site, random = ~ diag(site):rep + diag(site):clone + 
                  diag(site):core %in% clone, rcov = ~ at(site):units, data = threshold)
summary(moe.as2)

#$loglik
#[1] -148.5317

#$varcomp
#gamma    component  std.error   z.ratio constraint
#site:rep!site.Golden Downs.var           6.842800e-03 6.842800e-03 0.02322730 0.2946016   Positive
#site:rep!site.Waitarere Beach.var        5.856664e-01 5.856664e-01 0.45557702 1.2855487   Positive
#site:clone!site.Golden Downs.var         2.239028e-01 2.239028e-01 0.11642744 1.9231102   Positive
#site:clone!site.Waitarere Beach.var      4.439980e-01 4.439980e-01 0.21805546 2.0361701   Positive
#site:clone:core!site.Golden Downs.var    4.816808e-07 4.816808e-07         NA        NA   Boundary
#site:clone:core!site.Waitarere Beach.var 8.539546e-08 8.539546e-08         NA        NA   Boundary
#site_Golden Downs!variance               7.328308e-01 7.328308e-01 0.09192516 7.9720368   Positive
#site_Waitarere Beach!variance            1.073877e+00 1.073877e+00 0.14534288 7.3885782   Positive


# Third-model: full multivariate analysis
# There is no data to estimate rep or residual correlations, because we
# have two separate sites. However, we can estimate the correlation between
# clones
moe.as3 <- asreml(ring ~ site, random = ~ diag(site):rep + corgh(site):clone + 
                  diag(site):core %in% clone, rcov = ~ at(site):units, 
                  data = threshold)
summary(moe.as3)                  

# $loglik
# [1] -143.7327

# $varcomp
#                                                              gamma    component  std.error   z.ratio
# site:rep!site.Golden Downs.var                         6.224490e-03 6.224490e-03 0.02271525 0.2740225
# site:rep!site.Waitarere Beach.var                      5.329608e-01 5.329608e-01 0.41695869 1.2782102
# site:clone!site.Waitarere Beach:!site.Golden Downs.cor 9.589013e-01 9.589013e-01 0.14721162 6.5137612
# site:clone!site.Golden Downs                           2.201669e-01 2.201669e-01 0.11407334 1.9300466
# site:clone!site.Waitarere Beach                        4.278220e-01 4.278220e-01 0.21214974 2.0166039
# site:clone:core!site.Golden Downs.var                  4.490945e-07 4.490945e-07         NA        NA
# site:clone:core!site.Waitarere Beach.var               8.539546e-08 8.539546e-08         NA        NA
# site_Golden Downs!variance                             7.361216e-01 7.361216e-01 0.09249163 7.9587910
# site_Waitarere Beach!variance                          1.082456e+00 1.082456e+00 0.14706035 7.3606212

# constraint
# site:rep!site.Golden Downs.var                              Positive
# site:rep!site.Waitarere Beach.var                           Positive
# site:clone!site.Waitarere Beach:!site.Golden Downs.cor Unconstrained
# site:clone!site.Golden Downs                                Positive
# site:clone!site.Waitarere Beach                             Positive
# site:clone:core!site.Golden Downs.var                       Boundary
# site:clone:core!site.Waitarere Beach.var                    Boundary
# site_Golden Downs!variance                                  Positive
# site_Waitarere Beach!variance                               Positive

# Total genetic correlation between sites for threshold: 0.96
# Broad sense heritabilities for thresholds:
# H2 = sigma_g^2/(sigma_g^2 + sigma_rep^2 + sigma_e^2) Not using core, because is ~ 0
#  Golden Downs: 2.201669e-01/(2.201669e-01 + 6.224490e-03 + 7.361216e-01) = 0.23
#  Waitarere: 4.278220e-01/(4.278220e-01 + 5.329608e-01 + 1.082456e+00) = 0.21

}

if(drop.comp.80 == FALSE) {

# We have to get compartment back in the dataset
threshold$comp <- ifelse(threshold$site == 'Golden Downs', 63, 47)
threshold$comp <- ifelse(threshold$clone %in% c(253, 259, 272, 286), 80, threshold$comp)

# Full multivariate model
moe.as4 <- asreml(ring ~ site + at(site,2):comp, 
                  random = ~ at(site,1):rep + at(site,2):rep %in% comp + 
                             corgh(site):clone + diag(site):core %in% clone, 
                  rcov = ~ at(site):units, 
                  data = threshold)
summary(moe.as4)

# $loglik
# [1] -175.6418

# $varcomp
#                                                              gamma    component    std.error
# at(site, Golden Downs):rep!rep.var                     5.616531e-03 5.616531e-03 2.227641e-02
# rep:at(site, Waitarere Beach):comp!rep.var             3.711786e-05 3.711786e-05 3.459104e-05
# site:clone!site.Waitarere Beach:!site.Golden Downs.cor 9.770855e-01 9.770855e-01 1.337550e-01
# site:clone!site.Golden Downs                           2.253662e-01 2.253662e-01 1.129750e-01
# site:clone!site.Waitarere Beach                        4.753150e-01 4.753150e-01 2.088296e-01
# clone:site:core!site.Golden Downs.var                  4.402979e-07 4.402979e-07           NA
# clone:site:core!site.Waitarere Beach.var               8.285893e-08 8.285893e-08           NA
# site_Golden Downs!variance                             7.357379e-01 7.357379e-01 9.242360e-02
# site_Waitarere Beach!variance                          1.137019e+00 1.137019e+00 1.327555e-01
#                                                          z.ratio    constraint
# at(site, Golden Downs):rep!rep.var                     0.2521291      Positive
# rep:at(site, Waitarere Beach):comp!rep.var             1.0730484      Positive
# site:clone!site.Waitarere Beach:!site.Golden Downs.cor 7.3050420 Unconstrained
# site:clone!site.Golden Downs                           1.9948318      Positive
# site:clone!site.Waitarere Beach                        2.2760903      Positive
# clone:site:core!site.Golden Downs.var                         NA      Boundary
# clone:site:core!site.Waitarere Beach.var                      NA      Boundary
# site_Golden Downs!variance                             7.9604981      Positive
# site_Waitarere Beach!variance                          8.5647586      Positive

# Total genetic correlation between sites for threshold: 0.98
# Broad sense heritabilities for thresholds:
# H2 = sigma_g^2/(sigma_g^2 + sigma_rep^2 + sigma_e^2) Not using core, because is ~ 0
#  Golden Downs: 2.253662e-01/(2.253662e-01 + 5.616531e-03 + 7.357379e-01) = 0.23
#  Waitarere: 4.753150e-01/(4.753150e-01 + 3.711786e-05 + 1.137019e+00) = 0.29


}




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

clone.thr$site <- factor(substr(clone.thr$site, 6, 20))
clone.thr$clone <- factor(substr(clone.thr$clone, 7, 20))

if(drop.comp.80 == TRUE) clone.thr <- clone.thr[1:30, ]

str(clone.thr)


# Just in case, these results should add up to zero.
# Yes, they do.
sum(clone.thr$effect)
#[1] -2.592239e-12

# Saving genetic values
save(clone.thr, file = 'GeneticValues.Rdata')




# Plotting across sites
dotplot(clone ~ effect, groups = site, data = clone.thr, 
        xlab = 'Genetic value for first 7 GPa ring (lower is better)', 
        auto.key = TRUE)

xyplot(effect[1:15] ~ effect[16:30], data = clone.thr,
       xlab = 'GV for first 7 GPa ring - Waitarere Beach',
       ylab = 'GV for first 7 GPa ring - Golden Downs')


# Testing ggplot2 for some of the graphs
scat <- ggplot(clone.thr, aes(x = effect[16:30], y = effect[1:15], label = clone[1:15]))
scat + geom_point() + geom_text(vjust = 1.5, hjust = 1) + 
  scale_x_continuous('Genetic Value for first 7 GPa ring - Waitarere Beach') +
  scale_y_continuous('Genetic Value for first 7 GPa ring - Golden Downs')

