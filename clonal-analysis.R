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

# The following line produces a 'clean' file for Michael Hayes
write.csv(silv, file = 'clonal-data-for-Hayes.csv', col.names = TRUE, quote = FALSE, row.names = FALSE)



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

#### Reading code translation 2 ~ 3 digits
code <- read.csv('code-translation.csv', header = TRUE)
code$FGclone <- as.character(code$FGclone)
code$Harewood <- as.character(code$Harewood)
code$GoldenDowns00 <- as.character(code$GoldenDowns00)
names(code)[5] <- 'TENclone'

#### And here we start with the analyses ####

# 2012-02-09 Added this quick phenotypic analysis for ring 2 MoE
r2moe <- aggregate(moe ~ site + clone, data = silv[silv$new.age == 2,], FUN = mean, na.rm = TRUE)
r2moe[order(r2moe$site,-r2moe$moe),]

# 2012-02-13 Added a quick mean tabulation for plot in presentation


# 2012-02-13 Added a plot for presentation
silv.incognito <- silv[silv$site == 'Waitarere Beach',]
silv.incognito <- merge(silv.incognito, code[,c('Harewood', 'TENclone')], 
                        by.x = 'clone', by.y = 'TENclone')
silv.incognito$cloneAO <- factor(silv.incognito$Harewood, labels = LETTERS[1:15])
waip <- ggplot(silv.incognito, aes(new.age, mfa, groups = tree.core))
waip <- waip + geom_line(color = "#56B4E9") + facet_wrap(~cloneAO) +
        scale_x_continuous('Age (years)') + scale_y_continuous('MFA (degrees)') +
        opts(axis.title.x = theme_text(size = 12),
             axis.text.y = theme_text(size = 10, colour = 'black'),
             axis.text.x = theme_text(colour = 'black'))

pdf('waitarere-trends.pdf', width = 8, height = 5.3)
waip
dev.off()


# Threshold function (how early can we reach the critical GPa value?) 
moe.thres <- function(df) {
  ring <- min(df$new.age[df$moe >= critical])
}

# Total tree diameter function (sum of ring widths * 2, ignores missing rings)
dia.max <- function(df) {
  dia <- sum(df$rwidth)*2
}


threshold <- ddply(silv, .(site, rep, clone, core), moe.thres)
names(threshold)[5] <- 'ring'
threshold$ring <- with(threshold, ifelse(ring > 11, 12, ring))

diameters <- ddply(silv, .(site, rep, clone, core), dia.max)
names(diameters)[5] <- 'dbh'

dp <- ggplot(diameters, aes(clone, dbh))

dp + geom_point() + facet_grid(~site)

# We do have missing first ring sometimes. We'll ignore it
# for a quick analysis
summary(silv$rwidth[silv$age == 1] ~ silv$clone[silv$age == 1], FUN = mean)

diam.mns <- ddply(diameters, .(site, rep, clone), function(df) mean(df$dbh))
names(diam.mns)[4] <- 'dbhmn'


# Description of thresholds. There is a site difference
histogram(~ring | site, data = threshold)
xtabs(~ clone + site, data = threshold)
summary(ring ~ site + clone, data = threshold, na.rm = TRUE)
xyplot(jitter(ring) ~ clone | site*rep, group = core, data = threshold)

# Description of diameters
histogram(~ dbhmn | site, data = diam.mns)

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



# Quick and dirty correlation between diameters and threshold

# Site 1
dbh.as1s1 <- asreml(dbhmn ~ 1, random = ~ rep + clone, 
                    data = diam.mns, subset = site == 'Golden Downs')
summary(dbh.as1s1)

# $varcomp
# gamma    component    std.error  z.ratio constraint
# rep!rep.var     1.011929e-07 5.549250e-05 1.074390e-05 5.165026   Boundary
# clone!clone.var 2.047021e-01 1.122552e+02 9.107121e+01 1.232609   Positive
# R!variance      1.000000e+00 5.483834e+02 1.061724e+02 5.165026   Positive
# 
#   > 1.122552e+02/(1.122552e+02 + 5.483834e+02)
# [1] 0.1699192

dbhgv.s1 <- data.frame(coef(dbh.as1s1, pattern = 'clone'))
dbhgv.s1$clone = substr(rownames(dbhgv.s1), 7,9)

moe.as1s1 <- asreml(ring ~ 1, random = ~ rep + clone/core, 
                    data = threshold, subset = site == 'Golden Downs')

summary(moe.as1s1)
# 
# $varcomp
# gamma    component    std.error  z.ratio
# rep!rep.var          9.337511e-03 6.842824e-03 2.322790e-02 0.294595
# clone!clone.var      3.055313e-01 2.239031e-01 1.164247e-01 1.923158
# clone:core!clone.var 1.011929e-07 7.415736e-08 9.302187e-09 7.972035
# R!variance           1.000000e+00 7.328317e-01 9.192530e-02 7.972035
# constraint
# rep!rep.var            Positive
# clone!clone.var        Positive
# clone:core!clone.var   Boundary
# R!variance             Positive
# 
# > 2.239031e-01/(2.239031e-01 + 6.842824e-03 + 7.328317e-01)
# [1] 0.2323664

moegv.s1 <- data.frame(coef(moe.as1s1, pattern = 'clone'))
moegv.s1$clone <- substr(rownames(moegv.s1), 7,9)
moegv.s1 <- moegv.s1[1:15,]

gvsite1 <- merge(dbhgv.s1, moegv.s1, by = 'clone')
names(gvsite1) <- c('clone', 'gvdbh', 'gvmoe')
qplot(gvdbh, gvmoe, data = gvsite1)
with(gvsite1, cor.test(gvdbh, gvmoe))

# Pearson's product-moment correlation
# 
# data:  gvdbh and gvmoe 
# t = 0.4166, df = 13, p-value = 0.6838
# alternative hypothesis: true correlation is not equal to 0 
# 95 percent confidence interval:
#  -0.4223158  0.5922201 
# sample estimates:
#       cor 
# 0.1147777 


# Site 2
dbh.as1s2 <- asreml(dbhmn ~ 1, random = ~ rep + clone, 
                    data = diam.mns, subset = site == 'Waitarere Beach')
summary(dbh.as1s2)

# $loglik
# [1] -240.0974
# 
# $varcomp
# gamma    component    std.error   z.ratio constraint
# rep!rep.var     1.014576e-01 1.501101e+02 2.139468e+02 0.7016235   Positive
# clone!clone.var 4.733260e-07 7.003026e-04 1.366772e-04 5.1237700   Boundary
# R!variance      1.000000e+00 1.479536e+03 2.887592e+02 5.1237700   Positive

# H2 ~ NA (clonal variance is zero)

moe.as1s2 <- asreml(ring ~ 1, random = ~ rep + clone/core, 
                    data = threshold, subset = site == 'Waitarere Beach')

summary(moe.as1s2)

# $loglik
# [1] -86.10185
# 
# $varcomp
# gamma    component    std.error  z.ratio
# rep!rep.var          5.453764e-01 5.856674e-01 4.556011e-01 1.285483
# clone!clone.var      4.134532e-01 4.439981e-01 2.180546e-01 2.036179
# clone:core!clone.var 1.011929e-07 1.086688e-07 1.470765e-08 7.388590
# R!variance           1.000000e+00 1.073877e+00 1.453427e-01 7.388590
# constraint
# rep!rep.var            Positive
# clone!clone.var        Positive
# clone:core!clone.var   Boundary
# R!variance             Positive

# > 4.439981e-01/(4.439981e-01 + 5.856674e-01 + 1.073877e+00)
# [1] 0.2110716 H2

dbhgv.s2 <- data.frame(coef(dbh.as1s2, pattern = 'clone'))
dbhgv.s2$clone = substr(rownames(dbhgv.s2), 7,9)

moegv.s2 <- data.frame(coef(moe.as1s2, pattern = 'clone'))
moegv.s2$clone <- substr(rownames(moegv.s2), 7,9)
moegv.s2 <- moegv.s2[1:15,]

gvsite2 <- merge(dbhgv.s2, moegv.s2, by = 'clone')
names(gvsite2) <- c('clone', 'gvdbh', 'gvmoe')
qplot(gvdbh, gvmoe, data = gvsite2)
with(gvsite2, cor.test(gvdbh, gvmoe))

# Pearson's product-moment correlation
# 
# data:  gvdbh and gvmoe 
# t = 1.4754, df = 13, p-value = 0.1639
# alternative hypothesis: true correlation is not equal to 0 
# 95 percent confidence interval:
#  -0.1656922  0.7462117 
# sample estimates:
#       cor 
# 0.3787152 




# Second model: equivalent to two univariate analyses 
# at site level, allowing for different rep, clone and residual
# variances.
dbh.as2 <- asreml(dbhmn ~ site, random = ~ diag(site):rep + diag(site):clone,
                  rcov = ~ at(site):units, data = diam.mns)
summary(dbh.as2)

# $loglik
# [1] -491.5424
# 
# $varcomp
# gamma    component std.error
# site:rep!site.Golden Downs.var      1.379432e-04 1.379432e-04        NA
# site:rep!site.Waitarere Beach.var   1.501146e+02 1.501146e+02 213.97691
# site:clone!site.Golden Downs.var    1.122552e+02 1.122552e+02  91.07125
# site:clone!site.Waitarere Beach.var 6.589371e-04 6.589371e-04        NA
# site_Golden Downs!variance          5.483834e+02 5.483834e+02 106.17245
# site_Waitarere Beach!variance       1.479533e+03 1.479533e+03 288.75586
# z.ratio constraint
# site:rep!site.Golden Downs.var             NA   Boundary
# site:rep!site.Waitarere Beach.var   0.7015458   Positive
# site:clone!site.Golden Downs.var    1.2326085   Positive
# site:clone!site.Waitarere Beach.var        NA   Boundary
# site_Golden Downs!variance          5.1650256   Positive

# Third-model: full multivariate analysis
# There is no data to estimate rep or residual correlations, because we
# have two separate sites. However, we can estimate the correlation between
# clones
dbh.as3 <- asreml(dbhmn ~ site, random = ~ diag(site):rep + corgh(site):clone,
                  rcov = ~ at(site):units, 
                  data = diam.mns)
summary(dbh.as3)                  

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
if(drop.comp.80 == TRUE) query.data <- moe.as3 else query.data <- moe.as4


clone.thr <- data.frame(coef(query.data, pattern = 'clone'))
clone.thr$site <- apply(data.frame(rownames(clone.thr)), 1, 
                        FUN = function(x) unlist(strsplit(x, ':'))[1])
clone.thr$clone <- apply(data.frame(rownames(clone.thr)), 1, 
                         FUN = function(x) unlist(strsplit(x, ':'))[2])

if(drop.comp.80 == TRUE) {
  clone.thr <- clone.thr[1:30, ] 
} else {
  clone.thr <- clone.thr[1:38, ]
}

clone.thr$site <- factor(substr(clone.thr$site, 6, 20))
clone.thr$clone <- factor(substr(clone.thr$clone, 7, 20))



str(clone.thr)


# Just in case, these results should add up to zero.
# Yes, they do.
round(sum(clone.thr$effect), 6)
#[1] 0


# Saving genetic values
save(clone.thr, file = 'GeneticValues15Clones.Rdata')

# Let's add an overall mean, so graphs are more
# meaningful to audience
if(drop.comp.80 == TRUE) {
  intercept <- coef(query.data)$fixed[3]
} else {
  intercept <- coef(query.data)$fixed[4]
}

clone.thr$adj.effect <- clone.thr$effect + intercept


# Plotting across sites using ggplot
# Dotplot
dp <- ggplot(clone.thr, aes(adj.effect, clone, colour = site)) + geom_point(size = 2) +
      scale_x_continuous('Genetic value for first ring reaching 7 GPa (lower is better)') +
      scale_y_discrete('Clone code') + theme_bw() +
      opts(axis.title.x = theme_text(vjust = 0.2))

pdf('dotplot-for-clones.pdf', width = 8, height = 5.3)
dp
dev.off()


# Scatterplot
if(drop.comp.80 == TRUE) {
  scat <- ggplot(clone.thr, 
                 aes(x = adj.effect[16:30], y = adj.effect[1:15], label = clone[1:15]))
} else {
  scat <- ggplot(clone.thr, 
                 aes(x = adj.effect[20:38], y = adj.effect[1:19], label = clone[1:19]))
}

scat <- scat + geom_point(size = 1) + geom_text(size = 3, vjust = 1.5, hjust = 1) + 
        scale_x_continuous('Genetic Value for first 7 GPa ring - Waitarere Beach') +
        scale_y_continuous('Genetic Value for first 7 GPa ring - Golden Downs') +
        theme_bw() +
        opts(axis.title.x = theme_text(vjust = 0.2))

pdf('scatterplot-for-clones.pdf', width = 8, height = 5.3)
scat
dev.off()

