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
# Oh, the embarrassment! This should be done using something related
# to apply but I couldn't make it work with multiple argument functions
# (as finding the threshold). I have to try later

# Empty data-frame for results
size <- length(levels(silv$tree.core))
moe.thres <- data.frame(site = character(size), rep = character(size),
                        clone = character(size), tree.core = character(size),
                        threshold = numeric(size))
rec <- 0
for(co in levels(silv$tree.core)){
	rec <- rec + 1
	cdata <- subset(silv, tree.core == co)
	threshold <- min(cdata$new.age[cdata$moe >= critical])
	if(threshold > 11) threshold <- 12
	moe.thres[rec, 1:4] <- with(cdata, c(as.character(site[1]), as.character(rep[1]), 
	                            as.character(clone[1]), co))
	moe.thres[rec, 5] <- threshold
}

moe.thres$site <- factor(moe.thres$site)
moe.thres$rep <- factor(moe.thres$rep)
moe.thres$clone <- factor(moe.thres$clone)

# Description of thresholds. There doesn't seem to be a huge difference
# for thresholds (too rough a measurement?)
histogram(~threshold | site, data = moe.thres)
xtabs(~ clone + site, data = moe.thres)
summary(threshold ~ site + clone, data = moe.thres)

thres.lmer1 <- lmer(threshold ~ site + (1|rep %in% site) + (1|clone), data = moe.thres)
summary(thres.lmer1)

thres.lmer2 <- lmer(threshold ~ site + (1|rep %in% site) + (1|clone) + (1|clone:site), data = moe.thres)
summary(thres.lmer2)

anova(thres.lmer1, thres.lmer2)