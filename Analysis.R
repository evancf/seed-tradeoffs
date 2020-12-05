# Packages ---------------------------------------------------------------------

library("lme4")
library("phia")


# Read in data -----------------------------------------------------------------

dat.ps.traps <- read.csv("dat.ps.traps.csv")
dat.pr.traps <- read.csv("dat.pr.traps.csv")

dat.ps.nursery <- read.csv("dat.ps.nursery.csv")
dat.pr.nursery <- read.csv("dat.pr.nursery.csv")



# Seed density vs seed size ----------------------------------------------------

ps.counts.mod <- lmer(median.mass.g*1000 ~ seed.count.log10 + portion.handled + (1|trap.id),
                      data = dat.ps.traps)
ps.counts.mod
confint(ps.counts.mod, method = "Wald")


pr.counts.mod <- lmer(median.mass.g*1000 ~ seed.count.log10 + portion.handled + (1|trap.id),
                      data = dat.pr.traps)
pr.counts.mod
confint(pr.counts.mod, method = "Wald")



# Seed size vs density dependent germination -----------------------------------

dat.ps.nursery$treat.sh <- relevel(dat.ps.nursery$treat.sh, ref="pot")
ps.dd.mod <- glmer(tot.germ ~ scale(pred.mass) * treat.sh + (1|tree.id), data=dat.ps.nursery, family=binomial)
ps.dd.mod
confint(ps.dd.mod, method = "Wald")

# Non-zero slope in potting soil?
testInteractions(ps.dd.mod, custom=list(treat.sh=c(1,0,0)), slope="scale(pred.mass)")
# Non-zero slope in far soil?
testInteractions(ps.dd.mod, custom=list(treat.sh=c(0,1,0)), slope="scale(pred.mass)")
# Non-zero slope in near soil?
testInteractions(ps.dd.mod, custom=list(treat.sh=c(0,0,1)), slope="scale(pred.mass)")


dat.pr.nursery$treat.sh <- relevel(dat.pr.nursery$treat.sh, ref="pot")
pr.tot.germ.y <- cbind(dat.pr.nursery$tot.germ,4-dat.pr.nursery$tot.germ)
pr.dd.mod <- glmer(pr.tot.germ.y ~ scale(pred.mass) * treat.sh + (1|tree.id), data=dat.pr.nursery, family=binomial)
pr.dd.mod
confint(pr.dd.mod, method = "Wald")

# Non-zero slope in potting soil?
testInteractions(pr.dd.mod, custom=list(treat.sh=c(1,0,0)), slope="scale(pred.mass)")
# Non-zero slope in far soil?
testInteractions(pr.dd.mod, custom=list(treat.sh=c(0,1,0)), slope="scale(pred.mass)")
# Non-zero slope in near soil?
testInteractions(pr.dd.mod, custom=list(treat.sh=c(0,0,1)), slope="scale(pred.mass)")



# Seedling biomass -------------------------------------------------------------

ps.biomass.mod <- glmer(seedling.biomass ~ scale(pred.mass) + treat.sh + (1|tree.id), data=dat.ps.nursery,family=Gamma(link=log))
ps.biomass.mod
confint(ps.biomass.mod, method = "Wald")

pr.biomass.mod <- glmer(seedling.biomass ~ scale(pred.mass) + treat.sh + (1|tree.id), data=dat.pr.nursery,family=Gamma(link=log))
pr.biomass.mod
confint(pr.biomass.mod, method = "Wald")

