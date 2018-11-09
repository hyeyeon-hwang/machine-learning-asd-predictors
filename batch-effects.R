# example from: http://jtleek.com/genstats/inst/doc/02_13_batch-effects.html

library(devtools)
library(Biobase)
library(sva)
library(bladderbatch)
library(snpStats)

data(bladderdata)

pheno = pData(bladderEset)
edata = exprs(bladderEset)

# adjusting for batch effects with a linear model
mod = model.matrix(~as.factor(cancer) + as.factor(batch), data = pheno)
fit = lm.fit(mod, t(edata))
hist(fit$coefficients[2,], col = 2, breaks = 100)

table(pheno$cancer, pheno$batch)

# Combat returns a "cleaned" data matrix without batch effects
# adjusting for batch effects with Combat: use this
batch = pheno$batch
modcombat = model.matrix(~1, data = pheno)
modcancer = model.matrix(~cancer, data = pheno)
combat_edata = ComBat(dat = edata, 
                      batch = batch, 
                      mod = modcombat, 
                      par.prior = TRUE, 
                      prior.plots = FALSE)

combat_fit = lm.fit(modcancer, t(combat_edata))
hist(combat_fit$coefficients[2,], col = 2, breaks = 100)




