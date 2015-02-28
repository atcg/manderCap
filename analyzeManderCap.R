setwd("~/src/manderCap/")
mappingStats <- read.csv("mappingStats.csv", stringsAsFactors=FALSE)
mappingStatsNoDSN <- mappingStats[-c(18,21),]

summary(lm(mappingStatsNoDSN$percAveOverTargetOver10 ~ mappingStatsNoDSN$individual_dna_in_capture * mappingStatsNoDSN$cot.multiplier))
bothVar_discountMapRateMod <- (lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture * mappingStatsNoDSN$cot.multiplier))


# Generate the nested models
cot1Mod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$cot.multiplier)
inputDNAmod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture)
individualMod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual)
cot1AndInputDNAmod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual_dna_in_capture)
cot1AndIndividualMod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual)
inputDNAandIndividualMod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$individual)
all3Mod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual)


cot1AndInputDNAmod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual_dna_in_capture + G(mappingStatsNoDSN$individual))



totalDNAmod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$total_dna_in_capture)
all4Mod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual + mappingStatsNoDSN$total_dna_in_capture)
all3v2Mod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$total_dna_in_capture) 



# Are CTS/BTS/F1 distinguishable?
plot(mappingStatsNoDSN$condition, mappingStatsNoDSN$CTSonly6iterDiscountMapRate, col=mappingStatsNoDSN$individual, pch=20, ylab="Discounted capture rate", xlab="Reaction condition number", main="CTS-only assembly after six ARC iterations--looking for differences between CTS/BTS/F1")


# Looks at predictive power of our variables
plot(mappingStatsNoDSN$individual_dna_in_capture, mappingStatsNoDSN$CTSonly6iterDiscountMapRate, col=mappingStatsNoDSN$individual, pch=20, cex=1, ylab="Discounted capture rate", xlab="Individual DNA in capture", main="CTS-only assembly after six ARC iterations")
plot(mappingStatsNoDSN$individual_dna_in_capture, mappingStatsNoDSN$CTSonly6iterRawMapRate, col=mappingStatsNoDSN$individual, pch=20, cex=1, ylab="Discounted capture rate", xlab="Individual DNA in capture", main="RAW MAP")
plot(mappingStatsNoDSN$individual_dna_in_capture, mappingStatsNoDSN$CTSonly6iterPercDupes, col=mappingStatsNoDSN$individual, pch=20, cex=1, ylab="Discounted capture rate", xlab="Individual DNA in capture", main="CTS-only assembly after six ARC iterations")


indDNAmod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture)
summary(indDNAmod)
abline(indDNAmod, col="red", lwd=2)
indDNAmodResiduals <- resid(indDNAmod)
indDNAmodInd<- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$individual)
summary(indDNAmodInd)


all3mod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier + mappingStatsNoDSN$individual)
summary(all3mod)
anova(all3mod)

bothMod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier)
plot(bothMod)
plot(fitted(bothMod), mappingStatsNoDSN$CTSonly6iterDiscountMapRate)
abline(coef=c(0,1), col="red", lwd=2)

# Plot of the residuals from the individual DNA model
plot(mappingStatsNoDSN$individual_dna_in_capture, indDNAmodResiduals, col=mappingStatsNoDSN$individual, pch=20, cex=1, ylab="Residual", xlab="Individual DNA in capture", main="CTS-only assembly after six ARC iterations--residuals from individual input DNA model")




# cot-1 Analyses
plot(mappingStatsNoDSN$cot.multiplier, mappingStatsNoDSN$CTSonly6iterDiscountMapRate, col=mappingStatsNoDSN$individual, pch=20, cex=1, ylab="Discounted capture rate", xlab="cot-1 multiplier", main="CTS-only assembly after six ARC iterations")
cot1mod <- lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$cot.multiplier)
summary(cot1mod)
abline(cot1mod, col="red", lwd=2)
cot1modResiduals <- resid(cot1mod)

# Plot of the residuals from the cot-1 model
plot(mappingStatsNoDSN$individual_dna_in_capture, cot1modResiduals, col=mappingStatsNoDSN$individual, pch=20, cex=1, ylab="Residual", xlab="Individual DNA in capture", main="CTS-only assembly after six ARC iterations--residuals from cot-1 model")



# Predict individual DNA residuals based on cot-1
plot(mappingStatsNoDSN$cot.multiplier, indDNAmodResiduals, col=mappingStatsNoDSN$individual, pch=20)
summary(lm(indDNAmodResiduals ~ mappingStatsNoDSN$cot.multiplier))
abline(lm(indDNAmodResiduals ~ mappingStatsNoDSN$cot.multiplier), col="red", lwd=2)

# Predict cot-1 residuals based on input DNA
plot(mappingStatsNoDSN$individual_dna_in_capture, cot1modResiduals, col=mappingStatsNoDSN$individual, pch=20)
summary(lm(cot1modResiduals ~ mappingStatsNoDSN$individual_dna_in_capture))
abline(lm(cot1modResiduals ~ mappingStatsNoDSN$individual_dna_in_capture), col="red", lwd=2)




summary(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture + mappingStatsNoDSN$cot.multiplier))
summary(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$cot.multiplier))
summary(lm(mappingStatsNoDSN$CTSonly6iterDiscountMapRate ~ mappingStatsNoDSN$individual_dna_in_capture))
