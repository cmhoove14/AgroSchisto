#This work is licensed under a Creative Commons Attribution-NonCommercial 4.0 International License#########
#<http://creativecommons.org/licenses/by-nc/4.0/> by Christopher Hoover, Arathi Arakala, Manoj Gambhir 
#and Justin Remais. This work was supported in part by the National Institutes of Health/National Science 
#Foundation Ecology of Infectious Disease program funded by the Fogarty International Center 
#(grant R01TW010286), the National Institute of Allergy and Infectious Diseases (grant K01AI091864), 
#and the National Science Foundation Water Sustainability and Climate Program (grant 1360330).

#Per the terms of this license, if you are making derivative use of this work, you must identify that 
#your work is a derivative work, give credit to the original work, provide a link to the license, 
#and indicate changes that were made.###############

source("Agrochemical_Review/Response_Fxs/Baxter_Rohr2011_reanalysis_atrazine_snail_carrying_capacity_fit.R")

#Atrazine effect on phi_Nq (snail carrying capacity) from Baxter et al 2011 data with Rohr et al analysis####
  bax.phin.df = data.frame(atra = c(0:200),
                           logatra = log(c(0:200)+1),
                           Prediction = 0,
                           st.err = 0,
                           Prediction.no30 = 0,
                           st.err.no30 = 0)
  
  bax.phin.df[,3:4] <- predict(bax.mod, newdata = bax.phin.df, 
                               type = 'response', se.fit = T)[1:2]
  bax.phin.df[,5:6] <- predict(bax.mod.no30, newdata = bax.phin.df, 
                               type = 'response', se.fit = T)[1:2]

png("Agrochemical_Review/Response_Fxs/Plots/Baxter_Rohr_2011/Baxter_rohr_atrazine_growth_rate_fit.png")  
plot(atra.df$atra, atra.df$growthrate, pch = 16, ylim = c(0,0.3),
     xlab = 'atrazine (ppb)', ylab = 'growth rate', main = "model fits to growth rate data")
  segments(x0 = atra.df$atra, y0 = atra.df$growthrate + atra.df$st.err,
           x1 = atra.df$atra, y1 = atra.df$growthrate - atra.df$st.err)
  
  lines(bax.phin.df$atra, bax.phin.df$Prediction, lty=2)
    lines(bax.phin.df$atra, bax.phin.df$Prediction + 1.96*bax.phin.df$st.err, lty=3)
    lines(bax.phin.df$atra, bax.phin.df$Prediction - 1.96*bax.phin.df$st.err, lty=3)  
  lines(bax.phin.df$atra, bax.phin.df$Prediction.no30, lty=2, col=2)
    lines(bax.phin.df$atra, bax.phin.df$Prediction.no30 + 1.96*bax.phin.df$st.err.no30, lty=3, col=2)
    lines(bax.phin.df$atra, bax.phin.df$Prediction.no30 - 1.96*bax.phin.df$st.err.no30, lty=3, col=2)  
    legend('bottomright', legend = c('30ppb included', '30ppb excluded'), lty = rep(2,2), col = c(1,2), cex = 0.7, bty='n')
    
dev.off()    

png("Agrochemical_Review/Response_Fxs/Plots/Baxter_Rohr_2011/Baxter_rohr_atrazine_carrying_capacity_sim.png")  
plot(atra.df$atra, atra.df$growthrate / atra.df$growthrate[1], 
     pch = 16, ylim = c(0.8, 3.0),
     xlab = 'atrazine (ppb)', ylab = 'growth rate',
     main = 'atrazine:peak growth rate (~carrying capacity)')
#Add error bars
  segments(x0 = atra.df$atra, 
           y0 = (atra.df$growthrate + atra.df$st.err) / atra.df$growthrate[1],
           x1 = atra.df$atra, 
           y1 = (atra.df$growthrate - atra.df$st.err) / atra.df$growthrate[1])
#add model predictions  
    points(c(0:200), sapply(c(0:200), phi_Nq_atr_baxrohr), pch = 17, col=4, cex=0.6)
    points(c(0:200), sapply(c(0:200), phi_Nq_atr_baxrohr.no30), pch = 17, col=2, cex=0.6)

  legend('topleft', legend = c('30ppb included', '30ppb excluded'), pch = 17, col = c(4,2), cex = 0.7, bty='n')
  
  dev.off()    