##Load Packages
library(dplyr)
library(tidyr)
library(mgcv)
library(glmmTMB)
library(DHARMa)
library(grid)
library(gridExtra)
library(emmeans)
library(ggplot2)
library(viridis)

##Data Handling and model assumptions
#Load Copepod Survival Data CSV
cdata <- read.csv("LTEXP.csv")

#Create new proportions variables
#Create df with initial population values
initial_count <- cdata %>%
  filter(Experiment.Day == 0) %>%
  mutate(
    Initial_TF = Nb.of.F + Nb.of.PF,
    Initial_TC = Total.Copepods,
    Initial_M = Nb.of.M
  ) %>%
  select(Treatment, Replicate, Initial_TC, Initial_TF, Initial_M)

#Join initials to cdata
#Compute new variables using it.
cdata <- cdata %>%
  left_join(initial_count, by = c("Treatment", "Replicate")) %>%
  mutate(
    Current_TF = Nb.of.F + Nb.of.PF,
    Dead_TC = Initial_TC - Total.Copepods,
    Dead_TF = Initial_TF - Current_TF,
    Dead_M = Initial_M - Nb.of.M,
    Prop_Alive_TC = Total.Copepods/Initial_TC,
    Prop_Alive_TF = Current_TF/Initial_TF,
    Prop_Alive_M = Nb.of.M/Initial_M,
    Prop_PF = ifelse(Current_TF > 0,
                     Nb.of.PF/Current_TF, 0)
  )

#Check Data
summary(cdata)
head(cdata)
str(cdata)

#Convert Treatment and Replicate to Factor
cdata$Treatment <- as.factor(cdata$Treatment)
cdata$Replicate <- as.factor(cdata$Replicate)

#Check for Outliers
par(mfrow = c(2,2))

plot(Prop_Alive_M ~ Treatment, data = cdata)
plot(Prop_Alive_TC ~ Treatment, data = cdata)
plot(Prop_Alive_TF ~ Treatment, data = cdata)
plot(Prop_PF ~ Treatment, data = cdata)

boxplot(cdata$Prop_Alive_M,main = "Proportion of Alive Males")
boxplot(cdata$Prop_Alive_TC,main = "Proprtion of Alive Copepods")
boxplot(cdata$Prop_Alive_TF,main = "Proportion of Alive Females")
boxplot(cdata$Prop_PF,main = "Proportion of Pregant Females")
#Outliers visible but represent ecological variability so were not removed.

#Check Distribution, Zero Inflation
hist(cdata$Prop_Alive_TC, main = "Distribution of Prop_Alive_TC")
hist(cdata$Prop_Alive_TF, main = "Distribution of Prop_Alive_TF")
hist(cdata$Prop_Alive_M, main = "distribution of Prop_Alive_M")
hist(cdata$Prop_PF, main = "distribution of Prop_PF")
# Right-Skew for survival data, even distribution Prop_PF

#Colinearity among predictors
boxplot(Experiment.Day ~ Treatment, data = cdata)
#no colinearity

###Descriptive Stats:
##Full table of stats
Descriptive_Table <- cdata %>%
  group_by(Treatment) %>%
  summarise(
    n = n(),
    
    #Raw counts
    Total_Copepods_Mean = signif(mean(Total.Copepods), 3),
    Total_Copepods_SD = signif(sd(Total.Copepods), 3),
    
    Nb_of_Females_Mean = signif(mean(Current_TF), 3),
    Nb_of_Females_SD = signif(sd(Current_TF), 3),
    
    Nb_of_Males_Mean = signif(mean(Nb.of.M), 3),
    Nb_of_Males_SD = signif(sd(Nb.of.M), 3),
    
    Nb_of_PF_Mean = signif(mean(Nb.of.PF), 3),
    Nb_of_PF_SD = signif(sd(Nb.of.PF), 3),
    
    #Proportions
    Prop_Alive_TC_Mean = signif(mean(Prop_Alive_TC), 3),
    Prop_Alive_TC_SD = signif(sd(Prop_Alive_TC), 3),
    
    Prop_Alive_TF_Mean = signif(mean(Prop_Alive_TF), 3),
    Prop_Alive_TF_SD = signif(sd(Prop_Alive_TF), 3),
    
    Prop_Alive_M_Mean = signif(mean(Prop_Alive_M), 3),
    Prop_Alive_M_SD = signif(sd(Prop_Alive_M), 3),
    
    Prop_PF_Mean = signif(mean(Prop_PF), 3),
    Prop_PF_SD = signif(sd(Prop_PF), 3),
    
    .groups = 'drop'
  )
print(Descriptive_Table)

write.csv(Descriptive_Table, "descriptive_stats_summary.csv", row.names = FALSE)




##Models and Model Validation
##Proportion data - Binomial GAM Models
#Prop_Alive_TC Model
bgam_TC <- gam(cbind(Total.Copepods, Dead_TC) ~ Treatment + s(Experiment.Day, by = Treatment) + s(Replicate, bs = "re"),
               family = binomial(link = "logit"), data = cdata)

#Prop_Alive_TF Model
bgam_TF <- gam(cbind(Current_TF, Dead_TF) ~ Treatment + s(Experiment.Day, by = Treatment) + s(Replicate, bs = "re"),
               family = binomial(link = "logit"), data = cdata)

#Prop_Alive_M Model
bgam_M <- gam(cbind(Nb.of.M, Dead_M) ~ Treatment + s(Experiment.Day, by = Treatment) + s(Replicate, bs = "re"),
              family = binomial(link = "logit"), data = cdata)

#Prop_PF Model
bgam_Prop_PF <- gam(cbind(Nb.of.PF, Current_TF - Nb.of.PF) ~ Treatment + s(Experiment.Day, by = Treatment) + s(Replicate, bs = "re"),
                    family = binomial(link = "logit"), data = cdata)

##Model Validation
#Prop_Alive_TC Model
sim_tc <- simulateResiduals(bgam_TC)
plot(sim_tc)
testDispersion(sim_tc)
bgam_TC$deviance/bgam_TC$df.residual

#Prop_Alive_TF Model
sim_tf <- simulateResiduals(bgam_TF)
plot(sim_tf)
testDispersion(sim_tf)
bgam_TF$deviance/bgam_TF$df.residual

#Prop_Alive_M Model
sim_m <- simulateResiduals(bgam_M)
plot(sim_m)
testDispersion(sim_m)
bgam_M$deviance/bgam_M$df.residual

#Prop_PF Model
par(mfrow=c(2,2))
gam.check(bgam_Prop_PF)
bgam_Prop_PF$deviance/bgam_Prop_PF$df.residual

#Binomial GLMM - For Model Comparison
bglmm_TC <- glmmTMB(cbind(Total.Copepods, Dead_TC) ~ Experiment.Day * Treatment + (1|Replicate),
                    family = binomial(link = "logit"), data = cdata)

bglmm_TF <- glmmTMB(cbind(Current_TF, Dead_TF) ~ Experiment.Day * Treatment + (1|Replicate),
                    family = binomial(link = "logit"), data = cdata)

bglmm_M <- glmmTMB(cbind(Nb.of.M, Dead_M) ~ Experiment.Day * Treatment + (1|Replicate),
                   family = binomial(link = "logit"), data = cdata)

bglmm_Prop_PF <- glmmTMB(cbind(Nb.of.PF, Current_TF - Nb.of.PF) ~ Experiment.Day * Treatment + (1|Replicate),
                         family = binomial(link = "logit"), data = cdata)

#AIC Model Comparison to confirm final model
#TC
AIC(bgam_TC, bglmm_TC)
#TF
AIC(bgam_TF, bglmm_TF)
#M
AIC(bgam_M, bglmm_M)
#Preg Fem
AIC(bgam_Prop_PF, bglmm_Prop_PF)
##GAM model had lower AIC value in all cases so GAM model is final.

#Model Summaries
#Prop_Alive_TC Model
summary(bgam_TC)

#Prop_Alive_TF Model
summary(bgam_TF)

#Prop_Alive_M Model
summary(bgam_M)

#Prop_PF Model
summary(bgam_Prop_PF)

#Function to check models in more detail
check_model <- function(model, model_name) {
  cat(model_name)
  
  # Basic model check
  gam.check(model)
  
  # Plot smooth terms
  plot(model, pages = 1, main = paste("Smooth terms -", model_name))
  
  # Residuals vs fitted
  plot(fitted(model), residuals(model, type = "pearson"), 
       main = paste("Residuals vs Fitted -", model_name),
       xlab = "Fitted values", ylab = "Pearson residuals")
  abline(h = 0, col = "red")
}

#Check Models
par(mfrow = c(2,2))
check_model(bgam_TC, "Proportion of surviving copepods")
check_model(bgam_TF, "Proportion of surviving females") 
check_model(bgam_M, "Proportion of surviving males")
check_model(bgam_Prop_PF, "Proportion of pregnant females")

##look at variance components
gam.vcomp(bgam_TC)
gam.vcomp(bgam_TF)
gam.vcomp(bgam_M)
gam.vcomp(bgam_Prop_PF)

##Post Model Analysis 
##Pairwise comparison using emmeans at time points
#make function for each day so its not averaged over time
emmeans_by_day <- function (model, days = c(0, 5, 10, 15, 20), type = "response"){
  results_list <- list()
  
  for (day in days){
    em_result <- emmeans(model, pairwise ~ Treatment, 
                         at = list(Experiment.Day = day),
                         type = type)
    
    em_table <- as.data.frame(em_result$emmeans)
    em_table$Day <- rep(day, nrow(em_table))
    
    pairs_table <- as.data.frame(em_result$contrasts)
    pairs_table$Day <- rep(day, nrow(pairs_table))
    
    results_list[[paste0("Day_", day)]] <- list(
      emmeans = em_table,
      contrasts = pairs_table
    )
  }
  
  all_emmeans <- do.call(rbind, lapply(results_list, function(x) x$emmeans))
  all_contrasts <- do.call(rbind, lapply(results_list, function(x) x$contrasts))
  
  all_emmeans<- all_emmeans[,c("Day", setdiff(names(all_emmeans), "Day"))]
  all_contrasts <- all_contrasts[,c("Day", setdiff(names(all_contrasts), "Day"))]
  
  return(list(
    emmeans = all_emmeans,
    contrasts = all_contrasts,
    by_day = results_list
  ))
  
}

results_TC <- emmeans_by_day(bgam_TC)
results_TF <- emmeans_by_day(bgam_TF)
results_M <- emmeans_by_day(bgam_M)
results_PPF <- emmeans_by_day(bgam_Prop_PF)

write.csv(results_TC$emmeans, "TC_emmeans_by_day.csv", row.names = FALSE)
write.csv(results_TC$contrasts, "TC_contrasts_by_day.csv", row.names = FALSE)

write.csv(results_TF$emmeans, "TF_emmeans_by_day.csv", row.names = FALSE)
write.csv(results_TF$contrasts, "TF_contrasts_by_day.csv", row.names = FALSE)

write.csv(results_M$emmeans, "M_emmeans_by_day.csv", row.names = FALSE)
write.csv(results_M$contrasts, "M_contrasts_by_day.csv", row.names = FALSE)

write.csv(results_PPF$emmeans, "PPF_emmeans_by_day.csv", row.names = FALSE)
write.csv(results_PPF$contrasts, "PPF_contrasts_by_day.csv", row.names = FALSE)




##Model Plots
#Plots with Gam models with overlaying data points
#Custom colours for Treatments - colour blind friendly
custom_colours <- c(
  "Control" = magma(100)[3],
  "SS"     = plasma(100)[28],
  "DM"  = magma(100)[55],
  "SS+DM"    = magma(100)[75]
)


#Create function to plot model predictions and overlay data
#pgwd = plot gam with data
pgwd <- function(model, data, response_var, title, y_lab) {
  #Create all combinations of variables
  pred_data <- expand.grid(
    Experiment.Day = seq(min(data$Experiment.Day), max(data$Experiment.Day), length.out = 100),
    Treatment = unique(data$Treatment),
    Replicate = unique(data$Replicate)[1]
  )
  
  #Make predictions from models
  pred_data$fit <- predict(model, newdata = pred_data, type = "response")
  
  #Get SE for predictions
  pred_data$se <- predict(model, newdata = pred_data, type = "response", se.fit = TRUE)$se.fit
  
  #Confidence Intervals
  pred_data$lower <- pmax(0, pred_data$fit - 1.96 * pred_data$se)
  pred_data$upper <- pmin(1, pred_data$fit + 1.96 * pred_data$se)
  
  #Plot function
  p <- ggplot() +
    #Add predicted data
    geom_ribbon(data = pred_data,
                aes(x = Experiment.Day, ymin = lower, ymax = upper, fill = Treatment), alpha = 0.15) +
    geom_line(data = pred_data,
              aes(x = Experiment.Day, y = fit, colour = Treatment), linewidth = 0.75) +
    #Add OG data points
    geom_point(data = cdata,
               aes(x = Experiment.Day, y = !!sym(response_var), colour = Treatment), alpha = 0.7) +
    #Custom colours
    scale_colour_manual(values = custom_colours)+
    scale_fill_manual(values = custom_colours)+
    
    labs(
      title = title,
      x = "Experiment Day",
      y = y_lab,
      colour = "Treatment",
      fill = "Treatment"
    ) +
    theme_minimal()+
    theme(plot.title = element_text( hjust = 0.5, family = "serif", size = 5),
          legend.position = "bottom",
          panel.grid.minor = element_blank(),
          text = element_text(family = "serif", size = 11.5),
          axis.text = element_text(family = "serif", size = 11),
          axis.title = element_text(family = "serif", size = 11.5),
          legend.text = element_text(family = "serif", size = 12),
          legend.title = element_text(family = "serif", size = 12))+
    ylim(0,1)
  return(p)
}

P1 <- pgwd(bgam_TC, cdata, "Prop_Alive_TC", " ", "Proportion of Surviving Copepods")
P2 <- pgwd(bgam_TF, cdata, "Prop_Alive_TF", " ", "Proportion of Surviving Females")
P3 <- pgwd(bgam_M, cdata, "Prop_Alive_M", " ", "Proportion of Surviving Males")
P4 <- pgwd(bgam_Prop_PF, cdata, "Prop_PF", " ", "Proportion of Females Pregnant")

grid.arrange(P1,P4,P3,P2, ncol = 2)

##Combined Plot
shared_legend <- ggplot_gtable(ggplot_build(P1))
legend_index <- which(sapply(shared_legend$grobs, function(x) x$name) == "guide-box")
legend <- shared_legend$grobs[[legend_index]]

P1_nl <- P1 + theme(legend.position = "none")
P2_nl <- P2 + theme(legend.position = "none")
P3_nl <- P3 + theme(legend.position = "none")
P4_nl <- P4 + theme(legend.position = "none")

plots_grid <- arrangeGrob(P1_nl,P4_nl,P2_nl,P3_nl, ncol = 2)
combined_plot <- arrangeGrob(plots_grid, legend, ncol = 1, heights = c(10,1))
dev.off()
grid.draw(combined_plot)

#add labels for plots - (A),(B),(C),(D)
grid.text("(A)", x = 0.02, y = 0.985, gp = gpar(fontsize = 13, fontfamily = "serif"))
grid.text("(B)", x = 0.52, y = 0.985, gp = gpar(fontsize = 13, fontfamily = "serif"))
grid.text("(C)", x = 0.02, y = 0.53, gp = gpar(fontsize = 13, fontfamily = "serif"))
grid.text("(D)", x = 0.52, y = 0.53, gp = gpar(fontsize = 13, fontfamily = "serif"))


labelled_plot <- gTree(children = gList(
  combined_plot,
  textGrob("(A)", x = 0.02, y = 0.985, gp = gpar(fontsize = 13, fontfamily = "serif")),
  textGrob("(B)", x = 0.52, y = 0.985, gp = gpar(fontsize = 13, fontfamily = "serif")),
  textGrob("(C)", x = 0.02, y = 0.53, gp = gpar(fontsize = 13, fontfamily = "serif")),
  textGrob("(D)", x = 0.52, y = 0.53, gp = gpar(fontsize = 13, fontfamily = "serif"))
))

grid.draw(labelled_plot)  


#save as png with high resolution
png("labelled_plot.png", width = 8, height = 6.5, units = "in", res = 300)
grid.draw(labelled_plot)
dev.off()
