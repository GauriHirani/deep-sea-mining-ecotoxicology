# Copepod Survival Under Deep-Sea Mining Sediment Exposure

This project investigates the sex-specific effects of Clarion Clipperton Zone (CCZ) sediment on the survival and reproductive output of the marine copepod *Tigriopus californicus*, exposed via both the dissolved and particulate pathway. The findings inform understanding of the potential ecological impacts of deep-sea mining midwater discharges, which contain metal-rich CCZ sediment and associated trace metals.

Completed as part of an MSc research project at Imperial College London.

## Data

`LTEXP.csv` contains the survival data used in this analysis. Full metadata, including variable definitions, experimental methods, and data collection details, is provided in `Metadata_LTEXP.xlsx`.

## Setup

Ensure `LTEXP.csv` is in the same working directory as the script before running.

## Statistical Methods

Survival and reproductive data were explored to check for outliers, distribution shape, and collinearity between predictors, informing choice of statistical model. Binomial GAMs were fitted to model the proportion of survival (total copepods, females, males, pregnant females) across treatment and experimental day, with replicate included as a random effect. GAMs were compared against binomial GLMMs using AIC, with GAM providing a better fit in all cases. Model assumptions were validated using DHARMa residual diagnostics. Pairwise comparisons between treatments at specific experimental days were then carried out using estimated marginal means (emmeans).

## Key Findings

Sediment-containing treatments (SS and SS+DM) caused the highest copepod mortality, with the particulate pathway appearing to be the dominant route of toxicity. Survival effects were sex-specific, with males showing lower tolerance to stress than females. Reproductive output was significantly reduced in sediment-exposed treatments compared to Control and DM, suggesting suspended sediment impairs reproduction more than dissolved metal alone.

## Tools
R, dplyr, tidyr, mgcv, glmmTMB, DHARMa, emmeans, ggplot2, viridis
