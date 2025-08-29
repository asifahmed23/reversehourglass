set.seed(101)

# set working directory
setwd("C:/Users/sami001/OneDrive - Wageningen University & Research/PhD/Work Package 3 TAI/Analysis")


# load myTAI
library(myTAI)
library(tidyverse)


# load merged TPM counts for all species
load(file = "Github/TPM_counts_all_species.Rdata")


#-------------------------------------------------------------------------------
#---------------------------Arabidopsis thaliana--------------------------------
#-------------------------------------------------------------------------------

# calculate the TAI values
TAI(ath_expr)



# Visualize the TAI profile of correctly assigned Phylostrata
PlotSignature(ExpressionSet = ath_expr, 
              p.value       = TRUE, permutations = 50000, TestStatistic = "ReverseHourglassTest",
              modules = list(early = 1:7, mid = 8:17, late = 18:21)) +
              theme_bw()+
              theme(
                legend.position = "none",
                strip.text = element_text(size = 14),
                axis.text.x = element_text(size = 16, angle = 45, hjust=1), 
                axis.text.y = element_text(size = 16),
                axis.title.x = element_text(size = 16), # set x-axis label size
                axis.title.y = element_text(size = 16), # set y-axis label size
                legend.text = element_text(size = 20),
                legend.title = element_text(size = 20),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank()
              )

# check the stability of the reverse hourglass pattern to different data transformations
tfStability(ath_expr,
            TestStatistic = "ReverseHourglassTest",
            transforms = c("none", "sqrt", "log2", "rank", "squared"),
            modules = list(early = 1:7, mid = 8:17, late = 18:21),
            permutations = 50000,
            pseudocount = 1)

# Check the Cullen-Frey plot for FlatLineTest and ReverseHourglassTest
FlatLineTest(
  ath_expr,
  permutations = 50000,
  plotHistogram = T,
  runs = 5)
ReverseHourglassTest(
  ath_expr,
  permutations = 50000,
  modules = list(early = 1:7, mid = 8:17, late = 18:21),
  plotHistogram = TRUE,
  runs = 5,
  parallel = FALSE,
  lillie.test	= TRUE,
  custom.perm.matrix = NULL
)



# Plot relative expression of genes
PlotRE(ath_expr, Groups = list(1:4, 5:10),adjust.range = FALSE)

# plot the mean relative expression levels of phylostratum or divergence-stratum classes as barplot
PlotBarRE(ath_expr, Groups = list(1:4, 4:16),adjust.range = FALSE)+
  theme_bw()+
  theme(
    strip.text = element_text(size = 14),
    axis.text.x = element_text(size = 16, angle = 45, hjust=1), 
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16), # set x-axis label size
    axis.title.y = element_text(size = 16), # set y-axis label size
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )



# Plot the individual phylostrata contribution
PlotContribution(ath_expr,  legendName = "PS")



#-------------------------------------------------------------------------------
#---------------------------Brassica napus--------------------------------------
#-------------------------------------------------------------------------------


# calculate the TAI values
TAI(bna_expr)



# Visualize the TAI profile of correctly assigned Phylostrata
PlotSignature(ExpressionSet = bna_expr, 
              p.value       = TRUE, permutations = 50000, TestStatistic = "ReverseHourglassTest",
              modules = list(early = 1:7, mid = 8:15, late = 16:20)) +
  theme_bw()+
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14),
    axis.text.x = element_text(size = 16, angle = 45, hjust=1), 
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 16), # set x-axis label size
    axis.title.y = element_text(size = 16), # set y-axis label size
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank()
  )

# Plot the individual phylostrata contribution
PlotContribution(bna_expr,  legendName = "PS")


#-------------------------------------------------------------------------------
#---------------------------Solanum lycopersicum--------------------------------
#-------------------------------------------------------------------------------

# get the TAI value for each sample
TAI(sly_expr)


# Visualize the TAI profile of correctly assigned Phylostrata
PlotSignature( ExpressionSet = sly_expr, 
               p.value       = TRUE, permutations = 50000, TestStatistic = "ReverseHourglassTest",
               modules = list(early = 1:2, mid = 3:9, late = 10:14),
               ylab = "Transcriptome age index (TAI)\n",
               xlab = "\nDays after flowering (DAF)") +
               theme_bw()+
               theme(
                  legend.position = "none",
                  strip.text = element_text(size = 14),
                  axis.text.x = element_text(size = 20, angle = 45, hjust=1), 
                  axis.text.y = element_text(size = 20),
                  axis.title.x = element_text(size = 20), # set x-axis label size
                  axis.title.y = element_text(size = 20), # set y-axis label size
                  legend.text = element_text(size = 20),
                  legend.title = element_text(size = 20),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank()
                )

# Plot the individual phylostrata contribution
PlotContribution(sly_expr,  legendName = "PS")


#-------------------------------------------------------------------------------
#----------------------------------Zea mays-------------------------------------
#-------------------------------------------------------------------------------

TAI(zma_expr)


# Visualize the TAI profile of correctly assigned Phylostrata
PlotSignature( ExpressionSet = zma_expr, 
               p.value       = TRUE, permutations = 50000, TestStatistic = "ReverseHourglassTest",
               modules = list(early = 1:6, mid = 7:16, late = 17:21),
               ylab = "Transcriptome age index (TAI)\n",
               xlab = "\nDays after pollination (DAP)") +
               theme_bw()+
               theme(
                  legend.position = "none",
                  strip.text = element_text(size = 14),
                  axis.text.x = element_text(size = 20, angle = 45, hjust=1), 
                  axis.text.y = element_text(size = 20),
                  axis.title.x = element_text(size = 20), # set x-axis label size
                  axis.title.y = element_text(size = 20), # set y-axis label size
                  legend.text = element_text(size = 20),
                  legend.title = element_text(size = 20),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank()
                )

# Plot the individual phylostrata contribution
PlotContribution(zma_expr,  legendName = "PS")


