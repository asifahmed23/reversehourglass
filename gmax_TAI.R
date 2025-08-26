set.seed(101)

# set working directory
setwd("C:/Users/sami001/OneDrive - Wageningen University & Research/TAI_tpc/Resubmission_v.2 TPC/script Gmax")


# load myTAI
library(myTAI)
library(tidyverse)
library(viridis)

# load data from MNodine and Renake
load(file = "Gmax_TAI.Rdata")


# get the TAI value for each sample
TAI(tpmCountAvgGm)


# Visualize the TAI profile of correctly assigned Phylostrata
PlotSignature( ExpressionSet = tpmCountAvgGm, 
               p.value       = TRUE, permutations = 50000, TestStatistic = "ReverseHourglassTest",
               modules = list(early = 1:2, mid = 3:6, late = 7:9),
               ylab = "Transcriptome age index (TAI)\n",
               xlab = "\nOntogeny") +
              #ylim(0,6) +
              theme_bw()+
              theme(
                legend.position = "none",
                strip.text = element_text(size = 14),
                axis.text.x = element_text(size = 14, angle = 45, hjust=1, color = "black"), 
                axis.text.y = element_text(size = 14, color = "black"),
                axis.title.x = element_text(size = 16), # set x-axis label size
                axis.title.y = element_text(size = 16), # set y-axis label size
                legend.text = element_text(size = 14),
                legend.title = element_text(size = 16),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank()
              )


#ggsave("Gmax_TAI_value.png", device="png", dpi = 500, height = 12, width = 20)

#-------------------------------------------------------------------------------
# Stability of reverse hourglass pattern
#-------------------------------------------------------------------------------

# check the stability of the reverse hourglass pattern to different data transformations
tfStability(tpmCountAvgGm,
            TestStatistic = "ReverseHourglassTest",
            transforms = c("none", "sqrt", "log2", "rank", "squared"),
            modules = list(early = 1:2, mid = 3:6, late = 7:9),
            permutations = 20000,
            pseudocount = 1)

# check of the underlying statistical assumptions are met (also the Cullen-Frey plot)
ReverseHourglassTest(
  tpmCountAvgGm,
  permutations = 20000,
  modules = list(early = 1:2, mid = 3:6, late = 7:9),
  plotHistogram = TRUE,
  runs = 10,
  parallel = FALSE,
  lillie.test	= TRUE,
  custom.perm.matrix = NULL
)



#-------------------------------------------------------------------------------
# PS contribution
#-------------------------------------------------------------------------------

# Plot the individual phylostrata contribution
PlotContribution(tpmCountAvgGm,  legendName = "PS")


