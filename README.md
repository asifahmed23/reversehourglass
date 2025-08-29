# reversehourglass
All codes related to the article "The angiosperm seed life cycle follows a developmental reverse hourglass" are submitted here. 

1. The "TAI_calculation.R" script contains codes to reproduce the main TAI patterna and phylostrata contribution for Arabidopsis, B. napus, tomato, and maize. The script also checks for the effect of different data normalization methods on the TAI pattern seen in the Arabidopsis seed life cycle.

2. The "gmax_TAI.R" script contains the codes to produce the TAI pattern and phylostrata contribution of Glycine max seed life cycle from the rnaseq dataset published by Chen et al., 2024 (https://www.maxapress.com/article/id/675137defa6c5848635e86a1). The script also checks for data normalization effects and the Cullen-Frey plot to check if statistical assumptions are met.
3. The "Data check Lotharukpong etal. 2024" directory contains a script that performs data checks on the RNA-seq data used in the Lotharukpong etal. 2024 study (https://www.nature.com/articles/s41586-024-08059-8). The data checks for the Cullen-Frey plot using the FlatLineTest and ReductiveHourglassTest, and also checks the stability of the ReductiveHourglassTest test in Ficus distichus and F. serratus after different data transformations. Same parameters were used as the original study except that we increased the number of permuations to 20000.
