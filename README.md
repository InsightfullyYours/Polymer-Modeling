# Polymer-Modeling
MATLAB M-Files to model polymer scission and combination, intended to be compared to HPLC data.

This is a set of MATLAB m-files that I wrote between 2006 and 2012 to model the scission and combination reactions of arbitrary polymer distributions.  My thesis work (Photochemical Crosslinking Reactions in Polymers, Columbia University, 2012) was partially based upon it.  It takes arbitration polymer size distributions (as the number fraction of chains (yaxis) at each chain size (number of monomers), (xaxis)) and models scission and combination reactions based upon a number of theories.  This code/model is polymer-independent, in that the chemistry does not need to be understood to model.  Instead, the scission and combination events per unit time are determined from experimental data, so all mechanisms are included.  The experimental data was taken as HPLC data using RID, UV, and MALLS detectors. None of that data is included here.  The other repository (HPLC-Analysis-Code) includes a number of m-files and scripts that were used on the data itself in order to prepare it to the results of this modeling. There may be some overlap, duplication, and inter-dependencies between the two repositories.  I attempted to minimize this.

Hamielec.m is the main model file.  It changes the input data to number and weight fracitons, then approximates the degree of scission u and from there models the scission.  It then uses the modeled scission to calculate end-to-end combinations, the formation of 3-arm stars, and the formation of 4-arm stars.  Once the scission and combination have been modeled, it fits the results to the data.

FitShape.m scales the control peak by exp(-ru), fits it to the data (the left side of the Hamielec scission equation) and fits u to determine the degree of scission.

FitPeaks.m finds the fit of the combination peaks that best fit the data, ith or without g-factor incorporation. It minimizes the sum of squared errors between the sum of the 3 fitted peaks and the data.

MacroRadicalCombination.m calculates the macroradical recombination of any two input distributions and outputs its for further use.

AreaofPeaks.m calculates the area of the input peaks using trapezoidal Riemann sums.  It is used in the determination of peak height/area and relative reactions.

Probabilities.m ses the probabilities of each radical species to calculate the fraction reacted.

CombinationScissionwithRg.m is the heart of the simulation.  It calculates on either a user-inputted distribution or a user-defined one.  Uncomment the gaussian definition and modify it if you'd like a user-defined distribution.

RadiusOfGyration.m is a subfunction that implements Zimm and Stockmeyer's definition of radius of gyration from macromolecule structure.

Everyone is free to use them as the license dictates, but considering the specificity of the application and the wide variety of available software and analysis packages, they may not be relevant to your systems. They are, in general, not sufficiently explained in comments and I will not be updating the comments. I have uploaded them for nostalgic and demonstration reasons.
