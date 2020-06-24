# CleanCNA
Calling subclonal copy number aberrations (CNA) with Battenberg, DPClust, and iterative recalls to provide a high quality and reproducible set of copy number calls.

Subclonal copy number is called using Battenberg, DPClust takes Battenberg's output and clusters mutations by cancer cell fraction (CCF). Metrics are calculated to assess the accuracy of the copy number call. Adjusted purity (ρ) and ploidy (ψ) parameters are calculated if required, and the copy number profile is iteratively refitted where necessary. Detailed metrics for inclusion or exclusion of samples is provided.   
