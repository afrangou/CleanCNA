computed.expected.peak.loc <- function(multiplicity, purity, ploidy) {
  (multiplicity * purity) / (2 * (1 - purity) + ploidy * purity) # from ASCAT equation
}
