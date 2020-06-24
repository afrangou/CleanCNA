reestimate.purity <- function(loc, multiplicity, ploidy) {
  2 * loc / (multiplicity + loc * (2 - ploidy)) # from rearrangement of ASCAT equation
}
