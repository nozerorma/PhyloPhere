# Set of RERConverge functions used by CT Resample

function (namedvec, treewithbranchlengths) 
{
  vec = simulatevec(namedvec, treewithbranchlengths)
  simsorted = sort(vec)
  realsorted = sort(namedvec)
  l = length(simsorted)
  c = 1
  while (c <= l) {
    simsorted[c] = realsorted[c]
    c = c + 1
  }
  simsorted
}

function (namedvec, treewithbranchlengths) 
{
  library("geiger")
  rm = ratematrix(treewithbranchlengths, namedvec)
  sims = sim.char(treewithbranchlengths, rm, nsim = 1)
  nam = rownames(sims)
  s = as.data.frame(sims)
  simulatedvec = s[, 1]
  names(simulatedvec) = nam
  vec = simulatedvec
  vec
}
