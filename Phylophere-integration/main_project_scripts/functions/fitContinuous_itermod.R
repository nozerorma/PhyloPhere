fitCont_alt <- function (phy, dat, SE = 0, model = c("BM", "OU", "EB", "rate_trend",
  "lambda", "kappa", "delta", "mean_trend", "white"), bounds = list(),
  control = list(method = c("subplex", "L-BFGS-B"), niter = 1000,
    FAIL = 1e+200, hessian = FALSE, CI = 0.95), ncores = NULL,
  ...)
{
  td = treedata(phy, dat)
  if (nrow(td$data) != length(unique(rownames(td$data))))
    stop("Multiple records per tip label")
  phy = td$phy
  dat = td$data
  dd = dim(dat)
  trts = dd[2]
  if (trts > 1) {
    nm = colnames(dat)
    res = lapply(1:trts, function(idx) {
      fitContinuous(phy, dat[, idx], SE = SE, model = model,
        bounds = bounds, control = control)
    })
    names(res) = nm
    class(res) = c("gfits", class(res))
    return(res)
  }
  else {
    dat = dat[, 1]
  }
  ct = list(method = c("subplex", "L-BFGS-B"), niter = 1000,
    FAIL = 1e+200, hessian = FALSE, CI = 0.95)
  if (any(tmp <- !names(control) %in% names(ct)))
    warning("Unexpected 'control' parameters:\n\t", paste(names(control)[tmp],
      collapse = "\n\t"), sep = "")
  control = control[which(!tmp)]
  if ("method" %in% names(control))
    control$method = match.arg(control$method, c("Nelder-Mead",
      "BFGS", "CG", "L-BFGS-B", "SANN", "Brent", "subplex"),
      several.ok = TRUE)
  ct[names(control)] = control
  if (ct$niter < 2)
    stop("'niter' must be equal to or greater than 2")
  ct$hessian_P = 1 - ct$CI
  if (length(model) == 1) {
    if (model == "trend")
      model <- "rate_trend"
    else if (model == "drift")
      model <- "mean_trend"
  }
  model = match.arg(model, c("BM", "OU", "EB", "rate_trend",
    "lambda", "kappa", "delta", "mean_trend", "white"))
  if (model == "OU" & !is.ultrametric(phy)) {
    warning("Non-ultrametric tree with OU model, using VCV method.")
    lik = ou.lik(phy, dat, SE, model, ...)
  }
  else {
    con = list(method = "pruning", backend = "C")
    con[names(control)] = control
    lik = bm.lik(phy, dat, SE, model, ...)
  }
  attr(lik, "model") = model
  argn = argn(lik)
  mn = c(-500, -500, (log(10^(-5))/max(node.depth.edgelength(phy))),
    -100, -100, -500, -500, -500, -500)
  mx = c(100, 1, -1e-06, 100, 100, 0, 0, log(2.999999), 100)
  bnds = as.data.frame(cbind(mn, mx))
  bnds$typ = c("exp", "exp", "nat", "nat", "nat", "exp", "exp",
    "exp", "exp")
  rownames(bnds) = c("sigsq", "alpha", "a", "drift", "slope",
    "lambda", "kappa", "delta", "SE")
  bnds$model = c("BM", "OU", "EB", "mean_trend", "rate_trend",
    "lambda", "kappa", "delta", "SE")
  typs = bnds[argn, "typ"]
  if (length(bounds) > 0) {
    mm = match(names(bounds), rownames(bnds))
    if (any(is.na(mm))) {
      warning("Unexpected 'bounds' parameters:\n\t", paste(names(bounds)[is.na(mm)],
        collapse = "\n\t"), sep = "")
    }
    mm = mm[!is.na(mm)]
    if (length(mm)) {
      for (i in 1:length(mm)) {
        ww = mm[i]
        tmp = sort(bounds[[i]])
        if (bnds$typ[ww] == "exp") {
          if (any(tmp == 0))
            tmp[tmp == 0] = exp(-500)
          bnd = log(tmp)
        }
        else {
          bnd = tmp
        }
        bnds[ww, c("mn", "mx")] = bnd
      }
    }
  }
  if (any(!is.finite(as.matrix(bnds[, c("mn", "mx")])))) {
    stop("All bounds should be finite")
  }
  par = argn[1]
  xx = function(p) {
    pars = ifelse(typs == "exp", exp(p), p)
    tmp = -lik(pars, root = "max")
    if (is.infinite(tmp))
      tmp = ct$FAIL
    if (is.na(tmp))
      tmp = ct$FAIL
    tmp
  }
  boxconstrain = function(f, lower, upper, fail.value) {
    function(x) {
      if (any(x < lower | x > upper))
        fail.value
      else f(x)
    }
  }
  f = boxconstrain(xx, bnds[argn, "mn"], bnds[argn, "mx"],
    fail.value = ct$FAIL)
  if (par %in% c("alpha", "lambda", "delta", "kappa")) {
    bmstart = try(.bm.smartstart(phy, dat), silent = TRUE)
    if (inherits(bmstart, "try-error"))
      bmstart = 0.01
    bmstart = log(bmstart)
  }
  mm = matrix(NA, nrow = ct$niter, ncol = length(argn) + 2)
  mt = character(ct$niter)
  min = bnds[argn, "mn"]
  max = bnds[argn, "mx"]
  fxopt = .get.parallel(ncores)
  out = fxopt(1:ct$niter, function(i) {
    bnds$st = sapply(1:nrow(bnds), function(x) runif(1,
      bnds$mn[x], bnds$mx[x]))
    start = bnds[argn, "st"]
    if (par == "alpha") {
      oustart = log(.ou.smartstart(dat, unlist(exp(bnds["alpha",
        c("mx", "mn")]))))
      if (i == 1 | runif(1) < 0.25)
        start[match(c("sigsq", par), argn)] = c(bmstart,
          oustart)
      if (runif(1) < 0.5)
        start[match(c("sigsq", par), argn)] = c(0, oustart)
    }
    if (par %in% c("lambda", "delta", "kappa")) {
      ww = match(par, rownames(bnds))
      if (runif(1) < 0.5) {
        if (runif(1) < 0.5) {
          start[match(c("sigsq", par), argn)] = c(bmstart,
            bnds$mx[ww])
        }
        else {
          start[match(c("sigsq", par), argn)] = c(bmstart,
            bnds$mn[ww])
        }
      }
    }
    if (par == "white") {
      if (runif(1) < 0.5) {
        start[match("sigsq", argn)] = var(dat)
      }
    }
    names(start) = argn
    if (length(argn) == 1) {
      method = "Brent"
    }
    else {
      method = sample(ct$method, 1)
    }
    if (method == "subplex") {
      op <- try(suppressWarnings(subplex(par = start,
        fn = f, control = list(reltol = .Machine$double.eps^0.25,
          parscale = rep(0.1, length(argn))), hessian = ct$hessian)),
        silent = TRUE)
    }
    else {
      op <- try(suppressWarnings(optim(par = start, fn = f,
        upper = max, lower = min, method = method, hessian = ct$hessian)),
        silent = TRUE)
    }
    if (!inherits(op, "try-error")) {
      op$method = method
      op$value = -op$value
      names(op)[names(op) == "value"] = "lnL"
      names(op$par) = argn
      op$par = sapply(1:length(typs), function(x) if (typs[x] ==
        "exp")
        return(exp(op$par[x]))
      else return(op$par[x]))
      op$mm = c(op$par, op$lnL, op$convergence)
    }
    else {
      op = list(par = structure(rep(NA, length(argn)),
        names = argn), lnL = -Inf, convergence = 1,
        method = "FAIL")
    }
    op
  })
  for (i in 1:length(out)) {
    cur = out[[i]]
    if (cur$method != "FAIL")
      mm[i, ] = cur$mm
    mt[i] = cur$method
  }
  res = mm
  colnames(res) = c(argn, "lnL", "convergence")
  rownames(res) = mt
  colnames(mm) = c(argn, "lnL", "convergence")
  conv = mm[, "convergence"] == 0
  mm = mm[, -which(colnames(mm) == "convergence")]
  valid = apply(mm, 1, function(x) !any(is.na(x)))
  if (sum(valid & conv) >= 1) {
    mm = matrix(mm[valid, ], nrow = sum(valid), dimnames = dimnames(mm))
    mt = mt[valid]
    out = out[valid]
    mm = mm[z <- min(which(mm[, "lnL"] == max(mm[, "lnL"]))),
      ]
  }
  else {
    z = NA
    mm = c(rep(NA, length(argn)), -Inf)
    names(mm) = c(argn, "lnL")
  }
  zz = mm[-which(names(mm) %in% c("lnL"))]
  mm = as.list(mm)
  tmp = lik(unlist(mm[argn]))
  mm = c(mm[argn], z0 = attributes(tmp)$ROOT.MAX, mm[names(mm)[!names(mm) %in%
    argn]])
  mm$method = ifelse(is.na(z), NA, mt[z])
  mm$k = length(argn) + 1
  if (ct$hessian) {
    hessian = out[[z]]$hessian
    CI = .bnd.hessian(hessian, zz, typs, ct$hessian_P)
    if (!all(is.na(CI))) {
      if (is.constrained(lik)) {
        CI = rbind(lik(CI[1, ], pars.only = TRUE), rbind(lik(CI[2,
          ], pars.only = TRUE)))
      }
      dimnames(hessian) = NULL
      rownames(CI) = c("lb", "ub")
    }
  }
  else {
    hessian = NULL
    CI = NULL
  }
  range = as.data.frame(cbind(min, max))
  range$typ = typs
  range$mn = ifelse(range$typ == "exp", exp(range$min), range$min)
  range$mx = ifelse(range$typ == "exp", exp(range$max), range$max)
  par = mm[argn]
  rownames(range) = argn
  chk = sapply(1:length(par), function(idx) {
    p = par[[idx]]
    if (!is.na(p)) {
      return((p <= range$mn[idx] | p >= range$mx[idx]))
    }
    else {
      return(FALSE)
    }
  })
  if (any(chk)) {
    warning(paste("\nParameter estimates appear at bounds:\n\t",
      paste(names(par)[chk], collapse = "\n\t", sep = ""),
      sep = ""))
  }
  mm = .aic(mm, n = length(dat))
  mm$CI = CI
  mm$hessian = hessian
  res = list(lik = lik, bnd = range[, c("mn", "mx")], res = res,
    opt = mm)
  class(res) = c("gfit", class(res))
  return(res)
}
