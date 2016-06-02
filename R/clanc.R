## These are all the R functions required for using ClaNC.  Note the GUI interface has been 
## deprecated.  So all use must be done through the R command line.  See the 'example_session.txt' 
## file for an example of how to do this.

"cvClanc" <- function(data, id, prior = "equal", active = 1:10, gui = F, prntFrm = NULL) {
  ## data: expression data.  (m x n) matrix of class numeric.
  ## id: class id's.  n vector of class numeric.
  ## prior: either a vector of length p, "equal", or "class".  if "equal", then equal probabilities 
  ##   will be used for each class.  if "class", then the proportions of each class in the 
  ##   training data will be used as the prior.
  ## active: how many active features to consider?  can either be a single number or a vector 
  ##   containing a range of values to consider.
  ## gui: indicates whether call is coming from gui.  typical user should leave default.
  ## prntFrm: frame in which to print if called from gui.  typical user should leave default.

  cvIdx = balancedFolds(id, 5)

  m = nrow(data)
  n = ncol(data)
  p = length(unique(id))
  nn = as.numeric(table(id))
  folds = length(cvIdx)

  if(is.numeric(prior)) {
    if(length(prior) != p | sum(prior) != 1)
      stop("Invalid prior.")
    pi.k = prior
  } else {
    if(prior == "equal")
      pi.k = rep(1 / p, p)
    else if(prior == "class")
      pi.k = nn / n
  }

  if(is.matrix(active)) {
    d = nrow(active)
  } else {
    d = length(active)
  }

  if(gui) {
    printClanc = function(msg) {
      assign("shareMsg", msg, prntFrm)
      eval(expression(postMsg(shareMsg)), prntFrm)
    }
  } else {
    printClanc = cat
  }

  cv.error = array(rep(0, d * folds * p), dim = c(d, folds, p))

  cv.err.cnt.cls = matrix(NA, nrow = d, ncol = p)
  cv.err.prpn.cls = matrix(NA, nrow = d, ncol = p)
  n.features.cls = matrix(NA, nrow = d, ncol = p)
  n.features.ttl = rep(NA, d)

  ID = model.matrix(~ factor(id) - 1)
  dimnames(ID) = list(NULL, names(nn))

  ## cross validation
  printClanc("CV:")
  for(i in 1:folds) {
    printClanc(i)

    ## form initial statistics
    v = length(cvIdx[[i]])
    X = data[, -cvIdx[[i]]]
    Y = data[, cvIdx[[i]]]
    jd = id[-cvIdx[[i]]]
    truth = id[cvIdx[[i]]]
    JD = ID[-cvIdx[[i]], ]
        
    mm = table(jd)
    m.k = sqrt(1 / mm - 1 / (n - v))

    ## pooled standard deviations
    p.sd = pooledSDClanc(X, JD)
        
    ## class- and overall-centroids
    cntrd.k = scale(X %*% JD, center = F, scale = mm)
    cntrd.o = drop(X %*% rep(1 / (n - v), n - v))

    ## form statistics and order them
    d.k = scale((cntrd.k - cntrd.o) / p.sd, center = F, scale = m.k)
    d.k.ord = orderStatsClanc(d.k = d.k)

    ## select genes, update inactive centroid components
    for(j in 1:d) {
      if(is.matrix(active))
        aa = active[j, ]
      else
        aa = active[j]

      selected = selectClanc(d.k = d.k, d.k.ord = d.k.ord, active = aa)
      active.idx = (1:m)[drop(selected %*% rep(1, p)) != 0]

      cntrds = cntrd.k[active.idx, ]
      for(k in 1:p)
        cntrds[selected[active.idx, k] == 0, k] = cntrd.o[active.idx][selected[active.idx, k] == 0]

      ## classify test sample and assess error status
      for(k in 1:v) {
        dd = distClanc(data = Y[active.idx, k], cntrds = cntrds, sd = p.sd[active.idx], prior = pi.k)

        if(match(min(dd), dd) != truth[k])
          cv.error[j, i, truth[k]] = cv.error[j, i, truth[k]] + 1
      }
    }
  }

  ## record numbers and proportions of errors
  for(i in 1:p) {
    if(d > 1)
      cv.err.cnt.cls[, i] = apply(cv.error[, , i], 1, sum)
    else
      cv.err.cnt.cls[, i] = sum(cv.error[, , i])
    cv.err.prpn.cls[, i] = cv.err.cnt.cls[, i] / nn[i]
  }
  cv.err.cnt.ttl = apply(cv.err.cnt.cls, 1, sum)
  cv.err.prpn.ttl = cv.err.cnt.ttl / n

  printClanc("\n")
  return(list("classErrors" = cv.err.prpn.cls, "overallErrors" = cv.err.prpn.ttl, "prior" = pi.k))
}

"buildClanc" <- function(data, id, cNames, train, active) {
  m = nrow(data)
  n = ncol(data)
  p = length(unique(id))
  nn = as.numeric(table(id))
  m.k = sqrt(1 / nn - 1 / n)
  geneNames = train$geneNames
  classNames = cNames
  prior = train$prior

  if(is.numeric(prior)) {
    if(length(prior) != p | sum(prior) != 1 | any(prior <= 0))
      stop("Invalid prior.")
    pi.k = prior
  } else {
    if(prior == "equal")
      pi.k = rep(1 / p, p)
    else if(prior == "class")
      pi.k = nn / n
  }

  ## select genes, update inactive centroid components
  selected = selectClanc(d.k = train$tStats, d.k.ord = list("d.k.rnks" = train$tStatsOrderIndices, 
    "d.k.o" = train$tStatsOrdered, "no.ties" = train$tStatsNoTies), active = active)
  active.idx = (1:m)[selected %*% rep(1, p) != 0]

  cntrds = train$classCntrds[active.idx, ]
  for(k in 1:p)
    cntrds[selected[active.idx, k] == 0, k] = train$overallCntrd[active.idx][selected[active.idx, k] == 0]

  ## prepare output
  out = list("geneNames" = geneNames[active.idx], "pooledSD" = train$pooledSD[active.idx], "cntrds" = cntrds, 
    "prior" = pi.k, "classNames" = classNames, "classID" = seq(p))

  ## assess training error
  error = testClanc(data = data, id = id, geneNames = geneNames, fit = out)
  out$classError = error / nn
  out$overallError = sum(error) / n

  return(out)
}

"predictClanc" <- function(data, geneNames, fit) {
  n = ncol(data)
  cntrds = fit$cntrds

  active.idx = match(fit$geneNames, geneNames)
  sd = fit$pooledSD
  prior = fit$prior

  if(any(is.na(active.idx)))
    stop("Gene names do not match those in classifier.")

  cls = rep(NA, n)
  for(i in 1:n) {
    dd = distClanc(data = data[active.idx, i], cntrds = cntrds, sd = sd, prior = prior)
    cls[i] = match(min(dd), dd)
  }

  return(cls)
}

"selectClanc" <- function(d.k, d.k.ord, active) {
  m = nrow(d.k)
  p = ncol(d.k)
  d.k.rnks = d.k.ord$d.k.rnks
  d.k.o = d.k.ord$d.k.o
  no.ties = d.k.ord$no.ties

  if(length(active) != p)
      active = rep(active, p)
  if(any(active >= m))
      selected = matrix(1, nrow = m, ncol = p)
  else {
    if(any((avail <- t(no.ties) %*% rep(1, m)) < active)) {
      active[avail < active] = avail[avail < active]
            
      warning("Not enough unique genes in at least one class.")
    }

    selected = matrix(0, nrow = m, ncol = p)
 
    jdx = vector("list", p)
    delta = rep(NA, p)
    for(i in 1:p) {
      jdx[[i]] = d.k.rnks[no.ties[, i], i][1:active[i]]
      delta[i] = d.k.o[no.ties[, i], i][active[i] + 1]
    }
    cls = rep(1:p, active)

    while(any((kdx <- table(unlist(jdx))) > 1)) {
      nms = as.numeric(names(kdx))
      for(j in as.numeric(names(kdx))[kdx > 1]) {
        ldx = cls[unlist(jdx) == j]

        ii = match(max(abs(d.k[j, ldx])), abs(d.k[j, ldx]))

        for(k in ldx[-ii]) {
          active[k] = active[k] + 1
          jdx[[k]] = c(jdx[[k]][-match(j, jdx[[k]])], d.k.rnks[no.ties[, k], k][active[k]])
          delta[k] = d.k.o[no.ties[, k], k][active[k] + 1]
        }
      }
    }

    for(i in 1:p)
      selected[jdx[[i]], i] = 1
  }
    
  return(selected)
}

"testClanc" <- function(data, id, geneNames, fit) {
  m = nrow(data)
  n = ncol(data)
  p = ncol(fit$cntrds)
  cntrds = fit$cntrds

  active.idx = match(fit$geneNames, geneNames)
  sd = fit$pooledSD
  prior = fit$prior

  if(any(is.na(active.idx)))
    stop("Gene names do not match those in classifier.")

  error = matrix(0, nrow = n, ncol = p)
  cls = predictClanc(data = data, geneNames = geneNames, fit = fit)

  for(i in 1:n)
    if(cls[i] != id[i])
      error[i, id[i]] = 1

  return(drop(apply(error, 2, sum)))
}

"trainClanc" <- function(data, id, geneNames, prior = "equal") {
  ## data: expression data.  (m x n) matrix of class numeric.
  ## id: class id's.  n vector of class numeric.
  ## geneNames: gene names.  m vector of class character.

  m = nrow(data)
  n = ncol(data)
  p = length(unique(id))
  nn = as.numeric(table(id))
  m.k = sqrt(1 / nn - 1 / n)

  ID = model.matrix(~ factor(id) - 1)

  ## pooled standard deviations
  p.sd = pooledSDClanc(data, ID)

  ## class and overall centroids
  cntrd.k = scale(data %*% ID, center = F, scale = nn)
  cntrd.o = drop(data %*% rep(1 / n, n))

  ## form statistics and order them
  d.k = scale((cntrd.k - cntrd.o) / p.sd, center = F, scale = m.k)
  d.k.ord = orderStatsClanc(d.k = d.k)

  out = list("pooledSD" = p.sd, "classCntrds" = cntrd.k, "overallCntrd" = cntrd.o, "tStats" = d.k, 
    "tStatsOrdered" = d.k.ord$d.k.o, "tStatsOrderIndices" = d.k.ord$d.k.rnks, "tStatsNoTies" = d.k.ord$no.ties, 
    "geneNames" = geneNames, "prior" = prior)

  return(out)
}

"orderStatsClanc" <- function(d.k) {
  m = nrow(d.k)
  p = ncol(d.k)

  d.k.rnks = apply(d.k, 2, function(x) { order(abs(x), decreasing = T) })
  d.k.o = matrix(d.k[as.numeric(d.k.rnks) + m * rep(0:(p - 1), each = m)], nrow = m)
  no.ties = apply(d.k.o, 2, function(x) { !duplicated(x) })

  return(list("d.k.rnks" = d.k.rnks, "d.k.o" = d.k.o, "no.ties" = no.ties))
}

"pooledSDClanc" <- function(X, ID) {
  m = nrow(X)
  n = ncol(X)
  p = ncol(ID)
  nn = drop(t(ID) %*% rep(1, n))

  mn = t(t(X %*% ID) / nn)
  dif2 = (X - mn %*% t(ID)) ^ 2
  psd = sqrt(drop(dif2 %*% rep(1 / (n - p), n)))

  return(psd)
}

"balancedFolds" <- function(y, nfolds = 5) {
  permuteRows = function(x) {
    dd = dim(x)
    n = dd[1]
    p = dd[2]
    mm = runif(length(x)) + rep(seq(n) * 10, rep(p, n))
    matrix(t(x)[order(mm)], n, p, byrow = T)
  }

  totals = table(y)
  fmax = max(totals)
  nfolds = min(nfolds, fmax)
  folds = as.list(seq(nfolds))
  yids = split(seq(y), y)
  bigmat = matrix(NA, ceiling(fmax / nfolds) * nfolds, length(totals))
  for(i in seq(totals))
    bigmat[seq(totals[i]), i] = sample(yids[[i]])
  smallmat = matrix(bigmat, nrow = nfolds)
  smallmat = permuteRows(t(smallmat))
  res = vector("list", nfolds)
  for(j in 1:nfolds) {
    jj = !is.na(smallmat[, j])
    res[[j]] = smallmat[jj, j]
  }

  return(res)
}

"distClanc" <- function(data, cntrds, sd, prior) {
  vv = sd ^ 2
  pi.k = prior

  if(length(vv) > 1)
    dd = drop(t(cntrds ^ 2) %*% (1 / vv)) - 2 * drop(t(data * cntrds) %*% (1 / vv)) - 2 * log(pi.k)
  else
    dd = drop(cntrds ^ 2 / vv) - 2 * drop(data * cntrds / vv) - 2 * log(pi.k)

  return(dd)
}
