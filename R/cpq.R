
# backend for conditional probability queries.
conditional.probability.query = function(fitted, event, evidence, method,
    extra, probability = TRUE, cluster = NULL, debug = FALSE) {

  # consider only the upper closure of event and evidence to reduce the number
  # of variables in the Monte Carlo simulation.
  fitted = reduce.fitted(fitted = fitted, event = event, evidence = evidence,
             nodes = extra$query.nodes, method = method, debug = debug)

  if (method == "ls") {

    if (!is.null(cluster)) {

      # get the number of slaves.
      s = nSlaves(cluster)
      # divide the number of particles among the slaves.
      batch = n = ceiling(extra$n / s)

      if (probability) {

        results = parSapply(cluster, seq(s),
          function(x) {

            logic.sampling(fitted = fitted, event = event,
              evidence = evidence, n = n, batch = batch, debug = debug)

          })

        return(mean(results))

      }#THEN
      else {

        results = parLapply(cluster, seq(s),
          function(x) {

            logic.distribution(fitted = fitted, nodes = event,
              evidence = evidence, n = n, batch = batch, debug = debug)

          })

        return(do.call(rbind, results))

      }#ELSE

    }#THEN
    else {

      if (probability) {

        logic.sampling(fitted = fitted, event = event, evidence = evidence,
          n = extra$n, batch = extra$batch, debug = debug)

      }#THEN
      else {

        logic.distribution(fitted = fitted, nodes = event, evidence = evidence,
          n = extra$n, batch = extra$batch, debug = debug)

      }#ELSE

    }#ELSE

  }#THEN
  else if (method == "lw") {

    if (!is.null(cluster)) {

      # get the number of slaves.
      s = nSlaves(cluster)
      # divide the number of particles among the slaves.
      batch = n = ceiling(extra$n / s)

      if (probability) {

        results = parSapply(cluster, seq(s),
          function(x) {

            weighting.sampling(fitted = fitted, event = event,
              evidence = evidence, n = n, batch = batch, debug = debug)

          })

        return(mean(results))

      }#THEN
      else {

        results = parLapply(cluster, seq(s),
          function(x) {

            weighting.distribution(fitted = fitted, nodes = event,
              evidence = evidence, n = n, batch = batch, debug = debug)

          })

        return(do.call(rbind, results))

      }#ELSE

    }#THEN
    else {

      if (probability) {

        weighting.sampling(fitted = fitted, event = event, evidence = evidence,
          n = extra$n, batch = extra$batch, debug = debug)

      }#THEN
      else {

        weighting.distribution(fitted = fitted, nodes = event,
          evidence = evidence, n = extra$n, batch = extra$batch, debug = debug)

      }#ELSE

    }#ELSE

  }#THEN
  else if (method == "exact") {

    if (!is.fitted.discrete(fitted))
      stop("exact inference does only support discrete bns.")

    if (probability) {

      exact.prob(fitted = fitted, event = event, evidence = evidence, debug = debug)

    }#THEN
    else {

      exact.dist(fitted = fitted, event = event, evidence = evidence, debug = debug)

    }#ELSE

  }#THEN

}#CONDITIONAL.PROBABILITY.QUERY

# create an empty data frame from a bn.fit object.
fit.dummy.df = function(fitted, nodes) {

  dummy = sapply(nodes, function(x) {

    node = fitted[[x]]

    if (is(node, "bn.fit.dnode"))
      return(factor(character(0), levels = dimnames(node$prob)[[1]]))
    else if (is(node, "bn.fit.onode"))
      return(ordered(character(0), levels = dimnames(node$prob)[[1]]))
    else if (is(node, "bn.fit.gnode"))
      return(numeric(0))

  })

  return(minimal.data.frame(dummy))

}#FIT.DUMMY.DF

# reduce a bn.fit object to the upper closure of event and evidence nodes.
reduce.fitted = function(fitted, event, evidence, nodes, method, debug) {

  if (is.null(nodes)) {

    # find out which nodes are involved in the event and the evidence and
    # construct their upper closure.
    nodes = names(fitted)
    nodes.event = nodes[nodes %in% explode(event)]
    nodes.evidence = nodes[nodes %in% explode(evidence)]
    upper.closure = schedule(fitted, start = union(nodes.evidence, nodes.event),
                      reverse = TRUE)

    if (debug) {

      cat("* checking which nodes are needed.\n")
      cat("  > event involves the following nodes:", nodes.event, "\n")
      cat("  > evidence involves the following nodes:", nodes.evidence, "\n")
      cat("  > upper closure is '", upper.closure, "'\n")

    }#THEN

  }#THEN
  else {

    # construct the upper closure of the query nodes.
    upper.closure = schedule(fitted, start = nodes, reverse = TRUE)

    if (debug) {

      cat("* using specified query nodes.\n")
      cat("  > upper closure is '", upper.closure, "'\n")

    }#THEN

  }#ELSE

  # check whether the upper closure is correct: tricky expressions are not
  # always handled correctly by explode().
  dummy = fit.dummy.df(fitted, upper.closure)
  # testing evidence (TRUE or an expression in LS, a list in LW).
  if (!(is.language(evidence) || identical(evidence, TRUE)))
    try.evidence = TRUE
  else
    try.evidence = try(eval(evidence, dummy), silent = TRUE)
  # testing event (it's the label nodes in cpdist; TRUE or an expression 
  # in cpquery).
  if (!(is.language(event) || identical(event, TRUE)))
    try.event = TRUE
  else
    try.event = try(eval(event, dummy), silent = TRUE)

  # create the subgraph corresponding to the upper closure.
  if (is.logical(try.event) && is.logical(try.evidence)) {

    if (debug)
      cat("  > generating observations from", length(upper.closure), "/", 
        length(fitted), "nodes.\n")

    fitted = fitted[upper.closure]
    class(fitted) = "bn.fit"

  }#THEN
  else {

    if (debug)
      cat("  > unable use the upper closure, using the whole network.\n")

  }#ELSE

  return(fitted)

}#REDUCE.FITTED

# compute conditional probabilities with forward/logic sampling.
logic.sampling = function(fitted, event, evidence, n, batch, debug = FALSE) {

  cpxe = cpe = 0L
  filtered = logical(n)
  matching = logical(n)
  r = logical(n)

  # count how many complete batches we have to generate.
  nbatches = n %/% batch
  # count how many observations are in the last one.
  last.one = n %% batch

  for (m in c(rep(batch, nbatches), last.one)) {

    # do a hard reset of generated.data, so that the memory used by the data
    # set generated in the previous iteration can be garbage-collected and
    # reused _before_ rbn() returns.
    generated.data = NULL

    # generate random data from the bayesian network.
    if (m > 0) {

      generated.data = rbn.backend(x = fitted, n = m)

    }#THEN
    else
      break

    if (debug)
      cat("* generated", m, "samples from the bayesian network.\n")

    # evaluate the expression defining the evidence.
    if (identical(evidence, TRUE))
      r = rep(TRUE, m)
    else
      r = eval(evidence, generated.data, parent.frame())
    # double check that this is a logical vector.
    if (!is.logical(r))
      stop("evidence must evaluate to a logical vector.")
    # double check that it is of the right length.
    if (length(r) != m)
      stop("logical vector for evidence is of length ", length(r),
        " instead of ", m, ".")
    # filter out the samples not matching the evidence we assume.
    filtered = r & !is.na(r)

    # update the global counters.
    cpe = cpe + length(which(filtered))

    if (debug) {

      lwfilter = length(which(filtered))
      if (!identical(evidence, TRUE))
        cat("  > evidence matches ", lwfilter, " samples out of ", m,
          " (p = ", lwfilter/m, ").\n", sep = "")
      else
        cat("  > evidence matches ", m, " samples out of ", m,
          " (p = 1).\n", sep = "")

    }#THEN

    # evaluate the expression defining the event.
    if (identical(event, TRUE))
      r = rep(TRUE, m)
    else
      r = eval(event, generated.data, parent.frame())
    # double check that this is a logical vector.
    if (!is.logical(r))
      stop("event must evaluate to a logical vector.")
    # double check that it is of the right length.
    if (length(r) != m)
      stop("logical vector for event is of length ", length(r),
        " instead of ", m, ".")
    # filter out the samples not matching the event we are looking for.
    matching = filtered & r & !is.na(r)

    # update the global counters.
    cpxe = cpxe + length(which(matching))

    if (debug) {

      lwmatch = length(which(matching))
      lwratio = ifelse(lwfilter == 0, 0, lwmatch/lwfilter)
      if (!identical(event, TRUE))
        cat("  > event matches ", lwmatch, " samples out of ", lwfilter,
          " (p = ", lwratio, ").\n", sep = "")
      else
        cat("  > event matches ", lwfilter, " samples out of ", lwfilter,
          " (p = 1).\n", sep = "")

    }#THEN

  }#FOR

  # prevent divide-by-zero errors.
  result = ifelse(cpe == 0, 0, cpxe / cpe)

  if (debug && (nbatches > 1)) {

    cat("* generated a grand total of", n, "samples.\n")
    cat("  > event matches ", cpxe, " samples out of ", cpe,
      " (p = ", result, ").\n", sep = "")

  }#THEN

  return(result)

}#LOGIC.SAMPLING

# generate random observations from conditional distributions with forward/logic
# sampling.
logic.distribution = function(fitted, nodes, evidence, n, batch, debug = FALSE) {

  filtered = logical(n)
  result = NULL

  # count how many complete batches we have to generate.
  nbatches = n %/% batch
  # count how many observations are in the last one.
  last.one = n %% batch

  for (m in c(rep(batch, nbatches), last.one)) {

    # do a hard reset of generated.data, so that the memory used by the data
    # set generated in the previous iteration can be garbage-collected and
    # reused _before_ rbn() returns.
    generated.data = NULL

    # generate random data from the bayesian network.
    if (m > 0)
      generated.data = rbn.backend(x = fitted, n = m)
    else
      break

    if (debug)
      cat("* generated", m, "samples from the bayesian network.\n")

    # evaluate the expression defining the evidence.
    r = eval(evidence, generated.data, parent.frame())
    # double check that this is a logical vector.
    if (!is.logical(r))
      stop("evidence must evaluate to a logical vector.")
    # double check that it is of the right length.
    if ((length(r) != 1) && (length(r) != m))
      stop("logical vector for evidence is of length ", length(r),
        " instead of ", m, ".")
    # filter out the samples not matching the evidence we assume.
    filtered = r & !is.na(r)

    if (debug) {

      lwfilter = length(which(filtered))
      if (!identical(evidence, TRUE))
        cat("  > evidence matches ", lwfilter, " samples out of ", m,
          " (p = ", lwfilter/m, ").\n", sep = "")
      else
        cat("  > evidence matches ", m, " samples out of ", m,
          " (p = 1).\n", sep = "")

    }#THEN

    # update the return value.
    result = rbind(result, generated.data[filtered, nodes, drop = FALSE])

  }#FOR

  # reset the row names.
  rownames(result) = NULL

  if (debug && (nbatches > 1)) 
    cat("* generated a grand total of", n, "samples.\n")
 
  return(result)

}#LOGIC.DISTRIBUTION

# compute conditional probabilities with likelihood weighting.
weighting.sampling = function(fitted, event, evidence, n, batch, debug = FALSE) {

  cpxe = cpe = 0
  matching = logical(n)
  r = logical(n)

  # count how many complete batches we have to generate.
  nbatches = n %/% batch
  # count how many observations are in the last one.
  last.one = n %% batch

  weights = function(data) {

    if (isTRUE(evidence))
      return(rep(1, nrow(data)))
    else 
      exp(.Call("entropy_loss", 
                fitted = fitted[names(evidence)],
                data = data,
                by.sample = TRUE, 
                debug = FALSE, 
                PACKAGE = "bnlearn"))

  }#WEIGHTS

  for (m in c(rep(batch, nbatches), last.one)) {

    # do a hard reset of generated.data, so that the memory used by the data
    # set generated in the previous iteration can be garbage-collected and
    # reused _before_ rbn() returns.
    generated.data = NULL

    # generate random data from the bayesian network.
    if (m > 0)
      generated.data = rbn.backend(x = fitted, fix = evidence, n = m)
    else
      break

    if (debug)
      cat("* generated", m, "samples from the bayesian network.\n")

    # evaluate the expression defining the event.
    if (identical(event, TRUE))
      r = rep(TRUE, m)
    else
      r = eval(event, generated.data, parent.frame())
    # double check that this is a logical vector.
    if (!is.logical(r))
      stop("event must evaluate to a logical vector.")
    # double check that it is of the right length.
    if (length(r) != m)
      stop("logical vector for event is of length ", length(r),
        " instead of ", m, ".")
    # filter out the samples not matching the event we are looking for.
    matching = r & !is.na(r)

    # compute the probabilities and use them as weigths.
    w = weights(generated.data) 
    cpe = cpe + sum(w)
    cpxe = cpxe + sum(w[matching]) 

    if (debug)
      cat("  > event has a probability mass of ", cpxe, " out of ", cpe, ".\n", sep = "")

  }#FOR

  # compute the conditional probability.
  result = cpxe / cpe

  if (debug && (nbatches > 1)) {

    cat("* generated a grand total of", n, "samples.\n")
    cat("  > event has a probability mass of ", cpxe, " out of ", cpe, 
        " (p = ", result, ").\n", sep = "")

  }#THEN

  return(result)

}#WEIGHTING.SAMPLING

# generate random observations from conditional distributions with likelihood
# weighting.
weighting.distribution = function(fitted, nodes, evidence, n, batch, debug = FALSE) {

  result = NULL

  # count how many complete batches we have to generate.
  nbatches = n %/% batch
  # count how many observations are in the last one.
  last.one = n %% batch

  for (m in c(rep(batch, nbatches), last.one)) {

    # do a hard reset of generated.data, so that the memory used by the data
    # set generated in the previous iteration can be garbage-collected and
    # reused _before_ rbn() returns.
    generated.data = NULL

    # generate random data from the bayesian network.
    if (m > 0)
      generated.data = rbn.backend(x = fitted, fix = evidence, n = m)
    else
      break

    if (debug)
      cat("* generated", m, "samples from the bayesian network.\n")

    # update the return value.
    result = rbind(result, generated.data[, nodes, drop = FALSE])

  }#FOR

  if (debug && (nbatches > 1)) 
    cat("* generated a grand total of", n, "samples.\n")
 
  return(result)

}#WEIGHTING.DISTRIBUTION

exact.cpt = function(fitted, event, evidence, debug = FALSE){

  extract.names = function(call) {
    names = character(0)
    for (x in as.list(call)[-1]) {
      if (is.call(x)) {
        names = c(names, extract.names(x))
      }#THEN
      else if (is.name(x)) {
        names = c(names, as.character(x))
      }#ELSE
    }#FOR
    return(names)
  }#EXTRACT.NAMES

  if(is.call(event) | is.logical(event)) {
    nodes = extract.names(event)
  }#THEN
  else {
    nodes = event
  }#ELSE

  if(is.call(evidence) | is.logical(evidence)) {
    nodes = union(nodes, extract.names(evidence))
  }#THEN
  else {
    nodes = union(nodes, evidence)
  }#ELSE

  nodes = intersect(nodes, names(fitted))

  # Joint distribution of the target and conditional nodes knowing their
  # parents, their parent's parents, etc.
  to.check = nodes
  while(length(to.check) > 0) {
    for (node in to.check) {
      parents = fitted[[node]]$parents
      to.check = setdiff(to.check, node)
      to.check = union(to.check, setdiff(parents, nodes))
      nodes = union(nodes, parents)
    }#FOR
  }#WHILE

  nbnodes = length(nodes)

  cpt.table = expand.grid(lapply(nodes, function(x) {
    dimnames(fitted[[x]]$prob)[[1]]
  }))
  names(cpt.table) = nodes

  m = nrow(cpt.table)

  # evaluate the expression defining the evidence.
  if (identical(evidence, TRUE))
    r = rep(TRUE, m)
  else
    r = eval(evidence, cpt.table, parent.frame())

  # double check that it is a logical vector.
  if (!is.logical(r))
    stop("evidence must evaluate to a logical vector.")
  # double check that it is of the right length.
  if (length(r) != m)
    stop("logical vector for evidence is of length ", length(r), " instead of ", m, ".")

  cpt.table = cpt.table[r, , drop=FALSE]

  for (node in nodes) {
    cpt = as.data.frame(fitted[[node]]$prob, )
    names(cpt)[1] = node # FIX: with 0 parents variable name is missing
    names(cpt)[ncol(cpt)] = paste("p", node, sep=".")
    cpt.table = merge(cpt.table, cpt, by=names(cpt)[-length(cpt)])
  }#FOR

  joint.table = cbind(cpt.table[, nodes, drop=FALSE], "p"=1)
  for(i in 1:nbnodes) {
    joint.table[, nbnodes + 1] = joint.table[, nbnodes + 1] * cpt.table[, nbnodes + i]
  }#FOR

  return(joint.table)

}#EXACT.CPT

exact.prob = function(fitted, event, evidence, debug = FALSE){

  cpt = exact.cpt(fitted, event, evidence, debug)

  # evaluate the expression defining the event.
  if (identical(event, TRUE))
    r = rep(TRUE, nrow(cpt))
  else
    r = eval(event, cpt, parent.frame())

  # double check that this is a logical vector.
  if (!is.logical(r))
    stop("event must evaluate to a logical vector.")
  # double check that it is of the right length.
  if (length(r) != nrow(cpt))
    stop("logical vector for event is of length ", length(r), " instead of ", m, ".")

  p.AB = sum(cpt[r, ncol(cpt)])
  p.B = sum(cpt[, ncol(cpt)])

  if(p.AB == 0)
    return(0)

  return(p.AB/p.B)

}#EXACT.PROB

exact.dist = function(fitted, event, evidence, debug = FALSE) {

  probs = exact.cpt(fitted, event, evidence, debug)

  for(col in 1:(ncol(probs)-1)) {
    probs[, col] = factor(probs[, col], exclude=NULL)
    probs = probs[order(probs[, col]), ]
  }#FOR
  
  p.col = ncol(probs)
  cpt = table(probs[, -p.col], dnn=names(probs)[-p.col])
  cpt[1:length(cpt)] = probs[, p.col]
  cpt = prop.table(margin.table(cpt, 1:length(event)))

  return(cpt)

}#EXACT.DIST
