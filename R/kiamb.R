
kiamb.global = function(x, whitelist, blacklist, test, alpha, nbr.join, k,
                        test.args, strict, debug = FALSE) {
  
  nodes = names(x)
  
  # 1. [Compute Markov Blankets]
  mb = lapply(as.list(nodes), kiamb, data = x, whitelist = whitelist, blacklist = blacklist,
              nodes = nodes, alpha = alpha, test.args = test.args, k = k, test = test, debug = debug)
  names(mb) = nodes
  
  # check markov blankets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, mb = TRUE, debug = debug, filter = nbr.join)
  
  # 2. [Compute Graph Structure]
  mb = lapply(as.list(nodes), neighbour, mb = mb, data = x, alpha = alpha,
              test.args = test.args, test = test, debug = debug, whitelist = whitelist, blacklist = blacklist)
  names(mb) = nodes
  
  # check neighbourhood sets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, debug = debug, filter = nbr.join)
  
  return(mb)
  
}#KIAMB.GLOBAL

kiamb = function(x, data, nodes, alpha, test.args, k,
                 whitelist, blacklist, backtracking = NULL,
                 start = character(0), test, debug = FALSE) {
  
  nodes = setdiff(nodes, x)
  mb = start
  known.good = character(0)
  known.bad = character(0)
  
  if (debug) {
    
    cat("----------------------------------------------------------------\n")
    cat("* learning the markov blanket of", x, ".\n")
    if (length(start) > 0)
      cat("* initial set includes '", mb, "'.\n")
    
  }#THEN
  
  # whitelisted nodes are included by default (if there's a direct arc
  # between them of course they are in each other's markov blanket).
  # arc direction is irrelevant here
  known.good = nodes[sapply(nodes, function(y) {
    is.whitelisted(whitelist, c(x, y), either = TRUE) })]
  
  # vice versa with blacklisted nodes
  known.bad = nodes[sapply(nodes, function(y) {
    is.blacklisted(blacklist, c(x, y), either = TRUE) })]
  
  # use backtracking for a further screening of the nodes to be checked.
  # Nodes whose markov blanket includes this node are included, because
  # X \in MB(Y) <=> Y \in MB(X) and vice versa X \not\in MB(Y) <=> Y \not\in MB(X)
  if (!is.null(backtracking)) {
    
    known.good = union(known.good, names(backtracking[backtracking]))
    known.bad = union(known.bad, names(backtracking[!backtracking]))
    
    if (debug) {
      
      cat("    * known good (backtracking): '", known.good, "'.\n")
      cat("    * known bad (backtracking): '", known.bad, "'.\n")
      cat("    * nodes still to be tested for inclusion: '", setdiff(nodes, names(backtracking)), "'.\n")
      
    }#THEN
    
  }#THEN
  
  mb = union(mb, known.good)
  n.in.mb = nodes %in% mb
  
  # Add candidates nodes to Markov Boundary
  repeat {
    
    mb = nodes[n.in.mb]
    association = rep(1, length(nodes))
    
    # known bad nodes are not checked for inclusion
    can.mb = which((!n.in.mb) & (! nodes %in% known.bad))
    
    # randomisation
    if(length(can.mb) > 1) {
      can.mb = sample(can.mb)
    }
    
    for (i in can.mb) {
      
      y = nodes[i]
      association[i] = conditional.test(x, y, sx = setdiff(mb, y), test = test,
                                        data = data, test.args = test.args, alpha = alpha)
      
      # k=0 : early stop at the first dependent node
      if (k == 0 & association[i] <= alpha) {
        can.mb = i
        break
      }
    }
    association = association[can.mb]
    
    # filter candidates for inclusion
    s = association <= alpha
    can.mb = can.mb[s]
    association = association[s]
    
    # early stop if no candidates
    if (length(can.mb) == 0)
      break

    # draw a random (k-driven) subset
    s = 1:max(1, as.integer(length(can.mb)*k))
    can.mb = can.mb[s]
    association = association[s]
    
    # add the most dependent node
    best = which.min(association)
    n.in.mb[can.mb[best]] = TRUE
    
    if (debug) {
      
      cat("    * added: '", nodes[can.mb[best]], "'\n")
      
    }#THEN
    
  }#REPEAT
  
  # Remove wrongly added nodes
  repeat {
    
    mb = nodes[n.in.mb]
    association = rep(0, length(nodes))
    
    # known good nodes are not checked for exclusion
    can.not.mb = which(n.in.mb & ! nodes %in% known.good)
    
    # randomisation
    if(length(can.not.mb) > 1) {
      can.mb = sample(can.not.mb)
    }
    
    for (i in can.not.mb) {
      
      y = nodes[i]
      association[i] = conditional.test(x, y, sx = setdiff(mb, y), test = test,
                                        data = data, test.args = test.args, alpha = alpha)
      
      # k=0 : early stop at the first independent node
      if (k == 0 & association[i] > alpha) {
        can.not.mb = i
        break
      }
    }
    association = association[can.not.mb]
    
    # filter candidates for exclusion
    s = association > alpha
    can.not.mb = can.not.mb[s]
    association = association[s]
    
    # early stop if no candidates
    if (length(can.not.mb) == 0)
      break
    
    # draw a random (k-driven) subset
    s = 1:max(1, as.integer(length(can.not.mb)*k))
    can.not.mb = can.not.mb[s]
    association = association[s]
    
    # remove the least dependent node
    worst = which.max(association)
    n.in.mb[can.not.mb[worst]] = FALSE
    
    if (debug) {
      
      cat("    * removed: '", nodes[can.mb[worst]], "'\n")
      
    }#THEN
    
  }#REPEAT
  
  mb = nodes[n.in.mb]
  return(mb)
  
}#KIAMB
