
# this is to keep the old S3 behaviour inside the NAMESPACE.
is = function(x, class) {

  if (identical(class, "double"))
    is.double(x)
  else if ("double" %in% class)
    any(class(x) %in% class) || is.double(x)
  else
    any(class(x) %in% class)

}#IS

# make bnlearn's classes known to S4.
setClass("bn")
setClass("bn.fit")
setClass("bn.naive")
setClass("bn.tan")

# return the nodes in the graph.
.nodes = function(x) {

  # check x's class.
  check.bn.or.fit(x)

  if (inherits(x, "bn"))
    names(x$nodes)
  else
    names(x)

}#.NODES


# if no generic is present, create it.
if (!isGeneric("nodes"))
  setGeneric("nodes", function(object, ...) standardGeneric("nodes"))

setMethod("nodes", "bn", function(object) .nodes(object))
setMethod("nodes", "bn.fit", function(object) .nodes(object))
setMethod("nodes", "bn.naive", function(object) .nodes(object))
setMethod("nodes", "bn.tan", function(object) .nodes(object))

# get the degree of a node.
.degree = function(x, node) {

  # check x's class.
  check.bn.or.fit(x)
  # a valid node is needed.
  check.nodes(nodes = node, graph = x, max.nodes = 1)

  if (is(x, "bn"))
    length(x$nodes[[node]]$nbr)
  else
    length(x[[node]]$parents) + length(x[[node]]$children)

}#.DEGREE

# if no generic is present, create it.
if (!isGeneric("degree"))
  setGeneric("degree", function(object, Nodes, ...) standardGeneric("degree"))

setMethod("degree", "bn", function(object, Nodes) .degree(object, Nodes))
setMethod("degree", "bn.fit", function(object, Nodes) .degree(object, Nodes))
setMethod("degree", "bn.naive", function(object, Nodes) .degree(object, Nodes))
setMethod("degree", "bn.tan", function(object, Nodes) .degree(object, Nodes))

