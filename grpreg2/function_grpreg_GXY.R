setupG <- function(group, m = NULL, bilevel) {
  gf <- factor(group)
  if (any(levels(gf) == '0')) {
    g <- as.integer(gf) - 1
    lev <- levels(gf)[levels(gf) != '0']
  } else {
    g <- as.integer(gf)
    lev <- levels(gf)
  }
  if (is.numeric(group) | is.integer(group)) {
    lev <- paste0("G", lev)
  }
  # if (missing(m)) {
  if (is.null(m)) {
    m <- rep(NA, length(lev))
    names(m) <- lev
  } else {
    #if (all.equal(sort(names(m)), sort(group)))
    TRY <- try(as.integer(group) == g)
    if (inherits(TRY, 'try-error') ||
        any(!TRY))
      stop('Attempting to set group.multiplier is ambiguous if group is not a factor',
           call. = FALSE)
    if (length(m) != length(lev))
      stop("Length of group.multiplier must equal number of penalized groups",
           call. = FALSE)
    if (storage.mode(m) != "double")
      storage.mode(m) <- "double"
    if (any(m < 0))
      stop('group.multiplier cannot be negative', call. = FALSE)
  }
  structure(g, levels = lev, m = m)
}
######################################################################################
reorderG <- function(g, m, bilevel) {
  og <- g
  lev <- attr(g, 'levels')
  m <- attr(g, 'm')
  if (any(g == 0)) {
    g <- as.integer(relevel(factor(g), "0")) - 1
  }
  if (any(order(g) != 1:length(g))) {
    reorder <- TRUE
    gf <- factor(g)
    if (any(levels(gf) == "0")) {
      gf <- relevel(gf, "0")
      g <- as.integer(gf) - 1
    } else {
      g <- as.integer(gf)
    }
    ord <- order(g)
    ord.inv <- match(1:length(g), ord)
    g <- g[ord]
  } else {
    reorder <- FALSE
    ord <- ord.inv <- NULL
  }
  structure(
    g,
    levels = lev,
    m = m,
    ord = ord,
    ord.inv = ord.inv,
    reorder = reorder
  )
}
##################################################################################################################
newXG <- function(X, g, m = NULL, ncolY, bilevel) {
  # Coerce X to matrix
  if (!inherits(X, "matrix")) {
    tmp <- try(X <- model.matrix( ~ 0 + ., data = X), silent = TRUE)
    if (inherits(tmp, "try-error"))
      stop("X must be a matrix or able to be coerced to a matrix", call. = FALSE)
  }
  if (storage.mode(X) == "integer")
    storage.mode(X) <- "double"
  if (any(is.na(X)))
    stop(
      "Missing data (NA's) detected in X.  You must eliminate missing data (e.g., by removing cases, removing features, or imputation) before passing X to grpreg",
      call. = FALSE
    )
  if (length(g) != ncol(X))
    stop ("Dimensions of group is not compatible with X", call. = FALSE)
  xnames <-
    if (is.null(colnames(X)))
      paste("V", 1:ncol(X), sep = "")
  else
    colnames(X)
  
  # Setup group
  G <- setupG(g, m, bilevel)
  
  # Reconfigure for multiple outcomes, if necessary
  if (ncolY > 1) {
    X <- multiX(X, ncolY)
    G <- multiG(G, ncolY)
  }
  
  # Feature-level standardization
  std <- .Call("standardize", X)
  XX <- std[[1]]
  center <- std[[2]]
  scale <- std[[3]]
  ##
  XX = X
  center = center * 0
  scale = scale * 0 + 1
  
  nz <- which(scale > 1e-6)                # non-constant columns
  if (length(nz) != ncol(X)) {
    XX <- XX[, nz, drop = FALSE]
    G <- subsetG(G, nz)
  }
  
  # Reorder groups, if necessary
  G <- reorderG(G, attr(G, 'm'), bilevel)
  if (attr(G, 'reorder'))
    XX <- XX[, attr(G, 'ord')]
  
  # Group-level standardization
  if (!bilevel) {
    XX <- orthogonalize(XX, G)
    g <- attr(XX, "group")
  } else {
    g <- as.integer(G)
  }
  
  # Set group multiplier if missing
  m <- attr(G, 'm')
  if (all(is.na(m))) {
    m <- if (bilevel)
      rep(1, max(g))
    else
      sqrt(table(g[g != 0]))
  }
  
  # Return
  return(
    list(
      X = XX,
      g = g,
      m = m,
      reorder = attr(G, 'reorder'),
      ord.inv = attr(G, 'ord.inv'),
      names = xnames,
      center = center,
      scale = scale,
      nz = nz
    )
  )
}
###################################################################################################
newY <- function(y, family) {
  if (is.data.frame(y))
    y <- as.matrix(y)
  if (is.matrix(y)) {
    d <- dim(y)
    y <- t(y)
  } else {
    d <- c(length(y), 1)
  }
  
  # Convert fuzzy binomial data
  if (family == "binomial" && typeof(y) != "logical") {
    tab <- table(y)
    if (length(tab) > 2)
      stop("Attemping to use family='binomial' with non-binary data",
           call. = FALSE)
    if (!identical(names(tab), c("0", "1"))) {
      message(paste0("Logistic regression modeling Pr(y=", names(tab)[2], ")"))
      y <- as.double(as.character(y) == names(tab)[2])
      if (d[2] > 1)
        attr(y, "dim") <- d
    }
  }
  
  # Convert to double, if necessary
  if (typeof(y) != "double") {
    tryCatch(
      storage.mode(y) <-
        "double",
      warning = function(w) {
        stop("y must be numeric or able to be coerced to numeric", call. = FALSE)
      }
    )
  }
  if (any(is.na(y)))
    stop(
      "Missing data (NA's) detected in outcome y.  You must eliminate missing data (e.g., by removing cases or imputation) before passing y to grpreg",
      call. = FALSE
    )
  
  # Handle multi
  if (is.matrix(y)) {
    if (ncol(y) > 1) {
      if (is.null(colnames(y)))
        paste("Y", 1:ncol(y), sep = "")
    }
    attributes(y) <- NULL
  }
  
  if (family == "gaussian") {
    meanY <- mean(y) * 0
    # y <- y - meanY
    attr(y, "mean") <- meanY
  }
  attr(y, "m") <- d[2]
  y
}
