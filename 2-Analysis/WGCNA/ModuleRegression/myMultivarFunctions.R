my.drop1.mlm <- function (object, scope = NULL,
                       test = c("Wilks", "Pillai", "Hotelling-Lawley", "Roy"),
                       total=TRUE, add=FALSE, ...)
{
  ## add=TRUE   do add instead of drop1
  Pillai <- function(eig, q, df.res) {
    test <- sum(eig/(1 + eig))
    p <- length(eig)
    s <- min(p, q)
    n <- 0.5 * (df.res - p - 1)
    m <- 0.5 * (abs(p - q) - 1)
    tmp1 <- 2 * m + s + 1
    tmp2 <- 2 * n + s + 1
    c(test, (tmp2/tmp1 * test)/(s - test), s * tmp1, s *
        tmp2)
  }
  Wilks <- function(eig, q, df.res) {
    test <- prod(1/(1 + eig))
    p <- length(eig)
    tmp1 <- df.res - 0.5 * (p - q + 1)
    tmp2 <- (p * q - 2)/4
    tmp3 <- p^2 + q^2 - 5
    tmp3 <- if (tmp3 > 0)
      sqrt(((p * q)^2 - 4)/tmp3)
    else 1
    c(test, ((test^(-1/tmp3) - 1) * (tmp1 * tmp3 - 2 * tmp2))/p/q,
      p * q, tmp1 * tmp3 - 2 * tmp2)
  }
  HL <- function(eig, q, df.res) {
    test <- sum(eig)
    p <- length(eig)
    m <- 0.5 * (abs(p - q) - 1)
    n <- 0.5 * (df.res - p - 1)
    s <- min(p, q)
    tmp1 <- 2 * m + s + 1
    tmp2 <- 2 * (s * n + 1)
    c(test, (tmp2 * test)/s/s/tmp1, s * tmp1, tmp2)
  }
  Roy <- function(eig, q, df.res) {
    p <- length(eig)
    test <- max(eig)
    tmp1 <- max(p, q)
    tmp2 <- df.res - tmp1 + q
    c(test, (tmp2 * test)/tmp1, tmp1, tmp2)
  }
  ##-     if (!is.null(object$drop1)) return (object$drop1)
  if (!(inherits(object, "maov") || inherits(object, "mlm")))
    stop("object must be of class \"maov\" or \"mlm\"")
  test <- match.arg(test)
  asgn <- object$assign[object$qr$pivot[1:object$rank]]
  tl <- attr(object$terms, "term.labels")
  ## scope
  if (is.null(scope))
    scope <- if (add) attr(terms(update.formula(object, ~(.)^2)),
                           "term.labels") else
                             drop.scope(object)
  else {
    if (!is.character(scope))
      scope <- attr(terms(update.formula(object, scope)),
                    "term.labels")
    ##-         if (!all(match(scope, tl, FALSE)))
    if (!(add||all(match(scope, tl, FALSE))))
      stop("!drop1.mlm! scope is not a subset of term labels")
  }
  ns <- length(scope)
  rdf <- object$df.residual
  res <- object$resid
  ## ::: needed for finding the data later
  ldata <- eval(object$call$data,
                envir=environment(formula(object)))
  ladd <- 1
  if (add) {
    ladd <- -1
    lna <- i.add1na(object, scope)
    if (!is.null(lna)) res[lna,1] <- NA
  }
  res <- nainf.exclude(res)
  lna <- attr(res,"na.action")
  if (!is.null(lna)) ldata <- ldata[-lna,]
  ## full model
  rss <- crossprod(as.matrix(res))
  rss.qr <- qr(rss)
  if (rss.qr$rank < NCOL(res))
    stop(paste("!drop1.mlm! residuals have rank", rss.qr$rank, "<",
               ncol(res)))
  stats <- matrix(NA,length(scope),4)
  dimnames(stats) <- list(scope,c(test,"F.stat","dfnum","dfden"))
  tstfn <- switch(test, Pillai = Pillai, Wilks = Wilks,
                  "Hotelling-Lawley" = HL, HL = HL, Roy = Roy)
  object$call[[1]] <- as.name("lm")
  ## loop through scope
  for (lsc in scope) {
    lfo <- as.formula(paste(if (add) "~.+" else "~.-",lsc))
    lrg <- update(object, lfo, data=ldata, model=FALSE) # ,data=data
    dfj <- ladd * (lrg$df.residual - rdf)
    bss <- ladd * (crossprod(lrg$resid)-rss)
    eigs <- Re(eigen(qr.coef(rss.qr, bss),symmetric = FALSE)$values)
    stats[lsc,] <- tstfn(eigs, dfj, rdf)
  }
  ldf <- stats[1,3:4]
  names(ldf) <- c("numerator","denominator")
  if (total) {
    lpr <- predict(object)
    if (length(lna)) lpr <- lpr[-lna,] # drop rows with NA
    yy <- scale(res + lpr, scale=FALSE)
    bss <- crossprod(yy)-rss
    eigs <- Re(eigen(qr.coef(rss.qr, bss), symmetric = FALSE)$values)
    stats <- rbind(stats, "<total>"= tstfn(eigs, object$df[1], rdf))
  }
  data.frame(stats,
             p.value = pf(stats[,2],stats[,3],stats[,4], lower.tail = FALSE))
  ##    attr(stats,"df") <- ldf
  ##    stats
} ## {drop1.mlm}

my.add1.mlm <- #F
  function (object, scope=NULL,
            test = c("Wilks", "Pillai", "Hotelling-Lawley", "Roy"), ...)
  {
    ## Purpose:    add1 for regr objects
    ## ----------------------------------------------------------------------
    my.drop1.mlm(object, scope=scope, test=test, total=FALSE, add=TRUE, ...)
  }


# object=fit.mv1
# scope=NULL
# expand=FALSE
# scale = 0
# direction ="both"
# trace = TRUE
# keep = NULL
# steps = 1000
# k = 2



my.step.regr <- function (object, scope=NULL, expand=FALSE, scale = 0,
                       direction = c("both", "backward", "forward"), 
                       trace = FALSE, keep = NULL,
                       steps = 1000, k = 2, ...)
{
  ## ----------------------------------------------------------------------
  ## Author: Werner Stahel, Date: 15 May 2012, 07:58
  mydeviance <- function(x, ...) {
    dev <- deviance(x)
    if (!is.null(dev))
      dev
    else extractAIC(x, k = 0)[2L]
  }
  cut.string <- function(string) {
    if (length(string) > 1L)
      string[-1L] <- paste0("\n", string[-1L])
    string
  }
  re.arrange <- function(keep) {
    namr <- names(k1 <- keep[[1L]])
    namc <- names(keep)
    nc <- length(keep)
    nr <- length(k1)
    array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr,
                                                           namc))
  }
  step.results <- function(models, fit, object, usingCp = FALSE) {
    change <- sapply(models, "[[", "change")
    rd <- sapply(models, "[[", "deviance")
    if(is.matrix(rd)) rd <- rd[2,]
    dd <- c(NA, abs(diff(rd)))
    rdf <- sapply(models, "[[", "df.resid")
    ddf <- c(NA, diff(rdf))
    AIC <- sapply(models, "[[", "AIC")
    heading <- c("Stepwise Model Path \nAnalysis of Deviance Table",
                 "\nInitial Model:", deparse(formula(object)), "\nFinal Model:",
                 deparse(formula(fit)), "\n")
    aod <- data.frame(Step = I(change), Df = ddf, Deviance = dd,
                      `Resid. Df` = rdf, `Resid. Dev` = rd, AIC = AIC,
                      check.names = FALSE)
    if (usingCp) {
      cn <- colnames(aod)
      cn[cn == "AIC"] <- "Cp"
      colnames(aod) <- cn
    }
    attr(aod, "heading") <- heading
    fit$anova <- aod
    fit
  } ## end step.results
  ## ---
  Terms <- terms(object)
  object$call$formula <- object$formula <- Terms
  md <- missing(direction)
  direction <- match.arg(direction)
  backward <- direction == "both" | direction == "backward"
  forward <- direction == "both" | direction == "forward"
  ## scope
  lform <- formula(object)
  if (is.null(scope))  ## !! was missing(scope)
    scope <- if (expand) terms2order(object) else lform
  if (!is.list(scope)) {
    if (is.character(scope))
      scope <- as.formula(paste("~.+",paste(scope, collapse="+")))
    scope <- list(lower=NULL, upper=scope)
  }
  if (is.null(names(scope))&length(scope)==2)
    names(scope) <- c("lower","upper")
  fdrop <- if (length(fdrop <- scope$lower)) {
    attr(terms(update.formula(lform, fdrop)), "factors")
  } else numeric()
  if (length(fadd <- scope$upper)) {
    lform <- update.formula(lform, fadd)
    fadd <- attr(terms(lform), "factors")
  }
  ## data
  lcalldata <- object$call$data
  ldata <- object$allvars
  lvars <- all.vars(lform)
  if (any(lvars%nin%names(ldata)))
    ldata <- eval(object$call$data, envir=environment(formula(object)) )
  if (any(li <- lvars%nin%names(ldata)))
    stop("!step.regr! variable  ", paste(lvars[li], collapse=", "), "  not available")
  stepdata <- na.omit(ldata[,lvars])
  object$call$data <- stepdata
  if (length(object$funcall)) object$funcall$data <- stepdata
  models <- vector("list", steps)
  if (!is.null(keep))
    keep.list <- vector("list", steps)
  n <- nobs(object, use.fallback = TRUE)
  fit <- object
  if (inherits(fit, "survreg")) fit$residuals <- NULL
  # bAIC <- extractAIC(fit, scale, k = k, ...)
  bAIC <- extractAIC(fit, scale, k = k)
  
  edf <- bAIC[1L]
  bAIC <- bAIC[2L]
  if (is.na(bAIC))
    stop("AIC is not defined for this model, so `step` cannot proceed")
  nm <- 1
  if (trace) {
    cat("Start:  AIC=", format(round(bAIC, 2)), "\n",
        cut.string(deparse(formula(fit))), "\n\n", sep = "")
    utils::flush.console()
  }
  models[[nm]] <-
    list(deviance = mydeviance(fit), df.resid = n - edf, change = "",
         AIC = bAIC)
  if (!is.null(keep))
    keep.list[[nm]] <- keep(fit, bAIC)
  usingCp <- FALSE
  
  ## ------------------------
  lrgopt <- options(notices=FALSE)
  on.exit(options(lrgopt))
  while (steps > 0) {
    steps <- steps - 1
    AIC <- bAIC
    ffac <- attr(Terms, "factors")
    scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
    aod <- NULL
    change <- NULL
    ## backward
    if (backward && length(scope$drop)) {
      # aod <- drop1(fit, scope$drop, scale = scale, trace = trace,
                   # k = k, sorted=FALSE, ...)
      aod <- my.drop1.mlm(fit, scope$drop, scale = scale, trace = trace,
                   k = k, sorted=FALSE,test = "Pillai",total = FALSE)
      rn <- row.names(aod)
      # row.names(aod) <- c(rn[1L], paste("-", rn[-1L], sep = " "))
      row.names(aod) <- paste("-", rn, sep = " ")
      
      if (any(aod$dfnum == 0, na.rm = TRUE)) {
        zdf <- aod$dfnum == 0 & !is.na(aod$dfnum)
        change <- rev(rownames(aod)[zdf])[1L]
      }
    }
    ## forward
    if (is.null(change)) {
      if (forward && length(scope$add)) {
        # aodf <- add1(fit, scope$add, scale = scale, trace = trace,
        #              k = k, ...)
        aodf <- my.add1.mlm(fit, scope$add, scale = scale, trace = trace,
                     k = k)
        rn <- row.names(aodf)
        # row.names(aodf) <- c(rn[1L], paste("+", rn[-1L], sep = " "))
        row.names(aodf) <- paste("+", rn, sep = " ")
        aod <- if (is.null(aod)) aodf  else {
          names(aodf) <- names(aod)
          rbind(aod, aodf[-1, , drop = FALSE])
        }
      }
    }
    ## backward or forward
    attr(aod, "heading") <- NULL
    # nzdf <- if (!is.null(aod$Df))
    #   aod$Df != 0 | is.na(aod$Df)
    nzdf <- if (!is.null(aod$dfnum))
      aod$dfnum != 0 | is.na(aod$dfnum)
    aod <- aod[nzdf, ]
    if (is.null(aod) || ncol(aod) == 0)
      break
    # nc <- match(c("Cp", "AIC"), names(aod))
    nc <- match(c("Cp", "AIC","F.stat"), names(aod))
    nc <- nc[!is.na(nc)][1L]
    o <- order(aod[, nc])
    if (trace)
      print(aod[o, ])
    if (o[1L] == 1)  break
    change <- rownames(aod)[o[1L]]
    ## update
    usingCp <- match("Cp", names(aod), 0L) > 0L
    ##    if (is.null(change))  break  else {
    fit <- update(fit, paste("~ .", change), evaluate = FALSE)
    fit <- eval.parent(fit)
    if (inherits(fit, "survreg")) fit$residuals <- NULL
    ##-    nnew <- nobs(fit, use.fallback = TRUE)
    ##-     if (all(is.finite(c(n, nnew))) && nnew != n) {
    ##-       warning(":step.regr: number of rows in use has changed: \n  ",
    ##-               nnew," observations instead of ", n)
    ##-       n <- nnew
    ##-     }
    Terms <- terms(fit)
    # bAIC <- extractAIC(fit, scale, k = k, ...)
    bAIC <- extractAIC(fit, scale, k = k)
    
    edf <- bAIC[1L]
    bAIC <- bAIC[2L]
    ## output
    if (trace) {
      cat("\nStep:  AIC=", format(round(bAIC, 2)), "\n",
          cut.string(deparse(formula(fit))), "\n\n", sep = "")
      utils::flush.console()
    }
    ##        if (bAIC >= AIC + 1e-07)
    ##  else   break
    nm <- nm + 1
    models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n - edf,
                         change = change, AIC = bAIC)
    if (!is.null(keep))
      keep.list[[nm]] <- keep(fit, bAIC)
  }
  if (!is.null(keep))
    fit$keep <- re.arrange(keep.list[seq(nm)])
  lv <- all.vars(formula(fit))
  lnobs <- nrow(na.omit(ldata[,lvars]))
  lno <- NROW(na.omit(ldata[,lv]))
  if (lnobs<lno) {
    notice(" step.regr: Step worked only on ", lnobs,
           " observations, whereas result has ", lno,
           " -- due to missing values in dropped terms")
    fit <- update(fit, data=ldata)
  }
  fit$call$data <- lcalldata
  step.results(models = models[seq(nm)], fit, object, usingCp)
}

i.add1na <- function (object, scope)
{
  ##  determine rows with NA`s in model.frame for expanded model
  Terms <-
    terms(update.formula(object, paste("~.+",paste(scope, collapse = "+"))))
  fc <- object$call
  fc$formula <- Terms
  fob <- list(call = fc, terms = Terms)
  class(fob) <- oldClass(object)
  m <- model.frame(fob, xlev = object$xlevels)
  r <- cbind(object$resid) ## resid(object)
  if (nrow(r)!=nrow(m)) {
    notice(gettextf("add1! using the %d/%d rows from a combined fit",
                    nrow(m), nrow(r)), domain = NA)
    lna <- !row.names(r)%in%row.names(m)
  }
  else lan <- NULL
}
