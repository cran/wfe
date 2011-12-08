

pwfe <- function (formula, data, outcome, treat = "treat.name",
                  unit.index, time.index = NULL, method = "unit", within.unit = TRUE,
                  qoi = c("ate", "att"), estimator = NULL, C.it = NULL,
                  White.alpha = 0.05, hetero.se = TRUE, auto.se = TRUE) {


  pwfe.call <- match.call()
  tn.row <- nrow(data) # total number of rows in data


### Warnings
  ## Warning for missing unit & time index
  if (missing(unit.index))
    stop("'unit.index' or index for strata should be provided")

  if (is.null(time.index) & method == "time")
    stop("'time.index' should be provided")
  
  ## Warning for methods
  if(method=="time" && !is.null(estimator) && estimator == "fd")
    stop("First Difference is not compatible with 'time' method")
  
  if(is.null(time.index) && !is.null(estimator) && estimator  == "fd")
    stop("First Difference cannot calculate when 'time.index' is missing")
  
  
  ## Warning for C.it
  if (!is.null(C.it) && !is.numeric(C.it) && length(C.it)!= tn.row)
    stop("'C.it' must be a numeric vector with length equal to number of observations")
  
  
  
  ## cat("warnings done:\n")
  ## C.it
  ## Default for ATE
  if (is.null(C.it)){
    data$C.it <- as.integer(rep(1, nrow(data)))
  }

  ## White.alpha
  if (is.null(White.alpha)){
    White.alpha <- 0.05
  } else {
    White.alpha <- White.alpha
  }

  ## Warning for binary treatment
  ## treat should be 0,1 where 1 indicates treatment


### Unit and Time index
  ## Creating time index for strata fixed effect analysis

  numeric.u.index <- as.numeric(as.factor(data[,unit.index]))
  numeric.u.index[is.na(numeric.u.index)] <- 0
  ## handling missing unit index
  uniq.u <- unique(na.omit(numeric.u.index))
  uniq.u <- sort(uniq.u[!(uniq.u %in% 0)])
  J.u <- length(uniq.u)

  ## unit index
  data$u.index <- Index(numeric.u.index, uniq.u, J.u, tn.row)

  ## time index
  if (is.null(time.index)) {
    data$t.index <- GenTime(data$u.index, tn.row, length(uniq.u))
  } else {
    numeric.t.index <- as.numeric(as.factor(data[,time.index]))
    numeric.t.index[is.na(numeric.t.index)] <- 0
    ## handling missing time index
    uniq.t <- unique(na.omit(numeric.t.index))
    uniq.t <- sort(uniq.t[!(uniq.t %in% 0)])
    # needs to sort for unbalnced panel, See Index()
    J.t <- length(uniq.t)
    data$t.index <- Index(numeric.t.index, uniq.t, J.t, tn.row)
  }

  uniq.t <- unique(data$t.index)

  ## unique unit number and sorting data
  if (method == "unit"){
    unit.number <- length(uniq.u)
    data <- data[order(data$u.index, data$t.index),]   
  }
  if (method == "time"){
    unit.number <- length(uniq.t)
    data <- data[order(data$t.index, data$u.index),]   
  }
 
  ## treatment variable
  data$TR <- data[,treat]


  ## quantity of interest: qoi

  if (missing(qoi)){
    causal <- "ate"
    ate.n <- 1
  }

  ate.n <-  att.n <- 0
  if (qoi == "ate"){
    causal <- "ATE (Average Treatment Effect)"
    ate.n <- 1
  }
  if (qoi == "att"){
    causal <- "ATT (Average Treatment Effect for the Treated)"
    att.n <- 1
  }

  if (is.null(estimator)) {
    est <- "NULL"
  }
  if (!is.null(estimator) && estimator == "fd"){
    est <- "FD (First-Difference)"
  }
  if (!is.null(estimator) && estimator == "did"){
    est <- "DID (Difference-in-Differences)"
  } 


  ## change formula without intercept
  a <- unlist(strsplit(as.character(formula), "~"))
  formula.ni <- as.formula(paste(a[2], "~ -1 + ",  a[3]))

 
  data$index <- 1:nrow(data) # indexing for saving y.star


  ## transform Y_{it}
  if (method=="unit" && within.unit == TRUE){
    # regression for propensity score estimation is run within unit
    Data.ps <- c()
    for (j in 1:length(uniq.t)){
        #print(j)
        pool <- data[data$t.index == uniq.t[j], ] # subset data by time j
        n.pool <- nrow(pool) # gives total number of units for the year j
        n.treat <- sum(pool[,treat]) # number of treated units in a given year t
        n.cont <- sum(1-pool[,treat]) # number of controlled units in a given year t
        # estimate propensity score within each time t
        fit.ps <- bayesglm(formula, family=binomial(link="logit"),
                           data = pool)
        pool$fitted.ps <- fitted(fit.ps) # save the estimated propensity score
        pool.u.index <- unique(pool$u.index)
        pool$ystar <- Transform(pool[,outcome], pool[,treat], pool$fitted.ps)
        Data.ps <- rbind(Data.ps, pool)
    }
    data$y.star[Data.ps$index] <- Data.ps$ystar

  }

  if (method=="time" && within.unit == TRUE){
    # regression for propensity score estimation is run within unit
    Data.ps <- c()
    for (i in 1:length(uniq.u)){
        #print(j)
        pool <- data[data$u.index == uniq.u[i], ] # subset data by unit i
        n.pool <- nrow(pool) # gives total number of times for the unit i
        n.treat <- sum(pool[,treat]) # number of treated times in a given unit i
        n.cont <- sum(1-pool[,treat]) # number of controlled times in a given unit i
        # estimate propensity score within each time t
        fit.ps <- bayesglm(formula, family=binomial(link="logit"),
                           data = pool)
        pool$fitted.ps <- fitted(fit.ps) # save the estimated propensity score
        pool.t.index <- unique(pool$t.index)
        pool$ystar <- Transform(pool[,outcome], pool[,treat], pool$fitted.ps)
        Data.ps <- rbind(Data.ps, pool)
    }
    data$y.star[Data.ps$index] <- Data.ps$ystar

  }


  
  ## transform Y_{it}
  if (within.unit == FALSE){
    # regression for propensity score estimation is on entire data
    n.treat <- sum(data[,treat]) # number of treated units in a given year t
    n.cont <- sum(1-data[,treat]) # number of controlled units in a given year t
    # estimate propensity score within each time t
    fit.ps <- bayesglm(formula, family=binomial(link="logit"),
                       data = data)
    data$fitted.ps <- fitted(fit.ps) # save the estimated propensity score
    data$y.star <- Transform(data[,outcome], data[,treat], data$fitted.ps)
  }

  
  ## Unit fixed effects models
  if ( (method=="unit" & qoi=="ate" & is.null(estimator) ) | (method=="unit" & qoi=="att" & is.null(estimator)) ) {
    W <- GenWeights(data$u.index, data$t.index, data$TR, data$C.it, tn.row, length(uniq.u), length(uniq.t), ate.n, att.n, length(uniq.u)*length(uniq.t))
    W <- matrix(W, nrow=length(uniq.t), ncol=length(uniq.u), byrow=T)
    data$W.it <- Vectorize(as.matrix(W), data$t.index, data$u.index, tn.row)
  }

  ## Time fixed effects models
  if(( (method=="time") && (qoi == "ate") && (is.null(estimator)) ) | ( (method=="time") && (qoi == "att") && (is.null(estimator)) )) {  
    W <- GenWeights(data$t.index, data$u.index, data$TR, data$C.it, tn.row, length(uniq.t), length(uniq.u), ate.n, att.n, length(uniq.t)*length(uniq.u))
    W <- matrix(W, nrow=length(uniq.t), ncol=length(uniq.u), byrow=F)
    data$W.it <- Vectorize(as.matrix(W), data$t.index, data$u.index, tn.row)
  }


  ## Within Unit First Difference

  if(( (method=="unit") && (qoi == "ate") && (!is.null(estimator) && estimator == "fd")) | ((method == "unit") && (qoi =="att") && (!is.null(estimator) && estimator == "fd"))) {
    W <- GenWeightsFD(data$u.index, data$t.index, data$TR, data$C.it, tn.row, length(uniq.u), length(uniq.t), ate.n, att.n)
    W <- matrix(W, nrow=length(uniq.t), ncol=length(uniq.u), byrow=T)
    data$W.it <- Vectorize(as.matrix(W), data$t.index, data$u.index, tn.row)
  }

  Data.final <- as.data.frame(cbind(data$y.star, data$TR, data$W.it, data$u.index, data$t.index))
  colnames(Data.final) <- c("y.star", "TR", "W.it", "u.index", "t.index")

  
  ## Creating the matrix for final analysis
  Data.dm <- as.data.frame(matrix(NA, ncol = 2, nrow = nrow(Data.final)))
  Data.wdm <- as.data.frame(matrix(NA, ncol = 2, nrow = nrow(Data.final)))


  ## column names
  ## print(colnames(Data.final))
  colnames(Data.dm) <- colnames(Data.wdm) <- c("outcome", "treat")

  ## de-meaning for fast calculation
  for (k in 1:2){
    if (method == "unit"){
      w.wdemean <- WWDemean(Data.final[,k], Data.final$W.it, Data.final$u.index, unit.number, tn.row)
      demean <- Demean(Data.final[,k], Data.final$u.index, unit.number, tn.row)
    }
    if (method == "time"){
      w.wdemean <- WWDemean(Data.final[,k], Data.final$W.it, Data.final$t.index, unit.number, tn.row)
      demean <- Demean(Data.final[,k], Data.final$u.index, unit.number, tn.row)
    }
    Data.wdm[,k] <- as.vector(w.wdemean)
    Data.dm[,k] <- as.vector(demean)
  }
  
  Data.dm$unit <- Data.final$u.index
  Data.dm$time <- Data.final$t.index
  Data.wdm$unit <- Data.final$u.index
  Data.wdm$time <- Data.final$t.index
  
  ## Print de(weighted)-meaned data

  ## cat("#####", "\n","De-meaned Data:","\n")
  ## print(Data.dm)
  ## cat("#####", "\n","sqrt(W) x (weighted demeaned):","\n")
  ## print(Data.wdm)

  ## final regression on weighted demeaned data
  fit.final <- lm(outcome ~ -1 + treat, data = Data.wdm)
  fit.ols <- lm(outcome ~ -1 + treat, data = Data.dm)

  ## print(summary(fit.final))           

  ## set up data frame, with support for standard and modified responses
  
  mf.final <- model.frame(formula="outcome ~ -1 + treat", data=Data.wdm)
  X <- model.matrix(attr(mf.final, "terms"), data=mf.final)
  y <- model.response(mf.final, "numeric")
  p <- ncol(X)
  
  coef.wls <- fit.final$coef
  coef.ols <- fit.ols$coef
  d.f <- fit.final$df - unit.number

  sigma2 <- sum(resid(fit.final)^2)/d.f
  vcov.wls <- vcov(fit.final)*((fit.final$df)/(fit.final$df - unit.number))
  var.cov <- vcov.wls
  vcov.ols <- vcov(fit.ols)

  residual <-  resid(fit.final)*1/sqrt(Data.final$W.it)
  resid.ols <- resid(fit.ols)

  ## save residuals
  Data.wdm$u.tilde <- sqrt(Data.final$W.it)*resid(fit.final)
  Data.dm$u.hat <- resid(fit.ols)


  
### Robust Standard Errors
  
  ## contructing de(weighted)-meaned X matrix
  X.tilde <- as.matrix(Data.wdm[,2]) # NT x 1 (where p is number of covariates)
  X.hat <- as.matrix(Data.dm[,2]) # NT x 1 (where p is number of covariates)

  ## residuals
  u.tilde <- as.matrix(Data.wdm$u.tilde) # NT x 1
  u.hat <- as.matrix(Data.dm$u.hat)  # NT x 1
  
  ## constructing Omega.hat for WFE
  
  ## 1. independence across observations but heteroskedasticity (Eq 12)
  Omega.hat.heauto <- (1/unit.number)*( t(X.tilde) %*% diag(diag(u.tilde %*% t(u.tilde)), nrow = tn.row) %*% X.tilde )
  
  ## 2. arbitrary autocorrelation as well as heteroskedasticity (Eq 13)
  Omega.hat.hoauto <- (1/unit.number)*( t(X.tilde) %*% u.tilde %*% t(u.tilde) %*% X.tilde )


  ## by the same token, Omega.hat for FE
  ## 1. independence across observations but heteroskedasticity (Eq 12)
  Omega.hat.fe.heauto <- (1/unit.number)*( t(X.hat) %*% diag(diag(u.hat %*% t(u.hat)), nrow = tn.row) %*% X.hat)
  
  ## 2. arbitrary autocorrelation as well as heteroskedasticity (Eq 13)
  Omega.hat.fe.hoauto <- (1/unit.number)*( t(X.hat) %*% u.hat %*% t(u.hat) %*% X.hat )
  
  ## (Robust) standard errors

  if ((hetero.se == TRUE) & (auto.se == TRUE)){# Default is Arellano
    var.cov <- solve((t(X.tilde)%*% X.tilde)/unit.number) %*% Omega.hat.heauto %*% solve((t(X.tilde)%*% X.tilde)/unit.number)
    std.error <- "Heteroscedastic / Autocorrelation Robust Standard Error"
    Psi.hat.fe <- solve((t(X.hat)%*% X.hat)/unit.number) %*% Omega.hat.fe.heauto %*% solve((t(X.hat)%*% X.hat)/unit.number)
  }
  if ((hetero.se == TRUE) & (auto.se == FALSE)){# independence across observations but heteroskedasticity
    var.cov <- solve((t(X.tilde)%*% X.tilde)/unit.number) %*% Omega.hat.hoauto %*% solve((t(X.tilde)%*% X.tilde)/unit.number)
    std.error <- "Heteroscedastic Robust Standard Error"
    Psi.hat.fe <- solve((t(X.hat)%*% X.hat)/unit.number) %*% Omega.hat.fe.hoauto %*% solve((t(X.hat)%*% X.hat)/unit.number)
  }
  if ((hetero.se == FALSE) & (auto.se == FALSE)){# indepdence and homoskedasticity
    var.cov <- vcov.wls
    std.error <- "Homoskedastic Standard Error"
    Psi.hat.fe <- vcov.ols
  }
  if ((hetero.se == FALSE) & (auto.se == TRUE)){# Kiefer
    stop ("robust standard errors with autocorrelation and homoskedasiticy is not supported")
  }



### White (1980) Test: Theorem 4

  Gamma.hat <- 1/unit.number *( t(X.tilde) %*% u.tilde %*% t(u.hat) %*% X.hat )

  Phi.hat <- var.cov + Psi.hat.fe - solve((1/unit.number)*(t(X.hat)%*%
      X.hat)) %*% Gamma.hat %*% solve((1/unit.number)*(t(X.tilde) %*%
      X.tilde)) - solve((1/unit.number)*(t(X.tilde)%*%X.tilde)) %*%
      Gamma.hat %*% solve((1/unit.number)*(t(X.hat)%*%X.hat))
 
  ## White test: null hypothesis is ``no misspecification''

  white.stat <- unit.number * t(coef.ols - coef.wls) %*% solve(Phi.hat) %*% (coef.ols - coef.wls)
  test.null <- pchisq(as.numeric(white.stat), df=p, lower.tail=F) < White.alpha
  white.p <- pchisq(as.numeric(white.stat), df=p, lower.tail=F)


### Saving results
  z <- list(coefficients = coef.wls,
            x = X,                      
            y = y,
            call = pwfe.call,
            vcov = var.cov,
            se = sqrt(diag(var.cov)),
            sigma = sqrt(sigma2),
            df = d.f,
            residuals = residual,
            weights = W,
            uniq.n.units = unit.number,
            method = method,
            causal = causal,
            est = est,
            std.error = std.error,     
            White.pvalue = white.p,
            White.alpha = White.alpha,
            White.stat = white.stat,
            White.test = test.null)
  class(z) <- "pwfe"
  z



}


print.pwfe <- function(x,...){
    cat("Call:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print(x$coefficients)
    cat("\nStd.Err:\n")
    print(x$se)
}



summary.pwfe <- function(object, signif.stars = getOption("show.signif.stars"),...){
    se <- object$se
    sigma <- object$sigma
    df <- object$df
    tval <- coef(object) / se
    TAB <- cbind(Estimate = coef(object),
                 Std.Err = se,
                 t.value = tval,
                 p.value = 2*pt(-abs(tval), df = object$df))
    res <- list(call = object$call,
                coefficients = TAB,
                sigma = object$sigma,
                df = object$df,
                Weights = object$weights,
                residuals = object$residuals,
                method = object$method,
                causal = object$causal,
                estimator = object$est,
                std.error = object$std.error,
                White.pvalue = object$White.pvalue,
                White.alpha = object$White.alpha,
                White.stat = object$White.stat,
                White.test = object$White.test)
    class(res) <- "summary.pwfe"
    res
}


print.summary.pwfe <- function(x, ...){
    cat("\nMethod:", x$method, "Fixed Effects\n")
    cat("\nQuantity of Interest:", x$causal)
    cat("\nEstimator:", x$estimator)
    cat("\nStandard Error:", x$std.error)
    cat("\n")
    cat("\n")
    cat("Call:\n")
    print(x$call)
    cat("\n")
    printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
    cat("\nResidual standard error:", format(signif(x$sigma,
        4)), "on", x$df, "degrees of freedom")
    cat("\nWhite statistics for functional misspecification:", x$White.stat, "with Pvalue=", x$White.pvalue)
    cat("\nReject the null of NO misspecification:", x$White.test)
    cat("\n")
}
