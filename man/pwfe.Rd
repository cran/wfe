\name{pwfe}
\alias{pwfe}
\title{Fitting \code{pwfe} for Causal Inference with Propensity Score Weighting}
\description{
  \code{pwfe} is used to fit weighted fixed effects model for causal
  inference after transforming outcome variable based on estimated
  propensity score. \code{pwfe} also derives the regression weights for
  different causal quantity of interest.
}
\usage{
pwfe(formula, data, outcome, treat = "treat.name",
     unit.index, time.index = NULL, method = "unit", within.unit = TRUE,
     qoi = c("ate", "att"), estimator = NULL, 
     C.it = NULL, White.alpha = 0.05,
     hetero.se = TRUE, auto.se = TRUE)
}
\arguments{

  \item{formula}{a symbolic description of the model for estimating
	propensity score. The formula should not include dummmies for
	fixed effects. The details of model specifications are given
	under `Details'.}
      
  \item{data}{data frame containing the variables in the model.
  }

  \item{outcome}{a character string indicating the name of outcome
  variable.
  }

  \item{treat}{a character string indicating the name of treatment
    variable used in the models. The treatment should be binary
    indicator (integer with 0 for the control group and 1 for the
    treatment group).
  }
  \item{unit.index}{a character string indicating the name of unit
    variable used in the models. The index of unit should be factor.
  }
  \item{time.index}{a character string indicating the name of time
    variable used in the models. The index of time should be factor.
  }
  \item{method}{method for weighted fixed effects regression, either
    \code{unit} for unit fixed effects; \code{time} for time fixed
    effects. The default is \code{unit}.
  }

  \item{within.unit}{a logical value indicating whether propensity score
  is estimated within unit. The default is \code{TRUE}. 
  }

  \item{qoi}{one of \code{"ate"} or \code{"att"}. The default is
    \code{"ate"}.
  }

  \item{estimator}{an optional character string indicating the
  estimating method. One of \code{"fd"} or \code{"did"}. The default is
  \code{NULL}.
  }

  \item{C.it}{an optional non-negative numeric vector specifying
    relative weights for each unit of analysis.
  }
  \item{White.alpha}{level of functional specification test. See White
    (1980) and Imai . The default is \code{0.05}.
  }
  \item{hetero.se}{a logical value indicating whether heteroskedasticity
    across units is allowed in calculating standard errors. 
  }
  \item{auto.se}{a logical value indicating whether arbitrary
    autocorrelation is allowed in calculating standard errors.
  }

}
\value{
 \code{pwfe} returns an object of class "pwfe", a list that contains the
 components listed below.


 The function \code{summary} (i.e., \code{summary.pwfe}) can be used to
 obtain a table of the results.

 \item{coefficients}{a named vector of coefficients}
 \item{residuals}{the residuals, that is respons minus fitted values}
 \item{df}{the degree of freedom}
 \item{weights}{weights calculated from the model}
 \item{call}{the matched call}
 \item{causal}{causal quantity of interest}
 \item{estimator}{the estimating method}
 \item{method}{call of the method used}
 \item{vcov}{the variance covariance matrix}
 \item{x}{the design matrix}
 \item{y}{the response vector}
  
}
\details{
  
  To fit the weighted unit (time) fixed effects model with propensity
  score weighting, use the syntax for the formula, \code{treat ~ x1 +
  x2}, where \code{treat} is the name of treatment variable and
  \code{x1} and \code{x2} are unit (time) varying
  covariates. \code{bayesglm} will be used to estimate propensity
  scores. If \code{within.unit = TRUE}, propensity score is estimated
  within time (unit) separately if \code{method} is \code{unit}
  (\code{time}). Otherwise, propensity score is estimated on entire data
  at once.

  The estimated propensity scores will be used to transform the
  \code{outcome} variable as described in Imai and Kim (2011).

  \code{pwfe} calculates weights based on different underlying causal
  quantity of interest: Average Treatment Effect (\code{qoi = "ate"}) or
  Average Treatment Effect for the Treated (\code{qoi = "att"}).

  One can further set estimating methods: First-Difference
  (\code{estimator ="fd"}) or Difference-in-differences (\code{estimator
  = "did"}).

  To specify different ex-ante weights for each unit of analysis, use
  non-negative weights \code{C.it}. For instance, using the survey
  weights for \code{C.it} enables the estimation fo the average
  treatement effect for the target population.
  

}
\references{
  Imai, Kosuke and In Song Kim. (2011) ``Understanding and Improving
  Linear Fixed Effects Regression Models for Causal Inference.''
  Technical Report, Department of Politics, Princeton University.
  available at \url{http://imai.princeton.edu/research/FEmatch.html}
    
  White, Halbert. (1980) `Using Least Squares to Approximate Unknown
  Regression Functions.''  International Economic Review, 21, 1,
  149--170.
}
\seealso{
  \code{wfe} for fitting weighted fixed effect models.
}
\examples{
### NOTE: this example illustrates the use of wfe function with randomly
### generated panel data with arbitrary number of units and time.

## generate panel data with number of units = N, number of time = Time
N <- 10 # number of distinct units
Time <- 15 # number of distinct time

## generate treatment variable
treat <- matrix(rbinom(N*Time, size = 1, 0.25), ncol = N)
## make sure at least one observation is treated for each unit
while ((sum(apply(treat, 2, mean) == 0) > 0) | (sum(apply(treat, 2, mean) == 1) > 0) |
       (sum(apply(treat, 1, mean) == 0) > 0) | (sum(apply(treat, 1, mean) == 1) > 0)) {
  treat <- matrix(rbinom(N*Time, size = 1, 0.25), ncol = N)
}
treat.vec <- c(treat)

## unit fixed effects
alphai <- rnorm(N, mean = apply(treat, 2, mean))

## geneate two random covariates
x1 <- matrix(rnorm(N*Time, 0.5,1), ncol=N)
x2 <- matrix(rbeta(N*Time, 5,1), ncol=N)
x1.vec <- c(x1)
x2.vec <- c(x2)
## generate outcome variable
y <- matrix(NA, ncol = N, nrow = Time)
for (i in 1:N) {
    y[, i] <- alphai[i] + treat[, i] + x1[,i] + x2[,i] + rnorm(Time)
}
y.vec <- c(y)

## generate unit and time index
unit.index <- rep(1:N, each = Time)
time.index <- rep(1:Time, N)

Data.str <- as.data.frame(cbind(y.vec, treat.vec, unit.index, x1.vec, x2.vec))
colnames(Data.str) <- c("y", "tr", "strata.id", "x1", "x2")

Data.obs <- as.data.frame(cbind(y.vec, treat.vec, unit.index, time.index, x1.vec, x2.vec))
colnames(Data.obs) <- c("y", "tr", "unit", "time", "x1", "x2")


############################################################
# Example 1: Stratified Randomized Experiments
############################################################

## run the weighted fixed effect regression with strata fixed effect.
## Note: the quantity of interest is Average Treatment Effect ("ate")
## and the standard errors allow heteroskedasticity and arbitrary
## autocorrelation.


### Average Treatment Effect
ps.ate <- pwfe(tr~ x1+x2, data = Data.str, outcome = "y", treat = "tr",
               unit.index = "strata.id", method = "unit", within.unit = TRUE,
               qoi = "ate", hetero.se=TRUE, auto.se=TRUE)
## summarize the results
summary(ps.ate)

### Average Treatment Effect for the Treated
ps.att <- pwfe(tr~ x1+x2, data = Data.str, outcome = "y", treat = "tr",
               unit.index = "strata.id", method = "unit", within.unit = TRUE,
               qoi = "att", hetero.se=TRUE, auto.se=TRUE)
## summarize the results
summary(ps.att)


############################################################
# Example 2: Observational Studies with Unit Fixed-effects
############################################################

## run the weighted fixed effect regression with unit fixed effect.
## Note: the quantity of interest is Average Treatment Effect ("ate")
## and the standard errors allow heteroskedasticity and arbitrary
## autocorrelation.

### Average Treatment Effect
ps.obs <- pwfe(tr ~ x1+x2, data = Data.obs, outcome = "y", treat = "tr",
               unit.index = "unit", time.index = "time",
               method = "unit", within.unit = TRUE,
               qoi = "ate", hetero.se=TRUE, auto.se=TRUE)

## summarize the results
summary(ps.obs)

## extracting weigths
summary(ps.obs)$Weights

### Average Treatment Effect with First-difference

ps.fd <- pwfe(tr ~ x1+x2, data = Data.obs, outcome = "y", treat = "tr",
              unit.index = "unit", time.index = "time",
              method = "unit", within.unit = TRUE,
              qoi = "ate", estimator = "fd", hetero.se=TRUE, auto.se=TRUE)

## summarize the results
summary(ps.fd)
}

\author{Kosuke Imai, Princeton University, \email{kimai@princeton.edu}
  and In Song Kim, Princeton University, \email{insong@princeton.edu}}

\keyword{regression}



