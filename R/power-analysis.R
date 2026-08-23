#' Unified effect-size / power / sample-size calculator
#'
#' Dispatches to the \pkg{pwr}, \pkg{lsr}, \pkg{WebPower}, and \pkg{pwrss}
#' packages depending on \code{test}, covering the t-test/ANOVA/correlation/
#' logistic/published-summary-statistics power-analysis patterns commonly
#' used to justify a sample size. Only pass the arguments relevant to the
#' chosen \code{test}; leave the rest at their defaults.
#'
#' @param test One of \code{"t"}, \code{"anova"}, \code{"correlation"},
#'   \code{"logistic"}, or \code{"summary_means"}.
#' @param x,y Raw samples for \code{test \%in\% c("t", "anova")} when
#'   \code{d} is not supplied directly; \code{d} is then computed via
#'   \code{lsr::cohensD(x, y)}.
#' @param d Cohen's d for \code{test = "t"} or Cohen's f for
#'   \code{test = "anova"}. With two raw groups, d is converted to f as
#'   \code{abs(d) / 2} for an equal-sized two-group ANOVA.
#' @param k Number of groups (for \code{test = "anova"}).
#' @param r Correlation coefficient (for \code{test = "correlation"}).
#' @param mu1,mu2,sd1,sd2 Means and standard deviations of the two groups
#'   (for \code{test = "summary_means"}), as reported in a published
#'   summary-statistics table.
#' @param q1_1,q3_1,q1_2,q3_2 Optional first/third quartiles for group 1/2,
#'   used to approximate \code{sd1}/\code{sd2} via \code{(Q3-Q1)/1.35} when
#'   the SD itself is not reported (\code{test = "summary_means"}).
#' @param paired Logical; whether the design is paired (\code{test \%in\%
#'   c("t", "summary_means")}).
#' @param paired_r Assumed within-pair correlation (for
#'   \code{test = "summary_means"} with \code{paired = TRUE}).
#' @param p0,p1 Baseline and alternative probabilities (for
#'   \code{test = "logistic"}).
#' @param n Sample size (leave \code{NULL} to solve for it; supply to solve
#'   for power instead).
#' @param power Target power (leave \code{NULL} to solve for it).
#' @param sig.level Significance level (type I error rate).
#' @return The underlying power-analysis result object returned by the
#'   dispatched package function (e.g. an object of class \code{"power.htest"}
#'   or \code{"pwrss"}), printable for a formatted summary.
#' @examples
#' \dontrun{
#' data(dfMigrant)
#' PowerAnalysisR(test = "t", x = dfMigrant$Age[dfMigrant$Migrant == 1],
#'                y = dfMigrant$Age[dfMigrant$Migrant == 0], power = 0.80)
#' PowerAnalysisR(test = "summary_means", mu1 = 3.4, mu2 = 2.8,
#'                sd1 = 3.5, sd2 = 3.1, power = 0.80)
#' }
#' @export
PowerAnalysisR <- function(test = c("t", "anova", "correlation", "logistic", "summary_means"),
                           x = NULL, y = NULL, d = NULL, k = NULL, r = NULL,
                           mu1 = NULL, mu2 = NULL, sd1 = NULL, sd2 = NULL,
                           q1_1 = NULL, q3_1 = NULL, q1_2 = NULL, q3_2 = NULL,
                           paired = FALSE, paired_r = 0.5, p0 = NULL, p1 = NULL,
                           n = NULL, power = 0.80, sig.level = 0.05){
  test <- match.arg(test)
  if (!is.numeric(sig.level) || length(sig.level) != 1L || sig.level <= 0 || sig.level >= 1){
    stop("sig.level must be between 0 and 1.") }
  if (!is.null(power) && (!is.numeric(power) || length(power) != 1L || power <= 0 || power >= 1)){
    stop("power must be NULL or a value between 0 and 1.")}
  require_package <- function(package) {
    if(!requireNamespace(package, quietly = TRUE)){
      stop(sprintf("Package '%s' is required for test = '%s'.", package, test), call. = FALSE)}
  }
  #### Derive Cohen's d from raw samples if not supplied directly ####
  if(test %in% c("t", "anova") && is.null(d) && !is.null(x) && !is.null(y)){
    require_package("lsr")
    d <- lsr::cohensD(x = x, y = y, method = if(paired) "paired" else "unequal")
    if(test == "anova"){ d <- abs(d) / 2 }
  }
  if(test == "t" && is.null(d)){ stop("Supply d or both x and y for test = 't'.") }
  if(test == "anova" && (is.null(d) || is.null(k))){ stop("Supply Cohen's f (d) and k for test = 'anova'.") }
  if(test == "correlation" && is.null(r)){ stop("Supply r for test = 'correlation'.") }
  if(test == "logistic" && (is.null(p0) || is.null(p1))){ stop("Supply p0 and p1 for test = 'logistic'.") }
  if(test == "summary_means" && (is.null(mu1) || is.null(mu2))){ stop("Supply mu1 and mu2 for test = 'summary_means'.") }
  result <- switch(test,
    t = { require_package("pwr"); pwr::pwr.t.test(n = n, d = d, sig.level = sig.level, power = power,
                                                   type = if(paired) "paired" else "two.sample") },
    anova = { require_package("pwr"); pwr::pwr.anova.test(k = k, n = n, f = d, sig.level = sig.level, power = power) },
    correlation = { require_package("pwr"); pwr::pwr.r.test(n = n, r = r, sig.level = sig.level, power = power) },
    logistic = {
      require_package("WebPower")
      WebPower::wp.logistic(n = n, p0 = p0, p1 = p1, alpha = sig.level, power = power, family = "normal")
    },
    summary_means = {
      require_package("pwrss")
      if(is.null(sd1) && !is.null(q1_1) && !is.null(q3_1)){ sd1 <- (q3_1 - q1_1) / 1.35 }
      if(is.null(sd2) && !is.null(q1_2) && !is.null(q3_2)){ sd2 <- (q3_2 - q1_2) / 1.35 }
      pwrss::pwrss.t.2means(mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2, paired = paired,
                            paired.r = paired_r, power = power, alpha = sig.level, n = n)
    })
  return(result)
}
