# BINF 531 - Statistical Bioinformatics
# McGill University
# Final Project - Helper Functions
# Michael Shamash, Wen Da Lu, Garrie Peng

# Custom function for Hutcheson t-test, derived from ecolTest package
Hutcheson_t_test <- function(x, y, shannon.base = exp(1),
                             alternative = "two.sided",
                             difference = 0) {
  dname <-paste ((deparse(substitute(x))), ", ", (deparse(substitute(y))))
  x <- drop(as.matrix(x))
  y <- drop(as.matrix(y))
  if (!is.numeric(x) | !is.numeric(y)) {
    stop("x and y must be numeric")
  }
  if (any(c(x,y) < 0, na.rm = TRUE)) {
    stop("x and y must be non-negative")
  }
  if (any(c(length(x) < 2, length(y) < 2))) {
    stop("x and y must contain at least two elements")
  }
  if (any(c(sum(x, na.rm = TRUE) < 3, sum(y, na.rm = TRUE) < 3))) {
    stop("x and y abundance must be at least two")
  }
  if (!requireNamespace("stats", quietly = TRUE)) {
    stop('Package "stats" is needed')
  }
  
  if (any(is.na(c(x, y)))) {
    x[is.na(x)] <- 0
    y[is.na(y)] <- 0
    warning("missing values in x and y replaced with zeroes")
  }
  
  alternative <- char.expand(alternative, c("two.sided",
                                            "less", "greater", "auto"))
  if (length(alternative) > 1L || is.na(alternative)) {
    stop("alternative must be \"two.sided\", \"less\" or \"greater\"")
  }
  length_diff <- length(x)-length(y)
  if (length_diff > 0){
    y <- c(y, rep(0, length_diff))
  }
  else if(length_diff < 0){
    x <- c(x, rep(0, abs(length_diff)))
  }
  xy <- matrix(c(x, y), ncol=2)
  N <- colSums(xy)
  H <- (N*log(N, shannon.base)-colSums(xy*log(xy, shannon.base),
                                       na.rm = TRUE))/N
  S<-(colSums(xy*log(xy, shannon.base)**2, na.rm = TRUE) -
        ((colSums(xy*log(xy, shannon.base), na.rm = TRUE)**2)/N))/(N**2)
  Hutchesontstat <- (diff(H[c(2, 1)])-difference)/sqrt(sum(S))
  df <- (sum(S)**2)/(sum(S**2/N))
  estimate_dif <- diff(H[c(2,1)])
  if (alternative == "auto") {
    alternative <- if (estimate_dif < 0) {
      "less"
    }else{
      "greater"
    }
  }
  
  if (alternative == "less") {
    pval <- pt(Hutchesontstat, df)
  }
  else if (alternative == "greater") {
    pval <- pt(Hutchesontstat, df, lower.tail = FALSE)
  }
  else {
    pval <- 2 * pt(-abs(Hutchesontstat), df)
  }
  names(Hutchesontstat) <- "Hutcheson t-statistic"
  names(df) <- "df"
  names(H) <- c("x", "y")
  mu <- difference
  names(mu) <- "difference in H'"
  rval <- list(statistic = Hutchesontstat, parameter = df, p.value = pval,
               estimate = H, null.value = mu,
               method="Hutcheson t-test for two communities",
               alternative = alternative, data.name=dname)
  class(rval) <- "htest"
  return(rval)
}