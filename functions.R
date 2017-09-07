
fitpowbend2<-function (x, trunc, start.value, ...) {
	dots <- list(...)
	if (any(x <= 0) | any(!sads:::is.wholenumber(x))) stop("All x must be positive integers")
	if (!missing(trunc)) {
		if (min(x) <= trunc) stop("truncation point should be lower than the lowest data value")
	}
	if (missing(start.value)) {
		gamahat <- function(ga, xvec) eq <- -sum(log(xvec)) - VGAM::zeta(ga, deriv = 1) * length(xvec)/VGAM::zeta(ga)
		shat <- uniroot(gamahat, interval = c(1.01, 20), xvec = x)$root
		oMhat <- 3
	}
	else {
		shat <- start.value[1]
		oMhat <- start.value[2]
	}
	if (!"method" %in% names(dots)) {
		dots$method <- "L-BFGS-B"
		if (!"lower" %in% names(dots)) 
			dots$lower = c(s = 0.1, oM = 1)
		if (!"upper" %in% names(dots)) 
			dots$upper = c(s = 2.999, oM = 16)
	}
	if (missing(trunc)) {
		LL <- function(s, oM) -sum(dpowbend(x, s = s, oM = oM, log = TRUE))
	}
	else {
		LL <- function(s, oM) -sum(dtrunc2("powbend", x = x, coef = list(s = s, oM = oM), trunc = trunc, log = TRUE))
	}
	result <- lapply(c(.01,.1,1,10,100),function(mult)try(do.call("mle2", c(list(LL, start = list(s = shat, oM = oMhat*mult), data = list(x = x)), dots)),silent=TRUE))
  result<-result[[which.min(sapply(result,function(xx){if(inherits(xx,'try-error'))return(Inf);return(xx@min)}))]]
  if(inherits(result,'try-error'))stop(result)
	new("fitsad", result, sad = "powbend", trunc = ifelse(missing(trunc), NaN, trunc))
}




fitnbinom2<- function (x, trunc = 0, start.value, ...) {
  dots <- list(...)
  if (any(x <= 0) | any(!sads:::is.wholenumber(x))) 
    stop("All x must be positive integers")
  if (!is.null(trunc)) {
    if (min(x) <= trunc) stop("truncation point should be lower than the lowest data value")
  }
  if (missing(start.value)) {
    muhat <- length(x)/(length(x) + mean(x))
    sizehat <- muhat * mean(x)
  }
  else {
    sizehat <- start.value[[1]]
    muhat <- start.value[[2]]
  }
  if (is.null(trunc)) {
    LL <- function(size, mu) -sum(dnbinom(x, size = size, mu = mu, log = TRUE))
  }
  else {
    LL <- function(size, mu) -sum(dtrunc2("nbinom", x = x, coef = list(size = size, mu = mu), trunc = trunc, log = TRUE))
  }
  result <- lapply(c(.01,.1,1,10,100),function(mult)try(suppressWarnings(do.call("mle2", c(list(LL, start = list(size = mult, mu = muhat), data = list(x = x)),method='Nelder-Mead',skip.hessian=TRUE, dots))),silent=TRUE))
  result<-result[[which.min(sapply(result,function(xx){if(inherits(xx,'try-error'))return(Inf);return(xx@min)}))]]
  if(inherits(result,'try-error'))stop(result)
  new("fitsad", result, sad = "nbinom", trunc = ifelse(is.null(trunc), NaN, trunc))
}

fitls2<-function (x, trunc, start.value, upper = length(x), ...) {
    dots <- list(...)
    if (any(x <= 0) | any(!sads:::is.wholenumber(x))) stop("All x must be positive integers")
    S <- length(x)
    N <- sum(x)
    if (missing(start.value)) {
        f1 <- function(a) S + a * log((a/(a + N)))
        sol <- uniroot(f1, interval = c(1/N, N))
        alfa <- sol$root
        X <- N/(N + alfa)
    } else {
        alfa <- start.value
    }
    if (!missing(trunc)) {
        if (min(x) <= trunc) stop("truncation point should be lower than the lowest data value")
        else LL <- function(N, alpha) -sum(dtrunc2("ls", x = x, coef = list(N = N, alpha = alpha), trunc = trunc, log = TRUE))
    }
    if (missing(trunc)) LL <- function(N, alpha) -sum(dls2(x, N, alpha, log = TRUE))
    result <- do.call("mle2", c(list(LL, start = list(alpha = alfa), 
        data = list(x = x), fixed = list(N = N), method = "Brent", 
        lower = 0, upper = upper), dots))
    if (abs(as.numeric(result@coef) - upper) < 1e-07) warning("mle equal to upper bound provided. \\n Try new value for the 'upper' argument")
    new("fitsad", result, sad = "ls", trunc = ifelse(missing(trunc), NaN, trunc))
}



dls2<-function (x, N, alpha, log = FALSE) {
    N[!is.finite(N) | N <= 0] <- NaN
    alpha[!is.finite(alpha) | alpha <= 0] <- NaN
    X <- N/(N + alpha)
    gama <- function(y) 1/log(1/(1 - y))
    y <- log(gama(X)) + x*log(X)-log(x)
    if (any(is.nan(y))) warning("NaNs produced")
    if (any(!sads:::is.wholenumber(x))) warning("non integer values in x")
    y[!sads:::is.wholenumber(x) | x < 1] <- 0
    if (log) return(y)
    else return(exp(y))
}


fitbs2<-function (x, trunc, ...) {
    dots <- list(...)
    if (any(x <= 0)) stop("All x must be positive")
    s <- length(x)
    n <- sum(x)
    if (!missing(x)) {
        if (!missing(trunc)) {
            if (min(x) <= trunc) stop("truncation point should be lower than the lowest data value")
            else LL <- function(N, S) -sum(dtrunc2("bs", x = x, coef = list(N = N, S = S), trunc = trunc, log = TRUE))
        }
        if (missing(trunc)) LL <- function(N, S) -sum(dbs2(x = x, N = N, S = S, log = TRUE))
        result <- do.call("mle2", c(list(minuslogl = LL, data = list(x = x), fixed = list(N = n, S = s), eval.only = TRUE), dots))
        result@details$convergence = 0
        new("fitsad", result, sad = "bs", trunc = ifelse(missing(trunc), NaN, trunc))
    }
}


dbs2<-function(x,N,S,log=FALSE){
    N[!is.finite(N) | N <= 0] <- NaN
    S[!is.finite(S) | S <= 0] <- NaN
    y <- log(S - 1) + (S - 2)*log(1 - x/N)-log(N)
    if (any(is.nan(y))) warning("NaNs produced")
    y[x < 0 | x > N] <- 0
    if (log) return(y)
    else return(exp(y))
}

fitgamma2<- function (x, trunc, start.value, ...) {
  dots <- list(...)
  if (any(x <= 0)) stop("All x must be positive")
  if (!missing(trunc)) {
    if (min(x) <= trunc) stop("truncation point should, be lower than the lowest data value")
  }
  if (missing(start.value)) {
    if (missing(trunc)) {
      ka <- (mean(x)/sd(x))^2
      theta <- var(x)/mean(x)
      kahat <- function(k, dados) {
        eq <- length(dados) * (log(k) - log(mean(dados)) - digamma(k)) + sum(log(dados))
      }
      ka <- uniroot(kahat, interval = c(min(theta, ka), max(theta, ka)), dados = x)$root
      theta <- mean(x)/ka
    }
    else {
      xh <- hist(x, plot = FALSE)
      xbr <- xh$breaks
      Eh <- matrix(ncol = 2, nrow = length(xbr) - 1)
      for (i in 1:(length(xbr) - 1)) {
        m1 <- matrix(c(1, 1, -1, 1), ncol = 2)
        Eh[i, ] <- solve(m1, xbr[i:(i + 1)])
      }
      P <- xh$counts/sum(xh$counts)
      P[P == 0] <- min(P[P > 0])
      Y <- log(P[-length(P)]) - log(P[-1]) - (log(Eh[-nrow(Eh), 2]) - log(Eh[-1, 2]))
      X1 <- Eh[-1, 1] - Eh[-nrow(Eh), 1]
      X2 <- log(Eh[-nrow(Eh), 1]) - log(Eh[-1, 1])
      st1 <- unname(coef(lm(Y ~ X1 + X2 - 1)))
      ka <- st1[2] + 1
      theta <- 1/st1[1]
    }
  }
  else {
    ka <- start.value[1]
    theta <- start.value[2]
  }
  if (missing(trunc)) {
    LL <- function(shape, rate) -sum(dgamma(x, shape, rate, log = TRUE))
  }
  else {
    LL <- function(shape, rate) -sum(dtrunc2("gamma", x = x, coef = list(shape = shape, rate = rate), trunc = trunc, log = TRUE))
  }
  result <- lapply(c(.01,.1,1,10,100),function(theta)try(suppressWarnings(do.call("mle2", c(list(LL, start = list(shape = ka, rate = 1/theta), data = list(x = x)), dots))),silent=TRUE))
  result<-result[[which.min(sapply(result,function(xx){if(inherits(xx,'try-error'))return(Inf);return(xx@min)}))]]
  new("fitsad", result, sad = "gamma", trunc = ifelse(missing(trunc), NaN, trunc))
}

fitmzsm2 <-function (x, trunc, start.value, upper = length(x), ...) {
    dots <- list(...)
    if (any(x <= 0) | any(!sads:::is.wholenumber(x))) stop("All x must be positive integers")
    if (sum(x) < 100) warning("\\n small sample size (J<100); \\n mzsm may not be a good approximation")
    if (!missing(trunc)) {
        if (min(x) <= trunc) stop("truncation point should be lower than the lowest data value")
    }
    if (missing(start.value)) thetahat <- length(x)
    else thetahat <- start.value
    if (missing(trunc)) {
        LL <- function(J, theta) -sum(dmzsm2(x, J = J, theta = theta, log = TRUE))
    }
    else {
        LL <- function(J, theta) -sum(dtrunc2("mzsm", x = x, coef = list(J = J, theta = theta), trunc = trunc, log = TRUE))
    }
    result <- do.call("mle2", c(list(LL, start = list(theta = thetahat), 
        fixed = list(J = sum(x)), data = list(x = x), method = "Brent", 
        lower = 0.001, upper = upper), dots))
    if (abs(as.numeric(result@coef) - upper) < 1e-07) warning("mle equal to upper bound provided. \\n Try value for the 'upper' arguent")
    new("fitsad", result, sad = "mzsm", trunc = ifelse(missing(trunc), NaN, trunc))
}

dmzsm2<-function (x, J, theta, log = FALSE) {
    J[!is.finite(J) | J <= 0] <- NaN
    theta[!is.finite(theta) | theta <= 0] <- NaN
    mzsm <- function(y, J, theta) log(theta)-log(y) + (theta - 1)*log(1 - y/J)
    sn <- mzsm(y = x, J = J, theta = theta)
    mu <- mzsm(y = 1:J, J = J, theta = theta)
    lpn <- suppressWarnings(sn - log(sum(exp(mu))))
    if (any(!sads:::is.wholenumber(x))) warning("non integer values in x")
    lpn[x <= 0 | x > J | !sads:::is.wholenumber(x)] <- -Inf
    if (any(is.nan(lpn))) warning("NaNs produced")
    if (log) return(lpn)
    else return(exp(lpn))
}



fitgamma2<- function (x, trunc, start.value, ...) {
    dots <- list(...)
    if (any(x <= 0)) stop("All x must be positive")
    if (!missing(trunc)) {
        if (min(x) <= trunc) stop("truncation point should, be lower than the lowest data value")
    }
    if (missing(start.value)) {
        if (missing(trunc)) {
            ka <- (mean(x)/sd(x))^2
            theta <- var(x)/mean(x)
            kahat <- function(k, dados) {
                eq <- length(dados) * (log(k) - log(mean(dados)) - digamma(k)) + sum(log(dados))
            }
            ka <- uniroot(kahat, interval = c(min(theta, ka), max(theta, ka)), dados = x)$root
            theta <- mean(x)/ka
        }
        else {
					  xh <- hist(x, plot = FALSE)
            xbr <- xh$breaks
            Eh <- matrix(ncol = 2, nrow = length(xbr) - 1)
            for (i in 1:(length(xbr) - 1)) {
                m1 <- matrix(c(1, 1, -1, 1), ncol = 2)
                Eh[i, ] <- solve(m1, xbr[i:(i + 1)])
            }
            P <- xh$counts/sum(xh$counts)
            P[P == 0] <- min(P[P > 0])
            Y <- log(P[-length(P)]) - log(P[-1]) - (log(Eh[-nrow(Eh), 
                2]) - log(Eh[-1, 2]))
            X1 <- Eh[-1, 1] - Eh[-nrow(Eh), 1]
            X2 <- log(Eh[-nrow(Eh), 1]) - log(Eh[-1, 1])
            st1 <- unname(coef(lm(Y ~ X1 + X2 - 1)))
            ka <- st1[2] + 1
            theta <- 1/st1[1]
        }
    }
    else {
        ka <- start.value[1]
        theta <- start.value[2]
    }
    if (missing(trunc)) {
        LL <- function(shape, rate) -sum(dgamma(x, shape, rate, log = TRUE))
    }
    else {
        LL <- function(shape, rate) -sum(dtrunc2("gamma", x = x, coef = list(shape = shape, rate = rate), trunc = trunc, log = TRUE))
    }
    result <- do.call("mle2", c(list(LL, start = list(shape = ka, rate = 1/theta), data = list(x = x)),lower=1e-12,method='L-BFGS-B', dots))
    new("fitsad", result, sad = "gamma", trunc = ifelse(missing(trunc), NaN, trunc))
}

fitgeom2<-function (x, trunc = 0, start.value, ...) {
    dots <- list(...)
    if (any(x <= 0) | any(!sads:::is.wholenumber(x))) stop("All x must be positive integers")
    if (!is.null(trunc)) {
        if (min(x) <= trunc) stop("truncation point should be lower than the lowest data value")
    }
    if (missing(start.value)) phat <- 1/(mean(x))
    else phat <- start.value
    if (!"method" %in% names(dots)) {
        dots$method <- "Brent"
        if (!"lower" %in% names(dots)) dots$lower = max(c(phat/10, 1e-08))
        if (!"upper" %in% names(dots)) dots$upper = min(c(phat * 10, 0.99))
    }
    if (is.null(trunc)) {
        LL <- function(prob) -sum(dgeom(x, prob, log = TRUE))
    }
    else {
        LL <- function(prob) -sum(dtrunc2("geom", x = x, coef = prob, trunc = trunc, log = TRUE))
    }
    result <- do.call("mle2", c(list(LL, start = list(prob = phat), data = list(x = x)), dots))
    new("fitsad", result, sad = "geom", trunc = ifelse(is.null(trunc), NaN, trunc))
}


dtrunc2<-function (f, x, trunc, coef, log = FALSE) {
    pf <- get(paste("p", f, sep = ""), mode = "function")
    df <- get(paste("d", f, sep = ""), mode = "function")
    tt <- rep(0, length(x))
    if (!missing(trunc)) tt[x > trunc] <- do.call(df, c(list(x = x[x > trunc],log=TRUE), coef))-log(1 - do.call(pf, c(list(q = trunc), coef)))
    else tt <- do.call(df, c(list(x = x), coef,log=TRUE))
    if (log) return(tt)
    return(exp(tt))
}

