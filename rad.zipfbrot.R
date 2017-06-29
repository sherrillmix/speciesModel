library(vegan)
#https://github.com/vegandevs/vegan/blob/master/R/rad.zipfbrot.R
"rad.zipfbrot" <- function (x, family = poisson, ...) {
    mandelfun <- function(p, x, ...) {
        brnk <- log(rnk + exp(p))
        if(any(is.infinite(brnk)))return(Inf)
        sol <- glm(x ~ brnk + offset(off), family = family(link = "log"))
        -logLik(sol)
    }
    x <- as.rad(x)
    rnk <- seq(along = x)
    off <- rep(log(sum(x)), length(x))
    ps <- c(0,1,-1,10,-10,.1,-.1)
    fam <- family(link = "log")
    if (length(x) > 2){
      nls <- lapply(ps,function(p)try(nlm(mandelfun, p = p, x = x, rnk = rnk, off = off, 
                      family = fam, hessian = TRUE, ...),silent=TRUE))
      scores<-sapply(nls,function(xx){
        if(inherits(xx,'try-error'))return(Inf)
        xx$minimum
      })
      nl<-nls[[which.min(scores)]]
    }
    if (length(x) < 3) {
        aic <- NA
        dev <- rdf <-  0
        ln <- nl <- NA
        p <- rep(NA, 3)
        fit <- x
        res <- rep(0, length(x))
        wts <- rep(1, length(x))
    }
    else if (inherits(nl, "try-error")) {
        aic <- rdf <- ln <- nl <- dev <-  NA
        p <- rep(NA, 3)
        fit <- res <- wts <- rep(NA, length(x))
    }
    else {
        ln <- glm(x ~ log(rnk + exp(nl$estimate)) + offset(off), 
                  family = family(link = "log"))
        fit <- fitted(ln)
        p <- c(coef(ln), exp(nl$estimate))
        p[1] <- exp(p[1])
        aic <- AIC(ln) + 2
        rdf <- df.residual(ln) - 1
        dev <- deviance(ln)
        res <- ln$residuals
        wts <- weights(ln)
    }
    names(p) <- c("c", "gamma", "beta")
    out <- list(model = "Zipf-Mandelbrot", family = fam, 
                y = x, coefficients = p, fitted.values = fit, aic = aic, 
                rank = 3, df.residual = rdf, deviance = dev, 
                residuals = res, prior.weights = wts)
    class(out) <- c("radline", "glm")
    out
}

#http://r.789695.n4.nabble.com/how-to-override-replace-a-function-in-a-package-namespace-td866337.html
unlockBinding("rad.zipfbrot", as.environment('package:vegan'));
assignInNamespace('rad.zipfbrot',rad.zipfbrot,ns='vegan',envir='package:vegan')
assign('rad.zipfbrot',rad.zipfbrot,envir=as.environment('package:vegan'))
lockBinding("rad.zipfbrot", as.environment('package:vegan'))
