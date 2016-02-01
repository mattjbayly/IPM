#A single object matching ‘BCa’ was found
#It was found in the following places
#  package:bootBCa
#  namespace:bootBCa
#with value



BCa_validate_theta <- function (x) 
{
    if (length(x) != 1) 
        stop("BCa: theta returned a non-scalar value")
    if (!is.finite(x)) 
        stop("BCa: theta returned a non-finite value")
    if (!is.numeric(x)) 
        stop("BCa: theta returned a non-numeric value")
    x
}


BCa_fn <- function (thetastar, thetahat, acc, zalpha, qtype) 
{
    nboot <- length(thetastar)
    z0 <- qnorm(sum(thetastar < thetahat)/nboot)
    if (any(acc * (z0 + zalpha) >= 1)) 
        stop("BCa: alpha value too extreme for these data")
    tt <- pnorm(z0 + (z0 + zalpha)/(1 - acc * (z0 + zalpha)))
    c(mean(thetastar), quantile(x = thetastar, probs = tt, type = qtype))
}


BCa_ret <- function (v, alpha) 
{
    names(v) <- c("nboot", "prec", "est", sapply(alpha, function(x) sprintf("%0.3f", 
        x)))
    v
}



BCa <- function (x, delta, theta, ..., alpha = c(0.025, 0.975), verbose = F, 
    M = 25000, Mlimit = 1500000, qtype = 1) 
{
    if (is.function(x)) 
        stop("BCa: x must not be a function")
    if (is.expression(x)) 
        stop("BCa: x must not be an expression")
    if (is.factor(x)) 
        stop("BCa: x must not be a factor")
    if (length(x) == 0) 
        stop("BCa: x must have nonzero length")
    if (length(x) < 2 && is.na(x)) 
        stop("BCa: x must not be NA")
    if (is.function(delta)) 
        stop("BCa: delta must not be a function")
    if (is.expression(delta)) 
        stop("BCa: delta must not be an expression")
    if (length(delta) != 1) 
        stop("BCa: delta value must be scalar")
    if (!is.na(delta) && !is.numeric(delta)) 
        stop("BCa: delta value must be numeric")
    if (is.nan(delta)) 
        stop("BCa: delta value must not be NaN")
    if (!is.function(theta)) 
        stop("BCa: theta must be a function")
    if (is.function(alpha)) 
        stop("BCa: alpha must not be a function")
    if (is.expression(alpha)) 
        stop("BCa: alpha must not be an expression")
    if (length(alpha) == 0) 
        stop("BCa: alpha must have nonzero length")
    if (length(alpha) < 2 && is.na(alpha)) 
        stop("BCa: alpha must not be NA")
    if (any(!is.finite(alpha))) 
        stop("BCa: alpha values must be finite")
    if (any(!is.numeric(alpha))) 
        stop("BCa: alpha values must be numeric")
    if (any(alpha <= 0 | alpha >= 1)) 
        stop("BCa: alpha values must lie between 0 and 1 exclusive")
    if (is.function(verbose)) 
        stop("BCa: verbose must not be a function")
    if (is.expression(verbose)) 
        stop("BCa: verbose must not be an expression")
    if (length(verbose) != 1) 
        stop("BCa: verbose must be scalar")
    if (is.na(verbose)) 
        stop("BCa: verbose must not be NA")
    if (!is.logical(verbose)) 
        stop("BCa: verbose must be T or F")
    if (is.function(M)) 
        stop("BCa: M must not be a function")
    if (is.expression(M)) 
        stop("BCa: M must not be an expression")
    if (length(M) != 1) 
        stop("BCa: M must be scalar")
    if (!is.finite(M)) 
        stop("BCa: M must be finite")
    if (!is.numeric(M)) 
        stop("BCa: M must be numeric")
    if (M < 2) 
        stop("BCa: M must be greater than 1")
    if (is.function(Mlimit)) 
        stop("BCa: Mlimit must not be a function")
    if (is.expression(Mlimit)) 
        stop("BCa: Mlimit must not be an expression")
    if (length(Mlimit) != 1) 
        stop("BCa: Mlimit must be scalar")
    if (is.na(Mlimit)) 
        stop("BCa: Mlimit must not be NA")
    if (!is.numeric(Mlimit)) 
        stop("BCa: Mlimit must be numeric")
    if (!is.na(delta) && Mlimit < 2 * M) 
        stop("BCa: Mlimit cannot be less than 2M")
    if (is.function(qtype)) 
        stop("BCa: qtype must not be a function")
    if (is.expression(qtype)) 
        stop("BCa: qtype must not be an expression")
    if (length(qtype) != 1) 
        stop("BCa: qtype must be scalar")
    if (!is.finite(qtype)) 
        stop("BCa: qtype must be finite")
    if (!is.numeric(qtype)) 
        stop("BCa: qtype must be numeric")
    if (qtype < 1 || qtype > 9) 
        stop("BCa: qtype must be between 1 and 9 inclusive")
    n <- length(x)
    if (verbose) {
        if (is.na(delta)) 
            message(sprintf("BCa call n=%d, non-adaptive, fixed M=%d", 
                n, M))
        else message(sprintf("BCa call n=%d, delta=%f, M=%d, Mlimit=%d", 
            n, delta, M, Mlimit))
    }
    thetahat <- theta(x, ...)
    if (length(thetahat) != 1) 
        stop("BCa: theta returned a non-scalar value")
    if (length(unique(x)) == 1) {
        if (verbose) 
            warning("BCa: handling degenerate case (all values are the same).")
        return(BCa_ret(c(0, 0, rep(thetahat, 1 + length(alpha))), 
            alpha))
    }
    BCa_validate_theta(thetahat)
    u <- rep(0, n)
    for (i in 1:n) u[i] <- BCa_validate_theta(theta(x[-i], ...))
    uu <- mean(u) - u
    if (all(uu == 0)) 
        stop("BCa: acceleration is indeterminate for insensitive function")
    acc <- sum(uu * uu * uu)/(6 * (sum(uu * uu))^1.5)
    if (verbose) 
        message(sprintf("acceleration=%f", acc))
    zalpha <- qnorm(alpha)
    thetastar <- rep(NA, M)
    for (i in 1:M) thetastar[i] <- BCa_validate_theta(theta(sample(x, 
        size = n, replace = T), ...))
    answers <- BCa_fn(thetastar, thetahat, acc, zalpha, qtype)
    if (is.na(delta)) 
        return(BCa_ret(c(M, NA, answers), alpha))
    h <- 1
    u <- delta
    allthetastar <- thetastar
    while (u >= delta & h * M < Mlimit) {
        for (i in 1:M) thetastar[i] <- BCa_validate_theta(theta(sample(x, 
            size = n, replace = T), ...))
        answers <- rbind(answers, BCa_fn(thetastar, thetahat, 
            acc, zalpha, qtype))
        allthetastar <- rbind(allthetastar, thetastar)
        h <- h + 1
        u <- 2 * max(apply(answers, 2, sd))/sqrt(h)
        if (verbose) 
            message(sprintf("u=%f nboot=%d delta=%f", u, h * 
                M, delta))
    }
    BCa_ret(c(h * M, u, BCa_fn(allthetastar, thetahat, acc, zalpha, 
        qtype)), alpha)
}