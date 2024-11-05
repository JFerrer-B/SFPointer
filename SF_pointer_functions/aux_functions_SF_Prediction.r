#' @rdname InternalFunctions
getpij <- function(A){
    rSc <- rowSums(A)
    cSc <- colSums(A)
    rSums <- factor(rSc)
    cSums <- factor(cSc)
    
    Xbig2small <- sparseMatrix(i = 1:(nrow(A)+ncol(A)),
                               j = c(as.numeric(rSums), as.numeric(cSums)+length(levels(rSums))),
                               x =1)
    
    rM <- sparseMatrix(i = as.numeric(rSums),
                       j = 1:nrow(A),
                       x =1)
    cM <- sparseMatrix(i = 1:ncol(A),
                       j = as.numeric(cSums),
                       x =1)
    
    A1 <- rM %*% A %*% cM 
    A0 <- rM %*% (1-A) %*% cM 
    
    X1 <- sparseMatrix(i = 1:(nrow(A1)*ncol(A1)),
                       j = rep(1:nrow(A1), ncol(A1)),
                       x = 1,
                       dims = c((nrow(A1)*ncol(A1)), (nrow(A1)+ncol(A1))))
    
    X2 <- sparseMatrix(i = 1:(nrow(A1)*ncol(A1)),
                       j = nrow(A1)+rep(1:ncol(A1), each = nrow(A1)),
                       x = 1,
                       dims = c((nrow(A1)*ncol(A1)), (nrow(A1)+ncol(A1))))
    Xtxiki <- X1 + X2 # Matriz de diseño
    Xtxiki <- rbind(Xtxiki, c(rep(1,nrow(A1)), rep(0,ncol(A1))))
    Salida5 <- speedglm.wfit2(X=(Xtxiki),
                              y=cbind(c(as.vector(A1),round(mean(A1))), c(as.vector(A0),round(mean(A0)))),
                              family = binomial(), sparse=TRUE)
    
    cS5A <- coef(Salida5)
    cS5 <- cS5A %*% t(Xbig2small)
    elog5 <- matrix(exp(cS5[1:nrow(A)]), ncol = 1) %*% matrix(exp(cS5[nrow(A) + 1:ncol(A)]), nrow = 1)
    p5 <- elog5 /(1+elog5)
    
    ## TODO: include a warning if the colSums (rowSums) of p5 is
    # not similar to the colSums (rowSums) of A
    deltar <- rowSums(p5) - rSc
    if (max(abs(deltar))>1e-3) warning("Converge problems in the solution!!!")
    rownames(p5) <- rownames(A)
    colnames(p5) <- colnames(A)
    return(p5)
}

#' @rdname InternalFunctions
speedglm.wfit2 <- function (y, X, intercept = TRUE, weights = NULL, row.chunk = NULL, 
                            family = gaussian(), start = NULL, etastart = NULL, mustart = NULL, 
                            offset = NULL, acc = 1e-08, maxit = 25, k = 2, sparselim = 0.9, 
                            camp = 0.01, eigendec = TRUE, tol.values = 1e-07, tol.vectors = 1e-07, 
                            tol.solve = .Machine$double.eps, sparse = NULL, method = c("eigen", 
                                                                                       "Cholesky", "qr"), trace = FALSE, ...) 
{
    nobs <- NROW(y)
    nvar <- ncol(X)
    if (missing(y)) 
        stop("Argument y is missing")
    if (missing(X)) 
        stop("Argument X is missing")
    if (is.null(offset)) 
        offset <- rep.int(0, nobs)
    if (is.null(weights)) 
        weights <- rep(1, nobs)
    col.names <- dimnames(X)[[2]]
    method <- match.arg(method)
    fam <- family$family
    link <- family$link
    variance <- family$variance
    dev.resids <- family$dev.resids
    aic <- family$aic
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta
    if (is.null(sparse)) 
        sparse <- is.sparse(X = X, sparselim, camp)
    if (is.null(start)) {
        if (is.null(mustart)) 
            eval(family$initialize)
        eta <- if (is.null(etastart)) 
            family$linkfun(mustart)
        else etastart
        mu <- mustart
        start <- rep(0, nvar)
    }
    else {
        eta <- offset + as.vector(if (nvar == 1) 
            X * start
            else {
                if (sparse) 
                    X %*% start
                else tcrossprod(X, t(start))
            })
        mu <- linkinv(eta)
    }
    iter <- 0
    dev <- sum(dev.resids(y, mu, weights))
    tol <- 1
    if ((fam == "gaussian") & (link == "identity")) 
        maxit <- 1
    C_Cdqrls <- getNativeSymbolInfo("Cdqrls", PACKAGE = getLoadedDLLs()$stats)
    while ((tol > acc) & (iter < maxit)) {
        iter <- iter + 1
        beta <- start
        dev0 <- dev
        varmu <- variance(mu)
        mu.eta.val <- mu.eta(eta)
        z <- (eta - offset) + (y - mu)/mu.eta.val
        W <- (weights * mu.eta.val * mu.eta.val)/varmu
        X1 <- sqrt(W) * X
        XTX <- crossprod(X1)
        XTz <- t(crossprod((W * z), X))
        if (iter == 1 & method != "qr") {
            variable <- colnames(X)
            ris <- if (eigendec)
                control(XTX, , tol.values, tol.vectors, , method)
            else list(rank = nvar, pivot = 1:nvar)
            ok <- ris$pivot[1:ris$rank]
            if (eigendec) {
                XTX <- ris$XTX
                X <- X[, ok]
                XTz <- XTz[ok]
                start <- start[ok]
            }
            beta <- start
        }
        if (method == "qr") {
            ris <- .Call(C_Cdqrls, XTX, XTz, tol.values, FALSE)
            start <- if (ris$rank < nvar) 
                ris$coefficients[ris$pivot]
            else ris$coefficients
        }
        else {
            start <- solve(XTX, XTz, tol = tol.solve)
        }
        eta <-  drop(X %*% start)
        mu <- linkinv(eta <- eta + offset)
        dev <- sum(dev.resids(y, mu, weights))
        tol <- max(abs(dev0 - dev)/(abs(dev) + 0.1))
        if (trace) 
            cat("iter", iter, "tol", tol, "\n")
    }
    wt <- sum(weights)
    wtdmu <- if (intercept) 
        sum(weights * y)/wt
    else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- ris$rank
    dfr <- nobs - rank - sum(weights == 0)
    aic.model <- aic(y, nobs, mu, weights, dev) + k * rank
    #ll.nuovo <- speedglm:::ll.speedglm(fam, aic.model, rank)
    res <- (y - mu)/mu.eta(eta)
    resdf <- n.ok - rank
    RSS <- sum(W * res * res)
    var_res <- RSS/dfr
    dispersion <- if (fam %in% c("poisson", "binomial")) 
        1
    else var_res
    if (method == "qr") {
        coefficients <- start
        coefficients[coefficients == 0] = NA
        ok <- ris$pivot[1:rank]
    }
    else {
        coefficients <- rep(NA, nvar)
        start <- as(start, "numeric")
        coefficients[ok] <- start
    }
    names(coefficients) <- col.names
    rval <- list(coefficients = coefficients, 
                 iter = iter, tol = tol, family = family, link = link, 
                 df = dfr, XTX = XTX, dispersion = dispersion, ok = ok, 
                 rank = rank, RSS = RSS, method = method, aic = aic.model, 
                 sparse = sparse, deviance = dev, nulldf = nulldf, nulldev = nulldev, 
                 ngoodobs = n.ok, n = nobs, intercept = intercept, convergence = (!(tol > 
                                                                                        acc)))
    class(rval) <- "speedglm"
    rval
}


myphyper <- function(p, m, n, k,lower.tail=TRUE,log.p = FALSE) {
    #enrichment = p1 > p2, one - sided test, critical region right
    #this is the lower.tail = FALSE
    if(lower.tail==FALSE){
        pp <- phyper(p+1, m, n, k,lower.tail = FALSE) + .5 * dhyper(p, m, n, k)
    }
    
    #depletion = p1 < p2 one - sided test, critical region left
    #this is the lower.tail = TRUE
    if(lower.tail==TRUE){
        pp <- phyper(p-1, m, n, k) + .5 * dhyper(p, m, n, k)
    }
    
    if(log.p == TRUE){
        pp <- log(pp)
    }
    
    return(pp)
}