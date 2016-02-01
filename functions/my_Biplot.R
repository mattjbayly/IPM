myBiplot <- function (x, choices = 1L:2L, scale = 1,
                      pc.biplot = FALSE, var.axes = TRUE,
                      type = "t",
                      col,
                      col.arrows = "#FF0000",
                      col.text = "#000000",
                      cex = rep(par("cex"), 2),
                      expand = 1, 
                      xlabs = NULL, ylabs = NULL,
                      xlim = NULL, ylim = NULL, 
                      main = NULL, sub = NULL,
                      xlab = NULL, ylab = NULL, 
                      arrow.len = 0.1,
                      ...
                      )

{
	if (length(choices) != 2L) 
        stop("length of choices must be 2")
    if (!length(scores <- x$x)) 
        stop(gettextf("object '%s' has no scores", deparse(substitute(x))), 
            domain = NA)
    if (is.complex(scores)) 
        stop("biplots are not defined for complex PCA")
        
    lam <- x$sdev[choices]
    n <- NROW(scores)
    lam <- lam * sqrt(n)
    
    if (scale < 0 || scale > 1) 
        warning("'scale' is outside [0, 1]")
    if (scale != 0) 
        lam <- lam^scale
    else lam <- 1
    if (pc.biplot) 
        lam <- lam/sqrt(n)
        
    y <- t(t(x$rotation[, choices]) * lam)
    x <- t(t(scores[, choices])/lam)  # note that from here on
                                      # x is no longer the PC object
                                      # originally pased into the function
    n <- nrow(x)
    p <- nrow(y)
    
    if (missing(xlabs)) {
        xlabs <- dimnames(x)[[1L]]
        if (is.null(xlabs)) 
            xlabs <- 1L:n
    }
    xlabs <- as.character(xlabs)
    dimnames(x) <- list(xlabs, dimnames(x)[[2L]])
    
    if (missing(ylabs)) {
        ylabs <- dimnames(y)[[1L]]
        if (is.null(ylabs)) 
            ylabs <- paste("Var", 1L:p)
    }
    ylabs <- as.character(ylabs)
    dimnames(y) <- list(ylabs, dimnames(y)[[2L]])

    if (length(cex) == 1L) 
        cex <- c(cex, cex)

    unsigned.range <- function(x) c(-abs(min(x, na.rm = TRUE)), 
        abs(max(x, na.rm = TRUE)))
    rangx1 <- unsigned.range(x[, 1L])
    rangx2 <- unsigned.range(x[, 2L])
    rangy1 <- unsigned.range(y[, 1L])
    rangy2 <- unsigned.range(y[, 2L])

    if (missing(xlim) && missing(ylim)) 
        xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
    else if (missing(xlim)) 
        xlim <- rangx1
    else if (missing(ylim)) 
        ylim <- rangx2

    ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
    on.exit(par(op))
    op <- par(pty = "s")
    if (!is.null(main)) 
        op <- c(op, par(mar = par("mar") + c(0, 0, 1, 0)))

    # first, plot scores - either normally, or as row labels
    if (type == "p") {
	    plot(x, type = type, xlim = xlim, ylim = ylim, col = col, 
        xlab = xlab, ylab = ylab, sub = sub, main = main, ...)
    }
    else if (type == "t") {
        plot(x, type = "n", xlim = xlim, ylim = ylim, 
             xlab = xlab, ylab = ylab, sub = sub, main = main, ...)
        text(x, xlabs, cex = cex[1L], col = col.text, ...)
    	
    }
    else if (type == "n") {  # plot an empty frame
	    plot(x, type = type, xlim = xlim, ylim = ylim, 
        xlab = xlab, ylab = ylab, sub = sub, main = main, ...)
    }

    par(new = TRUE)
    dev.hold()
    on.exit(dev.flush(), add = TRUE)
    plot(y, axes = FALSE, type = "n", xlim = xlim * ratio, ylim = ylim * 
        ratio, xlab = "", ylab = "", col = col.arrows, ...)
    axis(3, col = col.arrows, ...)
    axis(4, col = col.arrows, ...)
    box(col = "#000000")
    text(y, labels = ylabs, cex = cex[2L], col = col.arrows, ...)
    if (var.axes) 
        arrows(0, 0, y[, 1L] * 0.8, y[, 2L] * 0.8, col = col.arrows, 
            length = arrow.len)
    # now replot into xlim, ylim scaled by lam, to reset par("usr")
    # to the  correct values needed for subsequent application of points(),
    # text() etc.
    par(new = TRUE)
    dev.hold()
    on.exit(dev.flush(), add = TRUE)
	plot(0, type = "n", xlim = xlim * lam[1], ylim = ylim * lam[2], 
    xlab = '', ylab = '', sub = '', main = '', xaxt='n', yaxt='n', axes=FALSE)

    invisible()

}
