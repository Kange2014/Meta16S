#! Plot Mothur's PCoA Data with an optional design file for groups in 3d
#!
#! input: axes: the axes file from mothur's pcoa function.
#!		  loading: the loading file from mothur's pcoa function.
#!		  design: the design file for the groups you which to color the points by. defaults to no design file.  
#!		  ax1, ax2, axes: the axes you wish to plot. Defaults to the 3 most relevant according to mothur.
#!		  Lab1, Lab2, Label2: titles for the 1st/x, 2nd/y and 3rd/z axis.
#!		  Pch: plotting "character", i.e. symbol to use. Defaults to 19 (solid circles)
#!		  Main: an overall title for the plot.
#!		  Sub: sub-title.
#!		  Scale.y: scale of y axis related to x- and z axis.
#!		  Angle: angle between x and y axis (Attention: result depends on scaling).
#!		  Axis: a logical value indicating whether axes should be drawn on the plot.
#!		  Tick.marks:a logical value indicating whether tick marks should be drawn on the plot (only if axis = TRUE).
#!		  Label.tick.marks: a logical value indicating whether tick marks should be labeled on the plot (only if axis = TRUE and tick.marks = TRUE).
#!		  X.ticklabs, Y.ticklabs, Z.ticklabs: vector of tick mark labels.
#!		  Y.margin.add: add additional space between tick mark labels and axis label of the y axis.
#!		  Grid: a logical value indicating whether a grid should be drawn on the plot.
#!		  Box: a logical value indicating whether a box should be drawn around the plot.
#!		  Lab: a numerical vector of the form c(x, y, len). The values of x and y give the (approximate) number of tickmarks on the x and y axes.
#!		  Lab.z: the same as lab, but for z axis.
#!		  Type: character indicating the type of plot: "p" for points, "l" for lines, "h" for vertical lines to x-y-plane, etc.
#!		  LX, LY: the x and y co-ordinates to be used to position the legend. They can be specified by keyword or in any way which is accepted by xy.coords.

# Check if required packages are already installed, and install if missing
#packages <-c("scatterplot3d") 

# Function to check whether the package is installed
#InsPack <- function(pack)
#{
#  if ((pack %in% installed.packages()) == FALSE) {
#    install.packages(pack, destdir="/results/plugins/Meta16S/scripts", lib=getwd())
#  } 
#}

scatterplot3d <-
function(x, y = NULL, z = NULL, color = par("col"), pch = par("pch"),
     main = NULL, sub = NULL, xlim = NULL, ylim = NULL, zlim = NULL,
     xlab = NULL, ylab = NULL, zlab = NULL, scale.y = 1, angle = 40,
     axis = TRUE, tick.marks = TRUE, label.tick.marks = TRUE,
     x.ticklabs = NULL, y.ticklabs = NULL, z.ticklabs = NULL,
     y.margin.add = 0, grid = TRUE, box = TRUE, lab = par("lab"),
     lab.z = mean(lab[1:2]), type = "p", highlight.3d = FALSE,
     mar = c(5, 3, 4, 3) + 0.1, bg = par("bg"), col.axis = par("col.axis"),
     col.grid = "grey", col.lab = par("col.lab"), cex.symbols = par("cex"),
     cex.axis = 0.8 * par("cex.axis"), cex.lab = par("cex.lab"),
     font.axis = par("font.axis"), font.lab = par("font.lab"),
     lty.axis = par("lty"), lty.grid = par("lty"), lty.hide=NULL,
     lty.hplot = par("lty"), log = "", asp = NA, ...)
     # log not yet implemented
{
    ## Uwe Ligges <ligges@statistik.tu-dortmund.de>,
    ## http://www.statistik.tu-dortmund.de/~ligges
    ##
    ## For MANY ideas and improvements thanks to Martin Maechler!!!
    ## Parts of the help files are stolen from the standard plotting functions in R.

    mem.par <- par(mar = mar)
#    on.exit(par(mem.par))
    x.scal <- y.scal <- z.scal <- 1
    xlabel <- if (!missing(x)) deparse(substitute(x))
    ylabel <- if (!missing(y)) deparse(substitute(y))
    zlabel <- if (!missing(z)) deparse(substitute(z))
    ## verification, init, ...
    if(highlight.3d && !missing(color))
        warning("color is ignored when highlight.3d = TRUE")


    ## color as part of `x' (data.frame or list):
    if(!is.null(d <- dim(x)) && (length(d) == 2) && (d[2] >= 4))
        color <- x[,4]
    else if(is.list(x) && !is.null(x$color))
        color <- x$color

    ## convert 'anything' -> vector
    xyz <- xyz.coords(x=x, y=y, z=z, xlab=xlabel, ylab=ylabel, zlab=zlabel,
                      log=log)
    if(is.null(xlab)) { xlab <- xyz$xlab; if(is.null(xlab)) xlab <- "" }
    if(is.null(ylab)) { ylab <- xyz$ylab; if(is.null(ylab)) ylab <- "" }
    if(is.null(zlab)) { zlab <- xyz$zlab; if(is.null(zlab)) zlab <- "" }

    if(length(color) == 1)
        color <- rep(color, length(xyz$x))
    else if(length(color) != length(xyz$x))
        stop("length(color) ", "must be equal length(x) or 1")

    if(length(pch) == 1)
        pch <- rep(pch, length(xyz$x))
    else if(length(pch) != length(xyz$x))
        stop("length(pch) ", "must be equal length(x) or 1")

    if(length(bg) == 1)
        bg <- rep(bg, length(xyz$x))
    else if(length(bg) != length(xyz$x))
        stop("length(bg) ", "must be equal length(x) or 1")


    angle <- (angle %% 360) / 90
    yz.f <- scale.y * abs(if(angle < 1) angle else if(angle > 3) angle - 4 else 2 - angle)
    yx.f <- scale.y * (if(angle < 2) 1 - angle else angle - 3)
    if(angle > 2) { ## switch y and x axis to ensure righthand oriented coord.
        temp <- xyz$x; xyz$x <- xyz$y; xyz$y <- temp
        temp <- xlab;  xlab <- ylab;   ylab <- temp
        temp <- xlim;  xlim <- ylim;   ylim <- temp
    }
    angle.1 <- (1 < angle && angle <= 2) || angle > 3
    angle.2 <- 1 < angle && angle <= 3
    dat <- data.frame(xyz[c("x","y","z")], col = color, pch = pch, bg = bg, stringsAsFactors = FALSE)
    ## xlim, ylim, zlim -- select the points inside the limits
    if(!is.null(xlim)) {
        xlim <- range(xlim)
        dat <- dat[ xlim[1] <= dat$x & dat$x <= xlim[2] , , drop = FALSE]
    }
    if(!is.null(ylim)) {
        ylim <- range(ylim)
        dat <- dat[ ylim[1] <= dat$y & dat$y <= ylim[2] , , drop = FALSE]
    }
    if(!is.null(zlim)) {
        zlim <- range(zlim)
        dat <- dat[ zlim[1] <= dat$z & dat$z <= zlim[2] , , drop = FALSE]
    }
    n <- nrow(dat)
    if(n < 1) stop("no data left within (x|y|z)lim")

    y.range <- range(dat$y[is.finite(dat$y)])

### 3D-highlighting / colors / sort by y
    if(type == "p" || type == "h") {
        y.ord <- rev(order(dat$y))
        dat <- dat[y.ord, ]
        if(length(cex.symbols) > 1)
            if(length(cex.symbols) != length(y.ord))
                stop("length(cex.symbols) ", "must be equal length(x) or 1")
            else cex.symbols <- cex.symbols[y.ord]
        daty <- dat$y
        daty[!is.finite(daty)] <- mean(daty[is.finite(daty)])
        if(highlight.3d && !(all(diff(daty) == 0)))
            dat$col <- rgb(red=seq(0, 1, length = n) * (y.range[2] - daty) / diff(y.range), green=0, blue=0)
    }

### optim. axis scaling
    p.lab <- par("lab")
    ## Y
    y.range <- y.range.fix <- range(dat$y[is.finite(dat$y)], ylim)
    y.prty <- pretty(y.range, n = lab[2],
        min.n = max(1, min(.5 * lab[2], p.lab[2])))
    y.scal <- round(diff(y.prty[1:2]), digits = 12)
    y.add <- min(y.prty)
    dat$y <- (dat$y - y.add) / y.scal
    y.max <- (max(y.prty) - y.add) / y.scal
    if(!is.null(ylim)) y.max <- max(y.max, ceiling((ylim[2] - y.add) / y.scal))
#    if(angle > 2) dat$y <- y.max - dat$y  ## turn y-values around
    ## X
    x.range <- x.range.fix <- range(dat$x[is.finite(dat$x)], xlim)
    x.prty <- pretty(x.range, n = lab[1],
        min.n = max(1, min(.5 * lab[1], p.lab[1])))
    x.scal <- round(diff(x.prty[1:2]), digits = 12)
    dat$x <- dat$x / x.scal
    x.range <- range(x.prty) / x.scal
    x.max <- ceiling(x.range[2])
    x.min <-   floor(x.range[1])
    if(!is.null(xlim)) {
        x.max <- max(x.max, ceiling(xlim[2] / x.scal))
        x.min <- min(x.min,   floor(xlim[1] / x.scal))
    }
    x.range <- range(x.min, x.max)
    ## Z
    z.range <- range(dat$z[is.finite(dat$z)], zlim)
    z.prty <- pretty(z.range, n = lab.z,
        min.n = max(1, min(.5 * lab.z, p.lab[2])))
    z.scal <- round(diff(z.prty[1:2]), digits = 12)
    dat$z <- dat$z / z.scal
    z.range <- range(z.prty) / z.scal
    z.max <- ceiling(z.range[2])
    z.min <-   floor(z.range[1])
    if(!is.null(zlim)) {
        z.max <- max(z.max, ceiling(zlim[2] / z.scal))
        z.min <- min(z.min,   floor(zlim[1] / z.scal))
    }
    z.range <- range(z.min, z.max)

### init graphics

### convert asp for plot (based on suggestions from Jari Oksanen)
    if(!is.na(asp)) {
        x.i <- x.min:x.max
        z.i <- z.min:z.max
        range.x <- abs(diff(range(x.i * x.scal)))
        range.z <- abs(diff(range(z.i * z.scal)))
        asp <- asp * (range.z / (length(z.i) - 1)) / (range.x / (length(x.i) - 1))
    }
    plot.new()
    if(angle.2) {x1 <- x.min + yx.f * y.max; x2 <- x.max}
    else        {x1 <- x.min; x2 <- x.max + yx.f * y.max}
    plot.window(c(x1, x2), c(z.min, z.max + yz.f * y.max), asp = asp)
    temp <- strwidth(format(rev(y.prty))[1], cex = cex.axis/par("cex"))

### lheight in usr units for numeric aspect is needed to locate
### side 2 and 4 axis annotation with fixes aspect.
    lheight <- (strheight("\n") - strheight("M")) * asp
    lheight2 <- (strheight("\n") - strheight("M"))

    if(angle.2) x1 <- x1 - temp - y.margin.add
    else        x2 <- x2 + temp + y.margin.add
    plot.window(c(x1, x2), c(z.min, z.max + yz.f * y.max), asp = asp)
    if(angle > 2) par("usr" = par("usr")[c(2, 1, 3:4)])
    usr <- par("usr") # we have to remind it for use in closures
    title(main, sub, ...)

### draw axis, tick marks, labels, grid, ...
    if(grid) {
        ## X
        i <- x.min:x.max
        segments(i, z.min, i + (yx.f * y.max), yz.f * y.max + z.min,
                 col = col.grid, lty = lty.grid)
        ## Y
        i <- 0:y.max
        segments(x.min + (i * yx.f), i * yz.f + z.min,
                 x.max + (i * yx.f), i * yz.f + z.min,
                 col = col.grid, lty = lty.grid)
    }
    if(axis) {
        xx <- if(angle.2) c(x.min, x.max) else c(x.max, x.min)
        if(tick.marks) { ## tick marks
            xtl <- (z.max - z.min) * (tcl <- -par("tcl")) / 50
            ztl <- (x.max - x.min) * tcl / 50
            mysegs <- function(x0,y0, x1,y1)
                segments(x0,y0, x1,y1, col=col.axis, lty=lty.axis)
            ## Y
            i.y <- 0:y.max
            mysegs(yx.f * i.y - ztl + xx[1], yz.f * i.y + z.min,
                   yx.f * i.y + ztl + xx[1], yz.f * i.y + z.min)
            ## X
            i.x <- x.min:x.max
            mysegs(i.x, -xtl + z.min, i.x, xtl + z.min)
            ## Z
            i.z <- z.min:z.max
            mysegs(-ztl + xx[2], i.z, ztl + xx[2], i.z)

            if(label.tick.marks) { ## label tick marks
                las <- par("las")
                mytext <- function(labels, side, at, line = -0.5, ...)
                    mtext(text = labels, side = side, at = at, line = line,
                          col=col.lab, cex=cex.axis, font=font.lab, ...)
                ## X
                if(is.null(x.ticklabs))
                    x.ticklabs <- format(i.x * x.scal)
                if(!is.na(asp)) {
                    linepad <- (usr[3] - z.min)/lheight2 + 0.5
                    mytext(x.ticklabs, side = 1, at = i.x, line = linepad)
                } else {
                    mytext(x.ticklabs, side = 1, at = i.x)
                }
                ## Z
                if(is.null(z.ticklabs))
                    z.ticklabs <- format(i.z * z.scal)
                if(!is.na(asp)) {
                    if(angle.1) {
                        if(angle > 2) {
                            linepad <- (x2 - usr[1])/lheight + 0.5
                        } else {
                            linepad <- (x2 - usr[2])/lheight + 0.5
                        }
                    } else {
                        if(angle > 2) {
                            linepad <- (usr[2] - x1)/lheight + 0.5
                        } else {
                            linepad <- (usr[1] - x1)/lheight + 0.5
                        }
                   }
                } else {
                    linepad = -0.5
                }
                mytext(z.ticklabs, side = if(angle.1) 4 else 2, at = i.z,
                       adj = if(0 < las && las < 3) 1 else NA, line = linepad)
                ## Y
                temp <- if(angle > 2) rev(i.y) else i.y ## turn y-labels around
                if(is.null(y.ticklabs))
                    y.ticklabs <- format(y.prty)
                else if (angle > 2)
                    y.ticklabs <- rev(y.ticklabs)
                text(i.y * yx.f + xx[1],
                     i.y * yz.f + z.min, y.ticklabs,
                     pos=if(angle.1) 2 else 4, offset=1,
                     col=col.lab, cex=cex.axis/par("cex"), font=font.lab)
            }
        }

        ## axis and labels

        ## determine position of labels
        if(!is.na(asp)) {
            if(angle.1) {
                if(angle > 2) {
                    linepad <- (x2 - usr[1])/lheight + 0.5
                } else {
                    linepad <- (x2 - usr[2])/lheight + 0.5
                }
            } else {
                if(angle > 2) {
                    linepad <- (usr[2] - x1)/lheight + 0.5
                } else {
                    linepad <- (usr[1] - x1)/lheight + 0.5
                }
            }
        } else {
            linepad = -0.5
        }
        
        mytext2 <- function(lab, side, line, at)
            mtext(lab, side = side, line = line, at = at, col = col.lab,
                  cex = cex.lab, font = font.axis, las = 0)
        ## X
        lines(c(x.min, x.max), c(z.min, z.min), col = col.axis, lty = lty.axis)
        if(!is.na(asp)) {
            mytext2(xlab, 1, line = (usr[3] - z.min)/lheight2 + 1.5, at = mean(x.range))
        } else {
            mytext2(xlab, 1, line = 1.5, at = mean(x.range))
        }
        ## Y
        lines(xx[1] + c(0, y.max * yx.f), c(z.min, y.max * yz.f + z.min),
              col = col.axis, lty = lty.axis)
        mytext2(ylab, if(angle.1) 2 else 4, line = linepad + 1, at = z.min + y.max * yz.f)
        ## Z
        lines(xx[c(2,2)], c(z.min, z.max), col = col.axis, lty = lty.axis)
        mytext2(zlab, if(angle.1) 4 else 2, line = linepad + 2, at = mean(z.range))
        if(box) {
            if(is.null(lty.hide)) lty.hide <- lty.axis
            ## X
            temp <- yx.f * y.max
            temp1 <- yz.f * y.max
            lines(c(x.min + temp, x.max + temp),
                  c(z.min + temp1, z.min + temp1), col = col.axis, lty = lty.hide)
            lines(c(x.min + temp, x.max + temp), c(temp1 + z.max, temp1 + z.max),
                  col = col.axis, lty = lty.axis)
            ## Y
            temp <- c(0, y.max * yx.f)
            temp1 <- c(0, y.max * yz.f)
            lines(temp + xx[2], temp1 + z.min, col = col.axis, lty = lty.hide)
            lines(temp + x.min, temp1 + z.max, col = col.axis, lty = lty.axis)
            ## Z
            temp <- yx.f * y.max
            temp1 <- yz.f * y.max
            lines(c(temp + x.min, temp + x.min), c(z.min + temp1, z.max + temp1),
                  col = col.axis, lty = if(!angle.2) lty.hide else lty.axis)
            lines(c(x.max + temp, x.max + temp), c(z.min + temp1, z.max + temp1),
                  col = col.axis, lty = if(angle.2) lty.hide else lty.axis)
        }
    }

### plot points
    x <- dat$x + (dat$y * yx.f)
    z <- dat$z + (dat$y * yz.f)
    col <- as.character(dat$col)
    if(type == "h") {
        z2 <- dat$y * yz.f + z.min
        segments(x, z, x, z2, col = col, cex = cex.symbols, lty = lty.hplot, ...)
        points(x, z, type = "p", col = col, pch = dat$pch, bg = dat$bg, cex = cex.symbols, ...)
    }
    else points(x, z, type = type, col = col, pch = dat$pch, bg = dat$bg, cex = cex.symbols, ...)

### box-lines in front of points (overlay)
    if(axis && box) {
        lines(c(x.min, x.max), c(z.max, z.max),
              col = col.axis, lty = lty.axis)
        lines(c(0, y.max * yx.f) + x.max, c(0, y.max * yz.f) + z.max,
              col = col.axis, lty = lty.axis)
        lines(xx[c(1,1)], c(z.min, z.max), col = col.axis, lty = lty.axis)
    }


    # par(mem.par) # we MUST NOT set the margins back
### Return Function Object
    ob <- ls() ## remove all unused objects from the result's enviroment:
    rm(list = ob[!ob %in% c("angle", "mar", "usr", "x.scal", "y.scal", "z.scal", "yx.f",
        "yz.f", "y.add", "z.min", "z.max", "x.min", "x.max", "y.max", "x.range.fix", "y.range.fix",
        "xlabel", "ylabel", "zlabel", "x.prty", "y.prty", "z.prty", "mem.par")])
    rm(ob)
    invisible(list(
        xyz.convert = function(x, y=NULL, z=NULL) {
            xyz <- xyz.coords(x, y, z)
            if(angle > 2) { ## switch y and x axis to ensure righthand oriented coord.
                temp <- xyz$x; xyz$x <- xyz$y; xyz$y <- temp
            }
            y <- (xyz$y - y.add) / y.scal
            return(list(x = xyz$x / x.scal + yx.f * y,
                y = xyz$z / z.scal + yz.f * y))
        },
        points3d = function(x, y = NULL, z = NULL, type = "p", ...) {
            xyz <- xyz.coords(x, y, z)
            if(angle > 2) { ## switch y and x axis to ensure righthand oriented coord.
                temp <- xyz$x; xyz$x <- xyz$y; xyz$y <- temp
            }
            y2 <- (xyz$y - y.add) / y.scal
            x <- xyz$x / x.scal + yx.f * y2
            y <- xyz$z / z.scal + yz.f * y2
            mem.par <- par(mar = mar, usr = usr)
            #on.exit(par(mem.par))
            if(type == "h") {
                y2 <- z.min + yz.f * y2
                segments(x, y, x, y2, ...)
                points(x, y, type = "p", ...)
            }
            else points(x, y, type = type, ...)
        },
        plane3d = function(Intercept, x.coef = NULL, y.coef = NULL, 
            lty = "dashed", lty.box = NULL, draw_lines = TRUE, draw_polygon = FALSE,
            polygon_args = list(border = NA, col = rgb(0,0,0,0.2)), 
            ...){
            if(!is.atomic(Intercept) && !is.null(coef(Intercept))){
                Intercept <- coef(Intercept)
                if(!("(Intercept)" %in% names(Intercept)))
                    Intercept <- c(0, Intercept)
            }
            if(is.null(lty.box)) lty.box <- lty
            if(is.null(x.coef) && length(Intercept) == 3){
                x.coef <- Intercept[if(angle > 2) 3 else 2]
                y.coef <- Intercept[if(angle > 2) 2 else 3]
                Intercept <- Intercept[1]
            }
            mem.par <- par(mar = mar, usr = usr)
            #on.exit(par(mem.par))
            x <- x.min:x.max
            y <- 0:y.max
            
            ltya <- c(lty.box, rep(lty, length(x)-2), lty.box)
            x.coef <- x.coef * x.scal
            z1 <- (Intercept + x * x.coef + y.add * y.coef) / z.scal
            z2 <- (Intercept + x * x.coef +
                (y.max * y.scal + y.add) * y.coef) / z.scal

            if(draw_polygon) 
                do.call("polygon", c(list(
                    c(x.min, x.min + y.max * yx.f, x.max + y.max * yx.f, x.max),
                    c(z1[1], z2[1] + yz.f * y.max, z2[length(z2)] + yz.f * y.max, z1[length(z1)])), 
                  polygon_args))
            if(draw_lines) 
                segments(x, z1, x + y.max * yx.f, z2 + yz.f * y.max, lty = ltya, ...)

            ltya <- c(lty.box, rep(lty, length(y)-2), lty.box)
            y.coef <- (y * y.scal + y.add) * y.coef
            z1 <- (Intercept + x.min * x.coef + y.coef) / z.scal
            z2 <- (Intercept + x.max * x.coef + y.coef) / z.scal
            if(draw_lines) 
                segments(x.min + y * yx.f, z1 + y * yz.f,
                  x.max + y * yx.f, z2 + y * yz.f, lty = ltya, ...)
        },
        box3d = function(...){
            mem.par <- par(mar = mar, usr = usr)
            #on.exit(par(mem.par))
            lines(c(x.min, x.max), c(z.max, z.max), ...)
            lines(c(0, y.max * yx.f) + x.max, c(0, y.max * yz.f) + z.max, ...)
            lines(c(0, y.max * yx.f) + x.min, c(0, y.max * yz.f) + z.max, ...)
            lines(c(x.max, x.max), c(z.min, z.max), ...)
            lines(c(x.min, x.min), c(z.min, z.max), ...)
            lines(c(x.min, x.max), c(z.min, z.min), ...)
        },
        contour3d = function(f, x.count = 10, y.count = 10, type = "l", lty = "24", 
            x.resolution = 50, y.resolution = 50, ...) {    
            if(class(f) == "lm"){
                #orig.vars <- c(xlabel, ylabel, zlabel)
                #orig.vars <- gsub(".*\\$", "", orig.vars)
                vars <- all.vars(formula(f))
            } else vars <- c("z", "x", "y")

            #vars.ordering <- names(sort(sapply(vars, function(v) grep(v, orig.vars)))[1:2])
            
            # x vor y in Formel!
            for(x1 in seq(x.range.fix[1], x.range.fix[2], length = x.count)){
                d <- data.frame(x1, seq(y.range.fix[1], y.range.fix[2], length = y.resolution))
                names(d) <- vars[-1]
                if(class(f) == "lm"){
                    d[vars[1]] <- predict(f, newdata=d)      
                } else d[vars[1]] <- f(d[[1]], d[[2]])
                xyz <- xyz.coords(d)
                if(angle > 2) { ## switch y and x axis to ensure righthand oriented coord.
                    temp <- xyz$x; xyz$x <- xyz$y; xyz$y <- temp
                }
                y2 <- (xyz$y - y.add) / y.scal
                x <- xyz$x / x.scal + yx.f * y2
                y <- xyz$z / z.scal + yz.f * y2
                mem.par <- par(mar = mar, usr = usr)
                if(type == "h") {
                    y2 <- z.min + yz.f * y2
                    segments(x, y, x, y2, ...)
                    points(x, y, type = "p", ...)
                }
                else points(x, y, type = type, lty = lty, ...)
            }       
            for(x2 in seq(y.range.fix[1], y.range.fix[2], length = y.count)){
                d <- data.frame(seq(x.range.fix[1], x.range.fix[2], length = x.resolution), x2)
                names(d) <- vars[-1]
                if(class(f) == "lm"){
                    d[vars[1]] <- predict(f, newdata=d)      
                } else d[vars[1]] <- f(d[[1]], d[[2]])
                xyz <- xyz.coords(d)
                if(angle > 2) { ## switch y and x axis to ensure righthand oriented coord.
                    temp <- xyz$x; xyz$x <- xyz$y; xyz$y <- temp
                }
                y2 <- (xyz$y - y.add) / y.scal
                x <- xyz$x / x.scal + yx.f * y2
                y <- xyz$z / z.scal + yz.f * y2
                mem.par <- par(mar = mar, usr = usr)
                if(type == "h") {
                    y2 <- z.min + yz.f * y2
                    segments(x, y, x, y2, ...)
                    points(x, y, type = "p", ...)
                }
                else points(x, y, type = type, lty = lty, ...)
            }
        },
        par.mar = mem.par
    ))
}


plotPCOA3d         <- function(axes,loading,design=FALSE,ax1=1, ax2=2, ax3=3, Lab1="PC1",Lab2="PC2",
                             Lab3="PC3", Pch = 19, Main = "PCoA plot in 3D", Sub=NULL,
                             Scale.y=1, Angle=40,
                             Axis=TRUE, Tick.marks=TRUE, Label.tick.marks=TRUE,
                             X.ticklabs=NULL, Y.ticklabs=NULL, Z.ticklabs=NULL,
                             Y.margin.add=0, Grid=TRUE, Box=TRUE, Lab=par("lab"),
                             Lab.z=mean(Lab[1:2]), Type="p", LX="topleft", LY= NULL){
  #require(scatterplot3d)
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  if (substrRight(axes,5) != ".axes") {
    print("warning file not .axes")
  }
  data           <- read.table(axes,header = T)
  names          <- data[,1]
  axis           <- data[,-1]
  if (ncol(axis) < 3) { return("No three axes!") }
  
  axi            <- axis[,c(ax1,ax2,ax3)]
  colnames(axi) <- c("axis1","axis2","axis3")
  
  loadings <- read.table(loading,header=T)
  Lab1 = paste(Lab1," (",loadings[1,2],"%)",sep="")
  Lab2 = paste(Lab2," (",loadings[2,2],"%)",sep="")
  Lab3 = paste(Lab3," (",loadings[3,2],"%)",sep="")
  
  png(file = "PCoA3d_plot.png",res=600,width=4800, height=4800)
  col=rainbow(length(names))
  col=sample(col)
  
  if (design==FALSE){
    scatterplot3d(x = axi$"axis1", y = axi$"axis2", z = axi$"axis3",xlab = Lab1, ylab = Lab2, zlab = Lab3,
                  color = rev(col),pch=Pch, main=Main, sub=Sub, scale.y=Scale.y,
                  angle=Angle, axis=Axis, tick.marks=Tick.marks,label.tick.marks=Label.tick.marks,
                  x.ticklabs=X.ticklabs, y.ticklabs=Y.ticklabs, z.ticklabs=Z.ticklabs, y.margin.add=Y.margin.add,
                  grid=Grid, box=Box, lab=Lab, lab.z=Lab.z, type=Type)
	legend(x=LX,y=LY, legend = c(levels(names)), col = col, pch = 16, cex=0.7, pt.cex = 1.2)
  }
  else{
    if (substrRight(design,7) != ".design") {
      warning("Warning message: design file is not .design")
    }
    designtable    <- suppressWarnings(read.table(design, header=F))
    tDesign        <- t(designtable)
    taxi           <- cbind(tDesign[,2],axi)
    names(taxi)[1] <- "partition" 
    
    scatterplot3d(x = taxi$"axis1", y = taxi$"axis2", z = taxi$"axis3", xlab = Lab1, ylab = Lab2, zlab = Lab3,
                  color = as.integer(factor(taxi$"partition")), pch = Pch, main=Main, sub=Sub, scale.y=Scale.y,
                  angle=Angle, axis=Axis, tick.marks=Tick.marks,label.tick.marks=Label.tick.marks,
                  x.ticklabs=X.ticklabs, y.ticklabs=Y.ticklabs, z.ticklabs=Z.ticklabs, y.margin.add=Y.margin.add,
                  grid=Grid, box=Box, lab=Lab, lab.z=Lab.z, type=Type)
    legend(x=LX,y=LY, legend = c(levels(taxi$partition)), col = as.integer(factor(levels(taxi$"partition"))), pch = Pch)
  }
  
  dev.off()
}

args<-commandArgs(T)
plotPCOA3d(axes=args[1],loading=args[2])
