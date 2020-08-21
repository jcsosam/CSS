adjacency.plot <- function (mat, show.grid = FALSE, cex.axis, tick, labs, col.axis, ...)
{ 
     JJ <- dim(mat)[1]
     colorscale <- c("white", rev(heat.colors(100)))
     if(missing(labs))     labs <- 1:JJ
     if(missing(col.axis)) col.axis <- rep("black", JJ)
     if(missing(cex.axis)) cex.axis <- 0.5
     if(missing(tick))     tick <- TRUE
     ## adjacency matrix
     image(seq(1, JJ), seq(1, JJ), mat, axes = FALSE, xlab = "", ylab = "",
           col = colorscale[seq(floor(100*min(mat)), floor(100*max(mat)))], ...)
     for(j in 1:JJ){
          axis(1, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick, col.axis = col.axis[j])
          axis(2, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick, col.axis = col.axis[j])
     }
     box()
     if(show.grid) grid(nx=JJ, ny=JJ)
}

heat.plot0 <- function (mat, show.grid = FALSE, cex.axis, tick, labs, col.axis,...)
{ 
     JJ <- dim(mat)[1]
     colorscale <- c("white", rev(heat.colors(100)))
     if(missing(labs))     labs <- 1:JJ
     if(missing(col.axis)) col.axis <- rep("black", JJ)
     if(missing(cex.axis)) cex.axis <- 0.5
     if(missing(tick))     tick <- TRUE
     ## adjacency matrix
     image(seq(1, JJ), seq(1, JJ), mat, axes = FALSE, xlab = "", ylab = "",
           col = colorscale[seq(floor(100*min(mat)), floor(100*max(mat)))], ...)
     for(j in 1:JJ){
          axis(1, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick, col.axis = col.axis[j])
          axis(2, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick, col.axis = col.axis[j])
     }
     box()
     if(show.grid) grid(nx = JJ, ny = JJ)
}

heat.plot <- function (mat, show.grid = FALSE, cex.axis, tick, labs, col.axis, ...)
{ 
     JJ <- dim(mat)[1]
     colorscale <- c("white", rev(heat.colors(100)))
     if(missing(labs))     labs <- 1:JJ
     if(missing(col.axis)) col.axis <- rep("black", JJ)
     if(missing(cex.axis)) cex.axis <- 0.5
     if(missing(tick))     tick <- TRUE
     ## adjacency matrix
     par(layout(mat = matrix(1:2, 1, 2), widths = c(9, 1)), mar = c(4, 2.5, 2, 1))
     image(seq(1, JJ), seq(1, JJ), mat, axes = FALSE, xlab = "", ylab = "",
           col = colorscale[seq(floor(100*min(mat)), floor(100*max(mat)))], ...)
     for(j in 1:JJ){
          axis(1, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick, col.axis = col.axis[j])
          axis(2, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick, col.axis = col.axis[j])
     }
     box()
     if(show.grid) grid(nx = JJ, ny = JJ)
     ## scale
     par(mar = c(3.1, 0.1, 1.3, 0.1))
     plot(1:100, 1:100, xlim = c(0,2), ylim = c(0, 100), axes = FALSE, type = "n", 
          xlab = "", ylab = "")
     yposr <- 1:100
     rect(0, yposr-.5, 0.5, yposr+.5, col = rev(heat.colors(100)), border = FALSE)
     rect(0, 0.5, 0.5, 100.5, col = "transparent")
     text(0.42, c(yposr[1], yposr[25], yposr[50], yposr[75], yposr[100]), 
          c("0", "0.25", "0.5", "0.75", "1"), pos = 4, cex = 0.6)
}

chain.plot <- function (xdom, x, ... )
{
     plot(xdom, x, cex = 0.25, type = 'p', ylim = range(x), 
          xlab = 'Iteration', ... )
}

chain.plot1 <- function (x, plot.mean = TRUE, xtruth = NULL, ... )
{
     N <- length(x)
     plot(1:N, x, xlab = 'Iteration', ylab = '', col = 'gray', ... )
     if (plot.mean) abline(h = mean(x), col = 'black', lwd = 2)
     if (!is.null(xtruth)) abline(h = xtruth, col = 'red', lty = 2)
}

post.plot <- function (x, adj = 1.15, n = 1000, col.border = 'gray', 
                       info = FALSE, ylim = NULL, col.shade = 'gray95', 
                       lwd.density = 1, tck = 0.015, summaries = TRUE, 
                       col.lines = 4, lty.lines = 3, lwd.lines = 1, main = NULL,
                       xtruth = NULL)
{
     den  <- density(x, adj = adj, n = n) 
     xx   <- den$x
     yy   <- den$y
     if(is.null(ylim)){ 
          ylim <- 1.05*c(0, max(range(yy)[2], hist(x, plot = FALSE)$density)) 
     } else { 
          ylim <- ylim 
     }
     
     info <- ifelse(info, paste('n =',length(x),'   ','bandwith =', round(den$bw, 2)),'')
     if(is.null(main)) main=''
     
     hist(x, freq = FALSE, xaxt = 'n', ann = FALSE, border = col.border, 
          ylim = ylim)
     polygon(c(xx,rev(xx)), c(rep(0,length(xx)), rev(yy)), lty = 1, col = col.shade, 
             xaxt = 'n', ann = FALSE)
     hist(x, freq = FALSE, xaxt = 'n', ann = FALSE, border = col.border, add = TRUE)
     if (!is.null(xtruth)) abline(v = xtruth, col = 'red', lty = 2)
     title(main)
     lines(den, lwd = lwd.density)
     
     axis(side = 1, at = x, labels = FALSE, tck = tck)
     m <- abs(max(x))
     r <- 1
     while (m < 1) {
          m <- m*10^r
          r <- r + 1
     }
     xticks <- round(seq(min(x), max(x), len=5), r)
     axis(side = 1, at = xticks, labels = xticks)
     
     if(summaries){
          closest <- function(x, x0) which(abs(x - x0) == min(abs(x - x0)))
          idx1    <- closest(xx, quantile(x, 0.025))
          idx2    <- closest(xx, quantile(x, 0.975))
          idx3    <- closest(xx, mean(x))
          idx     <- c(idx1,idx2,idx3)
          segments(x0 = xx[idx], y0 = 0, x1 = xx[idx], y1 = yy[idx], 
                   col = col.lines, lty = lty.lines, lwd = lwd.lines)
          axis(side = 1, at = xx[idx], col.ticks = col.lines, lwd.ticks = 1, 
               labels = FALSE)
     }
}

post.both <- function(x, xtruth = NULL, main.title = " ") 
{
     windows(width = 10, height = 5)
     par(mar = c(4, 4, 3, 2) - 1, mgp = c(2, 1, 0), oma = c(1, 1, 1, 1))
     layout(matrix(c(1, 2), 1, 2), heights = rep(1, 2), widths = c(1.75, 1))
     chain.plot1(x, xtruth = xtruth, type = 'p', cex = 0.5)
     post.plot(x, xtruth = xtruth)
     title(main.title, outer = T)
}

post.cis0 <- function(x, xtruth = NULL, plot.xaxis = TRUE, yrange = NULL, ...)
{
     n <- dim(x)[2]
     li <- apply(x, 2, quantile, probs = c(0.025))
     ul <- apply(x, 2, quantile, probs = c(0.975))
     pm <- colMeans(x)
     
     if (!is.null(xtruth)) {
          idx <- ((li <= xtruth) & (xtruth <= ul)) + 1
          ok  <- round(sum(idx - 1)/n, 3) * 100
          notok <- round(100 - ok, 1)
          if(is.null(yrange)) yrange <- 1.1 * range(li, ul, xtruth)
     } else {
          if(is.null(yrange)) yrange <- 1.1 * range(li, ul)
     }
     
     ntop <- 30
     if (n > ntop) {
          id.lab <- sort(sample(1:n, ntop))
          li <- li[id.lab]
          ul <- ul[id.lab]
          pm <- pm[id.lab]
          if (!is.null(xtruth)) { 
               xtruth <- xtruth[id.lab] 
               idx <- idx[id.lab]
               if(is.null(yrange)) yrange = 1.1 * range(li, ul, xtruth)
          } else {
               if(is.null(yrange)) yrange = 1.1 * range(li, ul)
          }
          n <- ntop
     } else {
          id.lab <- 1:n
     }
     
     par(mfrow = c(1, 1), mar = c(4, 4, 3, 2) - 0.5, mgp = c(2, 1, 0), oma = c(1, 1, 1, 1))
     plot(1:n, ylim = yrange, type = "n", xaxt = 'n', ...)
     abline(v = 1:n, lty = 3, col = 'lightgray')
     if(plot.xaxis) axis(side = 1, at = 1:n, labels = id.lab, cex.axis = 0.6)
     points(1:n, li, pch = "-", col = "darkorange", cex = 1.2)
     points(1:n, ul, pch = "-", col = "darkorange", cex = 1.2)
     segments(1:n, li, 1:n, ul, lwd = 2, col = "darkorange")
     points(1:n, pm, pch = 18, col = "dodgerblue2", cex = 1)
     if (!is.null(xtruth)) {
          col.points <- c("red", "green")[idx]
          pch.points <- c(4, 18)[idx] 
          lines(xtruth, type = 'p', col = col.points, pch = pch.points)
          legend("topleft", legend = c(paste(ok, "%", sep =''), paste(notok, "%", sep='')), 
                 col = c(3, 2), pch = c(18, 4), horiz = T)
     }
}

post.cis <- function(x, idlabs, yrange, hlines = T, h = 0, cex.idlabs = 0.9, 
                     cex.pm = 0.9, pch.pm = 18, h.lwd, ...)
{
     n  <- dim(x)[2]
     li <- apply(x, 2, quantile, probs = c(0.025))
     ul <- apply(x, 2, quantile, probs = c(0.975))
     pm <- colMeans(x)
     if (missing(yrange)) yrange <- 1.2 * range(li, ul)
     par(mfrow = c(1, 1), mar = c(4, 4, 1, 2) - 1, mgp = c(2, 1, 0), oma = c(1, 1, 1, 1))
     plot(1:n, ylim = yrange, type = "n", xaxt = "n", ...)
     abline(h = h, lty = 2, col = "lightgray")
     if (hlines) abline(v = 1:n, lty = 1, col = "lightgray")
     if (missing(idlabs)) idlabs <- 1:n
     #axis(side = 1, at = 1:n, labels = idlabs, cex.axis = cex.idlabs)
     axis(side = 1, at = 1:n, labels = NA)
     axis(side = 1, at = 1:n, cex.axis = cex.idlabs, labels = idlabs, line = -0.5, lwd = 0)
     points(1:n, li, pch = " ", col = "black", cex = 1.2)
     points(1:n, ul, pch = " ", col = "black", cex = 1.2)
     if (missing(h.lwd)) h.lwd <- c(2, 2)
     for (i in 1:n) {
          if ((li[i] < h) & (h < ul[i])) { 
               lwd <- h.lwd[1]
               cex <- cex.pm
          } else {
               lwd <- h.lwd[2]
               cex <- 1.1*cex.pm
          }
          segments(i, li[i], i, ul[i], lwd = lwd, col = "black")
          points(i, pm[i], pch = pch.pm, col = "black", cex = cex)
     }
}

post.beta <- function(x, xlabs, cex.axis = 1, ...)
{
        J <- nrow(x)
        xrange <- range(1:J) 
        yrange <- max(abs(range(x))) * c(-1, 1)
        
        plot(NA, NA, type = "n", xaxt = "n", xlim = xrange, ylim = yrange,
             xlab = 'Actor', cex.axis = cex.axis, ...)
        axis(side = 1, at = 0:J, labels = xlabs, cex.axis = cex.axis)
        abline(v = 1:J, col = 'gray', lty = 1)
        abline(h = 0,   col = 'gray', lty = 2)
        
        pch <- 18
        for (j in 1:J) {
                points(x = j, y = x[j, 1], cex = 1, pch = pch,  col = 'black')  # mean
                if ((x[j, 2] < 0) & (0 < x[j, 3])) { 
                        lwd <- 1
                } else {
                        lwd <- 2
                }
                segments(x0 = j, y0 = x[j, 2], x1 = j, y1 = x[j, 3], lwd = lwd, col = 'black')
        }
}
