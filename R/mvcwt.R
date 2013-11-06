#
# All code copyright 2013 Timothy H. Keitt
#

get.nscales = function(x) length(as.vector(unlist(x)))
get.min.scale = function(x) 2 * median(diff(as.vector(unlist(x))))
get.max.scale = function(x) diff(range(as.vector(unlist(x)))) / 2
  
mvcwt = function(x, y,
                 scale.trim = 1,
                 scale.exp = 0.5,
                 nscales = get.nscales(x),
                 min.scale = scale.trim * get.min.scale(x),
                 max.scale = get.max.scale(x) / scale.trim,
                 scales = log2Bins(min.scale, max.scale, nscales),
                 loc = regularize(x), wave.fun = "Morlet")
{
  s = 1 # this is a workaround for a bug in R's code checking
  wave.fun = match.fun(wave.fun)
  x = as.vector(unlist(x))
  lmat = lagMat(x, loc)
  y = matrix(unlist(y), nrow = length(x))
  w = foreach(s = scales) %dopar%
  {
    Conj(wave.fun(lmat / s) / s ^ scale.exp) %*% y
  }
  w = array(unlist(w), dim = c(length(loc), ncol(y), length(scales)))
  w = aperm(w, perm = c(1, 3, 2))
  structure(list(x = loc, y = scales, z = w), class = "mvcwt")
}

print.mvcwt = function(x, ...)
{
  print(str(x), ...)
  invisible(x)
}

wmr = function(w, smoothing = 1)
{
  with(w, {
    mods = Mod(rowSums(z, dims = 2))
    smod = rowSums(Mod(z), dims = 2)
    lmat = lagMat(x)
    exports = c("Gauss", "smoothing")
    flibs = c("mvcwt")
    modrat = foreach(i = 1:length(y),
                     .combine = cbind,
                     .export = exports,
                     .packages = flibs) %dopar%
    {
      kern = Gauss(lmat / y[i] / smoothing)
      modsv = kern %*% mods[,i]
      smodv = kern %*% smod[,i]
      modsv / smodv
    }
    dim(modrat) = c(length(x), length(y), 1)
    structure(list(x = x, y = y, z = modrat), class = "mvcwt")
  })
}

wmr.boot = function(w, smoothing = 1, reps = 1000, mr.func = "wmr")
{
  mr.func = match.fun(mr.func)
  mr.obs = mr.func(w, smoothing = smoothing)
  with(w, {
    nloc = length(x)
    nvars = dim(z)[3]
    nscales = length(y)
    exports = c("reps", "wmr", "lagMat", "regularize", "mr.func",
                "Gauss", "mr.obs", "nscales", "smoothing")
    flibs = c("mvcwt")
    mr.obs$z.boot = foreach(i = 1:nscales,
                            .combine = c,
                            .export = exports,
                            .packages = flibs) %dopar%
    {
      mr.boot = foreach(j = 1:reps,
                        .combine = cbind,
                        .inorder = FALSE) %dopar%
      {
        rphase = t(array(runif(nvars, -pi, pi), dim = c(nvars, nloc)))
        zp = z[, i,, drop = FALSE] * complex(argument = rphase)
        as.vector(mr.func(list(x = x, y = y[i], z = zp), smoothing = smoothing)$z)
      }
      res = foreach(j = 1:nloc, .combine = c) %dopar%
      {
         ecdf(mr.boot[j, ])(mr.obs$z[j, i,])
      }
      res
    }
    dim(mr.obs$z.boot) = c(length(x), length(y), 1)
    return(mr.obs)
  })
}

plot.mvcwt = function(x, var = 1, scale = 1, titles = TRUE, z.fun = "Re", ...)
{
  z.fun = match.fun(z.fun)
  opar = par(no.readonly = TRUE)
  on.exit(par(opar))
  with(x, {
    nvar = length(var)
    nscal = length(scale)
    par(mfrow = c(nvar * nscal, 1), mar = rep(0.6, 4), oma = rep(5, 4), xpd = NA)
    for ( s in scale )
    {
      scale.i = which.min(abs(y - s))
      for ( i in 1:nvar )
      {
        z.out = z.fun(z[, scale.i, var[i]])
        plot(x, z.out, xlab = NA, ylab = NA, axes = FALSE, type = "l", lwd = 2, ...)
        lines(range(x), rep(median(z.out), 2), lty = 2)
        if (titles )
        {
          tstr = paste("Scale =", signif(y[scale.i], 2), "Var =", var[i], sep = " ")
          title(main = tstr, line = 0.25)
        }
        if ( i %% 2 ) axis(2) else axis(4)
      }
    }
    axis(1)
    mtext("Location", 1, 3, outer = TRUE)
    mtext("Value", 2, 3, outer = TRUE)
  })  
  invisible(x)
}

image.mvcwt = function(x, z.fun = "Re", bound = 1, reset.par = TRUE, ...)
{
  z.fun = match.fun(z.fun)
  opar = par(no.readonly = TRUE)
  if ( reset.par ) on.exit(par(opar))
  pal = colorRampPalette(rev(brewer.pal(11, 'Spectral')))(1024)
  with(x, {
    nvar = ifelse(length(dim(z)) == 3, dim(z)[3], 1)
    par(mfrow = c(nvar, 1), mar = rep(0.2, 4), oma = rep(5, 4))
    for ( i in 1:nvar )
    {
      image(x, y, z.fun(z[,,i]), log = "y", col = pal, axes = FALSE, ...)
      if ( i %% 2 ) axis(2) else axis(4)
      if ( exists("z.boot") && !is.null(z.boot) )
      {
        z.boot = 1 - abs(1 - 2 * z.boot)
        contour(x, y, z.boot[,,i], levels = 0.05, lty = 3, add = TRUE, drawlabels = FALSE)
        zb = p.adjust(as.vector(z.boot), method = "BY")
        dim(zb) = dim(z.boot)
        contour(x, y, zb[,,i], levels = 0.05, lwd = 2, add = TRUE, drawlabels = FALSE)
      }
      if ( is.finite(bound) )
      {
        lines(min(x) + bound * y, y, lty = 2, lwd = 2, col = "darkgrey")
        lines(max(x) - bound * y, y, lty = 2, lwd = 2, col = "darkgrey")
      }
      box()
    }
    axis(1)
    mtext("Location", 1, 3, outer = TRUE)
    mtext("Scale", 2, 3, outer = TRUE)
  })
  return(invisible(x))     
}

contour.mvcwt = function(x, z.fun = "Re", bound = 1, reset.par = TRUE, ...)
{
  z.fun = match.fun(z.fun)
  opar = par(no.readonly = TRUE)
  if ( reset.par ) on.exit(par(opar))
  with(x, {
    nvar = ifelse(length(dim(z)) == 3, dim(z)[3], 1)
    par(mfrow = c(nvar, 1), mar = rep(0.2, 4), oma = rep(5, 4))
    for ( i in 1:nvar )
    {
      contour(x, y, z.fun(z[,,i]), log = "y", axes = FALSE, ...)
      if ( i %% 2 ) axis(2) else axis(4)
      if ( is.finite(bound) )
      {
        lines(min(x) + bound * y, y, lty = 2, lwd = 2, col = "darkgrey")
        lines(max(x) - bound * y, y, lty = 2, lwd = 2, col = "darkgrey")
      }
      box()
    }
    axis(1)
    mtext("Location", 1, 3, outer = TRUE)
    mtext("Scale", 2, 3, outer = TRUE)
  })
  return(invisible(x))     
}

midp <- function(x)
{
  (x[-1] + x[-length(x)]) / 2
}

log2Bins <- function(min, max, nbins)
{
  2 ^ midp(seq(log2(min), log2(max), length = nbins + 1))
}

regularize <- function(x, nsteps = length(as.vector(unlist(x))))
{
  seq(min(x), max(x), length = nsteps)
}
    
lagMat <- function(from, to = regularize(from))
{
  outer(to, from, "-")
}

Gauss <- function(lag)
{
  exp(-lag ^ 2 / 2) / sqrt(2 * pi)
}
  
Morlet <- function(lag)
{
  exp(-lag ^ 2 / 2 + 2i * pi * lag) / pi ^ 4
}

