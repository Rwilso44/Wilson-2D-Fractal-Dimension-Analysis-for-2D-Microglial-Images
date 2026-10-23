# Box-counting fractal dimension, lacunarity, and convex-hull morphometrics
# for a binary or grayscale image (base R + the 'tiff' package)
#
# Install once:   install.packages("tiff")
# Use:            source("boxcount.R")        # or open in RStudio and click Source
#                 res <- box_count_dimension("my_image.tif", threshold = 0.5)
#                 res$D; res$hull
#                 # white features instead of black:  invert = TRUE
#                 show_boxes("my_image.tif", size = 16)              # outline the 16 px boxes
#                 show_boxes("my_image.tif", size = 16, save = "boxes16.png")
#                 show_hull("my_image.tif")                          # draw the convex hull
#                 box_count_dimension("my_image.tif", csv = "results.csv")
#
# Measures reported
#   Fractal dimension D   slope of log(boxes) vs log(1/box size)
#   Lacunarity            (sd / mean)^2 of foreground pixels per box, per box size,
#                         plus the mean over all box sizes
#   Cell area             number of foreground pixels
#   Convex hull           area and perimeter (pixel units; each pixel is a unit square)
#   Density               cell area / convex hull area
#   Circularity           4*pi*hull area / hull perimeter^2   (1 = perfect circle)
#   Max span              longest distance between two points on the hull
#   Span ratio            major / minor axis of the ellipse with the same second
#                         moments as the hull polygon

library(tiff)

load_binary <- function(path, threshold = 0.5, invert = FALSE) {
  img <- readTIFF(path)                       # values scaled 0-1
  if (length(dim(img)) == 3) {                # color -> grayscale
    img <- 0.299 * img[, , 1] + 0.587 * img[, , 2] + 0.114 * img[, , 3]
  }
  features <- img < threshold                 # dark pixels = features (like Fractalyse)
  if (invert) features <- !features
  features
}

# ------------------------------------------------------------ box counting
box_count <- function(binary, sizes) {
  h <- nrow(binary); w <- ncol(binary)
  sapply(sizes, function(s) {
    H <- ceiling(h / s) * s; W <- ceiling(w / s) * s
    padded <- matrix(FALSE, H, W)
    padded[1:h, 1:w] <- binary
    bi <- (row(padded) - 1) %/% s
    bj <- (col(padded) - 1) %/% s
    length(unique(bi[padded] * (W / s) + bj[padded]))
  })
}

# Lacunarity per box size: (sd / mean)^2 of foreground pixels per box.
# Uses only complete boxes and a single fixed grid from the top-left corner.
lacunarity_by_size <- function(binary, sizes) {
  h <- nrow(binary); w <- ncol(binary)
  sapply(sizes, function(s) {
    H <- (h %/% s) * s; W <- (w %/% s) * s
    if (H == 0 || W == 0) return(NA_real_)
    sub <- binary[1:H, 1:W, drop = FALSE]
    nbc <- W / s
    ids <- ((row(sub) - 1) %/% s)[sub] * nbc + ((col(sub) - 1) %/% s)[sub] + 1
    mass <- tabulate(ids, nbins = (H / s) * nbc)   # empty boxes count as 0
    m <- mean(mass)
    if (m > 0) mean((mass - m)^2) / m^2 else NA_real_
  })
}

# ------------------------------------------------------------ convex hull
# Hull vertices (counter-clockwise), treating every pixel as a unit square.
# Only the leftmost and rightmost pixel of each row can be hull vertices.
hull_vertices <- function(binary) {
  rows <- which(rowSums(binary) > 0)
  pts <- do.call(rbind, lapply(rows, function(r) {
    cols <- which(binary[r, ])
    x0 <- min(cols) - 1; x1 <- max(cols)
    rbind(c(x0, r - 1), c(x0, r), c(x1, r - 1), c(x1, r))
  }))
  pts <- unique(pts)
  v <- pts[rev(chull(pts)), , drop = FALSE]    # chull is clockwise; reverse it
  v
}

hull_metrics <- function(binary) {
  cell_area <- sum(binary)
  out <- list(cell_area = cell_area, hull_area = NA_real_, hull_perimeter = NA_real_,
              density = NA_real_, circularity = NA_real_, max_span = NA_real_,
              span_ratio = NA_real_)
  v <- hull_vertices(binary)
  if (nrow(v) < 3) return(out)

  x <- v[, 1]; y <- v[, 2]
  nx <- c(x[-1], x[1]); ny <- c(y[-1], y[1])
  a <- x * ny - nx * y
  if (sum(a) < 0) {                            # make sure orientation is counter-clockwise
    v <- v[nrow(v):1, , drop = FALSE]
    x <- v[, 1]; y <- v[, 2]; nx <- c(x[-1], x[1]); ny <- c(y[-1], y[1])
    a <- x * ny - nx * y
  }
  area  <- sum(a) / 2
  perim <- sum(sqrt((nx - x)^2 + (ny - y)^2))
  max_span <- max(dist(v))

  # second moments of the hull polygon -> equivalent ellipse
  cx <- sum((x + nx) * a) / (6 * area)
  cy <- sum((y + ny) * a) / (6 * area)
  ixx <- sum((y^2 + y * ny + ny^2) * a) / 12 - area * cy^2
  iyy <- sum((x^2 + x * nx + nx^2) * a) / 12 - area * cx^2
  ixy <- sum((x * ny + 2 * x * y + 2 * nx * ny + nx * y) * a) / 24 - area * cx * cy
  cv <- matrix(c(iyy, ixy, ixy, ixx), 2) / area
  ev <- sort(eigen(cv, symmetric = TRUE)$values, decreasing = TRUE)
  span_ratio <- if (ev[2] > 1e-12) sqrt(ev[1] / ev[2]) else Inf

  list(cell_area = cell_area, hull_area = area, hull_perimeter = perim,
       density = cell_area / area, circularity = 4 * pi * area / perim^2,
       max_span = max_span, span_ratio = span_ratio)
}

# ---------------------------------------------------------------- drawing
# Draw the image with every occupied s x s box outlined in red
show_boxes <- function(path, size, threshold = 0.5, invert = FALSE, save = NULL) {
  binary <- load_binary(path, threshold, invert)
  h <- nrow(binary); w <- ncol(binary)
  if (!is.null(save)) png(save, width = 1200, height = round(1200 * h / w))
  op <- par(mar = c(1, 1, 3, 1)); on.exit(par(op))

  z <- t(!binary[h:1, , drop = FALSE]) + 0     # features black, background white
  image(1:w, 1:h, z, zlim = c(0, 1), col = c("black", "white"),
        axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)

  n <- 0
  for (i in seq(1, h, by = size)) {
    for (j in seq(1, w, by = size)) {
      block <- binary[i:min(i + size - 1, h), j:min(j + size - 1, w), drop = FALSE]
      if (any(block)) {
        n <- n + 1
        rect(j - 0.5, h - min(i + size - 1, h) + 0.5,
             min(j + size - 1, w) + 0.5, h - i + 1.5,
             border = "red", lwd = 1.5)
      }
    }
  }
  title(sprintf("Box size %d px: %d occupied boxes", size, n))
  if (!is.null(save)) { dev.off(); cat("Saved box picture to", save, "\n") }
  invisible(n)
}

# Draw the image with its convex hull outlined in blue
show_hull <- function(path, threshold = 0.5, invert = FALSE, save = NULL) {
  binary <- load_binary(path, threshold, invert)
  h <- nrow(binary); w <- ncol(binary)
  if (!is.null(save)) png(save, width = 1200, height = round(1200 * h / w))
  op <- par(mar = c(1, 1, 3, 1)); on.exit(par(op))

  z <- t(!binary[h:1, , drop = FALSE]) + 0
  image(1:w, 1:h, z, zlim = c(0, 1), col = c("black", "white"),
        axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
  v <- hull_vertices(binary)
  if (nrow(v) >= 3) {
    px <- c(v[, 1], v[1, 1]) + 0.5           # corner coordinates -> plot coordinates
    py <- h - c(v[, 2], v[1, 2]) + 0.5
    lines(px, py, col = "blue", lwd = 2)
  }
  title("Convex hull")
  if (!is.null(save)) { dev.off(); cat("Saved hull picture to", save, "\n") }
  invisible(v)
}

# ------------------------------------------------------------------- main
box_count_dimension <- function(path, threshold = 0.5, invert = FALSE,
                                min_size = 1, max_size = NULL, plot = TRUE,
                                csv = NULL) {
  binary <- load_binary(path, threshold, invert)
  if (!any(binary)) stop("No feature pixels found - try invert = TRUE or a different threshold.")
  if (is.null(max_size)) max_size <- min(dim(binary)) %/% 2

  sizes  <- 2^(floor(log2(min_size)):floor(log2(max_size)))   # 1, 2, 4, 8, ...
  counts <- box_count(binary, sizes)
  lac    <- lacunarity_by_size(binary, sizes)
  hull   <- hull_metrics(binary)

  x <- log(1 / sizes); y <- log(counts)
  fit <- lm(y ~ x)
  D  <- unname(coef(fit)[2])
  r2 <- summary(fit)$r.squared
  mean_lac <- mean(lac, na.rm = TRUE)

  cat(sprintf("Image: %d x %d px, %.1f%% feature pixels\n",
              ncol(binary), nrow(binary), 100 * mean(binary)))
  print(data.frame(box_size = sizes, boxes = counts, lacunarity = round(lac, 4)),
        row.names = FALSE)
  cat(sprintf("\nFractal dimension D   = %.4f  (R^2 = %.4f)\n", D, r2))
  cat(sprintf("Mean lacunarity       = %.4f\n", mean_lac))
  cat(sprintf("Cell area (pixels)    = %d\n", hull$cell_area))
  cat(sprintf("Convex hull area      = %.2f\n", hull$hull_area))
  cat(sprintf("Convex hull perimeter = %.2f\n", hull$hull_perimeter))
  cat(sprintf("Density               = %.4f\n", hull$density))
  cat(sprintf("Circularity           = %.4f\n", hull$circularity))
  cat(sprintf("Max span across hull  = %.2f\n", hull$max_span))
  cat(sprintf("Span ratio            = %.4f\n", hull$span_ratio))

  if (!is.null(csv)) {
    summary_df <- data.frame(
      measure = c("fractal_dimension", "r_squared", "mean_lacunarity", names(hull)),
      value   = c(D, r2, mean_lac, unlist(hull)))
    write.csv(summary_df, csv, row.names = FALSE)
    write.csv(data.frame(box_size = sizes, boxes = counts, lacunarity = lac),
              sub("\\.csv$", "_by_box_size.csv", csv), row.names = FALSE)
    cat("\nSaved results to", csv, "and", sub("\\.csv$", "_by_box_size.csv", csv), "\n")
  }

  if (plot) {
    plot(x, y, pch = 19, xlab = "log(1 / box size)", ylab = "log(number of boxes)")
    abline(fit, col = "red")
    legend("topleft", legend = sprintf("D = %.3f", D), bty = "n")
  }
  invisible(list(D = D, r2 = r2, sizes = sizes, counts = counts,
                 lacunarity = lac, mean_lacunarity = mean_lac, hull = hull))
}

# ------------------------------------------------------- one-call wrapper
# Give analyze_image() the image's file name (as text or as a variable) and it saves
# everything under that same name. For my_file <- "DFP Proj Veh Amygdala Microglia 1.tif":
#   DFP Proj Veh Amygdala Microglia 1.csv
#   DFP Proj Veh Amygdala Microglia 1_by_box_size.csv
#   DFP Proj Veh Amygdala Microglia 1_boxes16.png
#   DFP Proj Veh Amygdala Microglia 1_hull.png
# Options: box_size (for the box picture), threshold, invert, out_dir (folder for outputs),
#          and min_size / max_size (passed on to box_count_dimension).
analyze_image <- function(name, box_size = 5, threshold = 0.5, invert = FALSE,
                          out_dir = ".", ...) {
  path <- if (grepl("\\.tiff?$", name, ignore.case = TRUE)) name else paste0(name, ".tif")
  if (!file.exists(path)) stop("Can't find ", path, " in ", getwd())
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  base <- tools::file_path_sans_ext(basename(path))
  out  <- function(suffix) file.path(out_dir, paste0(base, suffix))

  res <- box_count_dimension(path, threshold = threshold, invert = invert,
                             csv = out(".csv"), plot = FALSE, ...)
  show_boxes(path, size = box_size, threshold = threshold, invert = invert,
             save = out(paste0("_boxes", box_size, ".png")))
  show_hull(path, threshold = threshold, invert = invert, save = out("_hull.png"))
  invisible(res)
}
