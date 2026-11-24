RAP_URL <- "https://rangeland.ntsg.umt.edu/data/rap/rap-vegetation-cover/v3/vegetation-cover-v3-%d.tif"
RAP_NODATA <- 255
RAP_YEARS <- 1986:2025
# RAP v3 grid (EPSG:4326, upper-left corner) and size, from the COG header
GRID_X0 <- -124.76979458874112
GRID_Y0 <- 49.393776065783484
RES <- 0.0002694945852358566
NCOL <- 214847L
NROW <- 92304L
TILE_PX <- 2048L
TILE_D <- TILE_PX * RES

RAP_BANDS <- c(
  AFG = "Annual forb and grass", BGR = "Bare ground", LTR = "Litter",
  PFG = "Perennial forb and grass", SHR = "Shrub", TRE = "Tree"
)

GDAL_CFG <- c(
  GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR",
  CPL_VSIL_CURL_ALLOWED_EXTENSIONS = ".tif",
  GDAL_HTTP_MULTIRANGE = "YES",
  GDAL_HTTP_MERGE_CONSECUTIVE_RANGES = "YES",
  GDAL_HTTP_MAX_RETRY = "5",
  GDAL_HTTP_RETRY_DELAY = "2",
  GDAL_HTTP_TIMEOUT = "120",
  VSI_CACHE = "TRUE"
)

read_geom <- function(x) {
  if (is.character(x) && length(x) == 1L) x <- sf::st_read(x, quiet = TRUE)
  if (inherits(x, "SpatVector")) x <- sf::st_as_sf(x)
  if (inherits(x, "sfc")) x <- sf::st_sf(geometry = x)
  if (!inherits(x, "sf")) stop("`geometry` must be sf, sfc, SpatVector, or a file path.", call. = FALSE)
  if (is.na(sf::st_crs(x))) stop("`geometry` has no CRS.", call. = FALSE)
  types <- unique(as.character(sf::st_geometry_type(x)))
  if (!all(types %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop("`geometry` must be POLYGON/MULTIPOLYGON.", call. = FALSE)
  }
  x
}

# Per-feature tile index ranges on the RAP pixel grid, clamped to the raster.
feature_bboxes <- function(g) {
  b <- vapply(g, function(p) {
    if (length(p) == 0L) return(rep(NA_real_, 4))
    as.numeric(sf::st_bbox(p))
  }, numeric(4))
  mc <- (NCOL - 1L) %/% TILE_PX
  mr <- (NROW - 1L) %/% TILE_PX
  cl <- function(v, hi) pmin(pmax(v, 0), hi)
  bb <- cbind(
    col0 = cl(floor((b[1, ] - GRID_X0) / TILE_D), mc),
    col1 = cl(floor((b[3, ] - GRID_X0) / TILE_D), mc),
    row0 = cl(floor((GRID_Y0 - b[4, ]) / TILE_D), mr),
    row1 = cl(floor((GRID_Y0 - b[2, ]) / TILE_D), mr)
  )
  out <- b[3, ] < GRID_X0 | b[1, ] > GRID_X0 + NCOL * RES |
    b[2, ] > GRID_Y0 | b[4, ] < GRID_Y0 - NROW * RES
  bb[out %in% TRUE, ] <- NA
  bb
}

tiles_for <- function(bb) {
  keys <- unlist(lapply(seq_len(nrow(bb)), function(i) {
    cc <- bb[i, "col0"]:bb[i, "col1"]
    rr <- bb[i, "row0"]:bb[i, "row1"]
    paste(rep(cc, times = length(rr)), rep(rr, each = length(cc)), sep = "_")
  }), use.names = FALSE)
  unique(keys)
}

is_tiff <- function(f) {
  if (!file.exists(f) || file.size(f) < 8) return(FALSE)
  m <- readBin(f, "raw", 4)
  identical(m, as.raw(c(0x49, 0x49, 0x2a, 0x00))) || identical(m, as.raw(c(0x4d, 0x4d, 0x00, 0x2a)))
}

# Reads one tile window from the remote COG with GDAL range requests.
fetch_one <- function(key, src, dst, cfg) {
  cr <- as.integer(strsplit(key, "_", fixed = TRUE)[[1]])
  x0 <- cr[1] * 2048L
  y0 <- cr[2] * 2048L
  w <- min(2048L, 214847L - x0)
  h <- min(2048L, 92304L - y0)
  tmp <- paste0(dst, ".part")
  ok <- tryCatch({
    sf::gdal_utils("translate", src, tmp, options = c(
      "-of", "GTiff", "-srcwin", x0, y0, w, h,
      "-co", "COMPRESS=LZW", "-co", "TILED=YES"
    ), quiet = TRUE, config_options = cfg)
    file.rename(tmp, dst)
  }, error = function(e) FALSE)
  if (!isTRUE(ok)) unlink(tmp)
  isTRUE(ok)
}

fetch_tiles <- function(keys, year, dir, workers = 4L, tries = 4L) {
  files <- file.path(dir, paste0(keys, ".tif"))
  src <- paste0("/vsicurl/", sprintf(RAP_URL, year))
  cl <- NULL
  on.exit(if (!is.null(cl)) parallel::stopCluster(cl))
  for (attempt in seq_len(tries)) {
    need <- which(!vapply(files, is_tiff, logical(1)))
    if (!length(need)) return(files)
    if (workers > 1L && length(need) > 1L) {
      if (is.null(cl)) cl <- parallel::makeCluster(min(workers, length(need)))
      parallel::clusterMap(cl, fetch_one, keys[need], files[need],
                           MoreArgs = list(src = src, cfg = GDAL_CFG))
    } else {
      for (i in need) fetch_one(keys[i], src, files[i], GDAL_CFG)
    }
    if (attempt < tries && !all(vapply(files, is_tiff, logical(1)))) {
      Sys.sleep(min(60, 2^(attempt + 1)))
    }
  }
  bad <- !vapply(files, is_tiff, logical(1))
  if (any(bad)) {
    stop(sum(bad), " RAP tile(s) failed to download after ", tries,
         " tries. Re-run to resume (downloaded tiles are cached).", call. = FALSE)
  }
  files
}

to_long <- function(res, ids, year, area, bands) {
  v <- as.vector(t(as.matrix(res)))
  data.frame(
    id = rep(ids, each = length(bands)),
    year = year,
    band = rep(bands, times = length(ids)),
    cover_type = rep(unname(RAP_BANDS[bands]), times = length(ids)),
    pct = v,
    area = v / 100 * rep(area, each = length(bands)),
    stringsAsFactors = FALSE
  )
}

to_wide <- function(df, ids, bands) {
  w <- data.frame(id = ids, year = if (nrow(df)) df$year[1] else NA_integer_)
  i <- match(df$id, ids)
  for (b in bands) {
    s <- df$band == b
    p <- a <- rep(NA_real_, length(ids))
    p[i[s]] <- df$pct[s]
    a[i[s]] <- df$area[s]
    w[[paste0("pct_", b)]] <- p
    w[[paste0("area_", b)]] <- a
  }
  w
}
