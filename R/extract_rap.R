#' Extract RAP vegetation cover for polygons
#'
#' Returns the mean percent cover of each Rangeland Analysis Platform (RAP)
#' vegetation cover v3 band inside every feature of `geometry`.
#'
#' Speed comes from four things:
#' * Only the 2048 x 2048 pixel windows of the national cloud-optimized
#'   GeoTIFF that touch your features are read, using HTTP range requests.
#'   Nothing national is downloaded.
#' * Windows are fetched by parallel workers (works on Windows too) and
#'   cached on disk, so repeat calls over the same area skip the network.
#' * Windows snap to the RAP grid, so no resampling happens.
#' * Extraction uses exact, cell-area-weighted coverage fractions
#'   (`exactextractr`, C++).
#'
#' Memory stays flat for large layers: features are sorted spatially and run
#' in chunks of `chunk_size`, and each chunk reads only its own tiles through a
#' lazy VRT.
#'
#' RAP was built for rangelands. Its authors note cover estimates may not suit
#' forests or agricultural land.
#'
#' @param year Integer, 1986 to 2025.
#' @param geometry An `sf`/`sfc` object, a `terra::SpatVector`, or a path to
#'   any vector file `sf::st_read()` can open. Polygons only. Must fall within
#'   CONUS.
#' @param cover Cover types to return, as RAP band codes: any of `"AFG"`
#'   (annual forb and grass), `"BGR"` (bare ground), `"LTR"` (litter), `"PFG"`
#'   (perennial forb and grass), `"SHR"` (shrub), `"TRE"` (tree). Case
#'   insensitive. Defaults to all 6.
#' @param units Area units: one of `"ha"`, `"km2"`, `"m2"`, `"acres"`, `"mi2"`.
#' @param id Optional column name in `geometry` to use as the feature id.
#'   Defaults to row number.
#' @param wide If `TRUE`, return one row per feature with `pct_<band>` and
#'   `area_<band>` columns. Default `FALSE` returns long (tidy) format.
#' @param cache_dir Directory for cached tiles. Defaults to a persistent user
#'   cache (`tools::R_user_dir("fastrap", "cache")`).
#' @param chunk_size Features per chunk. Lower it if memory is tight.
#' @param workers Parallel download workers. Raise it on a fast connection.
#' @param cores Worker processes for extraction (forked; Unix/macOS only,
#'   ignored on Windows).
#' @param quiet Suppress progress messages.
#'
#' @return A `data.frame`. Long format columns: `id`, `year`, `band` (RAP
#'   band code: `AFG`, `BGR`, `LTR`, `PFG`, `SHR`, `TRE`), `cover_type`,
#'   `pct` (area-weighted mean percent cover over the feature's cells with
#'   data), and `area` (`pct` / 100 times the feature's area, in `units`).
#'   Bands do not sum to 100: RAP predicts each cover type separately.
#' @export
#' @examples
#' \dontrun{
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"))
#' extract_rap(year = 2024, geometry = nc[1:5, ], units = "km2", id = "NAME")
#' extract_rap(2024, nc[1:5, ], cover = c("TRE", "SHR"), id = "NAME")
#' }
extract_rap <- function(year, geometry, cover = names(RAP_BANDS), units = c("ha", "km2", "m2", "acres", "mi2"),
                        id = NULL, wide = FALSE, cache_dir = NULL,
                        chunk_size = 500L, workers = 4L, cores = 1L, quiet = FALSE) {
  units <- match.arg(units)
  year <- as.integer(year)
  if (length(year) != 1L || is.na(year) || !year %in% RAP_YEARS) {
    stop("`year` must be a single year from 1986 to 2025.", call. = FALSE)
  }
  cover <- unique(toupper(cover))
  bad <- setdiff(cover, names(RAP_BANDS))
  if (length(bad) || !length(cover)) {
    stop("`cover` must be one or more of: ", paste(names(RAP_BANDS), collapse = ", "),
         ". Unknown: ", paste(bad, collapse = ", "), call. = FALSE)
  }
  cover <- names(RAP_BANDS)[names(RAP_BANDS) %in% cover]
  if (is.null(cache_dir)) cache_dir <- tools::R_user_dir("fastrap", "cache")
  ydir <- file.path(cache_dir, "v3", year)
  dir.create(ydir, recursive = TRUE, showWarnings = FALSE)

  g <- read_geom(geometry)
  ids <- if (is.null(id)) seq_len(nrow(g)) else {
    if (!id %in% names(g)) stop("`id` column '", id, "' not found.", call. = FALSE)
    g[[id]]
  }
  g <- sf::st_geometry(g)
  n <- length(g)
  if (n == 0L) stop("`geometry` has no features.", call. = FALSE)
  conv <- c(m2 = 1, ha = 1e-4, km2 = 1e-6, acres = 1 / 4046.8564224, mi2 = 1 / 2589988.110336)[[units]]
  farea <- as.numeric(sf::st_area(sf::st_transform(g, 5070))) * conv
  g <- sf::st_transform(g, 4326)

  bb <- feature_bboxes(g)
  ok <- stats::complete.cases(bb)
  if (!any(ok)) stop("No features overlap the RAP extent (CONUS).", call. = FALSE)
  ord <- order(bb[, "row0"], bb[, "col0"])
  ord <- ord[ok[ord]]
  chunks <- split(ord, ceiling(seq_along(ord) / chunk_size))
  out <- vector("list", length(chunks))
  t0 <- Sys.time()

  all_tiles <- tiles_for(bb[ord, , drop = FALSE])
  if (!quiet) message("[extract_rap] ", length(all_tiles), " tiles (",
                      sum(!file.exists(file.path(ydir, paste0(all_tiles, ".tif")))),
                      " to download)")
  fetch_tiles(all_tiles, year, ydir, workers = max(1L, as.integer(workers)))

  run_chunk <- function(k) {
    idx <- chunks[[k]]
    files <- file.path(ydir, paste0(tiles_for(bb[idx, , drop = FALSE]), ".tif"))
    vrt <- tempfile(fileext = ".vrt")
    on.exit(unlink(vrt))
    r <- terra::vrt(files, filename = vrt)
    terra::NAflag(r) <- RAP_NODATA
    names(r) <- names(RAP_BANDS)
    r <- r[[cover]]
    res <- exactextractr::exact_extract(
      r, sf::st_sf(geometry = g[idx]), fun = "weighted_mean", weights = "area",
      max_cells_in_memory = 3e7, progress = FALSE, full_colnames = TRUE
    )
    res <- res[paste0("weighted_mean.", cover)]
    to_long(res, ids[idx], year, farea[idx], cover)
  }

  cores <- max(1L, as.integer(cores))
  if (cores > 1L && .Platform$OS.type == "unix") {
    out <- parallel::mclapply(seq_along(chunks), run_chunk, mc.cores = cores,
                              mc.preschedule = FALSE)
    err <- vapply(out, inherits, logical(1), "try-error")
    if (any(err)) stop(out[[which(err)[1]]], call. = FALSE)
  } else {
    for (k in seq_along(chunks)) {
      out[[k]] <- run_chunk(k)
      if (!quiet) {
        message(sprintf("\r[extract_rap] %d/%d features | %.1fs",
                        sum(lengths(chunks[seq_len(k)])), length(ord),
                        as.numeric(Sys.time() - t0, units = "secs")), appendLF = FALSE)
      }
    }
  }
  if (!quiet) message("")
  if (any(!ok)) warning(sum(!ok), " empty or out-of-extent geometries skipped.", call. = FALSE)

  df <- do.call(rbind, out)
  df <- df[order(match(df$id, ids), match(df$band, names(RAP_BANDS))), , drop = FALSE]
  rownames(df) <- NULL
  if (wide) df <- to_wide(df, ids[ok], cover)
  attr(df, "units") <- units
  df
}
