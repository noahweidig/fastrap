# Extract RAP vegetation cover for polygons

Returns the mean percent cover of each Rangeland Analysis Platform (RAP)
vegetation cover v3 band inside every feature of `geometry`.

## Usage

``` r
extract_rap(
  year,
  geometry,
  cover = names(RAP_BANDS),
  units = c("ha", "km2", "m2", "acres", "mi2"),
  id = NULL,
  wide = FALSE,
  cache_dir = NULL,
  chunk_size = 500L,
  workers = 4L,
  cores = 1L,
  quiet = FALSE
)
```

## Arguments

- year:

  Integer, 1986 to 2025.

- geometry:

  An `sf`/`sfc` object, a
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html),
  or a path to any vector file
  [`sf::st_read()`](https://r-spatial.github.io/sf/reference/st_read.html)
  can open. Polygons only. Must fall within CONUS.

- cover:

  Cover types to return, as RAP band codes: any of `"AFG"` (annual forb
  and grass), `"BGR"` (bare ground), `"LTR"` (litter), `"PFG"`
  (perennial forb and grass), `"SHR"` (shrub), `"TRE"` (tree). Case
  insensitive. Defaults to all 6.

- units:

  Area units: one of `"ha"`, `"km2"`, `"m2"`, `"acres"`, `"mi2"`.

- id:

  Optional column name in `geometry` to use as the feature id. Defaults
  to row number.

- wide:

  If `TRUE`, return one row per feature with `pct_<band>` and
  `area_<band>` columns. Default `FALSE` returns long (tidy) format.

- cache_dir:

  Directory for cached tiles. Defaults to a persistent user cache
  (`tools::R_user_dir("fastrap", "cache")`).

- chunk_size:

  Features per chunk. Lower it if memory is tight.

- workers:

  Parallel download workers. Raise it on a fast connection.

- cores:

  Worker processes for extraction (forked; Unix/macOS only, ignored on
  Windows).

- quiet:

  Suppress progress messages.

## Value

A `data.frame`. Long format columns: `id`, `year`, `band` (RAP band
code: `AFG`, `BGR`, `LTR`, `PFG`, `SHR`, `TRE`), `cover_type`, `pct`
(area-weighted mean percent cover over the feature's cells with data),
and `area` (`pct` / 100 times the feature's area, in `units`). Bands do
not sum to 100: RAP predicts each cover type separately.

## Details

Speed comes from four things:

- Only the 2048 x 2048 pixel windows of the national cloud-optimized
  GeoTIFF that touch your features are read, using HTTP range requests.
  Nothing national is downloaded.

- Windows are fetched by parallel workers (works on Windows too) and
  cached on disk, so repeat calls over the same area skip the network.

- Windows snap to the RAP grid, so no resampling happens.

- Extraction uses exact, cell-area-weighted coverage fractions
  (`exactextractr`, C++).

Memory stays flat for large layers: features are sorted spatially and
run in chunks of `chunk_size`, and each chunk reads only its own tiles
through a lazy VRT.

RAP was built for rangelands. Its authors note cover estimates may not
suit forests or agricultural land.

## Examples

``` r
if (FALSE) { # \dontrun{
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"))
extract_rap(year = 2024, geometry = nc[1:5, ], units = "km2", id = "NAME")
extract_rap(2024, nc[1:5, ], cover = c("TRE", "SHR"), id = "NAME")
} # }
```
