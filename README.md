# fastrap <img src="man/figures/logo.png" align="right" height="139" alt="fastrap logo" />

Fast per-polygon vegetation cover from the **Rangeland Analysis Platform** (RAP vegetation cover v3, CONUS, 1986–2025, ~30 m). One function.

## Install

```r
# install.packages("remotes")
remotes::install_github("noahweidig/fastrap")
```

## Use

```r
library(fastrap)

nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"))
x <- extract_rap(year = 2023, geometry = nc, units = "km2", id = "NAME")

# only some cover types
extract_rap(2023, nc, cover = c("TRE", "SHR"), id = "NAME", wide = TRUE)
#>          id year   pct_SHR area_SHR  pct_TRE area_TRE
#> 1      Ashe 2023 0.6100578 6.938799 70.43895 801.1728
#> ...
```

`cover` takes any of the RAP band codes:

- `AFG` annual forb and grass
- `BGR` bare ground
- `LTR` litter
- `PFG` perennial forb and grass
- `SHR` shrub
- `TRE` tree

`geometry` takes an `sf`/`sfc` object, a `terra::SpatVector`, or a path to any vector file. `units` is one of `"ha"`, `"km2"`, `"m2"`, `"acres"`, `"mi2"`. Set `wide = TRUE` for one row per feature.

## Why it's fast

- Reads only the 2048 × 2048 pixel windows your features touch from the national cloud-optimized GeoTIFF, using HTTP range requests. Nothing national is downloaded.
- Windows download in parallel (`workers`, default 4, works on Windows) and are cached on disk. A second run over the same area skips the network.
- Windows snap to the RAP grid, so no resampling.
- Exact, cell-area-weighted coverage fractions in C++ (`exactextractr`).
- `cores = 4` splits extraction across forked workers (Unix/macOS).

## Output

Long format: `id`, `year`, `band`, `cover_type`, `pct`, `area`. `pct` is the area-weighted mean percent cover over the feature's cells with data. `area` is `pct / 100` times the feature's area. Bands don't sum to 100, because RAP predicts each cover type separately.

## Caveat

RAP was built for rangelands. Its authors say cover estimates may not suit forests or agricultural land.

## Data

[RAP vegetation cover v3](http://rangeland.ntsg.umt.edu/data/rap/rap-vegetation-cover/v3/README), CC0. Please cite:

Allred, B.W., et al. 2021. Improving Landsat predictions of rangeland fractional cover with multitask learning and uncertainty. *Methods in Ecology and Evolution*. <https://doi.org/10.1111/2041-210x.13564>
