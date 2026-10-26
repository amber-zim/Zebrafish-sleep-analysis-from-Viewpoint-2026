# Zebrafish Sleep Pipeline

An end-to-end R pipeline for larval zebrafish sleep/activity experiments recorded on a **ViewPoint** tracking system. It goes from the raw ViewPoint export to a per-fish summary table, genotype trace plots, lights-off arousal analysis, and manuscript-sized statistics figures, all in one script.

## What it does

| Part | Step | Output |
|------|------|--------|
| 2 | Convert raw long-format ViewPoint CSV to wide `FISH1..FISHn` + `CLOCK` | wide table |
| 3 | Crop to a clean N-hour window starting at first lights-on (ZT0) | `<gene>_wide_cropped.txt` |
| 4 | Load data, assign genotypes from the well map, flag sleep minutes | in-memory |
| 5 | Per-fish / per-day / per-phase summary (sleep, bouts, bout length, latency, activity, waking activity) | `R_summarytable.csv` |
| 6 | 48-hr genotype mean ± SEM traces (activity, sleep, waking activity) | `*_trace_48hr.pdf/png`, `hourly_genotype_summary.csv` |
| 7 | Lights-off dark-flash arousal response (Pre vs Post, delta, paired stats) | `*_lightsoff_*` figures, `*_lightsoff_delta_stats.csv` |
| 8 | Boxplots, SSMD, Z-score, t-test / Wilcoxon with BH-FDR | `*_boxplot_all_measures`, `*_SSMD_results_day_*.csv`, heatmaps |

**Sleep definition:** a minute is "asleep" if activity ≤ `sleep_cutoff` (default 0.5 s per 1-min bin). Sleep bouts start on a minute that is asleep after an awake minute. Sleep latency is the number of minutes from each light/dark transition to the first sleep minute.

**Lights-off detection (Part 7):** the transition is snapped to the closest 1-min bin so the "Pre" window (default 10 min) doesn't include the lights-off startle.

## Requirements

- R ≥ 4.1
- Packages (auto-installed on first run): `data.table`, `dplyr`, `tidyr`, `readr`, `ggplot2`, `stringr`, `pheatmap`

## Folder layout

Every experiment follows the same structure, so only the config block changes between genes:

```
<base_dir>/<gene_id>/
├── csv_files/            raw ViewPoint export
├── Genotype_files/       genotype / well map
└── R_analysis_output/
    ├── analysis_output/  all .csv / .txt written by the script
    └── figures/          all .pdf / .png written by the script
```

## Input formats

**Raw ViewPoint CSV** needs these columns: `location`, `start`, `end`, `sttime`, `middur`. Location strings look like `LocA01`; box A = wells 1-96, box B = wells 97-192. See `example_data/csv_files/`.

**Genotype file** is tab-delimited. The first line is skipped (title/notes). The second line has one genotype name per column, followed by the well numbers for that genotype. See `example_data/Genotype_files/`.

```
Genotype list (skipped)
control   dpydab
1         9
2         10
...
```

## Usage

1. Open `R/zebrafish_sleep_pipeline.R` and edit **PART 1: USER CONFIG** only:
   - `base_dir`, `gene_id`, `raw_csv_name`, `genotype_file_name`
   - `skip_raw_conversion`: `FALSE` for the first run on a raw CSV; `TRUE` to reuse an already-cropped wide file
   - `zt0_clock_hour`, `crop_hours`, `genotype_colors`, `control_genotype`, `mutant_genotype`
   - sleep cutoff, light schedule, and figure sizes if needed
2. Run it:

```r
source("R/zebrafish_sleep_pipeline.R")
```
or `Rscript R/zebrafish_sleep_pipeline.R`.

## Notes and limitations

- Designed for two genotypes (control vs. one mutant) in Part 8. Part 6/7 colors support more via `genotype_colors`.
- Part 8 compares Days 1 and 2 (a 48-hr window).
- Sleep length is `Inf`/`NaN` for fish with zero bouts in a phase (matching the original MATLAB behavior).

## Lineage

Consolidates earlier MATLAB scripts (`Viewpoint_to_Matlab.m`, `VT_analysis_2019aedited.m`) and `LightsOFF_ArousalResponse.R`.

## Citation

If you use this pipeline, please cite this repository (add a `CITATION.cff` or Zenodo DOI once released).

## License

MIT, see `LICENSE`.
