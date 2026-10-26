# ============================================================================
# Zebrafish Viewpoint Sleep / Activity / Lights-Off Arousal Pipeline (R)
# ----------------------------------------------------------------------------
# Combines, in one run:
#   PART 2   Raw Viewpoint export  -> wide FISH1..FISHn + CLOCK      (was Viewpoint_to_Matlab.m)
#   PART 3   Crop to a clean N-hour window
#   PART 4   Load wide/cropped data, assign genotype, flag sleep      (was VT_analysis_2019aedited.m)
#   PART 5   Per-fish / per-day / per-phase summary table
#   PART 6   48-hr genotype mean +/- SEM trace plots
#   PART 7   Lights-off dark-flash arousal response (ZT14 window)     (was LightsOFF_ArousalResponse.R)
#   PART 8   Stats on the summary table: boxplots, SSMD, heatmaps
#
# ----------------------------------------------------------------------------
# FOLDER STRUCTURE (same for every gene/experiment -- only PART 1 changes)
#
#   <base_dir>/<gene_id>/
#     csv_files/<raw_csv_name>                 raw Viewpoint export
#     Genotype_files/<genotype_file_name>       genotype/well map
#     R_analysis_output/
#       analysis_output/                       every .csv / .txt this script writes
#       figures/                                every .png / .pdf / .tif this script writes
#
# To move to a new gene/experiment, change ONLY the "USER CONFIG" block in
# PART 1 -- every path below is derived from base_dir + gene_id.
# ============================================================================


## ---- 0. Install & load required packages ----------------------------------

required_packages <- c("data.table", "dplyr", "tidyr", "readr", "ggplot2",
                       "stringr", "pheatmap")

new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)

invisible(lapply(required_packages, library, character.only = TRUE))


# ============================================================================
# PART 1 -- USER CONFIG  (the ONLY section you should need to edit)
# ============================================================================

## --- where this experiment lives -------------------------------------------
## Everything is built from these two lines. Point base_dir at the parent
## "Sleep Data" folder, and gene_id at the experiment subfolder -- e.g. for
## "/path/to/Sleep Data/dpydab/csv_files/dpydabSA120260623.csv" that's
## base_dir = ".../Sleep Data" and gene_id = "dpydab".
base_dir <- "path/to/Sleep Data"   # <-- EDIT: parent folder holding one subfolder per experiment
gene_id  <- "dpydab"

## exact filenames inside csv_files/ and Genotype_files/ (see folder
## structure above) -- these are the only filenames you should need to type
raw_csv_name        <- ""
genotype_file_name  <- "dpydabCombined_genotypeList.txt"

## --- raw-data conversion ----------------------------------------------------
## TRUE = skip PARTs 2-3 and read an already-cropped wide file (must live in
## analysis_output/, see wide_output_name below) straight into PART 4.
skip_raw_conversion <- TRUE

zt0_clock_hour  <- 9      # real wall-clock hour that corresponds to ZT0 (lights-on)
crop_hours      <- 48     # length of the window to crop to
clock_threshold <- 0.02   # a row counts as "CLOCK ~ 0" if CLOCK < this (~1 min tolerance)

wide_output_name <- paste0(gene_id, "_wide_cropped.txt")   # cropped wide file (written/read in analysis_output/)

## --- genotype colors ---------------------------------------------------
## Name each entry EXACTLY as it appears as a column header in the genotype
## list file. Add/remove entries if you have more than 2 genotypes. These
## colors are reused for BOTH the 48-hr trace plots (PART 6) and the
## lights-off arousal response plots (PART 7).
genotype_colors <- c(
  control = "gray50",
  dpydab    = "dodgerblue"
)



## --- genotypes compared in PART 8 (boxplots / SSMD) ------------------------
## Must match genotype column headers in the genotype file exactly.
control_genotype <- "control"
mutant_genotype  <- "dpydab"

## --- sleep / activity metric settings --------------------------------------
sleep_cutoff  <- 0.5   # sec of activity per 1-min bin; at/below this = "asleep"
bin_minutes   <- 60    # bin size for the PART 6 trace plots (60 = hourly; try 10 for 10-min bins)
light_on_zt   <- 0     # ZT hour lights turn ON
light_off_zt  <- 14    # ZT hour lights turn OFF (14/10 LD cycle) -- also used as the
# dark-flash transition ZT for PART 7 below
period_hours  <- 24    # length of one light/dark cycle
day_hours     <- light_off_zt - light_on_zt
night_hours   <- period_hours - day_hours

## --- lights-off dark-flash arousal response window (PART 7) ---------------
lightsoff_pre_min  <- 10   # minutes immediately BEFORE lights-off to average ("Pre")
lightsoff_post_min <- 2    # minutes immediately AT/AFTER lights-off to average ("Post")

## --- stats/measure order for PART 8 boxplot + heatmap ----------------------
custom_order <- c(
  "HourlyActivity_Day",
  "WakingActivity_Day",
  "HourlySleep_Day",
  "HourlySleepBout_Day",
  "SleepLength_Day",
  "HourlyActivity_Night",
  "WakingActivity_Night",
  "HourlySleep_Night",
  "HourlySleepBout_Night",
  "SleepLength_Night",
  "SleepLatency_Night"
)

## --- manuscript figure sizing -----------------------------------------------
## Every figure in this script is saved at one of these three sizes so the
## whole set drops cleanly into subpanels of a multi-panel manuscript figure.
## Tweak once here and every plot updates.
panel_dpi          <- 300
panel_base_size    <- 8      # ggplot theme base font size for panel-sized plots
panel_width_in     <- 3.25   # single trace / single bar-plot panel
panel_height_in    <- 2.5
facet_width_in     <- 7.5    # wider panel for multi-facet plots (boxplot grid, day-facets)
facet_height_in    <- 5
heatmap_width_in   <- 3.0    # narrow tall panel for the SSMD heatmaps
heatmap_height_in  <- 4.5

## ============================================================================
## Everything below this line is derived automatically -- no further path
## edits should be needed when switching genes/experiments.
## ============================================================================

project_dir   <- file.path(base_dir, gene_id)
raw_data_file <- file.path(project_dir, "csv_files", raw_csv_name)
genotype_file <- file.path(project_dir, "Genotype_files", genotype_file_name)

output_dir <- file.path(project_dir, "R_analysis_output")
csv_dir    <- file.path(output_dir, "analysis_output")   # every .csv / .txt this script writes
fig_dir    <- file.path(output_dir, "figures")            # every .png / .pdf / .tif this script writes

dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

wide_output_file <- file.path(csv_dir, wide_output_name)
data_file        <- wide_output_file   # PART 4 always reads from analysis_output/,
# whether PART 2-3 just wrote it or you're skipping to it

## ---- small helpers used throughout ----------------------------------------

## Manuscript-subpanel ggplot theme -- layer on top of any theme_* call.
theme_manuscript_panel <- function(base_size = panel_base_size) {
  theme(
    text              = element_text(size = base_size),
    axis.text         = element_text(size = base_size, color = "black"),
    axis.title        = element_text(size = base_size + 1),
    plot.title        = element_text(size = base_size + 1, face = "bold"),
    legend.title      = element_text(size = base_size),
    legend.text       = element_text(size = base_size - 1),
    legend.key.size   = unit(0.35, "cm"),
    strip.text        = element_text(size = base_size, face = "bold"),
    plot.margin       = margin(4, 8, 4, 4)
  )
}

## Save one plot as both .pdf (vector, for the manuscript) and .png (quick
## preview) into fig_dir, at a consistent subpanel size.
save_panel <- function(plot, filename_base, width = panel_width_in,
                       height = panel_height_in, dpi = panel_dpi) {
  ggsave(file.path(fig_dir, paste0(filename_base, ".pdf")), plot,
         width = width, height = height, units = "in", useDingbats = FALSE)
  ggsave(file.path(fig_dir, paste0(filename_base, ".png")), plot,
         width = width, height = height, units = "in", dpi = dpi)
}


# ============================================================================
# PART 2 -- Raw Viewpoint file -> wide format (FISH1...FISHn + CLOCK)
# ============================================================================

## 96-well location string -> global well index 1-192.
## The letter right after "Loc" identifies the BOX (A = box 1, wells 1-96;
## B = box 2, wells 97-192) -- both boxes are present in the same raw file.
## The trailing 2-digit number is already the well's position within its
## box (01-96), so e.g. "LocA01" -> FISH1, "LocB01" -> FISH97, "LocB96" -> FISH192.
well_index_from_location <- function(loc) {
  m <- regmatches(loc, regexec("Loc([AB])(\\d{2})", loc))
  box_letter <- vapply(m, function(x) x[2], character(1))
  well_num   <- as.integer(vapply(m, function(x) x[3], character(1)))
  box_offset <- ifelse(box_letter == "B", 96L, 0L)
  well_num + box_offset
}

## Robustly parse the raw "sttime" column (real wall-clock time-of-day) into h/m/s,
## whatever type data.table's fread happened to give it (character or ITime).
parse_sttime_to_hms <- function(x) {
  if (inherits(x, "ITime")) {
    hh <- as.integer(x) %/% 3600L
    mm <- (as.integer(x) %% 3600L) %/% 60L
    ss <- as.integer(x) %% 60L
    return(list(h = hh, m = mm, s = ss))
  }
  parts <- do.call(rbind, strsplit(as.character(x), ":"))
  list(h = as.integer(parts[, 1]), m = as.integer(parts[, 2]), s = as.numeric(parts[, 3]))
}

## Real wall-clock time -> ZT/CLOCK hour, treating zt0_clock_hour as ZT=0
## (identical formula to Viewpoint_to_Matlab.m's CLOCK calculation).
clock_from_sttime <- function(sttime, zt0_clock_hour = 9) {
  hms <- parse_sttime_to_hms(sttime)
  raw_hour <- hms$h + hms$m / 60 + hms$s / 3600
  clk <- raw_hour - zt0_clock_hour
  clk[clk < 0] <- clk[clk < 0] + 24
  clk
}

## Convert one raw, long-format Viewpoint CSV (both boxes combined) into a
## wide FISH1...FISH192 + CLOCK table. Uses data.table for a low memory
## footprint -- only the 5 needed columns are ever read into memory,
## regardless of how many extra columns the raw export has.
convert_viewpoint_raw_to_wide <- function(raw_file, zt0_clock_hour = 9) {
  
  needed_cols <- c("location", "start", "end", "sttime", "middur")
  dt <- fread(raw_file, select = needed_cols, showProgress = TRUE)
  
  dt[, well := well_index_from_location(location)]
  
  fish_ids <- sort(unique(dt$well))
  wide <- dcast(dt, start ~ well, value.var = "middur", fun.aggregate = sum)
  setnames(wide, as.character(fish_ids), paste0("FISH", fish_ids))
  
  endtimes <- dt[, .(end = end[1]), by = start]
  sttimes  <- dt[, .(sttime = sttime[1]), by = start]
  
  wide <- merge(wide, endtimes, by = "start")
  wide <- merge(wide, sttimes, by = "start")
  setorder(wide, start)
  
  wide[, CLOCK := clock_from_sttime(sttime, zt0_clock_hour)]
  wide[, sttime := NULL]
  
  fish_cols <- paste0("FISH", fish_ids)
  setcolorder(wide, c("start", "end", fish_cols, "CLOCK"))
  setnames(wide, c("start", "end"), c("starttimes", "endtimes"))
  wide[]
}


# ============================================================================
# PART 3 -- Crop the converted file to a clean N-hour window
# ============================================================================
# The converted file runs longer than the analysis window (extra minutes
# before the first lights-on, and/or extra data after the period of
# interest). This trims it to start at the first CLOCK ~ 0 (first
# lights-on) and keep exactly crop_hours (2880 rows for 48 hr), dropping
# everything before and after.

crop_to_window <- function(wide_dt, crop_hours = 48, clock_threshold = 0.02) {
  start_idx <- which(wide_dt$CLOCK < clock_threshold)[1]
  if (is.na(start_idx)) {
    stop("No row found with CLOCK below clock_threshold - check clock_threshold or zt0_clock_hour.")
  }
  n_rows_needed <- crop_hours * 60
  end_idx <- start_idx + n_rows_needed - 1
  if (end_idx > nrow(wide_dt)) {
    warning(sprintf(
      "Only %d rows available after the crop start (need %d for %g hours) - returning what is available.",
      nrow(wide_dt) - start_idx + 1, n_rows_needed, crop_hours
    ))
    end_idx <- nrow(wide_dt)
  }
  wide_dt[start_idx:end_idx]
}

if (!skip_raw_conversion) {
  
  message("Converting raw file: ", raw_data_file)
  wide_full <- convert_viewpoint_raw_to_wide(raw_data_file, zt0_clock_hour = zt0_clock_hour)
  
  cropped <- crop_to_window(wide_full, crop_hours = crop_hours, clock_threshold = clock_threshold)
  
  fwrite(cropped, wide_output_file, sep = "\t")
  message("Wrote cropped wide-format file (", nrow(cropped), " rows, ",
          length(grep('^FISH', names(cropped))), " wells) to: ",
          normalizePath(wide_output_file))
}


# ============================================================================
# PART 4 -- Load wide/cropped data, assign genotype, flag sleep minutes
# ============================================================================

read_genotype_list <- function(path) {
  raw <- read.delim(path, sep = "\t", skip = 1, header = TRUE,
                    check.names = FALSE, stringsAsFactors = FALSE)
  raw <- raw[, colSums(!is.na(raw)) > 0 & nzchar(trimws(names(raw))), drop = FALSE]
  names(raw) <- trimws(names(raw))
  
  geno_list <- lapply(names(raw), function(g) {
    wells <- suppressWarnings(as.integer(raw[[g]]))
    wells <- wells[!is.na(wells)]
    data.frame(genotype = g, well = wells)
  })
  df <- dplyr::bind_rows(geno_list)
  df$genotype <- factor(df$genotype, levels = names(raw))   # preserves column order
  df
}

geno_map <- read_genotype_list(genotype_file)

wide_data <- readr::read_delim(data_file, delim = "\t", col_types = readr::cols())
fish_cols <- grep("^FISH[0-9]+$", names(wide_data), value = TRUE)
n_min     <- nrow(wide_data)
CLOCK_vec <- wide_data$CLOCK

activity_long <- wide_data %>%
  dplyr::mutate(minute = dplyr::row_number()) %>%
  dplyr::select(minute, dplyr::all_of(fish_cols)) %>%
  tidyr::pivot_longer(dplyr::all_of(fish_cols), names_to = "well", values_to = "activity") %>%
  dplyr::mutate(well = as.integer(stringr::str_remove(well, "FISH"))) %>%
  dplyr::inner_join(geno_map, by = "well") %>%                # drop any unassigned wells
  dplyr::mutate(
    CLOCK     = CLOCK_vec[minute],
    sleep_min = as.integer(activity <= sleep_cutoff)
  )

setDT(activity_long)
setorder(activity_long, well, minute)

total_hours <- n_min / 60


# ============================================================================
# PART 5 -- Full per-fish, per-day, per-phase sleep/activity summary table
# ============================================================================

## phase: 0 = Day/light (CLOCK < light_off_zt), 1 = Night/dark
activity_long[, lightschedule := as.integer((CLOCK - light_on_zt) %% period_hours >= day_hours)]

## which 24-hr block (Day 1, Day 2, ...) each minute belongs to
first_clock <- CLOCK_vec[1]
activity_long[, daynumber := floor((first_clock + (minute - 1) / 60) / period_hours) + 1L]

## sleep bout onset = asleep this minute, awake (or start of recording) the minute before
activity_long[, sleep_prev := shift(sleep_min, 1, fill = 0L), by = well]
activity_long[, boutstart := as.integer(sleep_min == 1L & sleep_prev == 0L)]

## ---- sleep latency: minutes from each light/dark transition until first sleep ----
ls_full <- as.integer((CLOCK_vec - light_on_zt) %% period_hours >= day_hours)
boundary_idx <- c(1L, which(diff(ls_full) != 0) + 1L)

compute_latency_vec <- function(sleep_vec, boundaries) {
  n <- length(sleep_vec)
  lat <- rep(0L, n)
  for (k in boundaries) {
    lat[k] <- 1L
    if (k < n) {
      rest <- sleep_vec[(k + 1):n]
      p <- which(rest == 1L)[1]
      run_len <- if (is.na(p)) (n - k) else (p - 1L)   # minutes awake strictly before first sleep
      if (run_len > 0) {
        end_idx <- min(k + run_len, n)
        lat[(k + 1):end_idx] <- 1L
      }
    }
  }
  lat
}

wells <- sort(unique(activity_long$well))
latency_by_well <- lapply(wells, function(w) {
  sv <- activity_long[well == w][order(minute), sleep_min]
  compute_latency_vec(sv, boundary_idx)
})
names(latency_by_well) <- as.character(wells)

activity_long[, latency := latency_by_well[[as.character(well[1])]][minute], by = well]

## ---- aggregate to one row per fish, per day-block, per phase ----
n_hours_lookup <- c(`0` = day_hours, `1` = night_hours)

agg <- activity_long[, .(
  activity_sum  = sum(activity, na.rm = TRUE),
  sleep_sum     = sum(sleep_min, na.rm = TRUE),
  boutstart_sum = sum(boutstart, na.rm = TRUE),
  latency_sum   = sum(latency, na.rm = TRUE),
  waking_mean   = mean(activity[activity != 0], na.rm = TRUE),
  n_min         = .N
), by = .(genotype, well, daynumber, lightschedule)]

agg[, n_hours          := n_hours_lookup[as.character(lightschedule)]]
agg[, sleeplength      := sleep_sum / boutstart_sum]              # Inf if 0 bouts, matching MATLAB x/0
agg[, hourly_activity  := activity_sum / n_hours]
agg[, hourly_sleep     := sleep_sum / n_hours]
agg[, hourly_boutrate  := boutstart_sum / n_hours]
agg[, phase := ifelse(lightschedule == 0L, "Day", "Night")]

summary_table <- agg %>%
  dplyr::select(genotype, well, daynumber, phase,
                sleep_sum, boutstart_sum, sleeplength, latency_sum,
                hourly_activity, hourly_sleep, waking_mean, hourly_boutrate) %>%
  tidyr::pivot_wider(
    id_cols = c(genotype, well, daynumber),
    names_from = phase,
    values_from = c(sleep_sum, boutstart_sum, sleeplength, latency_sum,
                    hourly_activity, hourly_sleep, waking_mean, hourly_boutrate),
    names_glue = "{.value}_{phase}"
  ) %>%
  dplyr::rename(
    Genotype                = genotype,
    FishID                  = well,
    Day                     = daynumber,
    Sleep_Day               = sleep_sum_Day,
    SleepBout_Day           = boutstart_sum_Day,
    SleepLength_Day         = sleeplength_Day,
    SleepLatency_Day        = latency_sum_Day,
    HourlyActivity_Day      = hourly_activity_Day,
    HourlySleep_Day         = hourly_sleep_Day,
    WakingActivity_Day      = waking_mean_Day,
    HourlySleepBout_Day     = hourly_boutrate_Day,
    Sleep_Night             = sleep_sum_Night,
    SleepBout_Night         = boutstart_sum_Night,
    SleepLength_Night       = sleeplength_Night,
    SleepLatency_Night      = latency_sum_Night,
    HourlyActivity_Night    = hourly_activity_Night,
    HourlySleep_Night       = hourly_sleep_Night,
    WakingActivity_Night    = waking_mean_Night,
    HourlySleepBout_Night   = hourly_boutrate_Night
  ) %>%
  dplyr::arrange(Genotype, Day, FishID) %>%
  dplyr::select(Genotype, Day, FishID,
                Sleep_Day, SleepBout_Day, SleepLength_Day, SleepLatency_Day,
                HourlyActivity_Day, HourlySleep_Day, WakingActivity_Day, HourlySleepBout_Day,
                Sleep_Night, SleepBout_Night, SleepLength_Night, SleepLatency_Night,
                HourlyActivity_Night, HourlySleep_Night, WakingActivity_Night, HourlySleepBout_Night)

summary_table_file <- file.path(csv_dir, "R_summarytable.csv")
readr::write_csv(summary_table, summary_table_file)
message("Full per-fish, per-day summary table written to: ", normalizePath(summary_table_file))


# ============================================================================
# PART 6 -- 48-hr genotype mean +/- SEM trace plots (activity, sleep, waking activity)
# ============================================================================
# Each hourly point is plotted at the START of the hour it summarizes (e.g.
# the point for ZT 0-1 sits at x=0, not x=1), so the trace starts exactly at
# ZT0 instead of appearing to begin one bin late.

activity_long[, bin := ceiling(minute / bin_minutes)]

per_fish_bin <- activity_long[, .(
  bin_time        = min((minute - 1) / 60),   # ZT hour at the START of the bin
  activity_sum    = sum(activity, na.rm = TRUE),
  sleep_sum       = sum(sleep_min, na.rm = TRUE),
  waking_activity = mean(activity[activity != 0], na.rm = TRUE)   # mean activity per active (nonzero) minute
), by = .(genotype, well, bin)]

sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

geno_summary <- per_fish_bin %>%
  dplyr::group_by(genotype, bin, bin_time) %>%
  dplyr::summarise(
    activity_mean = mean(activity_sum, na.rm = TRUE),
    activity_sem  = sem(activity_sum),
    sleep_mean    = mean(sleep_sum, na.rm = TRUE),
    sleep_sem     = sem(sleep_sum),
    waking_mean   = mean(waking_activity, na.rm = TRUE),
    waking_sem    = sem(waking_activity),
    n_fish        = dplyr::n(),
    .groups = "drop"
  )

readr::write_csv(geno_summary, file.path(csv_dir, "hourly_genotype_summary.csv"))

build_dark_rects <- function(light_on, light_off, period, total_hours) {
  n_periods <- ceiling(total_hours / period)
  recs <- lapply(0:(n_periods - 1), function(p) {
    xmin <- p * period + light_off
    xmax <- min((p + 1) * period + light_on, total_hours)
    if (xmin >= total_hours) return(NULL)
    data.frame(xmin = xmin, xmax = xmax)
  })
  dplyr::bind_rows(recs)
}
dark_rects <- build_dark_rects(light_on_zt, light_off_zt, period_hours, total_hours)

plot_trace <- function(summary_df, mean_col, sem_col, y_lab, title, colors,
                       dark_rects, total_hours,
                       light_on_zt, light_off_zt, period_hours) {
  
  breaks_x <- sort(unique(c(light_on_zt, light_off_zt,
                            light_on_zt + period_hours,
                            light_off_zt + period_hours,
                            total_hours)))
  breaks_x <- breaks_x[breaks_x <= total_hours]
  
  ggplot(summary_df, aes(x = bin_time, y = .data[[mean_col]],
                         color = genotype, fill = genotype)) +
    geom_rect(data = dark_rects, inherit.aes = FALSE,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "grey85", alpha = 0.6) +
    geom_ribbon(aes(ymin = .data[[mean_col]] - .data[[sem_col]],
                    ymax = .data[[mean_col]] + .data[[sem_col]]),
                color = NA, alpha = 0.25) +
    geom_line(linewidth = 0.9) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    scale_x_continuous(breaks = breaks_x, limits = c(0, total_hours), expand = c(0, 0)) +
    labs(x = "Zeitgeber Time (h)", y = y_lab, title = title,
         color = "Genotype", fill = "Genotype") +
    theme_classic(base_size = panel_base_size) +
    theme_manuscript_panel() +
    theme(
      legend.position   = "top",
      axis.line         = element_line(linewidth = 0.6),
      axis.ticks.length = unit(0.15, "cm")
    )
}

activity_plot <- plot_trace(
  geno_summary, "activity_mean", "activity_sem",
  y_lab = paste0("Activity (sec / ", bin_minutes, " min)"),
  title = "Activity", colors = genotype_colors,
  dark_rects = dark_rects, total_hours = total_hours,
  light_on_zt = light_on_zt, light_off_zt = light_off_zt, period_hours = period_hours
)

sleep_plot <- plot_trace(
  geno_summary, "sleep_mean", "sleep_sem",
  y_lab = paste0("Sleep (min asleep / ", bin_minutes, " min)"),
  title = "Sleep", colors = genotype_colors,
  dark_rects = dark_rects, total_hours = total_hours,
  light_on_zt = light_on_zt, light_off_zt = light_off_zt, period_hours = period_hours
)

waking_plot <- plot_trace(
  geno_summary, "waking_mean", "waking_sem",
  y_lab = "Waking activity (sec / active min)",
  title = "Waking Activity", colors = genotype_colors,
  dark_rects = dark_rects, total_hours = total_hours,
  light_on_zt = light_on_zt, light_off_zt = light_off_zt, period_hours = period_hours
)

print(activity_plot)
print(sleep_plot)
print(waking_plot)

save_panel(activity_plot, paste0(gene_id, "_Activity_trace_48hr"))
save_panel(sleep_plot,    paste0(gene_id, "_Sleep_trace_48hr"))
save_panel(waking_plot,   paste0(gene_id, "_WakingActivity_trace_48hr"))

message("PART 6 done. Trace figures written to: ", normalizePath(fig_dir))


# ============================================================================
# PART 7 -- Lights-off dark-flash arousal response (was LightsOFF_ArousalResponse.R)
# ============================================================================
# Re-uses wide_data / geno_map / genotype_colors already loaded above -- no
# separate file paths or genotype file needed. Looks at 10
# minutes immediately before the light_off_zt transition ("Pre") and
# 2 minutes at/after it ("Post"), for every light/dark
# transition found in the cropped window. 
# Because the timing of the lights on/off could fall within the minute prior to ZT14,
# I made this look for the closest minute to determine the lights off 
# otherwise the pre window may contain the lights off arousal response.

# ============================================================================
# PART 7 -- Lights-off dark-flash arousal response 
# ============================================================================

## geno_map keyed by FISH-style id, to match wide_data's column names
geno_map_fish <- geno_map %>%
  dplyr::mutate(fish_id = paste0("FISH", well)) %>%
  dplyr::rename(treatment = genotype)

treatment_colors <- genotype_colors   # same palette as PART 6, keyed by the same genotype names
names(treatment_colors) <- names(genotype_colors)

## ---- find every lights-off transition shifted to 1 min bin BEFORE closest ZT value ----
# 1. Identify approximate transition points across cycles
raw_zt_crossings <- which(wide_data$CLOCK[-nrow(wide_data)] < light_off_zt &
                            wide_data$CLOCK[-1] >= light_off_zt) + 1

# 2. For each crossing region, find the index closest to light_off_zt and set idx_off to 1 bin prior
zt_off_indices <- sapply(raw_zt_crossings, function(idx) {
  # Search within a local 3-row window around the crossing to find the absolute minimum difference
  search_window <- max(1, idx - 1):min(nrow(wide_data), idx + 1)
  closest_idx <- search_window[which.min(abs(wide_data$CLOCK[search_window] - light_off_zt))]
  
  # Return the 1-min bin immediately before the value closest to ZT
  return(max(1, closest_idx - 1))
})

message("PART 7: shifted lights-off transition(s) to row index: ",
        paste(zt_off_indices, collapse = ", "))

## ---- Pre/Post mean-activity summary per fish, per transition --------------
summary_list <- list()

for (i in seq_along(zt_off_indices)) {
  idx_off <- zt_off_indices[i]
  
  pre_rows  <- max(1, idx_off - lightsoff_pre_min):(idx_off - 1)
  post_rows <- idx_off:min(nrow(wide_data), idx_off + lightsoff_post_min - 1)
  
  pre_df <- wide_data[pre_rows, ] %>%
    dplyr::select(dplyr::all_of(fish_cols)) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(.x, na.rm = TRUE))) %>%
    tidyr::pivot_longer(cols = dplyr::everything(), names_to = "fish_id", values_to = "activity") %>%
    dplyr::mutate(period = "Pre", Day = paste("Day", i))
  
  post_df <- wide_data[post_rows, ] %>%
    dplyr::select(dplyr::all_of(fish_cols)) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(.x, na.rm = TRUE))) %>%
    tidyr::pivot_longer(cols = dplyr::everything(), names_to = "fish_id", values_to = "activity") %>%
    dplyr::mutate(period = "Post", Day = paste("Day", i))
  
  summary_list[[i]] <- dplyr::bind_rows(pre_df, post_df)
}

plot_data <- dplyr::bind_rows(summary_list) %>%
  dplyr::inner_join(geno_map_fish, by = "fish_id") %>%
  dplyr::mutate(period = factor(period, levels = c("Pre", "Post")))

## ---- continuous minute-by-minute trace across the transition --------------
continuous_list <- list()

for (i in seq_along(zt_off_indices)) {
  idx_off <- zt_off_indices[i]
  
  window_rows <- max(1, idx_off - lightsoff_pre_min):min(nrow(wide_data), idx_off + lightsoff_post_min - 1)
  rel_minutes <- window_rows - idx_off   # 0 = exact minute lights turn off
  
  df_window <- wide_data[window_rows, ] %>%
    dplyr::select(dplyr::all_of(fish_cols)) %>%
    dplyr::mutate(minute_offset = rel_minutes, Day = paste("Day", i)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(fish_cols), names_to = "fish_id", values_to = "activity")
  
  continuous_list[[i]] <- df_window
}

continuous_summary <- dplyr::bind_rows(continuous_list) %>%
  dplyr::inner_join(geno_map_fish, by = "fish_id") %>%
  dplyr::mutate(treatment = factor(treatment, levels = names(treatment_colors))) %>%
  dplyr::group_by(Day, treatment, minute_offset) %>%
  dplyr::summarise(
    mean_act = mean(activity, na.rm = TRUE),
    sem_act  = sd(activity, na.rm = TRUE) / sqrt(sum(!is.na(activity))),
    .groups  = "drop"
  )

p_continuous <- ggplot(
  continuous_summary,
  aes(x = minute_offset, y = mean_act, color = treatment, group = treatment)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.4) +
  geom_errorbar(aes(ymin = mean_act - sem_act, ymax = mean_act + sem_act),
                width = 0.15, linewidth = 0.4) +
  scale_color_manual(values = treatment_colors, drop = FALSE) +
  scale_x_continuous(
    breaks = -lightsoff_pre_min:(lightsoff_post_min - 1),
    labels = c(paste0("-", lightsoff_pre_min:1), paste0("ZT", light_off_zt),
               if (lightsoff_post_min > 1) paste0("+", 1:(lightsoff_post_min - 1)) else NULL)
  ) +
  facet_wrap(~ Day) +
  labs(
    x = "Time Relative to Lights-Off (min)",
    y = "Mean Activity (sec/min)",
    color = "Genotype",
    title = "Lights-Off Arousal Response (Continuous Trace)"
  ) +
  theme_classic(base_size = panel_base_size) +
  theme_manuscript_panel() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    strip.background = element_blank()
  )

print(p_continuous)
save_panel(p_continuous, paste0(gene_id, "_lightsoff_continuous_trace"),
           width = facet_width_in, height = panel_height_in)

## ---- Pre vs Post summary line plot -----------------------------------------
plot_summary <- plot_data %>%
  dplyr::group_by(treatment, period, Day) %>%
  dplyr::summarise(
    mean_act = mean(activity, na.rm = TRUE),
    sem_act  = sd(activity, na.rm = TRUE) / sqrt(dplyr::n()),
    .groups  = "drop"
  )

p_prepost <- ggplot(plot_summary, aes(x = period, y = mean_act, color = treatment, group = treatment)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_act - sem_act, ymax = mean_act + sem_act), width = 0.1) +
  scale_color_manual(values = treatment_colors, drop = FALSE) +
  facet_wrap(~ Day) +
  labs(
    x = paste0("Lights-Off Transition (ZT", light_off_zt, ")"),
    y = "Mean Activity (sec/min)",
    color = "Genotype",
    title = "Lights-Off Response by Genotype"
  ) +
  theme_classic(base_size = panel_base_size) +
  theme_manuscript_panel()

print(p_prepost)
save_panel(p_prepost, paste0(gene_id, "_lightsoff_pre_vs_post"))

## ---- per-fish delta (Post - Pre) + paired stats ----------------------------
fish_deltas <- plot_data %>%
  tidyr::pivot_wider(names_from = period, values_from = activity) %>%
  dplyr::mutate(delta = Post - Pre)

stat_results <- fish_deltas %>%
  dplyr::group_by(Day, treatment) %>%
  dplyr::summarise(
    n_fish     = sum(!is.na(delta)),
    mean_pre   = mean(Pre, na.rm = TRUE),
    mean_post  = mean(Post, na.rm = TRUE),
    mean_delta = mean(delta, na.rm = TRUE),
    sem_delta  = sd(delta, na.rm = TRUE) / sqrt(n_fish),
    p_val      = t.test(Post, Pre, paired = TRUE)$p.value,
    .groups    = "drop"
  ) %>%
  dplyr::group_by(Day) %>%
  dplyr::mutate(
    padj_BH = p.adjust(p_val, method = "BH"),
    significance = dplyr::case_when(
      padj_BH < 0.001 ~ "***",
      padj_BH < 0.01  ~ "**",
      padj_BH < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  ) %>%
  dplyr::ungroup()

delta_summary_table <- stat_results %>%
  dplyr::select(
    Day,
    Genotype = treatment,
    N = n_fish,
    `Mean Pre` = mean_pre,
    `Mean Post` = mean_post,
    `Mean Delta (Post-Pre)` = mean_delta,
    `SEM Delta` = sem_delta,
    `p-value (Paired t-test)` = p_val,
    `Adjusted p-value (BH FDR)` = padj_BH,
    Significance = significance
  )

print(as.data.frame(delta_summary_table))
readr::write_csv(delta_summary_table, file.path(csv_dir, paste0(gene_id, "_lightsoff_delta_stats.csv")))

## ---- delta bar plot with individual fish points ----------------------------
max_y_per_day <- fish_deltas %>%
  dplyr::group_by(Day) %>%
  dplyr::summarise(max_val = max(delta, na.rm = TRUE), .groups = "drop")

plot_delta_summary <- stat_results %>%
  dplyr::left_join(max_y_per_day, by = "Day") %>%
  dplyr::mutate(y_pos = max_val + 0.1 * abs(max_val))

p_delta <- ggplot() +
  geom_col(
    data = plot_delta_summary,
    aes(x = treatment, y = mean_delta, fill = treatment),
    color = "black", width = 0.7, linewidth = 0.3, alpha = 0.6
  ) +
  geom_jitter(
    data = fish_deltas,
    aes(x = treatment, y = delta, color = treatment),
    width = 0.15, size = 1, alpha = 0.7
  ) +
  geom_errorbar(
    data = plot_delta_summary,
    aes(x = treatment, ymin = mean_delta - sem_delta, ymax = mean_delta + sem_delta),
    width = 0.2, linewidth = 0.5, color = "black"
  ) +
  geom_text(
    data = plot_delta_summary,
    aes(x = treatment, y = y_pos, label = significance),
    vjust = 0, size = 3, fontface = "bold"
  ) +
  scale_fill_manual(values = treatment_colors, drop = FALSE) +
  scale_color_manual(values = treatment_colors, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  facet_wrap(~ Day) +
  labs(
    x = "Genotype",
    y = expression(paste(Delta, " Activity (Post - Pre, sec/min)")),
    title = "Lights-Off Arousal Delta Response"
  ) +
  theme_classic(base_size = panel_base_size) +
  theme_manuscript_panel() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    legend.position = "none",
    strip.background = element_blank()
  )

print(p_delta)
save_panel(p_delta, paste0(gene_id, "_lightsoff_delta_bar_plot"),
           width = facet_width_in * 0.6, height = panel_height_in)

message("PART 7 done. Lights-off arousal figures/CSV written to: ",
        normalizePath(fig_dir), " / ", normalizePath(csv_dir))

## NOTE: the original LightsOFF_ArousalResponse.R also re-derived a full
## 48-hr continuous hourly-binned sleep/activity trace at the very end. That
## duplicates what PART 6 above already produces from the same data, so it
## was intentionally left out here rather than generated twice -- use the
## PART 6 Activity/Sleep trace figures for that.


# ============================================================================
# PART 8 -- Stats on the summary table: boxplots, SSMD, heatmaps
# ============================================================================
# Reads the summary table written by PART 5 (so this section can also be
# re-run on its own against an existing R_summarytable.csv in analysis_output/).

data8 <- readr::read_csv(summary_table_file, col_types = readr::cols())

# =========================
# REMOVE UNWANTED COLUMNS (day-only aggregate columns not used downstream)
# =========================
data8 <- data8 %>%
  dplyr::select(-c(
    Sleep_Day,
    SleepBout_Day,
    SleepLatency_Day,
    Sleep_Night,
    SleepBout_Night
  ))

# =========================
# LONG FORMAT
# =========================
long_data <- data8 %>%
  tidyr::pivot_longer(
    cols = -c(Genotype, Day, FishID),
    names_to = "Measure",
    values_to = "Value"
  ) %>%
  dplyr::mutate(Genotype_Day = paste(Genotype, paste0("Day", Day), sep = "_"))

long_data$Genotype_Day <- factor(
  long_data$Genotype_Day,
  levels = c(paste0(control_genotype, "_Day1"), paste0(mutant_genotype, "_Day1"),
             paste0(control_genotype, "_Day2"), paste0(mutant_genotype, "_Day2"))
)

# =========================
# BOXPLOT (all measures)
# =========================
box_plot <- ggplot(long_data, aes(x = Genotype_Day, y = Value, fill = Genotype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.6, size = 0.5) +
  facet_wrap(~ Measure, scales = "free_y") +
  scale_fill_manual(values = setNames(c("lightgray", "dodgerblue"), c(control_genotype, mutant_genotype))) +
  theme_bw(base_size = panel_base_size) +
  theme_manuscript_panel() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

print(box_plot)
save_panel(box_plot, paste0(gene_id, "_boxplot_all_measures"),
           width = facet_width_in, height = facet_height_in)

# =========================
# FUNCTION: CALCULATE SSMD TABLE
# =========================
calc_ssmd <- function(data, test_type = "t") {
  
  results <- list()
  measures <- unique(data$Measure)
  days <- c("1", "2")
  
  for (m in measures) {
    for (d in days) {
      
      df <- data %>% dplyr::filter(Measure == m, Day == d)
      if (length(unique(df$Genotype)) < 2) next
      
      control <- df$Value[df$Genotype == control_genotype]
      mutant  <- df$Value[df$Genotype == mutant_genotype]
      
      control_mean <- mean(control, na.rm = TRUE)
      mutant_mean  <- mean(mutant, na.rm = TRUE)
      control_sd   <- sd(control, na.rm = TRUE)
      mutant_sd    <- sd(mutant, na.rm = TRUE)
      
      # Z-SCORE (uses control mean and control SD)
      z_score <- (mutant_mean - control_mean) / control_sd
      
      # SSMD
      ssmd <- (mutant_mean - control_mean) / sqrt(mutant_sd^2 + control_sd^2)
      
      n_control <- length(control)
      n_mutant  <- length(mutant)
      
      se <- sqrt(
        (mutant_sd^2 / (n_mutant * (mutant_sd^2 + control_sd^2))) +
          (control_sd^2 / (n_control * (mutant_sd^2 + control_sd^2)))
      )
      
      ci_lower <- ssmd - 1.96 * se
      ci_upper <- ssmd + 1.96 * se
      
      if (test_type == "t") {
        pval <- t.test(Value ~ Genotype, data = df)$p.value
      } else {
        pval <- wilcox.test(Value ~ Genotype, data = df, exact = FALSE)$p.value
      }
      
      results[[paste(m, d, sep = "_")]] <- data.frame(
        Measure = m,
        Day = d,
        N_control = n_control,
        N_mutant = n_mutant,
        Control_mean = control_mean,
        Mutant_mean = mutant_mean,
        Z_score = z_score,
        SSMD = ssmd,
        SSMD_CI_lower = ci_lower,
        SSMD_CI_upper = ci_upper,
        p_value = pval
      )
    }
  }
  
  results <- dplyr::bind_rows(results)
  results$FDR <- p.adjust(results$p_value, method = "BH")
  
  results$Measure <- factor(results$Measure, levels = custom_order)
  results <- results %>% dplyr::arrange(Measure, Day)
  
  return(results)
}

# =========================
# RESULTS TABLES
# =========================
results_ttest <- calc_ssmd(long_data, "t")
results_wilc  <- calc_ssmd(long_data, "wilcox")

write.csv(results_ttest, file.path(csv_dir, paste0(gene_id, "_SSMD_results_day_ttest.csv")), row.names = FALSE)
write.csv(results_wilc,  file.path(csv_dir, paste0(gene_id, "_SSMD_results_day_wilcox.csv")), row.names = FALSE)

# =========================
# FUNCTION: HEATMAP
# =========================
make_heatmap <- function(results, file_out_base, width, height) {
  
  wide <- results %>%
    dplyr::select(Measure, Day, SSMD, FDR) %>%
    tidyr::pivot_wider(names_from = Day, values_from = c(SSMD, FDR))
  
  mat <- as.matrix(wide[, c("SSMD_1", "SSMD_2")])
  rownames(mat) <- wide$Measure
  mat <- mat[custom_order, , drop = FALSE]
  
  sig <- matrix("", nrow = nrow(mat), ncol = ncol(mat))
  rownames(sig) <- rownames(mat)
  colnames(sig) <- colnames(mat)
  
  for (i in 1:nrow(wide)) {
    if (!is.na(wide$FDR_1[i]) && wide$FDR_1[i] < 0.05) sig[i, 1] <- "*"
    if (!is.na(wide$FDR_2[i]) && wide$FDR_2[i] < 0.05) sig[i, 2] <- "*"
  }
  sig <- sig[custom_order, , drop = FALSE]
  
  # PDF (vector, for the manuscript) + PNG preview, both at panel size
  pdf(file.path(fig_dir, paste0(file_out_base, ".pdf")), width = width, height = height)
  pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           breaks = seq(-1, 1, length.out = 101),
           display_numbers = sig, fontsize_number = 10,
           fontsize = panel_base_size, border_color = "black")
  dev.off()
  
  png(file.path(fig_dir, paste0(file_out_base, ".png")), width = width, height = height,
      units = "in", res = panel_dpi)
  pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           breaks = seq(-1, 1, length.out = 101),
           display_numbers = sig, fontsize_number = 10,
           fontsize = panel_base_size, border_color = "black")
  dev.off()
}

make_heatmap(results_ttest, paste0(gene_id, "_SSMD_heatmap_ttest"), heatmap_width_in, heatmap_height_in)
make_heatmap(results_wilc,  paste0(gene_id, "_SSMD_heatmap_wilcox"), heatmap_width_in, heatmap_height_in)

message("PART 8 done. Stats CSVs + figures written to: ",
        normalizePath(csv_dir), " / ", normalizePath(fig_dir))

message("\nPipeline complete for gene_id = '", gene_id, "'.")
message("  CSV/TXT outputs -> ", normalizePath(csv_dir))
message("  Figures         -> ", normalizePath(fig_dir))