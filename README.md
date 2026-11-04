# ASCAT Ploidy Modifier

An interactive **R Shiny** app for reviewing and manually correcting
[ASCAT](https://github.com/VanLoo-lab/ascat) allele-specific copy-number solutions.

ASCAT picks a purity/ploidy solution automatically, but sometimes that solution isn't the right one,
for example in near-tetraploid samples, low-purity tumours or noisy arrays. ASCAT Ploidy Modifier
shows you the fitted profile next to its *sunrise plot* so you can:

- **Change the sample's purity and ploidy** by clicking on the sunrise plot
- **Edit the copy number of individual segments** (major and minor allele), by typing values or by clicking on the profile
- **Re-segment the data** with a different ASPCF penalty
- **Compare** the original and the modified profiles side by side
- **Save** the corrected result as an ASCAT `.Rda` object you can keep using in R

---

## Installation

### 1. Requirements

- R ≥ 4.1 (RStudio is recommended but not required)
- The R packages listed below

### 2. Install the R packages

```r
# CRAN packages
install.packages(c(
  "shiny", "shinybusy", "shinythemes", "shinyWidgets", "shinycssloaders",
  "shinyjs", "shinyalert", "shinydashboard", "shinyFiles",
  "spatstat.geom", "dipsaus", "dplyr", "RColorBrewer", "devtools"
))

# ASCAT and ASCAT.sc (GitHub)
devtools::install_github("VanLoo-lab/ascat/ASCAT")
devtools::install_github("VanLoo-lab/ASCAT.sc")
```

### 3. Get the app

```bash
git clone https://github.com/iovlaicu/ASCATPloidyModifier.git
```

or download the ZIP from GitHub (**Code → Download ZIP**) and unzip it.

### 4. Run it

From R:

```r
shiny::runApp("path/to/ASCATPloidyModifier")
```

Or open `app.R` in RStudio and click **Run App**. The app opens in your browser.

---

## Input data

The app loads a standard ASCAT R data file (`.Rda` / `.RData`) that contains these two objects:

| Object         | What it is                                                      |
|----------------|-----------------------------------------------------------------|
| `ascat.bc`     | The ASCAT data object (LogR/BAF, segmented data), e.g. after `ascat.aspcf()` |
| `ascat.output` | The result of `ascat.runAscat()`                                |

For example, at the end of a normal ASCAT run:

```r
ascat.bc     <- ascat.aspcf(ascat.bc, ...)
ascat.output <- ascat.runAscat(ascat.bc, ...)
save(ascat.bc, ascat.output, file = "mySample_ASCAT.Rda")
```

> **Notes**
> - Each file should contain **one sample**.
> - Plots use the **hg38** reference genome.

---

## How to use

### 1. Welcome tab: load your data
1. Click **Load data** and choose your `.Rda` file.
2. Click **Start**. Loading the object takes a few seconds, and then the app opens the **Modifier** tab.

Click **Help** on the Welcome tab for a short in-app guide.

### 2. Segmentation tab (optional): re-segment the profile
If the segmentation looks too fragmented or too coarse:
1. Choose a **Segmentation penalty (ASPCF)**: 25, 50, 70 (default), 100, 150 or 200.
   A higher penalty gives fewer, longer segments.
2. Click **Apply segmentation**. ASCAT re-segments the raw data and fits a new purity/ploidy solution, which is shown as a preview.
3. Click **Save segmentation** to keep it. It becomes the new working profile on the Modifier tab.
   Click **Discard** to drop the preview.

### 3. Modifier tab: fix the solution
The top panel (**Original**) shows the profile ASCAT produced. Changes are made in the
bottom panel, the **Modifier Working Station**, which shows a copy-number profile and its sunrise plot.

- **Change purity and ploidy:** click anywhere on the sunrise plot to refit the whole sample at that
  purity/ploidy. With **Automatic optima** turned on (the default), the app jumps to the nearest
  local optimum.
- **Change one segment:** open the modify menu, choose *Modify the copy number of 1 segment*,
  choose the chromosome, enter the **major** and **minor** allele copy numbers and click **Apply**.
- **Change a segment on the graph:** choose *Modify segment on graph*, click a segment on the
  Working Station profile and then click where it should be. Repeat for a second segment and click **Apply**.
  The app refits purity/ploidy so both segments match.
- **Display options** (plot settings menu): show raw (unrounded) segments, hide the main
  segments, and overlay the raw **LogR** or **BAF**.
- **Reset profile:** go back to the original ASCAT solution.

### 4. Save the result
Click **Save .Rda** to download `result_manualfitting.Rda`. It contains:

- `ascat.output$original`: the original ASCAT solution
- `ascat.output$manual`: your modified solution (purity, ploidy, segments, …)
- `ascat.bc`: the data object, including any new segmentation you saved

You can load the saved file back into the app to continue editing, or use it in R:

```r
load("result_manualfitting.Rda")
ascat.output$manual$purity
ascat.output$manual$ploidy
head(ascat.output$manual$segments)
```

---

## Repository contents

| File | Purpose |
|------|---------|
| `app.R` | The Shiny app (UI and server) |
| `PlotASCATAdjusted.R` | Plots the copy-number profile (with highlighting and raw-segment options) |
| `plotSunriseASCAT.R` | Plots the sunrise (purity × ploidy goodness-of-fit) and finds local optima |
| `PlotLogRASCAT.R` | Plots the raw LogR / BAF overlays |
| `plotSegmentedASCAT.R` | Segmented-data plotting helpers |
| `refitProfile_ASCAT.R` | Refits and merges segments after manual changes |
| `findMinima.R` | Finds local minima on the sunrise distance matrix |
| `www/` | Images and loading animations used by the interface |

---

## Troubleshooting

- **"there is no package called …"**: install the missing package (see *Installation*).
- **Nothing happens after Start**: check that the file contains objects named exactly
  `ascat.bc` and `ascat.output`, and that it holds a single sample.
- **"Cannot resegment profile with this penalty value"**: try a different penalty.
- **Plots look misplaced**: the app assumes the **hg38** reference genome.

Found a bug or have a suggestion? Please open an [issue](https://github.com/iovlaicu/ASCATPloidyModifier/issues).

---

## Citation

If you use this tool, please also cite ASCAT:

> Van Loo P, *et al.* Allele-specific copy number analysis of tumors. *PNAS* 107(39):16910–16915 (2010).
