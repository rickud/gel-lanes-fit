# Gel Lanes Fit user guide

Gel Lanes Fit is a Fiji (ImageJ) plugin for quantitative analysis of gel electrophoresis images. It takes a line profile from each lane, fits Gaussian peaks and a polynomial background to it, and sizes the bands against a ladder lane.

## Contents

- [Requirements](#requirements)
- [Installing](#installing)
- [Updating and uninstalling](#updating-and-uninstalling)
- [Checking the installation](#checking-the-installation)
- [Using the plugin](#using-the-plugin)
  - [1. Set up the lanes](#1-set-up-the-lanes)
  - [2. Choose the ladder](#2-choose-the-ladder)
  - [3. Choose the fit type and parameters](#3-choose-the-fit-type-and-parameters)
  - [4. Run the fit](#4-run-the-fit)
  - [5. Refine with custom peaks](#5-refine-with-custom-peaks)
- [The two fit types](#the-two-fit-types)
- [Fit parameters](#fit-parameters)
- [Reading the plots](#reading-the-plots)
- [Reading the results](#reading-the-results)
- [Where your results are saved](#where-your-results-are-saved)
- [Ladders and fragment distributions](#ladders-and-fragment-distributions)
- [Troubleshooting](#troubleshooting)

## Requirements

- **Fiji**, up to date. Download it from [fiji.sc](https://fiji.sc) if you don't have it, and run **Help › Update…** to update an existing installation. Gel Lanes Fit needs Java 11 or newer, which current Fiji includes. A very old Fiji still running on Java 8 can't run it (see [Troubleshooting](#fiji-reports-unsupportedclassversionerror)).
- **The Gel Lanes Fit jar**, `gel-lanes-fit-<version>.jar`, from the project's [Releases page](https://github.com/rickud/gel-lanes-fit/releases).

Everything else the plugin uses (JFreeChart, iText, JFreeSVG and the Apache Commons libraries) ships with Fiji, so there's nothing else to install.

## Installing

1. **Quit Fiji** if it's running.
2. **Find the Fiji installation folder.** It's the folder you unpacked Fiji into, usually called `Fiji` or `Fiji.app`, and it contains a `jars` and a `plugins` folder. On macOS, if you only see an app called `Fiji.app`, right-click it and choose **Show Package Contents** to get inside.
3. **Remove any older copy** of the plugin: look in both `jars` and `plugins` for files named `gel-lanes-fit-*.jar` or `GaussianFit-*.jar` (the plugin's name before version 1.0.5) and delete them. Two copies at once can make Fiji run the old one.
4. **Copy `gel-lanes-fit-<version>.jar` into the `jars` folder.**
5. **Start Fiji.** The plugin is in **Plugins › Gel Tools › Gel Lanes Fit**.

## Updating and uninstalling

- **To update**, quit Fiji, replace the old `gel-lanes-fit-*.jar` in `jars` with the new one, and start Fiji again. Your saved lanes, settings and custom peaks are kept, because they're stored next to your images, not in the plugin.
- **To uninstall**, quit Fiji and delete `gel-lanes-fit-*.jar` from `jars`. Your results next to your images are not affected.

## Checking the installation

1. Open a gel image with **File › Open…**. Grayscale images (8- or 16-bit) work best; convert a color image with **Image › Type › 8-bit** first.
2. Choose **Plugins › Gel Tools › Gel Lanes Fit**. The main window's title shows the version, for example `[1.0.8] Gel Lanes Fit - my gel.tif`, so you can confirm Fiji is running the version you installed.
3. The plugin opens its main window and a **Profiles** window with one plot per lane.

Fiji's **Console** (**Window › Console**) shows the plugin's messages, such as the data folder it uses. It's the first place to look when something doesn't work.

## Using the plugin

Open your gel image, then choose **Plugins › Gel Tools › Gel Lanes Fit**. The plugin works on one image at a time and opens two windows:

- **The main window**, with the lane, ladder and fit settings, and the **Fit**, **Open Data Folder** and **Close** buttons.
- **The Profiles window**, with one plot per lane, in tabs of four lanes each.

The lanes are drawn on the image as rectangles. For each lane, the plugin averages the intensity across the lane's width, row by row, which gives the lane's **profile**: intensity against distance down the gel. Distances are in pixels, measured from the top of the image.

If you've analysed the image before, the plugin restores your lanes, ladder, settings and custom peaks, and repeats the last fit (see [Where your results are saved](#where-your-results-are-saved)).

A typical analysis has five steps.

### 1. Set up the lanes

Choose how lanes are placed at the top of the main window:

- **Automatic Rectangle Selection** places equally sized, evenly spaced lanes. Set **Number of Lanes**, then adjust the sliders until the rectangles cover your lanes:
  - **Width** and **Height** of each lane, in pixels;
  - **Space** between lanes;
  - **Horizontal Offset** and **Vertical Offset** of the first lane from the image's top-left corner.
- **Manual Rectangle Selection** lets you draw each lane yourself, which suits bent or unevenly spaced lanes:
  - **Add a lane:** draw a rectangle on the image with ImageJ's rectangle tool.
  - **Move or resize a lane:** drag it or its handles.
  - **Delete a lane:** click inside it (not on a handle) and confirm.

  Lanes are numbered from left to right. Your manual lanes are kept when you switch to Automatic and back.

Make each lane a little narrower than the band so the profile isn't diluted by the gaps between lanes, and tall enough to cover every band you want to fit.

Moving the mouse over a lane highlights its plot in the Profiles window and draws a line at the same position on it.

### 2. Choose the ladder

1. In **Select Ladder Lane**, pick the lane with the size ladder. Its plot turns light blue.
2. In the ladder type list, pick your ladder: **Hi-Lo**, **100bp**, **Quick-Load**, **Tapestation**, or **Custom Ladder** to load your own from a text file (see [Ladders and fragment distributions](#ladders-and-fragment-distributions)).
3. In the **Ladder Range** dialog, pick the **First Band** and **Last Band** that are visible in the ladder lane, in either order. The fit expects to find exactly these bands, so leave out bands that ran off the gel or aren't resolved.

The list then shows the chosen range, for example `Hi-Lo [1.4 kbp - 750 bp]`. To change the range, pick the ladder type again.

### 3. Choose the fit type and parameters

- **Fit Type:** **Banded** for distinct bands, or **Continuum** for smears such as fragmented or tagmented DNA. For Continuum, also choose a **fragment distribution** in the list below the ladder. [The two fit types](#the-two-fit-types) explains the difference.
- **Parameters:** the defaults are a good start. [Fit parameters](#fit-parameters) explains what each one does and when to change it.

### 4. Run the fit

Click **Fit**. The plugin:

1. **Fits the ladder lane first**, always with the Banded method. It must find exactly as many bands as the ladder range contains. If it doesn't, it stops and shows how many bands it detected and how many the range has. To fix it:
   - lower **Peak Tolerance** to detect weaker bands, or raise it to ignore noise;
   - narrow the ladder range to the bands you can actually see;
   - or add the missing bands as custom peaks (step 5);

   then click **Fit** again.
2. **Labels the ladder bands**, with their sizes, as dashed vertical lines on every lane's plot.
3. **Fits all other lanes** with the chosen fit type.
4. **Shows and saves the results:** the curves in the plots, the **Results Display** table, the **LOG** summary, and the files in the data folder.

If you change the lanes or re-run a fit, the plugin warns that the current fit will be replaced. Tick **Don't show this again** to skip that warning in future.

Tick **Show Bands** to mark the bands on the image as short ticks at each lane's edges: magenta for fitted bands, blue for the initial guesses, green for custom peaks. The setting is remembered.

### 5. Refine with custom peaks

When the automatic detection misses a band or finds one where there isn't any, add your own **custom peaks**:

1. Click **Edit Custom Peaks**. It's available after the first fit. In edit mode, the plot of the lane under the mouse on the image has a light green background.
2. **Add a peak:** click on a lane's plot at the band's position. A dialog shows the position and intensity, and asks for the band's width (**FWHM**, the full width at half maximum, in pixels). A custom peak within 2 pixels of a detected band replaces it; otherwise it's added.
3. **Remove a peak:** click within 2 pixels of it.
4. Custom peaks show as green dots. Click **Fit** again to use them.

To remove the custom peaks from some lanes, click **Reset Custom Peaks** and select the lanes. Custom peaks are saved with the image, keep their position on the gel when you adjust the lanes, and are kept when you switch between automatic and manual lanes.

In **Continuum** fits, each custom peak replaces the closest fragment position in the starting guess.

## The two fit types

Both fit types model each lane's profile as a **background** (a polynomial, see **Polynomial Degree**) plus a sum of **Gaussian peaks**, and adjust them until they match the profile as closely as possible. They differ in where the peaks come from.

### Banded

For gels with **distinct bands**.

- **Starting guess:** the plugin finds the bands in each profile. A band is a local maximum that stands out from the dips on either side by more than the **Peak Tolerance**. Each band gets one Gaussian, with a starting width taken from where the band falls to half its height. Custom peaks are added or replace nearby bands.
- **Fitting:** each peak's position, height and width are adjusted. A peak may move less than 80 % of the distance to its neighbouring band. Its height stays between 1 % and 200 % of the profile's height above the background, and its width between 0.4 and 2 times its starting width.
- **Results:** position, height, width and area of every band. Banded fits don't convert positions to sizes; compare the band positions with the ladder labels in the plots.

### Continuum

For lanes where DNA runs as a **smear**, a continuous range of fragment sizes, such as fragmented, digested or tagmented DNA. It estimates the size distribution of the fragments in each lane.

- **Fragment distribution:** you choose a list of fragment lengths with their relative frequencies, the sizes you expect in the sample (see [Fragment distributions](#fragment-distributions)).
- **Starting guess:** each fragment length gets one Gaussian:
  - **Position** is predicted from the ladder. Between ladder bands, migration distance is interpolated linearly against the logarithm of the fragment's molecular weight, and extrapolated beyond the first and last band. Fragments predicted up to 20 % beyond either end of the lane are included, so that bands at the edges fit properly.
  - **Width** is interpolated from the widths of the ladder bands.
  - **Height** comes from the profile at that position, scaled by the fragment's frequency in the distribution.
- **Fitting:** positions may shift by up to 10 % of the spacing between neighbouring fragments. Widths stay within the range set by **SD Drift**. The areas may depart from the distribution's proportions only as far as **Area Drift** allows.
- **Results:** one row per fragment with its length (**BP**), molecular weight (**MW**) and frequency in the distribution. The **LOG** adds each lane's **Average Fragment Size** ± its standard deviation, and a plot of the fitted size distribution per lane.

Molecular weights are computed for double-stranded DNA as MW = 607.4 × bp + 157.9 Da.

## Fit parameters

The parameters are in the main window. They're saved with each image, and the last values you used become the defaults for new images.

| Parameter | Range (default) | Used in | What it does |
|---|---|---|---|
| **Polynomial Degree** | −1 to 15 (2) | both | Shape of the background under each lane. **−1**: no background. **0**: a constant offset. **1**: a straight slope. **2**: a gentle curve. Higher degrees follow uneven backgrounds more closely but can absorb broad bands, so raise it only when the background clearly isn't smooth. The background is always kept below the lane's lowest intensity, so it can't take over signal. |
| **Max Polynomial Derivative** (grayvalue/px) | 0 to 10 (10) | both | Limit on how steep the background may be, as its average slope in gray values per pixel. If the fitted background is steeper, its slope is reduced until it's within the limit. Lower it when the background should be nearly flat; **0** forces a flat background. |
| **Peak Tolerance** (fraction) | 0.01 to 1 (0.1) | Banded, and the ladder lane in both | How much a band must stand out to be detected, as a fraction of the lane's intensity range: a band's peak must rise above the dips on both sides by more than this fraction. **0.05** means 5 % of the range. Lower it to catch faint bands; raise it if noise is detected as bands. |
| **Area Drift** | 0.001 to 1 (0.1) | Continuum | How far the fragments' fitted areas may depart from the proportions in the fragment distribution. It limits the spread of the ratio *fitted area / expected area* across the fragments, measured as the standard deviation of its logarithm. Small values keep the result close to the distribution's shape; larger values let the fit follow the data. When a fit exceeds the limit, the areas are re-drawn at random around the expected proportions, so Continuum results can differ slightly between runs. |
| **SD Drift** | 1 to 5 (1) | Continuum | How far each fragment's width may depart from the width expected from the ladder bands: between *expected ÷ SD Drift* and *expected × SD Drift*. **1** keeps the widths fixed at the ladder-based values; **2** allows half to double. Raise it if bands in the sample lanes are clearly wider or narrower than in the ladder. |

The **LOG** window lists the parameters used for each fit, next to the results.

## Reading the plots

Each lane's plot shows:

| Curve | Colour | Meaning |
|---|---|---|
| Profile | black | the lane's measured intensity |
| Background | blue | the fitted polynomial background |
| Peaks | red | each fitted Gaussian, drawn on top of the background |
| Fit | orange | background plus all peaks; it should follow the profile |
| Custom peaks | green dots | your custom peaks |

Dashed vertical lines mark the ladder bands, labelled with their sizes. The ladder lane's plot has a light blue background.

- **Zoom:** drag a box over the region to enlarge.
- **Reset the view:** drag a box up and to the left, or right-click and choose **Auto Range › Both Axes**. The reset view shows the whole lane, from the lowest to the highest curve, with enough room above for the ladder labels.
- **More options:** right-click for zoom, copy, save and print.

Plots keep your zoom while you fit, edit custom peaks or toggle edit mode. Changing the lanes rebuilds the profiles, which resets the view. The plot images are saved to the data folder as PNG, PDF and SVG after each fit.

## Reading the results

After a fit, three places show the results:

- **Results Display** (a table window), also saved as `Fit of <image>.xls`, a tab-separated text file that Excel and other spreadsheets open. One row per band or fragment within the lane:

  | Column | Meaning |
  |---|---|
  | Lane, Band | lane number, and band number within the lane from the top |
  | Distance | fitted band position, in pixels from the top of the image |
  | Amplitude | fitted peak height above the background, in gray values |
  | FWHM | fitted band width (full width at half maximum), in pixels |
  | Area | band area, the peak integrated over the lane; it's proportional to the amount of material in the band |
  | Dist. G., Amp. G., FWHM G. | the starting guesses for Distance, Amplitude and FWHM, for comparison |
  | Frequency, BP, MW | Continuum only: the fragment's relative frequency in the distribution, its length in base pairs and its molecular weight. They're shown as "-" for the ladder lane. |

  For Continuum fits, the file also lists each lane's mean fragment size and standard deviation at the end.
- **LOG**, also saved as `<image>_log.html`:
  - the fit type and parameters;
  - each lane's **RMS**, the typical difference between the fit and the profile as a fraction of the lane's intensity range (lower is better);
  - for Continuum fits, each lane's **Average Fragment Size ± standard deviation** in base pairs. It's weighted by each fragment's fitted area times its length. The value in parentheses is the average over the distribution itself, for comparison.
  - the data folder's path at the end.

  For Continuum fits, the LOG window also has one tab per lane with the fitted size distribution.
- **The plots**, see [Reading the plots](#reading-the-plots).

## Where your results are saved

For each image, the plugin creates a folder **next to the image**, named after it:

```
my gel.tif                    ← your image
my gel - Gel Lanes Fit/       ← its results
    saved-state.bak           lanes, ladder, fit settings and custom peaks
    Fit of my gel.xls         results table (tab-separated)
    my gel_log.html           fit summary
    my gel.tif                copy of the image, with the lanes and bands as an overlay
    Lane 1.png / .pdf / .svg  one plot per lane
    …
```

The folder is created when you start the plugin on the image. `saved-state.bak` is updated whenever you change the lanes, the ladder or the custom peaks, run a fit, or close the plugin. The results table, summary, image copy and plots are written each time you run a fit.

- Click **Open Data Folder** in the main window to open this folder in Finder or Explorer. Hovering over the button shows the full path.
- The next time you open the same image, the plugin restores the lanes, ladder, settings and custom peaks from `saved-state.bak` and repeats the last fit.
- Results from versions before 1.0.8 were saved in `gel-lanes-fit/<image name without spaces>/` next to the image. The plugin moves that folder to the new name the first time you open the image, so your saved state carries over.
- If the image has never been saved to disk, the folder is created in the Fiji installation folder instead.

## Ladders and fragment distributions

### Built-in ladders

Bands are listed in the order they appear in the lane, from the top (largest) down.

| Ladder | Bands |
|---|---|
| **Hi-Lo** | 10, 8, 6, 4, 3, 2, 1.55, 1.4, 1 kbp; 750, 500, 400, 300, 200, 100, 50 bp |
| **100bp** | 1.5, 1.2, 1 kbp; 900, 800, 700, 600, 500, 400, 300, 200, 100 bp |
| **Quick-Load** | 48.5, 20, 15, 10, 8, 6, 5, 4, 3, 2, 1.5, 1 kbp; 500 bp |
| **Tapestation** | 1.5, 1 kbp; 700, 500, 400, 300, 200, 100, 50, 25 bp |

The 100bp ladder's top band is 1,517 bp, shown as 1.5 kbp.

### Custom ladder files

Choose **Custom Ladder** in the ladder type list to load your own ladder from a plain text file:

- one band size per line, as a whole number of base pairs;
- in the order the bands appear in the lane, from the top (largest) down;
- empty lines are ignored.

For example:

```
1500
1000
700
500
400
300
```

If the file can't be read or has no sizes, the plugin says so and keeps the previous ladder. The ladder you load is saved with the image, so you don't need the file again for that image.

### Fragment distributions

Continuum fits need a **fragment distribution**: the fragment lengths expected in the sample and how common each one is. Choose it in the list below the ladder:

| Distribution | Description |
|---|---|
| **AciI-Lambda**, **AciI-Lambda4**, **AciI-Lambda3**, **AciI-Lambda2** | four approximations of the same fragment distribution, an AciI digest of lambda DNA (2 to 1,086 bp), made with decreasing numbers of fragment lengths: 206, 155, 138 and 103. Each length becomes one peak in the fit, so fewer lengths mean fewer peaks to fit. Fitting the same lanes with each variant shows whether a coarser approximation, with fewer fragment lengths, gives a poorer fit or not; compare the lanes' RMS in the LOG. |
| **Uniform** | every length between a **Lower** and an **Upper** size, in steps of **Every** bp, all equally frequent; the plugin asks for the three values when you choose it |
| **Ladder** | a fixed list of 22 sizes from 50 bp to 10 kbp (the Hi-Lo and 100bp band sizes), each once; it doesn't follow the ladder you selected |

The distributions are tab-separated text files: a header line, then one line per fragment with its length in base pairs and its number of copies. Further columns are ignored. They're built into the plugin, so adding your own requires rebuilding it.

## Troubleshooting

### "The number of peaks detected does not match the ladder range"

The fit found a different number of bands in the ladder lane than the ladder range contains, so it stopped before fitting the other lanes. The message shows both numbers. Check the ladder lane's plot, then:

- **Too few bands detected:** lower **Peak Tolerance**, or add the missing bands as custom peaks.
- **Too many detected:** raise **Peak Tolerance**, or check the lane covers only the ladder.
- **Range doesn't match what you see:** pick the ladder type again and choose the first and last bands that are actually visible.
- Make sure **Select Ladder Lane** points at the lane with the ladder.

Then click **Fit** again.

### The plugin is not in the Plugins menu

- Check that the jar is in the `jars` folder of the Fiji you're starting. It's easy to have more than one Fiji installed.
- Quit Fiji completely and start it again. Fiji only looks for new plugins when it starts.
- Look in **Window › Console** for errors mentioning `gellanesfit`.

### An old version still runs after updating

Fiji found more than one copy of the plugin. Check the version in the main window's title bar, then remove every `gel-lanes-fit-*.jar` and `GaussianFit-*.jar` from both `jars` and `plugins`, except the one you want, and restart Fiji.

### Fiji reports `UnsupportedClassVersionError`

A message like this means your Fiji runs on a Java version that's too old:

```
gellanesfit/GelLanesFit has been compiled by a more recent version of the Java Runtime (class file version 55.0), this version of the Java Runtime only recognizes class file versions up to 52.0
```

Install a current Fiji from [fiji.sc](https://fiji.sc). Updating a very old installation may not replace its Java.

### "There are no images open"

The plugin works on the image that is open and selected. Open your gel image first, click its window, then start Gel Lanes Fit.

### Messages in the Console that you can ignore

These appear when Fiji starts and don't affect Gel Lanes Fit:

```
SLF4J(W): No SLF4J providers were found.
SLF4J(W): Defaulting to no-operation (NOP) logger implementation
```

A library inside Fiji has no logging back-end installed, so its own log messages are dropped. The plugin's messages still appear.

```
[ERROR] Cannot create plugin: org.scijava.plugins.scripting.javascript.JavaScriptScriptLanguage
```

Fiji's JavaScript scripting support needs a component that newer Java versions no longer include. It only matters if you run JavaScript macros or scripts; Gel Lanes Fit doesn't use it.

### I can't find the results folder

- Click **Open Data Folder** in the main window. It opens the folder directly.
- The **Console** shows the path on a line starting `Data folder:`.
- On macOS, Finder may list the folder away from the image if the window is grouped or sorted by date, for example under **Today**. Choose **View › Use Groups** to turn grouping off, or sort by name.
- Finder sometimes doesn't show a folder that was just created or renamed by another program until it refreshes. Go up a folder and back, or press ⇧⌘G and paste the path. If that doesn't help, relaunch Finder: hold ⌥ (Option), right-click the Finder icon in the Dock and choose **Relaunch**.

### "Cannot create the data folder" in the Console

The plugin can't write next to the image, for example because the image is on a read-only drive or network share, or in a folder you don't have write access to. Copy the image to a folder you can write to, such as your Documents folder, and open it from there.
