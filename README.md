# Gel Lanes Fit

Gel Lanes Fit is a [Fiji](https://fiji.sc) (ImageJ) plugin for quantitative analysis of gel electrophoresis images.

It takes an intensity profile along each lane and fits it with Gaussian peaks on a polynomial background, with a size ladder as reference:

- **Lanes:** placed automatically, equally sized and spaced, or drawn one by one for bent or uneven lanes.
- **Two fit types:** **Banded**, for gels with distinct bands, reports each band's position, height, width and area. **Continuum**, for smears such as fragmented or tagmented DNA, estimates each lane's fragment size distribution and average fragment size against the ladder.
- **Ladders:** built-in Hi-Lo, 100bp, Quick-Load and Tapestation ladders, or your own from a text file.
- **Custom peaks:** add or remove bands by clicking on a lane's plot, then fit again.
- **Results:** a results table, a fit summary, and each lane's plot as PNG, PDF and SVG, saved in a folder next to the image.
- **Saved analysis:** lanes, ladder, settings and custom peaks are saved with each image, and the last fit is repeated when you open the image again.

## Installing

1. Install or update [Fiji](https://fiji.sc).
2. Download `gel-lanes-fit-<version>.jar` from the [Releases page](https://github.com/rickud/gel-lanes-fit/releases).
3. Quit Fiji, copy the jar into Fiji's `jars` folder (removing any older `gel-lanes-fit-*.jar`), and start Fiji again.

The plugin is in **Plugins › Gel Tools › Gel Lanes Fit**. The [user guide](docs/user-guide.md#installing) has the details and a troubleshooting section.

## Quick start

1. Open a gel image in Fiji and choose **Plugins › Gel Tools › Gel Lanes Fit**.
2. Place the lanes: **Automatic Rectangle Selection** with the sliders, or **Manual Rectangle Selection** to draw them.
3. Pick the ladder lane, the ladder type and the range of ladder bands visible in it.
4. Choose **Banded** or **Continuum** (and a fragment distribution for Continuum), then click **Fit**.
5. Check the fit in the Profiles window, add custom peaks where needed and fit again. Click **Open Data Folder** to see the saved results.

## Documentation

The **[user guide](docs/user-guide.md)** covers:

- [installing and troubleshooting](docs/user-guide.md#installing);
- [using the plugin, step by step](docs/user-guide.md#using-the-plugin);
- [the two fit types](docs/user-guide.md#the-two-fit-types) and [what each fit parameter does](docs/user-guide.md#fit-parameters);
- [reading the plots](docs/user-guide.md#reading-the-plots) and [the results](docs/user-guide.md#reading-the-results);
- [where results are saved](docs/user-guide.md#where-your-results-are-saved);
- [ladders, custom ladder files and fragment distributions](docs/user-guide.md#ladders-and-fragment-distributions).

The same guide is published in the [wiki](https://github.com/rickud/gel-lanes-fit/wiki).

Questions and bug reports are welcome in [Issues](https://github.com/rickud/gel-lanes-fit/issues).

## Development

Build with Maven and install into a local Fiji:

```
mvn -Dscijava.app.directory=/path/to/Fiji install
```

To avoid passing the path every time, set `scijava.app.directory` in a profile in `~/.m2/settings.xml`.

`GelLanesFit.main` starts ImageJ with a sample image, so you can run the plugin from an IDE. On Java 17+ add these JVM options to the run configuration (the second one only matters on macOS):

```
--add-opens=java.base/java.lang=ALL-UNNAMED --add-exports=java.desktop/com.apple.eawt=ALL-UNNAMED
```

Without the first option, ImageJ fails at startup with `No _hooks field found in ij.IJ`.

## Authors

Developed at [The University of Texas at Dallas](https://www.utdallas.edu) by Rick Ziraldo ([@rickud](https://github.com/rickud)), with contributions from Massa Shoura and Stephen Levene.

## License

Gel Lanes Fit is free software under the [GNU Affero General Public License v3.0](LICENSE) (AGPL-3.0).
