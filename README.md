# Gel Lanes Fit

This plugin for ImageJ provides some functionality to set up and perform quatitative analysys of gel electrophoresis images.
The plugin can draw same-size or custom ROIs (regions of interest).
The areas selected by the ROIs are used to compute respective line profiles for each gel lane. The profiles obtained by averaging pixel intensities along horizontal lines in the ROI.
Gaussian peaks can be fitted to the line profiles. A polynomial function can be used to represent the background signal in each profile.

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
