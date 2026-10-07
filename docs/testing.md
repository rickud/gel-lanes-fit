# Testing Gel Lanes Fit

Every change is tested in two ways before it's merged into `master`: the automated tests, and the manual checklist below, which covers what only a person using the plugin in Fiji can check.

## Automated tests

The tests are in `src/test/java` and use JUnit 4. They run without opening any windows.

- **In Eclipse:** right-click `src/test/java` and choose **Run As › JUnit Test**.
- **With Maven:** `mvn test`. They also run as part of `mvn install`, so **Build + Fiji** stops if a test fails.

| Test | What it checks |
|---|---|
| `LadderTest` | band names and molecular weights of the built-in ladders, ranges, serialization |
| `BandDetectionTest` | band detection and starting widths on synthetic profiles, and the Peak Tolerance |
| `BandedFitTest` | Banded fits of synthetic lanes, custom peaks, the results table |
| `ContinuumFitTest` | the bundled distributions, a Continuum fit against a ladder, reproducible results, the average fragment size |
| `SavedStateTest` | `saved-state.bak` files from current and older versions still load |
| `SampleImageFitTest` | lane profiles and a Banded fit of four lanes of the sample image `Long_5s.tif` |

### Reference files

The fit tests compare their results with reference files in `src/test/resources/gellanesfit/reference`, so that any change in the results shows up, not just large errors. When a change is meant to change the results, for example an improvement to the fitting:

1. Run the tests once with the system property `-Dglf.record=true`. In Eclipse, add it under **Run Configurations › Arguments › VM arguments**; with Maven, run `mvn test -DargLine=-Dglf.record=true`. The reference files are rewritten.
2. Review the changes to the reference files with `git diff` before committing them, and say in the commit message why the results changed.

The fixtures in `src/test/resources/gellanesfit/saved-state` are saved-state files written by earlier versions. Don't regenerate them: they make sure users' existing saved analyses keep loading.

## Manual checklist

Run **Build + Fiji** with the sample image `src/main/resources/sample-images/tagment-test/gel-camera-1/Long_5s.tif`: the Fiji launch configuration's arguments should open that image and run the plugin. Watch **Window › Console** for errors throughout. Delete or rename the image's data folder (`Long_5s - Gel Lanes Fit`, next to the image) first, so you start without a saved state.

### Lanes

- [ ] **Automatic lanes:** changing **Number of Lanes** and each slider moves the rectangles on the image and updates the plots.
- [ ] **Manual lanes:** switch to **Manual Rectangle Selection**. Draw a lane, move it, resize it; the plots follow. Click inside a lane and confirm: it's deleted and the lanes are renumbered from left to right.
- [ ] **Switching modes:** switching back to Automatic and then to Manual again restores the manual lanes.
- [ ] **Lane under the mouse:** moving the mouse over a lane highlights its plot and draws a line at the same position.

### Ladder and fit

- [ ] **Ladder:** pick the ladder lane (its plot turns light blue), a ladder type, and a range in the Ladder Range dialog, also with the bands picked in reverse order.
- [ ] **Banded fit:** click **Fit**. The ladder bands are labelled in every plot, the curves appear, and the Results Display, LOG and data folder files are created.
- [ ] **Band count mismatch:** narrow the ladder range so it doesn't match and fit again. A warning shows both counts and nothing else is fitted.
- [ ] **Continuum fit:** choose **Continuum**, then each distribution in turn (the AciI-Lambda ones, Ladder, and Uniform with its range dialog), and fit. The results table has Frequency, BP and MW columns, and the LOG shows Average Fragment Size and a tab per lane.
- [ ] **Uniform dialog:** typing Lower, Upper and Every updates the number of fragment lengths under the fields, with a note above 100.
- [ ] **Slow fit warning:** a Continuum fit with a Uniform distribution of a few hundred lengths asks before fitting; **Cancel** stops after the ladder lane, **Fit Anyway** fits.
- [ ] **Fit warning:** with a fit shown, changing the lanes or fitting again shows the warning; **Don't show this again** stops it.

### Custom peaks and display

- [ ] **Edit Custom Peaks:** clicking in a plot asks for the FWHM and adds a green dot; clicking near a dot removes it; clicking in the plot margin does nothing.
- [ ] **Refit with custom peaks:** **Fit** uses the custom peaks. **Reset Custom Peaks** removes them from the chosen lanes.
- [ ] **Show Bands:** marks the bands on the image (magenta fitted, blue guesses, green custom) and keeps its state.
- [ ] **Plots:** zoom with a dragged box; toggling edit mode, fitting and editing peaks keep the zoom. Dragging up and to the left resets the view, with the ladder labels clear of the curves.
- [ ] **Open Data Folder:** opens the image's data folder.

### Saved state

- [ ] **Restore:** close the plugin and start it again on the same image. The lanes, ladder, settings, Show Bands, and custom peaks come back, and the last fit runs on its own.
- [ ] **Interrupted fit:** start a slow fit (for example Continuum with a Uniform distribution of a few hundred lengths) and quit Fiji while it runs. On the next start the plugin restores the lanes and settings, doesn't repeat the fit, and explains why in a message.
- [ ] **Old data folder:** an old-style `gel-lanes-fit/<name>/` folder next to an image is moved to `<image name> - Gel Lanes Fit` on first start.

Note anything unexpected in the pull request or merge commit, with the Console output if there was an error.
