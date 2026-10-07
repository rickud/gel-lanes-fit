package gellanesfit;

import static gellanesfit.TestProfiles.band;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;
import org.scijava.Context;
import org.scijava.app.StatusService;
import org.scijava.display.DisplayService;
import org.scijava.log.LogService;

/** Banded fits through the Fitter, as the Fit button runs them. */
public class BandedFitTest {

	private static Context context;

	@Rule
	public TemporaryFolder folder = new TemporaryFolder();

	@BeforeClass
	public static void createContext() {
		context = new Context(LogService.class, StatusService.class,
			DisplayService.class);
	}

	@AfterClass
	public static void disposeContext() {
		context.dispose();
	}

	/** A Fitter with the plugin's default parameters, in Banded mode */
	static Fitter bandedFitter(final Context ctx, final DataSeries... lanes) {
		final Fitter fitter = new Fitter(ctx, "test");
		fitter.setDegBG(2);
		fitter.setPolyDerivative(10);
		fitter.setTolPK(0.1);
		fitter.setAreaDrift(0.1);
		fitter.setSDDrift(1.0);
		fitter.setFitMode(Fitter.bandMode);
		fitter.setInputData(new ArrayList<>(Arrays.asList(lanes)));
		return fitter;
	}

	private static final double[][] BANDS = { band(100, 180, 6), band(60, 300,
		8), band(80, 420, 5) };

	@Test
	public void recoversKnownBandsOnASlopedBackground() {
		final double[] y = TestProfiles.values(100, 500, 20, 0.02, BANDS);
		final Fitter fitter = bandedFitter(context, TestProfiles.profile(1, 100,
			y));
		fitter.doFit(1);

		final List<Peak> peaks = fitter.getFittedPeaks(1);
		assertEquals(3, peaks.size());
		for (int i = 0; i < BANDS.length; i++) {
			final Peak p = peaks.get(i);
			assertEquals("position of band " + i, BANDS[i][1], p.getMean(), 1.0);
			// The narrowest band comes out about 19 % too wide (5.95 for 5)
			assertEquals("width of band " + i, BANDS[i][2], p.getSigma(),
				0.25 * BANDS[i][2]);
			assertEquals("height of band " + i, BANDS[i][0], p.getNorm(),
				0.15 * BANDS[i][0]);
		}
		Reference.check("banded-sloped-background", Reference.peaks(peaks));
		Reference.check("banded-sloped-background-guess", Reference.peaks(fitter
			.getGuessPeaks(1)));
	}

	@Test
	public void fitsWithoutABackground() {
		// Polynomial Degree -1: no background (used to throw in doGuess)
		final double[] y = TestProfiles.values(100, 500, 0, 0, BANDS);
		final Fitter fitter = bandedFitter(context, TestProfiles.profile(1, 100,
			y));
		fitter.setDegBG(-1);
		final List<DataSeries> curves = fitter.doFit(1);

		final List<Peak> peaks = fitter.getFittedPeaks(1);
		assertEquals(3, peaks.size());
		for (int i = 0; i < BANDS.length; i++)
			assertEquals("position of band " + i, BANDS[i][1], peaks.get(i)
				.getMean(), 1.0);
		for (final DataSeries d : curves)
			if (d.getType() == DataSeries.BACKGROUND) assertEquals(
				"the background is zero", 0, d.getY().getLInfNorm(), 0);
	}

	@Test
	public void customPeakAddsABandTheDetectionMissed() {
		// The middle band is too weak for the 0.1 tolerance
		final double[] y = TestProfiles.values(100, 500, 20, 0, band(100, 180, 6),
			band(8, 300, 6), band(100, 420, 6));
		Fitter fitter = bandedFitter(context, TestProfiles.profile(1, 100, y));
		fitter.doFit(1);
		assertEquals(2, fitter.getFittedPeaks(1).size());

		fitter = bandedFitter(context, TestProfiles.profile(1, 100, y));
		fitter.addCustomPeak(new Peak(1, 8, 300, 6));
		fitter.doFit(1);
		final List<Peak> peaks = fitter.getFittedPeaks(1);
		assertEquals(3, peaks.size());
		assertEquals(300, peaks.get(1).getMean(), 2.0);
		Reference.check("banded-custom-peak", Reference.peaks(peaks));
	}

	@Test
	public void resultsTableListsFittedLanesAndSkipsTheOthers()
		throws Exception
	{
		final double[] y = TestProfiles.values(100, 500, 20, 0.02, BANDS);
		final Fitter fitter = bandedFitter(context, TestProfiles.profile(1, 100,
			y), TestProfiles.profile(2, 100, y));
		fitter.doFit(1); // like the ladder lane, fitted before the others

		final String dir = folder.getRoot().getAbsolutePath() + File.separator;
		fitter.updateResultsTable(dir); // used to fail on lane 2 (issue #7)

		final List<String> rows = Files.readAllLines(new File(dir +
			"Fit of test.xls").toPath());
		assertTrue(rows.get(0).startsWith("Lane\tBand\tDistance"));
		assertEquals("header and 3 bands of lane 1", 4, rows.size());
		assertTrue(rows.get(1).startsWith("1\t1\t"));
		Reference.checkText("banded-results-file", rows);
	}
}
