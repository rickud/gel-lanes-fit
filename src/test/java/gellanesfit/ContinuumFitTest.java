package gellanesfit;

import static gellanesfit.TestProfiles.band;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;
import org.scijava.Context;
import org.scijava.app.StatusService;
import org.scijava.display.DisplayService;
import org.scijava.log.LogService;

/** Continuum fits: a ladder lane fitted first, then a smear. */
public class ContinuumFitTest {

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

	/** The smear lane's profile, set by fitContinuum */
	private double[] smear;

	/** Hi-Lo bands 1 kbp to 100 bp: indexes 8 to 14 */
	private static final int[] LADDER_RANGE = { 8, 14 };

	/** Migration distance of a band of the given molecular weight */
	private static double position(final double mw) {
		return 100 + 250 * (Math.log10(1000 * 607.4 + 157.9) - Math.log10(mw));
	}

	/** A bundled fragment distribution, read as the plugin reads it */
	static double[][] distribution(final String name) throws Exception {
		final InputStream in = ContinuumFitTest.class.getClassLoader()
			.getResourceAsStream(FragmentDistribution.FOLDER + name + ".txt");
		assertNotNull(name + " is packaged", in);
		try (BufferedReader r = new BufferedReader(new InputStreamReader(in))) {
			return FragmentDistribution.read(r, w -> fail(name + ": " + w));
		}
	}

	private Fitter fitContinuum() throws Exception {
		return fitContinuum(2);
	}

	private Fitter fitContinuum(final int degBG) throws Exception {
		final Ladder ladder = Ladder.create(Ladder.HILO);
		ladder.setRange(LADDER_RANGE);
		final RealVector mw = ladder.getMolecularWeights();

		final double[][] ladderBands = new double[mw.getDimension()][];
		for (int i = 0; i < ladderBands.length; i++)
			ladderBands[i] = band(100, position(mw.getEntry(i)), 4);
		final double[] ladderLane = TestProfiles.values(50, 450, 10, 0,
			ladderBands);
		// A smear centred around 300 bp
		smear = TestProfiles.values(50, 450, 10, 0, band(120,
			position(300 * 607.4 + 157.9), 40));

		final Fitter fitter = BandedFitTest.bandedFitter(context, TestProfiles
			.profile(1, 50, ladderLane), TestProfiles.profile(2, 50, smear));
		fitter.setDegBG(degBG);
		fitter.setReferenceLane(1);
		fitter.doFit(1);
		assertEquals("all ladder bands found", mw.getDimension(), fitter
			.getFittedPeaks(1).size());

		fitter.setLadder(mw);
		fitter.setFitMode(Fitter.continuumMode);
		fitter.setFragmentDistribution(distribution("AciI-Lambda2"));
		fitter.doFit(2);
		return fitter;
	}

	@Test
	public void bundledDistributionsLoad() throws Exception {
		final String[] names = { "AciI-Lambda", "AciI-Lambda2", "AciI-Lambda3",
			"AciI-Lambda4", "Ladder" };
		final int[] lengths = { 206, 103, 138, 155, 22 };
		for (int i = 0; i < names.length; i++) {
			final double[][] d = distribution(names[i]);
			assertEquals(names[i], lengths[i], d.length);
			double sum = 0;
			for (final double[] row : d)
				sum += row[0];
			assertEquals(names[i] + " frequencies add up to 1", 1, sum, 1e-9);
		}
	}

	@Test
	public void fitsTheSmearWithOnePeakPerFragmentLength() throws Exception {
		final Fitter fitter = fitContinuum();
		final List<Peak> peaks = fitter.getFittedPeaks(2);
		assertTrue("fragments in the lane", peaks.size() > 10);
		Reference.check("continuum-smear", Reference.peaks(peaks));
	}

	@Test
	public void fitsWithoutABackground() throws Exception {
		// Polynomial Degree -1: no background (used to throw in doGuess)
		final Fitter fitter = fitContinuum(-1);
		assertTrue("fragments in the lane", fitter.getFittedPeaks(2).size() > 10);
	}

	@Test
	public void usesAFlatBackgroundBelowTheLanesMinimum() throws Exception {
		final Fitter fitter = fitContinuum();
		final List<DataSeries> curves = fitter.doFit(2);
		double laneMin = Double.MAX_VALUE;
		for (final double v : smear)
			laneMin = Math.min(laneMin, v);
		for (final DataSeries d : curves)
			if (d.getType() == DataSeries.BACKGROUND) {
				final RealVector bg = d.getY();
				assertEquals("98 % of the minimum", 0.98 * laneMin, bg.getMinValue(),
					1e-9);
				assertEquals("flat", bg.getMinValue(), bg.getMaxValue(), 0);
			}
	}

	@Test
	public void uniformCoversTheLaneFromBottomToTop() throws Exception {
		final Fitter fitter = fitContinuum();
		// The synthetic ladder: 1 kbp at x = 100, 250 px per decade of MW
		final double[] range = fitter.sizeRange(new ArrayRealVector(new double[] {
			50, 450 }));
		assertEquals("top of the lane", basePairsAt(50), range[1], 0.01 *
			range[1]);
		assertEquals("bottom of the lane", basePairsAt(450), range[0], 0.01 *
			range[0]);

		fitter.setUniformDistribution(25, 0);
		fitter.resetFit(2); // as before every fit in the plugin
		fitter.doFit(2);
		final List<Peak> guess = fitter.getGuessPeaks(2);
		assertEquals(25, guess.size());
		assertEquals("first peak at the top", 50, guess.get(0).getMean(), 1);
		assertEquals("last peak at the bottom", 450, guess.get(24).getMean(), 1);

		fitter.setUniformDistribution(0, 100);
		assertEquals("at most 100 bp apart", FragmentDistribution.uniformCount(
			range[0], range[1], 100), fitter.fragmentsToFit(2));
	}

	@Test
	public void widthsStayPositiveBeyondTheLadder() throws Exception {
		// Bands widening steeply down the gel: the line through their widths
		// is below zero at the top of the lane
		final Ladder ladder = Ladder.create(Ladder.HILO);
		ladder.setRange(LADDER_RANGE);
		final RealVector mw = ladder.getMolecularWeights();
		final double[][] ladderBands = new double[mw.getDimension()][];
		for (int i = 0; i < ladderBands.length; i++) {
			final double x = position(mw.getEntry(i));
			ladderBands[i] = band(100, x, 1.5 + 0.025 * (x - 100));
		}
		final Fitter fitter = BandedFitTest.bandedFitter(context, TestProfiles
			.profile(1, 50, TestProfiles.values(50, 450, 10, 0, ladderBands)),
			TestProfiles.profile(2, 50, TestProfiles.values(50, 450, 10, 0, band(
				120, 250, 60))));
		fitter.setReferenceLane(1);
		fitter.doFit(1);
		assertEquals(mw.getDimension(), fitter.getFittedPeaks(1).size());
		fitter.setLadder(mw);
		fitter.setFitMode(Fitter.continuumMode);
		fitter.setUniformDistribution(25, 0);
		fitter.doFit(2); // used to throw: a negative width at the top
		for (final Peak p : fitter.getGuessPeaks(2))
			assertTrue("width at " + p.getMean(), p.getSigma() > 1);
	}

	@Test
	public void logSpacedPlacesFragmentsEvenlyAlongTheLane() throws Exception {
		final Fitter fitter = fitContinuum();
		fitter.setLogDistribution(41);
		fitter.resetFit(2);
		fitter.doFit(2);
		final List<Peak> guess = fitter.getGuessPeaks(2);
		assertEquals(41, guess.size());
		assertEquals("first at the top", 50, guess.get(0).getMean(), 1);
		assertEquals("last at the bottom", 450, guess.get(40).getMean(), 1);
		// Log of the size runs with the position: about 10 px apart throughout
		for (int i = 1; i < guess.size(); i++)
			assertEquals("spacing before fragment " + i, 10, guess.get(i)
				.getMean() - guess.get(i - 1).getMean(), 1);
	}

	@Test
	public void guessHeightsAreTheFitsStartingPoint() throws Exception {
		// Amp. G. shows where the fit starts: with closely spaced peaks, the
		// guessed peaks add up to the profile above the background
		final Fitter fitter = fitContinuum();
		fitter.setLogDistribution(81);
		fitter.resetFit(2);
		fitter.doFit(2);
		double laneMin = Double.MAX_VALUE;
		for (final double v : smear)
			laneMin = Math.min(laneMin, v);
		for (final int x : new int[] { 200, 250, 300 }) {
			double sum = 0;
			for (final Peak p : fitter.getGuessPeaks(2))
				sum += new org.apache.commons.math3.analysis.function.Gaussian(p
					.getNorm(), p.getMean(), p.getSigma()).value(x);
			final double above = smear[x - 50] - 0.98 * laneMin;
			assertEquals("guess at " + x, above, sum, 0.15 * above);
		}
	}

	/** The size of a fragment at a position of the synthetic ladder, in bp */
	private static double basePairsAt(final double x) {
		return Ladder.basePairs(Math.pow(10, Math.log10(1000 * 607.4 + 157.9) -
			(x - 100) / 250));
	}

	@Test
	public void resultsLabelEachFragmentWithItsOwnLength() throws Exception {
		// Fragments beyond the lane's ends are fitted too; the rows must still
		// pair each fitted peak with its own length
		final Fitter fitter = fitContinuum();
		fitter.setUniformDistribution(25, 0);
		fitter.resetFit(2); // as before every fit in the plugin
		fitter.doFit(2);
		final String dir = folder.getRoot().getAbsolutePath() + File.separator;
		fitter.updateResultsTable(dir);
		final List<String> rows = Files.readAllLines(new File(dir +
			"Fit of test.xls").toPath());
		int checked = 0;
		for (final String row : rows) {
			final String[] cols = row.split("\t");
			if (cols.length < 11 || !cols[10].trim().matches("[0-9.]+")) continue;
			final double guessPosition = Double.parseDouble(cols[3]);
			final double bp = Double.parseDouble(cols[10]);
			assertEquals("length of the fragment guessed at " + guessPosition,
				basePairsAt(guessPosition), bp, 0.02 * bp + 1);
			checked++;
		}
		assertTrue("rows checked", checked > 10);
	}

	@Test
	public void predictsTheNumberOfFragmentsToFit() throws Exception {
		final Fitter fitter = fitContinuum();
		assertEquals(fitter.getFittedPeaks(2).size(), fitter.fragmentsToFit(2));
	}

	@Test
	public void repeatedFitsGiveTheSameResult() throws Exception {
		assertEquals(Reference.peaks(fitContinuum().getFittedPeaks(2)), Reference
			.peaks(fitContinuum().getFittedPeaks(2)));
	}

	@Test
	public void reportsTheAverageFragmentSize() throws Exception {
		final Fitter fitter = fitContinuum();
		final String dir = folder.getRoot().getAbsolutePath() + File.separator;
		fitter.updateResultsTable(dir);

		final List<String> rows = Files.readAllLines(new File(dir +
			"Fit of test.xls").toPath());
		assertTrue(rows.get(0).endsWith("Frequency\tBP\tMW\t"));
		assertTrue(fitter.getSummary().contains("Average Fragment Size"));
		final List<String> stats = new ArrayList<>();
		for (final String row : rows)
			if (row.startsWith("Lane ")) stats.add(row.replace("Lane ", "")
				.replace("\t", " "));
		Reference.check("continuum-average-size", stats);
		Reference.checkText("continuum-results-file", rows);
	}
}
