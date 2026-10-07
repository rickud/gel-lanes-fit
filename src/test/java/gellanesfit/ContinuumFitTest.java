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
		final Ladder ladder = Ladder.create(Ladder.HILO);
		ladder.setRange(LADDER_RANGE);
		final RealVector mw = ladder.getMolecularWeights();

		final double[][] ladderBands = new double[mw.getDimension()][];
		for (int i = 0; i < ladderBands.length; i++)
			ladderBands[i] = band(100, position(mw.getEntry(i)), 4);
		final double[] ladderLane = TestProfiles.values(50, 450, 10, 0,
			ladderBands);
		// A smear centred around 300 bp
		final double[] smear = TestProfiles.values(50, 450, 10, 0, band(120,
			position(300 * 607.4 + 157.9), 40));

		final Fitter fitter = BandedFitTest.bandedFitter(context, TestProfiles
			.profile(1, 50, ladderLane), TestProfiles.profile(2, 50, smear));
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
