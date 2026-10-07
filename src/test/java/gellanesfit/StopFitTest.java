package gellanesfit;

import static gellanesfit.TestProfiles.band;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.Arrays;
import java.util.concurrent.CancellationException;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.apache.commons.math3.linear.RealVector;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import org.scijava.Context;
import org.scijava.app.StatusService;
import org.scijava.display.DisplayService;
import org.scijava.log.LogService;

/** A slow fit running in another thread stops when asked to. */
public class StopFitTest {

	private static Context context;

	@BeforeClass
	public static void createContext() {
		context = new Context(LogService.class, StatusService.class,
			DisplayService.class);
	}

	@AfterClass
	public static void disposeContext() {
		context.dispose();
	}

	@Test(timeout = 60000)
	public void stopsASlowContinuumFitPartWay() throws Exception {
		final Ladder ladder = Ladder.create(Ladder.HILO);
		ladder.setRange(new int[] { 8, 14 }); // 1 kbp to 100 bp
		final RealVector mw = ladder.getMolecularWeights();
		final double[][] bands = new double[mw.getDimension()][];
		for (int i = 0; i < bands.length; i++)
			bands[i] = band(100, 100 + 250 * (Math.log10(Ladder.molecularWeight(
				1000)) - Math.log10(mw.getEntry(i))), 4);
		final double[] smear = TestProfiles.values(50, 450, 10, 0, band(120, 330,
			60));
		final Fitter fitter = BandedFitTest.bandedFitter(context, TestProfiles
			.profile(1, 50, TestProfiles.values(50, 450, 10, 0, bands)),
			TestProfiles.profile(2, 50, smear), TestProfiles.profile(3, 50, smear));
		fitter.setReferenceLane(1);
		fitter.doFit(1);
		fitter.setLadder(mw);
		fitter.setFitMode(Fitter.continuumMode);
		// Many fragment lengths: a fit that takes a long time
		fitter.setFragmentDistribution(FragmentDistribution.uniform(50, 1500, 3));
		assertTrue(fitter.fragmentsToFit(2) > 300);

		final ExecutorService executor = Executors.newSingleThreadExecutor();
		try {
			final Future<?> fit = executor.submit(() -> fitter.doFit(Arrays.asList(
				2, 3)));
			Thread.sleep(500);
			final long asked = System.nanoTime();
			fitter.requestStop();
			try {
				fit.get(20, TimeUnit.SECONDS);
				fail("The fit should have been stopped");
			}
			catch (final ExecutionException e) {
				assertTrue(e.getCause() instanceof CancellationException);
			}
			System.out.printf("Stopped %.1f s after being asked%n", (System
				.nanoTime() - asked) / 1e9);
		}
		finally {
			executor.shutdownNow();
		}
		assertEquals("lane 3 wasn't started", 0, fitter.getFittedPeaks(3).size());

		// The stop doesn't carry over to the next fit
		fitter.setFragmentDistribution(FragmentDistribution.uniform(200, 500,
			50));
		fitter.doFit(Arrays.asList(3));
		assertTrue(fitter.getFittedPeaks(3).size() > 0);
	}
}
