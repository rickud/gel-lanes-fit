package gellanesfit;

import static gellanesfit.TestProfiles.band;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

/** The band fitter on its own, on profiles made from known bands. */
public class BandFitTest {

	@Test
	public void envelopeFollowsASlopeUnderBands() {
		// A straight slope with narrow bands: the envelope is the slope
		final double[] y = TestProfiles.values(0, 400, 50, 0.2, band(100, 100, 4),
			band(150, 250, 5));
		final double[] env = BandFit.lowerEnvelope(y, 40);
		// Wherever the profile is on the slope (away from the bands' tails)
		int checked = 0;
		for (int j = 30; j < 370; j++) {
			if (y[j] - (50 + 0.2 * j) > 0.01) continue;
			assertEquals("envelope at " + j, 50 + 0.2 * j, env[j], 0.2);
			checked++;
		}
		assertTrue(checked > 250);
	}

	@Test
	public void envelopeNeverRisesAboveTheProfile() {
		final double[] y = TestProfiles.values(0, 200, 10, 0, band(-8, 100, 3));
		final double[] env = BandFit.lowerEnvelope(y, 40);
		for (int j = 0; j < y.length; j++)
			assertTrue("at " + j, env[j] <= y[j]);
	}

	/** Fits three known bands on a curved background */
	private static double[] fit(final int deg, final double maxSlope,
		final double[] background)
	{
		final double[][] bands = { band(100, 120, 5), band(60, 200, 6), band(80,
			300, 4) };
		final double[] x = new double[401];
		final double[] y = new double[x.length];
		for (int j = 0; j < x.length; j++) {
			x[j] = j;
			y[j] = background[j];
			for (final double[] b : bands)
				y[j] += new org.apache.commons.math3.analysis.function.Gaussian(b[0],
					b[1], b[2]).value(j);
		}
		// Starting guesses a little off
		final BandFit.Result r = new BandFit(x, y, deg, new double[] { 122, 197,
			303 }, new double[] { 6, 5, 5 }, maxSlope, () -> false).fit();
		return r.parameters;
	}

	@Test
	public void recoversBandsOnACurvedBackground() {
		final double[] bg = new double[401];
		for (int j = 0; j < bg.length; j++)
			bg[j] = 40 + 0.1 * j - 0.0002 * j * j;
		final double[] p = fit(2, 10, bg);
		final double[][] expected = { { 100, 120, 5 }, { 60, 200, 6 }, { 80, 300,
			4 } };
		for (int i = 0; i < 3; i++) {
			assertEquals("height " + i, expected[i][0], p[4 + 3 * i], 1);
			assertEquals("position " + i, expected[i][1], p[5 + 3 * i], 0.1);
			assertEquals("width " + i, expected[i][2], p[6 + 3 * i], 0.05);
		}
		// The background, in x: 40 + 0.1 x - 0.0002 x^2
		assertEquals(40, p[1], 0.5);
		assertEquals(0.1, p[2], 0.005);
		assertEquals(-0.0002, p[3], 0.00002);
	}

	@Test
	public void maxPolynomialDerivativeLimitsTheSlope() {
		final double[] bg = new double[401];
		for (int j = 0; j < bg.length; j++)
			bg[j] = 40 + 0.1 * j;
		// Average slope limited to 0.02 per px, against 0.1 in the profile
		final double[] p = fit(1, 0.02, bg);
		assertEquals("slope held at the limit", 0.02, p[2], 1e-6);
	}

	@Test
	public void bandsKeepTheirOrder() {
		final double[] bg = new double[401];
		final double[] p = fit(0, 10, bg);
		assertTrue(p[3] < p[6] && p[6] < p[9]);
	}
}
