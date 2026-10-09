package gellanesfit;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Random;

import org.apache.commons.math3.analysis.function.Gaussian;
import org.junit.Test;

/** The smear fitter on its own, on profiles made from known peaks. */
public class ContinuumFitSolverTest {

	@Test
	public void nonNegativeLeastSquaresKeepsSolutionsAtOrAboveZero() {
		// Unconstrained, z = (2, -1); with z >= 0 the best is (1.5, 0)
		final double[][] m = { { 1, 0.5 }, { 0.5, 1 } };
		final double[] q = { 1.5, 0 };
		assertArrayEquals(new double[] { 1.5, 0 }, ContinuumFit.nnls(m, q),
			1e-9);
		// Without a bound in the way it's the plain solution
		assertArrayEquals(new double[] { 1, 1 }, ContinuumFit.nnls(m,
			new double[] { 1.5, 1.5 }), 1e-9);
	}

	@Test
	public void noiseLevelOfAProfile() {
		final Random random = new Random(1);
		final double[] y = new double[2000];
		for (int j = 0; j < y.length; j++)
			y[j] = 100 + 3 * random.nextGaussian();
		assertEquals(3, ContinuumFit.noiseLevel(y), 0.3);
	}

	/** A smear of dense peaks with smoothly varying areas, plus noise */
	private static double[] smear(final double[] x, final double[] means,
		final double[] areas, final double sd, final double noise)
	{
		final Random random = new Random(2);
		final double[] y = new double[x.length];
		for (int j = 0; j < x.length; j++) {
			for (int i = 0; i < means.length; i++)
				y[j] += areas[i] * new Gaussian(means[i], sd).value(x[j]);
			y[j] += noise * random.nextGaussian();
		}
		return y;
	}

	@Test
	public void recoversASmoothSmearWithinTheNoise() {
		final double[] x = new double[401];
		for (int j = 0; j < x.length; j++)
			x[j] = 100 + j;
		final int n = 100;
		final double[] means = new double[n];
		final double[] sds = new double[n];
		final double[] areas = new double[n];
		final double[] expected = new double[n];
		for (int i = 0; i < n; i++) {
			means[i] = 100 + 4 * i;
			sds[i] = 4;
			areas[i] = 400 * (1 + Math.sin(i / 15.0)) + 100;
			expected[i] = 1;
		}
		final double[] y = smear(x, means, areas, 4, 2);
		final ContinuumFit.Result r = new ContinuumFit(x, y, means, sds, expected,
			1, 10, () -> false).fit();
		assertEquals("fits to the noise level", 2, r.rms, 0.5);
		assertEquals("no offset", 0, r.offset, 0.5);
		for (int i = 5; i < n - 5; i++)
			assertEquals("area of peak " + i, areas[i], r.areas[i], 0.1 *
				areas[i]);
	}

	@Test
	public void smoothnessKeepsNeighbouringAreasAlike() {
		final double[] x = new double[401];
		for (int j = 0; j < x.length; j++)
			x[j] = 100 + j;
		final int n = 100;
		final double[] means = new double[n];
		final double[] sds = new double[n];
		final double[] areas = new double[n];
		final double[] expected = new double[n];
		for (int i = 0; i < n; i++) {
			means[i] = 100 + 4 * i;
			sds[i] = 4;
			areas[i] = 500;
			expected[i] = 1;
		}
		final double[] y = smear(x, means, areas, 4, 5);
		final double rough = roughness(new ContinuumFit(x, y, means, sds,
			expected, 1, 0.01, () -> false).fit().areas);
		final double smooth = roughness(new ContinuumFit(x, y, means, sds,
			expected, 1, 10, () -> false).fit().areas);
		assertTrue("smoother: " + smooth + " vs " + rough, smooth < rough / 3);
	}

	/** Mean absolute second difference, relative to the mean */
	private static double roughness(final double[] a) {
		double s = 0;
		double m = 0;
		for (int i = 1; i < a.length - 1; i++)
			s += Math.abs(a[i - 1] - 2 * a[i] + a[i + 1]);
		for (final double v : a)
			m += v;
		return s / (a.length - 2) / (m / a.length);
	}
}
