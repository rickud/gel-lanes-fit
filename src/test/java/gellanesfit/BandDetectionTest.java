package gellanesfit;

import static gellanesfit.TestProfiles.band;
import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import org.junit.Test;

/** The starting guess of a Banded fit: one Gaussian per detected band. */
public class BandDetectionTest {

	private static final int X0 = 100;
	private static final int X1 = 500;

	private static SortedParameters guess(final double[] y,
		final double absoluteTolerance)
	{
		return new GaussianArrayCurveFitter.ParameterGuesser(TestProfiles.points(
			X0, y), 2, absoluteTolerance, 0.98).guess();
	}

	@Test
	public void findsSeparatedBandsAtTheirMaxima() {
		final double[] y = TestProfiles.values(X0, X1, 20, 0, band(100, 180, 6),
			band(60, 300, 8), band(80, 420, 5));
		final SortedParameters g = guess(y, 0.1 * 140);
		assertArrayEquals(new double[] { 180, 300, 420 }, g.getMean().toArray(),
			0.5);
	}

	@Test
	public void startingWidthComesFromTheHalfMaximum() {
		final double[] y = TestProfiles.values(X0, X1, 20, 0, band(100, 300, 8));
		final SortedParameters g = guess(y, 10);
		assertEquals(1, g.getMean().getDimension());
		assertEquals(8, g.getSD().getEntry(0), 1.5);
	}

	@Test
	public void toleranceDecidesWhichBandsCount() {
		// A weak band (height 10) next to two strong ones
		final double[] y = TestProfiles.values(X0, X1, 20, 0, band(100, 180, 6),
			band(10, 300, 6), band(100, 420, 6));
		assertEquals(3, guess(y, 5).getMean().getDimension());
		assertEquals(2, guess(y, 20).getMean().getDimension());
	}
}
