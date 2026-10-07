package gellanesfit;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

/** Synthetic lane profiles with known bands, for the tests. */
final class TestProfiles {

	private TestProfiles() {}

	/** A band: height above the background, position and standard deviation */
	static double[] band(final double height, final double mean,
		final double sd)
	{
		return new double[] { height, mean, sd };
	}

	/**
	 * Profile values at x = x0, x0 + 1, ..., x1: the bands plus a linear
	 * background a + b * x.
	 */
	static double[] values(final int x0, final int x1, final double a,
		final double b, final double[]... bands)
	{
		final double[] y = new double[x1 - x0 + 1];
		for (int i = 0; i < y.length; i++) {
			final double x = x0 + i;
			y[i] = a + b * x;
			for (final double[] g : bands)
				y[i] += new Gaussian(g[0], g[1], g[2]).value(x);
		}
		return y;
	}

	static double[] xs(final int x0, final int x1) {
		final double[] x = new double[x1 - x0 + 1];
		for (int i = 0; i < x.length; i++)
			x[i] = x0 + i;
		return x;
	}

	static List<WeightedObservedPoint> points(final int x0, final double[] y) {
		final List<WeightedObservedPoint> obs = new ArrayList<>();
		for (int i = 0; i < y.length; i++)
			obs.add(new WeightedObservedPoint(1, x0 + i, y[i]));
		return obs;
	}

	/** A lane profile as the plugin builds it, x being the distance in px */
	static DataSeries profile(final int lane, final int x0, final double[] y) {
		final RealVector x = new ArrayRealVector(xs(x0, x0 + y.length - 1));
		return new DataSeries("Lane " + lane, lane, DataSeries.PROFILE, x,
			new ArrayRealVector(y), Color.BLACK);
	}
}
