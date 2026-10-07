package gellanesfit;

import static org.junit.Assert.assertEquals;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.junit.Test;

/** The fitted curve: Gaussian peaks on a polynomial background. */
public class GaussianArrayTest {

	private static GaussianArray curve(final double... poly) {
		return new GaussianArray(new SortedParameters(new ArrayRealVector(poly),
			new ArrayRealVector(new double[] { 100, 50 }), new ArrayRealVector(
				new double[] { 200, 300 }), new ArrayRealVector(new double[] { 5,
					8 })));
	}

	@Test
	public void valueIsPeaksPlusBackground() {
		final GaussianArray g = curve(20, 0.02);
		assertEquals(20 + 0.02 * 200 + 100 + 50 * Math.exp(-0.5 * 10000.0 / 64), g
			.value(200), 1e-9);
	}

	@Test
	public void derivativesMatchTheValueAndItsSlope() {
		final GaussianArray g = curve(20, 0.02, -1e-5);
		for (final double x : new double[] { 150, 197, 200, 260, 310 }) {
			final DerivativeStructure t = new DerivativeStructure(1, 2, 0, x);
			final DerivativeStructure v = g.value(t);
			assertEquals(g.value(x), v.getValue(), 1e-9);
			final double h = 1e-4;
			assertEquals("first derivative at " + x, (g.value(x + h) - g.value(x -
				h)) / (2 * h), v.getPartialDerivative(1), 1e-5);
		}
	}

	@Test
	public void noBackground() {
		final GaussianArray g = curve(); // Polynomial Degree -1
		final DerivativeStructure v = g.value(new DerivativeStructure(1, 1, 0,
			200));
		assertEquals(100 + 50 * Math.exp(-0.5 * 10000.0 / 64), v.getValue(), 1e-9);
		assertEquals(v.getValue(), g.value(200), 1e-9);
	}
}
