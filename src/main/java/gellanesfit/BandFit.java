/*
 * Gel Lanes Fit - BandFit.java
 * Author: Rick Ziraldo, 2017
 * The University of Texas at Dallas, Richardson, TX
 *
 * Licensed under the GNU Affero General Public License v3.0; see LICENSE.
 * Source: https://github.com/rickud/gel-lanes-fit
 */

package gellanesfit;

import java.util.Arrays;
import java.util.concurrent.CancellationException;
import java.util.function.BooleanSupplier;

import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.optim.ConvergenceChecker;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Pair;

/**
 * Fits distinct bands: a Gaussian peak per band on a polynomial background.
 * <p>
 * The background must stay below the profile's lower envelope: the lowest
 * values within a few band widths, smoothed. It can follow a sloped or hazy
 * background, or a smear under the bands, but can't rise into the bands.
 * That ceiling is a penalty in the least-squares sum, so the optimizer can
 * follow it. The background's average slope (Max Polynomial Derivative) is a
 * parameter of its own, bounded by tanh, so it never exceeds the limit. Each band's height, position and width are kept
 * within bounds of their starting values by mapping them through tanh: a
 * band moves at most half way to its nearest neighbour, so bands keep their
 * order; its width stays between 0.4 and 2 times its start; its height
 * between 0 and twice the profile above the envelope at its start.
 * </p>
 */
final class BandFit {

	/** A band moves at most this fraction of the gap to its nearest neighbour */
	static final double MAX_SHIFT = 0.5;
	static final double MIN_WIDTH = 0.4;
	static final double MAX_WIDTH = 2.0;
	/** Highest band, as a multiple of the profile above the envelope */
	static final double MAX_HEIGHT = 2.0;
	/** The envelope's window, in starting band widths (standard deviations) */
	static final double ENVELOPE_WINDOW = 10;
	/** Weight of the penalty terms, against misfits of one noise level */
	private static final double PENALTY = 30;
	private static final double CUTOFF = 8;
	/** The fit stops when an iteration improves the cost by less than this */
	static final double COST_TOLERANCE = 1e-6;
	static final int MAX_ITERATIONS = 200;

	/** The fit's result */
	static final class Result {

		/** Degree, background coefficients in x (constant first), then height,
		 * position and width per band, as Fitter's parameter arrays */
		double[] parameters;
		/** Root mean square of the profile minus the fit */
		double rms;
		int iterations;
	}

	private final double[] x, y;
	private final int deg, n;
	private final double[] mean0, sd0, maxShift, maxHeight;
	private final double[] envelope;
	private final double maxSlope;
	private final double noise, centre, half;
	private final BooleanSupplier stopCheck;

	/**
	 * @param x positions of the profile, ascending
	 * @param y the profile
	 * @param deg the background's degree, -1 for none
	 * @param means starting band positions, in order
	 * @param sds starting band widths
	 * @param maxSlope Max Polynomial Derivative: the background's average slope
	 *          limit, in gray values per px
	 */
	BandFit(final double[] x, final double[] y, final int deg,
		final double[] means, final double[] sds, final double maxSlope,
		final BooleanSupplier stopCheck)
	{
		this.x = x;
		this.y = y;
		this.deg = deg;
		this.n = means.length;
		this.mean0 = means;
		this.sd0 = sds;
		this.maxSlope = maxSlope;
		this.stopCheck = stopCheck;
		centre = (x[0] + x[x.length - 1]) / 2;
		half = FastMath.max(1e-9, (x[x.length - 1] - x[0]) / 2);
		noise = ContinuumFit.noiseLevel(y);
		envelope = lowerEnvelope(y, (int) FastMath.round(ENVELOPE_WINDOW *
			median(sds)));
		maxShift = new double[n];
		maxHeight = new double[n];
		for (int i = 0; i < n; i++) {
			double gap = Double.POSITIVE_INFINITY;
			if (i > 0) gap = FastMath.min(gap, means[i] - means[i - 1]);
			if (i < n - 1) gap = FastMath.min(gap, means[i + 1] - means[i]);
			if (Double.isInfinite(gap)) gap = 2 * half;
			maxShift[i] = MAX_SHIFT * FastMath.abs(gap);
			maxHeight[i] = MAX_HEIGHT * FastMath.max(valueAt(y, means[i]) -
				valueAt(envelope, means[i]), noise);
		}
	}

	private static double median(final double[] v) {
		if (v.length == 0) return 1;
		final double[] s = v.clone();
		Arrays.sort(s);
		return FastMath.max(1, s[s.length / 2]);
	}

	/**
	 * The profile's lower envelope: a morphological opening (the lowest value
	 * within half a window on either side, then the highest of those), which
	 * removes anything narrower than the window, such as bands, and follows a
	 * straight slope exactly; then averaged over the same window, and never
	 * above the profile
	 */
	static double[] lowerEnvelope(final double[] y, final int window) {
		final int w = FastMath.max(1, window / 2);
		final double[] min = slide(y, w, true);
		final double[] open = slide(min, w, false);
		extendEnds(open, w);
		final double[] env = new double[y.length];
		for (int j = 0; j < y.length; j++) {
			double s = 0;
			int c = 0;
			for (int k = FastMath.max(0, j - w); k <= FastMath.min(y.length - 1, j +
				w); k++)
			{
				s += open[k];
				c++;
			}
			env[j] = FastMath.min(s / c, y[j]);
		}
		return env;
	}

	/**
	 * Within w points of either end, the opening flattens a slope (there's
	 * nothing beyond the end to take the highest value from): there, it's
	 * replaced by a straight line through the next w points inwards
	 */
	private static void extendEnds(final double[] v, final int w) {
		final int n = v.length;
		if (n < 3 * w + 2) return;
		// Left end: line through points w .. 2w
		double slope = (v[2 * w] - v[w]) / w;
		for (int j = 0; j < w; j++)
			v[j] = v[w] - slope * (w - j);
		// Right end: line through points n-1-2w .. n-1-w
		slope = (v[n - 1 - w] - v[n - 1 - 2 * w]) / w;
		for (int j = n - w; j < n; j++)
			v[j] = v[n - 1 - w] + slope * (j - (n - 1 - w));
	}

	/** The lowest (or highest) value within w points on either side */
	private static double[] slide(final double[] v, final int w,
		final boolean lowest)
	{
		final double[] out = new double[v.length];
		for (int j = 0; j < v.length; j++) {
			double m = lowest ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
			for (int k = FastMath.max(0, j - w); k <= FastMath.min(v.length - 1, j +
				w); k++)
				m = lowest ? FastMath.min(m, v[k]) : FastMath.max(m, v[k]);
			out[j] = m;
		}
		return out;
	}

	private double valueAt(final double[] v, final double pos) {
		int j = Arrays.binarySearch(x, pos);
		if (j >= 0) return v[j];
		j = -j - 1;
		if (j <= 0) return v[0];
		if (j >= x.length) return v[x.length - 1];
		final double t = (pos - x[j - 1]) / (x[j] - x[j - 1]);
		return v[j - 1] + t * (v[j] - v[j - 1]);
	}

	private int bgCount() {
		return deg + 1;
	}

	/** The background's starting coefficients (in t): a fit to the envelope,
	 * lowered until it's under it */
	private double[] startingBackground() {
		final int m = bgCount();
		if (m == 0) return new double[0];
		final WeightedObservedPoints obs = new WeightedObservedPoints();
		for (int j = 0; j < x.length; j++)
			obs.add((x[j] - centre) / half, envelope[j]);
		final double[] c = PolynomialCurveFitter.create(deg).fit(obs.toList());
		double over = 0;
		for (int j = 0; j < x.length; j++)
			over = FastMath.max(over, poly(c, (x[j] - centre) / half) - envelope[j]);
		c[0] -= over;
		// Within the slope limit
		final double slope = (poly(c, 1) - poly(c, -1)) / (2 * half);
		if (FastMath.abs(slope) > maxSlope && deg >= 1) {
			final double f = maxSlope / FastMath.abs(slope);
			for (int k = 1; k < c.length; k++)
				c[k] *= f;
			double o = 0;
			for (int j = 0; j < x.length; j++)
				o = FastMath.max(o, poly(c, (x[j] - centre) / half) - envelope[j]);
			c[0] -= o;
		}
		return c;
	}

	private static double poly(final double[] c, final double t) {
		double v = 0;
		for (int k = c.length - 1; k >= 0; k--)
			v = v * t + c[k];
		return v;
	}

	Result fit() {
		final int m = bgCount();
		final double[] start = new double[m + 3 * n];
		System.arraycopy(toParameters(startingBackground()), 0, start, 0, m);
		for (int i = 0; i < n; i++) {
			// Height: start at the profile above the envelope, half the upper bound
			start[m + 3 * i] = 0;
			start[m + 3 * i + 1] = 0;
			start[m + 3 * i + 2] = FastMath.atanh(FastMath.max(-0.999, FastMath.min(
				0.999, -widthOffset() / widthScale())));
		}
		final int rows = x.length + x.length;
		final double[] target = new double[rows];
		for (int j = 0; j < x.length; j++)
			target[j] = y[j] / noise;
		// Stops when an iteration barely improves the fit, or after
		// MAX_ITERATIONS, keeping the point reached: a band the profile doesn't
		// support only approaches zero height, and the library's own criteria
		// may then never be met
		final ConvergenceChecker<LeastSquaresProblem.Evaluation> checker = (
			iteration, previous, current) -> iteration >= MAX_ITERATIONS || FastMath
				.abs(previous.getCost() - current.getCost()) <= COST_TOLERANCE *
					current.getCost();
		final LeastSquaresOptimizer.Optimum opt = new LevenbergMarquardtOptimizer()
			.optimize(new LeastSquaresBuilder().start(start).model(
				this::valueAndJacobian).target(target).checker(checker)
				.maxEvaluations(Integer.MAX_VALUE).maxIterations(Integer.MAX_VALUE)
				.lazyEvaluation(false).build());
		final double[] p = opt.getPoint().toArray();

		final Result r = new Result();
		r.iterations = opt.getIterations();
		final double[] coeffsT = coefficients(p);
		final double[] coeffsX = toX(coeffsT);
		final double[] out = new double[1 + m + 3 * n];
		out[0] = deg;
		System.arraycopy(coeffsX, 0, out, 1, m);
		double ss = 0;
		final double[] h = new double[n], mu = new double[n], sd = new double[n];
		for (int i = 0; i < n; i++) {
			h[i] = height(i, p[m + 3 * i]);
			mu[i] = mean(i, p[m + 3 * i + 1]);
			sd[i] = width(i, p[m + 3 * i + 2]);
			out[1 + m + 3 * i] = h[i];
			out[2 + m + 3 * i] = mu[i];
			out[3 + m + 3 * i] = sd[i];
		}
		for (int j = 0; j < x.length; j++) {
			double v = poly(coeffsT, (x[j] - centre) / half);
			for (int i = 0; i < n; i++) {
				final double z = (x[j] - mu[i]) / sd[i];
				if (FastMath.abs(z) < CUTOFF) v += h[i] * FastMath.exp(-0.5 * z * z);
			}
			ss += (y[j] - v) * (y[j] - v);
		}
		r.rms = FastMath.sqrt(ss / x.length);
		r.parameters = out;
		return r;
	}

	// Bounded parameters

	/**
	 * The background's parameters: its coefficients in t, except that the
	 * linear one is replaced by v, with the average slope maxSlope * tanh(v).
	 * The average slope over the lane is (c1 + c3 + c5 + ...) / half.
	 */
	private double[] toParameters(final double[] c) {
		final double[] p = c.clone();
		if (c.length > 1) {
			final double slope = oddSum(c) / half;
			p[1] = maxSlope > 0 ? FastMath.atanh(FastMath.max(-0.999, FastMath.min(
				0.999, slope / maxSlope))) : 0;
		}
		return p;
	}

	/** The background's coefficients in t from the parameters */
	private double[] coefficients(final double[] p) {
		final double[] c = Arrays.copyOf(p, bgCount());
		if (c.length > 1) {
			c[1] = 0;
			c[1] = half * maxSlope * FastMath.tanh(p[1]) - oddSum(c);
		}
		return c;
	}

	private static double oddSum(final double[] c) {
		double s = 0;
		for (int k = 1; k < c.length; k += 2)
			s += c[k];
		return s;
	}

	/**
	 * d(background at t)/d(parameter k): t^k for the constant and even terms;
	 * for v, t times d(c1)/dv; for the other odd terms, t^k - t, since c1
	 * makes up for them
	 */
	private double[] backgroundGradient(final double[] p, final double t) {
		final int m = bgCount();
		final double[] g = new double[m];
		double tk = 1;
		for (int k = 0; k < m; k++) {
			g[k] = tk;
			tk *= t;
		}
		if (m > 1) {
			final double th = FastMath.tanh(p[1]);
			g[1] = t * half * maxSlope * (1 - th * th);
			for (int k = 3; k < m; k += 2)
				g[k] -= t;
		}
		return g;
	}
	private double height(final int i, final double u) {
		return maxHeight[i] * (1 + FastMath.tanh(u)) / 2;
	}

	private double dHeight(final int i, final double u) {
		final double t = FastMath.tanh(u);
		return maxHeight[i] * (1 - t * t) / 2;
	}

	private double mean(final int i, final double u) {
		return mean0[i] + maxShift[i] * FastMath.tanh(u);
	}

	private double dMean(final int i, final double u) {
		final double t = FastMath.tanh(u);
		return maxShift[i] * (1 - t * t);
	}

	/** log(width / start width) = offset + scale * tanh(u) */
	private static double widthOffset() {
		return (FastMath.log(MAX_WIDTH) + FastMath.log(MIN_WIDTH)) / 2;
	}

	private static double widthScale() {
		return (FastMath.log(MAX_WIDTH) - FastMath.log(MIN_WIDTH)) / 2;
	}

	private double width(final int i, final double u) {
		return sd0[i] * FastMath.exp(widthOffset() + widthScale() * FastMath.tanh(
			u));
	}

	/** d(width)/du */
	private double dWidth(final int i, final double u) {
		final double t = FastMath.tanh(u);
		return width(i, u) * widthScale() * (1 - t * t);
	}

	/**
	 * Rows: the profile (over the noise), and the envelope penalty at each
	 * point
	 */
	private Pair<RealVector, RealMatrix> valueAndJacobian(final RealVector point) {
		if (stopCheck.getAsBoolean()) throw new CancellationException(
			"The fit was stopped");
		final double[] p = point.toArray();
		final int m = bgCount();
		final int cols = p.length;
		final int nx = x.length;
		final double[] value = new double[2 * nx];
		final double[][] jac = new double[2 * nx][cols];
		final double[] c = coefficients(p);
		for (int j = 0; j < nx; j++) {
			final double t = (x[j] - centre) / half;
			final double bg = poly(c, t);
			final double[] g = backgroundGradient(p, t);
			for (int k = 0; k < m; k++)
				jac[j][k] = g[k] / noise;
			value[j] = bg / noise;
			// Envelope penalty: only where the background is above it
			final double over = bg - envelope[j];
			if (over > 0) {
				value[nx + j] = PENALTY * over / noise;
				for (int k = 0; k < m; k++)
					jac[nx + j][k] = PENALTY * g[k] / noise;
			}
		}
		for (int i = 0; i < n; i++) {
			final double uh = p[m + 3 * i], um = p[m + 3 * i + 1], us = p[m + 3 * i +
				2];
			final double h = height(i, uh), mu = mean(i, um), sd = width(i, us);
			final double dh = dHeight(i, uh), dm = dMean(i, um), ds = dWidth(i, us);
			final double reach = CUTOFF * sd;
			int j = Arrays.binarySearch(x, mu - reach);
			if (j < 0) j = -j - 1;
			for (; j < nx && x[j] <= mu + reach; j++) {
				final double z = (x[j] - mu) / sd;
				final double e = FastMath.exp(-0.5 * z * z) / noise;
				value[j] += h * e;
				jac[j][m + 3 * i] = e * dh;
				jac[j][m + 3 * i + 1] = h * e * z / sd * dm;
				jac[j][m + 3 * i + 2] = h * e * z * z / sd * ds;
			}
		}
		return new Pair<>(new ArrayRealVector(value, false),
			new Array2DRowRealMatrix(jac, false));
	}

	/** Coefficients in t = (x - centre) / half, rewritten in x */
	private double[] toX(final double[] ct) {
		final int m = ct.length;
		final double[] cx = new double[m];
		// sum_k ct[k] ((x - centre) / half)^k, expanded binomially
		for (int k = 0; k < m; k++) {
			final double a = ct[k] / FastMath.pow(half, k);
			double binom = 1;
			for (int r = 0; r <= k; r++) {
				// term: a * C(k, r) * x^r * (-centre)^(k - r)
				cx[r] += a * binom * FastMath.pow(-centre, k - r);
				binom = binom * (k - r) / (r + 1);
			}
		}
		return cx;
	}
}
