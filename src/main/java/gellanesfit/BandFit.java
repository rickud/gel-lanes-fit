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
 * The background may also never bend upward: its second derivative must be
 * zero or negative everywhere (a straight line, or a cap shape), so that a
 * high degree can't add wiggles. The ceiling and the curvature limit are
 * penalties in the least-squares sum, so the optimizer can follow them. The background's average slope (Max Polynomial Derivative) is a
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
	/** Iterations in a row that must barely improve the fit to stop it */
	static final int CALM_ITERATIONS = 5;

	/** The fit's result */
	static final class Result {

		/** Degree, background coefficients in x (constant first), then height,
		 * position and width per band, as Fitter's parameter arrays */
		double[] parameters;
		/** Root mean square of the profile minus the fit */
		double rms;
		int iterations;
		/** The starting background's coefficients in x, constant first */
		double[] startBackground;
	}

	private final double[] x, y;
	private final int deg, n;
	private final double[] mean0, sd0, maxShift, maxHeight;
	private final double[] envelope;
	/** Half the envelope's window, in points */
	private final int envelopeHalf;
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
		final int window = (int) FastMath.round(ENVELOPE_WINDOW * median(sds));
		envelope = lowerEnvelope(y, window);
		envelopeHalf = FastMath.max(1, window / 2);
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
	 * straight slope exactly; then averaged over the same window (centred, so
	 * narrower near the ends), and never above the profile. Within half a
	 * window of either end, where an opening can't follow a slope, it's
	 * continued in a straight line with the overall slope.
	 */
	static double[] lowerEnvelope(final double[] y, final int window) {
		final int w = FastMath.max(1, window / 2);
		final double[] min = slide(y, w, true);
		final double[] open = slide(min, w, false);
		// Within half a window of either end, the opening flattens a slope
		// (there's nothing beyond the end to take the highest value from):
		// it's used from w to n-1-w only, averaged over a window centred on
		// each point and kept within that range, so that a straight slope
		// averages to itself
		final int n = y.length;
		final int first = w, last = n - 1 - w;
		final double[] env = new double[n];
		if (last - first < 2) { // Too short a lane for the window
			double low = Double.POSITIVE_INFINITY;
			for (final double v : y)
				low = FastMath.min(low, v);
			Arrays.fill(env, low);
			return env;
		}
		for (int j = first; j <= last; j++) {
			final int h = FastMath.min(w, FastMath.min(j - first, last - j));
			double s = 0;
			for (int k = j - h; k <= j + h; k++)
				s += open[k];
			env[j] = FastMath.min(s / (2 * h + 1), y[j]);
		}
		// Near the ends: a straight line with the slope of the whole reliable
		// part, from its first and last points, under the profile. Bands near
		// an end stay above it, and a straight background is followed exactly.
		double sx = 0, sy = 0, sxx = 0, sxy = 0;
		for (int j = first; j <= last; j++) {
			sx += j;
			sy += env[j];
			sxx += (double) j * j;
			sxy += j * env[j];
		}
		final double count = last - first + 1;
		final double slope = (count * sxy - sx * sy) / (count * sxx - sx * sx);
		for (int j = 0; j < first; j++)
			env[j] = FastMath.min(y[j], env[first] + slope * (j - first));
		for (int j = last + 1; j < n; j++)
			env[j] = FastMath.min(y[j], env[last] + slope * (j - last));
		return env;
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
		// Fitted where the envelope is an opening, away from the ends
		final WeightedObservedPoints obs = new WeightedObservedPoints();
		final int skip = x.length > 4 * envelopeHalf ? envelopeHalf : 0;
		for (int j = skip; j < x.length - skip; j++)
			obs.add((x[j] - centre) / half, envelope[j]);
		double[] c = PolynomialCurveFitter.create(deg).fit(obs.toList());
		if (bendsUpward(c)) { // Not allowed: start from a straight line
			c = Arrays.copyOf(PolynomialCurveFitter.create(FastMath.min(deg, 1))
				.fit(obs.toList()), m);
		}
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

	/** Whether the background's second derivative is positive anywhere */
	private boolean bendsUpward(final double[] c) {
		for (int j = 0; j < x.length; j++)
			if (curvature(c, (x[j] - centre) / half) > 1e-9 * noise) return true;
		return false;
	}

	/** The second derivative in t of the polynomial with these coefficients */
	private static double curvature(final double[] c, final double t) {
		double v = 0;
		for (int k = c.length - 1; k >= 2; k--)
			v = v * t + k * (k - 1) * c[k];
		return v;
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
		final double[] startT = startingBackground();
		System.arraycopy(toParameters(startT), 0, start, 0, m);
		for (int i = 0; i < n; i++) {
			// Height: start at the profile above the envelope, half the upper bound
			start[m + 3 * i] = 0;
			start[m + 3 * i + 1] = 0;
			start[m + 3 * i + 2] = FastMath.atanh(FastMath.max(-0.999, FastMath.min(
				0.999, -widthOffset() / widthScale())));
		}
		final int rows = 3 * x.length;
		final double[] target = new double[rows];
		for (int j = 0; j < x.length; j++)
			target[j] = y[j] / noise;
		// Stops when CALM_ITERATIONS iterations in a row barely improve the fit,
		// or after MAX_ITERATIONS, keeping the point reached: a band the profile
		// doesn't support only approaches zero height, and the library's own
		// criteria may then never be met. One small step isn't enough: the
		// first steps can be small while the optimizer finds its scale.
		final int[] calm = { 0 };
		final ConvergenceChecker<LeastSquaresProblem.Evaluation> checker = (
			iteration, previous, current) -> {
			if (iteration >= MAX_ITERATIONS) return true;
			if (FastMath.abs(previous.getCost() - current.getCost()) <=
				COST_TOLERANCE * current.getCost()) calm[0]++;
			else calm[0] = 0;
			return calm[0] >= CALM_ITERATIONS;
		};
		final LeastSquaresOptimizer.Optimum opt = new LevenbergMarquardtOptimizer()
			.optimize(new LeastSquaresBuilder().start(start).model(
				this::valueAndJacobian).target(target).checker(checker)
				.maxEvaluations(Integer.MAX_VALUE).maxIterations(Integer.MAX_VALUE)
				.lazyEvaluation(false).build());
		final double[] p = opt.getPoint().toArray();

		final Result r = new Result();
		r.iterations = opt.getIterations();
		r.startBackground = toX(startT);
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
	 * Rows: the profile (over the noise), the envelope penalty and the
	 * curvature penalty at each point
	 */
	private Pair<RealVector, RealMatrix> valueAndJacobian(final RealVector point) {
		if (stopCheck.getAsBoolean()) throw new CancellationException(
			"The fit was stopped");
		final double[] p = point.toArray();
		final int m = bgCount();
		final int cols = p.length;
		final int nx = x.length;
		final double[] value = new double[3 * nx];
		final double[][] jac = new double[3 * nx][cols];
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
			// Curvature penalty: only where the background bends upward. The
			// linear term doesn't enter the curvature, so the slope parameter
			// and the odd terms' share of c1 don't either.
			final double bend = curvature(c, t);
			if (bend > 0) {
				value[2 * nx + j] = PENALTY * bend / noise;
				double tk = 1;
				for (int k = 2; k < m; k++) {
					jac[2 * nx + j][k] = PENALTY * k * (k - 1) * tk / noise;
					tk *= t;
				}
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
