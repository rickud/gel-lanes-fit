/*
 * Gel Lanes Fit - ContinuumFit.java
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

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.CholeskyDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.SimplexOptimizer;
import org.apache.commons.math3.util.FastMath;

/**
 * Fits a smear: a sum of Gaussian peaks, one per fragment length, whose
 * positions and widths come from the ladder.
 * <p>
 * The peaks keep their order and spacing: only a lane-wide offset and stretch
 * of the positions, and a lane-wide scale of the widths (within SD Drift),
 * are fitted, plus each peak's area. With the positions and widths set, the
 * areas enter the fit linearly, so they're solved directly, as non-negative
 * least squares; the offset, stretch and width scale are searched around
 * them.
 * </p>
 * <p>
 * A smoothness penalty keeps neighbouring peaks alike unless the profile says
 * otherwise: it adds the curvature, across neighbouring peaks, of each area
 * relative to its expected and starting values to the least-squares sum, so
 * a sharp change between peaks only survives where the profile has one.
 * </p>
 * <p>
 * The profile is expected with its background already taken off.
 * </p>
 */
final class ContinuumFit {

	/** Largest lane offset, as a fraction of the lane length */
	static final double MAX_OFFSET = 0.05;
	/** Largest stretch of the positions about the lane's centre */
	static final double MAX_STRETCH = 0.05;
	/** Gaussians are taken as zero beyond this many standard deviations */
	private static final double CUTOFF = 6;
	private static final double SQRT_2PI = FastMath.sqrt(2 * FastMath.PI);

	/** The fitted peaks and how the fit went */
	static final class Result {

		double[] areas, means, sds;
		/**
		 * The heights the fit starts from: the profile's level at each peak,
		 * shared out among overlapping peaks (see {@link ContinuumFit#startAreas})
		 */
		double[] startHeights;
		double offset, stretch, widthScale;
		/** Root mean square of the profile minus the fit */
		double rms;
		/** Number of area solves the search took */
		int iterations;

		/** A peak's height: its area over its width times the square root of 2 pi */
		double height(final int i) {
			return areas[i] / (sds[i] * SQRT_2PI);
		}
	}

	private final double[] x, y;
	private final double[] mean0, sd0;
	/** Each area is penalized relative to this: expected times starting area */
	private final double[] reference;
	/**
	 * The starting areas: the profile's level at each peak times the smaller of
	 * its width (times the square root of 2 pi) and the spacing to its
	 * neighbours, so that overlapping peaks add up to the profile
	 */
	private final double[] startAreas;
	private final double maxOffset, centre, logSdDrift;
	private final double noise, weight;
	private final BooleanSupplier stopCheck;
	private final int n;
	private int solves;

	/**
	 * @param x positions of the profile, in px, ascending
	 * @param y the profile, background off
	 * @param mean0 each peak's position from the ladder, in order
	 * @param sd0 each peak's width from the ladder
	 * @param expected each peak's expected area, in any unit: only ratios count
	 * @param sdDrift widths may be scaled by 1 / sdDrift to sdDrift
	 * @param smoothness strength of the smoothness penalty; 0 for none
	 * @param stopCheck checked on every area solve; true stops the fit
	 */
	ContinuumFit(final double[] x, final double[] y, final double[] mean0,
		final double[] sd0, final double[] expected, final double sdDrift,
		final double smoothness, final BooleanSupplier stopCheck)
	{
		this.x = x;
		this.y = y;
		this.mean0 = mean0;
		this.sd0 = sd0;
		this.stopCheck = stopCheck;
		n = mean0.length;
		final double xMin = x[0];
		final double xMax = x[x.length - 1];
		maxOffset = MAX_OFFSET * (xMax - xMin);
		centre = (xMin + xMax) / 2;
		logSdDrift = FastMath.log(FastMath.max(1, sdDrift));
		noise = noiseLevel(y);
		// At smoothness 1, a relative curvature of 1 at every peak costs as much
		// as a misfit of one noise level at every point of the profile
		weight = n > 2 ? smoothness * FastMath.sqrt(x.length / (double) (n -
			2)) : 0;
		startAreas = new double[n];
		reference = referenceAreas(expected);
	}

	/**
	 * The profile's noise: the robust spread of the differences between
	 * neighbouring points, over the square root of 2. Never 0.
	 */
	static double noiseLevel(final double[] y) {
		final double[] d = new double[y.length - 1];
		double range = 0;
		for (int j = 0; j < d.length; j++) {
			d[j] = FastMath.abs(y[j + 1] - y[j]);
			range = FastMath.max(range, FastMath.abs(y[j]));
		}
		Arrays.sort(d);
		final double mad = d.length == 0 ? 0 : d[d.length / 2];
		return FastMath.max(1.4826 * mad / FastMath.sqrt(2), 1e-6 * FastMath.max(
			range, 1));
	}

	/**
	 * The areas the smoothness penalty compares each area with: the expected
	 * area, rescaled smoothly along the lane to the profile's level. Without
	 * the rescaling, the penalty would favour the distribution's shape over
	 * the smear's; with it, it only discourages changes from one peak to the
	 * next that the profile doesn't show.
	 */
	private double[] referenceAreas(final double[] expected) {
		double max = 0;
		for (final double v : y)
			max = FastMath.max(max, v);
		final double[] ref = new double[n];
		for (int i = 0; i < n; i++) {
			final double width = FastMath.min(spacing(i), sd0[i] * SQRT_2PI);
			// The profile smoothed over the peak's width, so noise doesn't count
			final double level = FastMath.max(smoothedProfileAt(mean0[i], 2 *
				sd0[i]), 1e-3 * max) + 1e-300;
			ref[i] = level * width;
			startAreas[i] = ref[i];
		}
		// The expected areas' shape, without their scale
		double meanLogRatio = 0;
		for (int i = 0; i < n; i++)
			meanLogRatio += FastMath.log(FastMath.max(expected[i], 1e-300) /
				ref[i]);
		meanLogRatio /= n;
		for (int i = 0; i < n; i++) {
			// Half way, in log, between the profile's level and the expected
			// area: the penalty then pulls gently towards the distribution
			final double e = FastMath.max(expected[i], 1e-300) / FastMath.exp(
				meanLogRatio);
			ref[i] = FastMath.sqrt(ref[i] * e);
		}
		return ref;
	}

	private double spacing(final int i) {
		final double left = i > 0 ? mean0[i] - mean0[i - 1] : Double.NaN;
		final double right = i < n - 1 ? mean0[i + 1] - mean0[i] : Double.NaN;
		final double s = Double.isNaN(left) ? right : Double.isNaN(right) ? left
			: (left + right) / 2;
		return Double.isNaN(s) ? sd0[i] * SQRT_2PI : FastMath.abs(s);
	}

	/** The mean of the profile within half a width of a position, 0 beyond it */
	private double smoothedProfileAt(final double pos, final double width) {
		double sum = 0;
		int count = 0;
		for (int j = lowerIndex(pos - width / 2); j < x.length && x[j] <= pos +
			width / 2; j++)
		{
			sum += y[j];
			count++;
		}
		return count == 0 ? 0 : FastMath.max(0, sum / count);
	}

	/** Index of the first profile position at or after pos */
	private int lowerIndex(final double pos) {
		int j = Arrays.binarySearch(x, pos);
		if (j < 0) j = -j - 1;
		return FastMath.max(0, j);
	}

	Result fit() {
		// Searched parameters, mapped into their bounds by tanh: offset and
		// stretch, and the width scale when SD Drift allows any
		final boolean widths = logSdDrift > 0;
		final int dims = widths ? 3 : 2;
		final MultivariateFunction cost = p -> solve(p)[0][0];
		final SimplexOptimizer optimizer = new SimplexOptimizer(1e-7, 1e-10);
		final PointValuePair best = optimizer.optimize(new MaxEval(400),
			new ObjectiveFunction(cost), GoalType.MINIMIZE, new InitialGuess(
				new double[dims]), new NelderMeadSimplex(dims, 0.3));
		final double[] p = best.getPoint();
		final double[][] solution = solve(p);

		final Result r = new Result();
		r.offset = offset(p);
		r.stretch = stretch(p);
		r.widthScale = widthScale(p);
		r.areas = solution[1];
		r.means = new double[n];
		r.sds = new double[n];
		r.startHeights = new double[n];
		for (int i = 0; i < n; i++) {
			r.startHeights[i] = startAreas[i] / (sd0[i] * SQRT_2PI);
			r.means[i] = mean0[i] + r.offset + r.stretch * (mean0[i] - centre);
			r.sds[i] = sd0[i] * r.widthScale;
		}
		double ss = 0;
		for (int j = 0; j < x.length; j++) {
			double m = 0;
			for (int i = 0; i < n; i++)
				m += gaussian(r.areas[i], r.means[i], r.sds[i], x[j]);
			ss += (y[j] - m) * (y[j] - m);
		}
		r.rms = FastMath.sqrt(ss / x.length);
		r.iterations = solves;
		return r;
	}

	private double offset(final double[] p) {
		return maxOffset * FastMath.tanh(p[0]);
	}

	private double stretch(final double[] p) {
		return MAX_STRETCH * FastMath.tanh(p[1]);
	}

	private double widthScale(final double[] p) {
		return p.length > 2 ? FastMath.exp(logSdDrift * FastMath.tanh(p[2])) : 1;
	}

	private static double gaussian(final double area, final double mean,
		final double sd, final double x)
	{
		final double z = (x - mean) / sd;
		if (FastMath.abs(z) > CUTOFF) return 0;
		return area / (sd * SQRT_2PI) * FastMath.exp(-0.5 * z * z);
	}

	/**
	 * The best areas for these offset, stretch and width scale.
	 *
	 * @return the cost (squared misfit over the noise, plus the penalty) as
	 *         {{cost}, areas}
	 */
	private double[][] solve(final double[] p) {
		if (stopCheck.getAsBoolean()) throw new CancellationException(
			"The fit was stopped");
		solves++;
		final double offset = offset(p);
		final double stretch = stretch(p);
		final double scale = widthScale(p);

		// Normal equations of the data rows (Gaussians of unit area, over the
		// noise) and the penalty rows (curvature of area / reference)
		final double[][] basis = new double[n][];
		final int[] first = new int[n];
		for (int i = 0; i < n; i++) {
			final double mean = mean0[i] + offset + stretch * (mean0[i] - centre);
			final double sd = sd0[i] * scale;
			first[i] = lowerIndex(mean - CUTOFF * sd);
			int last = first[i];
			while (last < x.length && x[last] <= mean + CUTOFF * sd)
				last++;
			basis[i] = new double[last - first[i]];
			for (int j = first[i]; j < last; j++)
				basis[i][j - first[i]] = gaussian(1, mean, sd, x[j]) / noise;
		}
		final double[][] ata = new double[n][n];
		final double[] atb = new double[n];
		double btb = 0;
		for (int j = 0; j < x.length; j++)
			btb += (y[j] / noise) * (y[j] / noise);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < basis[i].length; j++)
				atb[i] += basis[i][j] * y[first[i] + j] / noise;
			for (int k = i; k < n; k++) {
				final int lo = FastMath.max(first[i], first[k]);
				final int hi = FastMath.min(first[i] + basis[i].length, first[k] +
					basis[k].length);
				double s = 0;
				for (int j = lo; j < hi; j++)
					s += basis[i][j - first[i]] * basis[k][j - first[k]];
				ata[i][k] = s;
				ata[k][i] = s;
			}
		}
		for (int i = 1; i < n - 1; i++) {
			final int[] idx = { i - 1, i, i + 1 };
			final double[] c = { weight / reference[i - 1], -2 * weight /
				reference[i], weight / reference[i + 1] };
			for (int a = 0; a < 3; a++)
				for (int b = 0; b < 3; b++)
					ata[idx[a]][idx[b]] += c[a] * c[b];
		}
		final double[] areas = nnls(ata, atb);
		// cost = |A z - b|^2 = z'A'Az - 2 z'A'b + b'b
		double cost = btb;
		for (int i = 0; i < n; i++) {
			if (areas[i] == 0) continue;
			cost -= 2 * areas[i] * atb[i];
			for (int k = 0; k < n; k++)
				cost += areas[i] * ata[i][k] * areas[k];
		}
		return new double[][] { { cost }, areas };
	}

	/**
	 * Non-negative least squares from the normal equations (Lawson and
	 * Hanson's active set method): minimizes z'Mz - 2z'q for z >= 0.
	 */
	static double[] nnls(final double[][] m, final double[] q) {
		final int n = q.length;
		final double[] z = new double[n];
		final boolean[] passive = new boolean[n];
		for (int outer = 0; outer < 3 * n + 10; outer++) {
			// Gradient of the cost's negative half: q - Mz
			int best = -1;
			double bestW = 1e-12 * (1 + norm(q));
			for (int i = 0; i < n; i++) {
				if (passive[i]) continue;
				double w = q[i];
				for (int k = 0; k < n; k++)
					w -= m[i][k] * z[k];
				if (w > bestW) {
					bestW = w;
					best = i;
				}
			}
			if (best < 0) break;
			passive[best] = true;
			for (int inner = 0; inner < 3 * n + 10; inner++) {
				final double[] s = solvePassive(m, q, passive);
				double alpha = 1;
				int blocking = -1;
				for (int i = 0; i < n; i++)
					if (passive[i] && s[i] <= 0) {
						final double a = z[i] / (z[i] - s[i]);
						if (a < alpha) {
							alpha = a;
							blocking = i;
						}
					}
				if (blocking < 0) {
					System.arraycopy(s, 0, z, 0, n);
					break;
				}
				for (int i = 0; i < n; i++) {
					if (!passive[i]) continue;
					z[i] += alpha * (s[i] - z[i]);
					if (z[i] <= 1e-15 * (1 + FastMath.abs(s[i])) || i == blocking) {
						z[i] = 0;
						passive[i] = false;
					}
				}
			}
		}
		return z;
	}

	private static double norm(final double[] v) {
		double s = 0;
		for (final double e : v)
			s += e * e;
		return FastMath.sqrt(s);
	}

	/** Solves the normal equations on the passive set; 0 elsewhere */
	private static double[] solvePassive(final double[][] m, final double[] q,
		final boolean[] passive)
	{
		final int n = q.length;
		int k = 0;
		final int[] idx = new int[n];
		for (int i = 0; i < n; i++)
			if (passive[i]) idx[k++] = i;
		final double[][] sub = new double[k][k];
		final double[] rhs = new double[k];
		double diag = 0;
		for (int a = 0; a < k; a++) {
			rhs[a] = q[idx[a]];
			for (int b = 0; b < k; b++)
				sub[a][b] = m[idx[a]][idx[b]];
			diag = FastMath.max(diag, sub[a][a]);
		}
		// A touch of damping keeps nearly identical peaks solvable
		for (int a = 0; a < k; a++)
			sub[a][a] += 1e-10 * diag + 1e-300;
		final RealMatrix mat = new Array2DRowRealMatrix(sub, false);
		final double[] sol = new CholeskyDecomposition(mat, 1e-14, 1e-300)
			.getSolver().solve(new org.apache.commons.math3.linear.ArrayRealVector(
				rhs, false)).toArray();
		final double[] s = new double[n];
		for (int a = 0; a < k; a++)
			s[idx[a]] = sol[a];
		return s;
	}
}
