
package gausscurvefit;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.analysis.differentiation.UnivariateDifferentiableFunction;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.exception.NumberIsTooSmallException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.fitting.AbstractCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.ParameterValidator;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;

class GaussianArrayCurveFitter extends AbstractCurveFitter {

	/** Parametric function to be fitted. */
	private static final GaussianArray.Parametric FUNCTION =
		new GaussianArray.Parametric();

	/** Initial guess. */
	private final double[] initialGuess;
	/** Maximum number of iterations of the optimization algorithm. */
	private final int maxIter;
	private final double peakTol;
	private final int deg;

	/**
	 * Constructor used by the factory methods.
	 *
	 * @param initialGuess Initial guess. If set to {@code null}, the initial
	 *          guess will be estimated using the {@link ParameterGuesser}.
	 * @param maxIter Maximum number of iterations of the optimization algorithm.
	 */
	private GaussianArrayCurveFitter(final double[] initialGuess,
		final int maxIter, final double peakTol, final int deg)
	{
		this.initialGuess = initialGuess;
		this.peakTol = peakTol;
		this.maxIter = maxIter;
		this.deg = deg;
	}

	/**
	 * Creates a default curve fitter. The initial guess for the parameters will
	 * be {@link ParameterGuesser} computed automatically, and the maximum number
	 * of iterations of the optimization algorithm is set to
	 * {@link Integer#MAX_VALUE}.
	 *
	 * @return a curve fitter.
	 * @see #withStartPoint(double[])
	 * @see #withMaxIterations(int)
	 */
	static GaussianArrayCurveFitter create(final double peakTol,
		final int deg)
	{
		return new GaussianArrayCurveFitter(null, Integer.MAX_VALUE, peakTol, deg);
	}

	/**
	 * Configure the start point (initial guess).
	 *
	 * @param newStart new start point (initial guess)
	 * @return a new instance.
	 */
	GaussianArrayCurveFitter withStartPoint(final double[] newStart) {
		return new GaussianArrayCurveFitter(newStart.clone(), maxIter, peakTol,
			deg);
	}

	/**
	 * Configure the maximum number of iterations.
	 *
	 * @param newMaxIter maximum number of iterations
	 * @return a new instance.
	 */
	public GaussianArrayCurveFitter withMaxIterations(final int newMaxIter) {
		return new GaussianArrayCurveFitter(initialGuess, newMaxIter, peakTol, deg);
	}

	/** {@inheritDoc} */
	@Override
	protected LeastSquaresProblem getProblem(
		final Collection<WeightedObservedPoint> observations)
	{

		// Prepare least-squares problem.
		final int len = observations.size();
		final double[] xx = new double[len];
		final double[] target = new double[len];
		final double[] weights = new double[len];

		int i = 0;
		for (final WeightedObservedPoint obs : observations) {
			xx[i] = obs.getX();
			target[i] = obs.getY();
			weights[i] = obs.getWeight();
			++i;
		}

		final AbstractCurveFitter.TheoreticalValuesFunction model =
			new AbstractCurveFitter.TheoreticalValuesFunction(FUNCTION, observations);

		final double[] startPoint = initialGuess != null ? initialGuess
			: new ParameterGuesser(observations, peakTol, deg).guess();

		final GaussianArrayParameterValidator parValid =
			new GaussianArrayParameterValidator(startPoint, xx, target);

		return new LeastSquaresBuilder().parameterValidator(parValid)
			.maxEvaluations(Integer.MAX_VALUE).maxIterations(maxIter).lazyEvaluation(
				false).start(startPoint).target(target).weight(new DiagonalMatrix(
					weights)).model(model.getModelFunction(), model
						.getModelFunctionJacobian()).build();
	}

	/** {@inheritDoc} */
	@Override
	protected LeastSquaresOptimizer getOptimizer() {
		return new LevenbergMarquardtOptimizer().withCostRelativeTolerance(1e-12)
			.withOrthoTolerance(1e-12).withParameterRelativeTolerance(1e-10);
	}

	/**
	 * Guesses the parameters {@code norm}, {@code mean}, and {@code sigma} of a
	 * {@link GaussianArray.Parametric} based on the specified observed points.
	 */
	static class ParameterGuesser {

		/** Normalization factor. */
		private RealVector norm = new ArrayRealVector();
		/** Mean. */
		private RealVector mean = new ArrayRealVector();
		/** Standard deviation. */
		private RealVector sigma = new ArrayRealVector();

		/**
		 * Polynomial background. First element is {@code deg}, followed by all
		 * coefficients.
		 */
		private RealVector poly = new ArrayRealVector();

		/**
		 * Constructs instance with the specified observed points.
		 *
		 * @param observations Observed points from which to guess the parameters of
		 *          the train of Gaussians.
		 * @throws NullArgumentException if {@code observations} is {@code null}.
		 * @throws NumberIsTooSmallException if there are less than 3 observations.
		 */
		ParameterGuesser(
			final Collection<WeightedObservedPoint> observations,
			final double peakTol, final int deg)
		{
			if (observations == null) {
				throw new NullArgumentException(LocalizedFormats.INPUT_ARRAY);
			}
			if (observations.size() < 3) {
				throw new NumberIsTooSmallException(observations.size(), 3, true);
			}

			final List<WeightedObservedPoint> sorted = sortObservations(observations);
			final RealVector params = basicGuess(sorted.toArray(
				new WeightedObservedPoint[0]), peakTol, deg);

			int gaussStart = 1;
			poly = poly.append(params.getEntry(0));
			if (deg != -1) {
				gaussStart = deg + 2;
				poly = poly.append(params.getSubVector(1, deg + 1));
			}
			for (int p = gaussStart; p < params.getDimension(); p += 3) {
				norm = norm.append(params.getEntry(p));
				mean = mean.append(params.getEntry(p + 1));
				sigma = sigma.append(params.getEntry(p + 2));
			}
		}

		/**
		 * Gets an estimation of the parameters.
		 *
		 * @return the guessed parameters, in the following order:
		 *         <ul>
		 *         <li>Normalization factor</li>
		 *         <li>Mean</li>
		 *         <li>Standard deviation</li>
		 *         </ul>
		 */
		double[] guess() {
			RealVector guess = new ArrayRealVector();
			guess = guess.append(poly);
			for (int vv = 0; vv < norm.getDimension(); vv++) {
				guess = guess.append(norm.getEntry(vv)).append(mean.getEntry(vv))
					.append(sigma.getEntry(vv));
			}
			return guess.toArray();
		}

		/**
		 * Sort the observations.
		 *
		 * @param unsorted Input observations.
		 * @return the input observations, sorted.
		 */
		private List<WeightedObservedPoint> sortObservations(
			final Collection<WeightedObservedPoint> unsorted)
		{
			final List<WeightedObservedPoint> observations = new ArrayList<>(
				unsorted);

			final Comparator<WeightedObservedPoint> cmp =
				new Comparator<WeightedObservedPoint>()
			{

					/** {@inheritDoc} */
					@Override
					public int compare(final WeightedObservedPoint p1,
						final WeightedObservedPoint p2)
				{
						if (p1 == null && p2 == null) {
							return 0;
						}
						if (p1 == null) {
							return -1;
						}
						if (p2 == null) {
							return 1;
						}
						int comp = Double.compare(p1.getX(), p2.getX());
						if (comp != 0) {
							return comp;
						}
						comp = Double.compare(p1.getY(), p2.getY());
						if (comp != 0) {
							return comp;
						}
						comp = Double.compare(p1.getWeight(), p2.getWeight());
						if (comp != 0) {
							return comp;
						}
						return 0;
					}
				};
			Collections.sort(observations, cmp);
			return observations;
		}

		/**
		 * Guesses the parameters based on the specified observed points.
		 *
		 * @param points Observed points, sorted.
		 * @return the guessed parameters (normalization factor, mean and sigma).
		 */
		private RealVector basicGuess(final WeightedObservedPoint[] points,
			final double tolpk, final int deg)
		{
			RealVector xvals = new ArrayRealVector();
			RealVector yvals = new ArrayRealVector();

			for (int o = 0; o < points.length; o++) {
				xvals = xvals.append(points[o].getX());
				yvals = yvals.append(points[o].getY());
			}
			final double minY = yvals.getMinValue();
			RealVector polyGuess = new ArrayRealVector();
			RealVector gaussGuess = new ArrayRealVector();
			polyGuess = polyGuess.append(deg);
			if (deg >= 0) polyGuess = polyGuess.append(minY / 2).append(
				new ArrayRealVector(new double[deg]));

			// Local maxima where the peaks are
			final int[] maximaIdx = findMaxima(points, tolpk, true);

			// Define a Gaussian at each peak
			final double[] means = new double[maximaIdx.length];
			final double[] sds = new double[maximaIdx.length];
			final double[] norms = new double[maximaIdx.length];

			for (int m = 0; m < maximaIdx.length; m++) {
				final double yMin = yvals.getMinValue();
				norms[m] = yvals.getEntry(maximaIdx[m]) - yMin;
				means[m] = xvals.getEntry(maximaIdx[m]);

				// estimate sds from maxima and mean
				boolean foundFWHM = false; // Full width at half maximum
				boolean foundRWHM = false;
				boolean foundLWHM = false;
				double LWHM = 0.0;
				double RWHM = 0.0;
				double FWHM = (xvals.getMaxValue() - xvals.getMinValue()) / 2;

				final double yRange = yvals.getEntry(maximaIdx[m]) - yMin;
				double hm = yRange / 2.0; // Actual profile value; incledes offset
				final int pkPos = maximaIdx[m];

				final int range = 3;
				final double inc = 1.05;
				double peakDistance = 0.0;

				int count = 0;
				for (int p = m - range; p < m + range; p++) {
					if (p >= 0 && p < maximaIdx.length - 1) {
						peakDistance += xvals.getEntry(maximaIdx[p + 1]) - xvals.getEntry(
							maximaIdx[p]);
						count++;
					}
				}
				peakDistance = peakDistance / count;

				while (!(foundFWHM && FWHM < peakDistance) && hm * inc < 0.9 * yRange) {
					foundFWHM = false;
					foundRWHM = false;
					foundLWHM = false;
					// Right side, check 3 consecutive points for smoothing
					int p = 0;
					while (pkPos + p + 2 < xvals.getDimension() && !foundRWHM) {
						if (yvals.getEntry(pkPos + p) < (hm + yMin) && yvals.getEntry(
							pkPos + p + 1) < (hm + yMin) && yvals.getEntry(pkPos + p +
								2) < (hm + yMin))
						{
							foundRWHM = true;
							RWHM = xvals.getEntry(pkPos + p) - means[m];
						}
						p++;
					}

					// Left side, check 3 consecutive points for smoothing
					p = 0;
					while (pkPos - p - 2 >= 0 && !foundLWHM) {
						if (yvals.getEntry(pkPos - p) < (hm + yMin) && yvals.getEntry(
							pkPos - p - 1) < (hm + yMin) && yvals.getEntry(pkPos - p -
								2) < (hm + yMin))
						{
							foundLWHM = true;
							LWHM = means[m] - xvals.getEntry(pkPos - p);
						}
						p++;
					}

					if (foundLWHM && foundRWHM) {
						foundFWHM = true;
						FWHM = 2 * FastMath.min(LWHM, RWHM);
					}

					// Do another round with larger hm
					hm = hm * inc;
				}
//				System.out.println("PeakDist: " + peakDistance + " (" + LWHM + ":" +
//					RWHM + "); " + hm / yRange);
				sds[m] = FWHM / (2 * FastMath.sqrt(2 * FastMath.log(2)));
				gaussGuess = gaussGuess.append(norms[m]);
				gaussGuess = gaussGuess.append(means[m]);
				gaussGuess = gaussGuess.append(sds[m]);
			}
			// Array of parameters for Fitter
			return polyGuess.append(gaussGuess);
		}

		/**
		 * Adapted From: Calculates peak positions of 1D array N.Vischer,
		 * 13-sep-2013
		 *
		 * @param x Array containing peaks.
		 * @param tolerance Depth of a qualified valley must exceed tolerance.
		 *          Tolerance must be >= 0. Flat tops are marked at their centers.
		 * @param includeEnds If 'false', a peak is only accepted if it is separated
		 *          by two qualified valleys. If 'true', a peak is also accepted if
		 *          separated by one qualified valley and by a border.
		 * @return Positions of peaks, sorted with decreasing amplitude
		 */
		private int[] findMaxima(final WeightedObservedPoint[] x, double tolerance,
			final boolean includeEnds)
		{
			final int len = x.length;
			if (len < 2) return new int[0];
			if (tolerance < 0) tolerance = 0;
			int[] maxPositions = new int[len];
			double max = x[0].getY();
			double min = x[0].getY();
			int maxPos = 0;
			int lastMaxPos = -1;
			boolean leftValleyFound = includeEnds;
			int maxCount = 0;
			for (int j = 1; j < len; j++) {
				final double val = x[j].getY();
				if (val > min + tolerance) leftValleyFound = true;
				if (val > max && leftValleyFound) {
					max = val;
					maxPos = j;
				}
				if (leftValleyFound) lastMaxPos = maxPos;
				if (val < max - tolerance && leftValleyFound) {
					maxPositions[maxCount] = maxPos;
					maxCount++;
					leftValleyFound = false;
					min = val;
					max = val;
				}
				if (val < min) {
					min = val;
					if (!leftValleyFound) max = val;
				}
			}
			if (includeEnds) {
				if (maxCount > 0 && maxPositions[maxCount - 1] != lastMaxPos)
					maxPositions[maxCount++] = lastMaxPos;
				if (maxCount == 0 && max - min >= tolerance) maxPositions[maxCount++] =
					lastMaxPos;
			}
			final int[] cropped = new int[maxCount];
			System.arraycopy(maxPositions, 0, cropped, 0, maxCount);
			maxPositions = cropped;
			final double[] maxValues = new double[maxCount];
			for (int m = 0; m < maxCount; m++) {
				int pos = maxPositions[m];
				double midPos = pos;
				while (pos < len - 1 && x[pos] == x[pos + 1]) {
					midPos += 0.5;
					pos++;
				}
				maxPositions[m] = (int) midPos;
				maxValues[m] = x[maxPositions[m]].getY();
			}
			// int[] rankPositions = Tools.rank(maxValues);
			// int[] returnArr = new int[maxCount];
			// for (int n = 0; n < maxCount; n++) {
			// int pos = maxPositions[rankPositions[n]];
			// returnArr[maxCount - n - 1] = pos; //in descending order
			// }
			// return returnArr;
			Arrays.sort(maxPositions);
			return maxPositions;
		}
	}

	/**
	 * Checks the polynomial coefficients {@code poly} and the parameter triplets
	 * {@code norm}, {@code mean}, and {@code sigma} of a
	 * {@code GaussianArray.Parametric} for applied constraints, based on the
	 * specified points and initial guess.
	 */
	private static class GaussianArrayParameterValidator implements
		ParameterValidator
	{

		private final double[] initialParameterSet;
		private final double maxX;
		private final double minX;
		private final double minY;
		private final double minSD;
		private final double maxMeanDiff;
		private final int gNumber;

		private GaussianArrayParameterValidator(final double[] initialGuess,
			final double[] xtarget, final double[] ytarget)
		{
			this.gNumber = (initialGuess.length - (int) initialGuess[0] - 2) / 3; // #
																																						// //
																																						// peaks
			this.initialParameterSet = initialGuess;
			this.minX = new ArrayRealVector(xtarget).getMinValue();
			this.maxX = new ArrayRealVector(xtarget).getMaxValue();
			this.minY = new ArrayRealVector(ytarget).getMinValue();
			this.maxMeanDiff = FastMath.min((maxX - minX) / gNumber, (maxX - minX) /
				500);
			this.minSD = 2.0;
		}

		private double peakDistance(final int i) {
			final int m = i - 1; // means positioned right before SDs
			final int range = 3;
			double peakDistance = 0.0;

			int count = 0;
			for (int p = m - 3 * range; p < m + 3 * range; p += 3) {
				if (p >= (int) initialParameterSet[0] + 2 &&
					p < initialParameterSet.length - 4)
				{
					peakDistance += initialParameterSet[p + 3] - initialParameterSet[p];
					count++;
				}
			}
			return peakDistance / count;
		}

		@Override
		// Set parameter constraints here
		public RealVector validate(final RealVector param) {
			// Do not fit the order of the polynomial
			if (param.getEntry(0) != initialParameterSet[0]) {
				param.setEntry(0, initialParameterSet[0]);
			}

			int gaussStart = 1;
			// Polynomyal parameters
			if (gaussStart != -1) {
				gaussStart = (int) initialParameterSet[0] + 2;
				if (param.getEntry(1) > 1.2 * minY) param.setEntry(1, 1.2 * minY);
				if (param.getEntry(1) < 0.0) param.setEntry(1, 0.0);
			}

			// Gaussian Parameters
			for (int i = gaussStart; i < param.getDimension(); i++) {
				// {0=Norm, 1=Mean, 2=Sigma}
				// Lower bound of 0 for all 3
				if ((i - gaussStart) % 3 == 0 && param.getEntry(i) < 0.0) // Norm
					param.setEntry(i, 0.0);
				if ((i - gaussStart) % 3 == 1 && param.getEntry(i) < 1e-1) // Mean
					param.setEntry(i, 1e-1);
				if ((i - gaussStart) % 3 == 2 && param.getEntry(i) < minSD) // SD
					param.setEntry(i, minSD);
				// Upper bounds for each parameter
				if (((i - gaussStart) % 3 == 0) && param.getEntry(
					i) > initialParameterSet[i]) param.setEntry(i,
						initialParameterSet[i]);
				if ((i - gaussStart) % 3 == 1) {
					if (param.getEntry(i) > maxX) param.setEntry(i, maxX);
					// Keep means close to maxima
					final double diff = param.getEntry(i) - initialParameterSet[i];
					final double sign = diff / Math.abs(diff);
					if (Math.abs(diff) > maxMeanDiff) param.setEntry(i,
						initialParameterSet[i] + sign * maxMeanDiff);
				}
				if ((i - gaussStart) % 3 == 2) {
					final double maxSD = peakDistance(i);
					if (param.getEntry(i) > maxSD) param.setEntry(i, maxSD);
				}
			}
			return param;
		}

	}
}

class GaussianArray implements UnivariateDifferentiableFunction {
	// Implements a train of Gaussian peaks
	// with an optional polynomial background function

	private final RealVector means;
	private final RealVector norms;
	private final RealVector sds;

	private final RealVector poly; // Polynomial background curve

	public GaussianArray(final RealVector norms, final RealVector means,
		final RealVector sds, final RealVector poly)
		throws NotStrictlyPositiveException
	{
		this.means = means;
		this.norms = norms;
		this.sds = sds;
		this.poly = poly;
	}

	@Override
	public double value(final double x) {
		double output = 0;
		final int degPoly = (int) poly.getEntry(0);
		for (int i = 0; i < norms.getDimension(); i++) {
			output += new Gaussian(norms.getEntry(i), means.getEntry(i), sds.getEntry(
				i)).value(x);
		}
		output += new PolynomialFunction(poly.getSubVector(1, degPoly + 1)
			.toArray()).value(x);
		return output;
	}

	@Override
	public DerivativeStructure value(final DerivativeStructure t)
		throws DimensionMismatchException
	{
		// TODO Write Derivative Structure for future use
		// final double[] u = new double[]
		// {is.multiplyToSelf(means.subtractToself(t.getValue()).multiplyToSelf(-1)).toArray()};

		// final double[] f = new double[t.getOrder() + 1];

		// the nth order derivative of the Gaussian has the form:
		// dn(g(x)/dxn = (norm / s^n) P_n(u) exp(-u^2/2) with u=(x-m)/s
		// where P_n(u) is a degree n polynomial with same parity as n
		// P_0(u) = 1, P_1(u) = -u, P_2(u) = u^2 - 1, P_3(u) = -u^3 + 3 u...
		// the general recurrence relation for P_n is:
		// P_n(u) = P_(n-1)'(u) - u P_(n-1)(u)
		// as per polynomial parity, we can store coefficients of both P_(n-1)
		// and
		// P_n in the same array
		final PolynomialFunction p = new PolynomialFunction(poly.toArray());
		return p.value(t);
	}

	/**
	 * Parametric function where the input array contains the parameters of the
	 * Gaussian, ordered as follows:
	 * <ul>
	 * <li>Norm</li>
	 * <li>Mean</li>
	 * <li>Standard deviation</li>
	 * </ul>
	 */
	public static class Parametric implements ParametricUnivariateFunction {

		/**
		 * Computes the value of the Gaussian Array at {@code x}.
		 *
		 * @param x Value for which the function must be computed.
		 * @param param Values of norm, mean and standard deviation, passed in
		 *          subsequent triplets.
		 * @return the value of the function.
		 * @throws NullArgumentException if {@code param} is {@code null}.
		 * @throws DimensionMismatchException if the size of {@code param} is not a
		 *           multiple of 3 + number of polynomial coefficients + 1.
		 * @throws NotStrictlyPositiveException if {@code param[2]} is negative.
		 */
		@Override
		public double value(final double x, final double... param) {
			// validateParameters(param);
			final RealVector allParam = new ArrayRealVector(param);
			final int degPoly = (int) param[0];
			RealVector means = new ArrayRealVector();
			RealVector norms = new ArrayRealVector();
			RealVector sds = new ArrayRealVector();
			RealVector poly = new ArrayRealVector();

			int gaussStart = 1;
			if (degPoly != -1) { // No background, param[0] == -1
				poly = allParam.getSubVector(0, degPoly + 2);
				gaussStart = degPoly + 2;
			}

			for (int i = gaussStart; i < param.length; i += 3) {
				means = means.append(param[i + 1]);
				norms = norms.append(param[i]);
				sds = sds.append(param[i + 2]);
			}
			return new GaussianArray(norms, means, sds, poly).value(x);
		}

		/**
		 * Computes the value of the gradient at {@code x}. The components of the
		 * gradient vector are the partial derivatives of the function with respect
		 * to each of the <em>parameters</em> (norm, mean and standard deviation,
		 * for each Gaussian in the array).
		 *
		 * @param x Value at which the gradient must be computed.
		 * @param param Values of poly coeffs, norm, mean and standard deviation
		 *          triplets.
		 * @return the gradient vector at {@code x}.
		 * @throws NullArgumentException if {@code param} is {@code null}.
		 * @throws DimensionMismatchException if the size of {@code param} minus the
		 *           number of polyomial coefficients, minus 1, is not divisible by
		 *           3.
		 * @throws NotStrictlyPositiveException if any gaussian SD is negative.
		 */
		@Override
		public double[] gradient(final double x, final double... param)
			throws NullArgumentException, DimensionMismatchException,
			NotStrictlyPositiveException
		{
			// validateParameters(param);

			final int degPoly = (int) param[0];
			RealVector out = new ArrayRealVector();
			RealVector poly = new ArrayRealVector();

			final int gaussStart = degPoly + 2;
			out = out.append(1.0);
			for (int pp = 1; pp < gaussStart; pp++)
				poly = poly.append(param[pp]);
			poly = new ArrayRealVector(new PolynomialFunction.Parametric().gradient(x,
				poly.toArray()));
			out = out.append(poly);

			for (int gg = gaussStart; gg < param.length; gg += 3) {
				out = out.append(new ArrayRealVector(new Gaussian.Parametric().gradient(
					x, param[gg], param[gg + 1], param[gg + 2])));
			}
			return out.toArray();
		}
	}
}
