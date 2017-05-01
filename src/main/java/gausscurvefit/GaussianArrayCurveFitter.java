package gausscurvefit;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.math3.analysis.FunctionUtils;
import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.analysis.UnivariateFunction;
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

public class GaussianArrayCurveFitter extends AbstractCurveFitter {

	/** Parametric function to be fitted. */
//	ParametricUnivariateFunction a = new Gaussian.Parametric();
//	ParametricUnivariateFunction b = new Gaussian.Parametric();
//	private final ParametricUnivariateFunction FUNCTION = FunctionUtils.add(a, b);
		private static final GaussianArray.Parametric FUNCTION = new GaussianArray.Parametric();

//		/** {@inheritDoc} */
//		@Override
//		public double value(final double x, final double... p) {
//			double v = Double.POSITIVE_INFINITY;
//			try {
//				v = super.value(x, p);
//			} catch (final NotStrictlyPositiveException e) { // NOPMD
//				// Do nothing.
//			}
//			return v;
//		}
//
//		/** {@inheritDoc} */
//		@Override
//		public double[] gradient(final double x, final double... p) {
//			RealVector v = new ArrayRealVector();
//			for (int i = 0; i < p.length; i++) {
//				v = v.append(Double.POSITIVE_INFINITY);
//			}
//			try {
//				v = new ArrayRealVector(super.gradient(x, p));
//			} catch (final NotStrictlyPositiveException e) { // NOPMD
//				System.out.println(e.getMessage());
//			}
//			return v.toArray();
//		}

	/** Initial guess. */
	private final double[] initialGuess;
	/** Maximum number of iterations of the optimization algorithm. */
	private final int maxIter;
	private final double peakTol;
	private final int deg;

	/**
	 * Contructor used by the factory methods.
	 *
	 * @param initialGuess
	 *            Initial guess. If set to {@code null}, the initial guess will
	 *            be estimated using the {@link ParameterGuesser}.
	 * @param maxIter
	 *            Maximum number of iterations of the optimization algorithm.
	 */
	private GaussianArrayCurveFitter(final double[] initialGuess, final int maxIter,
	        final double peakTol, final int deg) {
		this.initialGuess = initialGuess;
		this.peakTol = peakTol;
		this.maxIter = maxIter;
		this.deg = deg;
	}

	/**
	 * Creates a default curve fitter. The initial guess for the parameters will
	 * be {@link ParameterGuesser} computed automatically, and the maximum
	 * number of iterations of the optimization algorithm is set to
	 * {@link Integer#MAX_VALUE}.
	 *
	 * @return a curve fitter.
	 * @see #withStartPoint(double[])
	 * @see #withMaxIterations(int)
	 */
	public static GaussianArrayCurveFitter create(final double peakTol, final int deg) {
		return new GaussianArrayCurveFitter(null, Integer.MAX_VALUE, peakTol, deg);
	}

	/**
	 * Configure the start point (initial guess).
	 *
	 * @param newStart
	 *            new start point (initial guess)
	 * @return a new instance.
	 */
	public GaussianArrayCurveFitter withStartPoint(final double[] newStart) {
		return new GaussianArrayCurveFitter(newStart.clone(), maxIter, peakTol, deg);
	}

	/**
	 * Configure the maximum number of iterations.
	 *
	 * @param newMaxIter
	 *            maximum number of iterations
	 * @return a new instance.
	 */
	public GaussianArrayCurveFitter withMaxIterations(final int newMaxIter) {
		return new GaussianArrayCurveFitter(initialGuess, newMaxIter, peakTol, deg);
	}

	/** {@inheritDoc} */
	@Override
	protected LeastSquaresProblem getProblem(final Collection<WeightedObservedPoint> observations) {

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

		final AbstractCurveFitter.TheoreticalValuesFunction model = new AbstractCurveFitter.TheoreticalValuesFunction(
		        FUNCTION, observations);

		final double[] startPoint = initialGuess != null ? initialGuess
		        : new ParameterGuesser(observations, peakTol, deg).guess();

		final GaussianArrayParameterValidator parValid = new GaussianArrayParameterValidator(
		        startPoint, xx, target);

		// Return a new least squares problem set up to fit a Gaussian curve to
		// the observed points
		return new LeastSquaresBuilder().parameterValidator(parValid)
		        .maxEvaluations(Integer.MAX_VALUE).maxIterations(maxIter).lazyEvaluation(false)
		        .start(startPoint).target(target).weight(new DiagonalMatrix(weights))
		        .model(model.getModelFunction(), model.getModelFunctionJacobian()).build();
	}

	/** {@inheritDoc} */
	@Override
	protected LeastSquaresOptimizer getOptimizer() {
		// Return a new least squares problem set up to fit a univariate curve
		// to the observed points.
		return new LevenbergMarquardtOptimizer().withCostRelativeTolerance(1e-12)
		        .withOrthoTolerance(1e-12).withParameterRelativeTolerance(1e-10);
	}

	/**
	 * Guesses the parameters {@code norm}, {@code mean}, and {@code sigma} of a
	 * {@link GaussianArrayBG.Parametric} based on the specified observed
	 * points.
	 */
	public static class ParameterGuesser {

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
		 * @param observations
		 *            Observed points from which to guess the parameters of the
		 *            train of Gaussians.
		 * @throws NullArgumentException
		 *             if {@code observations} is {@code null}.
		 * @throws NumberIsTooSmallException
		 *             if there are less than 3 observations.
		 */
		public ParameterGuesser(final Collection<WeightedObservedPoint> observations,
		        final double peakTol, final int deg) {
			if (observations == null) {
				throw new NullArgumentException(LocalizedFormats.INPUT_ARRAY);
			}
			if (observations.size() < 3) {
				throw new NumberIsTooSmallException(observations.size(), 3, true);
			}

			final List<WeightedObservedPoint> sorted = sortObservations(observations);
			final RealVector params = basicGuess(sorted.toArray(new WeightedObservedPoint[0]),
			        peakTol, deg);

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
		public double[] guess() {
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
		 * @param unsorted
		 *            Input observations.
		 * @return the input observations, sorted.
		 */
		private List<WeightedObservedPoint> sortObservations(
		        final Collection<WeightedObservedPoint> unsorted) {
			final List<WeightedObservedPoint> observations = new ArrayList<>(unsorted);

			final Comparator<WeightedObservedPoint> cmp = new Comparator<WeightedObservedPoint>() {

				/** {@inheritDoc} */
				@Override
				public int compare(final WeightedObservedPoint p1, final WeightedObservedPoint p2) {
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
		 * @param points
		 *            Observed points, sorted.
		 * @return the guessed parameters (normalization factor, mean and
		 *         sigma).
		 */
		private RealVector basicGuess(final WeightedObservedPoint[] points, final double tolpk,
		        final int deg) {
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
			if (deg >= 0)
				polyGuess = polyGuess.append(minY / 2)
			        .append(new ArrayRealVector(new double[deg]));
			
			// Local maxima where the peaks are
			final int[] maximaIdx = findMaxima(points, tolpk, true);

			// Define a Gaussian at each peak
			final double[] means = new double[maximaIdx.length];
			final double[] sds = new double[maximaIdx.length];
			final double[] norms = new double[maximaIdx.length];

			for (int m = 0; m < maximaIdx.length; m++) {
				norms[m] = yvals.getEntry(maximaIdx[m]);
				means[m] = xvals.getEntry(maximaIdx[m]);

				// estimate sds from maxima and mean
				boolean foundHWHM = false;
				double leftHWHM = 0.0;
				double rightHWHM = 0.0;
				final double yMin = yvals.getMinValue();
				final double maxSD = xvals.getEntry(xvals.getMaxIndex()) / maximaIdx.length;
				final double yRange = norms[m] - yMin;
				double hm = yRange / 2.0 + yMin;
				final int pkPos = maximaIdx[m];
				int p = 0;

				while (!foundHWHM) {
					if (pkPos + p + 2 > xvals.getMaxIndex()) {
						rightHWHM = xvals.getEntry(xvals.getMaxIndex()) - means[m];
						foundHWHM = true;
						// Right side check 3 consecutive points for smoothing
					} else if (yvals.getEntry(pkPos + p) < hm && yvals.getEntry(pkPos + p + 1) < hm
					        && yvals.getEntry(pkPos + p + 2) < hm) {
						foundHWHM = true;
						rightHWHM = xvals.getEntry(pkPos + p) - means[m] > 0.0
						        ? xvals.getEntry(pkPos + p) - means[m] : 0.1;
					}
					if (pkPos - p - 2 < 0) {
						leftHWHM = xvals.getEntry(pkPos) - xvals.getEntry(0);
						foundHWHM = true;
						// Left side check 3 consecutive points for smoothing
					} else if (yvals.getEntry(pkPos - p) < hm && yvals.getEntry(pkPos - p - 1) < hm
					        && yvals.getEntry(pkPos - p - 2) < hm) {
						foundHWHM = true;
						leftHWHM = means[m] - xvals.getEntry(pkPos - p);
					}
					if (foundHWHM) {
						sds[m] = FastMath.max(rightHWHM, leftHWHM)
						        / (2.0 * FastMath.sqrt(2.0 * FastMath.log(2)));
						if (sds[m] <= 0.0 || sds[m] > maxSD) {
							if (hm < 0.95 * yRange) {
								hm *= 1.1;
								p = 0; // Do another round with smaller hm
							} else {
								// If exiting without finding anything
								sds[m] = maxSD;
								break;
							}
						}
					}
					p++;
				}
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
		 * @param x
		 *            Array containing peaks.
		 * @param tolerance
		 *            Depth of a qualified valley must exceed tolerance.
		 *            Tolerance must be >= 0. Flat tops are marked at their
		 *            centers.
		 * @param includeEnds
		 *            If 'false', a peak is only accepted if it is separated by
		 *            two qualified valleys. If 'true', a peak is also accepted
		 *            if separated by one qualified valley and by a border.
		 * @return Positions of peaks, sorted with decreasing amplitude
		 */
		private int[] findMaxima(final WeightedObservedPoint[] x, double tolerance,
		        final boolean includeEnds) {
			final int len = x.length;
			if (len < 2)
				return new int[0];
			if (tolerance < 0)
				tolerance = 0;
			int[] maxPositions = new int[len];
			double max = x[0].getY();
			double min = x[0].getY();
			int maxPos = 0;
			int lastMaxPos = -1;
			boolean leftValleyFound = includeEnds;
			int maxCount = 0;
			for (int j = 1; j < len; j++) {
				final double val = x[j].getY();
				if (val > min + tolerance)
					leftValleyFound = true;
				if (val > max && leftValleyFound) {
					max = val;
					maxPos = j;
				}
				if (leftValleyFound)
					lastMaxPos = maxPos;
				if (val < max - tolerance && leftValleyFound) {
					maxPositions[maxCount] = maxPos;
					maxCount++;
					leftValleyFound = false;
					min = val;
					max = val;
				}
				if (val < min) {
					min = val;
					if (!leftValleyFound)
						max = val;
				}
			}
			if (includeEnds) {
				if (maxCount > 0 && maxPositions[maxCount - 1] != lastMaxPos)
					maxPositions[maxCount++] = lastMaxPos;
				if (maxCount == 0 && max - min >= tolerance)
					maxPositions[maxCount++] = lastMaxPos;
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
	 * Checks the polynomial coefficients {@code poly} and the parameter
	 * triplets {@code norm}, {@code mean}, and {@code sigma} of a
	 * {@code GaussianArray.Parametric} for applied constraints, based on the
	 * specified points and initial guess.
	 */
	public static class GaussianArrayParameterValidator implements ParameterValidator {

		private final double[] initialParameterSet;
		private final double maxXvalue;
		private final double minYvalue;
		private final double maxMeandiff;
		private final int gNumber;

		public GaussianArrayParameterValidator(final double[] initialGuess, final double[] xtarget,
		        final double[] ytarget) {
			this.gNumber = (initialGuess.length - (int) initialGuess[0] - 2) / 3; // #
			                                                                      // peaks
			this.initialParameterSet = initialGuess;
			this.maxXvalue = xtarget[xtarget.length - 1];
			this.minYvalue = new ArrayRealVector(ytarget).getMinValue();
			this.maxMeandiff = maxXvalue / gNumber / 20;
		}

		@Override
		// Set parameter constraints here
		public RealVector validate(final RealVector param) {
			if (param.getEntry(0) != initialParameterSet[0]) {
				param.setEntry(0, initialParameterSet[0]);
			}

			int gaussStart = 1;
			// Polynomyal parameters
			if (gaussStart != -1) {
				gaussStart = (int) initialParameterSet[0] + 2;
				if (param.getEntry(1) > 1.2 * minYvalue)
					param.setEntry(1, 1.2 * minYvalue);
				if (param.getEntry(1) < 0.0)
					param.setEntry(1, 0.0);
			}

			// Gaussian Parameters
			for (int i = gaussStart; i < param.getDimension(); i++) {
				// {0=Norm, 1=Mean, 2=Sigma}
				// Lower bound of 0 for all 3
				if ((i - gaussStart) % 3 == 0 && param.getEntry(i) < 0.0) // Norm
					param.setEntry(i, 0.0);
				if (((i - gaussStart) % 3 == 1 || (i - gaussStart) % 3 == 2)
				        && param.getEntry(i) < 1e-1)
					param.setEntry(i, 1e-1);
				// Upper bounds for each parameter
				if (((i - gaussStart) % 3 == 0) && param.getEntry(i) > initialParameterSet[i])
					param.setEntry(i, initialParameterSet[i]);
				if ((i - gaussStart) % 3 == 1 && param.getEntry(i) > maxXvalue)
					param.setEntry(i, maxXvalue);
				if ((i - gaussStart) % 3 == 2) {
					double meanPD = 0.0;
					int nPDsamples = 0;
					// Average Peak Distance 2 neighbor peaks
					for (int pp = 0; pp < 9; pp += 3) {
						if (i - pp - 3 > gaussStart) {
							meanPD += FastMath
							        .abs(param.getEntry(i - pp - 1) - param.getEntry(i - pp - 4));
							nPDsamples++;
						}
						if (i + pp + 3 < param.getMaxIndex()) {
							meanPD += FastMath
							        .abs(param.getEntry(i + pp - 1) - param.getEntry(i + pp + 2));
							nPDsamples++;
						}
						if (meanPD < 0.0)
							System.out.println("meanPD: " + meanPD);
					}
					meanPD = meanPD / nPDsamples / (2 * FastMath.sqrt(2 * FastMath.log(2)));
					if (meanPD == 0.0)
						meanPD = maxXvalue / gNumber / (2 * FastMath.sqrt(2 * FastMath.log(2)));
					if (param.getEntry(i) > meanPD)
						param.setEntry(i, meanPD);
				}
				// Keep means close to maxima
				final double diff = param.getEntry(i) - initialParameterSet[i];
				if ((i - gaussStart) % 3 == 1 && (Math.abs(diff) > maxMeandiff))
					param.setEntry(i, initialParameterSet[i] + maxMeandiff * diff / Math.abs(diff));
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

	public GaussianArray(final RealVector norms, final RealVector means, final RealVector sds,
	        final RealVector poly) throws NotStrictlyPositiveException {
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
			output += new Gaussian(norms.getEntry(i), means.getEntry(i), sds.getEntry(i)).value(x);
		}
		output += new PolynomialFunction(poly.getSubVector(1, degPoly + 1).toArray()).value(x);
		return output;
	}


	@Override
	public DerivativeStructure value(final DerivativeStructure t)
	        throws DimensionMismatchException {
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
		PolynomialFunction p = new PolynomialFunction(poly.toArray());
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
		 * Computes the value of the Gaussian at {@code x}.
		 *
		 * @param x
		 *            Value for which the function must be computed.
		 * @param param
		 *            Values of norm, mean and standard deviation, passed in
		 *            subsequent triplets.
		 * @return the value of the function.
		 * @throws NullArgumentException
		 *             if {@code param} is {@code null}.
		 * @throws DimensionMismatchException
		 *             if the size of {@code param} is not a multiple of 3 +
		 *             number of polynomial coefficients + 1.
		 * @throws NotStrictlyPositiveException
		 *             if {@code param[2]} is negative.
		 */
		@Override
		public double value(final double x, final double... param) {
			//validateParameters(param);
			RealVector allParam = new ArrayRealVector(param);
			final int degPoly = (int) param[0];
			RealVector means = new ArrayRealVector();
			RealVector norms = new ArrayRealVector();
			RealVector sds = new ArrayRealVector();
			RealVector poly = new ArrayRealVector();

			int gaussStart = 1;
			if (degPoly != -1) { // No background, param[0] == -1
				poly = allParam.getSubVector(0, degPoly+2);
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
		 * Computes the value of the gradient at {@code x}. The components of
		 * the gradient vector are the partial derivatives of the function with
		 * respect to each of the <em>parameters</em> (norm, mean and standard
		 * deviation, for each Gaussian in the array).
		 *
		 * @param x
		 *            Value at which the gradient must be computed.
		 * @param param
		 *            Values of poly coeffs, norm, mean and standard deviation
		 *            triplets.
		 * @return the gradient vector at {@code x}.
		 * @throws NullArgumentException
		 *             if {@code param} is {@code null}.
		 * @throws DimensionMismatchException
		 *             if the size of {@code param} minus the number of
		 *             polyomial coefficients, minus 1, is not divisible by 3.
		 * @throws NotStrictlyPositiveException
		 *             if any gaussian SD is negative.
		 */
		@Override
		public double[] gradient(final double x, final double... param)
		        throws NullArgumentException, DimensionMismatchException,
		        NotStrictlyPositiveException {
			//validateParameters(param);

			final int degPoly = (int) param[0];
			RealVector out = new ArrayRealVector();
			RealVector poly = new ArrayRealVector();

			final int gaussStart = degPoly + 2;
			out = out.append(1.0);
			for (int pp = 1; pp < gaussStart; pp++)
				poly = poly.append(param[pp]);
			poly = new ArrayRealVector(
			        new PolynomialFunction.Parametric().gradient(x, poly.toArray()));
			out = out.append(poly);

			for (int gg = gaussStart; gg < param.length; gg += 3) {
				out = out.append(new ArrayRealVector(new Gaussian.Parametric().gradient(x,
				        param[gg], param[gg + 1], param[gg + 2])));
			}
			return out.toArray();
		}
	}
}
