import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.math4.analysis.ParametricUnivariateFunction;
import org.apache.commons.math4.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math4.analysis.differentiation.UnivariateDifferentiableFunction;
import org.apache.commons.math4.analysis.function.Gaussian;
import org.apache.commons.math4.exception.DimensionMismatchException;
import org.apache.commons.math4.exception.NotStrictlyPositiveException;
import org.apache.commons.math4.exception.NullArgumentException;
import org.apache.commons.math4.exception.NumberIsTooSmallException;
import org.apache.commons.math4.exception.OutOfRangeException;
import org.apache.commons.math4.exception.ZeroException;
import org.apache.commons.math4.exception.util.LocalizedFormats;
import org.apache.commons.math4.fitting.AbstractCurveFitter;
import org.apache.commons.math4.fitting.WeightedObservedPoint;
import org.apache.commons.math4.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math4.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math4.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math4.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math4.fitting.leastsquares.ParameterValidator;
import org.apache.commons.math4.linear.ArrayRealVector;
import org.apache.commons.math4.linear.DiagonalMatrix;
import org.apache.commons.math4.linear.RealVector;
import org.apache.commons.math4.util.FastMath;



public class GaussianArrayCurveFitter extends AbstractCurveFitter {
	/** Parametric function to be fitted. */
	private static final GaussianArray.Parametric FUNCTION = new GaussianArray.Parametric() {
		/** {@inheritDoc} */
		@Override
		public double value(double x, double ... p) {
			double v = Double.POSITIVE_INFINITY;
			try {
				v = super.value(x, p);
			} catch (NotStrictlyPositiveException e) { // NOPMD
				// Do nothing.
			}
			return v;
		}

		/** {@inheritDoc} */
		@Override
		public double[] gradient(double x, double ... p) {
			double[] v = { Double.POSITIVE_INFINITY,
					Double.POSITIVE_INFINITY,
					Double.POSITIVE_INFINITY };
			try {
				v = super.gradient(x, p);
			} catch (NotStrictlyPositiveException e) { // NOPMD
				// Do nothing.
			}
			return v;
		}
	};
	/** Initial guess. */
	private double[] initialGuess;
	/** Maximum number of iterations of the optimization algorithm. */
	private final int maxIter;

	/**
	 * Contructor used by the factory methods.
	 *
	 * @param initialGuess Initial guess. If set to {@code null}, the initial guess
	 * will be estimated using the {@link ParameterGuesser}.
	 * @param maxIter Maximum number of iterations of the optimization algorithm.
	 */
	private GaussianArrayCurveFitter(double[] initialGuess,
			int maxIter) {
		this.initialGuess = initialGuess;
		this.maxIter = maxIter;
	}

	/**
	 * Creates a default curve fitter.
	 * The initial guess for the parameters will be {@link ParameterGuesser}
	 * computed automatically, and the maximum number of iterations of the
	 * optimization algorithm is set to {@link Integer#MAX_VALUE}.
	 *
	 * @return a curve fitter.
	 *
	 * @see #withStartPoint(double[])
	 * @see #withMaxIterations(int)
	 */
	public static GaussianArrayCurveFitter create() {
		return new GaussianArrayCurveFitter(null, Integer.MAX_VALUE);
	}

	/**
	 * Configure the start point (initial guess).
	 * @param newStart new start point (initial guess)
	 * @return a new instance.
	 */
	public GaussianArrayCurveFitter withStartPoint(double[] newStart) {
		return new GaussianArrayCurveFitter(newStart.clone(),
				maxIter);
	}

	/**
	 * Configure the maximum number of iterations.
	 * @param newMaxIter maximum number of iterations
	 * @return a new instance.
	 */
	public GaussianArrayCurveFitter withMaxIterations(int newMaxIter) {
		return new GaussianArrayCurveFitter(initialGuess,
				newMaxIter);
	}

	/** {@inheritDoc} */
	@Override
	protected LeastSquaresProblem getProblem(Collection<WeightedObservedPoint> observations) {

		// Prepare least-squares problem.
		final int len = observations.size();
		final double[] xx      = new double[len];
		final double[] target  = new double[len];
		final double[] weights = new double[len];

		int i = 0;
		for (WeightedObservedPoint obs : observations) {
			xx[i]      = obs.getX();
			target[i]  = obs.getY();
			weights[i] = obs.getWeight();
			++i;
		}

		final AbstractCurveFitter.TheoreticalValuesFunction model =
				new AbstractCurveFitter.TheoreticalValuesFunction(FUNCTION, observations);
		
		final double[] startPoint = initialGuess != null ?
				initialGuess :
					// Compute estimation.
					// TODO: Check consistency with guessing function
					new ParameterGuesser(observations, initialGuess.length/3).guess();
		
		final GaussianArrayParameterValidator parValid = 
				new GaussianArrayParameterValidator(initialGuess, xx, target);

		// Return a new least squares problem set up to fit a Gaussian curve to the
		// observed points
		return new LeastSquaresBuilder().
				parameterValidator(parValid).
				maxEvaluations(Integer.MAX_VALUE).
				maxIterations(maxIter).
				lazyEvaluation(false).
				start(startPoint).
				target(target).
				weight(new DiagonalMatrix(weights)).
				model(model.getModelFunction(), model.getModelFunctionJacobian()).
				build();

	}
	
	/** {@inheritDoc} */
	@Override
	protected LeastSquaresOptimizer getOptimizer() {
		// Return a new least squares problem set up to fit a Gaussian curve to the
		// observed points.
		return new LevenbergMarquardtOptimizer()
				.withCostRelativeTolerance(1e-15)
				.withOrthoTolerance(1e-15)
				.withParameterRelativeTolerance(1e-15);
	}
	/**
	 * Guesses the parameters {@code norm}, {@code mean}, and {@code sigma}
	 * of a {@link org.apache.commons.math4.analysis.function.GaussianArray.Parametric}
	 * based on the specified observed points.
	 */
	public static class ParameterGuesser {
		/** Normalization factor. */
		private final double[] norm;
		/** Mean. */
		private final double[] mean;
		/** Standard deviation. */
		private final double[] sigma;

		/**
		 * Constructs instance with the specified observed points.
		 *
		 * @param observations Observed points from which to guess the
		 * parameters of the Gaussian.
		 * @throws NullArgumentException if {@code observations} is
		 * {@code null}.
		 * @throws NumberIsTooSmallException if there are less than 3
		 * observations.
		 */
		public ParameterGuesser(Collection<WeightedObservedPoint> observations, int gaussN) {
			if (observations == null) {
				throw new NullArgumentException(LocalizedFormats.INPUT_ARRAY);
			}
			if (observations.size() < 3) {
				throw new NumberIsTooSmallException(observations.size(), 3, true);
			}

			final List<WeightedObservedPoint> sorted = sortObservations(observations);
			final double[] params = basicGuess(sorted.toArray(new WeightedObservedPoint[0]));
			norm  = new double[gaussN];
			mean  = new double[gaussN];
			sigma = new double[gaussN];
			for (int p=0; p<params.length; p++ ){
				norm [p/3] = params[p];
				mean [p/3] = params[p+1];
				sigma[p/3] = params[p+2];
			}
		}

		/**
		 * Gets an estimation of the parameters.
		 *
		 * @return the guessed parameters, in the following order:
		 * <ul>
		 *  <li>Normalization factor</li>
		 *  <li>Mean</li>
		 *  <li>Standard deviation</li>
		 * </ul>
		 */
		public double[] guess() {
			// TODO:
			return null;
			//return new double[] { norm, mean, sigma };
		}

		/**
		 * Sort the observations.
		 *
		 * @param unsorted Input observations.
		 * @return the input observations, sorted.
		 */
		private List<WeightedObservedPoint> sortObservations(Collection<WeightedObservedPoint> unsorted) {
			final List<WeightedObservedPoint> observations = new ArrayList<WeightedObservedPoint>(unsorted);

			final Comparator<WeightedObservedPoint> cmp = new Comparator<WeightedObservedPoint>() {
				/** {@inheritDoc} */
				@Override
				public int compare(WeightedObservedPoint p1,
						WeightedObservedPoint p2) {
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
		 * @return the guessed parameters (normalization factor, mean and
		 * sigma).
		 */
		private double[] basicGuess(WeightedObservedPoint[] points) {
			// TODO: Guess multiple peaks
			final int maxYIdx = findMaxY(points);
			final double n = points[maxYIdx].getY();
			final double m = points[maxYIdx].getX();

			double fwhmApprox;
			try {
				final double halfY = n + ((m - n) / 2);
				final double fwhmX1 = interpolateXAtY(points, maxYIdx, -1, halfY);
				final double fwhmX2 = interpolateXAtY(points, maxYIdx, 1, halfY);
				fwhmApprox = fwhmX2 - fwhmX1;
			} catch (OutOfRangeException e) {
				// TODO: Exceptions should not be used for flow control.
				fwhmApprox = points[points.length - 1].getX() - points[0].getX();
			}
			final double s = fwhmApprox / (2 * FastMath.sqrt(2 * FastMath.log(2)));

			return new double[] { n, m, s };
		}

		/**
		 * Finds index of point in specified points with the largest Y.
		 *
		 * @param points Points to search.
		 * @return the index in specified points array.
		 */
		private int findMaxY(WeightedObservedPoint[] points) {
			int maxYIdx = 0;
			for (int i = 1; i < points.length; i++) {
				if (points[i].getY() > points[maxYIdx].getY()) {
					maxYIdx = i;
				}
			}
			return maxYIdx;
		}

		/**
		 * Interpolates using the specified points to determine X at the
		 * specified Y.
		 *
		 * @param points Points to use for interpolation.
		 * @param startIdx Index within points from which to start the search for
		 * interpolation bounds points.
		 * @param idxStep Index step for searching interpolation bounds points.
		 * @param y Y value for which X should be determined.
		 * @return the value of X for the specified Y.
		 * @throws ZeroException if {@code idxStep} is 0.
		 * @throws OutOfRangeException if specified {@code y} is not within the
		 * range of the specified {@code points}.
		 */
		private double interpolateXAtY(WeightedObservedPoint[] points,
				int startIdx,
				int idxStep,
				double y)
						throws OutOfRangeException {
			if (idxStep == 0) {
				throw new ZeroException();
			}
			final WeightedObservedPoint[] twoPoints
			= getInterpolationPointsForY(points, startIdx, idxStep, y);
			final WeightedObservedPoint p1 = twoPoints[0];
			final WeightedObservedPoint p2 = twoPoints[1];
			if (p1.getY() == y) {
				return p1.getX();
			}
			if (p2.getY() == y) {
				return p2.getX();
			}
			return p1.getX() + (((y - p1.getY()) * (p2.getX() - p1.getX())) /
					(p2.getY() - p1.getY()));
		}

		/**
		 * Gets the two bounding interpolation points from the specified points
		 * suitable for determining X at the specified Y.
		 *
		 * @param points Points to use for interpolation.
		 * @param startIdx Index within points from which to start search for
		 * interpolation bounds points.
		 * @param idxStep Index step for search for interpolation bounds points.
		 * @param y Y value for which X should be determined.
		 * @return the array containing two points suitable for determining X at
		 * the specified Y.
		 * @throws ZeroException if {@code idxStep} is 0.
		 * @throws OutOfRangeException if specified {@code y} is not within the
		 * range of the specified {@code points}.
		 */
		private WeightedObservedPoint[] getInterpolationPointsForY(WeightedObservedPoint[] points,
				int startIdx,
				int idxStep,
				double y)
						throws OutOfRangeException {
			if (idxStep == 0) {
				throw new ZeroException();
			}
			for (int i = startIdx;
					idxStep < 0 ? i + idxStep >= 0 : i + idxStep < points.length;
					i += idxStep) {
				final WeightedObservedPoint p1 = points[i];
				final WeightedObservedPoint p2 = points[i + idxStep];
				if (isBetween(y, p1.getY(), p2.getY())) {
					if (idxStep < 0) {
						return new WeightedObservedPoint[] { p2, p1 };
					} else {
						return new WeightedObservedPoint[] { p1, p2 };
					}
				}
			}

			// Boundaries are replaced by dummy values because the raised
			// exception is caught and the message never displayed.
			// TODO: Exceptions should not be used for flow control.
			throw new OutOfRangeException(y,
					Double.NEGATIVE_INFINITY,
					Double.POSITIVE_INFINITY);
		}

		/**
		 * Determines whether a value is between two other values.
		 *
		 * @param value Value to test whether it is between {@code boundary1}
		 * and {@code boundary2}.
		 * @param boundary1 One end of the range.
		 * @param boundary2 Other end of the range.
		 * @return {@code true} if {@code value} is between {@code boundary1} and
		 * {@code boundary2} (inclusive), {@code false} otherwise.
		 */
		private boolean isBetween(double value,
				double boundary1,
				double boundary2) {
			return (value >= boundary1 && value <= boundary2) ||
					(value >= boundary2 && value <= boundary1);
		}
	}
	
	/**
	 * Checks the parameter triplets {@code norm}, {@code mean}, and {@code sigma}
	 * of a {@code GaussianArray.Parametric} for applied constraints,
	 *  based on the specified points and initial guess.
	 */
	public static class GaussianArrayParameterValidator 
						implements ParameterValidator {
		private double[] initialParameterSet;
		private double   maxXvalue;
		private double   maxMeandiff;
		public GaussianArrayParameterValidator(double[] initialGuess, double[] xtarget, double[] ytarget) {
			this.initialParameterSet = initialGuess;
			this.maxXvalue = xtarget[xtarget.length-1];
			this.maxMeandiff = maxXvalue/(initialGuess.length/3)/50;
//			for (int p = 1; p<initialGuess.length-4; p+=3){
//				if (initialGuess[p+3]-initialGuess[p]<maxMeandiff) 
//					maxMeandiff = initialGuess[p+3]-initialGuess[p];
//			} // any mean not further than this from initial guess
		}
		
		@Override
		// Set parameter constraints here
		public RealVector validate(RealVector pars) {
			for (int i = 0; i<pars.getDimension(); i++) {
				// {0=Norm, 1=Mean, 2=Sigma}
				// Lower bound of 0 for all 3
			    if ((i % 3 == 0 || i % 3 == 1 || i % 3 == 2) && pars.getEntry(i) < 1e-1)  
			    	pars.setEntry(i, 0.0);
			    // Upper bounds for each parameter
			    if ((i % 3 == 0) && pars.getEntry(i) > initialParameterSet[i])  
			    	pars.setEntry(i, initialParameterSet[i]);
			    if ((i % 3 == 1) && 
			    		(pars.getEntry(i) > maxXvalue ))  
			    	pars.setEntry(i, maxXvalue);
			    if ((i % 3 == 2) && pars.getEntry(i) > 1.0)  
			    	pars.setEntry(i, 1.0);
			    // Keep means close to maxima
			    double diff = pars.getEntry(i)-initialParameterSet[i];
			    if ((i % 3 == 1) && (Math.abs(diff) > maxMeandiff))
			    	pars.setEntry(i,initialParameterSet[i]+maxMeandiff*diff/Math.abs(diff));
			}
			return pars;
		}
		
	}
}

class GaussianArray implements UnivariateDifferentiableFunction {
	private double[] means;
	private double[] norms;
	
	private double[] is;  // 1*std
	private double[] i2s2;// 1/(2*std^2)
	
	public GaussianArray(	double[] norms,
							double[] means,
							double[] stds)
					throws NotStrictlyPositiveException {
	
		this.means = means;
		this.norms = norms;
		for (int q = 0; q< means.length; q++) {
			this.is  [q] = 1 / stds[q];
			this.i2s2[q] = 0.5*is[q]*is[q];
		}
	}
	
	@Override
	public double value(double x) {
		RealVector diffs = new ArrayRealVector(means);
		diffs.mapAddToSelf(-x).mapMultiplyToSelf(-1);
		return value(diffs.toArray(), norms, i2s2);
	}
	
	private static double value(double[] xMinusMean, 
								double[] norms, 
								double[] i2s2) {
		double output = 0;
		for (int i=0;i<xMinusMean.length;i++) {
			output += norms[i] * FastMath.exp(-xMinusMean[i] * xMinusMean[i] * i2s2[i]);
		}
		return output;
	} 
	
	@Override
	public DerivativeStructure value(DerivativeStructure t) throws DimensionMismatchException {
		// TODO Auto-generated method stub
		//final double[] u = new double[] {is.multiplyToSelf(means.subtractToself(t.getValue()).multiplyToSelf(-1)).toArray()};
		
        double[] f = new double[t.getOrder() + 1];

        // the nth order derivative of the Gaussian has the form:
        // dn(g(x)/dxn = (norm / s^n) P_n(u) exp(-u^2/2) with u=(x-m)/s
        // where P_n(u) is a degree n polynomial with same parity as n
        // P_0(u) = 1, P_1(u) = -u, P_2(u) = u^2 - 1, P_3(u) = -u^3 + 3 u...
        // the general recurrence relation for P_n is:
        // P_n(u) = P_(n-1)'(u) - u P_(n-1)(u)
        // as per polynomial parity, we can store coefficients of both P_(n-1) and P_n in the same array
		return null;
	}
	
	/**
	 * Parametric function where the input array contains the parameters of
	 * the Gaussian, ordered as follows:
	 * <ul>
	 *  <li>Norm</li>
	 *  <li>Mean</li>
	 *  <li>Standard deviation</li>
	 * </ul>
	 */
	public static class Parametric implements ParametricUnivariateFunction {
		/**
		 * Computes the value of the Gaussian at {@code x}.
		 *
		 * @param x Value for which the function must be computed.
		 * @param param Values of norm, mean and standard deviation,
		 * 		passed in subsequent triplets.
		 * @return the value of the function.
		 * @throws NullArgumentException if {@code param} is {@code null}.
		 * @throws DimensionMismatchException if the size of {@code param} is
		 * not a multiple of 3.
		 * @throws NotStrictlyPositiveException if {@code param[2]} is negative.
		 */
		@Override
		public double value(double x, double... pars) {
			validateParameters(pars);
			final double[] diffs = new double[pars.length / 3];
			final double[] norms = new double[pars.length / 3];
			final double[] i2s2  = new double[pars.length / 3];
			
			for (int i = 0; i<pars.length; i+=3) {
				diffs[i/3] = x - pars[i+1];
				norms[1/3] = pars[i];
				i2s2 [i/3] = 1/(2*pars[i+2] * pars[i+2]);
			}
			return GaussianArray.value(diffs,norms,i2s2);
		}
		
        /**
         * Computes the value of the gradient at {@code x}.
         * The components of the gradient vector are the partial
         * derivatives of the function with respect to each of the
         * <em>parameters</em> (norm, mean and standard deviation,
         * for each Gaussian in the array).
         *
         * @param x Value at which the gradient must be computed.
         * @param param Values of norm, mean and standard deviation triplets.
         * @return the gradient vector at {@code x}.
         * @throws NullArgumentException if {@code param} is {@code null}.
         * @throws DimensionMismatchException if the size of {@code param} is
         * not 3.
         * @throws NotStrictlyPositiveException if {@code param[2]} is negative.
         */
		@Override
		public double[] gradient(double x, double ... param)
                throws NullArgumentException,
                       DimensionMismatchException,
                       NotStrictlyPositiveException {
                validateParameters(param);
                double[] out = new double[param.length];
                
            	final double[] norm  = new double[param.length/3];
            	final double[] diff  = new double[param.length/3];
            	final double[] sigma = new double[param.length/3];
            	final double[] i2s2  = new double[param.length/3];
                
                for (int gg=0; gg<param.length; gg+=3) {
                	norm [gg/3] =     param[gg]  ;
                	diff [gg/3] = x - param[gg+1];
                	sigma[gg/3] =     param[gg+2];
                	i2s2 [gg/3] = 1 / (2 * sigma[gg/3] * sigma[gg/3]);
                	
                	// Only 1 Gaussian function at a time contributes to the gradient
                	final double n = new Gaussian(1.0, param[gg+1],param[gg+2]).value(x);
                	final double m = norm[gg/3] * n * 2 * i2s2[gg/3] * diff[gg/3];
                	final double s = m * diff[gg/3] / sigma[gg/3];
                	out[gg]   = n;
                	out[gg+1] = m;
                	out[gg+2] = s;
                }
                return out;
            }

		/**
		 * Validates parameters to ensure they are appropriate for the evaluation of
		 * the {@link #value(double,double[])} and {@link #gradient(double,double[])}
		 * methods.
		 *
		 * @param param Values of norm, mean and standard deviation.
		 * @throws NullArgumentException if {@code param} is {@code null}.
		 * @throws DimensionMismatchException if the size of {@code param} is
		 * not 3.
		 * @throws NotStrictlyPositiveException if {@code param[2]} is negative.
		 */
		private void validateParameters(double[] param)
				throws NullArgumentException,
				DimensionMismatchException,
				NotStrictlyPositiveException {
			if (param == null) {
				throw new NullArgumentException();
			}
			if ((param.length % 3) != 0) {
				throw new DimensionMismatchException(param.length % 3, 0);
			}
			
			for (int p=0; p<param.length;p++)
				if ((p % 3 ==2) && param[p] <= 0) {
					throw new NotStrictlyPositiveException(param[p]);
				}
			}
		}
}

