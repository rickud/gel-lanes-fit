/**
 * Gel Lanes Fit
 * GelLanesFit.java
 * author: Rick Ziraldo, 2017
 * The /University of Texas at Dallas, Richardson, TX
 * http://www.utdallas.edu
 *
 * The source code is maintained and made available on GitHub
 * https://github.com/rickud/gauss-curve-fit
 *
 */

package gellanesfit;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.analysis.ParametricUnivariateFunction;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.analysis.differentiation.UnivariateDifferentiableFunction;
import org.apache.commons.math3.analysis.function.Abs;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.analysis.function.Log;
import org.apache.commons.math3.analysis.function.Log10;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.exception.NumberIsTooSmallException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.fitting.AbstractCurveFitter;
import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresBuilder;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.fitting.leastsquares.ParameterValidator;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DiagonalMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.util.FastMath;

class GaussianArrayCurveFitter extends AbstractCurveFitter {

	/** Parametric function to be fitted. */
	private static final GaussianArray.Parametric FUNCTION =
		new GaussianArray.Parametric();
	private static final int bandMode = 0;
	private static final int continuumMode = 1;
	private static final double sd2FWHM = 2 * FastMath.sqrt(2 * FastMath.log(2));
	

	/** Initial guess. */
	private final SortedParameters initialGuess;
	/** Maximum number of iterations of the optimization algorithm. */
	private final int maxIter;
	private final int fitMode;
	private final int deg;
	private final double peakTol, areaDrift, sdDrift, polyDerivative, polyOffset;

	/**
	 * Constructor used by the factory methods.
	 *
	 * @param initialGuess Initial guess. If set to {@code null}, the initial
	 *          guess will be estimated using the {@link ParameterGuesser}.
	 * @param maxIter Maximum number of iterations of the optimization algorithm.
	 */
	private GaussianArrayCurveFitter(final SortedParameters initialGuess,
		final int maxIter, final int fitMode, final int deg,
		final double polyDerivative, final double peakTol,
		final double polyOffset, final double areaDrift, final double sdDrift)
	{
		this.initialGuess = initialGuess;
		this.maxIter = maxIter;
		this.fitMode = fitMode;
		this.deg = deg;
		this.polyDerivative = polyDerivative;
		this.peakTol = peakTol;
		this.polyOffset = polyOffset;
		this.areaDrift = areaDrift;
		this.sdDrift = sdDrift;

//		f.getContentPane().add(new ChartPanel(chart));
	}

	/**
	 * Creates a default curve fitter. The initial guess for the parameters will
	 * be {@link ParameterGuesser} computed automatically, and the maximum number
	 * of iterations of the optimization algorithm is set to
	 * {@link Integer#MAX_VALUE}.
	 *
	 * @param fitMode
	 * @param sdDrift 
	 * @return a curve fitter.
	 * @see #withStartPoint(final SortedParameters newStart)
	 */
	static GaussianArrayCurveFitter create(final int fitMode, final int deg,
		final double polyDerivative, final double polyOffset, final double peakTol, final double areaDrift, final double sdDrift)
	{
		return new GaussianArrayCurveFitter(null, Integer.MAX_VALUE, fitMode, deg,
			polyDerivative, polyOffset, peakTol, areaDrift, sdDrift);
	}

	/**
	 * Configure the start point (initial guess).
	 *
	 * @param newStart new start point (initial guess)
	 * @return a new instance.
	 */
	GaussianArrayCurveFitter withStartPoint(final SortedParameters newStart) {
		return new GaussianArrayCurveFitter(newStart, maxIter, fitMode, deg,
			polyDerivative, polyOffset, peakTol, areaDrift, sdDrift);
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

		final SortedParameters startPoint = initialGuess != null ? initialGuess
			: new ParameterGuesser(observations, deg, peakTol, 1.0).guess();

		final GaussianArrayParameterValidator parValid =
			new GaussianArrayParameterValidator(fitMode, startPoint, xx, target,
				polyDerivative, polyOffset, areaDrift, sdDrift);

		return new LeastSquaresBuilder().parameterValidator(parValid)
			.maxEvaluations(Integer.MAX_VALUE).maxIterations(maxIter).lazyEvaluation(
				false).start(startPoint.getParameters()).target(target).weight(new DiagonalMatrix(
					weights)).model(model.getModelFunction(), model
						.getModelFunctionJacobian()).build();
	}

	/** {@inheritDoc} */
	@Override
	protected LeastSquaresOptimizer getOptimizer() {
		return new LevenbergMarquardtOptimizer()
			.withCostRelativeTolerance(1e-18)
			.withOrthoTolerance(1e-18)
			.withParameterRelativeTolerance(1e-15);
	}

	/**
	 * Guesses the parameters {@code norm}, {@code mean}, and {@code sigma} of a
	 * {@link GaussianArray.Parametric} based on the specified observed points.
	 */
	static class ParameterGuesser {
		private final int fitMode;
		private final List<WeightedObservedPoint> sorted;
		private final int deg;
		private final double peakTol, polyOffset;
		private final RealMatrix distMatrix;
		private final List<Peak> ladderPeaks;
		private final double[] ladderMW;
		
		/**
		 * Constructs instance with the specified observed points.
		 *
		 * @param observations Observed points from which to guess the parameters of
		 *          the train of Gaussians.
		 * @throws NullArgumentException if {@code observations} is {@code null}.
		 * @throws NumberIsTooSmallException if there are less than 3 observations.
		 */
		ParameterGuesser(final Collection<WeightedObservedPoint> observations,
			final int deg, final double peakTol, final double polyOffset)
		{
			fitMode = GaussianArrayCurveFitter.bandMode;
			this.deg = deg;
			this.peakTol = peakTol;
			this.polyOffset = polyOffset;
			this.distMatrix = null;
			this.ladderPeaks = null;
			this.ladderMW = null;
			
			if (observations == null) {
				throw new NullArgumentException(LocalizedFormats.INPUT_ARRAY);
			}
			if (observations.size() < 3) {
				throw new NumberIsTooSmallException(observations.size(), 3, true);
			}
		
			this.sorted = sortObservations(observations);
		}

		ParameterGuesser(final Collection<WeightedObservedPoint> observations,
			final int deg, final double peakTol, final double polyOffset,
			final RealMatrix distMatrix, final List<Peak> ladderPeaks,
			final double[] ladderMW)
		{
			fitMode = GaussianArrayCurveFitter.continuumMode;
			this.deg = deg;
			this.peakTol = peakTol;
			this.polyOffset = polyOffset;
			this.distMatrix = distMatrix;
			this.ladderPeaks = ladderPeaks;
			this.ladderMW = ladderMW;
			
			if (observations == null) {
				throw new NullArgumentException(LocalizedFormats.INPUT_ARRAY);
			}
			if (observations.size() < 3) {
				throw new NumberIsTooSmallException(observations.size(), 3, true);
			}
		
			this.sorted = sortObservations(observations);
		}
		
		/**
		 * Gets an estimation of the parameters.
		 *
		 * @return the guessed parameters, in the following order:
		 *         <ul>
		 *         <li>Degree of the polynomial background</li>
		 *         <li>Coefficients of the polynomial</li>
		 *         <li>Normalization factor</li>
		 *         <li>Mean</li>
		 *         <li>Standard deviation</li>
		 *         </ul>
		 */
		SortedParameters guess() {
			return basicGuess(sorted.toArray(new WeightedObservedPoint[0]), peakTol);
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
		 * @return the guessed parameters (amplitude, mean and sigma).
		 */
		private SortedParameters basicGuess(final WeightedObservedPoint[] points,
			final double tolpk)
		{
			RealVector xvals = new ArrayRealVector();
			RealVector yvals = new ArrayRealVector();

			for (int o = 0; o < points.length; o++) {
				xvals = xvals.append(points[o].getX());
				yvals = yvals.append(points[o].getY());
			}
			
			final double minX = xvals.getMinValue();
			final double maxX = xvals.getMaxValue();
			final double minY = yvals.getMinValue();
			final double maxY = yvals.getMinValue();
			final double slope = (maxY - minY) / (maxX - minX);
			final double p0 = minY - slope * minX;
			
			
			RealVector normG = new ArrayRealVector();
			RealVector meanG = new ArrayRealVector();
			RealVector sdG = new ArrayRealVector();
			RealVector polyG = new ArrayRealVector(deg + 1);
			if (deg >= 0) 
				polyG.setEntry(0,p0);
			
			if (fitMode == GaussianArrayCurveFitter.bandMode) {
				// Local maxima where the peaks are
				final int[] maximaIdx = findMaxima(points, tolpk, true);
	
				// Define a Gaussian at each peak
				final double[] means = new double[maximaIdx.length];
				final double[] sds = new double[maximaIdx.length];
				final double[] norms = new double[maximaIdx.length];
	
				for (int m = 0; m < maximaIdx.length; m++) {
					norms[m] = yvals.getEntry(maximaIdx[m]) - minY;
					means[m] = xvals.getEntry(maximaIdx[m]);
				}
	
				// estimate sds from maxima and mean
				for (int m = 0; m < maximaIdx.length; m++) {
					boolean foundFWHM = false; // Full width at half maximum
					boolean foundRWHM = false;
					boolean foundLWHM = false;
					double LWHM = 0.0;
					double RWHM = 0.0;
					double FWHM = (xvals.getMaxValue() - xvals.getMinValue()) / 2;
	
					final double yRange = yvals.getEntry(maximaIdx[m]) - minY;
					double hm = yRange / 2.0; // Actual profile value; incledes offset
					final int pkPos = maximaIdx[m];
	
					final double inc = 1.05;
					double peakDistance;
					if (maximaIdx.length == 1) {
						peakDistance = FastMath.min(means[0] - minX, maxX - means[0]);
					}
					else if (m == 0) {
						peakDistance = FastMath.min(means[m + 1] - means[m], means[m] - minX);
					}
					else if (m + 1 == maximaIdx.length) {
						peakDistance = FastMath.min(maxX - means[m], means[m] - means[m - 1]);
					}
					else {
						peakDistance = FastMath.min(means[m + 1] - means[m], means[m] -
							means[m - 1]);
					}
	
					while ((!foundFWHM || FWHM > peakDistance) && hm * inc < 0.9 * yRange) {
						foundFWHM = false;
						foundRWHM = false;
						foundLWHM = false;
						// Right side, check 3 consecutive points for smoothing
						int p = 0;
						while (pkPos + p + 2 < xvals.getDimension() && !foundRWHM) {
							if (yvals.getEntry(pkPos + p) < (hm + minY) && yvals.getEntry(
								pkPos + p + 1) < (hm + minY) && yvals.getEntry(pkPos + p +
									2) < (hm + minY))
							{
								foundRWHM = true;
								RWHM = xvals.getEntry(pkPos + p) - means[m];
							}
							p++;
						}
	
						// Left side, check 3 consecutive points for smoothing
						p = 0;
						while (pkPos - p - 2 >= 0 && !foundLWHM) {
							if (yvals.getEntry(pkPos - p) < (hm + minY) && yvals.getEntry(
								pkPos - p - 1) < (hm + minY) && yvals.getEntry(pkPos - p -
									2) < (hm + minY))
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
					normG = normG.append(norms[m]);
					meanG = meanG.append(means[m]);
					sdG   = sdG.append(sds[m]);
				}
			}
			else if (fitMode == continuumMode) {
				final RealVector profile = yvals.mapSubtractToSelf(yvals.getMinValue()*polyOffset);
				
				PolynomialSplineFunction pr = new LinearInterpolator()
						.interpolate(xvals.toArray(), profile.toArray());

				// Use the stored distribution as a guess
				// fragmentDistribution[:][0] = Frequency
				// fragmentDistribution[:][1] = Length (bp)
				// fragmentDistribution[:][2] = MW
	
				final double[] meanLadder = new double[ladderPeaks.size()];
				final double[] sdLadder = new double[ladderPeaks.size()];
				for (int p = 0; p < ladderPeaks.size(); p++) {
					meanLadder[p] = ladderPeaks.get(p).getMean();
					sdLadder[p] = ladderPeaks.get(p).getSigma();
				}
				final double[] means = interpolateDisplacement(meanLadder, ladderMW, distMatrix.getColumnVector(2));
				final double[] sds = interpolateSD(meanLadder, sdLadder, means);
				final RealVector scaledFrequency = distMatrix.getColumnVector(0);
				for (int s = 0; s < scaledFrequency.getDimension(); s++) {
					double pv = 0.0;
					if (means[s] < xvals.getMinValue()) {
						double pv0 = pr.getPolynomials()[0].value(means[s] - pr.getKnots()[0]);
						pv = FastMath.max(0.0, pv0);
					} else if (means[s] > xvals.getMaxValue()) {
						int n = pr.getPolynomials().length - 1;
						double pv0 = pr.getPolynomials()[n].value(means[s] - pr.getKnots()[n]);
						pv = FastMath.max(0.0, pv0);
					}	else 
						pv = pr.value(means[s]);
					
					scaledFrequency.setEntry(s, scaledFrequency.getEntry(s) * pv);
				}
				final RealVector mwArray = distMatrix.getColumnVector(2);
				RealVector scale = scaledFrequency.ebeMultiply(mwArray.map(new Log()));
				scale = scale.mapDivide(scale.getMaxValue());
			
				List<Integer> fragmentSubset = new ArrayList<>();
				double margin = (xvals.getMaxValue() - xvals.getMinValue() )* 0.2;
				for (int i = 0; i < means.length; i++) {
					if (means[i] > xvals.getMinValue() - margin &&
							means[i] < xvals.getMaxValue() + margin) {
						fragmentSubset.add(i);
						meanG = meanG.append(means[i]);
					}
				}
				for (int i : fragmentSubset) {
					sdG = sdG.append(sds[i]);
					normG = normG.append(scale.getEntry(i));
				}
				double scale2  = profile.getMaxValue() 
						/ FastMath.log(normG.getDimension());
				normG = normG.mapMultiply(scale2);
			}
			return new SortedParameters(polyG, normG, meanG, sdG);
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
			Arrays.sort(maxPositions);
			return maxPositions;
		}
	
	
		private double[] interpolateDisplacement(final double[] y,
			final double[] ladder, final RealVector dist) {
			final RealVector logs = dist.map(new Log10());
			RealVector yi = new ArrayRealVector();
			// Linear interpolator between ladder points
			final LinearInterpolator li = new LinearInterpolator();
			final double[] l = new ArrayRealVector(ladder).map(new Log10()).toArray();
			ArrayUtils.reverse(l);
			ArrayUtils.reverse(y);
			final PolynomialSplineFunction f = li.interpolate(l, y);
			final double[] kn = f.getKnots();
			final PolynomialFunction[] fi = f.getPolynomials();
	
			for (int i = 0; i < logs.getDimension(); i++) {
				final double logsi = logs.getEntry(i);
				if (logsi < kn[0]) {
					yi = yi.append(fi[0].value(logsi - kn[0]));
				}
				else if (logsi > kn[kn.length - 1]) {
					yi = yi.append(fi[fi.length - 1].value(logsi - kn[kn.length - 2]));
				}
				else {
					yi = yi.append(f.value(logsi));
				}
			}
			return yi.toArray();
		}
	
		private double[] interpolateSD(final double[] y, final double[] sd,
			final double[] yi)
		{
			final WeightedObservedPoints obs = new WeightedObservedPoints();
			for (int l = 0; l < y.length; l++)
				obs.add(y[y.length - 1 - l], sd[l]);
	
			// First-degree polynomial fitter (line)
			final PolynomialCurveFitter linfit = PolynomialCurveFitter.create(1);
			final double[] coeffs = linfit.fit(obs.toList());
			final UnivariateFunction f = new PolynomialFunction(coeffs);
			final RealVector sdi = new ArrayRealVector(yi).map(f);
			return sdi.toArray();
		}

		public List<Integer> getUsedFragments(double[] xrange) {
			final double[] meanLadder = new double[ladderPeaks.size()];
			for (int p = 0; p < ladderPeaks.size(); p++) {
				meanLadder[p] = ladderPeaks.get(p).getMean();
			}
			List<Integer> out = new ArrayList<>(); 
			final double[] means = interpolateDisplacement(meanLadder, ladderMW, distMatrix.getColumnVector(2));
//			double margin = (xrange[1] - xrange[0]) * 0.2;
			double margin = 0.0;
			for (int i = 0; i < means.length; i++) {
				if (means[i] > xrange[0] - margin && means[i] < xrange[1] + margin) out.add(i);
			}
			return out;
		}
	}
	/**
	 * Checks on the polynomial coefficients {@code poly} and the parameter
	 * triplets {@code norm}, {@code mean}, and {@code sigma} of a
	 * {@code GaussianArray.Parametric} for applied constraints, based on the
	 * specified points and initial guess.
	 */
	private static class GaussianArrayParameterValidator implements
		ParameterValidator
	{

		private final int fitMode;
		private final SortedParameters iniSP;
		private final RealVector xtarget;
		private final RealVector ytarget;

		private final int deg;

		private final double areaDrift;
		private RealVector maxMeanDiff = new ArrayRealVector();

		private final double maxX, minX, maxY, minY, margin;
		private final double minN;
		private double minSD;
		private double maxSD;
		private final double minD1;
		private final double maxD1;
		private final double polyOffset;
		private final PolynomialSplineFunction profile;

		private GaussianArrayParameterValidator(final int fitMode,
			final SortedParameters iniSP, final double[] xtarget,
			final double[] ytarget, final double polyDerivative, final double polyOffset,
			final double areaDrift, final double sdDrift)
		{

			this.fitMode = fitMode;
			this.iniSP = iniSP;
			this.xtarget = new ArrayRealVector(xtarget);
			this.ytarget = new ArrayRealVector(ytarget);
			this.deg = iniSP.getDeg();
			this.areaDrift = areaDrift;

			this.polyOffset = polyOffset; // proportion of the profile value
			if (fitMode == GaussianArrayCurveFitter.continuumMode)
				this.margin = (this.xtarget.getMaxValue() - this.xtarget.getMinValue()) * 0.2;
			else this.margin = 0.0;
			minX = this.xtarget.getMinValue() - margin;
			maxX = this.xtarget.getMaxValue() + margin;
			minY = this.ytarget.getMinValue();
			maxY = this.ytarget.getMaxValue();

			maxD1 = polyDerivative;
			minD1 = -polyDerivative;
			this.profile = new LinearInterpolator().interpolate(xtarget, ytarget);
			double mds = 0.8; // Distance from guess peak mean, as fraction of initial inter-peak distance
			if (fitMode == continuumMode) {
				mds = 0.1;
			}
			if (iniSP.getMean().getDimension() > 1) {
				maxMeanDiff = iniSP.getMean().getSubVector(1, iniSP.getMean().getDimension() - 1).subtract(
					iniSP.getMean().getSubVector(0, iniSP.getMean().getDimension() - 1)).map(new Abs())
					.mapMultiply(mds);
				maxMeanDiff = maxMeanDiff.append(maxMeanDiff.getEntry(maxMeanDiff
					.getDimension() - 1));
			}
			else {
				maxMeanDiff = maxMeanDiff.append((maxX - minX) / 2.0);
			}

			minN = 0.01; // proportion of the profile-bg difference

			minSD = 0.4; // proportion of sd0[i]
			maxSD = 2.0;
			if (fitMode == continuumMode) { // controllable from interface
				minSD = 1/sdDrift;
				maxSD = sdDrift;
			}
		}

		@Override
		// Set parameter constraints here
		public RealVector validate(final RealVector param) {
			// Sort the parameter array the same way as the initial array
			RealVector poly = new ArrayRealVector();
			RealVector norm = new ArrayRealVector();
			RealVector mean = new ArrayRealVector();
			RealVector sd = new ArrayRealVector();
			RealVector area = new ArrayRealVector();
			PolynomialFunction p = new PolynomialFunction(new double[] { 0.0 });
			if (deg != -1) { // no polynomial
				poly = param.getSubVector(1, deg + 1);
			}
			// Gaussian parameters in order {Norm, Mean, Sigma}
			for (int i = deg + 2; i < param.getDimension(); i += 3) {
				norm = norm.append(param.getEntry(i));
				mean = mean.append(param.getEntry(i + 1));
				sd = sd.append(param.getEntry(i + 2));
			}

			// Polynomyal parameters
			if (poly.getDimension() > 0) {
				final double tolHigh = ytarget.getMinValue() * polyOffset;
				final double tolLow = 0.4 * tolHigh;
				p = new PolynomialFunction(poly.toArray());
				PolynomialFunction p1 = p.polynomialDerivative();
				double d1bg = new Mean().evaluate(xtarget.map(p1).toArray());
				final RealVector coeffs0 = poly.getSubVector(0, 1);
				RealVector coeffs1end = poly.getSubVector(1, poly.getDimension() - 1);
				RealVector bg = xtarget.map(p);

				boolean tooSteep = (d1bg <= minD1 || d1bg > maxD1) || bg.getMaxValue() -
					bg.getMinValue() > tolHigh - tolLow;
				while (tooSteep) {
					coeffs1end = coeffs1end.mapMultiplyToSelf(0.9);
					poly = coeffs0.append(coeffs1end);
					p = new PolynomialFunction(poly.toArray());
					p1 = p.polynomialDerivative();
					d1bg = new Mean().evaluate(xtarget.map(p1).toArray());
					bg = xtarget.map(p);
					tooSteep = (d1bg <= minD1 || d1bg > maxD1) || bg.getMaxValue() - bg
						.getMinValue() > tolHigh - tolLow;
				}

				final boolean tooHigh = bg.getMaxValue() > tolHigh;
				if (tooHigh) {
					final double marginUpper = bg.getMaxValue() - tolHigh;
					coeffs0.setEntry(0, coeffs0.getEntry(0) - (marginUpper < 1e-3
						? 1e-3 : marginUpper));
					poly = coeffs0.append(coeffs1end);
				}

				final boolean tooLow = bg.getMinValue() < tolLow;
				if (tooLow) {
					final double marginLower = tolLow - bg.getMinValue();
					coeffs0.setEntry(0, coeffs0.getEntry(0) + (marginLower < 1e-3
						? 1e-3 : marginLower));
					poly = coeffs0.append(coeffs1end);
				}
			}

			// Gaussian Parameters
			int peakCount = iniSP.getMean().getDimension();
			for (int i = 0; i < peakCount; i++) {
				// Keep means close to maxima or original mean guess`
				if (peakCount == 1) {
					// Do not restrict
				}
//				else if (i == 0) {
//					if (mean.getEntry(i) > mean.getEntry(i + 1))
//						mean.setEntry(i, mean.getEntry(i + 1));
//				}
//				else if (i + 1 == peakCount) {
//					if (mean.getEntry(i) < mean.getEntry(i - 1))
//						mean.setEntry(i, mean.getEntry(i - 1));
//				}
//				else {
//					if (mean.getEntry(i) > mean.getEntry(i + 1))
//						mean.setEntry(i, mean.getEntry(i + 1));
//					if (mean.getEntry(i) < mean.getEntry(i - 1))
//						mean.setEntry(i, mean.getEntry(i - 1));
//				}

				final double diff = mean.getEntry(i) - iniSP.getMean().getEntry(i);
				double sign = 0.0;
				if (diff != 0.0) 
					sign = diff / Math.abs(diff);
				if (Math.abs(diff) > maxMeanDiff.getEntry(i)) 
					mean.setEntry(i, iniSP.getMean().getEntry(i) + sign * maxMeanDiff.getEntry(i));

				// Upper/Lower bound for each parameter
				double ni = norm.getEntry(i);
				double mi = mean.getEntry(i);
				if (mi < minX) {
					mi = minX;
					mean.setEntry(i, mi);
				}
				if (mi > maxX) {
					mi = maxX;
					mean.setEntry(i, mi);
				}
				
				double profv = 0.0;
				if (mi < xtarget.getMinValue()) {
					double profv0 = profile.getPolynomials()[0].value(mi - profile.getKnots()[0]);
					profv = FastMath.max(0.0, profv0);
				} else if (mi > xtarget.getMaxValue()) {
					int n = profile.getPolynomials().length - 1;
					double profv0 = profile.getPolynomials()[n].value(mi - profile.getKnots()[n]);
					profv = FastMath.max(0.0, profv0);
				}	else 
					profv = profile.value(mi);

				double minNi = FastMath.max((profv - p.value(mi)) * minN, 0.0);
				double maxNi = FastMath.max(2.0*minNi, (profv - p.value(mi))*2.0);
				if (ni < minNi)
					norm.setEntry(i, minNi);
				if (ni > maxNi) 
					norm.setEntry(i, maxNi);

				if (sd.getEntry(i) < minSD * iniSP.getSD().getEntry(i))
					sd.setEntry(i, minSD * iniSP.getSD().getEntry(i));
				if (sd.getEntry(i) > maxSD * iniSP.getSD().getEntry(i))
					sd.setEntry(i, maxSD * iniSP.getSD().getEntry(i));
			}

			// Maintain initial AREA proportions between peaks
			if (fitMode == continuumMode) {
				final Variance varCalculator = new Variance();
				final Mean meanCalculator = new Mean();
				area = norm.ebeMultiply(sd).mapMultiply(FastMath.sqrt(2 * FastMath.PI));
				final RealVector ratio = area.ebeDivide(iniSP.getArea());
				double meanRatio = meanCalculator.evaluate(ratio.toArray());
				double varRatio = varCalculator.evaluate(ratio.toArray());
				double lowNorm = meanCalculator.evaluate(
					ytarget.mapSubtract(ytarget.getMinValue()*polyOffset).toArray())*minN;
				double sigma = FastMath.sqrt(FastMath.log(varRatio
									/(meanRatio*meanRatio) + 1));
				if (sigma > areaDrift) {
					// Use a log-normal distribution with parameters mu, sigma
					double mu = FastMath.log(meanRatio /
						FastMath.sqrt(1 + varRatio/(meanRatio*meanRatio)));
					final String outStr = String.format("%1$.4f; %2$.4f; %3$.4f; %4$.4f; %5$.4f; %6$.4f",
						mu, sigma,
						meanRatio, FastMath.exp(mu + 0.5*sigma*sigma), 
						varRatio, (FastMath.exp(sigma*sigma)-1)*FastMath.exp(2*mu+sigma*sigma));
					System.out.println(outStr);
//					double mu = meanRatio;
					sigma = areaDrift;
					LogNormalDistribution logNormal = new LogNormalDistribution(mu, sigma);
					norm = new ArrayRealVector();
					for (int i = 0; i < area.getDimension(); i++) {
						norm = norm.append(logNormal.sample());
					}
					norm = norm.ebeMultiply(iniSP.getNorm());
				}
				double normMin = norm.getMinValue();
				if (normMin < lowNorm) 
					norm.mapMultiply(lowNorm/normMin);

			}

			// Repackage parameter array
			RealVector out = new ArrayRealVector();
			out = out.append(deg).append(poly);
			for (int i = 0; i < iniSP.getMean().getDimension(); i++)
				out = out.append(norm.getEntry(i)).append(mean.getEntry(i)).append(sd
					.getEntry(i));
			return out;
		}
	}
}

class SortedParameters {
	private double[] parameters;
	private RealVector mean;
	private RealVector norm;
	private RealVector sd;
	private RealVector area;
	private final int deg;
	private final RealVector poly; // Polynomial background curve
	
	public SortedParameters(double[] parameters) throws DimensionMismatchException
	{
		this.parameters = parameters;
		this.deg = (int) parameters[0];
		RealVector par = new ArrayRealVector(parameters);
		
		this.mean = new ArrayRealVector();
		this.norm = new ArrayRealVector();
		this.sd   = new ArrayRealVector();
		this.poly = par.getSubVector(1, deg + 1);
		this.area = new ArrayRealVector();
		
		for (int i = deg + 2; i < par.getDimension(); i += 3) {
			norm = norm.append(par.getEntry(i));
			mean = mean.append(par.getEntry(i + 1));
			sd = sd.append(par.getEntry(i + 2));
		}
		area = norm.ebeMultiply(sd).mapMultiply(FastMath.sqrt(2 *
			FastMath.PI));
	}
	
	public SortedParameters(RealVector poly, RealVector norm,
			RealVector mean, RealVector sd)
	{
		if (norm.getDimension() != mean.getDimension())
			throw new DimensionMismatchException(norm.getDimension(), mean.getDimension());
		if (mean.getDimension() != sd.getDimension())
			throw new DimensionMismatchException(mean.getDimension(), sd.getDimension());
		
		this.norm = norm;
		this.mean = mean;
		this.sd   = sd;
		this.poly = poly;
		deg = poly.getDimension() - 1;
		RealVector parArray = new ArrayRealVector().append(deg).append(poly);
		for (int i = 0; i< norm.getDimension(); i++) {
			parArray = parArray.append(norm.getEntry(i));
			parArray = parArray.append(mean.getEntry(i));
			parArray = parArray.append(sd.getEntry(i));
		}
		this.parameters = parArray.toArray();
	}
	
	public int getDeg() {
		return deg;
	}
	public RealVector getMean() {
		return mean;
	}
	public RealVector getNorm() {
		return norm;
	}
	public RealVector getSD() {
		return sd;
	}
	public RealVector getPoly() {
		return poly;
	}
	public RealVector getArea() {
		return area;
	}
	public double[] getParameters() {
		return parameters;
	}
}

class GaussianArray implements UnivariateDifferentiableFunction {
	// Implements a train of Gaussian peaks
	// with an optional polynomial background function

	private final SortedParameters sp;

	public GaussianArray(final SortedParameters sp)
		throws NotStrictlyPositiveException {
		this.sp = sp;
	}

	@Override
	public double value(final double x) {
		double output = 0;
		final int degPoly = (int) sp.getPoly().getEntry(0);
		for (int i = 0; i < sp.getNorm().getDimension(); i++) {
			output += new Gaussian( sp.getNorm().getEntry(i),
															sp.getMean().getEntry(i), 
															sp.getSD().getEntry(i)).value(x);
		}
		output += new PolynomialFunction(sp.getPoly().toArray()).value(x);
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
		System.out.println("Need derivative structure!");
		final PolynomialFunction p = new PolynomialFunction(sp.getPoly().toArray());
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
			SortedParameters sp = new SortedParameters(param);
			return new GaussianArray(sp).value(x);
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
