/*
 * Gel Lanes Fit - ParameterGuesser.java
 * Author: Rick Ziraldo, 2017
 * The University of Texas at Dallas, Richardson, TX
 *
 * Licensed under the GNU Affero General Public License v3.0; see LICENSE.
 * Source: https://github.com/rickud/gel-lanes-fit
 */

package gellanesfit;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.function.Log;
import org.apache.commons.math3.analysis.function.Log10;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.exception.NumberIsTooSmallException;
import org.apache.commons.math3.exception.util.LocalizedFormats;
import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;

/**
 * Starting guesses for the fits: the bands found in a lane's profile (Banded),
 * or one peak per fragment length of a distribution, placed and sized from
 * the fitted ladder lane (Continuum).
 * <p>
 * A guess is a {@link SortedParameters}: the background polynomial's degree and
 * coefficients, then a (height, position, standard deviation) triplet per peak.
 * </p>
 */
final class ParameterGuesser {

	static final int bandMode = 0;
	static final int continuumMode = 1;

	/**
	 * Continuum: fragments predicted up to this fraction of the lane length
	 * beyond either end are fitted too, so edge bands fit properly
	 */
	private static final double CONTINUUM_MARGIN = 0.2;
	/**
	 * Starting band widths: the half-maximum level is raised by this factor
	 * while the width found is wider than the gap to the next band, but not
	 * above HALF_MAX_LIMIT of the band's height
	 */
	private static final double HALF_MAX_STEP = 1.05;
	private static final double HALF_MAX_LIMIT = 0.9;

	/**
	 * A profile's value at position x, interpolated linearly. Beyond the lane,
	 * it's extrapolated from the first or last segment, and never negative.
	 */
	static double profileAt(final PolynomialSplineFunction profile,
		final double x)
	{
		final double[] knots = profile.getKnots();
		final PolynomialFunction[] segments = profile.getPolynomials();
		if (x < knots[0]) return FastMath.max(0.0, segments[0].value(x -
			knots[0]));
		if (x > knots[knots.length - 1]) {
			final int n = segments.length - 1;
			return FastMath.max(0.0, segments[n].value(x - knots[n]));
		}
		return profile.value(x);
	}
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
		fitMode = bandMode;
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
		fitMode = continuumMode;
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
			new Comparator<>()
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
		
		// A flat background at the lane's minimum. Fitter replaces it with its
		// own; it's only used when a fit is started without a start point
		final RealVector polyG = new ArrayRealVector(deg + 1);
		if (deg >= 0) polyG.setEntry(0, yvals.getMinValue());

		if (fitMode == bandMode) return guessBands(points, xvals, yvals, polyG,
			tolpk);
		if (fitMode == continuumMode) return guessFragments(xvals, yvals, polyG);
		return new SortedParameters(polyG, new ArrayRealVector(),
			new ArrayRealVector(), new ArrayRealVector());
	}

	/**
	 * Banded starting guess: a Gaussian at each band found by
	 * {@link #findMaxima}, its height measured from the lane's minimum and
	 * its width from where the band falls to half its height.
	 */
	private SortedParameters guessBands(final WeightedObservedPoint[] points,
		final RealVector xvals, final RealVector yvals, final RealVector polyG,
		final double tolpk)
	{
		final double minX = xvals.getMinValue();
		final double maxX = xvals.getMaxValue();
		final double minY = yvals.getMinValue();
		RealVector normG = new ArrayRealVector();
		RealVector meanG = new ArrayRealVector();
		RealVector sdG = new ArrayRealVector();
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

			final double inc = HALF_MAX_STEP;
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

			while ((!foundFWHM || FWHM > peakDistance) && hm * inc < HALF_MAX_LIMIT * yRange) {
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
			sds[m] = FWHM / (2 * FastMath.sqrt(2 * FastMath.log(2)));
			normG = normG.append(norms[m]);
			meanG = meanG.append(means[m]);
			sdG   = sdG.append(sds[m]);
		}
		return new SortedParameters(polyG, normG, meanG, sdG);
	}

	/**
	 * Continuum starting guess: a Gaussian for each fragment of the
	 * distribution predicted to run within the lane (or within
	 * {@link #CONTINUUM_MARGIN} of it), at the position and width predicted
	 * from the ladder, with a height from the profile at that position scaled
	 * by the fragment's frequency.
	 */
	private SortedParameters guessFragments(final RealVector xvals,
		final RealVector yvals, final RealVector polyG)
	{
		RealVector normG = new ArrayRealVector();
		RealVector meanG = new ArrayRealVector();
		RealVector sdG = new ArrayRealVector();
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
			final double pv = profileAt(pr, means[s]);
			
			scaledFrequency.setEntry(s, scaledFrequency.getEntry(s) * pv);
		}
		final RealVector mwArray = distMatrix.getColumnVector(2);
		RealVector scale = scaledFrequency.ebeMultiply(mwArray.map(new Log()));
		scale = scale.mapDivide(scale.getMaxValue());
	
		List<Integer> fragmentSubset = new ArrayList<>();
		double margin = (xvals.getMaxValue() - xvals.getMinValue()) *
			CONTINUUM_MARGIN;
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
		return new SortedParameters(polyG, normG, meanG, sdG);
	}


	/**
	 * Adapted From: Calculates peak positions of 1D array N.Vischer,
	 * 13-sep-2013
	 *
	 * @param x Array containing peaks.
	 * @param depth Depth of a qualified valley must exceed depth.
	 *          Tolerance must be >= 0. Flat tops are marked at their centers.
	 * @param includeEnds If 'false', a peak is only accepted if it is separated
	 *          by two qualified valleys. If 'true', a peak is also accepted if
	 *          separated by one qualified valley and by a border.
	 * @return Positions of peaks, sorted with decreasing amplitude
	 */
	private int[] findMaxima(final WeightedObservedPoint[] x, final double depth,
		final boolean includeEnds)
	{
		final int len = x.length;
		if (len < 2) return new int[0];
		final double tolerance = FastMath.max(depth, 0);
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


	/**
	 * Predicts where fragments run: migration distance is interpolated linearly
	 * against log10 of the molecular weight between the ladder bands, and
	 * extrapolated from the first and last segments beyond them.
	 * <p>
	 * Note: reverses {@code y} in place. {@link #interpolateSD} relies on this
	 * when it is called after it with the same array.
	 * </p>
	 * @param y positions of the ladder bands, in pixels, top band first
	 * @param ladder molecular weights of the ladder bands, in Da
	 * @param dist molecular weights of the fragments, in Da
	 * @return the predicted position of each fragment, in pixels
	 */
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

	/**
	 * Predicts fragment widths from a straight line fitted to the ladder bands'
	 * widths against their positions. Beyond the ladder, the width stays at
	 * the line's value at the nearest band: extended further, the line can
	 * reach zero or below.
	 * @param y ladder band positions, as left reversed by
	 *          {@link #interpolateDisplacement}
	 * @param sd ladder band standard deviations, top band first
	 * @param yi positions of the fragments
	 * @return the predicted standard deviation of each fragment
	 */
	private double[] interpolateSD(final double[] y, final double[] sd,
		final double[] yi)
	{
		final WeightedObservedPoints obs = new WeightedObservedPoints();
		for (int l = 0; l < y.length; l++)
			obs.add(y[y.length - 1 - l], sd[l]);

		// First-degree polynomial fitter (line)
		final PolynomialCurveFitter linfit = PolynomialCurveFitter.create(1);
		final double[] coeffs = linfit.fit(obs.toList());
		final UnivariateFunction line = new PolynomialFunction(coeffs);
		final double first = FastMath.min(y[0], y[y.length - 1]);
		final double last = FastMath.max(y[0], y[y.length - 1]);
		final UnivariateFunction f = x -> line.value(FastMath.max(first, FastMath
			.min(last, x)));
		final RealVector sdi = new ArrayRealVector(yi).map(f);
		return sdi.toArray();
	}

	/**
	 * Indexes of the distribution's fragments in the guess for a lane with the
	 * given range of positions: those predicted to run within it, or within
	 * {@link #CONTINUUM_MARGIN} of its length beyond either end, as in
	 * {@link #guessFragments}. The fitted peaks are in the same order.
	 */
	public List<Integer> getUsedFragments(double[] xrange) {
		final double[] meanLadder = new double[ladderPeaks.size()];
		for (int p = 0; p < ladderPeaks.size(); p++) {
			meanLadder[p] = ladderPeaks.get(p).getMean();
		}
		List<Integer> out = new ArrayList<>(); 
		final double[] means = interpolateDisplacement(meanLadder, ladderMW, distMatrix.getColumnVector(2));
		final double margin = (xrange[1] - xrange[0]) * CONTINUUM_MARGIN;
		for (int i = 0; i < means.length; i++) {
			if (means[i] > xrange[0] - margin && means[i] < xrange[1] + margin) out.add(i);
		}
		return out;
	}
}
