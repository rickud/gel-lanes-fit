/**
 * Gauss Fit
 * GelLanesFit.java
 * author: Rick Ziraldo, 2017
 * The /University of Texas at Dallas, Richardson, TX
 * http://www.utdallas.edu
 *
 * Feature: Fitting of multiple Gaussian functions to intensity profiles along the gel lanes
 * Gauss Fit is a tool for fitting gaussian profiles and estimating
 * the profile parameters on selected lanes in gel electrophoresis images.
 *
 *    The GaussianArrayCurveFitter class is implemented using
 *    Abstract classes from Apache Commons project
 *
 *    The source code is maintained and made available on GitHub
 *    https://github.com/rickud/gauss-curve-fit
 */

package gellanesfit;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import net.imagej.table.DefaultGenericTable;
import net.imagej.table.DefaultTableDisplay;
import net.imagej.table.GenericColumn;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.analysis.integration.gauss.GaussIntegrator;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;
import org.scijava.Context;
import org.scijava.app.StatusService;
import org.scijava.display.Display;
import org.scijava.display.DisplayService;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;

import gellanesfit.GaussianArrayCurveFitter.ParameterGuesser;

class Fitter {

	@Parameter
	private LogService log;
	@Parameter
	private DisplayService displayServ;
	@Parameter
	private StatusService statusServ;

	static final int regMode = 0;
	static final int fragmentMode = 1;
	private static final int noLadderLane = -1;

	public static final double peakDistanceTol = 2;
	public static final double sd2FWHM = 2 * FastMath.sqrt(2 * FastMath.log(2));

	private int degBG;
	private double tolPK;
	private int fitMode = regMode;
	private int ladderLane = noLadderLane;
	private double[] ladder; // MW, Dalton
	private double[][] fragmentDistribution; // MW, Dalton; relative frequency
	private final String title; // for the datafile, main image's title
	private List<DataSeries> inputData;

	private List<Peak> allGuessList;
	private List<Peak> allFittedList;
	private List<Peak> allCustomList;

	public Fitter(final Context context, final String title, final int degBG, final double tolPK) {
		context.inject(this);
		this.title = title;
		this.degBG = degBG;
		this.tolPK = tolPK;
		this.inputData = new ArrayList<>();

		this.allGuessList = new ArrayList<>();
		this.allFittedList = new ArrayList<>();
		this.allCustomList = new ArrayList<>();
	}

	private List<Peak> arrayToPeaks(final int ln, final double[] param) {
		final int gaussStart = (int) param[0] + 2;
		final List<Peak> peaks = new ArrayList<>();
		for (int p = gaussStart; p < param.length; p += 3) {
			peaks.add(new Peak(ln, param[p], param[p + 1], param[p + 2]));
		}
		return peaks;
	}

	private RealVector peaksToArray(final List<Peak> peaks) {
		RealVector guessArray = new ArrayRealVector();
		for (final Peak p : peaks) {
			guessArray = guessArray.append(p.getNorm());
			guessArray = guessArray.append(p.getMean());
			guessArray = guessArray.append(p.getSigma());
		}
		return guessArray;
	}

	public void updateResultsTable() {
		// Results Table Columns
		String[] headersB = { "", "Band", "Distance", "Dist. G.", "Amplitude", "Amp. G.", "FWHM", "FWHM G.", "Area"};
		String[] headersF = { "", "Band", "Distance", "Dist. G.", "Amplitude", "Amp. G.", "FWHM", "FWHM G.", "Area", "Frequency", "BP", "MW" };
		String[] headers = null;
		if (fitMode == regMode) headers = headersB;
		else headers = headersF;
		
		final GenericColumn[] tableCol = new GenericColumn[headers.length];
		final DefaultGenericTable rt = new DefaultGenericTable();
		for (int cc = 0; cc < headers.length; cc++)
			tableCol[cc] = new GenericColumn(headers[cc]);

		for (final DataSeries d : inputData) {
			final int lane = d.getLane();
			final List<Peak> guess = getGuessPeaks(lane);
			final List<Peak> fitted = getFittedPeaks(lane);
			if (guess.size() != fitted.size()) {
				log.error("Size mismatch: " + guess.size() + ", " + fitted.size());
			}
			int band = 1;
			for (int p = 0; p < guess.size(); p++) {
				final double n = fitted.get(p).getNorm();
				final double m = fitted.get(p).getMean();
				final double s = fitted.get(p).getSigma();
				final double a = doIntegrate(d.getX(), n, m, s);

				final double n0 = guess.get(p).getNorm();
				final double m0 = guess.get(p).getMean();
				final double s0 = guess.get(p).getSigma();

				// Prepare columns for Results Table
				if (band == 1) {
					tableCol[0].add("Lane " + lane);
				} else {
					tableCol[0].add("");
				}
				tableCol[1].add(band);
				tableCol[2].add(String.format("%1$.1f", m));
				tableCol[3].add(String.format("%1$.1f", m0));
				tableCol[4].add(String.format("%1$.1f", n));
				tableCol[5].add(String.format("%1$.1f", n0));
				tableCol[6].add(String.format("%1$.2f", s * sd2FWHM));
				tableCol[7].add(String.format("%1$.2f", s0 * sd2FWHM));
				tableCol[8].add(String.format("%1$.1f", a));

				if (fitMode == fragmentMode) {
					if (lane == ladderLane) { 
						tableCol[9].add(" - ");
						tableCol[10].add(" - ");
						tableCol[11].add(" - ");
					} else {
						tableCol[9].add(String.format("%1$.3f", fragmentDistribution[band-1][0]));
						tableCol[10].add(String.format("%1$.1f", fragmentDistribution[band-1][1]));
						tableCol[11].add(String.format("%1$.3f", fragmentDistribution[band-1][2]));
					}	
				}
				band++;
			}
		}

		for (int cc = 0; cc < headers.length; cc++)
			rt.add(tableCol[cc]);

		// Update results table
		DefaultTableDisplay rtDisplay = null;
		for (final Display<?> d : displayServ.getDisplays()) {
			//log.info(d.getClass().toString());
			if (d instanceof DefaultTableDisplay) {
				rtDisplay = (DefaultTableDisplay) d;
				rtDisplay.clear();
				rtDisplay.add(rt);
				rtDisplay.update();
			}
		}
		if (rtDisplay == null) {
			rtDisplay = (DefaultTableDisplay) displayServ.createDisplay("Results Display", rt);
		}
		displayServ.setActiveDisplay(rtDisplay);
		
		// Save Results Table
		final String path = "gel-lanes-fit/data/";
		final String file = "Fit of " + title + ".xls";
		new File(path).mkdirs();
		log.info("Saving to " + path + file + " ...");
		try (BufferedWriter out = new BufferedWriter(new FileWriter(path + file))) {
			String outText = "";
			for (int cc = 0; cc < headers.length; cc++)
				outText += headers[cc] + "\t";
			outText += "\n";
			for (int rr = 0; rr < tableCol[0].size(); rr++) {
				for (int cc = 0; cc < headers.length; cc++) {
					outText += tableCol[cc].get(rr) + "\t";
				}
				outText = outText.substring(0, outText.length() - 1) + "\n";
			}
			out.write(outText);
			out.close();
		} catch (final IOException e) {
			log.info("Exception", e);
		}
	}

	private double[] interpolateDisplacement(final double[] y) {
		final WeightedObservedPoints obs = new WeightedObservedPoints();
		final double[] logs = new double[y.length];
		final double distLogs[] = new double[fragmentDistribution.length];
		for (int l = 0; l < y.length; l++) {
			logs[l] = Math.log10(ladder[l]);
			obs.add(Math.log10(ladder[l]), y[l]);
		}
		// First-degree polynomial fitter (line).
		final PolynomialCurveFitter fitter = PolynomialCurveFitter.create(1);
		final double[] coeff = fitter.fit(obs.toList());
		final double[] yi = new double[fragmentDistribution.length];
		final UnivariateFunction f = new PolynomialFunction(coeff);
		for (int i = 0; i < fragmentDistribution.length; i++) {
			distLogs[i] = Math.log10(fragmentDistribution[i][2]);
			yi[i] = f.value(Math.log10(fragmentDistribution[i][2]));
		}
		return yi;
	}

	private double[] interpolateSD(final double[] y, final double[] sd, final double[] yi) {
		final double[] sdi = new double[yi.length];
		final PolynomialSplineFunction f = new LinearInterpolator().interpolate(y, sd);
		for (int i = 0; i < yi.length; i++) {
			if (yi[i] < y[0])
				sdi[i] = sd[0];
			else if (yi[i] > y[y.length - 1])
				sdi[i] = sd[sd.length - 1];
			else
				sdi[i] = f.value(yi[i]);
		}
		return sdi;
	}

	private RealVector doGuess(final int lane, final ParameterGuesser pg) {
		// Remove guess and fit for this lane
		final Iterator<Peak> peakIter = allGuessList.iterator();
		while (peakIter.hasNext()) {
			if (peakIter.next().getLane() == lane)
				peakIter.remove();
		}

		final Iterator<Peak> peakIter2 = allFittedList.iterator();
		while (peakIter2.hasNext()) {
			if (peakIter2.next().getLane() == lane)
				peakIter2.remove();
		}

		List<Peak> peaks = new ArrayList<>();
		RealVector guess = new ArrayRealVector();
		if (fitMode == regMode || (fitMode == fragmentMode && lane == ladderLane)) {
			// Guess a set of peaks a set of peaks automatically
			guess = new ArrayRealVector(pg.guess());
			peaks = arrayToPeaks(lane, guess.toArray());

		} else if (fitMode == fragmentMode) {
			if (ladder == null || fragmentDistribution == null) {
				log.info("Missing ladder/distribution");
				return null;
			}
			// Use the stored distribution as a guess 
			// fragmentDistribution[:][0] = Frequency
			// fragmentDistribution[:][1] = Length
			// fragmentDistribution[:][2] = MW
			
			final List<Peak> ladderPeaks = getFittedPeaks(ladderLane);
			final double[] meanLadder = new double[ladderPeaks.size()];
			final double[] sdLadder = new double[ladderPeaks.size()];
			for (int p = 0; p < ladderPeaks.size(); p++) {
				meanLadder[p] = ladderPeaks.get(p).getMean();
				sdLadder[p] = ladderPeaks.get(p).getSigma();
			}
			final double[] means = interpolateDisplacement(meanLadder);
			final double[] sds = interpolateSD(meanLadder, sdLadder, means);
			final double[] norms = new double[means.length]; // range
			RealVector scaledFrequency = new ArrayRealVector(new Array2DRowRealMatrix(fragmentDistribution).getColumn(0));
			scaledFrequency = scaledFrequency.mapMultiplyToSelf(1 / scaledFrequency.getMaxValue());
			RealVector mwArray = new ArrayRealVector(new Array2DRowRealMatrix(fragmentDistribution).getColumn(2));
			for (final DataSeries d : inputData) {
				if (d.getLane() == lane) {
					RealVector profile = d.getY().mapSubtractToSelf(d.getMinY());
//					final double[] y = d.getX().toArray();
//					final LinearInterpolator li = new LinearInterpolator();
//					final UnivariateFunction f = li.interpolate(y, profile);
					for (int i = 0; i < means.length; i++) {
						// exclude peaks outside the profile domain
						if (means[i] > d.getMinX() && means[i] < d.getMaxX()) {
							double mw = mwArray.getEntry(i);
							double frequency = scaledFrequency.getEntry(i);
							double scale = profile.getMaxValue()/(FastMath.log10(mwArray.getMaxValue())*scaledFrequency.getMaxValue());
							norms[i] = scale * FastMath.log10(mw) * frequency;
							peaks.add(new Peak(lane, norms[i], means[i], sds[i]));
						}
					}
				}
			}
		}

		// Check if custom peaks are needed
		for (final Peak c : getCustomPeaks(lane)) {
			// Check if guessList contains peak and update
			boolean found = false;
			for (final Peak g : peaks) {
				if (FastMath.abs(c.getMean() - g.getMean()) <= peakDistanceTol) {
					found = true;
					log.info("Peak guess:    " + g.getLane() + ", "+ g.getNorm() + ", " + g.getMean() + ", "
					        + g.getSigma());
					g.setMean(c.getMean());
					g.setNorm(c.getNorm());
					g.setSigma(c.getSigma());
					log.info("\tReplaced with: " + c.getNorm() + ", " + c.getMean() + ", "
					        + c.getSigma());
					break;
				}
			}
			if (!found) {
				peaks.add(c);
				Collections.sort(peaks);
			}
		}
		allGuessList.addAll(peaks);
		RealVector poly = new ArrayRealVector(); 
		poly = poly.append(degBG).append(new ArrayRealVector(degBG + 1));
		if (fitMode == regMode) {
			poly.setEntry(1, guess.getEntry(1));
		} else if (fitMode == fragmentMode) {
			for (DataSeries d : inputData) {
				if (d.getLane() == lane) 
					poly.setEntry(1, d.getMinY());
			}
		}
		return poly.append(peaksToArray(peaks));
	}

	private double doIntegrate(final RealVector xvals, final double n, final double m,
	        final double s) {
		double area;
		RealVector weights = new ArrayRealVector(xvals.getDimension()).mapAddToSelf(1.0);
		final GaussIntegrator ti = new GaussIntegrator(xvals.toArray(), weights.toArray());
		final Gaussian gauss = new Gaussian(n, m, s);
		area = ti.integrate(gauss);
		return area;
	}

	public List<DataSeries> doFit() {
		int progress = 1;
		final ArrayList<DataSeries> out = new ArrayList<>();
		for (final DataSeries d : inputData) {
			out.addAll(doFit(d.getLane()));
		}
		statusServ.showProgress(++progress, inputData.size());
		return out;
	}

	public List<DataSeries> doFit(final int lane) {
		final List<DataSeries> out = new ArrayList<>();
		for (final DataSeries d : inputData) {
			if (d.getLane() == lane) {
				out.addAll(doFit(d));
			}
		}
		return out;
	}

	private List<DataSeries> doFit(final DataSeries in) {
		final int lane = in.getLane();
		
		int oldFitMode = fitMode;
		if (lane == ladderLane) fitMode = Fitter.regMode;
		
		final RealVector xvals = in.getX();
		final RealVector yvals = in.getY();

		final List<Peak> fittedPeaks = new ArrayList<>();
		final List<DataSeries> output = new ArrayList<>();

		// Tolerance as percentage of the range
		final double tolpk = tolPK * (yvals.getMaxValue() - yvals.getMinValue());

		final WeightedObservedPoints obs = new WeightedObservedPoints();
		for (int o = 0; o < xvals.getDimension(); o++)
			obs.add(xvals.getEntry(o), yvals.getEntry(o));

		final ParameterGuesser pg = new GaussianArrayCurveFitter.ParameterGuesser(obs.toList(),
		        tolpk, degBG);

		final double[] firstGuess = doGuess(lane, pg).toArray();
		final LeastSquaresProblem problem = GaussianArrayCurveFitter.create(tolpk, degBG, fitMode)
		        .withStartPoint(firstGuess).getProblem(obs.toList());

		final LeastSquaresOptimizer.Optimum optimum = new LevenbergMarquardtOptimizer()
		        .optimize(problem);

		final RealVector fitted = new ArrayRealVector(optimum.getPoint());

		// After fitting
		RealVector norms = new ArrayRealVector();
		RealVector means = new ArrayRealVector();
		RealVector sds = new ArrayRealVector();
		final RealVector poly = fitted.getSubVector(0, degBG + 2);

		for (int b = degBG + 2; b < fitted.getDimension(); b += 3) {
			// Fitted Peaks
			final double n = fitted.getEntry(b);
			final double m = fitted.getEntry(b + 1);
			final double s = fitted.getEntry(b + 2);

			norms = norms.append(n);
			means = means.append(m);
			sds = sds.append(s);
			fittedPeaks.add(new Peak(lane, n, m, s));
		}
		allFittedList.addAll(fittedPeaks);

		// Create new DataSeries for Plotter
		final PolynomialFunction bg = new PolynomialFunction(
		        poly.getSubVector(1, degBG + 1).toArray());
		output.add(new DataSeries("Background", lane, DataSeries.BACKGROUND, xvals, bg,
		        Plotter.bgColor));

		for (final Peak p : fittedPeaks) {
			final Gaussian gauss = new Gaussian(p.getNorm(), p.getMean(), p.getSigma());
			final UnivariateFunction[] functs = { bg, gauss };
			output.add(new DataSeries("Band " + 1, lane, DataSeries.GAUSS_BG, xvals, functs,
			        Plotter.gaussColor));
		}
		final GaussianArray fittedCurve = new GaussianArray(norms, means, sds, poly);
		output.add(new DataSeries("Fit", lane, DataSeries.FITTED, xvals, fittedCurve,
		        Plotter.fittedColor));

		// Print RMS to console
		final String outStr = String.format("Lane " + lane + ", RMS: %1$.2f; ", optimum.getRMS());
		log.info(outStr);
		fitMode = oldFitMode;
		return output;
	}

	public void addCustomPeak(final Peak peak) {
		boolean found = false;
		// If close to existing replace
		for (final Peak p : allCustomList) {
			if (p.getLane() == peak.getLane()) {
				if (FastMath.abs(peak.getMean() - p.getMean()) <= Fitter.peakDistanceTol) {
					p.setSigma(peak.getSigma());
					found = true;
					break;
				}
			}
		}
		if (!found)
			allCustomList.add(peak);
		Collections.sort(allCustomList);
	}

	public boolean removeCustomPeak(final Peak peak) {
		// Remove from custom list
		final Iterator<Peak> peakIter = allCustomList.iterator();
		while (peakIter.hasNext()) {
			final Peak p = peakIter.next();
			if (p.getLane() == peak.getLane()) {
				final double diff = FastMath.abs(peak.getMean() - p.getMean());
				if (diff <= Fitter.peakDistanceTol) {
					peakIter.remove();
					return true;
				}
			}
		}
		return false;
	}

	public void resetCustomPeaks(final int lane) {
		final Iterator<Peak> peakIter = allCustomList.iterator();
		while (peakIter.hasNext()) {
			if (peakIter.next().getLane() == lane)
				peakIter.remove();
		}
	}

	public void resetFit(final int lane) {
		final Iterator<Peak> itFitted = allFittedList.iterator();
		while (itFitted.hasNext()) {
			if (itFitted.next().getLane() == lane)
				itFitted.remove();
		}
		final Iterator<Peak> itGuess = allGuessList.iterator();
		while (itGuess.hasNext()) {
			if (itGuess.next().getLane() == lane)
				itGuess.remove();
		}
	}

	public void resetAllFitter() {
		inputData = new ArrayList<>();
		allGuessList = new ArrayList<>();
		allFittedList = new ArrayList<>();
		allCustomList = new ArrayList<>();
	}

	public List<Peak> getGuessPeaks(final int lane) {
		final List<Peak> g = new ArrayList<>();
		for (final Peak p : allGuessList) {
			if (p.getLane() == lane)
				g.add(p);
		}
		return g;
	}

	public List<Peak> getCustomPeaks(final int lane) {
		final List<Peak> c = new ArrayList<>();
		for (final Peak p : allCustomList) {
			if (p.getLane() == lane)
				c.add(p);
		}
		return c;
	}

	public List<Peak> getFittedPeaks(final int lane) {
		final List<Peak> f = new ArrayList<>();
		for (final Peak p : allFittedList) {
			if (p.getLane() == lane)
				f.add(p);
		}
		return f;
	}

	/**
	 * @param degBG
	 */
	public void setDegBG(final int degBG) {
		this.degBG = degBG;
	}

	public void setInputData(final ArrayList<DataSeries> inputData) {
		this.inputData = inputData;
	}

	public void setFitMode(final int fitMode) {
		this.fitMode = fitMode;
	}

	public void setReferenceLane(final int ladderLane) {
		this.ladderLane = ladderLane;
	}

	public void setFragmentDistribution(final double[][] fragmentDistribution) {
		this.fragmentDistribution = fragmentDistribution;
	}

	/**
	 * @param tolPK
	 */
	public void setTolPK(final double tolPK) {
		this.tolPK = tolPK;
	}

	public void setLadder(final RealVector ladder) {
		this.ladder = ladder.toArray();
	}
}

class Peak implements Comparable<Peak> {

	private String name = "";
	private final int lane;
	private double mean;
	private double norm;
	private double sd;

	public Peak(final int lane, final double norm, final double mean) {
		this.lane = lane;
		this.norm = norm;
		this.mean = mean;
		this.sd = 0.0;
	}

	public Peak(final int lane, final double norm, final double mean, final double sd) {
		this.lane = lane;
		this.norm = norm;
		this.mean = mean;
		if (sd > 0) {
			this.sd = sd;
		} else {
			throw new NotStrictlyPositiveException(sd);
		}
	}

	public double getNorm() {
		return norm;
	}

	public double getMean() {
		return mean;
	}

	public double getSigma() {
		return sd;
	}

	public int getLane() {
		return lane;
	}

	public String getName() {
		return name;
	}

	public void setMean(final double mean) {
		this.mean = mean;
	}

	public void setName(final String name) {
		this.name = name;
	}

	public void setNorm(final double norm) {
		this.norm = norm;
	}

	public void setSigma(final double sd) {
		if (sd > 0) {
			this.sd = sd;
		} else {
			throw new NotStrictlyPositiveException(sd);
		}
	}

	@Override
	public int compareTo(final Peak p) {
		final double m = mean - p.getMean();
		final int l = lane - p.getLane();

		if (l == 0) {
			return (int) m;
		}
		return l;
	}
}
