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

import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.analysis.function.Log10;
import org.apache.commons.math3.analysis.integration.gauss.GaussIntegrator;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
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
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.util.FastMath;
import org.scijava.Context;
import org.scijava.app.StatusService;
import org.scijava.display.DisplayService;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;

import gellanesfit.GaussianArrayCurveFitter.ParameterGuesser;
import ij.measure.ResultsTable;


class Fitter {

	@Parameter
	private LogService log;
	@Parameter
	private DisplayService displayServ;
	@Parameter
	private StatusService statusServ;

	static final int bandMode = 0;
	static final int continuumMode = 1;

	public static final double peakDistanceTol = 2;
	public static final double sd2FWHM = 2 * FastMath.sqrt(2 * FastMath.log(2));

	private int degBG;
	private double tolPK;
	private double areaDrift;
	private double polyDerivative;
	private int fitMode = bandMode;
	private int ladderLane = MainDialog.noLadderLane;
	private double[] ladder; // MW, Dalton
	private double[][] fragmentDistribution; // MW, Dalton; relative frequency
	private List<List<Integer>> selectedFragments;
	private final String title; // for the datafile, main image's title
	private List<DataSeries> inputData;

	private List<Peak> allGuessList;
	private List<Peak> allFittedList;
	private List<Peak> allCustomList;
	
	// Fields used in the summary
	private RealVector rms;
	private Array2DRowRealMatrix statMatrix;
	private List<List<Double>> fittedDistributions;
	
	public Fitter(final Context context, final String title ) {
		context.inject(this);
		this.title = title;
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

	public void updateResultsTable(String savePath) {
		// Results Table Columns
		String[] headersB = { "Lane", "Band", "Distance", "Dist. G.", "Amplitude", "Amp. G.", "FWHM", "FWHM G.", "Area"};
		String[] headersF = { "Lane", "Band", "Distance", "Dist. G.", "Amplitude", "Amp. G.", "FWHM", "FWHM G.", "Area", "Frequency", "BP", "MW" };
		String[] headers = null;
		if (fitMode == bandMode) headers = headersB;
		else headers = headersF;
		final ResultsTable rt = new ResultsTable();
		
		for (final DataSeries d : inputData) {
			final int lane = d.getLane();
			final List<Peak> guess = getGuessPeaks(lane);
			final List<Peak> fitted = getFittedPeaks(lane);
			RealVector areas = new ArrayRealVector();
			if (guess.size() != fitted.size()) {
				log.error("Size mismatch: " + guess.size() + ", " + fitted.size());
			}
			
			int band = 1;
			double bpamax = Double.NEGATIVE_INFINITY;
			double[] scaledCount = new double[guess.size()];
			for (int p = 0; p < guess.size(); p++) {
				final double n = fitted.get(p).getNorm();
				final double m = fitted.get(p).getMean();
				final double s = fitted.get(p).getSigma();
				final double a = doIntegrate(d.getX(), n, m, s);
				
				final double n0 = guess.get(p).getNorm();
				final double m0 = guess.get(p).getMean();
				final double s0 = guess.get(p).getSigma();
				areas = areas.append(a);
				
				// Add values to Results Table
				rt.incrementCounter();
				if (band == 1) {
					rt.addValue(headers[0], "" + lane);
				} else {
					rt.addValue(headers[0], "");
				}
				rt.addValue(headers[1], band);
				rt.addValue(headers[2], String.format("%1$.1f", m));
				rt.addValue(headers[3], String.format("%1$.1f", m0));
				rt.addValue(headers[4], String.format("%1$.1f", n));
				rt.addValue(headers[5], String.format("%1$.1f", n0));
				rt.addValue(headers[6], String.format("%1$.2f", s * sd2FWHM));
				rt.addValue(headers[7], String.format("%1$.2f", s0 * sd2FWHM));
				rt.addValue(headers[8], String.format("%1$.1f", a));
				
				if (fitMode == continuumMode) {
					if (lane == ladderLane) { 
						rt.addValue(headers[9], " - ");
						rt.addValue(headers[10], " - ");
						rt.addValue(headers[11], " - ");
					} else {
						double freq = fragmentDistribution[selectedFragments.get(lane-1).get(p)][0];
						double bp = fragmentDistribution[selectedFragments.get(lane-1).get(p)][1];
						double mw = fragmentDistribution[selectedFragments.get(lane-1).get(p)][2];
						rt.addValue(headers[9],  String.format("%1$.3f", freq));
						rt.addValue(headers[10], String.format("%1$.0f", bp));
						rt.addValue(headers[11], String.format("%1$.3f", mw));
						
						scaledCount[p] = bp * a * 1000;
						if (bp * a > bpamax)
							bpamax = bp * a;
					}
				}
				band++;
			}
			
			if (fitMode == continuumMode && lane != ladderLane) {
				RealVector bpSubarray = new ArrayRealVector();
				RealMatrix distMatrix = new Array2DRowRealMatrix(fragmentDistribution);
				RealVector bpArray = distMatrix.getColumnVector(1);
				bpSubarray = new ArrayRealVector();
				for (int i : selectedFragments.get(lane-1)) 
					bpSubarray = bpSubarray.append(bpArray.getEntry(i));
				double meanLength = new Mean().evaluate(bpSubarray.toArray());
				RealVector scaleFactor = areas.ebeMultiply(bpSubarray);
				double weighedMeanLength = new Mean().evaluate(bpSubarray.toArray(), scaleFactor.toArray());
				double sdLength = FastMath.sqrt(
									new Variance()
									.evaluate(bpSubarray.toArray(), scaleFactor.toArray()));
				statMatrix.setRow(lane-1, new double[] {weighedMeanLength, sdLength});
				String info = String.format("Lane %1$d, Average length: %2$.2f(%3$.2f) %4$.2f", 
					lane, weighedMeanLength, sdLength, meanLength);
				
				log.info(info);
				fittedDistributions.set(lane-1, new ArrayList<>());
				for (int i = 0; i < bpSubarray.getDimension(); i++) {
					int reps = (int) (1000 * scaleFactor.getEntry(i) / scaleFactor.getL1Norm());
					for (int j = 0; j < reps; j++)
						fittedDistributions.get(lane-1).add(bpSubarray.getEntry(i));
				}
			}
		}	
		rt.show("Results Display");
		
		// Save Results Table
		String outText = "";
		for (int cc = 0; cc < headers.length; cc++)
			outText += headers[cc] + "\t";
		outText += "\n";
		for (int rr = 0; rr < rt.getCounter(); rr++) {
			for (int cc = 0; cc < headers.length; cc++) {
				outText += rt.getStringValue(cc, rr) + "\t";
			}
			outText = outText.substring(0, outText.length() - 1) + "\n";
		}
		if (fitMode == continuumMode) {
			String line = "\n\tMean\tStandard Deviation\n";
			for (int i = 0; i < inputData.size(); i++) {
				if (i !=  ladderLane - 1) {
					line = String.format("Lane %1$d\t%2$.2f \t %3$.2f\n", 
						(i + 1), statMatrix.getEntry(i, 0), statMatrix.getEntry(i, 1));
					outText += line;
				}
			}
		}
		
		final String file = "Fit of " + title + ".xls";
		new File(savePath).mkdirs();
		log.info("Saving to " + savePath + file + " ...");
		try (BufferedWriter out = new BufferedWriter(new FileWriter(savePath + file))) {
			out.write(outText);
			out.close();
		} catch (final IOException e) {
			log.info("Exception", e);
		}
	}

	private double[] interpolateDisplacement(final double[] y) {
		final WeightedObservedPoints obs = new WeightedObservedPoints();
		RealVector logs = new ArrayRealVector(new Array2DRowRealMatrix(fragmentDistribution).getColumn(2));
		for (int l = 0; l < y.length; l++)
			obs.add(Math.log10(ladder[l]), y[l]);

		// First-degree polynomial fitter (line)
		final PolynomialCurveFitter linfit = PolynomialCurveFitter.create(1);
		final double[] coeffs = linfit.fit(obs.toList());
		final UnivariateFunction f = new PolynomialFunction(coeffs);
		RealVector yi = logs.map(new Log10()).map(f); 
		return yi.toArray();
	}

	private double[] interpolateSD(final double[] y, final double[] sd, final double[] yi) {
		final WeightedObservedPoints obs = new WeightedObservedPoints();
		for (int l = 0; l < y.length; l++)
			obs.add(y[l], sd[l]);
		
		// First-degree polynomial fitter (line)
		final PolynomialCurveFitter linfit = PolynomialCurveFitter.create(1);
		final double[] coeffs = linfit.fit(obs.toList());
		final UnivariateFunction f = new PolynomialFunction(coeffs);
		RealVector sdi = new ArrayRealVector(yi).map(f); 
		return sdi.toArray();
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
		if (fitMode == bandMode || (fitMode == continuumMode && lane == ladderLane)) {
			// Guess a set of peaks a set of peaks automatically
			guess = new ArrayRealVector(pg.guess());
			peaks = arrayToPeaks(lane, guess.toArray());

		} else if (fitMode == continuumMode) {
			if (ladder == null || fragmentDistribution == null) {
				log.info("Missing ladder/distribution");
				return null;
			}
			// Use the stored distribution as a guess 
			// fragmentDistribution[:][0] = Frequency
			// fragmentDistribution[:][1] = Length (bp)
			// fragmentDistribution[:][2] = MW
			RealMatrix distMatrix = new Array2DRowRealMatrix(fragmentDistribution);
			final List<Peak> ladderPeaks = getFittedPeaks(ladderLane);
			final double[] meanLadder = new double[ladderPeaks.size()];
			final double[] sdLadder = new double[ladderPeaks.size()];
			for (int p = 0; p < ladderPeaks.size(); p++) {
				meanLadder[p] = ladderPeaks.get(p).getMean();
				sdLadder[p] = ladderPeaks.get(p).getSigma();
			}
			final double[] means = interpolateDisplacement(meanLadder);
			final double[] sds = interpolateSD(meanLadder, sdLadder, means);

			RealVector scaledFrequency = distMatrix.getColumnVector(0);
			scaledFrequency = scaledFrequency.mapDivide(scaledFrequency.getMaxValue());
			RealVector mwArray = distMatrix.getColumnVector(2);
			mwArray = mwArray.map(new Log10());
			RealVector scaledMW = mwArray.mapDivide(mwArray.getMaxValue());
			for (final DataSeries d : inputData) {
				if (d.getLane() == lane) {
					selectedFragments.set(lane-1, new ArrayList<>());
					RealVector profile = d.getY().mapSubtractToSelf(d.getMinY());
					for (int i = 0; i < means.length; i++) {
						// exclude peaks outside the profile domain
						double m = means[i];
						if (m > d.getMinX() && m < d.getMaxX()) {
							selectedFragments.get(lane-1).add(i);
							double frequency = scaledFrequency.getEntry(i);
							double mw = scaledMW.getEntry(i);
							double n = profile.getMaxValue() * mw * frequency;
							double s = sds[i];
							peaks.add(new Peak(lane, n, m, s));
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
					log.info("[LANE " + g.getLane() + "] Peak guess: " + g.getNorm() + ", " + g.getMean() + ", "
					        + g.getSigma());
					g.setMean(c.getMean());
					g.setNorm(c.getNorm());
					g.setSigma(c.getSigma());
					log.info("\t\tReplaced with: " + c.getNorm() + ", " + c.getMean() + ", "
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
		if (fitMode == bandMode) {
			poly.setEntry(1, guess.getEntry(1));
		} else if (fitMode == continuumMode) {
			for (DataSeries d : inputData) {
				if (d.getLane() == lane) 
					poly.setEntry(1, d.getMinY());
			}
		}
		return poly.append(peaksToArray(peaks));
	}

	private double doIntegrate(final RealVector xvals, final double n, final double m,
	        final double s) {
		RealVector weights = new ArrayRealVector(xvals.getDimension()).mapAddToSelf(1.0);
		final GaussIntegrator ti = new GaussIntegrator(xvals.toArray(), weights.toArray());
		final Gaussian gauss = new Gaussian(n, m, s);
		double area = ti.integrate(gauss);
		return area;
	}

	public List<DataSeries> doFit(List<Integer> lanes) {
		int progress = 1;
		ArrayList<DataSeries> in = new ArrayList<>();
		for (int i : lanes) {
			for (DataSeries d : inputData) {
				if (d.getLane() == i) in.add(d);
			}
		}
		ArrayList<DataSeries> out = new ArrayList<>();
		final StopWatch sw = new StopWatch();
		sw.start();
		try {
			in.parallelStream().forEach((d) -> {
				out.addAll(doFit(d.getLane()));
			});
		} catch (NullPointerException np) {
			log.info(in.size());
		}
		String t = String.format("Time elapsed: %1$.1f s\n", sw.getTime()/1000.0);
		log.info(t);
		statusServ.showProgress(++progress, in.size());
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
		        degBG, tolpk);

		final double[] firstGuess = doGuess(lane, pg).toArray();
		final LeastSquaresProblem problem = GaussianArrayCurveFitter.create(fitMode, degBG, polyDerivative, tolpk, areaDrift)
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
		rms.setEntry(lane-1, optimum.getRMS());
		final String outStr = String.format("Lane " + lane + ", RMS: %1$.2f; ", rms.getEntry(lane - 1));
		log.info(outStr);
		return output;
	}

	public String getSummary() {
		String s = "" + "<h1>FIT SUMMARY</h1>";
		s += "<h2>PARAMETERS</h2>";
		s += "<table>";
		if (fitMode == bandMode)
			s += "<tr> <td>Fitting Mode</td> <td>Banded</td>   </tr>";
		else if (fitMode == continuumMode)
			s += "<tr> <td>Fitting Mode</td> <td>Continuum</td></tr>";
		s += "<tr> <td>Peak Height Tolerance</td>        <td>" + tolPK          + "</td></tr>";
		s += "<tr> <td>Polynomial Degree</td>            <td>" + degBG          + "</td></tr>";
		s += "<tr> <td>Maximum Polynomial Derivative</td><td>" + polyDerivative + "</td></tr>";
		s += "<tr> <td>Area Drift Limit</td>             <td>" + areaDrift      + "</td></tr></table>";
		s += "<h2>RESULTS</h2>";
		s += "<table>";
		for (DataSeries d : inputData) {
			int l = d.getLane();
			s += String.format("<tr><td>Lane %1$d;</td>", l);
			s += String.format("    <td>RMS = %1$.2f</td>", rms.getEntry(l-1));
			if (fitMode == bandMode || l == ladderLane)
				s += "<td></td><td></td></tr>";
			else { //(fitMode == continuumMode)
				double mean = statMatrix.getEntry(l-1, 0);
				double standardDeviation = statMatrix.getEntry(l-1, 1);
				s += String.format("<td>Average Fragment Size = %1$.2f</td> <td>(%2$.2f)</td></tr>",
					mean,standardDeviation);
			}
		}
		s += "</table>";
		return s;
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

	public List<Peak> getCustomPeaks(final int lane) {
		final List<Peak> c = new ArrayList<>();
		for (final Peak p : allCustomList) {
			if (p.getLane() == lane)
				c.add(p);
		}
		return c;
	}

	public List<Peak> getGuessPeaks(final int lane) {
		final List<Peak> g = new ArrayList<>();
		Iterator<Peak> it = allGuessList.iterator();
		while (it.hasNext()) {
			Peak p = it.next();
			if (p.getLane() == lane)
				g.add(p);
		}
		return g;
	}

	public List<Peak> getFittedPeaks(final int lane) {
		final List<Peak> f = new ArrayList<>();
		for (final Peak p : allFittedList) {
			if (p.getLane() == lane)
				f.add(p);
		}
		return f;
	}

	public double[] getFittedDistribution(int l) {
		List<Double> list = fittedDistributions.get(l-1);
		double[] array = new double[list.size()];
		for(int i = 0; i < array.length; i++) 
			array[i] = list.get(i);
		return array;
	}

	public void setDegBG(final int degBG) {
		this.degBG = degBG;
	}

	public void setInputData(final ArrayList<DataSeries> inputData) {
		this.inputData = inputData;
		this.selectedFragments = new ArrayList<>();
		this.fittedDistributions = new ArrayList<>();
		this.rms = new ArrayRealVector(inputData.size());
		this.statMatrix = new Array2DRowRealMatrix(inputData.size(), 2);
		for (int i = 0; i< inputData.size(); i++) {
			selectedFragments.add(new ArrayList<>());
			fittedDistributions.add(new ArrayList<>());
		}
	}

	public void setFitMode(final int fitMode) {
		this.fitMode = fitMode;
	}

	public void setLadder(final RealVector ladder) {
		this.ladder = ladder.toArray();
	}

	public void setReferenceLane(final int ladderLane) {
		this.ladderLane = ladderLane;
	}

	public void setFragmentDistribution(final double[][] fragmentDistribution) {
		this.fragmentDistribution = fragmentDistribution;
	}

	public void setPolyDerivative(final double polyDerivative) {
		this.polyDerivative = polyDerivative;
	}
	
	public void setTolPK(final double tolPK) {
		this.tolPK = tolPK;
	}

	public void setAreaDrift(double areaDrift) {
		this.areaDrift = areaDrift;
	}
}

/** Class to generate Gaussian Peak objects
 * Could be expanded to represent other types of peaks
 **/
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
