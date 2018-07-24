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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

import org.apache.commons.lang3.time.StopWatch;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.analysis.function.Pow;
import org.apache.commons.math3.analysis.integration.gauss.GaussIntegrator;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.fitting.PolynomialCurveFitter;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
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
	private double tolPK, areaDrift, sdDrift, polyDerivative, polyOffset;
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
	private RealVector rms; // normalized to y-range
	private Array2DRowRealMatrix statMatrix;
	private List<DataSeries> fittedDistributions;

	public Fitter(final Context context, final String title) {
		context.inject(this);
		this.title = title;
		this.inputData = new ArrayList<>();

		this.allGuessList = new ArrayList<>();
		this.allFittedList = new ArrayList<>();
		this.allCustomList = new ArrayList<>();
		this.polyOffset = 0.98;
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

	public void updateResultsTable(final String savePath) {
		// Results Table Columns
		final String[] headersB = { "Lane", "Band", "Distance", "Dist. G.",
			"Amplitude", "Amp. G.", "FWHM", "FWHM G.", "Area" };
		final String[] headersF = { "Lane", "Band", "Distance", "Dist. G.",
			"Amplitude", "Amp. G.", "FWHM", "FWHM G.", "Area", "Frequency", "BP",
			"MW" };
		String[] headers = null;
		if (fitMode == bandMode) headers = headersB;
		else headers = headersF;
		final ResultsTable rt = new ResultsTable();

		for (final DataSeries d : inputData) {
			final int lane = d.getLane();
			final List<Peak> guess = getGuessPeaks(lane);
			final List<Peak> fitted = getFittedPeaks(lane);
			RealVector areas = new ArrayRealVector();
			if (guess.size() != fitted.size() || 
					guess.size() == 0 || fitted.size() == 0) {
				log.info("Data Size " + lane + " : " + guess.size() + ", " + fitted.size());
			}

			int band = 1;
			RealVector means = new ArrayRealVector();
			int listSize = (fitMode == continuumMode && lane != ladderLane) ? 
				selectedFragments.get(lane - 1).size() : guess.size();
			for (int p = 0; p < listSize; p++) {
				final double n = fitted.get(p).getNorm();
				final double m = fitted.get(p).getMean();
				final double s = fitted.get(p).getSigma();
				final double a = doIntegrate(d.getX(), n, m, s);
				means = means.append(m);

				final double n0 = guess.get(p).getNorm();
				final double m0 = guess.get(p).getMean();
				final double s0 = guess.get(p).getSigma();
				areas = areas.append(a);

				// Add values to Results Table
				if (m>d.getMinX() && m< d.getMaxX()) {
					rt.incrementCounter();
					if (band == 1) {
						rt.addValue(headers[0], "" + lane);
					}
					else {
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
						}
						else {
							final double freq = fragmentDistribution[selectedFragments.get(
								lane - 1).get(p)][0];
							final double bp = fragmentDistribution[selectedFragments.get(lane -
								1).get(p)][1];
							final double mw = fragmentDistribution[selectedFragments.get(lane -
								1).get(p)][2];
							rt.addValue(headers[9], String.format("%1$.3f", freq));
							rt.addValue(headers[10], String.format("%1$.0f", bp));
							rt.addValue(headers[11], String.format("%1$.3f", mw));
						}
					}
					band++;
				}
			}

			if (fitMode == continuumMode && lane != ladderLane) {
				RealVector bpSubarray = new ArrayRealVector();
				RealVector wSubarray  = new ArrayRealVector();
				final RealMatrix distMatrix = new Array2DRowRealMatrix(
					fragmentDistribution);
				final RealVector bpArray = distMatrix.getColumnVector(1);
				final RealVector wArray  = distMatrix.getColumnVector(2);

				for (final int i : selectedFragments.get(lane - 1)) {
					bpSubarray = bpSubarray.append(bpArray.getEntry(i));
					wSubarray  =  wSubarray.append( wArray.getEntry(i));
				}
				final double meanLength = new Mean().evaluate(
					bpSubarray.toArray(), wSubarray.toArray());
				final RealVector scaleFactor = areas.ebeMultiply(bpSubarray);
				final double weighedMeanLength = new Mean().evaluate(bpSubarray
					.toArray(), scaleFactor.toArray());
				final double sdLength = FastMath.sqrt(new Variance().evaluate(bpSubarray
					.toArray(), scaleFactor.toArray()));
				statMatrix.setRow(lane - 1, new double[] { weighedMeanLength,
					sdLength, meanLength });
				final String info = String.format(
					"Lane %1$d, Average length: %2$.2f(%3$.2f) %4$.2f", lane,
					weighedMeanLength, sdLength, meanLength);

				log.info(info);

				// Update the fitted distribution list with appropriate scale
				RealVector x = new ArrayRealVector();
				RealVector y = new ArrayRealVector();
				final Iterator<DataSeries> it = fittedDistributions.iterator();
				while (it.hasNext()) {
					final DataSeries dfit = it.next();
					if (dfit.getLane() == lane) {
						x = dfit.getX();
						y = dfit.getY();
						it.remove();
					}
				}

				final WeightedObservedPoints obs = new WeightedObservedPoints();
				for (int l = 0; l < bpSubarray.getDimension(); l++)
					obs.add(means.getEntry(l), Math.log10(bpSubarray.getEntry(l)));

				// First-degree polynomial fitter (line)
				final PolynomialCurveFitter linfit = PolynomialCurveFitter.create(1);
				final double[] coeffs = linfit.fit(obs.toList());
				final UnivariateFunction f = new PolynomialFunction(coeffs);
				x = x.map(f);
				for (int i = 0; i < x.getDimension(); i++)
					x.setEntry(i, new Pow().value(10, x.getEntry(i)));
				y = y.ebeMultiply(x);
				y = y.mapDivide(y.getMaxValue());
				fittedDistributions.add(new DataSeries("Fit", lane, DataSeries.FITTED,
					x, y, Plotter.fittedColor));
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
				if (i != ladderLane - 1) {
					line = String.format("Lane %1$d\t%2$.2f \t %3$.2f\n", (i + 1),
						statMatrix.getEntry(i, 0), statMatrix.getEntry(i, 1));
					outText += line;
				}
			}
		}

		final String file = "Fit of " + title + ".xls";
		new File(savePath).mkdirs();
		log.info("Saving to " + savePath + file + " ...");
		try (BufferedWriter out = new BufferedWriter(new FileWriter(savePath +
			file)))
		{
			out.write(outText);
			out.close();
		}
		catch (final IOException e) {
			log.info("Exception", e);
		}
	}

	

	private SortedParameters doGuess(final int lane, final ParameterGuesser pg) {
		// Guess a set of peaks a set of peaks
		SortedParameters guess = pg.guess();
		List<Peak> peaks = arrayToPeaks(lane, guess.getParameters());

		// Check if custom peaks are needed
		for (final Peak c : getCustomPeaks(lane)) {
			// Check if guessList contains peak and update
			double minDistance = Double.POSITIVE_INFINITY;
			Peak closest = c;
			for (final Peak g : peaks) {
				final double distance = FastMath.abs(c.getMean() - g.getMean());
				if (distance < minDistance) {
					minDistance = distance;
					closest = g;
				}
			}
			if ((fitMode == bandMode && minDistance <= peakDistanceTol) ||
				fitMode == continuumMode)
			{ // replace closest
				log.info("[LANE " + closest.getLane() + "] Peak guess: " + closest
					.getNorm() + ", " + closest.getMean() + ", " + closest.getSigma());
				closest.setMean(c.getMean());
				closest.setNorm(c.getNorm());
				closest.setSigma(c.getSigma());
				log.info("\t\tReplaced with: " + c.getNorm() + ", " + c.getMean() +
					", " + c.getSigma());
				continue;
			}
			peaks.add(c);
			Collections.sort(peaks);
		}
		allGuessList.addAll(peaks);
		RealVector poly = new ArrayRealVector();
		poly = poly.append(degBG).append(new ArrayRealVector(degBG + 1));

		for (final DataSeries d : inputData) {
			if (d.getLane() == lane) poly.setEntry(1, d.getMinY());
		}
		
		return new SortedParameters(poly.append(peaksToArray(peaks)).toArray());
	}

	private double doIntegrate(final RealVector xvals, final double n,
		final double m, final double s)
	{
		final RealVector weights = new ArrayRealVector(xvals.getDimension())
			.mapAddToSelf(1.0);
		final GaussIntegrator ti = new GaussIntegrator(xvals.toArray(), weights
			.toArray());
		final Gaussian gauss = new Gaussian(n, m, s);
		final double area = ti.integrate(gauss);
		return area;
	}

	public List<DataSeries> doFit(final List<Integer> lanes) {
		
		final ArrayList<DataSeries> in = new ArrayList<>();
		for (final int i : lanes) {
			for (final DataSeries d : inputData) {
				if (d.getLane() == i) in.add(d);
			}
		}
		final ArrayList<DataSeries> out = new ArrayList<>();
		final StopWatch sw = new StopWatch();
		sw.start();
		try {
			AtomicInteger progress = new AtomicInteger();
			in.parallelStream().forEach((d) -> {
				out.addAll(doFit(d.getLane()));
				statusServ.showProgress(progress.incrementAndGet(), in.size());
			});
		}
		catch (final NullPointerException np) {
			log.info(in.size());
		}
		final String t = String.format("Time elapsed: %1$.1f s\n", sw.getTime() /
			1000.0);
		log.info(t);
		
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
		ParameterGuesser pg = null;
		if (fitMode == bandMode) {
			pg = new GaussianArrayCurveFitter.ParameterGuesser(
				obs.toList(), degBG, tolpk, polyOffset);
		} else {
			final RealMatrix distMatrix = new Array2DRowRealMatrix(
				fragmentDistribution);
			pg = new GaussianArrayCurveFitter.ParameterGuesser(
				obs.toList(), degBG, tolpk, polyOffset, distMatrix, 
				getFittedPeaks(ladderLane), ladder);
			selectedFragments.set(lane - 1, 
				pg.getUsedFragments(new double[] {xvals.getMinValue(), xvals.getMaxValue()}));
		}
		final SortedParameters firstGuess = doGuess(lane, pg);
		final GaussianArrayCurveFitter cf = GaussianArrayCurveFitter.create(fitMode,
			degBG, polyDerivative, polyOffset, tolpk, areaDrift, sdDrift).withStartPoint(firstGuess);
		final LeastSquaresProblem problem = cf.getProblem(obs.toList());

		final LeastSquaresOptimizer.Optimum optimum = cf.getOptimizer().optimize(problem);

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
		final PolynomialFunction bg = new PolynomialFunction(poly.getSubVector(1,
			degBG + 1).toArray());
		output.add(new DataSeries("Background", lane, DataSeries.BACKGROUND, xvals,
			bg, Plotter.bgColor));

		for (final Peak p : fittedPeaks) {
			final Gaussian gauss = new Gaussian(p.getNorm(), p.getMean(), p
				.getSigma());
			final UnivariateFunction[] functs = { bg, gauss };
			output.add(new DataSeries("Band " + 1, lane, DataSeries.GAUSS_BG, xvals,
				functs, Plotter.gaussColor));
		}
		double[] par = poly.append(peaksToArray(fittedPeaks)).toArray();
		final GaussianArray fittedCurve = new GaussianArray(new SortedParameters(par));
		final DataSeries fit = new DataSeries("Fit", lane, DataSeries.FITTED, xvals,
			fittedCurve, Plotter.fittedColor);
		boolean contains = false;
		for (DataSeries d : fittedDistributions) {
			if (d.getLane() == lane) {
				contains = true;
				d = fit;
			}
		}
		if (!contains) fittedDistributions.add(fit);
		output.add(fit);

		// Print RMS to console
		rms.setEntry(lane - 1, optimum.getRMS()/
			(inputData.get(lane - 1).getMaxY() - inputData.get(lane - 1).getMinY()));
		final String outStr = String.format("Lane " + lane + ", RMS: %1$.4f; ", rms
			.getEntry(lane - 1));
		log.info(outStr);
		return output;
	}

	public String getSummary() {
		String s = "" + "<h1>FIT SUMMARY</h1>";
		s += "<h2>PARAMETERS</h2>";
		s += "<table>";
		if (fitMode == bandMode) s +=
			"<tr> <td>Fitting Mode</td> <td>Banded</td>   </tr>";
		else if (fitMode == continuumMode) s +=
			"<tr> <td>Fitting Mode</td> <td>Continuum</td></tr>";
		s += "<tr> <td>Peak Height Tolerance</td>        <td>" + tolPK +
			"</td></tr>";
		s += "<tr> <td>Polynomial Degree</td>            <td>" + degBG +
			"</td></tr>";
		s += "<tr> <td>Maximum Polynomial Derivative</td><td>" + polyDerivative +
			"</td></tr>";
		s += "<tr> <td>Area Drift Limit</td>             <td>" + areaDrift +
			"</td></tr>";
		s += "<tr> <td>SD Drift Limit</td>             <td>" + sdDrift +
				"</td></tr></table>";
		s += "<h2>RESULTS</h2>";
		s += "<table>";
		for (final DataSeries d : inputData) {
			final int l = d.getLane();
			s += String.format("<tr><td>Lane %1$d;</td>", l);
			s += String.format("    <td>RMS = %1$.4f</td>", rms.getEntry(l - 1));
			if (fitMode == bandMode || l == ladderLane) s +=
				"<td></td><td></td></tr>";
			else { // (fitMode == continuumMode)
				final double mean = statMatrix.getEntry(l - 1, 0);
				final double standardDeviation = statMatrix.getEntry(l - 1, 1);
				final double aciMean = statMatrix.getEntry(l - 1, 2);
				s += String.format(
					"<td>Average Fragment Size = %1$.0f</td> <td>± %2$.0f (%3$.0f)</td></tr>",
					mean, standardDeviation, aciMean);
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
				if (FastMath.abs(peak.getMean() - p
					.getMean()) <= Fitter.peakDistanceTol)
				{
					p.setSigma(peak.getSigma());
					found = true;
					break;
				}
			}
		}
		if (!found) allCustomList.add(peak);
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
			if (peakIter.next().getLane() == lane) peakIter.remove();
		}
	}

	public void resetFit(final int lane) {
		final Iterator<Peak> itFitted = allFittedList.iterator();
		while (itFitted.hasNext()) {
			final Peak f = itFitted.next();
			if (f == null || f.getLane() == lane) itFitted.remove();
		}
		Collections.sort(allFittedList);

		final Iterator<Peak> itGuess = allGuessList.iterator();
		while (itGuess.hasNext()) {
			final Peak g = itGuess.next();
			if (g == null || g.getLane() == lane) itGuess.remove();
		}
		Collections.sort(allGuessList);
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
			if (p.getLane() == lane) c.add(p);
		}
		return c;
	}

	public List<Peak> getGuessPeaks(final int lane) {
		final List<Peak> g = new ArrayList<>();
		final Iterator<Peak> it = allGuessList.iterator();
		while (it.hasNext()) {
			final Peak p = it.next();
			if (p == null) continue;	
			if (p.getLane() == lane) g.add(p);
		}
		return g;
	}

	public List<Peak> getFittedPeaks(final int lane) {
		final List<Peak> f = new ArrayList<>();
		for (final Peak p : allFittedList) {
			if (p.getLane() == lane) f.add(p);
		}
		return f;
	}

	public DataSeries getFittedDistribution(final int l) {
		return fittedDistributions.get(l - 1);
	}

	public void setDegBG(final int degBG) {
		this.degBG = degBG;
	}

	public void setInputData(final ArrayList<DataSeries> inputData) {
		for (DataSeries d : this.inputData) 
			resetFit(d.getLane());
		
		this.inputData = inputData;

		for (DataSeries d :inputData) {
			double minX = d.getMinX();
			double maxX = d.getMaxX();
			Iterator<Peak> pIter = allCustomList.iterator();
			while (pIter.hasNext()) {
				Peak c = pIter.next();
				if ( c.getLane() == d.getLane() && (c.getMean() < minX || c.getMean() > maxX))
					pIter.remove();
			}
		}
		this.selectedFragments = new ArrayList<>();
		this.fittedDistributions = new ArrayList<>();
		this.rms = new ArrayRealVector(inputData.size());
		this.statMatrix = new Array2DRowRealMatrix(inputData.size(), 3);
		for (int i = 0; i < inputData.size(); i++) {
			selectedFragments.add(new ArrayList<>());
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

	public void setAreaDrift(final double areaDrift) {
		this.areaDrift = areaDrift;
	}
	public void setSDDrift(final double sdDrift) {
		this.sdDrift = sdDrift;
	}
}

/**
 * Class to generate Gaussian Peak objects Could be expanded to represent other
 * types of peaks
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

	public Peak(final int lane, final double norm, final double mean,
		final double sd)
	{
		this.lane = lane;
		this.norm = norm;
		this.mean = mean;

		if (sd > 0) {
			this.sd = sd;
		}
		else {
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
		}
		else {
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
