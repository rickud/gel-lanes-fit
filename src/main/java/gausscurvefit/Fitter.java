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

package gausscurvefit;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

import net.imagej.table.DefaultGenericTable;
import net.imagej.table.DefaultTableDisplay;
import net.imagej.table.GenericColumn;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.analysis.integration.TrapezoidIntegrator;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;
import org.scijava.Context;
import org.scijava.app.StatusService;
import org.scijava.display.DisplayService;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;

import gausscurvefit.GaussianArrayCurveFitter.ParameterGuesser;

class Fitter {

	@Parameter
	private LogService log;
	@Parameter
	private DisplayService displayServ;
	@Parameter
	private StatusService statusServ;

	public static final double peakDistanceTol = 2;
	private static final double sd2FWHM = 2 * FastMath.sqrt(2 *
		FastMath.log(2));
	private int degBG;
	private double tolPK;
	private final String title; // for the datafile, main image's title
	private ArrayList<DataSeries> inputData;
	private ArrayList<PeaksList> allPeaksLists;

	public Fitter(final Context context, final String title, final int degBG,
		final double tolPK)
	{
		context.inject(this);
		this.title = title;
		this.degBG = degBG;
		this.tolPK = tolPK;
		this.inputData = new ArrayList<>();
		this.allPeaksLists = new ArrayList<>();
	}

	private double[] peaksToArray(final ArrayList<Peak> peakList) {
		RealVector guessArray = new ArrayRealVector();
		guessArray = guessArray.append(degBG).append(new ArrayRealVector(degBG +
			1));
		for (final Peak p : peakList) {
			guessArray = guessArray.append(p.getNorm());
			guessArray = guessArray.append(p.getMean());
			guessArray = guessArray.append(p.getSigma());
		}
		return guessArray.toArray();
	}

	private void updateResultsTable() {
		// Results Table Columns
		final String[] headers = { "", "Band", "Distance", "Amplitude", "FWHM",
			"Area", "Dist. G.", "Amp. G.", "FWHM G." };
		final GenericColumn[] tableCol = new GenericColumn[headers.length];
		final DefaultGenericTable rt = new DefaultGenericTable();

		for (int cc = 0; cc < headers.length; cc++)
			tableCol[cc] = new GenericColumn(headers[cc]);

		for (final PeaksList pl : allPeaksLists) {
			final int lane = pl.getLane();
			final ArrayList<Peak> guess = pl.getGuessList();
			final ArrayList<Peak> fitted = pl.getFittedList();
			int band = 1;
			for (final Peak p : fitted) {
				final double n = p.getNorm();
				final double m = p.getMean();
				final double s = p.getSigma();
				double a = 0.0;
				final double n0 = guess.get(band - 1).getNorm();
				final double m0 = guess.get(band - 1).getMean();
				final double s0 = guess.get(band - 1).getSigma();

				for (final DataSeries d : inputData) {
					if (d.getLane() == lane) {
						a = doIntegrate(new ArrayRealVector(d.getX()), n, m, s);
					}
				}
				// Prepare columns for Results Table
				if (band == 1) {
					tableCol[0].add("Lane " + lane);
				}
				else {
					tableCol[0].add("");
				}
				tableCol[1].add(band);
				tableCol[2].add(String.format("%1$.1f", m));
				tableCol[3].add(String.format("%1$.1f", n));
				tableCol[4].add(String.format("%1$.2f", s * sd2FWHM));
				tableCol[5].add(String.format("%1$.1f", a));
				tableCol[6].add(String.format("%1$.1f", m0));
				tableCol[7].add(String.format("%1$.1f", n0));
				tableCol[8].add(String.format("%1$.2f", s0 * sd2FWHM));
				band++;
			}
		}

		for (int cc = 0; cc < headers.length; cc++)
			rt.add(tableCol[cc]);

		final DefaultTableDisplay tableDisplay = (DefaultTableDisplay) displayServ
			.createDisplay("Results Display", rt);
		displayServ.setActiveDisplay(tableDisplay);

		final String saveFile = "Fit of " + title + ".xls";
		try (BufferedWriter out = new BufferedWriter(new FileWriter(saveFile))) {
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
		}
		catch (final IOException e) {
			log.info("Exception", e);
		}
	}

	private double doIntegrate(final RealVector xvals, final double n,
		final double m, final double s)
	{
		double area;
		final TrapezoidIntegrator ti = new TrapezoidIntegrator();
		final Gaussian gauss = new Gaussian(n, m, s);
		area = ti.integrate(Integer.MAX_VALUE, gauss, xvals.getMinValue(), xvals
			.getMaxValue());
		return area;
	}

	public ArrayList<DataSeries> doFit() {
		int progress = 1;
		final ArrayList<DataSeries> outputAll = new ArrayList<>();
		for (final DataSeries d : inputData) {
			outputAll.addAll(doFit(d));
			statusServ.showProgress(++progress, inputData.size());
		}
		updateResultsTable();
		return outputAll;
	}

	private ArrayList<DataSeries> doFit(final DataSeries in) {
		final int lane = in.getLane();
		final RealVector xvals = in.getX();
		final RealVector yvals = in.getY();

		final ArrayList<Peak> fittedPeaks = new ArrayList<>();
		ArrayList<DataSeries> output = new ArrayList<>();
		
		// Tolerance as percentage of the range
		final double tolpk = tolPK * (yvals.getMaxValue() - yvals.getMinValue());

		final WeightedObservedPoints obs = new WeightedObservedPoints();
		for (int o = 0; o < xvals.getDimension(); o++) {
			obs.add(xvals.getEntry(o), yvals.getEntry(o));
		}
		
		double[] firstGuess = peaksToArray(getGuessPeaks(lane));
		boolean foundList = false;
		for (PeaksList pl : allPeaksLists) {
			if (pl.getLane() == lane) foundList = true;
		}
		
		if (!foundList) {
			final ParameterGuesser pg = new GaussianArrayCurveFitter.ParameterGuesser(
				obs.toList(), tolpk, degBG);
			PeaksList pl = new PeaksList(lane, pg);
			allPeaksLists.add(pl);
			firstGuess = peaksToArray(pl.getGuessList());
		} else {
			firstGuess = peaksToArray(getGuessPeaks(lane));
		}
		
		final LeastSquaresProblem problem = GaussianArrayCurveFitter.create(tolpk,
			degBG).withStartPoint(firstGuess).getProblem(obs.toList());

		final LeastSquaresOptimizer.Optimum optimum =
			new LevenbergMarquardtOptimizer().optimize(problem);

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
		for (final PeaksList pl : allPeaksLists) {
			if (pl.getLane() == lane) pl.setFittedList(fittedPeaks);
		}
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
		final GaussianArray fittedCurve = new GaussianArray(norms, means, sds, poly);
		output.add(new DataSeries("Fit", lane, DataSeries.FITTED, xvals, fittedCurve,
			Plotter.fittedColor));

		// Print RMS to console
		final String outStr = String.format("Lane " + lane + ", RMS: %1$.2f; ",
			optimum.getRMS());
		log.info(outStr);
		return output;
	}

	public void addCustomPeak(final int lane, final Peak peak) {
		for (final PeaksList pl : allPeaksLists) {
			if (pl.getLane() == lane) {
				pl.addToList(peak);
			}
		}
	}

	public void removeCustomPeak(final int lane, final Peak peak)
	{
		for (final PeaksList pl : allPeaksLists) {
			if (pl.getLane() == lane) {
				pl.removeFromList(peak);
			}
		}
	}

	public void resetCustomPeaks(final int lane) {
		for (final PeaksList pl : allPeaksLists) {
			if (pl.getLane() == lane) pl.setAddList(new ArrayList<Peak>());
			pl.setGuessList(new ArrayList<Peak>());
		}
	}

	public void resetFit() {
		for (final PeaksList pl : allPeaksLists) {
			pl.setFittedList(new ArrayList<>());
		}
	}

	public void resetAllFitter() {
		inputData = new ArrayList<>();
		allPeaksLists = new ArrayList<>();
	}

	/**
	 * return guessList
	 */
	public ArrayList<Peak> getGuessPeaks(final int lane) {
		for (final PeaksList pl : allPeaksLists) {
			if (pl.getLane() == lane) return pl.getGuessList();
		}
		return new ArrayList<>();
	}

	/**
	 * return addList
	 */
	public ArrayList<Peak> getCustomPeaks(final int lane) {
		for (final PeaksList pl : allPeaksLists) {
			if (pl.getLane() == lane) {
				return pl.getAddList();
			}
		}
		return new ArrayList<>();
	}

	/**
	 * return fittedList
	 */
	public ArrayList<Peak> getFittedPeaks(final int lane) {
		for (final PeaksList pl : allPeaksLists) {
			if (pl.getLane() == lane) {
				return pl.getFittedList();
			}
		}
		return new ArrayList<>();
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

	/**
	 * @param tolPK
	 */
	public void setTolPK(final double tolPK) {
		this.tolPK = tolPK;

	}
}

class PeaksList {

	private final int lane;
	private ParameterGuesser pg;
	private ArrayList<Peak> guessList;
	private ArrayList<Peak> fittedList;
	private ArrayList<Peak> addList;

	public PeaksList(final int lane, final ParameterGuesser pg) {
		this.lane = lane;
		this.pg = pg;
		this.guessList = new ArrayList<>();
		this.fittedList = new ArrayList<>();
		this.addList = new ArrayList<>();
		reGuess();
	}

	private void reGuess() {
		guessList = new ArrayList<>();
		final RealVector guess = new ArrayRealVector(pg.guess());
		int gaussStart = 0;
		if (guess.getEntry(0) != -1) 
			gaussStart = (int) guess.getEntry(0)+2;
		
		RealVector norms0 = new ArrayRealVector();
		RealVector means0 = new ArrayRealVector();
		RealVector sds0 = new ArrayRealVector();
		for (int b = gaussStart; b < guess.getDimension(); b += 3) {
			// Initial Guess
			final double n = guess.getEntry(b);
			final double m = guess.getEntry(b + 1);
			final double s = guess.getEntry(b + 2);
	
			norms0 = norms0.append(n);
			means0 = means0.append(m);
			sds0 = sds0.append(s);
			guessList.add(new Peak(lane, n, m, s));
		}
		
		for (Peak al : addList) {
			// Check if guessList contains peak and update
			boolean found = false;
			for (Peak gl : guessList) {
				if (FastMath.abs(al.getMean() - gl
					.getMean()) <= Fitter.peakDistanceTol)
				{
					gl.setSigma(al.getSigma());
					found = true;
					break;
				}
			}
			if (!found) {
				guessList.add(al);
				Collections.sort(guessList);
			}
		}
	}

	public void addToList(final Peak peak) {
		boolean found = false;
		for (final Peak p : addList) {
			// If close to existing replace
			if (FastMath.abs(peak.getMean() - p
				.getMean()) <= Fitter.peakDistanceTol)
			{
				p.setSigma(peak.getSigma());
				found = true;
				break;
			}
		}
		if (!found) addList.add(peak);
		reGuess();
	}

	public void removeFromList(final Peak peak) {
		// Remove from custom list
		final Iterator<Peak> addPeakIter = addList.iterator();
		while (addPeakIter.hasNext()) {
			final double diff = FastMath.abs(peak.getMean() - addPeakIter.next()
				.getMean());
			if (diff <= Fitter.peakDistanceTol) addPeakIter.remove();
		}
		reGuess();
	}

	public int getLane() {
		return lane;
	}

	public ArrayList<Peak> getGuessList() {
		return guessList;
	}

	public ArrayList<Peak> getAddList() {
		return addList;
	}

	public ArrayList<Peak> getFittedList() {
		return fittedList;
	}

	public void setGuessList(final ArrayList<Peak> guessList) {
		this.guessList = guessList;
		reGuess();
	}

	public void setAddList(final ArrayList<Peak> addList) {
		this.addList = addList;
		reGuess();
	}

	public void setFittedList(final ArrayList<Peak> fittedList) {
		this.fittedList = fittedList;
	}
}

class Peak implements Comparable<Peak> {

	private final int lane;
	private final double mean;
	private final double norm;
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

	public double getLane() {
		return lane;
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
		return Double.compare(mean, p.getMean());
	}
}
