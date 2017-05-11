
package gausscurvefit;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

import net.imagej.table.Column;
import net.imagej.table.DefaultGenericTable;
import net.imagej.table.DefaultTableDisplay;
import net.imagej.table.GenericColumn;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.function.Abs;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.analysis.integration.TrapezoidIntegrator;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
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

import gausscurvefit.Peak;
import gausscurvefit.GaussianArrayCurveFitter.ParameterGuesser;

public class Fitter {

	@Parameter
	private LogService log;
	@Parameter
	private DisplayService displayServ;
	@Parameter
	private StatusService statusServ;

	private int degBG;
	private double tolPK;
	private final String title; // for the datafile, main image's title
	public static final double peakDistanceTol = 2;
	
	private ArrayList<DataSeries> inputData;
	private ArrayList<PeaksList> allPeaksLists;
	//private ArrayList<ArrayList<Peak>> fittedPeaks;

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

	public ArrayList<ArrayList<DataSeries>> doFit() {
		int progress = 1;
		ArrayList<ArrayList<DataSeries>> outputAll = new ArrayList<>();
		for (DataSeries d : inputData) {
			outputAll.add(doFit(d));
			statusServ.showProgress(++progress, inputData.size());
		}
		return outputAll;
	}
	
	public ArrayList<DataSeries> doFit(final DataSeries in) {	
		final int lane = in.getLane();
		final RealVector xvals = new ArrayRealVector(in.getX());
		final RealVector yvals = new ArrayRealVector(in.getY());
		final ArrayList<DataSeries> funout = new ArrayList<>();
		final ArrayList<Peak> fittedPeaks = new ArrayList<>();
		ArrayList<DataSeries> output = new ArrayList<>();

		// Tolerance as percentage of the range
		final double tolpk = tolPK * (yvals.getMaxValue() - yvals.getMinValue());

		final WeightedObservedPoints obs = new WeightedObservedPoints();
		for (int o = 0; o < xvals.getDimension(); o++) {
			obs.add(xvals.getEntry(o), yvals.getEntry(o));
		}
		
		if (getGuessPeaks(lane) == null) {
			final ParameterGuesser pg = new GaussianArrayCurveFitter.ParameterGuesser(
				obs.toList(), tolpk, degBG);
			RealVector firstGuess = new ArrayRealVector(pg.guess());
		 
			// Initial Guess
			RealVector norms0 = new ArrayRealVector();
			RealVector means0 = new ArrayRealVector();
			RealVector sds0 = new ArrayRealVector();
			ArrayList<Peak> guessPeaks = new ArrayList<>(); 
			for (int b = degBG + 2; b < firstGuess.getDimension(); b += 3) {
				// Initial Guess
				double n = firstGuess.getEntry(b);
				double m = firstGuess.getEntry(b + 1);
				double s = firstGuess.getEntry(b + 2);
				
				norms0 = norms0.append(n);
				means0 = means0.append(m);
				sds0 = sds0.append(s);
				guessPeaks.add(new Peak(lane, n, m, s));
			}
			PeaksList pl = new PeaksList(lane);
			pl.setGuessList(guessPeaks);
		}

		final RealVector firstGuess = peaksToArray(getGuessPeaks(lane));
		final LeastSquaresProblem problem = GaussianArrayCurveFitter.create(tolpk,
			degBG).withStartPoint(firstGuess.toArray()).getProblem(obs.toList());

		final LeastSquaresOptimizer.Optimum optimum =
			new LevenbergMarquardtOptimizer().optimize(problem);

		final RealVector pars = new ArrayRealVector(optimum.getPoint());

		// After fitting
		RealVector norms = new ArrayRealVector();
		RealVector means = new ArrayRealVector();
		RealVector sds = new ArrayRealVector();
		final RealVector poly = pars.getSubVector(0, degBG + 2);

		for (int b = degBG + 2; b < pars.getDimension(); b += 3) {
			norms = norms.append(pars.getEntry(b));
			means = means.append(pars.getEntry(b + 1));
			sds = sds.append(pars.getEntry(b + 2));
			fittedPeaks.add(new Peak(lane,pars.getEntry(b),pars.getEntry(b + 1),pars.getEntry(b + 2)));
		}
		final PolynomialFunction bg = new PolynomialFunction(poly.getSubVector(1,
			degBG + 1).toArray());
		funout.add(new DataSeries("Background", lane, DataSeries.BACKGROUND,
			xvals, bg, Plotter.bgColor));

		for (int gg = 0; gg < norms.getDimension(); gg++) {
			final Gaussian gauss = new Gaussian(norms.getEntry(gg), means.getEntry(
				gg), sds.getEntry(gg));
			final UnivariateFunction[] functs = { bg, gauss };
			funout.add(new DataSeries("Band " + gg+1, lane, DataSeries.GAUSS_BG,
				xvals, functs, Plotter.gaussColor));
		}
		final GaussianArray fitted = new GaussianArray(norms, means, sds,
			poly);
		output.add(new DataSeries("Fit", lane, DataSeries.FITTED, xvals, fitted,
			Plotter.fittedColor));

		final String outStr = String.format("Lane " + lane +
			", RMS: %1$.2f; ", optimum.getRMS());
		log.info(outStr);
		updateResultsTable();
		return output;
	}

	private RealVector peaksToArray(ArrayList<Peak> peakList) {
		RealVector guessArray = new ArrayRealVector();
		guessArray = guessArray.append(degBG).append(new ArrayRealVector(degBG+1));
		for (Peak p : peakList) {
			guessArray = guessArray.append(p.getIntensity());
			guessArray = guessArray.append(p.getDistance());
			guessArray = guessArray.append(p.getFWHM());
		}
		return guessArray;
	}

	private void updateResultsTable() {
		// Results Table Columns
		final String[] headers = { "", "Band", "Distance", "Amplitude", "FWHM",
			"Area", "Dist. G.", "Amp. G.", "FWHM G." };
		final GenericColumn[] tableCol = new GenericColumn[headers.length];
		final DefaultGenericTable rt = new DefaultGenericTable();

		for (int cc = 0; cc < headers.length; cc++)
			tableCol[cc] = new GenericColumn(headers[cc]);
		
		for (PeaksList pl : allPeaksLists) {
			int lane = pl.getLane();
			ArrayList<Peak> guess  = pl.getGuessList();
			ArrayList<Peak> fitted = pl.getFittedList();
			int band = 1;
			for (Peak p : fitted){
				double n = p.getIntensity();
				double m = p.getDistance();
				double s = p.getFWHM();
				double a = 0.0;
				double n0 = guess.get(band).getIntensity();
				double m0 = guess.get(band).getDistance();
				double s0 = guess.get(band).getFWHM();
				
				for (DataSeries d : inputData) {
					if (d.getLane() == lane){
						a = doIntegrate(new ArrayRealVector(d.getX()), n, m, s);
					}
				}
				// Prepare columns for Results Table
				if (band == 1) {
					tableCol[0].add("Lane " + lane);
				} else {
					tableCol[0].add("");
				}
				tableCol[1].add(band);
				tableCol[2].add(String.format("%1$.1f", m));
				tableCol[3].add(String.format("%1$.1f", n));
				tableCol[4].add(String.format("%1$.2f", s));
				tableCol[5].add(String.format("%1$.1f", a));
				tableCol[6].add(String.format("%1$.1f", m0));
				tableCol[7].add(String.format("%1$.1f", n0));
				tableCol[8].add(String.format("%1$.2f", s0));
				//s0 * 2 * FastMath.sqrt(2 * FastMath.log(2))));
			}
			band++;
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
			System.out.println("Exception");
		}
	}
	
	private double doIntegrate(RealVector xvals, double n, double m, double s)
	{
		double area;
		final TrapezoidIntegrator ti = new TrapezoidIntegrator();
		final Gaussian gauss = new Gaussian(n, m, s);
		area = ti.integrate(Integer.MAX_VALUE, gauss, xvals.getMinValue(), xvals.getMaxValue());
		return area;
	}

	public void addCustomPeak(int lane, Peak peak) {
		for (PeaksList pl : allPeaksLists) {
			if (pl.getLane() == lane) {
				pl.addToList(peak);
			}
		}
	}

	public void removeCustomPeak(int lane, Peak peak) {
		for (PeaksList pl : allPeaksLists) {
			if (pl.getLane() == lane) {
				pl.removeFromList(peak);
			}
		}
	}

	public void resetCustomPeaks(final int lane) {
		for (PeaksList pl : allPeaksLists) {
			if (pl.getLane() == lane) 
				pl.setAddList(new ArrayList<Peak>());
		}
	}

	public ArrayList<Peak> getGuessPeaks(final int lane) {
		for (PeaksList pl : allPeaksLists) {
			if (pl.getLane() == lane) return pl.getGuessList();
		}
		return null;
	}

	public ArrayList<Peak> getCustomPeaks(final int lane) {
		for (PeaksList pl : allPeaksLists){
			if (pl.getLane() == lane) {
				return pl.getAddList();
			}
		}
		return null;
	}

	/**
	 * return fittedPeaks
	 */
	public ArrayList<Peak> getFittedPeaks(int lane) {
		for (PeaksList pl : allPeaksLists) {
			if (pl.getLane() == lane) {
				return pl.getFittedList();
			}
		}
		return null;
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
	private ArrayList<Peak> guessList;
	private ArrayList<Peak> fittedList;
	private ArrayList<Peak> addList;

	public PeaksList(final int lane) {
		this.lane = lane;
		this.guessList = new ArrayList<>();
		this.fittedList = new ArrayList<>();
		this.addList = new ArrayList<>();
	}

	public void addToList(Peak peak) {
		boolean found = false;
		for (Peak p : addList) {
			// If close to existing replace
			if (FastMath.abs(peak.getDistance()-p.getDistance()) <= 2) {
				p.setFWHM(peak.getFWHM());
				found = true;break;
			}
		}
		if(!found) addList.add(peak);
		
		// Check if guessList contains peak and update
		found = false;
		for (final Peak gp : guessList) {
			if (FastMath.abs(peak.getDistance()-gp.getDistance()) <= Fitter.peakDistanceTol) {
				gp.setFWHM(peak.getFWHM());
				found = true;
			}
			if(!found) guessList.add(peak);
			Collections.sort(guessList);
			}
		return;
	}

	public void removeFromList(final Peak peak)
	{
		Iterator<Peak> addPeakIter = addList.iterator(); 
		while (addPeakIter.hasNext()) {
			double diff = FastMath.abs(peak.getDistance()-addPeakIter.next().getDistance());
			if (diff <= Fitter.peakDistanceTol) addPeakIter.remove();
		}
		
		Iterator<Peak> guessPeakIter = guessList.iterator(); 
		while (guessPeakIter.hasNext()) {
			double diff = FastMath.abs(peak.getDistance()-guessPeakIter.next().getDistance());
			if (diff <= Fitter.peakDistanceTol) guessPeakIter.remove();
		}
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
	
	public void setGuessList(ArrayList<Peak> guessList) {
		this.guessList = guessList;
	}
	
	public void setAddList(ArrayList<Peak> addList) {
		this.addList = addList;
	}
	
	public void setFittedList(ArrayList<Peak> fittedList) {
		this.fittedList = fittedList;
	}
}

class Peak implements Comparable<Peak> {

		private double distance;
		private double intensity;
		private double fwhm;
		private final int lane;

		public Peak(int lane, final double intensity, final double distance) {
			this.intensity = intensity;
			this.distance = distance;
			this.fwhm = 0.0;
			this.lane = lane;
		}

		public Peak(int lane, final double intensity, final double distance, final double fwhm) {
			this.intensity = intensity;
			this.distance = distance;
			this.fwhm = fwhm;
			this.lane = lane;
		}

		public double getDistance() {
			return distance;

		}

		public double getIntensity() {
			return intensity;
		}

		public double getFWHM() {
			return fwhm;
		}
		
		public double getLane() {
			return lane;
		}
		
		public void setFWHM(double fwhm) {
			this.fwhm = fwhm;
		}

		@Override
		public int compareTo(final Peak p) {
			return Double.compare(distance, p.getDistance());
		}

}

