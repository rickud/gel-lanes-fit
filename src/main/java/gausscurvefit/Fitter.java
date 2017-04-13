package gausscurvefit;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;

import net.imagej.table.DefaultGenericTable;
import net.imagej.table.DefaultTableDisplay;
import net.imagej.table.GenericColumn;

import org.apache.commons.math3.analysis.UnivariateFunction;
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
	private String title; // for the datafile, typically the main imag's title

	private ArrayList<DataSeries> inputData; 
	private ArrayList<CustomPoints> customPoints;
	
	public Fitter(Context context, String title, int degBG, double tolPK) {
		context.inject(this);
		this.title = title;
		this.degBG = degBG;
		this.tolPK = tolPK;
		this.inputData= new ArrayList<>();
		this.customPoints = new ArrayList<>();
	}

	public ArrayList<ArrayList<DataSeries>> doFit() {
		// Results Table Columns
		final ArrayList<String> colLane = new ArrayList<>();
		final ArrayList<Integer> colBand = new ArrayList<>();
		RealVector colDistance = new ArrayRealVector();
		RealVector colAmplitude = new ArrayRealVector();
		RealVector colFWHM = new ArrayRealVector();
		RealVector colArea = new ArrayRealVector();
		RealVector colDistance_g = new ArrayRealVector();
		RealVector colAmplitude_g = new ArrayRealVector();
		RealVector colFWHM_g = new ArrayRealVector();
	
		final ArrayList<ArrayList<DataSeries>> outputData = new ArrayList<>();
		final Iterator<DataSeries> dataIter = inputData.iterator();
		int progress = 1;
		int datasetNumber = 0;
		while (dataIter.hasNext()) {
			final DataSeries profile = dataIter.next();
			final RealVector xvals = new ArrayRealVector(profile.getX());
			final RealVector yvals = new ArrayRealVector(profile.getY());
			final ArrayList<DataSeries> funout = new ArrayList<>();

			// Tolerance as percentage of the range
			final double tolpk = tolPK * (yvals.getMaxValue() - yvals
				.getMinValue());
	
			final WeightedObservedPoints obs = new WeightedObservedPoints();
			for (int o = 0; o < xvals.getDimension(); o++) {
				obs.add(xvals.getEntry(o), yvals.getEntry(o));
			}
	
			final ParameterGuesser pg = new GaussianArrayCurveFitter.ParameterGuesser(
				obs.toList(), tolpk, degBG);
			final RealVector firstGuess = new ArrayRealVector(pg.guess());
			final LeastSquaresProblem problem = GaussianArrayCurveFitter.create(tolpk,
				degBG).getProblem(obs.toList());
	
			final LeastSquaresOptimizer.Optimum optimum =
				new LevenbergMarquardtOptimizer().optimize(problem);
	
			final RealVector pars = new ArrayRealVector(optimum.getPoint());
	
			// Initial Guess
			RealVector norms0 = new ArrayRealVector();
			RealVector means0 = new ArrayRealVector();
			RealVector sds0 = new ArrayRealVector();
	
			// After fitting
			RealVector norms = new ArrayRealVector();
			RealVector means = new ArrayRealVector();
			RealVector sds = new ArrayRealVector();
			final RealVector poly = pars.getSubVector(0, degBG + 2);
	
			final PolynomialFunction bg = new PolynomialFunction(poly.getSubVector(1,
				degBG + 1).toArray());
			funout.add(new DataSeries("Background", DataSeries.BACKGROUND, xvals, bg,
				Plotter.bgColor));
			for (int b = degBG + 2; b < pars.getDimension(); b += 3) {
				// Initial Guess
				norms0 = norms0.append(firstGuess.getEntry(b));
				means0 = means0.append(firstGuess.getEntry(b + 1));
				sds0 = sds0.append(firstGuess.getEntry(b + 2));
	
				// After fitting
				norms = norms.append(pars.getEntry(b));
				means = means.append(pars.getEntry(b + 1));
				sds = sds.append(pars.getEntry(b + 2));
			}
	
			for (int gg = 1; gg < norms.getDimension(); gg++) {
				final Gaussian gauss = new Gaussian(norms.getEntry(gg), means.getEntry(
					gg), sds.getEntry(gg));
				final UnivariateFunction[] functs = { bg, gauss };
				funout.add(new DataSeries("Band " + gg, DataSeries.GAUSS_BG, xvals,
					functs, Plotter.gaussColor));
			}
			final GaussianArrayBG fitted = new GaussianArrayBG(norms, means, sds,
				poly);
			funout.add(new DataSeries("Fit", DataSeries.FITTED, xvals, fitted,
				Plotter.fittedColor));
			outputData.add(funout);
	
			final RealVector peakAreas = doIntegrate(xvals, norms, means, sds);
	
			// Prepare columns for Results Table
			colLane.add("Lane " + datasetNumber++);
			colBand.add(1);
			for (int rr = 1; rr < peakAreas.getDimension(); rr++) {
				colLane.add("");
				colBand.add(rr + 1);
			}
			colDistance = colDistance.append(means);
			colAmplitude = colAmplitude.append(norms);
			colFWHM = colFWHM.append(sds.mapMultiplyToSelf(2 * FastMath.sqrt(2 *
				FastMath.log(2))));
			colArea = colArea.append(peakAreas);
			colDistance_g = colDistance_g.append(means0);
			colAmplitude_g = colAmplitude_g.append(norms0);
			colFWHM_g = colFWHM_g.append(sds0.mapMultiplyToSelf(2 * FastMath.sqrt(2 *
				FastMath.log(2))));
	
			final String output = String.format("Lane " + datasetNumber +", RMS: %1$.2f; ",
				optimum.getRMS());
			log.info(output);
			statusServ.showProgress(++progress, inputData.size());
		}
	
		final String[] headers = { "", "Band", "Distance", "Amplitude", "FWHM",
			"Area", "Dist. G.", "Amp. G.", "FWHM G." };
		final GenericColumn[] tableCol = new GenericColumn[headers.length];
		final DefaultGenericTable rt = new DefaultGenericTable();
	
		for (int cc = 0; cc < headers.length; cc++)
			tableCol[cc] = new GenericColumn(headers[cc]);
		tableCol[0].addAll(colLane);
		tableCol[1].addAll(colBand);
		for (int rr = 0; rr < colLane.size(); rr++) {
			tableCol[2].add(String.format("%1$.1f", colDistance.getEntry(rr)));
			tableCol[3].add(String.format("%1$.1f", colAmplitude.getEntry(rr)));
			tableCol[4].add(String.format("%1$.2f", colFWHM.getEntry(rr)));
			tableCol[5].add(String.format("%1$.1f", colArea.getEntry(rr)));
			tableCol[6].add(String.format("%1$.1f", colDistance_g.getEntry(rr)));
			tableCol[7].add(String.format("%1$.1f", colAmplitude_g.getEntry(rr)));
			tableCol[8].add(String.format("%1$.2f", colFWHM_g.getEntry(rr)));
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
		return outputData;
	}

	private RealVector doIntegrate(final RealVector xvals, final RealVector norms,
		final RealVector means, final RealVector sds)
	{
		RealVector areas = new ArrayRealVector();
		for (int i = 0; i < norms.getDimension(); i++) {
			final TrapezoidIntegrator ti = new TrapezoidIntegrator();
			areas = areas.append(ti.integrate(Integer.MAX_VALUE, new Gaussian(norms
				.getEntry(i), means.getEntry(i), sds.getEntry(i)), xvals.getMinValue(),
				xvals.getMaxValue()));
		}
		return areas;
	}
	
	public ArrayList<CustomPoints> getCustomPoints(){
		return customPoints;
	}
	
	public void setCustomPoints(ArrayList<CustomPoints> customPoints){
		this.customPoints = customPoints;
	}
	
	public void setInputData(ArrayList<DataSeries> inputData) {
		this.inputData = inputData;
	}
	
	public void resetCustomPoints(int laneNumber) {
		final Iterator<CustomPoints> pointsListsIter = customPoints.iterator();
		while (pointsListsIter.hasNext()) {
			final CustomPoints pointsList = pointsListsIter.next();
			if (pointsList.getLane() == laneNumber) 
				pointsListsIter.remove();
		}
		customPoints.add(new CustomPoints(laneNumber));
	}
}

class CustomPoints {

	private final int lane;
	private final ArrayList<FitPoint> removeList;
	private final ArrayList<FitPoint> addList;

	public CustomPoints(final int lane) {
		this.lane = lane;
		removeList = new ArrayList<>();
		addList = new ArrayList<>();
	}

	public void addToAddList(final double distance, final double intensity,
		final double sigma)
	{
		addList.add(new FitPoint(distance, intensity, sigma));
	}

	public void addToRemoveList(final double distance, final double intensity) {
		removeList.add(new FitPoint(distance, intensity));
	}

	public int getLane() {
		return lane;
	}

	public ArrayList<FitPoint> getAddList() {
		return addList;
	}

	public ArrayList<FitPoint> getRemoveList() {
		return removeList;
	}

	class FitPoint {

		private final double distance;
		private final double intensity;
		private final double sigma;

		public FitPoint(final double x, final double y) {
			this.distance = x;
			this.intensity = y;
			this.sigma = 0.0;
		}

		public FitPoint(final double x, final double y, final double sigma) {
			this.distance = x;
			this.intensity = y;
			this.sigma = sigma;
		}

		public double getDistance() {
			return distance;

		}

		public double getIntensity() {
			return intensity;
		}

		public double getSigma() {
			return sigma;
		}
	}
}
