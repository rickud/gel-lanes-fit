
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
	private final String title; // for the datafile, typically the main image's
															// title

	private ArrayList<DataSeries> inputData;
	private ArrayList<CustomPeaks> customPeaks;
	private ArrayList<ArrayList<Peak>> fittedPeaks;

	public Fitter(final Context context, final String title, final int degBG,
		final double tolPK)
	{
		context.inject(this);
		this.title = title;
		this.degBG = degBG;
		this.tolPK = tolPK;
		this.inputData = new ArrayList<>();
		this.customPeaks = new ArrayList<>();
		this.fittedPeaks = new ArrayList<>();
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
		fittedPeaks = new ArrayList<>();
		int progress = 1;
		int datasetNumber = 0;
		final Iterator<DataSeries> dataInputIter = inputData.iterator();
		while (dataInputIter.hasNext()) {
			datasetNumber++;
			final DataSeries profile = dataInputIter.next();
			final int lane = profile.getLane();
			final RealVector xvals = new ArrayRealVector(profile.getX());
			final RealVector yvals = new ArrayRealVector(profile.getY());
			final ArrayList<DataSeries> funout = new ArrayList<>();
			final ArrayList<Peak> fittedPeaksThisLane = new ArrayList<>();

			// Tolerance as percentage of the range
			final double tolpk = tolPK * (yvals.getMaxValue() - yvals.getMinValue());

			final WeightedObservedPoints obs = new WeightedObservedPoints();
			for (int o = 0; o < xvals.getDimension(); o++) {
				obs.add(xvals.getEntry(o), yvals.getEntry(o));
			}

			final ParameterGuesser pg = new GaussianArrayCurveFitter.ParameterGuesser(
				obs.toList(), tolpk, degBG);
			RealVector firstGuess = new ArrayRealVector(pg.guess());

			// Initial Guess
			RealVector norms0 = new ArrayRealVector();
			RealVector means0 = new ArrayRealVector();
			RealVector sds0 = new ArrayRealVector();
			for (int b = degBG + 2; b < firstGuess.getDimension(); b += 3) {
				// Initial Guess
				norms0 = norms0.append(firstGuess.getEntry(b));
				means0 = means0.append(firstGuess.getEntry(b + 1));
				sds0 = sds0.append(firstGuess.getEntry(b + 2));
			}

			// Check CustomPoint for the lane and adjust initial guess
			for (final CustomPeaks cp : customPeaks) {
				if (cp.getLane() == lane) {
					for (final Peak fp : cp.getList()) {
						final double mean = fp.getDistance();
						int idx = firstGuess.mapSubtract(mean).mapToSelf(new Abs())
							.getMinIndex();
						if (FastMath.abs(firstGuess.getEntry(idx) - mean) <= 5){
							// the point is already detected, ONLY replace the SD
							firstGuess.setEntry(idx + 1, fp.getSigma());
							int idx2 = (idx - degBG - 1)/3;
							sds0.setEntry(idx2, fp.getSigma());
						}
						else {
							// Place the new peak in order of mean, distance y from top 
							for (int m  = 0; m < means0.getDimension(); m++) {
								if (fp.getDistance() < means0.getEntry(m)) {
									norms0 = norms0.getSubVector(0,m).append(fp.getIntensity()).append(norms0.getSubVector(m,norms0.getDimension()-m));
									means0 = means0.getSubVector(0,m).append(fp.getDistance()).append(means0.getSubVector(m,means0.getDimension()-m));
									sds0 = sds0.getSubVector(0,m).append(fp.getSigma()).append(sds0.getSubVector(m,sds0.getDimension()-m));
									firstGuess = firstGuess.getSubVector(0,degBG+2+3*m).append(fp.getIntensity()).append(fp
										.getDistance()).append(fp.getSigma()).append(firstGuess.getSubVector(degBG+2+3*m,(firstGuess.getDimension()-(degBG+2+3*m))));
									break;
								}
								// mean is larger than all existing ones
								if (m == means0.getDimension()-1) {
									norms0 = norms0.append(fp.getIntensity());
									means0 = means0.append(fp.getDistance());
									sds0 = sds0.append(fp.getSigma());
									firstGuess = firstGuess.append(fp.getIntensity()).append(fp.getDistance()).append(fp.getSigma());
									break;
								}
							}
						}
					}

					// FIXME Remove points previously added MANUALLY
					// Does not make sense to remove points that where detected automatically.
					// Better to increase threshold to avoid detecting noise.
					// Use remove feature to remove a previously added custom peak.
				}
			}

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
				fittedPeaksThisLane.add(new Peak(lane,pars.getEntry(b),pars.getEntry(b + 1),pars.getEntry(b + 2)));
			}
			fittedPeaks.add(fittedPeaksThisLane);
			final PolynomialFunction bg = new PolynomialFunction(poly.getSubVector(1,
				degBG + 1).toArray());
			funout.add(new DataSeries("Background", lane, DataSeries.BACKGROUND,
				xvals, bg, Plotter.bgColor));

			for (int gg = 1; gg < norms.getDimension(); gg++) {
				final Gaussian gauss = new Gaussian(norms.getEntry(gg), means.getEntry(
					gg), sds.getEntry(gg));
				final UnivariateFunction[] functs = { bg, gauss };
				funout.add(new DataSeries("Band " + gg, lane, DataSeries.GAUSS_BG,
					xvals, functs, Plotter.gaussColor));
			}
			final GaussianArray fitted = new GaussianArray(norms, means, sds,
				poly);
			funout.add(new DataSeries("Fit", lane, DataSeries.FITTED, xvals, fitted,
				Plotter.fittedColor));
			outputData.add(funout);

			final RealVector peakAreas = doIntegrate(xvals, norms, means, sds);

			// Prepare columns for Results Table
			colLane.add("Lane " + datasetNumber);
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

			final String output = String.format("Lane " + datasetNumber +
				", RMS: %1$.2f; ", optimum.getRMS());
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
			final Gaussian gauss = new Gaussian(norms.getEntry(i), means.getEntry(i),
				sds.getEntry(i));
			areas = areas.append(ti.integrate(Integer.MAX_VALUE, gauss, xvals
				.getMinValue(), xvals.getMaxValue()));
		}
		return areas;
	}

	public void addCustomPeak(int lane, Peak peak) {
		for (final CustomPeaks cp : customPeaks) {
			if (cp.getLane() == lane) {
				cp.addToList(peak);
				return;
			}
		}
		CustomPeaks cp = new CustomPeaks(lane);
		cp.addToList(peak);
		customPeaks.add(cp);
		return;
	}

	public CustomPeaks getCustomPeaks(int lane) {
		Iterator<CustomPeaks> cpIter = customPeaks.iterator();
		while (cpIter.hasNext()){
			CustomPeaks cp = cpIter.next();
			if (cp.getLane() == lane) {
				return cp;
			}
		}
		return null;
	}

	/**
	 * return fittedPeaks
	 */
	public ArrayList<Peak> getFittedPeaks(int lane) {
		for (ArrayList<Peak> f : fittedPeaks) {
			if (!f.isEmpty() && f.get(0).getLane() == lane) {
				return f;
			}
		}
		return new ArrayList<>();
	}

	public void removeCustomPeak(int lane, Peak peak) {
		for (final CustomPeaks cp : customPeaks) {
			if (cp.getLane() == lane) {
					cp.removeFromList(peak);
				}
		}
	}

	public void resetCustomPeaks(final int laneNumber) {
		final Iterator<CustomPeaks> pointsListsIter = customPeaks.iterator();
		while (pointsListsIter.hasNext()) {
			final CustomPeaks pointsList = pointsListsIter.next();
			if (pointsList.getLane() == laneNumber) pointsListsIter.remove();
		}
		customPeaks.add(new CustomPeaks(laneNumber));
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

class CustomPeaks {

	private final int lane;
	private final ArrayList<Peak> addList;

	public CustomPeaks(final int lane) {
		this.lane = lane;
		addList = new ArrayList<>();
	}

	public void addToList(Peak peak) {
		addList.add(peak);
	}

	public void removeFromList(final Peak peak)
	{
		Iterator<Peak> peakIter = addList.iterator(); 
		while (peakIter.hasNext()) {
			double diff = FastMath.abs(peak.getDistance()-peakIter.next().getDistance());
			if (diff <= 5) peakIter.remove();
		}
	}

	public int getLane() {
		return lane;
	}

	public ArrayList<Peak> getList() {
		return addList;
	}
}

class Peak {

		private final double distance;
		private final double intensity;
		private final double sigma;
		private final int lane;

		public Peak(int lane, final double intensity, final double distance) {
			this.intensity = intensity;
			this.distance = distance;
			this.sigma = 0.0;
			this.lane = lane;
		}

		public Peak(int lane, final double intensity, final double distance, final double sigma) {
			this.intensity = intensity;
			this.distance = distance;
			this.sigma = sigma;
			this.lane = lane;
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
		
		public double getLane() {
			return lane;
		}
}

