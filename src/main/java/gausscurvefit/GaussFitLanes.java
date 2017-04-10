
package gausscurvefit;
/**
 * Gauss Fit
 * GaussFitLanes.java
 * author: Rick Ziraldo, 2017
 * The /University of Texas at Dallas, Richardson, TX
 * http://www.utdallas.edu
 *
 * Feature:   Fitting of Gaussian profiles along gel lanes
 * v is a tool for fitting gaussian profiles and estimating
 * the profile parameters on selected lanes in gel electrophoresis images.
 *
 *    The GaussianArrayCurveFitter class is implemented using
 *    Abstract classes from Apache Commons project
 *
 *    The source code is maintained and made available on GitHub
 *    https://github.com/rickud/gauss-curve-fit
 */

import io.scif.services.DatasetIOService;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Properties;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

import net.imagej.ImageJ;
import net.imagej.display.ImageDisplayService;
import net.imagej.ops.OpService;
import net.imagej.table.DefaultGenericTable;
import net.imagej.table.DefaultTableDisplay;
import net.imagej.table.GenericColumn;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math.analysis.UnivariateRealFunction;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.analysis.integration.TrapezoidIntegrator;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.fitting.WeightedObservedPoints;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math3.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math3.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;
import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.convert.ConvertService;
import org.scijava.display.Display;
import org.scijava.display.DisplayService;
import org.scijava.event.EventService;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import gausscurvefit.GaussianArrayCurveFitter.ParameterGuesser;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.ProfilePlot;
import ij.gui.Roi;
import ij.gui.TextRoi;
import ij.io.Opener;
import ij.measure.Calibration;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;
import ij.util.Tools;

@Plugin(type = Command.class, headless = true,
	menuPath = "Plugins>Gel Tools>Gauss Fit")
public class GaussFitLanes implements Command {

	@Parameter
	private static LogService log;

	@Parameter
	private EventService eventServ;

	@Parameter
	private StatusService statusServ;

	@Parameter
	private static DatasetIOService datasetIOServ;

	@Parameter
	private ConvertService convertServ;

	@Parameter
	private ImageDisplayService imgDisplayServ;

	@Parameter
	private static DisplayService displayServ;

	@Parameter
	private static OpService ops;

	// Default Parameters
	private Thread mainWindowThread; // thread for the main window
	private Thread plotThread; 			 // thread for plotting

	private boolean setup = true;
	private boolean doFit; // tells the background thread to update
	private boolean fitDone = false;

	private final ImagePlus plotImage = new ImagePlus();
	private ArrayList<MyPlot> plots = new ArrayList<>();
	private Color vLineColor = Color.RED;
	private Color profileColor = Color.BLACK;
	private Color gaussColor = Color.RED;
	private Color bgColor = Color.BLUE;
	private Color fittedColor = new Color(0, 180, 0);
	
	private ImagePlus imp;
	private MainDialog md;
	private String version;

	private final int rows = 2; // Number of plot Rows in display
	private final int cols = 2; // Number of plot Rows in display
	private int nLanes = 5; 
	private int degBG = 3; // Order of Background Polynomial
	private double tolPK = 0.05; // Peak detection tolerance as % of range

	
	
	public void init() {
		final double SW = IJ.getScreenSize().getWidth();
		final double SH = IJ.getScreenSize().getHeight();

		imp = IJ.getImage();
		about();
		md = new MainDialog("[v" + version + "] Gel Lanes Gauss Fitting: " + imp
			.getTitle());
		IJ.wait(50); // delay to make sure ROIs have updated
		plotImage.setTitle("Profiles of " + imp.getShortTitle());
		for (int p : getAllLaneNumbers()) updateProfile(p);
		plotsMontage();
		
		plotImage.show();
		imp.getCanvas().requestFocus();
		final ImageWindow iwin = imp.getWindow();
		final ImageWindow pwin = plotImage.getWindow();
		if (iwin == null || pwin == null || md == null) return;

		final Dimension imageSize = iwin.getSize();
		final Dimension plotSize = pwin.getSize();
		final Point imageLoc = iwin.getLocation();

		iwin.setLocation(0, 100);
		pwin.setLocation(imageSize.width + 30, imageLoc.y);
		iwin.getCanvas().requestFocus();

		// thread for plotting in the background
		plotThread = new Thread(this, "Dynamic Plots");
		plotThread.setPriority(Math.max(plotThread.getPriority() - 3,
			Thread.MIN_PRIORITY));
		plotThread.start();
	}

	@Override
	public void run() {
		if (setup) {
			init();
			setup = false;
		}

		synchronized (this) {
			// if (doUpdate) {
			// doUpdate = false; //and loop again
			// } else {
			// try {wait();} //notify wakes up the thread
			// catch(InterruptedException e) { //interrupted tells the thread to exit
			// return;
			// }
			// }
		}
	}

	private int[] getAllLaneNumbers() {
		ArrayList<Roi> rois =   md.getRois();
		Iterator<Roi> roisIter = rois.iterator();
		int[] laneNumbers = new int[rois.size()];
		int i = 0;
		while (roisIter.hasNext()) {
			Roi roi = roisIter.next();
			if (roi != null) {
				String name  = roi.getName();
				laneNumbers[i] = Integer.parseInt(name.substring(5));
				i++;
			}
		}
		return ArrayUtils.subarray(laneNumbers, 0, i);
	}

	/**
	 * Profile data from Roi, ready for fitting (null if not possible)
	 *
	 * @param roi
	 */
	private DataSeries getLaneProfile(final Roi roi) {
		if (md.getRois().size() == 0) return null;
		imp.setRoi(roi);
		final ProfilePlot profileP = new ProfilePlot(imp, true); // get the profile
		final RealVector profile = new ArrayRealVector(profileP.getProfile());
		if (profile.getDimension() < 2) return null;
		imp.killRoi();
		double y0 = roi.getBounds().getMinY();
		double[] y = new double[profile.getDimension()];
		for (int p = 0;  p < y.length; p++) y[p] = y0 + p;
		return new DataSeries("Lane Profile", DataSeries.PROFILE, new ArrayRealVector(y), profile, Color.BLACK);
	}

	
	private void updateProfile(int laneNumber) {
		final Roi roi = md.getRoi("Lane " + laneNumber);
		if (roi != null) {
			//Plot the profiles
			final DataSeries profile = getLaneProfile(roi);
			
			final String xLabel = "Distance (px)";
			final String yLabel = "Grayscale Value";
			MyPlot newPlot = new MyPlot(roi.getName(), xLabel, yLabel);
			newPlot.addDataSeries(profile);
			newPlot.updatePlot();
			
			Iterator<MyPlot> plotIter = plots.iterator();
			while (plotIter.hasNext()) {
				MyPlot plot = plotIter.next();
				if (plot.getNumber() == newPlot.getNumber()) plotIter.remove();
			}
			plots.add(newPlot);
		}
	}
	
	private void setVLine(int laneNumber, double x) {
		Iterator<MyPlot> plotIter = plots.iterator();
		while (plotIter.hasNext()) {
			MyPlot plot = plotIter.next();
			boolean vline = false;
			if (plot.getNumber() == laneNumber) {
				double[] minmax = plot.getDataMinMax();
				RealVector xvals = new ArrayRealVector(new double[] {x, x});
				RealVector yvals = new ArrayRealVector(new double[] {minmax[2], minmax[3]});
				DataSeries l = new DataSeries("", DataSeries.VLINE, xvals, yvals, vLineColor);
				for (DataSeries d : plot.getDataSeries()) {
					if (d.getType() == DataSeries.VLINE) {
						d.setX(xvals.toArray()); 
						d.setY(yvals.toArray());
						vline = true; break;
					}
				}
				if (!vline) {
					plot.setSelected(true);
					plot.addDataSeries(l); 
				}
				plot.updatePlot();
				break;
			} 
		}
	}
	
	private void removeVLine(int laneNumber) {
		Iterator<MyPlot> plotIter = plots.iterator();
		while (plotIter.hasNext()) {
			MyPlot plot = plotIter.next();
			if (plot.getNumber() == laneNumber) {
				for (DataSeries d : plot.getDataSeries()) {
					if (d.getType() == DataSeries.VLINE) {
						 plot.removeDataSeries(d);
						 plot.setSelected(false);
						 plot.updatePlot();
						 break;
					}
				}
				break;
			} 
		}
	}
	
			
	
	
	/**
	 * Method to output plot collage to plotImage
	 */
	private void plotsMontage(){
		// Construct the plot image
		final int plotSpacing = 5; // black border
		ImageProcessor plot0 = plots.get(0).getPlot().getProcessor(); 
		final int plotW = plot0.getWidth();
		final int plotH = plot0.getHeight();
		final int plotsWidth = plotSpacing + cols * (plotW + plotSpacing);
		final int plotsHeight = plotSpacing + rows * (plotH + plotSpacing);
		final ImageProcessor blank = plot0.duplicate();
		IJ.setForegroundColor(255, 255, 255);
		blank.fill();

		if (plots.size() != md.getRois().size() || plots == null) {
			System.out.println(plots.size() + " plots, " + md.getRois().size() +
				" rois");
			return;
		}
		final ArrayList<ImageProcessor> pageMontages = new ArrayList<>();
		Iterator<MyPlot> plotIter = plots.iterator();
		int psel = 0; // selected page
		while (plotIter.hasNext()) {
			final ImageProcessor plotsMontage = new ColorProcessor(plotsWidth,
				plotsHeight);
			Plot subplot;
			ImageProcessor subplotProcessor;
			String subplotTitle;
			for (int r = 0; r < rows; r++) {
				for (int c = 0; c < cols; c++) {
					if (plotIter.hasNext()) {
						subplot = plotIter.next().getPlot();
						if (subplot.getProcessor() == null) 
							System.out.println("Stop!");
						subplotProcessor = subplot.getProcessor();
						subplotTitle = subplot.getTitle();
					}
					else {
						subplotProcessor = blank;
						subplotTitle = "";
					}
					plotsMontage.insert(subplotProcessor, plotSpacing + c * (plotSpacing +
						plotW), plotSpacing + r * (plotSpacing + plotH));

					plotsMontage.drawString(subplotTitle, plotSpacing + plotW / 2 + c *
						(plotSpacing + plotW), plotSpacing + plotH / 5 + r * (plotSpacing +
							plotH), Color.LIGHT_GRAY);

					if (subplotTitle.equals(md.getRoiSelected())) psel = pageMontages
						.size() + 1;
				}
			}
			pageMontages.add(plotsMontage);
		}
		plotImage.setProcessor(pageMontages.get(0));
		final ImageStack plotStack = plotImage.createEmptyStack();
		final Iterator<ImageProcessor> montageIter = pageMontages.iterator();
		int p = 0;
		while (montageIter.hasNext()) {
			plotStack.addSlice("Lanes " + (rows * cols * p + 1) + "-" + (rows * cols *
				p + rows * cols), montageIter.next(), p);
			p++;
		}
		plotImage.setStack(plotStack);
		if (!md.getRoiSelected().equals("none")) plotImage.setSlice(psel);
		plotImage.updateImage();
		plotImage.show();
		if (plotImage.getRoi() != null) plotImage.killRoi();
	}
	
	

	private ArrayList<ArrayList<DataSeries>> doFit() {
		final String warning =
			"The current plots will be reset and the current fitting data will be lost.";

		if (fitDone) {
			if (md.askUser(warning)) {
				// Reset the plots window
				if (displayServ.getDisplay("Results Display") != null) displayServ
					.getDisplay("Results Display").close();
			}
			else return null;
		}

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
		
		ArrayList<ArrayList<DataSeries>> outputData = new ArrayList<>();
		ArrayList<Roi> rois = md.getRois();
		final Iterator<Roi> roiIter = rois.iterator();
		int progress = 1;
		while (roiIter.hasNext()) {
			final Roi roi = roiIter.next();
			DataSeries profile = getLaneProfile(roi);
			final RealVector xvals = new ArrayRealVector(profile.getX());
			final RealVector yvals = new ArrayRealVector(profile.getY());
			ArrayList<DataSeries> funout = new ArrayList<>();
			final int degbg = md.getDegBG();
			// Tolerance as percentage of the range
			final double tolpk = md.getTolPK() * (yvals.getMaxValue() - yvals
				.getMinValue());

			final WeightedObservedPoints obs = new WeightedObservedPoints();
			for (int o = 0; o < xvals.getDimension(); o++) {
				obs.add(xvals.getEntry(o), yvals.getEntry(o));
			}

			final ParameterGuesser pg = new GaussianArrayCurveFitter.ParameterGuesser(
				obs.toList(), tolpk, degBG);
			final RealVector firstGuess = new ArrayRealVector(pg.guess());
			final LeastSquaresProblem problem = GaussianArrayCurveFitter.create(tolpk,
				degbg).getProblem(obs.toList());

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
			RealVector poly = pars.getSubVector(0, degBG + 2);
					
			PolynomialFunction bg = new PolynomialFunction(poly.getSubVector(1, degBG + 1).toArray());
			funout.add(new DataSeries("Background", DataSeries.BACKGROUND, xvals, bg, Color.BLUE));
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
				Gaussian gauss = new Gaussian(norms.getEntry(gg), means.getEntry(gg), sds.getEntry(gg));
				UnivariateFunction[] functs = {bg, gauss};
				funout.add(new DataSeries("Band " + gg, DataSeries.GAUSS_BG, xvals, functs, Color.RED));
			}
			GaussianArrayBG fitted = new GaussianArrayBG(norms, means, sds, poly);
			funout.add(new DataSeries("Fit", DataSeries.FITTED, xvals, fitted, new Color(100, 100, 10)) );
			outputData.add(funout);

			final RealVector peakAreas = doIntegrate(xvals, norms, means, sds);

			// Prepare columns for Results Table
			colLane.add(roi.getName());
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

			final String output = String.format(roi.getName() + ", RMS: %1$.2f; ",
				optimum.getRMS());
			log.info(output);
			statusServ.showProgress(++progress, rois.size());
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

		final String saveFile = "Fit of " + imp.getTitle() + ".xls";
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
		fitDone = true;
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

	
	/**
	 *
	 */
	public void about() {
		final ClassLoader loader = Thread.currentThread().getContextClassLoader();
		final Properties props = new Properties();
		try (InputStream input = loader.getResourceAsStream("about.properties")) {
			props.load(input);
			version = props.getProperty("version");
			log.info("Gauss Fit - v" + version);
			input.close();
		}
		catch (final IOException e) {
			e.printStackTrace();
			log.info("Gauss Fit - v[unknown]");
		}
	}

	
	public static void main(final String... args) throws Exception {
		// create the ImageJ application context with all available services
		final ImageJ ij = net.imagej.Main.launch(args);
		final ImagePlus iPlus = new Opener().openImage(
			"src//main//resources//sample//All[01-17-2017].tif");
		// display it via ImageJ
		iPlus.show();
		iPlus.getCanvas().setLocation(0, 100);
		// wrap it into an ImgLib image (no copying)
		// final Img image = ImagePlusAdapter.wrap(imp);
		// display it via ImgLib using ImageJ
		// ImageJFunctions.show(image);
		// Launch the "OpenImage" command.
		ij.command().run(GaussFitLanes.class, true);
	}
	
	
	class MyPlot implements Comparable<MyPlot> {
		private Plot plot;
		private final Color selectedBG = new Color(240, 240, 240);
		private final Color unselectedBG = Color.white;
		private boolean selected = false;
		private int number;
		private String title;
		private String xLabel;
		private String yLabel;
		private ArrayList<DataSeries> data = new ArrayList<>();
		double xmin;
		double xmax;
		double ymin;
		double ymax;
		
		public MyPlot(String title, String xLabel, String yLabel) {
			this.title = title;
			this.xLabel = xLabel;
			this.yLabel = yLabel;
			this.number = Integer.parseInt(title.substring(5));
		}
	
		public void updatePlot() {
			Plot newPlot = new Plot(title, xLabel, yLabel);
			newPlot.useTemplate(plot, Plot.COPY_LABELS+Plot.COPY_LEGEND+Plot.COPY_SIZE);
			if (plot != null) plot.dispose();
			plot = newPlot;
			if (selected) plot.setBackgroundColor(selectedBG);
			else plot.setBackgroundColor(unselectedBG);
			for (DataSeries d : data) {
				plot.setColor(d.getColor());
				plot.addPoints(d.getX(), d.getY(), Plot.LINE);
			}
			plot.setLimitsToFit(true);
		}
		
		private void sortDataSeries() {
			Collections.sort(data);
			Iterator<DataSeries> dataIter = data.iterator();
			int nullIndex = 0;
			xmin = Double.MAX_VALUE; xmax = Double.MIN_VALUE;
			ymin = Double.MAX_VALUE; ymax = Double.MIN_VALUE;
			while (dataIter.hasNext()) {
				DataSeries d = dataIter.next();
				if (d == null) break;
				double xmin2, xmax2, ymin2, ymax2;
				RealVector x = new ArrayRealVector(d.getX());
				RealVector y = new ArrayRealVector(d.getY());
				xmin2 = x.getMinValue(); xmax2=x.getMaxValue();
				ymin2 = y.getMinValue(); ymax2=y.getMaxValue();
				if (xmin2 < xmin) xmin = xmin2; if (xmax2 > xmax) xmax = xmax2;
				if (ymin2 < ymin) ymin = ymin2; if (ymax2 > ymax) ymax = ymax2;			
				nullIndex++;
			}
			data = new ArrayList<>(data.subList(0, nullIndex));
		}
		
		public void addDataSeries(DataSeries d) {
			data.add(d);
			sortDataSeries();
		}
		
		public void addDataSeries(ArrayList<DataSeries> d) {
			data.addAll(d);
			sortDataSeries();
		}
		
		public void removeDataSeries(DataSeries d) {
			data.remove(d);
			sortDataSeries();
		}
		
		public ArrayList<DataSeries> getDataSeries() {
			return data;
		}
		
		public Plot getPlot() {
			return plot;
		}
		
		public int getNumber() {
			return number;
		}
		public double[] getDataMinMax() {
			return new double[] {xmin, xmax, ymin, ymax};
		}
		public void setSelected(boolean s) {
			selected = s;
		}
		
		@Override
		public int compareTo(MyPlot p) {
			return (number - p.getNumber());
		}	
	}
	
	class DataSeries implements Comparable<DataSeries> {
		private Color color; // Plot color
		private RealVector x; 
		private RealVector y;
		private String name; // Name for Legend
		private int type; // Type of function
		
		final static int VLINE = 0;
		final static int PROFILE = 1;
		final static int BACKGROUND = 2;
		final static int GAUSS_BG = 3;
		final static int FITTED = 4;
		
		 
		public DataSeries(final String name, final int type, final RealVector x, final RealVector y, final Color color) {
			this.name = name;
			this.type = type;
			this.x=x;
			this.y=y;
			this.color=color;
			if (x.getDimension()==2 && x.getEntry(0)==x.getEntry(1)){
				this.type = VLINE;
			} else
				this.type = PROFILE;
		}
		
		public DataSeries(final String name, final int type, final RealVector x, final UnivariateFunction[] function, final Color color) {
			this.name = name;
			this.type = type;
			this.x = x;
			this.color = color;
			for (int i = 0; i < function.length; i++) {
				if (i == 0) {
					this.y = new ArrayRealVector(x.map(function[i]));
				} else {
					this.y = y.add(x.map(function[i]));
				}
			}
		}
		
		public DataSeries(final String name, final int type, final RealVector x, final UnivariateFunction function,
			final Color color) {
			this.name = name;
			this.type = type;
			this.x = x;
			this.color = color;
			this.y = new ArrayRealVector(x.map(function));
		}
		
		public int getType() {
			return type;
		}

		public double[] getX() {
			return x.toArray();
		}
		
		public double[] getY() {
			return y.toArray();
		}
		
		public Color getColor() {
			return color;
		}
		
		public String getName() {
			return name;
		}
		
		public void setX(double[] x){
			this.x = new ArrayRealVector(x);
		}
		
		public void setY(double[] y){
			this.y = new ArrayRealVector(y);
		}
		
		@Override
		public int compareTo(DataSeries d) {
			return (type - d.getType());
		}
	}
	
	@SuppressWarnings("serial")
	class MainDialog extends JFrame implements ActionListener, ChangeListener,
		DocumentListener, ItemListener, MouseMotionListener, MouseListener,
		MouseWheelListener
	{

		private boolean auto = true; // AUTO Lane size mode is ON
		private boolean selectionUpdate = false; // active updating is off
		
		
		private String roiSelected = "none";
		private String roiPreviouslySelected = "none";

		private final int IW = imp.getWidth();
		private final int IH = imp.getHeight();

		// Default lane size/offset (Just center 4 lanes in the image)
		private int LW = (int) Math.round(0.8 * IW / nLanes);
		private int LH = (int) Math.round(IH * 0.8);
		private int LSp = Math.round((IW - LW * nLanes) / (nLanes + 1));
		private int LHOff = LSp / 2;
		private int LVOff = (IH - LH) / 2;

		private ArrayList<Roi> rois = new ArrayList<>();

		private JPanel roiButtonsPanel;
		private JRadioButton buttonAuto;
		private JRadioButton buttonManual;
		private ButtonGroup roiButtons;

		private JPanel textPanel;
		private JLabel labelNLanes;
		private JTextField textNLanes;

		private JPanel sliderPanel;
		private JSlider sliderW;
		private JSlider sliderH;
		private JSlider sliderSp;
		private JSlider sliderHOff;
		private JSlider sliderVOff;

		private JPanel buttonPanel;
		private JButton buttonMeasure;
		private JButton buttonClose;
		private JCheckBox chkBoxBands;

		private JPanel settingsPanel;
		private JPanel degPanel;
		private JPanel tolPanel;
		private JLabel labelDegBG;
		private JLabel labelTolPK;
		private JTextField textDegBG;
		private JTextField textTolPK;

		private JPanel dialogPanel;
		private final JFrame frame;

		public MainDialog(final String string) {
			frame = new JFrame(string);
			setupMainDialog();
		}

		private void setupMainDialog() {
			imp.getCanvas().addMouseMotionListener(this);
			imp.getCanvas().addMouseListener(this);
			imp.getCanvas().addMouseWheelListener(this);
			imp.getWindow().addMouseListener(this);
			imp.getWindow().addMouseMotionListener(this);

			roiButtonsPanel = new JPanel();
			buttonAuto = new JRadioButton("Automatic Rectangle Selection");
			buttonAuto.setActionCommand("Auto");
			buttonAuto.setSelected(true);
			buttonAuto.addActionListener(this);
			buttonManual = new JRadioButton("Manual Rectangle Selection");
			buttonManual.setActionCommand("Manual");
			buttonManual.addActionListener(this);
			roiButtons = new ButtonGroup();
			roiButtonsPanel.setLayout(new BoxLayout(roiButtonsPanel,
				BoxLayout.Y_AXIS));
			roiButtons.add(buttonAuto);
			roiButtons.add(buttonManual);
			
			roiButtonsPanel.add(buttonAuto);
			roiButtonsPanel.add(buttonManual);
			roiButtonsPanel.setBorder(new TitledBorder(BorderFactory
				.createEtchedBorder(), "Lane Selection", TitledBorder.LEADING,
				TitledBorder.BELOW_TOP, new Font("Sans", Font.PLAIN, 11)));

			sliderPanel = new JPanel();
			textPanel = new JPanel();
			textPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
			textNLanes = new JTextField(3);
			textNLanes.setText("" + nLanes);
			textNLanes.setVisible(true);
			textNLanes.getDocument().addDocumentListener(this);
			textPanel.add(textNLanes);
			labelNLanes = new JLabel("Number of Lanes");
			textPanel.add(labelNLanes);
			sliderPanel.add(textPanel);

			sliderPanel.setLayout(new BoxLayout(sliderPanel, BoxLayout.Y_AXIS));
			sliderW = makeTitledSlider("Width ( " + LW + " px )", Color.black, 1, IW /
				4, LW);
			sliderPanel.add(sliderW);
			sliderH = makeTitledSlider("Height ( " + LH + " px )", Color.black, IH /
				10, IH, LH);
			sliderPanel.add(sliderH);
			sliderSp = makeTitledSlider("Space ( " + LSp + " px )", Color.black, 1,
				IW / nLanes, LSp);
			sliderPanel.add(sliderSp);
			sliderHOff = makeTitledSlider("Horizontal Offset ( " + LHOff + " px )",
				Color.black, 0, (int) Math.round(IW * 0.9), LHOff);
			sliderPanel.add(sliderHOff);
			sliderVOff = makeTitledSlider("Vertical Offset ( " + LVOff + " px )",
				Color.black, 0, (int) Math.round(IH * 0.9), LVOff);
			sliderPanel.add(sliderVOff);

			buttonPanel = new JPanel();
			buttonPanel.setLayout(new FlowLayout(FlowLayout.TRAILING));
			buttonMeasure = new JButton("Measure");
			buttonMeasure.addActionListener(this);
			buttonPanel.add(buttonMeasure);
			buttonClose = new JButton("Close");
			buttonClose.addActionListener(this);
			buttonPanel.add(buttonClose);

			settingsPanel = new JPanel();
			degPanel = new JPanel();
			tolPanel = new JPanel();
			labelDegBG = new JLabel("Deg BG");
			labelTolPK = new JLabel("Tol PK");
			textDegBG = new JTextField(3);
			textDegBG.setText("" + degBG);
			textTolPK = new JTextField(3);
			textTolPK.setText("" + tolPK);
			chkBoxBands = new JCheckBox("Show Bands");

			degPanel.setLayout(new FlowLayout(FlowLayout.TRAILING));
			degPanel.add(labelDegBG);
			degPanel.add(textDegBG);
			tolPanel.setLayout(new FlowLayout(FlowLayout.TRAILING));
			tolPanel.add(labelTolPK);
			tolPanel.add(textTolPK);
			chkBoxBands.addItemListener(this);
			chkBoxBands.setSelected(false);
			chkBoxBands.setEnabled(false);

			settingsPanel.setLayout(new GridLayout(10, 1));
			settingsPanel.add(degPanel);
			settingsPanel.add(tolPanel);
			settingsPanel.add(chkBoxBands);

			dialogPanel = new JPanel();
			dialogPanel.setBackground(Color.darkGray);
			dialogPanel.setLayout(new BorderLayout());
			dialogPanel.add(roiButtonsPanel, BorderLayout.NORTH);
			dialogPanel.add(sliderPanel, BorderLayout.CENTER);
			dialogPanel.add(settingsPanel, BorderLayout.EAST);
			dialogPanel.add(buttonPanel, BorderLayout.SOUTH);

			frame.addWindowListener(new WindowAdapter() {

				@Override
				public void windowClosing(final WindowEvent e) {
					cleanupAndClose();
					//frame.dispatchEvent(new WindowEvent(frame, WindowEvent.WINDOW_CLOSING));
					log.info("Gauss Fit terminated.");
				}
			});
			frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
			frame.getContentPane().add(dialogPanel);
			frame.setResizable(true);
			frame.validate();
			frame.pack();
			frame.setVisible(true);

			resetAutoROIs(LW, LH, LSp, LHOff, LVOff, nLanes);
			reDrawROIs(imp, "none");
			imp.killRoi();
		}

		private void cleanupAndClose() {
			// Remove Listeners on imp
			if (imp != null) {
				imp.getCanvas().removeMouseListener(this);
				imp.getCanvas().removeMouseMotionListener(this);
				imp.getCanvas().removeMouseWheelListener(this);
				imp.getWindow().removeMouseListener(this);
				imp.getWindow().removeMouseMotionListener(this);
				imp.setOverlay(null);
			}
			if (plotImage.getWindow() != null) plotImage.getWindow().close();
			final List<Display<?>> displays = displayServ.getDisplays();
			for (final Display<?> d : displays)
				d.close();
		}

		private boolean askUser(final String question) {
			final GenericDialog gd = new GenericDialog("WARNING!");
			gd.addMessage(question);
			gd.showDialog();
			if (gd.wasOKed()) return true;
			return false;
		}

		private JSlider makeTitledSlider(final String string, final Color color,
			final int minVal, final int maxVal, final int val)
		{
			final JSlider slider = new JSlider(SwingConstants.HORIZONTAL, minVal,
				maxVal, val);
			final TitledBorder tb = new TitledBorder(BorderFactory
				.createEtchedBorder(),
				// empty,
				"", TitledBorder.CENTER, TitledBorder.BELOW_TOP, new Font("Sans",
					Font.PLAIN, 11));
			tb.setTitle(string);
			tb.setTitleJustification(TitledBorder.LEFT);
			tb.setTitleColor(color);
			slider.setBorder(tb);
			slider.setMajorTickSpacing((maxVal - minVal) / 10);
			slider.setPaintTicks(true);
			slider.addChangeListener(this);
			return slider;
		}

		private void setSliderTitle(final JSlider slider, final String str) {
			final TitledBorder tb = (TitledBorder) slider.getBorder();
			tb.setTitle(str);
			slider.setBorder(tb);
		}

		private void setSliderPanelEnabled(final boolean enabled) {
			if (enabled) {
				textNLanes.setEnabled(true);
				sliderW.setEnabled(true);
				sliderH.setEnabled(true);
				sliderSp.setEnabled(true);
				sliderHOff.setEnabled(true);
				sliderVOff.setEnabled(true);
			}
			else {
				textNLanes.setEnabled(false);
				sliderW.setEnabled(false);
				sliderH.setEnabled(false);
				sliderSp.setEnabled(false);
				sliderHOff.setEnabled(false);
				sliderVOff.setEnabled(false);
			}
		}

		private void reDrawROIs(final ImagePlus imgPlus, final String roiName) {
			if (!auto) {
				// HouseKeeping: Sort the rois based on name and remove null elements
				final Comparator<Roi> nameComparator = new Comparator<Roi>() {

					@Override
					public int compare(final Roi r1, final Roi r2) {
						if (r1 == null && r2 == null) {
							return 0;
						}
						if (r1 == null) {
							return -1;
						}
						if (r2 == null) {
							return 1;
						}
						final int n1 = Integer.parseInt(r1.getName().substring(5));
						final int n2 = Integer.parseInt(r2.getName().substring(5));
						final int comp = Integer.compare(n1, n2);
						if (comp != 0) {
							return comp;
						}
						return 0;
					}
				};
				Collections.sort(rois, nameComparator);
				final Iterator<Roi> roisIter = rois.iterator();
				int nullIndex = 0;
				while (roisIter.hasNext()) {
					if (roisIter.next() == null) break;
					nullIndex++;
				}
				rois = new ArrayList<>(rois.subList(0, nullIndex));
			}
			final Overlay overlay = new Overlay();
			final Iterator<Roi> roiIter = rois.iterator();
			while (roiIter.hasNext()) {
				final Roi roi = roiIter.next();
				final String label = roi.getName().substring(5);
				final Font font = new Font("SansSerif", Font.BOLD, 20);
				final TextRoi labelRoi = new TextRoi(0, 0, label, font);
				labelRoi.setAntialiased(true);
				final double lh = font.getSize() * 1.4;
				final double lw = font.getSize() * 1.4;
				final double rh = roi.getBounds().getHeight();
				final double rw = roi.getBounds().getWidth();
				final double x0 = roi.getBounds().getX();
				final double y0 = roi.getBounds().getY();
				if (lh > rh || lw > rw) labelRoi.setLocation(x0 + 0.1 * lw, y0 + rh +
					1.1 * lh);
				else labelRoi.setLocation(x0 + 0.2 * lw, y0 + rh - 1.1 * lh);
				if (roi.getName().equals(roiName)) {
					roi.setStrokeColor(Color.YELLOW);
					roi.setStrokeWidth(3);
					if (buttonManual.isSelected() && !roiName.equals("none")) {
						imgPlus.setRoi(roi);
					}
				}
				else {
					roi.setStrokeWidth(1);
					roi.setStrokeColor(Color.RED);
					labelRoi.setFillColor(Color.RED);
				}
				overlay.add(roi);
				overlay.add(labelRoi);
				if (buttonAuto.isSelected()) imgPlus.killRoi();
			}
			imgPlus.setOverlay(overlay);
		}

		private void resetAutoROIs(final int lw, final int lh, final int lsp,
			final int lhoff, final int lvoff, final int nlanes)
		{
			rois = new ArrayList<>();
			for (int i = 0; i < nlanes; i++) {
				final Roi roi = new Roi(lhoff + lw * i + lsp * i, lvoff, lw, lh);
				roi.setName("Lane " + (i + 1));
				rois.add(roi);
			}
		}
		
		private void changeLaneNumber() {
			nLanes = getNLanes();
			if (nLanes == -1) return;
			resetAutoROIs(LW, LH, LSp, LHOff, LVOff, nLanes);
			reDrawROIs(imp, "none");
			plots = new ArrayList<>();
			for (int i = 1; i <= nLanes; i++) updateProfile(i);
			plotsMontage();
		}

		public int getNLanes() {
			try {
				return Integer.parseInt(textNLanes.getText().trim());
			}
			catch (final NumberFormatException e1) {
				return -1;
			}
		}

		public int getDegBG() {
			try {
				degBG = Integer.parseInt(textDegBG.getText().trim());
				return degBG;
			}
			catch (final NumberFormatException e1) {
				return -1;
			}
		}

		/**
		 * @return the currently selected @param tolpk, the tolerance on the peak
		 *         selection
		 */
		public double getTolPK() {
			try {
				tolPK = Double.parseDouble(textTolPK.getText().trim());
				return tolPK;
			}
			catch (final NumberFormatException e1) {
				return 0.0;
			}
		}

		/**
		 * @return the current ROI set
		 */
		public ArrayList<Roi> getRois() {
			return rois;
		}

		public Roi getRoi(final String title) {
			final Iterator<Roi> roiIter = rois.iterator();
			while (roiIter.hasNext()) {
				final Roi roi = roiIter.next();
				if (roi.getName().equals(title)) return roi;
			}
			return null;
		}

		public String getRoiSelected() {
			return roiSelected;
		}

		public String getRoiPreviouslySelected() {
			return roiPreviouslySelected;
		}

		@Override
		public synchronized void stateChanged(final ChangeEvent e) {
			final JSlider slider = (JSlider) e.getSource();
			if (slider == sliderW) {
				LW = sliderW.getValue();
				final String str = "Width ( " + LW + " px )";
				setSliderTitle(sliderW, str);
			}
			else if (slider == sliderH) {
				LH = sliderH.getValue();
				final String str = "Height ( " + LH + " px )";
				setSliderTitle(sliderH, str);
			}
			else if (slider == sliderSp) {
				LSp = sliderSp.getValue();
				final String str = "Spacing ( " + LSp + " px )";
				setSliderTitle(sliderSp, str);
			}
			else if (slider == sliderHOff) {
				LHOff = sliderHOff.getValue();
				final String str = "Horizontal Offset ( " + LHOff + " px )";
				setSliderTitle(sliderHOff, str);
			}
			else if (slider == sliderVOff) {
				LVOff = sliderVOff.getValue();
				final String str = "Vertical Offset ( " + LVOff + " px )";
				setSliderTitle(sliderVOff, str);
			}
			resetAutoROIs(LW, LH, LSp, LHOff, LVOff, nLanes);
			reDrawROIs(imp, "none");
			plots = new ArrayList<>();
			for (int i = 1; i <= nLanes; i++) updateProfile(i);
			plotsMontage();
		}

		// TextField Listeners
		@Override
		public void changedUpdate(final DocumentEvent e) {
			changeLaneNumber();
		}

		@Override
		public void removeUpdate(final DocumentEvent e) {
			changeLaneNumber();
		}

		@Override
		public void insertUpdate(final DocumentEvent e) {
			changeLaneNumber();
		}

		@Override
		public void actionPerformed(final ActionEvent e) {
			if (e.getSource().equals(buttonMeasure)) {
				ArrayList<ArrayList<DataSeries>> fitted = doFit();
				for (ArrayList<DataSeries> f : fitted) {
					MyPlot plot = plots.get(fitted.indexOf(f));
					plot.addDataSeries(f);
					plot.updatePlot();
				}
				plotsMontage();
				chkBoxBands.setEnabled(true);
			}

			if (e.getSource().equals(buttonClose)) {
				if (askUser("Would you like to exit?")) {
					frame.dispatchEvent(new WindowEvent(frame, WindowEvent.WINDOW_CLOSING));
				}
			}

			if (e.getSource().equals(buttonAuto)) {
				if (!auto) { // NOT already AUTO
					if (askUser(
						"When switching to AUTO mode, the current lane selections will be reset!"))
					{
						auto = true;
						setSliderPanelEnabled(true);
					}
					else buttonManual.setSelected(true);
				}
			}

			if (e.getSource().equals(buttonManual)) {
				auto = false;
				setSliderPanelEnabled(false);
			}

		}

		@Override
		public void itemStateChanged(final ItemEvent e) {
			if (e.getItemSelectable() == chkBoxBands) {
				// TODO Add feature to show/hide horizontal bars for detected peaks in
				// original image
			}
		}

		@Override
		public void mouseDragged(final MouseEvent e) {

		}

		@Override
		public void mouseMoved(final MouseEvent e) {
			if (e.getSource() == imp.getCanvas()) {
				final int x = ((ImageCanvas) e.getSource()).offScreenX(e.getX());
				final int y = ((ImageCanvas) e.getSource()).offScreenX(e.getY());
				statusServ.showStatus("[" + x + ":" + y + "]");
				if (!selectionUpdate) {
					String roiCurrent = "none"; // None selected
					int roiN = 0;
					final Iterator<Roi> roisIter = rois.iterator();
					while (roisIter.hasNext()) {
						final Roi roi = roisIter.next();
						if (roi.contains(x, y)) { 
							roiCurrent = roi.getName(); // This selected
							roiN = Integer.parseInt(roiCurrent.substring(5));
						}
					}
							
					if (!(roiCurrent.equals(roiSelected))) {
						roiPreviouslySelected = roiSelected;
						roiSelected = roiCurrent;
						reDrawROIs(imp, roiCurrent);
						if (!roiPreviouslySelected.equals("none")) {
							int roiPN = Integer.parseInt(roiPreviouslySelected.substring(5));
							removeVLine(roiPN);
						}
						plotsMontage();
						
						if (roiCurrent.equals("none")) {
							if (imp.getRoi() != null) imp.killRoi();						
						}
					} else {
						if (!roiSelected.equals("none")) {
							setVLine(roiN, y);
							plotsMontage();
						}
					}				
				}
			}
		}

		@Override
		public void mouseClicked(final MouseEvent e) {
			if (e.getClickCount() == 2) {
				if (!auto && !roiSelected.equals("none") && !selectionUpdate) {
					// Remove Roi
					final Roi roiRemove = md.getRoi(roiSelected);
					if (roiRemove == null) return;
					if (roiRemove.isHandle(e.getX(), e.getY()) != -1) {
						// Might be dragging existing roi
						return;
					}

					if (askUser("Delete " + roiSelected + "?")) {
						final int x = ((ImageCanvas) e.getSource()).offScreenX(e.getX());
						final int y = ((ImageCanvas) e.getSource()).offScreenX(e.getY());
						final Iterator<Roi> roiIter = rois.iterator();
						while (roiIter.hasNext()) {
							Roi roi = roiIter.next();
							if (roi.contains(x, y)) {
								roiIter.remove();
								// System.out.println(imp.getTitle() + ": " + roiSelected +
								// " - REMOVED");
							}
						}
						reDrawROIs(imp, "none");
						final Iterator<Roi> roiIter2 = rois.iterator();
						plots = new ArrayList<>();
						while (roiIter2.hasNext()) {
							Roi roi2 = roiIter2.next();
							int laneNumber = Integer.parseInt(roi2.getName().substring(5));
							updateProfile(laneNumber);
						}
						plotsMontage();
					}
				}
			}
		}

		@Override
		public void mousePressed(final MouseEvent e) {
			if (!auto) selectionUpdate = true;
		}

		@Override
		public void mouseReleased(final MouseEvent e) {
			// Add or modify Roi
			if (buttonManual.isSelected()) {
				if (selectionUpdate) {
					final Roi roiNew = imp.getRoi();
					if (roiNew != null && roiNew.getBounds().getWidth() * roiNew
						.getBounds().getHeight() > 100)
					{

						if (roiSelected.equals("none")) {
							// Creating New ROI
							// Find an available integer for lane name
							int i = 0;
							while (++i <= rois.size() + 1) {
								final Iterator<Roi> roisIter = rois.iterator();
								boolean contains = false;
								while (roisIter.hasNext()) {
									if (i == Integer.parseInt(roisIter.next().getName().substring(
										5)))
									{
										contains = true;
										break;
									}
								}
								if (!contains) {
									roiNew.setName("Lane " + i);
									roiSelected = roiNew.getName();
									break;
								}
							}
							rois.add(roiNew);
							// System.out.println(imp.getTitle() + ": " + roiSelected +
							// " - ADDED");
						}
						else {
							// System.out.println("Modify");
							// Moving/Resizing a specific ROI
							final Iterator<Roi> roisIter = rois.iterator();
							while (roisIter.hasNext()) {
								final Roi roi = roisIter.next();
								if (roi != null) {
									if (roi.getName().equals(roiSelected)) {
										roiNew.setName(roi.getName());
										rois.set(rois.indexOf(roi), roiNew);
										break;
									}
								}
							}
						}
						reDrawROIs(imp, roiSelected);
						final Iterator<Roi> roiIter2 = rois.iterator();
						plots = new ArrayList<>();
						while (roiIter2.hasNext()) {
							Roi roi = roiIter2.next();					
							updateProfile(Integer.parseInt(roi.getName().substring(5)));
						}
						plotsMontage();
					}
				}
			}
			selectionUpdate = false;
		}

		@Override
		public void mouseEntered(final MouseEvent e) {
			// TODO Auto-generated method stub

		}

		@Override
		public void mouseExited(final MouseEvent e) {
			// TODO Auto-generated method stub

		}

		@Override
		public void mouseWheelMoved(final MouseWheelEvent e) {
			if (e.getSource() == imp.getCanvas() && imp.getImageStackSize() > 1) {
				final Iterator<Roi> roiIter = rois.iterator();
				while (roiIter.hasNext()) {
					Roi roi = roiIter.next();					
					updateProfile(Integer.parseInt(roi.getName().substring(5)));
				}
				plotsMontage();
			}
		}
	}
}
