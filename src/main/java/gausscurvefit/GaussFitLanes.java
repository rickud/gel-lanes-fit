
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
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.FileInputStream;
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
import ij.gui.Line;
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
	private Thread mainThread; // thread for the main window
	private Thread plotThread; // thread for plotting
	private boolean setup = true;
	private boolean doFit; // tells the background thread to update
	private final ImagePlus plotImage = new ImagePlus();

	private ArrayList<Plot> plots;
	private ImagePlus imp;
	private MainDialog md;

	private final int rows = 2; // Number of plot Rows in display
	private final int cols = 2; // Number of plot Rows in display
	private int nLanes = 8;

	private int degBG = 3; // Order of Background Polynomial
	private double tolPK = 0.05; // Peak detection tolerance as % of range

	public void init() {
		imp = IJ.getImage();
		about();
		md = new MainDialog("Gel Lanes Gauss Fitting: " + imp.getTitle());
		IJ.wait(50); // delay to make sure ROIs have updated
		plotImage.setTitle("Profiles of " + imp.getShortTitle());
		updatePlots();
		plotImage.show();
		imp.getCanvas().requestFocus();
		final ImageWindow iwin = imp.getWindow();
		final ImageWindow pwin = plotImage.getWindow();
		if (iwin == null || pwin == null || md == null) return;

		final Dimension imageSize = iwin.getSize();
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

	/**
	 * Profile data from Roi, ready for fitting (null if not possible)
	 *
	 * @param roi
	 */
	private RealMatrix getLaneProfile(final Roi roi) {
		if (md.getRois().size() == 0) return null;
		imp.setRoi(roi);
		final ProfilePlot profileP = new ProfilePlot(imp, true); // get the profile
		final RealVector profile = new ArrayRealVector(profileP.getProfile());
		if (profile.getDimension() < 2) return null;
		imp.killRoi();
		// the following code is mainly for x calibration
		final Calibration cal = imp.getCalibration();
		final RealVector x = calibrateX(profile, roi, cal);
		final RealMatrix output = MatrixUtils.createRealMatrix(new double[][] { x
			.toArray(), profile.toArray() });
		return output;
	}

	/**
	 * Method for x calibration
	 **/
	private RealVector calibrateX(final RealVector y, final Roi roi,
		final Calibration cal)
	{
		double xInc = 1;
		if (roi.getType() == Roi.LINE) {
			final Line line = (Line) roi;
			if (cal != null) {
				final double dx = cal.pixelWidth * (line.x2 - line.x1);
				final double dy = cal.pixelHeight * (line.y2 - line.y1);
				final double length = Math.sqrt(dx * dx + dy * dy);
				xInc = length / (y.getMaxIndex());
			}
		}
		else if (roi.getType() == Roi.RECTANGLE) {
			if (cal != null) {
				xInc = roi.getBounds().getHeight() * cal.pixelHeight / (y
					.getDimension());
			}
		}
		else return null;

		final double[] x = new double[y.getDimension()];// create the x axis
		for (int i = 0; i < y.getDimension(); i++)
			x[i] = i * xInc;
		return new ArrayRealVector(x);
	}

	/**
	 * Method for Output plot collage
	 */
	private void updatePlots() {
		if (md.getRois().size() == 0) return;

		if (plots == null) { // Plots have been reset
			plots = new ArrayList<>();

			Iterator<Roi> roisIter1 = md.getRois().iterator();
			while (roisIter1.hasNext()) {
				final Roi roi = roisIter1.next();
				if (roi != null) {
					final RealVector profile = new ArrayRealVector(getLaneProfile(roi)
						.getRow(1));
					if (profile.getDimension() < 2) return;

					final Calibration cal = imp.getCalibration();
					String xUnit;
					if (cal.getUnit() == null) xUnit = "pixels";
					else xUnit = cal.getUnit();
					final RealVector x = new ArrayRealVector(getLaneProfile(roi).getRow(
						0));
					final String xLabel = "Distance (" + xUnit + ")";
					final String yLabel = (cal.getValueUnit() != null && !cal
						.getValueUnit().equals("Gray Value")) ? "Value (" + cal
							.getValueUnit() + ")" : "Value";
					Plot plot = new Plot(roi.getName(), xLabel, yLabel, x.toArray(),
						profile.toArray());

					final double fixedMin = ProfilePlot.getFixedMin();
					final double fixedMax = ProfilePlot.getFixedMax();
					if (fixedMin != 0 || fixedMax != 0) {
						final double[] a = Tools.getMinMax(x.toArray());
						plot.setLimits(a[0], a[1], fixedMin, fixedMax);
					}
					plots.add(plot);
				}
			}
		}
		// Construct the plot image
		final Color selectedBG = new Color(255, 218, 185);
		final Color unselectedBG = Color.white;
		Iterator<Plot> plotIter = plots.iterator();
		while (plotIter.hasNext()) {
			Plot plot = plotIter.next();
			if (plot.getTitle().equals(md.getPlotSelected())) {
				plot.setBackgroundColor(selectedBG);
			}
			else if (plot.getTitle().equals(md.getPlotPreviouslySelected())) {
				plot.setBackgroundColor(unselectedBG);
			}
			plot.updateImage();
		}

		final int plotSpacing = 5; // black border
		final int plotW = plots.get(0).getProcessor().getWidth();
		final int plotH = plots.get(0).getProcessor().getHeight();
		final int plotsWidth = plotSpacing + cols * (plotW + plotSpacing);
		final int plotsHeight = plotSpacing + rows * (plotH + plotSpacing);
		final ImageProcessor blank = plots.get(0).getProcessor().duplicate();
		IJ.setForegroundColor(255, 255, 255);
		blank.fill();

		if (plots.size() != md.getRois().size() || plots == null) {
			System.out.println(plots.size() + " plots, " + md.getRois().size() +
				" rois");
			return;
		}
		final ArrayList<ImageProcessor> pageMontages = new ArrayList<>();
		plotIter = plots.iterator();
		int psel = 0; // selected page
		while (plotIter.hasNext()) {
			final ImageProcessor plotsMontage = new ColorProcessor(plotsWidth,
				plotsHeight);
			Plot subplot;
			ImageProcessor subplotProcessor;
			String subplotTitle;
			for (int c = 0; c < cols; c++) {
				for (int r = 0; r < rows; r++) {
					if (plotIter.hasNext()) {
						subplot = plotIter.next();
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

					if (subplotTitle.equals(md.getPlotSelected())) psel = pageMontages
						.size() + 1;
				}
			}
			pageMontages.add(plotsMontage);
		}
		plotImage.setProcessor(pageMontages.get(0));
		final ImageStack plotStack = plotImage.createEmptyStack();
		Iterator<ImageProcessor> montageIter = pageMontages.iterator();
		int p = 0;
		while (montageIter.hasNext()) {
			plotStack.addSlice("Lanes " + (rows * cols * p + 1) + "-" + (rows * cols *
				p + rows * cols), montageIter.next(), p);
			p++;
		}
		plotImage.setStack(plotStack);
		if (!md.getPlotSelected().equals("none")) plotImage.setSlice(psel);
		plotImage.updateImage();
		if (plotImage.getRoi() != null) plotImage.killRoi();
	}

	private void doFit() {
		// Reset the plots window
		plots = null;
		updatePlots();

		// Results Table Columns
		final ArrayList<String> colLane = new ArrayList<>();
		final ArrayList<Integer> colBand = new ArrayList<>();
		RealVector colDistance = new ArrayRealVector();
		RealVector colAmplitude = new ArrayRealVector();
		RealVector colFWHM = new ArrayRealVector();
		RealVector colArea = new ArrayRealVector();

		Iterator<Plot> plotIter = plots.iterator();
		int progress = 1;
		while (plotIter.hasNext()) {
			Plot plot = plotIter.next();
			if (log.getLevel() == 0) log.initialize();

			Roi roi = md.getRoi(plot.getTitle());
			final RealVector xvals = new ArrayRealVector(getLaneProfile(roi).getRow(
				0));
			final RealVector yvals = new ArrayRealVector(getLaneProfile(roi).getRow(
				1));
			RealVector bg = new ArrayRealVector();

			final int degbg = md.getDegBG();
			// Tolerance as percentage of the range
			final double tolpk = md.getTolPK() * (yvals.getMaxValue() - yvals
				.getMinValue());

			final WeightedObservedPoints obs = new WeightedObservedPoints(); // All Y
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
			RealVector poly0 = new ArrayRealVector();

			// After fitting
			RealVector norms = new ArrayRealVector();
			RealVector means = new ArrayRealVector();
			RealVector sds = new ArrayRealVector();
			RealVector gauss = new ArrayRealVector(); // Each Gaussian for plotting
			RealVector poly = new ArrayRealVector(); // Polynomial for plotting
			RealVector fitted = new ArrayRealVector(); // Train of Gaussians with BG
			RealVector guessed = new ArrayRealVector();

			poly = pars.getSubVector(0, degBG + 2);
			poly0 = firstGuess.getSubVector(0, degBG + 2);

			bg = xvals.map(new PolynomialFunction(poly.getSubVector(1, degBG + 1)
				.toArray()));
			plot.setColor(Color.blue);
			plot.addPoints(xvals.toArray(), bg.toArray(), PlotWindow.LINE);

			for (int b = degBG + 2; b < pars.getDimension(); b += 3) {
				// Initial Guess
				norms0 = norms0.append(firstGuess.getEntry(b));
				means0 = means0.append(firstGuess.getEntry(b + 1));
				sds0 = sds0.append(firstGuess.getEntry(b + 2));

				// After fitting
				norms = norms.append(pars.getEntry(b));
				means = means.append(pars.getEntry(b + 1));
				sds = sds.append(pars.getEntry(b + 2));
				gauss = xvals.map(new Gaussian(pars.getEntry(b), pars.getEntry(b + 1),
					pars.getEntry(b + 2)));
				gauss = gauss.add(bg); // add the background Polynomial
				plot.setColor(Color.red);
				plot.addPoints(xvals.toArray(), gauss.toArray(), PlotWindow.LINE);
			}
			fitted = xvals.map(new GaussianArrayBG(norms, means, sds, poly));
			guessed = xvals.map(new GaussianArrayBG(norms0, means0, sds0, poly0));

			plot.setColor(Color.blue);
			plot.addPoints(means0.toArray(), norms0.toArray(), PlotWindow.CROSS);
			plot.setColor(new Color(0, 128, 0));
			plot.addPoints(xvals.toArray(), fitted.toArray(), PlotWindow.LINE);
			plot.setLimitsToFit(true);

			final RealVector peakAreas = doIntegrate(xvals, norms, means, sds);

			// Prepare columns for Results Table
			colLane.add(plot.getTitle());
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

			final String output = String.format(plot.getTitle() + ", RMS: %1$.2f; ",
				optimum.getRMS());
			log.info(output);
			statusServ.showProgress(++progress, plots.size());
		}

		updatePlots(); // without resetting or re-reading lanes
		final String[] headers = { "", "Band", "Distance", "Amplitude", "FWHM",
			"Area" };
		GenericColumn[] tableCol = new GenericColumn[headers.length];
		final DefaultGenericTable rt = new DefaultGenericTable();

		for (int cc = 0; cc < headers.length; cc++)
			tableCol[cc] = new GenericColumn(headers[cc]);
		tableCol[0].addAll(colLane);
		tableCol[1].addAll(colBand);
		for (int rr = 0; rr < colLane.size() - 1; rr++) {
			tableCol[2].add(String.format("%1$.1f", colDistance.getEntry(rr)));
			tableCol[3].add(String.format("%1$.1f", colAmplitude.getEntry(rr)));
			tableCol[4].add(String.format("%1$.2f", colFWHM.getEntry(rr)));
			tableCol[5].add(String.format("%1$.1f", colArea.getEntry(rr)));
		}
		for (int cc = 0; cc < headers.length; cc++)
			rt.add(tableCol[cc]);
		final DefaultTableDisplay tableDisplay = (DefaultTableDisplay) displayServ
			.createDisplay("Results Display", rt);
		displayServ.setActiveDisplay(tableDisplay);
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

	public void about() {
		Properties prop = new Properties();
		InputStream input = null;

		try {
			input = new FileInputStream("src//main//resources//about.properties");
			prop.load(input);
			System.out.println("Gauss Fit - v" + prop.getProperty("version"));
		}
		catch (IOException ex) {
			ex.printStackTrace();
		}
		finally {
			if (input != null) {
				try {
					input.close();
				}
				catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}

	public static void main(final String... args) throws Exception {
		// create the ImageJ application context with all available services
		ImageJ ij = net.imagej.Main.launch(args);
		ImagePlus iPlus = new Opener().openImage(
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

	@SuppressWarnings("serial")
	class MainDialog extends JFrame implements ActionListener, ChangeListener,
		DocumentListener, ItemListener, MouseMotionListener, MouseListener
	{

		private boolean auto = true; // AUTO Lane size mode is ON
		private boolean selectionUpdate = false; // active updating is off
		private String roiSelected = "none";
		private String plotSelected = "none";
		private String plotPreviouslySelected = "none";

		private int IW = imp.getWidth();
		private int IH = imp.getWidth();

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
		private JButton buttonCancel;
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

		public MainDialog(String string) {
			frame = new JFrame(string);
			setupMainDialog();
		}

		private void setupMainDialog() {
			imp.getCanvas().addMouseMotionListener(this);
			imp.getCanvas().addMouseListener(this);

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
			buttonCancel = new JButton("Cancel");
			buttonCancel.addActionListener(this);
			buttonPanel.add(buttonCancel);

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
					cleanup();
					frame.dispose();
				}
			});

			frame.getContentPane().add(dialogPanel);
			frame.setResizable(true);
			frame.validate();
			frame.pack();
			frame.setVisible(true);

			resetAutoROIs(LW, LH, LSp, LHOff, LVOff, nLanes);
			reDrawROIs(imp, "none");
			imp.killRoi();
		}

		public void cleanup() {
			// Remove Listeners on imp
			imp.getCanvas().removeMouseListener(this);
			imp.getCanvas().removeMouseMotionListener(this);
		}

		private boolean askUser(String question) {
			GenericDialog gd = new GenericDialog("WARNING!");
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
			TitledBorder tb = (TitledBorder) slider.getBorder();
			tb.setTitle(str);
			slider.setBorder(tb);
		}

		private void setSliderPanelEnabled(boolean enabled) {
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

		private void reDrawROIs(ImagePlus imgPlus, String roiName) {
			if (!auto) {
				// HouseKeeping: Sort the rois based on name and remove null elements
				final Comparator<Roi> nameComparator = new Comparator<Roi>() {

					@Override
					public int compare(Roi r1, Roi r2) {
						if (r1 == null && r2 == null) {
							return 0;
						}
						if (r1 == null) {
							return -1;
						}
						if (r2 == null) {
							return 1;
						}
						int n1 = Integer.parseInt(r1.getName().substring(5));
						int n2 = Integer.parseInt(r2.getName().substring(5));
						int comp = Integer.compare(n1, n2);
						if (comp != 0) {
							return comp;
						}
						return 0;
					}
				};
				Collections.sort(rois, nameComparator);
				Iterator<Roi> roisIter = rois.iterator();
				int nullIndex = 0;
				while (roisIter.hasNext()) {
					if (roisIter.next() == null) break;
					nullIndex++;
				}
				rois = new ArrayList<>(rois.subList(0, nullIndex));
			}
			Overlay overlay = new Overlay();
			Iterator<Roi> roiIter = rois.iterator();
			while (roiIter.hasNext()) {
				Roi roi = roiIter.next();
				String label = roi.getName().substring(5);
				Font font = new Font("SansSerif", Font.BOLD, 20);
				TextRoi labelRoi = new TextRoi(0, 0, label, font);
				labelRoi.setAntialiased(true);
				double lh = font.getSize() * 1.4;
				double lw = font.getSize() * 1.4;
				double rh = roi.getBounds().getHeight();
				double rw = roi.getBounds().getWidth();
				double x0 = roi.getBounds().getX();
				double y0 = roi.getBounds().getY();
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

		private void resetAutoROIs(int lw, int lh, int lsp, int lhoff, int lvoff,
			int nlanes)
		{
			rois = new ArrayList<>();
			for (int i = 0; i < nlanes; i++) {
				Roi roi = new Roi(lhoff + lw * i + lsp * i, lvoff, lw, lh);
				roi.setName("Lane " + (i + 1));
				rois.add(roi);
			}
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

		public double getTolPK() {
			try {
				tolPK = Double.parseDouble(textTolPK.getText().trim());
				return tolPK;
			}
			catch (final NumberFormatException e1) {
				return 0.0;
			}
		}

		public ArrayList<Roi> getRois() {
			return rois;
		}

		public Roi getRoi(String title) {
			Iterator<Roi> roiIter = rois.iterator();
			while (roiIter.hasNext()) {
				Roi roi = roiIter.next();
				if (roi.getName().equals(title)) return roi;
			}
			return null;
		}

		public String getPlotSelected() {
			return plotSelected;
		}

		public String getPlotPreviouslySelected() {
			return plotPreviouslySelected;
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
			plots = null; // need to recalculate the profiles
			updatePlots();
		}

		// TextField Listeners
		@Override
		public void changedUpdate(final DocumentEvent e) {
			nLanes = getNLanes();
			resetAutoROIs(LW, LH, LSp, LHOff, LVOff, nLanes);
			reDrawROIs(imp, "none");
			plots = null; // need to recalculate the profiles
			updatePlots();
		}

		@Override
		public void removeUpdate(final DocumentEvent e) {
			nLanes = getNLanes();
			resetAutoROIs(LW, LH, LSp, LHOff, LVOff, nLanes);
			reDrawROIs(imp, "none");
			plots = null; // need to recalculate the profiles
			updatePlots();
		}

		@Override
		public void insertUpdate(final DocumentEvent e) {
			nLanes = getNLanes();
			resetAutoROIs(LW, LH, LSp, LHOff, LVOff, nLanes);
			reDrawROIs(imp, "none");
			plots = null; // need to recalculate the profiles
			updatePlots();
		}

		// Buttons
		@Override
		public void actionPerformed(final ActionEvent e) {
			if (e.getSource().equals(buttonMeasure)) {
				doFit();
				chkBoxBands.setEnabled(true);
			}

			if (e.getSource().equals(buttonCancel)) {
				cleanup();
				plotImage.getWindow().close();
				List<Display<?>> displays = displayServ.getDisplays();
				for (Display<?> d : displays)
					d.close();
				this.dispose();
				log.info("Gauss Fit terminated.");
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
				// TODO Auto-generated method stub
			}
		}

		@Override
		public void mouseDragged(MouseEvent e) {

		}

		@Override
		public void mouseMoved(MouseEvent e) {
			statusServ.showStatus("[" + e.getX() + ":" + e.getY() + "]");
			if (!selectionUpdate) {
				String roiCurrent = "none"; // None selected
				// Local copy for iteration
				Iterator<Roi> roisIter = rois.iterator();
				while (roisIter.hasNext()) {
					Roi roi = roisIter.next();
					int x = ((ImageCanvas) e.getSource()).offScreenX(e.getX());
					int y = ((ImageCanvas) e.getSource()).offScreenX(e.getY());
					if (roi.contains(x, y)) roiCurrent = roi.getName();
				}
				if (roiCurrent.equals("none") && imp.getRoi() != null) imp.killRoi();
				if (!roiCurrent.equals(roiSelected)) {
					roiSelected = roiCurrent;
					plotPreviouslySelected = plotSelected;
					plotSelected = roiCurrent;
					if (buttonAuto.isSelected()) resetAutoROIs(LW, LH, LSp, LHOff, LVOff,
						nLanes);
					reDrawROIs(imp, roiCurrent);
					// This does not change the profiles, so do not set plots=null
					updatePlots();
				}
			}
		}

		@Override
		public void mouseClicked(MouseEvent e) {
			if (e.getClickCount() == 2) {
				System.out.println("M 2-clicked");
				if (!auto && !roiSelected.equals("none") && !selectionUpdate) {
					// Remove Roi
					Roi roi = md.getRoi(roiSelected);
					if (roi == null) return;
					if (roi.isHandle(e.getX(), e.getY()) != -1) {
						// might be dragging existing roi
						return;
					}

					if (askUser("Delete " + roiSelected + "?")) {
						int x = ((ImageCanvas) e.getSource()).offScreenX(e.getX());
						int y = ((ImageCanvas) e.getSource()).offScreenX(e.getY());
						Iterator<Roi> roiIter = rois.iterator();
						while (roiIter.hasNext()) {
							roi = roiIter.next();
							if (roi.contains(x, y)) {
								roiIter.remove();
								System.out.println(imp.getTitle() + ": sel(" + roiSelected +
									"); REMOVED");
							}
						}
						plots = null;
						reDrawROIs(imp, "none");
						updatePlots();
					}
				}
			}
		}

		@Override
		public void mousePressed(MouseEvent e) {
			System.out.println("M pressed");
			if (!auto) selectionUpdate = true;
		}

		@Override
		public void mouseReleased(MouseEvent e) {
			// Add or modify Roi
			if (buttonManual.isSelected()) {
				if (selectionUpdate) {
					Roi roiNew = imp.getRoi();
					if (roiNew != null && roiNew.getBounds().getWidth() * roiNew
						.getBounds().getHeight() > 100)
					{

						if (roiSelected.equals("none")) {
							System.out.println("Add");
							// Creating New ROI
							// Find an available integer for lane name
							int i = 0;
							while (++i <= rois.size() + 1) {
								Iterator<Roi> roisIter = rois.iterator();
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
						}
						else {
							System.out.println("Modify");
							// Moving/Resizing a specific ROI
							Iterator<Roi> roisIter = rois.iterator();
							while (roisIter.hasNext()) {
								Roi roi = roisIter.next();
								if (roi != null) {
									if (roi.getName().equals(roiSelected)) {
										roiNew.setName(roi.getName());
										rois.set(rois.indexOf(roi), roiNew);
										break;
									}
								}
							}
						}
						plots = null;
						reDrawROIs(imp, roiSelected);
						updatePlots();
					}
				}
			}
			selectionUpdate = false;
		}

		@Override
		public void mouseEntered(MouseEvent e) {
			// TODO Auto-generated method stub

		}

		@Override
		public void mouseExited(MouseEvent e) {
			// TODO Auto-generated method stub

		}
	}
}
