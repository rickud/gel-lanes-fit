
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
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

import net.imagej.Dataset;
import net.imagej.DatasetService;
import net.imagej.ImageJ;
import net.imagej.ImgPlus;
import net.imagej.display.ImageDisplay;
import net.imagej.display.ImageDisplayService;
import net.imagej.display.OverlayService;
import net.imagej.ops.OpService;
import net.imagej.overlay.RectangleOverlay;
import net.imagej.table.DefaultGenericTable;
import net.imagej.table.DefaultTableDisplay;
import net.imagej.table.GenericColumn;
import net.imglib2.type.numeric.RealType;

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
import org.scijava.display.DisplayService;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.util.Colors;

import gausscurvefit.GaussianArrayCurveFitter.ParameterGuesser;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;
import ij.gui.Line;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.ProfilePlot;
import ij.gui.Roi;
import ij.io.Opener;
import ij.measure.Calibration;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;
import ij.util.Tools;

@Plugin(type = Command.class, headless = true,
	menuPath = "Plugins>Gel Tools>Gauss Fit")
public class GaussFitLanes implements Command {

	@Parameter
	private LogService log;

	@Parameter
	private StatusService statusServ;

	@Parameter
	private static DatasetService datasetServ;

	@Parameter
	private static DatasetIOService datasetIOServ;

	@Parameter
	private ConvertService convertServ;

	@Parameter
	private ImageDisplayService imgDisplayServ;

	@Parameter
	private OverlayService overlayServ;

	@Parameter
	private static DisplayService displayServ;

	@Parameter
	private static OpService ops;

	// Default Parameters
	private Thread mainThread; // thread for plotting (in the background)
	private Thread plotThread; // thread for plotting (in the background)
	private boolean setup = true;
	private boolean doFit; // tells the background thread to update
	private final ImagePlus plotImage = new ImagePlus();
	private ArrayList<Roi> rois = new ArrayList<>();
	private Plot[] plots;
	private static ImgPlus<? extends RealType<?>> img;
	private ImageDisplay imageDisplay;
	private ImagePlus imp;
	private CustomDialog cd;

	private final int rows = 2; // Number of plot Rows in display
	private final int cols = 2; // Number of plot Rows in display
	private int nLanes = 4;

	private int degBG = 3; // Order of Background Polynomial
	private double tolPK = 0.05; // Peak detection tolerance as % of range

	public void init() {

		imp = IJ.getImage();

		try {
			Dataset dataset = datasetIOServ.open(
				"src//main//resources//sample//All[01-17-2017].tif");
			displayServ.createDisplay(dataset.getName(), dataset);
			imageDisplay = (ImageDisplay) displayServ.getDisplay(dataset.getName());
			img = dataset.getImgPlus();
		}
		catch (IOException exc) {
			// TODO Auto-generated catch block
			exc.printStackTrace();
		}

		cd = new CustomDialog("Gel Lanes Gauss Fitting:" + imp.getTitle());
		IJ.wait(50); // delay to make sure ROIs have updated
		plotImage.setTitle("Profiles of " + imp.getShortTitle());
		updatePlots();
		plotImage.show();
		if (plotImage.getRoi() != null) plotImage.deleteRoi();
		imp.getCanvas().requestFocus();
		final ImageWindow iwin = imp.getWindow();
		final ImageWindow pwin = plotImage.getWindow();
		if (iwin == null || pwin == null) return;
		final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
		final Dimension imageSize = iwin.getSize();
		final Dimension dialogSize = pwin.getSize();
		final Point imageLoc = iwin.getLocation();
		int x = imageLoc.x + imageSize.width + 30;
		if (x + dialogSize.width > screen.width) x = screen.width -
			dialogSize.width;
		pwin.setLocation(x, imageLoc.y);
		final ImageCanvas canvas = iwin.getCanvas();
		canvas.requestFocus();
		pwin.setVisible(true);

		// thread for plotting in the background
		plotThread = new Thread(this, "Dynamic Plots");
		plotThread.setPriority(Math.max(plotThread.getPriority() - 3,
			Thread.MIN_PRIORITY));
		plotThread.start();
	}

	// the background thread for plotting.
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
	 * @param laneRoi
	 */
	private RealMatrix getLaneProfile(final int laneRoi) {
		if (rois.size() == 0) return null;

		final Roi roi = rois.get(laneRoi);
		imp.setRoi(roi);
		final ProfilePlot profileP = new ProfilePlot(imp, true); // get the profile
		final RealVector profile = new ArrayRealVector(profileP.getProfile());
		if (profile.getDimension() < 2) return null;

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
		if (rois.size() == 0) return;

		final ImageProcessor ip = imp.getProcessor();
		if (plots == null) { // Plots have been reset
			plots = new Plot[rois.size()];
			for (int p = 0; p < rois.size(); p++) {
				final Roi roi = rois.get(p);;
				if (ip == null || roi == null) return; // these may change
																								// asynchronously
				final RealVector profile = new ArrayRealVector(getLaneProfile(p).getRow(
					1));
				if (profile.getDimension() < 2) return;

				final Calibration cal = imp.getCalibration();
				String xUnit;
				if (cal.getUnit() == null) xUnit = "pixels";
				else xUnit = cal.getUnit();
				final RealVector x = new ArrayRealVector(getLaneProfile(p).getRow(0));
				final String xLabel = "Distance (" + xUnit + ")";
				final String yLabel = (cal.getValueUnit() != null && !cal.getValueUnit()
					.equals("Gray Value")) ? "Value (" + cal.getValueUnit() + ")"
						: "Value";
				plots[p] = new Plot("Lane " + (p + 1), xLabel, yLabel, x.toArray(),
					profile.toArray());

				final double fixedMin = ProfilePlot.getFixedMin();
				final double fixedMax = ProfilePlot.getFixedMax();
				if (fixedMin != 0 || fixedMax != 0) {
					final double[] a = Tools.getMinMax(x.toArray());
					plots[p].setLimits(a[0], a[1], fixedMin, fixedMax);
				}
			}
		}
		// Construct the plot image
		int pages = Math.floorDiv(plots.length, rows * cols);
		if (Math.floorMod(plots.length, rows * cols) != 0) pages++;

		final int plotSpacing = 5; // black border
		final int plotW = plots[0].getProcessor().getWidth();
		final int plotH = plots[0].getProcessor().getHeight();
		final int plotsWidth = plotSpacing + cols * (plotW + plotSpacing);
		final int plotsHeight = plotSpacing + rows * (plotH + plotSpacing);
		final ImageProcessor blank = plots[0].getProcessor().duplicate();
		IJ.setForegroundColor(255, 255, 255);
		blank.fill();

		final ImageProcessor[] pageMontages = new ImageProcessor[pages];
		for (int pg = 0; pg < pages; pg++) {
			final ImageProcessor plotsMontage = new ColorProcessor(plotsWidth,
				plotsHeight);
			for (int c = 0; c < cols; c++) {
				for (int r = 0; r < rows; r++) {
					if (c * rows + r < plots.length) {
						plotsMontage.insert(plots[pg * cols * rows + (c * rows + r)]
							.getProcessor(), // here NP
							plotSpacing + c * (plotSpacing + plotW), plotSpacing + r *
								(plotSpacing + plotH));
						plotsMontage.drawString(plots[pg * cols * rows + (c * rows + r)]
							.getTitle(), plotSpacing + plotW / 2 + c * (plotSpacing + plotW),
							plotSpacing + plotH / 5 + r * (plotSpacing + plotH),
							Color.LIGHT_GRAY);
					}
					else {
						plotsMontage.insert(blank, plotSpacing + c * (plotSpacing + plotW),
							plotSpacing + r * (plotSpacing + plotH));
					}
				}
			}
			pageMontages[pg] = plotsMontage;
		}
		plotImage.setProcessor(pageMontages[0]);
		final ImageStack plotStack = plotImage.createEmptyStack();

		for (int i = 0; i < pageMontages.length; i++) {
			if (pageMontages[i] != null) plotStack.addSlice("Lanes " + (rows * cols *
				i + 1) + "-" + (rows * cols * i + rows * cols), pageMontages[i], i);
		}
		plotImage.setStack(plotStack);
		plotImage.updateImage();
		if (plotImage.getRoi() != null) plotImage.deleteRoi();
	}

	private void doFit() {
		// Reset the plots window
		if (plots == null) updatePlots();

		// Results Table Colums
		final ArrayList<String> colLane = new ArrayList<>();
		final ArrayList<Integer> colBand = new ArrayList<>();
		RealVector colDistance = new ArrayRealVector();
		RealVector colAmplitude = new ArrayRealVector();
		RealVector colFWHM = new ArrayRealVector();
		RealVector colArea = new ArrayRealVector();

		for (int i = 0; i < plots.length; i++) {
			if (log.getLevel() == 0) log.initialize();

			final RealVector xvals = new ArrayRealVector(getLaneProfile(i).getRow(0));
			final RealVector yvals = new ArrayRealVector(getLaneProfile(i).getRow(1));
			RealVector bg = new ArrayRealVector();

			final int degbg = cd.getDegBG();
			// Tolerance as percentage of the range
			final double tolpk = cd.getTolPK() * (yvals.getMaxValue() - yvals
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
			final String output = String.format("Lane %1$d, RMS: %2$.2f; ", i + 1,
				optimum.getRMS());
			log.info(output);
			bg = xvals.map(new PolynomialFunction(poly.getSubVector(1, degBG + 1)
				.toArray()));
			plots[i].setColor(Color.blue);
			plots[i].addPoints(xvals.toArray(), bg.toArray(), PlotWindow.LINE);

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
				plots[i].setColor(Color.red);
				plots[i].addPoints(xvals.toArray(), gauss.toArray(), PlotWindow.LINE);
			}
			fitted = xvals.map(new GaussianArrayBG(norms, means, sds, poly));
			guessed = xvals.map(new GaussianArrayBG(norms0, means0, sds0, poly0));

			plots[i].setColor(Color.blue);
			plots[i].addPoints(means0.toArray(), norms0.toArray(), PlotWindow.CROSS);
			plots[i].setColor(new Color(0, 128, 0));
			plots[i].addPoints(xvals.toArray(), fitted.toArray(), PlotWindow.LINE);
			plots[i].setLimitsToFit(true);

			final RealVector peakAreas = doIntegrate(xvals, norms, means, sds);

			// Prepare columns for Results Table
			colLane.add("Lane " + (i + 1));
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

	public static void main(final String... args) throws Exception {
		// create the ImageJ application context with all available services
		ImageJ ij = net.imagej.Main.launch(args);

		final ImagePlus imp = new Opener().openImage(
			"src//main//resources//sample//All[01-17-2017].tif");
		// display it via ImageJ
		imp.show();
		// wrap it into an ImgLib image (no copying)
		// final Img image = ImagePlusAdapter.wrap(imp);
		// display it via ImgLib using ImageJ
		// ImageJFunctions.show(image);

	}

	@SuppressWarnings("serial")
	class CustomDialog extends JFrame implements ActionListener, ChangeListener,
		DocumentListener, ItemListener
	{

		int IW = imp.getWidth();
		int IH = imp.getWidth();

		// Default lane size/offset (Just center 4 lanes in the image)
		int LW = (int) Math.round(0.9 * IW / nLanes);
		int LH = (int) Math.round(IH * 0.9);
		int LSp = Math.round((IW - LW * nLanes) / (nLanes + 1));
		int LHOff = LSp;
		int LVOff = (IH - LH) / 2;

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

		public CustomDialog(String string) {
			frame = new JFrame(string);
			showMainDialog();
		}

		private void showMainDialog() {
			textPanel = new JPanel();
			textPanel.setLayout(new FlowLayout(FlowLayout.LEADING));

			textNLanes = new JTextField(3);
			textNLanes.setText("" + nLanes);
			textNLanes.setVisible(true);
			textNLanes.getDocument().addDocumentListener(this);
			textPanel.add(textNLanes);

			labelNLanes = new JLabel("Number of Lanes");
			textPanel.add(labelNLanes);

			sliderPanel = new JPanel();
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
			dialogPanel.setBackground(Color.lightGray);
			dialogPanel.setLayout(new BorderLayout());
			dialogPanel.add(textPanel, BorderLayout.NORTH);
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

			final ImageWindow iwin = imp.getWindow();
			if (iwin == null) return;
			final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
			final Dimension imageSize = iwin.getSize();
			final Dimension dialogSize = frame.getSize();
			final Point imageLoc = iwin.getLocation();
			int x = imageLoc.x + imageSize.width + 10;
			if (x + dialogSize.width > screen.width) x = screen.width -
				dialogSize.width;
			frame.setLocation(x, imageLoc.y);
			final ImageCanvas canvas = iwin.getCanvas();
			canvas.requestFocus();
			frame.setVisible(true);
			reDrawROIs(LW, LH, LSp, LHOff, LVOff);
		}

		public void cleanup() {
			// Remove Listeners on imp
			// imp.removeImageListener(listener);

		}

		private JSlider makeTitledSlider(final String string, final Color color,
			final int minVal, final int maxVal, final int val)
		{
			// Border empty = BorderFactory.createTitledBorder(
			// BorderFactory.createEmptyBorder() );

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

		private void setSliderTitle(final JSlider slider, final Color color,
			final String str)
		{
			// Border empty = BorderFactory.createTitledBorder(
			// BorderFactory.createEmptyBorder() );
			final TitledBorder tb = new TitledBorder(BorderFactory
				.createEtchedBorder(), // empty,
				"", TitledBorder.CENTER, TitledBorder.BELOW_TOP, new Font("Sans",
					Font.PLAIN, 11));
			tb.setTitleJustification(TitledBorder.LEFT);
			tb.setTitle(str);
			tb.setTitleColor(color);
			slider.setBorder(tb);
		}

		@Override
		public synchronized void stateChanged(final ChangeEvent e) {
			final JSlider slider = (JSlider) e.getSource();
			if (slider == sliderW) {
				LW = sliderW.getValue();
				final String str = "Width ( " + LW + " px )";
				setSliderTitle(sliderW, Color.black, str);
			}
			else if (slider == sliderH) {
				LH = sliderH.getValue();
				final String str = "Height ( " + LH + " px )";
				setSliderTitle(sliderH, Color.black, str);
			}
			else if (slider == sliderSp) {
				LSp = sliderSp.getValue();
				final String str = "Spacing ( " + LSp + " px )";
				setSliderTitle(sliderSp, Color.black, str);
			}
			else if (slider == sliderHOff) {
				LHOff = sliderHOff.getValue();
				final String str = "Horizontal Offset ( " + LHOff + " px )";
				setSliderTitle(sliderHOff, Color.black, str);
			}
			else if (slider == sliderVOff) {
				LVOff = sliderVOff.getValue();
				final String str = "Vertical Offset ( " + LVOff + " px )";
				setSliderTitle(sliderVOff, Color.black, str);
			}
			reDrawROIs(LW, LH, LSp, LHOff, LVOff);
			// IJ.wait(50);
			// delay to make sure the roi has been updated
		}

		private void reDrawROIs(int lw, int lh, int lsp, int lhoff, int lvoff) {
//			if (!WindowManager.getCurrentImage().equals(image)) {
//				displayServ.setActiveDisplay(imageDisplay);
//			}
//			if (imageDisplay.getInfo() != null || rois.size() != 0) {
//				rois = null;
//				imp.deleteRoi();
//			}
			nLanes = getNLanes();
			ArrayList<RectangleOverlay> overlays = new ArrayList<>();

			for (int i = 0; i < nLanes; i++) {
				Roi roi = new Roi(lhoff + lw * i + lsp * i, lvoff, lw, lh);
				rois.add(roi);
				RectangleOverlay ol = new RectangleOverlay(imageDisplay.getContext());
				ol.setOrigin(lhoff + lw * i + lsp * i, 0);
				ol.setOrigin(lvoff, 1);
				ol.setExtent(lw, 0);
				ol.setExtent(lh, 1);
				ol.setLineColor(Colors.YELLOW);
				ol.setLineWidth(2);
				ol.setFillColor(Colors.YELLOW);
				ol.setAlpha(128);
				overlays.add(ol);
			}
			overlayServ.addOverlays(imageDisplay, overlays);
			plots = null; // Reset the profile plots
			updatePlots();
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

//		boolean isRoi() {
//			if (imp == null) return false;
//			final Roi roi = imp.getRoi();
//			if (roi == null) return false;
//			return roi.getType() == Roi.LINE || roi.getType() == Roi.RECTANGLE;
//		}

		// TextField Listeners
		@Override
		public void changedUpdate(final DocumentEvent e) {
			reDrawROIs(LW, LH, LSp, LHOff, LVOff);
		}

		@Override
		public void removeUpdate(final DocumentEvent e) {
			reDrawROIs(LW, LH, LSp, LHOff, LVOff);
		}

		@Override
		public void insertUpdate(final DocumentEvent e) {
			reDrawROIs(LW, LH, LSp, LHOff, LVOff);
		}

		// Measure Button
		@Override
		public void actionPerformed(final ActionEvent e) {
			if (e.getSource().equals(buttonMeasure)) {
				doFit();
				chkBoxBands.setEnabled(true);
			}
			if (e.getSource().equals(buttonCancel)) {
				cleanup();
				plotImage.getWindow().close();
				this.dispose();
				log.info("Gauss Fit terminated.");
			}
		}

		@Override
		public void itemStateChanged(final ItemEvent e) {
			if (e.getItemSelectable() == chkBoxBands) {
				// TODO Auto-generated method stub
			}
		}
	}
}
