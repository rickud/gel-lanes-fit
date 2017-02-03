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

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

import org.apache.commons.collections4.CollectionUtils;
import org.apache.commons.collections4.IterableMap;
import org.apache.commons.collections4.MapIterator;
import org.apache.commons.collections4.collection.IndexedCollection;
import org.apache.commons.collections4.iterators.ArrayListIterator;
import org.apache.commons.collections4.list.FixedSizeList;
import org.apache.commons.collections4.map.HashedMap;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.math.NumberUtils;
import org.apache.commons.math4.analysis.function.Gaussian;
import org.apache.commons.math4.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math4.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math4.exception.ConvergenceException;
import org.apache.commons.math4.exception.DimensionMismatchException;
import org.apache.commons.math4.fitting.WeightedObservedPoints;
import org.apache.commons.math4.fitting.leastsquares.LeastSquaresOptimizer;
import org.apache.commons.math4.fitting.leastsquares.LeastSquaresProblem;
import org.apache.commons.math4.fitting.leastsquares.LevenbergMarquardtOptimizer;
import org.apache.commons.math4.linear.ArrayRealVector;
import org.apache.commons.math4.linear.MatrixUtils;
import org.apache.commons.math4.linear.RealMatrix;
import org.apache.commons.math4.linear.RealVector;
import org.apache.commons.math4.util.OpenIntToDoubleHashMap.Iterator;
import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.command.Previewable;
import org.scijava.convert.ConvertService;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;
import ij.gui.Line;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.ProfilePlot;
import ij.gui.Roi;
import ij.io.Opener;
import ij.measure.Calibration;
import ij.plugin.frame.RoiManager;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;
import ij.util.Tools;
import io.scif.services.DatasetIOService;
import net.imagej.DatasetService;
import net.imagej.ImageJ;
import net.imagej.display.ImageDisplayService;
import net.imagej.ops.OpService;
import scala.collection.mutable.HashMap;


@Plugin(type = Command.class, headless = true,
menuPath = "Plugins>Gauss Fit")
public class GaussFitLanes implements 
Command, Previewable, Runnable {

	@Parameter
	private LogService log;

	@Parameter
	private StatusService statusServ;

	@Parameter
	private DatasetService datasetServ;

	@Parameter
	private DatasetIOService datasetIOServ;

	@Parameter
	private ConvertService convertServ;

	@Parameter
	private ImageDisplayService imgServ;

	@Parameter
	private static OpService ops;

	// Default Parameters
	private Thread mainThread;           //thread for plotting (in the background)
	private Thread plotThread;           //thread for plotting (in the background)
	private boolean setup = true;
	private boolean doFit;               //tells the background thread to update
	ImagePlus imp;
	private Plot[] plots; 
	private boolean subMin = true;  // subtract min y value from all profile plots

	private ImagePlus plotImage = new ImagePlus();
	private RoiManager roiMan = new RoiManager();
	private CustomDialog cd;
	private int rows = 2; // Number of plot Rows in display
	private int cols = 2; // Number of plot Rows in display

	public void init() {
		imp = IJ.getImage();
		roiMan = RoiManager.getInstance(); 
		if(roiMan==null) { 
			roiMan = new RoiManager(true); 
		} 
		roiMan.runCommand("Show All");
		roiMan.runCommand("Labels");

		cd = new CustomDialog();		
		cd.showMainDialog();

		IJ.wait(50); //delay to make sure ROIs have updated
		plotImage.setTitle("Profiles of " + imp.getShortTitle()); 

		updatePlots(rows, cols, subMin);
		
		plotImage.show();
		if (plotImage.getRoi()!=null) plotImage.deleteRoi();
		imp.getCanvas().requestFocus();

		IJ.wait(50);
		ImageWindow iwin = imp.getWindow();
		ImageWindow pwin = plotImage.getWindow();

		if (iwin == null || pwin==null) return;
		Dimension screen     = Toolkit.getDefaultToolkit().getScreenSize();
		Dimension imageSize  = iwin.getSize();
		Dimension dialogSize = pwin.getSize();
		Point imageLoc = iwin.getLocation();
		int x = imageLoc.x+imageSize.width+30;
		if (x+dialogSize.width>screen.width)
			x = screen.width-dialogSize.width;
		pwin.setLocation(x, imageLoc.y);
		ImageCanvas canvas = iwin.getCanvas();
		canvas.requestFocus();
		pwin.setVisible(true);

		// thread for plotting in the background
		plotThread = new Thread(this, "Dynamic Plots");
		plotThread.setPriority(Math.max(plotThread.getPriority()-3, Thread.MIN_PRIORITY));
		plotThread.start();
	}

	// the background thread for plotting.
	public void run() {
		if (setup) { init(); setup = false;}		

		synchronized(this) {
			//				if (doUpdate) {
			//					doUpdate = false;               //and loop again
			//				} else {
			//					try {wait();}                   //notify wakes up the thread
			//					catch(InterruptedException e) { //interrupted tells the thread to exit
			//						return;
			//					}
			//				}
		};
	}

	/** Profile data from Roi, ready for fitting 
	 *  (null if not possible) */
	private RealMatrix getLaneProfile(int laneRoi, boolean minValueZero) {
		if (roiMan.getCount() == 0) return null;
		roiMan.select(laneRoi);
		Roi roi = roiMan.getRoi(laneRoi);	
		ProfilePlot profileP = new ProfilePlot(imp, true);//get the profile
		RealVector profile = new ArrayRealVector(profileP.getProfile());
		if (profile==null || profile.getDimension()<2)
			return null;
		if (minValueZero)
			profile.mapSubtractToSelf(profile.getMinValue());

		//the following code is mainly for x calibration                   
		Calibration cal = imp.getCalibration();
		RealVector x = calibrateX(profile,roi,cal);
		RealMatrix output = MatrixUtils.createRealMatrix(
				new double[][]{x.toArray(),profile.toArray()});
		return output;
	}
	
	/**
	 * Method for x calibration
	 **/
	private RealVector calibrateX(RealVector y, Roi roi, Calibration cal) {                    
		double xInc = 1;
		if (roi.getType() == Roi.LINE) {
			Line line = (Line)roi;
			if (cal != null) {
				double dx = cal.pixelWidth*(line.x2 - line.x1);
				double dy = cal.pixelHeight*(line.y2 - line.y1);
				double length = Math.sqrt(dx*dx + dy*dy);
				xInc = length/(y.getMaxIndex());
			}
		} else if (roi.getType() == Roi.RECTANGLE) {
			if (cal != null) {
				xInc = roi.getBounds().getHeight()*cal.pixelHeight/(y.getMaxIndex());
			}
		} else return null;
		

		double[] x = new double[y.getDimension()];// create the x axis
		for (int i=0; i<y.getDimension(); i++)
			x[i] = i*xInc;
		return new ArrayRealVector(x);
	}
	
	/**
	 * Method for Output plot collage
	 * */
	private void updatePlots(int rows, int cols, boolean minValueZero){
		if (roiMan.getCount() == 0) return;
		ImageProcessor ip = imp.getProcessor();
		if (plots == null){
			plots = new Plot[roiMan.getCount()]; 
			for (int p = 0; p<roiMan.getCount(); p++) {
				roiMan.select(p);
				Roi roi = roiMan.getRoi(p);
				if (ip == null || roi == null) return; //these may change asynchronously
				if (roi.getType() == Roi.LINE)
					ip.setInterpolate(PlotWindow.interpolate);
				else
					ip.setInterpolate(false);
				ProfilePlot profileP = new ProfilePlot(imp, true);//get the profile
				RealVector profile = new ArrayRealVector(profileP.getProfile());
				if (minValueZero) {
					profile.mapSubtractToSelf(profile.getMinValue());
				}
				if (profile==null || profile.getDimension()<2)
					return;

				Calibration cal = imp.getCalibration();
				String xUnit;
				if (cal.getUnit() == null)
					xUnit = "pixels";
				else
					xUnit = cal.getUnit();
				RealVector x = calibrateX(profile, roi, cal);
				String xLabel = "Distance (" + xUnit + ")";
				String yLabel = (cal !=null && cal.getValueUnit()!=null && !cal.getValueUnit().equals("Gray Value")) ?
						"Value ("+cal.getValueUnit()+")" : "Value";
				plots[p] = new Plot("Lane "+(p+1), 
						xLabel, yLabel, x.toArray(), profile.toArray());
				//plots[p].addText(plots[p].getTitle(), plots[p].getSize().getWidth()/5, plots[p].getSize().getHeight()/2);
				double fixedMin = ProfilePlot.getFixedMin();
				double fixedMax = ProfilePlot.getFixedMax();
				if (fixedMin!=0 || fixedMax!=0) {
					double[] a = Tools.getMinMax(x.toArray());
					plots[p].setLimits(a[0],a[1], fixedMin, fixedMax);
				}
			}
		}
		// Construct the plot image
		int pages = Math.floorDiv(plots.length, rows*cols);
		if (Math.floorMod(plots.length, rows*cols)!=0) pages++;

		int plotSpacing = 5; // black border
		int plotW = plots[0].getProcessor().getWidth(); int plotH = plots[0].getProcessor().getHeight();
		int plotsWidth  = plotSpacing + cols*( plotW + plotSpacing );
		int plotsHeight = plotSpacing + rows*( plotH + plotSpacing );
		ImageProcessor blank = plots[0].getProcessor().duplicate();
		IJ.setForegroundColor(255,255,255);
		blank.fill();

		ImageProcessor[] pageMontages = new ImageProcessor[pages];
		for (int pg = 0; pg<pages; pg++){
			ImageProcessor plotsMontage = new ColorProcessor(plotsWidth, plotsHeight);
			for (int c = 0; c<cols; c++) {
				for (int r = 0; r<rows; r++) {
					if (c*rows+r<plots.length){
						plotsMontage.insert(plots[pg*cols*rows+(c*rows+r)].getProcessor(), //here NP
								plotSpacing+c*(plotSpacing + plotW), 
								plotSpacing+r*(plotSpacing + plotH));
						plotsMontage.drawString(plots[pg*cols*rows+(c*rows+r)].getTitle(),
								plotSpacing+plotW/2+c*(plotSpacing + plotW),
								plotSpacing+plotH/5+r*(plotSpacing + plotH),
								Color.LIGHT_GRAY);
					} else {
						plotsMontage.insert(blank, 
								plotSpacing+c*(plotSpacing + plotW), 
								plotSpacing+r*(plotSpacing + plotH));
					}
				}
			}
			pageMontages[pg] = plotsMontage;
		}
		plotImage.setProcessor(pageMontages[0]);
		ImageStack plotStack = plotImage.createEmptyStack();

		for (int i=0; i<pageMontages.length; i++){
			if (pageMontages[i] != null)
				plotStack.addSlice("Lanes "+ (rows*cols*i+1) + "-" + (rows*cols*i+rows*cols),
						pageMontages[i], i);
		}
		plotImage.setStack(plotStack);
		plotImage.updateImage();
		if (plotImage.getRoi()!=null) plotImage.deleteRoi();
	}

	private void doFit() throws IOException {
		//Reset the plots window
		for (int i = 0; i<plots.length; i++) {
			RealVector xvals = new ArrayRealVector(getLaneProfile(i,subMin).getRow(0));
			RealVector yvals = new ArrayRealVector(getLaneProfile(i,subMin).getRow(1));
			RealVector bg    = new ArrayRealVector(yvals.getDimension(), yvals.getMinValue()); 
			
			WeightedObservedPoints obs   = new WeightedObservedPoints(); // All Y
			BufferedWriter br = new BufferedWriter(new FileWriter("output//Plot"+(i+1)+".csv"));
			StringBuilder  sb = new StringBuilder();
			yvals = new ArrayRealVector(yvals).mapSubtractToSelf(yvals.getMinValue());
			for (int o = 0 ; o<xvals.getDimension(); o++) {
				obs.add(xvals.getEntry(o), yvals.getEntry(o));
				sb.append(xvals.getEntry(o) +"\t"+yvals.getEntry(o)+"\n");
			}
			br.write(sb.toString());
			br.close();
			
//			plots[i].setColor(Color.black);
//			plots[i].addPoints(xvals.toArray(), yvals.toArray(), PlotWindow.LINE);

			// Tolerances as percentage of the range
			double tolbg = cd.getTolBG()*(yvals.getMaxValue()-yvals.getMinValue());
			double tolpk = cd.getTolPK()*(yvals.getMaxValue()-yvals.getMinValue());

			

			//double[] pars = GaussianCurveFitter.create().withStartPoint(ArrayUtils.subarray(guessStart, 0, 3)).fit(obs.toList());
			LeastSquaresProblem problem  = GaussianArrayCurveFitter.
					create(tolpk).
					getProblem(obs.toList());

			LeastSquaresOptimizer.Optimum optimum;
			RealVector pars = new ArrayRealVector();
			
				optimum = new LevenbergMarquardtOptimizer()
						//.withCostRelativeTolerance(1e-100)
						//.withOrthoTolerance(1e-100)
						//.withParameterRelativeTolerance(1e-15)
						//.withInitialStepBoundFactor(0.1)
						.optimize(problem);
				System.out.println("RMS: "           + optimum.getRMS());
				System.out.println("evaluations: "   + optimum.getEvaluations());
				System.out.println("iterations: "    + optimum.getIterations());


				pars  = new ArrayRealVector(
						GaussianArrayCurveFitter.
						create(tolpk).
						fit(obs.toList()));
//			} catch (ConvergenceException e) {
//				System.out.println("No Convergence: "+ e.getCause());
//				pars= new ArrayRealVector();
//			}
			RealVector diff = new ArrayRealVector(pars).subtract(optimum.getPoint());
			RealVector gauss  = new ArrayRealVector();
			RealVector gauss0 = new ArrayRealVector();
			for (int b=0; b<pars.getDimension()/3; b++){ // plot one gaussian for each peak
				gauss = xvals.map(
						new Gaussian(pars.getEntry(b*3), pars.getEntry(b*3+1), pars.getEntry(b*3+2)));
	
				plots[i].setColor(Color.red);
				plots[i].addPoints(xvals.toArray(), gauss.toArray(), PlotWindow.LINE);
				gauss0 = xvals.map(
						new Gaussian(norms[b],means[b],stds[b]));					
				plots[i].setColor(Color.blue);
				plots[i].addPoints(xvals.toArray(), gauss0.toArray(), PlotWindow.LINE);
			}
			plots[i].setColor(Color.blue);
			plots[i].addPoints(means, norms, PlotWindow.CIRCLE);


			String parStrN = String.format("Lane %1$d-parN: ", i+1);
			String parStrM = String.format("Lane %1$d-parM: ", i+1);
			String parStrS = String.format("Lane %1$d-parS: ", i+1);
			for (int pp = 0; pp<pars.getDimension(); pp+=3){
				parStrN += String.format("%1$.2f, ", pars.getEntry(pp  ));
				parStrM += String.format("%1$.2f, ", pars.getEntry(pp+1));
				parStrS += String.format("%1$.2f, ", pars.getEntry(pp+2));
			}

			System.out.println(parStrN.substring(0,parStrN.length()-2));
			System.out.println(parStrM.substring(0,parStrM.length()-2));
			System.out.println(parStrS.substring(0,parStrS.length()-2));

			plots[i].setColor(Color.blue);
			plots[i].addPoints(xvals.toArray(), bg.toArray(), PlotWindow.LINE);
			plots[i].setLimitsToFit(true);
		}	
	}

	
	public static void main(final String... args) throws Exception {
		// create the ImageJ application context with all available services
		final ImageJ ij = net.imagej.Main.launch(args);
		ImagePlus imp = new Opener().openImage("src//main//resources//sample//All[01-17-2017].tif");
		// display it via ImageJ
		imp.show();
		// wrap it into an ImgLib image (no copying)
		// final Img image = ImagePlusAdapter.wrap(imp);
		// display it via ImgLib using ImageJ
		// ImageJFunctions.show( image );

	}

	public void cancel() {
		log.info("Process canceled");
		return;
	}

	public void preview() {
		// TODO Auto-generated method stub

	}

	@SuppressWarnings("serial")
	class CustomDialog extends JFrame implements
	ActionListener, ChangeListener, DocumentListener, ItemListener {
		int IW = imp.getWidth();
		int IH = imp.getWidth();

		// Default lane size/offset (Just center 4 lanes in the image)
		int nLanes = 8;
		int LW     = (int) Math.round(0.9*IW/nLanes);
		int LH     = (int) Math.round(IH*0.9);
		int LSp    = (int) Math.round((IW-LW*nLanes)/(nLanes+1));
		int LHOff  = LSp;
		int LVOff  = (IH-LH)/2;

		double tolBG = 0.3; 
		double tolPK = 0.1;

		private JPanel     textPanel;
		private JLabel     labelNLanes;
		private JTextField textNLanes; 

		private JPanel    sliderPanel;
		private JSlider   sliderW;
		private JSlider   sliderH;
		private JSlider   sliderSp;
		private JSlider   sliderHOff;
		private JSlider   sliderVOff;

		private JPanel    buttonPanel;
		private JButton   buttonMeasure;
		private JCheckBox chkBoxBands;
		private JLabel     labelTolBG;
		private JLabel     labelTolPK;
		private JTextField textTolBG;
		private JTextField textTolPK;

		private JPanel    dialogPanel;
		private JFrame    frame = new JFrame("Gel Lanes Gauss Fitting:" + imp.getTitle());

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
			sliderW  = makeTitledSlider("Width ( "+ LW +" px )",                  Color.black, 1, IW/4, LW);
			sliderPanel.add(sliderW);
			sliderH  = makeTitledSlider("Height ( "+ LH +" px )",                 Color.black, IH/10, IH, LH);
			sliderPanel.add(sliderH);
			sliderSp = makeTitledSlider("Space ( "+ LSp +" px )",                 Color.black, 1, IW/nLanes, LSp);
			sliderPanel.add(sliderSp);
			sliderHOff = makeTitledSlider("Horizontal Offset ( "+ LHOff +" px )", Color.black, 0, (int) Math.round(IW*0.9), LHOff);
			sliderPanel.add(sliderHOff);
			sliderVOff = makeTitledSlider("Vertical Offset ( "+ LVOff +" px )",   Color.black, 0, (int) Math.round(IH*0.9), LVOff);
			sliderPanel.add(sliderVOff);

			buttonPanel = new JPanel();
			buttonPanel.setLayout(new BoxLayout(buttonPanel, BoxLayout.X_AXIS));
			buttonMeasure = new JButton("Measure");
			buttonMeasure.addActionListener(this);
			buttonPanel.add(buttonMeasure);
			chkBoxBands = new JCheckBox("Show Bands");
			chkBoxBands.addItemListener(this);
			chkBoxBands.setSelected(false);
			//chkBoxBands.setVisible(false);
			chkBoxBands.setEnabled(false);
			buttonPanel.add(chkBoxBands);
			labelTolBG = new JLabel("Tol BG");
			labelTolPK = new JLabel("Tol PK");
			textTolBG  = new JTextField(3); textTolBG.setText("" + tolBG);
			textTolPK  = new JTextField(3); textTolPK.setText("" + tolPK);
			buttonPanel.add(labelTolBG);
			buttonPanel.add(textTolBG);
			buttonPanel.add(textTolPK);
			buttonPanel.add(labelTolPK);

			dialogPanel = new JPanel();
			dialogPanel.setBackground(Color.lightGray);
			dialogPanel.setLayout(new BorderLayout());
			dialogPanel.add(textPanel,   BorderLayout.NORTH);
			dialogPanel.add(sliderPanel, BorderLayout.CENTER);
			dialogPanel.add(buttonPanel, BorderLayout.SOUTH);

			frame.addWindowListener(new WindowAdapter() {
				public void windowClosing(WindowEvent e) {
					cleanup();
					frame.dispose();
				}
			});

			frame.setLocation(0,0);
			//frame.setJMenuBar(new Menu());
			frame.getContentPane().add(dialogPanel);
			frame.setResizable(true);
			frame.validate();
			frame.pack();

			ImageWindow iwin = imp.getWindow();
			if (iwin==null) return;
			Dimension screen     = Toolkit.getDefaultToolkit().getScreenSize();
			Dimension imageSize  = iwin.getSize();
			Dimension dialogSize = frame.getSize();
			Point imageLoc = iwin.getLocation();
			int x = imageLoc.x+imageSize.width+10;
			if (x+dialogSize.width>screen.width)
				x = screen.width-dialogSize.width;
			frame.setLocation(x, imageLoc.y);
			ImageCanvas canvas = iwin.getCanvas();
			canvas.requestFocus();
			frame.setVisible(true);
			reDrawROIs(imp,LW,LH,LSp,LHOff,LVOff);
		}


		public void cleanup() {
			//imp.removeImageListener(listener);

		}


		private JSlider makeTitledSlider(String string, Color color, int minVal, int maxVal, int val) {
			//Border empty = BorderFactory.createTitledBorder( BorderFactory.createEmptyBorder() );

			JSlider slider = new JSlider(JSlider.HORIZONTAL, minVal, maxVal, val );
			TitledBorder tb = new TitledBorder(BorderFactory.createEtchedBorder(), 
					//empty,
					"", TitledBorder.CENTER, TitledBorder.BELOW_TOP,
					new Font("Sans", Font.PLAIN, 11));
			tb.setTitle(string);
			tb.setTitleJustification(TitledBorder.LEFT);
			tb.setTitleColor(color);
			slider.setBorder(tb);
			slider.setMajorTickSpacing((maxVal - minVal)/10 );
			slider.setPaintTicks(true);
			slider.addChangeListener(this);
			return slider;
		}

		private void setSliderTitle(JSlider slider, Color color, String str) {
			//Border empty = BorderFactory.createTitledBorder( BorderFactory.createEmptyBorder() );
			TitledBorder tb = new TitledBorder(BorderFactory.createEtchedBorder(), //empty,
					"", TitledBorder.CENTER, TitledBorder.BELOW_TOP,
					new Font("Sans", Font.PLAIN, 11));
			tb.setTitleJustification(TitledBorder.LEFT);
			tb.setTitle(str);
			tb.setTitleColor(color);
			slider.setBorder(tb);
		}

		@Override
		public synchronized void stateChanged(ChangeEvent e) {
			JSlider slider = (JSlider) e.getSource();
			if (slider == sliderW ) {
				LW = sliderW.getValue();
				String str = "Width ( "+ LW +" px )"; 
				setSliderTitle(sliderW, Color.black, str );
			}
			else if (slider == sliderH) {
				LH = sliderH.getValue();
				String str = "Height ( "+ LH +" px )"; 
				setSliderTitle(sliderH, Color.black, str );
			}
			else if (slider == sliderSp) {
				LSp = sliderSp.getValue();
				String str = "Spacing ( "+ LSp +" px )"; 
				setSliderTitle(sliderSp, Color.black, str );
			}
			else if (slider == sliderHOff) {
				LHOff = sliderHOff.getValue();
				String str = "Horizontal Offset ( "+ LHOff +" px )"; 
				setSliderTitle(sliderHOff, Color.black, str );
			}
			else if (slider == sliderVOff) {
				LVOff = sliderVOff.getValue();
				String str = "Vertical Offset ( "+ LVOff +" px )"; 
				setSliderTitle(sliderVOff, Color.black, str );
			}
			reDrawROIs(imp,LW,LH,LSp,LHOff,LVOff);
			//IJ.wait(50);
			//delay to make sure the roi has been updated
			plots = null;
			updatePlots(rows,cols,subMin);
		}

		private void reDrawROIs(ImagePlus image, int lw, int lh, int lsp, int lhoff, int lvoff) {
			if (!WindowManager.getCurrentImage().equals(image)) {
				WindowManager.setTempCurrentImage(image);
			}
			if (imp.getRoi() != null || roiMan.getCount() !=0) {
				roiMan.reset();
				imp.deleteRoi();
			}

			try {
				//IJ.wait(50); //make sure the user has finished entering a number 
				nLanes = (int) Integer.parseInt(textNLanes.getText().trim());
			} catch (NumberFormatException e1) {
				//log.error("Cannot parse Number of Lanes: \""+textNLanes.getText().trim()+ "\"");
			}
			for (int i=0;i<nLanes;i++) {
				Roi roi = new Roi(lhoff+lw*i+lsp*i, lvoff, lw, lh);
				//roi.setFillColor(new Color(r/255,g/255,b/255,alpha));
				roiMan.addRoi(roi);
			}
			//			for (int i=0;i<roiMan.getCount();i++) {
			//				//roiMan.select(i, true, false);
			//				roiMan.add(imp,roiMan.getRoi(i),i);
			//			}
		}

		public int getNLanes(){
			return nLanes;
		}

		public double getTolBG(){
			try {
				tolBG = (double) Double.parseDouble(textTolBG.getText().trim());
				return tolBG;
			} catch (NumberFormatException e1) {
				return 1.0;
			}
		}

		public double getTolPK(){
			try {
				tolPK = (double) Double.parseDouble(textTolPK.getText().trim());
				return tolPK;
			} catch (NumberFormatException e1) {
				//log.error("Cannot parse Number of Lanes: \""+textNLanes.getText().trim()+ "\"");
				return 1.0;
			}
		}

		boolean isRoi() {
			if (imp==null)
				return false;
			Roi roi = imp.getRoi();
			if (roi==null)
				return false;
			return roi.getType()==Roi.LINE || roi.getType()==Roi.RECTANGLE;
		}


		// TextField Listeners
		@Override
		public void changedUpdate(DocumentEvent e) {reDrawROIs(imp,LW,LH,LSp,LHOff,LVOff);}
		@Override
		public void removeUpdate (DocumentEvent e) {reDrawROIs(imp,LW,LH,LSp,LHOff,LVOff);}
		@Override
		public void insertUpdate (DocumentEvent e) {reDrawROIs(imp,LW,LH,LSp,LHOff,LVOff);}

		// Measure Button
		@Override
		public void actionPerformed(ActionEvent e) {
			if (e.getSource().equals(buttonMeasure)){
				try {
					doFit();
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
				updatePlots(rows,cols,subMin);
				chkBoxBands.setEnabled(true);
			}
		}

		@Override
		public void itemStateChanged(ItemEvent e) {
			if (e.getItemSelectable() == chkBoxBands){
				// TODO Auto-generated method stub
			}
		}
	}
}
