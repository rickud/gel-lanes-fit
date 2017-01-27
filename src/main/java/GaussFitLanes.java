/**
 * 
 * GaussFit_OnSpot.java
 * author: by Peter Haub, Nov 2013
 * PD Dr. Tobias Meckel, University Darmstadt, www.bio.tu-darmstadt.de/meckel
 * Peter Haub, phaub@dipsystems.de
 * 
 * 
 * Feature:   Fitting of Gaussian profiles at manually selected PointROIs
 * 
 * GaussFit_OnSpot is a tool for fitting gaussian profiles and estimating
 * the profile parameters on manually selected positions in
 * sub-diffraction-limited microscopic images.
 * 
 * 
 * The code is based on the Gaussian Fitting package developed by Nico Stuurman.
 * Original code and license information can be found on:
 * 		http://valelab.ucsf.edu/~nstuurman/IJplugins/GaussianFit.html
 * 		http://valelab.ucsf.edu/~nstuurman/IJplugins/license.html
 * 
 * Original files from the Gaussian Fitting package included:
 *    GaussianSpotData.java
 *    GaussianFit.java
 *    MultiVariateGaussianFunction.java
 *    ParametricGaussianFunction.java
 * 
 * Files based on the Gaussian Fitting package:
 *    GFUtils.java (contains parts of the GaussianUtils.java)
 * 
 * 
 * History :
 * 
 * Modifications thanks to Wayne Rasband, July 2014
 * 
 * Correction by Peter Haub, Jan. 2015
 * Bug: Identical CenterX and CenterY value
 * Update: setHideLabels() changed to setShowLabels() ( for IJ.version >= 1.49m )
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

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.math.NumberUtils;
import org.apache.commons.math4.analysis.function.Gaussian;
import org.apache.commons.math4.fitting.PolynomialCurveFitter;
import org.apache.commons.math4.fitting.WeightedObservedPoints;
import org.apache.commons.math4.linear.ArrayRealVector;
import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.command.Previewable;
import org.scijava.convert.ConvertService;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import ij.IJ;
import ij.ImagePlus;
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
    private Thread mainThread;                //thread for plotting (in the background)
    private Thread plotThread;                //thread for plotting (in the background)
    private boolean setup = true;
    private boolean doFit;               //tells the background thread to update
	ImagePlus imp;

	private ImagePlus plotImage;
	private RoiManager roiMan = new RoiManager();
	private CustomDialog cd;
	private int rows = 4; // Number of plot Rows in display
	private Plot[] plots;

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
		getProfilePlots();
		ImageProcessor ip = makePlotsMontage(rows);  
        plotImage = new ImagePlus("Profile of "+imp.getShortTitle(),ip);

        if (ip != null) plotImage.setProcessor(null, ip);
		if (plotImage.getRoi()!=null) plotImage.deleteRoi();
		plotImage.show();
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
		}
	}
	
    /** get a profile, analyze it and return a plot (or null if not possible) */
	private Plot[] getProfilePlots() {
		if (roiMan.getCount() == 0) return null;
		ImageProcessor ip = imp.getProcessor();
		plots = new Plot[roiMan.getCount()];
		
		for (int p = 0; p<roiMan.getCount(); p++) {
			roiMan.select(p);
			Roi roi = roiMan.getRoi(p);
			if (ip == null || roi == null) return null; //these may change asynchronously
			if (roi.getType() == Roi.LINE)
				ip.setInterpolate(PlotWindow.interpolate);
			else
				ip.setInterpolate(false);
			ProfilePlot profileP = new ProfilePlot(imp, true);//get the profile
			double[] profile = profileP.getProfile();
			if (profile==null || profile.length<2)
				return null;
			String xUnit = "pixels";                    //the following code is mainly for x calibration
			double xInc = 1;
			Calibration cal = imp.getCalibration(); 
			if (roi.getType() == Roi.LINE) {
				Line line = (Line)roi;
				if (cal != null) {
					double dx = cal.pixelWidth*(line.x2 - line.x1);
					double dy = cal.pixelHeight*(line.y2 - line.y1);
					double length = Math.sqrt(dx*dx + dy*dy);
					xInc = length/(profile.length-1);
					xUnit = cal.getUnits();
				}
			} else if (roi.getType() == Roi.RECTANGLE) {
				if (cal != null) {
					xInc = roi.getBounds().getWidth()*cal.pixelWidth/(profile.length-1);
					xUnit = cal.getUnits();
				}
			} else return null;
			String xLabel = "Distance (" + xUnit + ")";
			String yLabel = (cal !=null && cal.getValueUnit()!=null && !cal.getValueUnit().equals("Gray Value")) ?
					"Value ("+cal.getValueUnit()+")" : "Value";

			int n = profile.length;                 // create the x axis
			double[] x = new double[n];
			for (int i=0; i<n; i++)
				x[i] = i*xInc;

			Plot plot = new Plot("Lane "+(p+1), xLabel, yLabel, x, profile);
			plot.addText(plot.getTitle(), plot.getSize().getWidth()/2, plot.getSize().getHeight()/2);
			double fixedMin = ProfilePlot.getFixedMin();
			double fixedMax = ProfilePlot.getFixedMax();
			if (fixedMin!=0 || fixedMax!=0) {
				double[] a = Tools.getMinMax(x);
				plot.setLimits(a[0],a[1], fixedMin, fixedMax);
			}
			plots[p] = plot;
		}
		return plots;
	}
		
	private ImageProcessor makePlotsMontage(int rows){
		// Construct the plot image
        if (ArrayUtils.contains(plots,null)) {                          // no profile?
            IJ.error("Gauss Fit","No Profiles Obtained"); return null;
        }
		
		int cols = Math.floorDiv(plots.length, rows);
		if (Math.floorMod(plots.length, rows)!=0) cols++;
		
		int plotSpacing = 5; // black border
		int plotW = plots[0].getProcessor().getWidth(); int plotH = plots[0].getProcessor().getHeight();
		int plotsWidth  = plotSpacing + cols*( plotW + plotSpacing );
		int plotsHeight = plotSpacing + rows*( plotH + plotSpacing );
		ImageProcessor blank = plots[0].getProcessor().duplicate();
		IJ.setForegroundColor(255,255,255);
		blank.fill();

		ImageProcessor plotsMontage = new ColorProcessor(plotsWidth, plotsHeight);
		for (int c = 0; c<cols; c++) {
			for (int r = 0; r<rows; r++) {
				if (c*rows+r<plots.length){
					plotsMontage.insert(plots[c*rows+r].getProcessor(), //here NP
										plotSpacing+c*(plotSpacing + plotW), 
										plotSpacing+r*(plotSpacing + plotH));
					plotsMontage.drawString(plots[c*rows+r].getTitle(),
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
		return plotsMontage;
	}
    
	private void doFit() {
		//Reset the plots window
		getProfilePlots();
		makePlotsMontage(rows);
		
		for (int i = 0; i<plots.length; i++) {
			float[] xvalsf   = plots[i].getXValues();
			float[] yvalsf   = plots[i].getYValues();
			double[] xvals   = new double[xvalsf.length];
			double[] yvals   = new double[yvalsf.length];
					
			WeightedObservedPoints obs   = new WeightedObservedPoints(); // All Y
			WeightedObservedPoints obsbg = new WeightedObservedPoints(); // Only minima for BG
			for (int o = 0 ; o<xvals.length; o++) {
				xvals[o] = (double) xvalsf[o];
				yvals[o] = (double) yvalsf[o];
				obs.add(xvals[o], yvals[o]);
			}
			plots[i].setColor(Color.black);
			plots[i].addPoints(xvals, yvals, PlotWindow.LINE);
			double tolbg = cd.getTolBG()*(NumberUtils.max(yvals)-NumberUtils.min(yvals));
			double tolpk = cd.getTolPK()*(NumberUtils.max(yvals)-NumberUtils.min(yvals));
			
			int[] minimaIdx = findMaxima(new ArrayRealVector(yvals).mapMultiplyToSelf(-1).toArray(), 
											tolbg, true);
			
			double[] xvalsbg = new double[minimaIdx.length]; // Valleys to determine BG
			double[] yvalsbg = new double[minimaIdx.length];
			
			for (int o = 0 ; o<minimaIdx.length; o++) {
				xvalsbg[o] = xvals[minimaIdx[o]];
				yvalsbg[o] = yvals[minimaIdx[o]];
				obsbg.add(xvalsbg[o], yvalsbg[o]);
			}
			
			// fit background poly
			double[] bg = new double[xvals.length];
			if (minimaIdx.length>3) {
				double[] pars_bg = PolynomialCurveFitter.create(3).fit(obsbg.toList());
				for (int b=0; b<xvals.length; b++) {
					bg[b] = pars_bg[0] + xvals[b]*(pars_bg[1] + xvals[b]*(pars_bg[2] + xvals[b]*pars_bg[3]));
				}
			} else if (minimaIdx.length>2) {
				double[] pars_bg = PolynomialCurveFitter.create(2).fit(obsbg.toList());
				for (int b=0; b<xvals.length; b++) {
					bg[b] = pars_bg[0] + xvals[b]*(pars_bg[1] + xvals[b]*pars_bg[2]);
				}
			} else if (minimaIdx.length>1) {
				double[] pars_bg = PolynomialCurveFitter.create(1).fit(obsbg.toList());
				for (int b=0; b<xvals.length; b++) {
					bg[b] = pars_bg[0] + xvals[b]*pars_bg[1];
				}
			} else {
				for (int b=0; b<xvals.length; b++) {
					bg[b] = NumberUtils.min(yvals);
				}
			}

			for (int b=0; b<xvals.length; b++) {
				Math.max(yvals[b] -= bg[b], 0);
			}

			// Local maxima where the peaks are
			int[] maximaIdx = findMaxima(yvals, tolpk, true);

			// Define a Gaussian at each peak
			double[] means      = new double[maximaIdx.length];
			double[] stds       = new double[maximaIdx.length];
			double[] norms      = new double[maximaIdx.length];
			double[] guessStart = new double[3*maximaIdx.length]; // Array of parameters for Fitter
			
			// Initial parameters based on the position of the peaks
			for (int m=0; m<maximaIdx.length; m++) {
				norms[m] = yvals[maximaIdx[m]];
				means[m] = xvals[maximaIdx[m]];
				stds [m] = 0.1;
				guessStart[3*m]   = norms[m];
				guessStart[3*m+1] = means[m];
				guessStart[3*m+2] = stds [m];
			}

			//double[] pars = GaussianCurveFitter.create().withStartPoint(ArrayUtils.subarray(guessStart, 0, 3)).fit(obs.toList());
			double[] pars  = GaussianArrayCurveFitter.create().withStartPoint(guessStart).fit(obs.toList());
			double[] gauss = new double[xvals.length];
			if (pars.length != maximaIdx.length*3) {
				System.out.println("Par number mismatch: " +
									pars.length+ " parameters, " +
									maximaIdx.length + " peaks!");
				return;
			}
			
			for (int b=0; b<maximaIdx.length; b++){ // plot one gaussian for each peak
				Gaussian g = new Gaussian(pars[b*3], pars[b*3+1], pars[b*3+2]);
				for (int c=0; c<xvals.length; c++) {
					gauss[c] = g.value(xvals[c])+bg[c];
				}
				plots[i].setColor(Color.red);
				plots[i].addPoints(xvals, gauss, PlotWindow.LINE);
			}
			
			plots[i].setColor(Color.blue);
			plots[i].addPoints(xvals, bg, PlotWindow.LINE);		
			plots[i].setLimitsToFit(true);
			
			String parStr = String.format("Lane %1$d-par: ", i);
			String maxStr = String.format("Lane %1$d-max: ", i);;
			String minStr = String.format("Lane %1$d-min: ", i);;
			for (int pp = 0; pp<pars.length; pp++){
				parStr += String.format("%1$.2f, ", pars[pp]);
			}
			for (int pp = 0; pp<maximaIdx.length; pp++){
				maxStr += String.format("[%1$.2f, %2$.2f] ", xvals[maximaIdx[pp]], yvals[maximaIdx[pp]]);
			}
			for (int pp = 0; pp<minimaIdx.length; pp++){
				minStr += String.format("[%1$.2f, %2$.2f] ", xvals[minimaIdx[pp]], yvals[minimaIdx[pp]]);
			}
			System.out.println(parStr.substring(0,parStr.length()-2));
			System.out.println(maxStr.substring(0,maxStr.length()));
			System.out.println(minStr.substring(0,minStr.length()));
		}	
	}
	
	/**
	 * Adapted From:
	 * Calculates peak positions of 1D array N.Vischer, 13-sep-2013
	 *
	 * @param x Array containing peaks.
	 * @param tolerance Depth of a qualified valley must exceed tolerance.
	 * Tolerance must be >= 0. Flat tops are marked at their centers.
	 * @param  excludeOnEdges If 'true', a peak is only
	 * accepted if it is separated by two qualified valleys. If 'false', a peak
	 * is also accepted if separated by one qualified valley and by a border.
	 * @return Positions of peaks, sorted with decreasing amplitude
	 */
	public static int[] findMaxima(double[] x, double tolerance, boolean includeEnds) {
		int len = x.length;
		if (len<2)
			return new int[0];
		if (tolerance < 0)
			tolerance = 0;
		int[] maxPositions = new int[len];
		double max = x[0];
		double min = x[0];
		int maxPos = 0;
		int lastMaxPos = -1;
		boolean leftValleyFound = includeEnds;
		int maxCount = 0;
		for (int j = 1; j < len; j++) {
			double val = x[j];
			if (val > min + tolerance)
				leftValleyFound = true;
			if (val > max && leftValleyFound) {
				max = val;
				maxPos = j;
			}
			if (leftValleyFound)
				lastMaxPos = maxPos;
			if (val < max - tolerance && leftValleyFound) {
				maxPositions[maxCount] = maxPos;
				maxCount++;
				leftValleyFound = false;
				min = val;
				max = val;
			}
			if (val < min) {
				min = val;
				if (!leftValleyFound)
					max = val;
			}
		}
		if (includeEnds) {
			if (maxCount > 0 && maxPositions[maxCount - 1] != lastMaxPos)
				maxPositions[maxCount++] = lastMaxPos;
			if (maxCount == 0 && max - min >= tolerance)
				maxPositions[maxCount++] = lastMaxPos;
		}
		int[] cropped = new int[maxCount];
		System.arraycopy(maxPositions, 0, cropped, 0, maxCount);
		maxPositions = cropped;
		double[] maxValues = new double[maxCount];
		for (int m = 0; m < maxCount; m++) {
			int pos = maxPositions[m];
			double midPos = pos;
			while (pos < len - 1 && x[pos] == x[pos + 1]) {
				midPos += 0.5;
				pos++;
			}
			maxPositions[m] = (int) midPos;
			maxValues[m] = x[maxPositions[m]];
		}
		int[] rankPositions = Tools.rank(maxValues);
		int[] returnArr = new int[maxCount];
		for (int n = 0; n < maxCount; n++) {
			int pos = maxPositions[rankPositions[n]];
			returnArr[maxCount - n - 1] = pos; //in descending order
		}
		return returnArr;
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
		
		double tolBG = 0.5; 
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
			getProfilePlots();
			ImageProcessor ip = makePlotsMontage(rows);
			if (ip != null) plotImage.setProcessor(ip);
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
				doFit();
				plotImage.setProcessor(null, makePlotsMontage(rows));
				if (plotImage.getRoi()!=null) plotImage.deleteRoi();
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
