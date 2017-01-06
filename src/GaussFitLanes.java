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
import java.awt.GridLayout;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.lang.Math;
import java.util.ArrayList;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JButton;
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

import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.command.Previewable;
import org.scijava.convert.ConvertService;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import gauss.GaussianFit;
import gauss.GaussianSpotData;
import ij.IJ;
import ij.ImageListener;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.*;
import ij.io.Opener;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;
import ij.util.Tools;
import io.scif.img.ImgIOException;
import io.scif.img.ImgOpener;
import io.scif.services.DatasetIOService;
import net.imagej.DatasetService;
import net.imagej.ImageJ;
import net.imagej.display.ImageDisplayService;
import net.imagej.ops.OpService;
import net.imglib2.img.ImagePlusAdapter;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;

//import net.imglib2.type.Type;


@Plugin(type = Command.class, headless = true,
menuPath = "Plugins>Gauss Fit")
public class GaussFitLanes implements Command, Previewable{
	
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
    private Thread bgThread;                //thread for plotting (in the background)
    private boolean setup = true;
    private boolean doUpdate;               //tells the background thread to update
	ImagePlus imp;

	private ImagePlus plotImage;
	private RoiManager roiMan = new RoiManager();
	

	public void init() {
		imp = IJ.getImage();
        
		CustomDialog cd = new CustomDialog();		
		cd.showMainDialog();

		ImageProcessor ip = getProfilePlot();  // get a profile
        if (ip==null) {                     // no profile?
            IJ.error("Dynamic Profiler","No Profile Obtained"); return;
        }

        plotImage = new ImagePlus("Profile of "+imp.getShortTitle(), ip);
        plotImage.show();
	}
	
    // the background thread for plotting.
	public void run() {
		if (setup) { init(); setup = false;}
		while (true) {
			IJ.wait(50);                            //delay to make sure the roi has been updated
			for (int i=0; i<roiMan.getCount(); i++) {
				ImageProcessor ip = getProfilePlot();
				if (ip != null) plotImage.setProcessor(null, ip);
				synchronized(this) {
					if (doUpdate) {
						doUpdate = false;               //and loop again
					} else {
						try {wait();}                   //notify wakes up the thread
						catch(InterruptedException e) { //interrupted tells the thread to exit
							return;
						}
					}
				}
			}
		}
	}

	
    /** get a profile, analyze it and return a plot (or null if not possible) */
	ImageProcessor getProfilePlot() {
		if (roiMan.getCount() == 0) return null;
		ImageProcessor ip = imp.getProcessor();
		for (int p = 0; p<roiMan.getCount(); p++) {
			roiMan.select(p);
			Roi roi = roiMan.getRoi(p);
			if (ip == null || roi == null) return null; //these may change asynchronously
			if (roi.getType() == Roi.LINE)
				ip.setInterpolate(PlotWindow.interpolate);
			else
				ip.setInterpolate(false);
			ProfilePlot profileP = new ProfilePlot(imp, Prefs.verticalProfile);//get the profile
			if (profileP == null) return null;
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

			Plot plot = new Plot("profile", xLabel, yLabel, x, profile);
			double fixedMin = ProfilePlot.getFixedMin();
			double fixedMax = ProfilePlot.getFixedMax();
			if (fixedMin!=0 || fixedMax!=0) {
				double[] a = Tools.getMinMax(x);
				plot.setLimits(a[0],a[1], fixedMin, fixedMax);
			}
		}
		return ip;
	}
    
	private void doFit(ImagePlus imp2) {
		// TODO Auto-generated method stub
		
	}



    
    public static void main(final String... args) throws Exception {
		// create the ImageJ application context with all available services
		final ImageJ ij = net.imagej.Main.launch(args);
		ImagePlus imp = new Opener().openImage("target//01_05_16.tiff");
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
			ActionListener, ChangeListener, DocumentListener, RoiListener, ImageListener{
		
		int IW = imp.getWidth();
		int IH = imp.getWidth();
		
		// Default lane size/offset (Just center 4 lanes in the image)
		int nLanes = 4;
		int LW     = IW/10;
		int LH     = (int) Math.round(IH*0.9);
		int LSp    = IW/10;
		int LHOff  = (IW-(LW+LSp)*nLanes)/2;
		int LVOff  = (IH-LH)/2;
		
		private JPanel textPanel;
		private JLabel labelNLanes = new JLabel("");
		private JTextField textNLanes= new JTextField(3); 
		
		private JPanel sliderPanel;
		private JSlider sliderW;
		private JSlider sliderH;
		private JSlider sliderSp;
		private JSlider sliderHOff;
		private JSlider sliderVOff;
		
		private JPanel buttonPanel;
		private JButton buttonPeaks;
		private JPanel dialogPanel;
		private JFrame frame = new JFrame("Gel Lanes Gauss Fitting:" + imp.getTitle());
		
		
		private void showMainDialog() {
			textPanel = new JPanel();
			textPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
			
			textNLanes.setText("" + nLanes);
			textNLanes.setVisible(true);
			textPanel.add(textNLanes);
			textNLanes.getDocument().addDocumentListener(this);

			labelNLanes.setText("Number of Lanes");
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
			buttonPeaks = new JButton("Peaks");
			buttonPeaks.addActionListener(this);
			buttonPanel.add(buttonPeaks);
			
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

			reDrawROIs(LW,LH,LSp,LHOff,LVOff);
			imp.addImageListener(this);
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
		public void stateChanged(ChangeEvent e) {
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
			reDrawROIs(LW,LH,LSp,LHOff,LVOff);
		}

		private void reDrawROIs(int lw, int lh, int lsp, int lhoff, int lvoff) {
			if (imp.getRoi() != null ) {
				roiMan.reset();
				imp.deleteRoi();
			}
			
			for (int i=0;i<nLanes;i++) {
				Roi roi = new Roi(lhoff+lw*i+lsp*i, lvoff, lw, lh);
				roiMan.addRoi(roi);
			}
			for (int i=0;i<roiMan.getCount();i++) {
				roiMan.select(i, true, false);
				//roiMan.
			}

		}



		private void readNLanes(){
			try {
				nLanes = (int) Integer.parseInt(textNLanes.getText().trim());
				reDrawROIs(LW,LH,LSp,LHOff,LVOff);
			} catch (NumberFormatException e1) {
				log.error("Cannot parse Number of Lanes: \""+textNLanes.getText().trim()+ "\"");
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
		@Override
		public void actionPerformed(ActionEvent e) {
			// TODO Auto-generated method stub
			
		}
	    public void changedUpdate(DocumentEvent e) {readNLanes();}
		@Override
		public void removeUpdate(DocumentEvent e)  {readNLanes();}
		@Override
	    public void insertUpdate(DocumentEvent e)  {readNLanes();}
		@Override
		public void roiModified(ImagePlus arg0, int arg1) {log.info("" + arg0 + " " + arg1);}


		@Override
		public void imageClosed(ImagePlus impl) {
			// if any image we are listening to is closed, exit
			//if (impl == this.imp || imp == plotImage) {
	            //removeListeners();
	            //closePlotImage(); //also terminates the background thread
	        //}
		}

		@Override
		public void imageOpened(ImagePlus arg0) {
			// TODO Auto-generated method stub
			
		}

		@Override
		public synchronized void imageUpdated(ImagePlus arg0) {
			// TODO Auto-generated method stub
			
		}
	}

}
