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
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.Rectangle;
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
import ij.ImagePlus;
import ij.gui.*;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.gui.RoiListener;
import ij.io.Opener;
import ij.measure.ResultsTable;
import ij.plugin.frame.RoiManager;
import ij.process.ImageProcessor;
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

@SuppressWarnings({ "rawtypes" })
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
	private JFrame frame;
	static int halfSizeDefault = 4;
	static double pixelSizeDefault = 207.0; // "nm"
    // @param shape - fit circle (1) ellipse(2), or ellipse with varying angle (3)
	static int shapeDefault = 1;
    // @param fitmode - algorithm use: NelderMead (1), or Levenberg Marquard (2)
	static int fitModeDefault = 1;
	static int maxIterationsDefault = 500;
	static double cPCFDefault = 1.0;
	static double baseLevelDefault = 100;

	int halfSize;
	int rectSize; 
	double pixelSize; 
	int shape; 
	int fitMode; 
	int maxIterations; 
	double cPCF; 
	double baseLevel; 
	ImagePlus imp;

	public void run() {
		imp = IJ.getImage();
		preprocess(imp);
		doFit(imp);
	}
	
	private void preprocess(ImagePlus impIn) {
		CustomDialog cd = new CustomDialog();		
		cd.showMainDialog();
		frame = new JFrame("Gel Lanes Gauss Fitting:" + impIn.getTitle());
		frame.addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				//cd.cleanup();
				frame.dispose();
			}
		});
		frame.setLocation(0,0);
		//frame.setJMenuBar(new Menu());
		frame.getContentPane().add(cd);
		frame.setResizable(true);
		frame.validate();
		frame.pack();
		frame.setVisible(true);
	}
	

	private void doFit(ImagePlus impIn){
		// Get list of Point ROIs
		Roi roi = impIn.getRoi();
		if (roi instanceof PointRoi){			
			int nPoints = ((PointRoi) roi).getNCoordinates();	
			Rectangle rect = roi.getBounds();
			int[] xP = ((PointRoi) roi).getXCoordinates();
			int[] yP = ((PointRoi) roi).getYCoordinates();
			for (int i=0; i<nPoints; i++){
				xP[i] += rect.x;
				yP[i] += rect.y;
			}
			
//			************ code derived from Nico Stuurman ********************
			
		   /**
		    * Gaussian fit can be run by estimating parameter c (width of Gaussian)
		    * as 1 (circle), 2 (width varies in x and y), or 3 (ellipse) parameters
		    * 
		    * @param shape - fit circle (1) ellipse(2), or ellipse with varying angle (3)
		    * @param fitmode - algorithm use: NelderMead (1), or Levenberg Marquard (2)
		    */
		    GaussianFit gs = new GaussianFit(shape, fitMode);
		    List<GaussianSpotData> resultList = new ArrayList<GaussianSpotData>();
  
		    for (int i = 0; i < nPoints; i++) {
		        GaussianSpotData spot;

		        // Give user feedback
		        ij.IJ.showStatus("Fitting...");
		        ij.IJ.showProgress(i, nPoints);

		        Roi spotRoi = new Roi(xP[i] - halfSize, yP[i] - halfSize, 2 * halfSize + 1, 2 * halfSize + 1);
		        impIn.setRoi(spotRoi, false);

		        ImageProcessor ipRoi = impIn.getProcessor().crop();
		        spot = new GaussianSpotData(ipRoi, 1, 1, 1, 1, 1, xP[i], yP[i]);
		        double[] paramsOut = gs.doGaussianFit(ipRoi, maxIterations);
		        double sx = 0;
		        double sy = 0;
		        double a = 1.0;
		        double theta = 0.0;
		        if (paramsOut.length >= 4) {
		           // normalize the intensity from the Gaussian fit
		           double N = cPCF * paramsOut[GaussianFit.INT] * (2 * Math.PI * paramsOut[GaussianFit.S] * paramsOut[GaussianFit.S]);           
		           double xpc = paramsOut[GaussianFit.XC];
		           double ypc = paramsOut[GaussianFit.YC];
		           double x = (xpc - halfSize + xP[i]) * pixelSize;
		           double y = (ypc - halfSize + yP[i]) * pixelSize;
	              
		           double s = paramsOut[GaussianFit.S] * pixelSize;
		           // express background in photons after base level correction
		           double bgr = cPCF * (paramsOut[GaussianFit.BGR] - baseLevel);
		           // calculate error using formular from Thompson et al (2002)
		           // (dx)2 = (s*s + (a*a/12)) / N + (8*pi*s*s*s*s * b*b) / (a*a*N*N)
		           double sigma = (s * s + (pixelSize * pixelSize) / 12) / 
		                   N + (8 * Math.PI * s * s * s * s * bgr * bgr) / (pixelSize * pixelSize * N * N);
		           sigma = Math.sqrt(sigma);

		           if (paramsOut.length >= 6) {
		               sx = paramsOut[GaussianFit.S1] * pixelSize;
		               sy = paramsOut[GaussianFit.S2] * pixelSize;
		               a = sx/sy;
		            }

		            if (paramsOut.length >= 7) {
		               theta = paramsOut[GaussianFit.S3];
		            }
	                spot.setData(N, bgr, x, y, 0.0, 2 * s, a, theta, sigma);
	                resultList.add(spot);
		        }
		     }
//			END ************ code derived from Nico Stuurman ********************

			// Log results
		    ResultsTable rt = new ResultsTable();	
		    for (int i=0; i<resultList.size(); i++){
		    	GaussianSpotData spot = resultList.get(i);
		    	rt.incrementCounter();
		    	rt.addValue("X", spot.getX());
		    	rt.addValue("Y", spot.getY());
		    	rt.addValue("XC", spot.getXCenter());
		    	rt.addValue("YC", spot.getYCenter());
		    	rt.addValue("Sigma", spot.getSigma());
		    	rt.addValue("Theta", spot.getTheta());
		    	rt.addValue("Width", spot.getWidth());
		    	rt.addValue("A", spot.getA());
		    	rt.addValue("BGrd", spot.getBackground());
		    	rt.addValue("Intens", spot.getIntensity());
		    }	
		    rt.show("Results");
		    
		    // Restore PointROIs
            Roi points = new PointRoi(xP, yP, nPoints);
            ((PointRoi)points).setShowLabels(false);
            impIn.setRoi(points);
		}
		//else{
		//	IJ.showMessage("New Point-ROI selection required");
		//}		
	}

    private boolean showParameterDialog() {
        String[] shapes = {"Circle", "Ellipse", "Ellipse with varying angle"};
        String[] modes = {"NelderMead", "Levenberg Marquard"};
        GenericDialog gd = new GenericDialog("GaussFitOnSpot Parameters");
        gd.addChoice("Shape ", shapes, shapes[shapeDefault-1]);
        gd.addChoice("FitMode ", modes, modes[fitModeDefault-1]);
        gd.addNumericField("Rectangle HalfSize ", halfSizeDefault, 0, 6, "pixels");
        gd.addNumericField("Pixel size ", pixelSizeDefault, 1, 6, "nm");
        gd.addNumericField("Max iterations ", maxIterationsDefault, 0, 6, "");
        gd.addNumericField("cPCF ", cPCFDefault, 1, 6, "(pixel correction)");
        gd.addNumericField("Base level ", baseLevelDefault, 0, 6, "photons");
        gd.showDialog();          //input by the user (or macro) happens here
        
        if (gd.wasCanceled())
             return false;
        shape = gd.getNextChoiceIndex() + 1;
        fitMode = gd.getNextChoiceIndex() + 1;
        halfSize = (int) gd.getNextNumber();
    	rectSize = 2 * halfSize + 1;
        pixelSize = (int) gd.getNextNumber();
        maxIterations = (int) gd.getNextNumber();
        cPCF = (int) gd.getNextNumber();
        baseLevel = (int) gd.getNextNumber();  
        
        return true;
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
	
	
	class CustomDialog extends JPanel implements 
			ChangeListener, DocumentListener, RoiListener{
		RoiManager roiMan = new RoiManager();
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
		private JButton buttonProfile;
		
		private void showMainDialog() {			
			setBackground(Color.lightGray);
			setLayout(new BorderLayout());
			
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
			buttonProfile = new JButton("Profile");
			buttonPanel.add(buttonProfile);
			
			add(textPanel,   BorderLayout.NORTH);
			add(sliderPanel, BorderLayout.CENTER);
			add(buttonPanel, BorderLayout.SOUTH);
			reDrawROIs(LW,LH,LSp,LHOff,LVOff);
		}


		private JSlider makeTitledSlider(String string, Color color, int minVal, int maxVal, int val) {
			//Border empty = BorderFactory.createTitledBorder( BorderFactory.createEmptyBorder() );

			JSlider slider = new JSlider(JSlider.HORIZONTAL, minVal, maxVal, val );
			TitledBorder tb = new TitledBorder(BorderFactory.createEtchedBorder(), 
					//empty,
					"", TitledBorder.CENTER, TitledBorder.ABOVE_BOTTOM,
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
					"", TitledBorder.CENTER, TitledBorder.ABOVE_BOTTOM,
					new Font("Sans", Font.PLAIN, 11));
			tb.setTitleJustification(TitledBorder.LEFT);
			tb.setTitle(str);
			tb.setTitleColor(color);
			slider.setBorder(tb);
		}


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
		@Override
	    public void changedUpdate(DocumentEvent e) {readNLanes();}
		@Override
		public void removeUpdate(DocumentEvent e)  {readNLanes();}
		@Override
	    public void insertUpdate(DocumentEvent e)  {readNLanes();}
		private void readNLanes(){
			try {
				nLanes = (int) Integer.parseInt(textNLanes.getText().trim());
				reDrawROIs(LW,LH,LSp,LHOff,LVOff);
			} catch (NumberFormatException e1) {
				log.error("Cannot parse Number of Lanes: \""+textNLanes.getText().trim()+ "\"");
			}
		}
		@Override
		public void roiModified(ImagePlus arg0, int arg1) {
			// TODO Auto-generated method stub
			
		}
	}
}
