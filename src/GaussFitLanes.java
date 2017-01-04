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


import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.List;

import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.command.Previewable;
import org.scijava.convert.ConvertService;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import gauss.GaussianFit;
import gauss.GaussianSpotData;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.io.Opener;
import ij.measure.ResultsTable;
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
		imp = preprocess();
		doFit(imp);
	}
	
	private ImagePlus preprocess() {

		imp = new Opener().openImage("/Users/Rick/Box Sync/Gels/01_05_16.tiff");
		// display it via ImageJ
		imp.show();
		// wrap it into an ImgLib image (no copying)
		final Img image = ImagePlusAdapter.wrap( imp );
		// display it via ImgLib using ImageJ
		ImageJFunctions.show( image );

		showMainDialog();
		return null;
	}
	
    private boolean showMainDialog() {
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
		
	}

	public void cancel() {
		log.info("Process canceled");
		return;
	}

	public void preview() {
		// TODO Auto-generated method stub
		
	}
}
