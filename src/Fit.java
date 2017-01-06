import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.List;

import gauss.GaussianFit;
import gauss.GaussianSpotData;
import ij.ImagePlus;
import ij.gui.PointRoi;
import ij.gui.Roi;
import ij.measure.ResultsTable;
import ij.process.ImageProcessor;

public class Fit {
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
}
