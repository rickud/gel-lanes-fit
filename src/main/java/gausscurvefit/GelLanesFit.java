/**
 * Gauss Fit
 * GelLanesFit.java
 * author: Rick Ziraldo, 2017
 * The /University of Texas at Dallas, Richardson, TX
 * http://www.utdallas.edu
 *
 * Feature: Fitting of multiple Gaussian functions to intensity profiles along the gel lanes
 * Gauss Fit is a tool for fitting gaussian profiles and estimating
 * the profile parameters on selected lanes in gel electrophoresis images.
 *
 *    The GaussianArrayCurveFitter class is implemented using
 *    Abstract classes from Apache Commons project
 *
 *    The source code is maintained and made available on GitHub
 *    https://github.com/rickud/gauss-curve-fit
 */

package gausscurvefit;

import java.io.IOException;
import java.io.InputStream;
import java.util.Properties;

import net.imagej.ImageJ;

import org.scijava.Context;
import org.scijava.app.StatusService;
import org.scijava.command.Command;
import org.scijava.display.DisplayService;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.ImageWindow;
import ij.io.Opener;

@Plugin(type = Command.class, headless = true,
	menuPath = "Plugins>Gel Tools>Gel Lanes Fit")
public class GelLanesFit implements Command {

	/**
	 *
	 */
	@Parameter
	private LogService log;
	@Parameter
	private DisplayService displayServ;
	@Parameter
	private StatusService statusServ;
	@Parameter
	private static Context context;

	// Default Parameters
	// TODO: private Thread mainWindowThread; // thread for the main window
	private Thread plotThread; // thread for plotting

	private boolean setup = true;
	// TODO: private boolean doPlot; // tells the background thread to update

	private ImagePlus imp;
	private String version;

	/**
	 * Initialization method
	 */
	public void init() {
		final double SW = IJ.getScreenSize().getWidth();
		final double SH = IJ.getScreenSize().getHeight();

		imp = IJ.getImage();
		final String impName = imp.getTitle().substring(0, imp.getTitle().indexOf(
			"."));
		about();
		final String title = "[v" + version + "] Gel Lanes Fit: " + imp
			.getTitle();
		final MainDialog md = new MainDialog(context, title, imp, SH, SW);

		final Fitter fitter = new Fitter(context, impName, md.getDegBG(), md
			.getTolPK());
		final Plotter plotter = new Plotter(context, imp, md.getRois(), SH, SW);
		md.setPlotter(plotter);
		md.setFitter(fitter);

		imp.getCanvas().requestFocus();
		final ImageWindow iwin = imp.getWindow();
		if (iwin == null) return;
		
		iwin.setLocation(0, (int) SH/2);
		iwin.setSize((int) SW/2, (int) SH/2);
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
	 * General info About the Software
	 */
	public void about() {
		final ClassLoader loader = Thread.currentThread().getContextClassLoader();
		final Properties props = new Properties();
		try (InputStream input = loader.getResourceAsStream("about.properties")) {
			props.load(input);
			version = props.getProperty("version");
			log.info("Gauss Fit - v" + version + "\n");
			input.close();
		}
		catch (final IOException e) {
			e.printStackTrace();
			log.info("Gel Lanes Fit - v[unknown]");
		}
	}

	/**
	 * Main method to execute the plugin in Eclipse
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(final String... args) throws Exception {
		// create the ImageJ application context with all available services
		final ImageJ ij = net.imagej.Main.launch(args);
		final ImagePlus iPlus = new Opener().openImage(
			"src//main//resources//sample//SumTimeSeriesBright[100-40ms].tif");
		iPlus.show();
		ij.command().run(GelLanesFit.class, true);
	}
}
