/**
 * Gauss Fit
 * GelLanesFit.java
 * author: Rick Ziraldo, 2017
 * The /University of Texas at Dallas, Richardson, TX
 * http://www.utdallas.edu
 *
 * Feature: Fitting of multiple Gaussian functions to intensity profiles along the gel lanes
 * Gel Lanes Fit is a tool for fitting gaussian profiles and estimating
 * the profile parameters on selected lanes in gel electrophoresis images.
 *
 * The GaussianArrayCurveFitter class is implemented using
 * Abstract classes from Apache Commons project
 *
 * The source code is maintained and made available on GitHub
 * https://github.com/rickud/gauss-curve-fit
 *
 */

package gellanesfit;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.Enumeration;
import java.util.jar.Attributes;
import java.util.jar.Manifest;
import java.util.prefs.Preferences;

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
import ij.process.ImageConverter;
import ij.process.LUT;
import net.imagej.ImageJ;

@Plugin(type = Command.class, headless = true,
	menuPath = "Plugins>Gel Tools>Gel Lanes Fit")
@SuppressWarnings("ucd")
public class GelLanesFit implements Command {

	@Parameter
	private LogService log;
	@Parameter
	private DisplayService displayServ;
	@Parameter
	private StatusService statusServ;
	@Parameter
	private static Context context;

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
		final Preferences prefs = Preferences.userRoot().node(this.getClass()
			.getName());
		final double SW = IJ.getScreenSize().getWidth();
		final double SH = IJ.getScreenSize().getHeight();

		imp = IJ.getImage();

		final String impName = imp.getTitle().substring(0, imp.getTitle().indexOf(
			"."));
		about();
		final String title = "[" + version + "] Gel Lanes Fit - " + imp.getTitle();

		final Fitter fitter = new Fitter(context, impName);
		final Plotter plotter = new Plotter(context, imp);
		new MainDialog(context, title, imp, prefs, plotter, fitter);

		imp.getCanvas().requestFocus();
		final ImageWindow iwin = imp.getWindow();
		if (iwin == null) return;

		iwin.setLocation(0, (int) SH / 2);
		iwin.setSize((int) SW / 2, (int) SH / 2);
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
	}

	/**
	 * General info About the Software
	 */
	public void about() {
		try {
			final Enumeration<URL> resources = getClass().getClassLoader()
				.getResources("META-INF/MANIFEST.MF");
			while (resources.hasMoreElements()) {
				try {
					final Manifest manifest = new Manifest(resources.nextElement()
						.openStream());
					// check that this is your manifest and do what you need or get the
					// next one
					final Attributes a = manifest.getMainAttributes();
					final String name = a.getValue("Implementation-Title");
					if (name == null) continue;
					if (name.equals("Gel Lanes Fit")) {
						log.info(name);
						version = a.getValue("Implementation-Version");
						log.info(name + " " + version);
					}
				} catch (final IOException e) {
					log.info("Manifest not found");
				}
			}
		} catch (final IOException e) {
			log.info("Manifest not found");
		}
	}

	
	/**
	 * Main method to execute the plugin in Eclipse
	 *
	 * @param args
	 * @throws Exception
	 * 
	 */
	public static void main(final String... args) throws Exception {
		// create the ImageJ application context with all available services
		final ImageJ ij = new ImageJ();
		ij.launch(args);
		final String sep = File.separator;
		final String folder = "src" + sep + "main" + sep + "resources" + sep +
			"sample" + sep + "tape" + sep;
		String file = "042021_16-02-23_unanalyzed.tif";
		
		final ImagePlus iPlus = new Opener().openImage(folder + sep + file);
		if (iPlus.getType() != ImagePlus.GRAY8 && iPlus.getType() != ImagePlus.GRAY16) {
			ImageConverter ic = new ImageConverter(iPlus);
			ic.convertToGray8();
			iPlus.getProcessor().invert();
			//iPlus.updateAndDraw();
		}
		final LUT[] lut = iPlus.getLuts();
		iPlus.setLut(lut[0].createInvertedLut());
		iPlus.show();
		ij.command().run(GelLanesFit.class, true);
	}
}
