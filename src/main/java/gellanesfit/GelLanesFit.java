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
import java.io.InputStream;
import java.lang.reflect.InvocationTargetException;
import java.net.URL;
import java.util.Enumeration;
import java.util.jar.Attributes;
import java.util.jar.Manifest;
import java.util.prefs.Preferences;

import javax.swing.SwingUtilities;

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
public class GelLanesFit implements Command {

	@Parameter
	private LogService log;
	@Parameter
	private DisplayService displayServ;
	@Parameter
	private StatusService statusServ;
	@Parameter
	private static Context context;

	private boolean setup = true;

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
	}

	@Override
	public void run() {
		if (!setup) return;
		setup = false;
		// Swing is not thread safe: build and show the windows on the event
		// thread, not on the thread SciJava runs the command on
		try {
			if (SwingUtilities.isEventDispatchThread()) init();
			else SwingUtilities.invokeAndWait(this::init);
		}
		catch (final InterruptedException e) {
			Thread.currentThread().interrupt();
		}
		catch (final InvocationTargetException e) {
			log.error("Gel Lanes Fit could not start", e.getCause());
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
				try (InputStream in = resources.nextElement().openStream()) {
					final Manifest manifest = new Manifest(in);
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
			"sample-images" + sep;
		String file = "tagment-test" + sep + "gel-camera-1" + sep + "Long_5s.tif";
		
		final ImagePlus iPlus = new Opener().openImage(folder + sep + file);
		if (iPlus.getType() != ImagePlus.GRAY8 && iPlus.getType() != ImagePlus.GRAY16) {
			ImageConverter ic = new ImageConverter(iPlus);
			ic.convertToGray8();
			iPlus.getProcessor().invert();
		}
		final LUT[] lut = iPlus.getLuts();
		iPlus.setLut(lut[0].createInvertedLut());
		iPlus.show();
		ij.command().run(GelLanesFit.class, true);
	}
}
