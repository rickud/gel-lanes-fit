package gausscurvefit;

import java.awt.Color;
import java.awt.event.WindowEvent;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.gui.ProfilePlot;
import ij.gui.Roi;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;

public class Plotter {
	private final ImagePlus plotImage = new ImagePlus();
	private final ImagePlus imp;
	private ArrayList<MyPlot> plots;
	
	private final int rows = 2; // Number of plot Rows in display
	private final int cols = 2; // Number of plot Rows in display
	
	// Color are listed here for consistent, easy modification
	public static final Color vLineRegColor = new Color(175, 175, 175);
	public static final Color vLineAddPointColor = new Color(0, 185, 19);
	public static final Color vLineRemovePointColor = new Color(185, 7, 0);
	// Colors for DataSeries in plots
	public static final Color profileColor = Color.BLACK;
	public static final Color gaussColor = Color.RED;
	public static final Color bgColor = Color.BLUE;
	public static final Color fittedColor = new Color(255, 153, 0);
	// Background Color of selected plot
	public static final Color plotSelColor = new Color(240, 240, 240);
	public static final Color plotAddSelColor = new Color(192, 255, 185);
	public static final Color plotRemoveSelColor = new Color(255, 188, 185);
	
	public Plotter(ImagePlus imp) {
		this.imp = imp;
		this.plots = new ArrayList<>();
		plotImage.setTitle("Profiles of " + imp.getShortTitle());
	}
	
	public void create(ArrayList<Roi> rois) {
		for (final Roi r : rois)
			updateProfile(r);
		plotsMontage();
		plotImage.show();
	}
	
	/**
	 * Profile data from Roi, ready for fitting (null if not possible)
	 *
	 * @param roi
	 */
	private DataSeries getLaneProfile(final Roi roi) {
		imp.setRoi(roi);
		final ProfilePlot profileP = new ProfilePlot(imp, true); // get the profile
		final RealVector profile = new ArrayRealVector(profileP.getProfile());
		if (profile.getDimension() < 2) return null;
		imp.killRoi();
		final double y0 = roi.getBounds().getMinY();
		final double[] y = new double[profile.getDimension()];
		for (int p = 0; p < y.length; p++)
			y[p] = y0 + p;
		return new DataSeries("Lane Profile", DataSeries.PROFILE,
			new ArrayRealVector(y), profile, Color.BLACK);
	}

	public void setVLine(final int laneNumber, final double x,
		final Color color)
	{
		final Iterator<MyPlot> plotIter = plots.iterator();
		while (plotIter.hasNext()) {
			final MyPlot plot = plotIter.next();
			boolean vline = false;
			if (plot.getNumber() == laneNumber) {
				final double[] minmax = plot.getDataMinMax();
				final RealVector xvals = new ArrayRealVector(new double[] { x, x });
				final RealVector yvals = new ArrayRealVector(new double[] { minmax[2],
					minmax[3] });
				final DataSeries l = new DataSeries("", DataSeries.VLINE, xvals, yvals,
					color);
				for (final DataSeries d : plot.getDataSeries()) {
					if (d.getType() == DataSeries.VLINE) {
						d.setX(xvals.toArray());
						d.setY(yvals.toArray());
						vline = true;
						break;
					}
				}
				if (!vline) {
					plot.setSelected(true);
					plot.addDataSeries(l);
				}
				plot.updatePlot();
				break;
			}
		}
	}

	public void removeVLine(final int laneNumber) {
		final Iterator<MyPlot> plotIter = plots.iterator();
		while (plotIter.hasNext()) {
			final MyPlot plot = plotIter.next();
			if (plot.getNumber() == laneNumber) {
				for (final DataSeries d : plot.getDataSeries()) {
					if (d.getType() == DataSeries.VLINE) {
						plot.removeDataSeries(d);
						plot.setSelected(false);
						plot.updatePlot();
						break;
					}
				}
				break;
			}
		}
	}

	/**
	 * Method to output plot collage to plotImage
	 */
	public void plotsMontage() {
		// Construct the plot image
		final int plotSpacing = 5; // black border
		final ImageProcessor plot0 = plots.get(0).getPlot().getProcessor();
		final int plotW = plot0.getWidth();
		final int plotH = plot0.getHeight();
		final int plotsWidth = plotSpacing + cols * (plotW + plotSpacing);
		final int plotsHeight = plotSpacing + rows * (plotH + plotSpacing);
		final ImageProcessor blank = plot0.duplicate();
		IJ.setForegroundColor(255, 255, 255);
		blank.fill();
	
		final ArrayList<ImageProcessor> pageMontages = new ArrayList<>();
		final Iterator<MyPlot> plotIter = plots.iterator();
		int psel = 0; // selected page
		while (plotIter.hasNext()) {
			final ImageProcessor plotsMontage = new ColorProcessor(plotsWidth,
				plotsHeight);
			MyPlot subplot;
			ImageProcessor subplotProcessor;
			String subplotTitle;
			for (int r = 0; r < rows; r++) {
				for (int c = 0; c < cols; c++) {
					if (plotIter.hasNext()) {
						subplot = plotIter.next();
						if (subplot.getPlot().getProcessor() == null) System.out.println("Stop!");
						if (subplot.isSelected()) psel = pageMontages.size() + 1;
						subplotProcessor = subplot.getPlot().getProcessor();
						subplotTitle = subplot.getPlot().getTitle();
					}
					else {
						subplotProcessor = blank;
						subplotTitle = "";
					}
					plotsMontage.insert(subplotProcessor, plotSpacing + c * (plotSpacing +
						plotW), plotSpacing + r * (plotSpacing + plotH));
	
					plotsMontage.drawString(subplotTitle, plotSpacing + plotW / 2 + c *
						(plotSpacing + plotW), plotSpacing + plotH / 5 + r * (plotSpacing +
							plotH), Color.LIGHT_GRAY);
				}
			}
			pageMontages.add(plotsMontage);
		}
		plotImage.setProcessor(pageMontages.get(0));
		final ImageStack plotStack = plotImage.createEmptyStack();
		final Iterator<ImageProcessor> montageIter = pageMontages.iterator();
		int p = 0;
		while (montageIter.hasNext()) {
			plotStack.addSlice("Lanes " + (rows * cols * p + 1) + "-" + (rows * cols *
				p + rows * cols), montageIter.next(), p);
			p++;
		}
		plotImage.setStack(plotStack);
		if (psel != 0) plotImage.setSlice(psel);
		plotImage.updateImage();
		plotImage.show();
		if (plotImage.getRoi() != null) plotImage.killRoi();
	}

	public void updateProfile(final Roi roi) {
		final int roiNumber = Integer.parseInt(roi.getName().substring(5));
		final DataSeries profile = getLaneProfile(roi);
		profile.setColor(profileColor);
		final String xLabel = "Distance (px)";
		final String yLabel = "Grayscale Value";
		final MyPlot newPlot = new MyPlot(roi.getName(), xLabel, yLabel);
		newPlot.setSelectedBGColor(plotSelColor);
		newPlot.addDataSeries(profile);
		newPlot.updatePlot();

		final Iterator<MyPlot> plotIter = plots.iterator();
		while (plotIter.hasNext()) {
			final MyPlot plot = plotIter.next();
			if (plot.getNumber() == (roiNumber)) plotIter.remove();
		}
		plots.add(newPlot);
	}

	public void resetPlots() {
		plots = new ArrayList<>();
	}
	
	public ImagePlus getPlotImage() {
		return plotImage;
	}
	
	public ArrayList<MyPlot> getPlots() {
		return plots;
	}
	
	public ArrayList<DataSeries> getProfiles() {
		ArrayList<DataSeries> profiles = new ArrayList<>();
		for (MyPlot p : plots) {
			ArrayList<DataSeries> dlist = p.getDataSeries();
			for (DataSeries d : dlist) {
				if (d.getType() == DataSeries.PROFILE) profiles.add(d);
			}
		}
		return profiles;
	}
	
	public void closePlot(){ 
		if (plotImage.getWindow() != null) plotImage.getWindow().close();
	}
}

class DataSeries implements Comparable<DataSeries> {

	private Color color; // Plot color
	private RealVector x;
	private RealVector y;
	private final String name; // Name for Legend
	private int type; // Type of function

	final static int VLINE = 0;
	final static int PROFILE = 1;
	final static int BACKGROUND = 2;
	final static int GAUSS_BG = 3;
	final static int FITTED = 4;

	public DataSeries(final String name, final int type, final RealVector x,
		final RealVector y, final Color color)
	{
		this.name = name;
		this.type = type;
		this.x = x;
		this.y = y;
		this.color = color;
		if (x.getDimension() == 2 && x.getEntry(0) == x.getEntry(1)) {
			this.type = VLINE;
		}
		else this.type = PROFILE;
	}

	public DataSeries(final String name, final int type, final RealVector x,
		final UnivariateFunction[] function, final Color color)
	{
		this.name = name;
		this.type = type;
		this.x = x;
		this.color = color;
		for (int i = 0; i < function.length; i++) {
			if (i == 0) {
				this.y = new ArrayRealVector(x.map(function[i]));
			}
			else {
				this.y = y.add(x.map(function[i]));
			}
		}
	}

	public DataSeries(final String name, final int type, final RealVector x,
		final UnivariateFunction function, final Color color)
	{
		this.name = name;
		this.type = type;
		this.x = x;
		this.color = color;
		this.y = new ArrayRealVector(x.map(function));
	}

	public int getType() {
		return type;
	}

	public double[] getX() {
		return x.toArray();
	}

	public double[] getY() {
		return y.toArray();
	}

	public Color getColor() {
		return color;
	}

	public String getName() {
		return name;
	}

	/**
	 * @param color
	 */
	public void setColor(final Color color) {
		this.color = color;
	}

	public void setX(final double[] x) {
		this.x = new ArrayRealVector(x);
	}

	public void setY(final double[] y) {
		this.y = new ArrayRealVector(y);
	}

	@Override
	public int compareTo(final DataSeries d) {
		return (type - d.getType());
	}
}

class MyPlot implements Comparable<MyPlot> {

	private Plot plot;

	private boolean selected = false;
	private final int number;

	private final String title;
	private final String xLabel;
	private final String yLabel;

	private ArrayList<DataSeries> data = new ArrayList<>();
	double xmin;
	double xmax;
	double ymin;
	double ymax;

	private Color selectedBG = Color.white;
	private final Color unselectedBG = Color.white;

	public MyPlot(final String title, final String xLabel,
		final String yLabel)
	{
		this.title = title;
		this.xLabel = xLabel;
		this.yLabel = yLabel;
		this.number = Integer.parseInt(title.substring(5));
	}

	public void updatePlot() {
		final Plot newPlot = new Plot(title, xLabel, yLabel);
		newPlot.useTemplate(plot, Plot.COPY_LABELS + Plot.COPY_LEGEND +
			Plot.COPY_SIZE);
		if (plot != null) plot.dispose();
		plot = newPlot;
		if (selected) plot.setBackgroundColor(selectedBG);
		else plot.setBackgroundColor(unselectedBG);
		for (final DataSeries d : data) {
			plot.setColor(d.getColor());
			plot.addPoints(d.getX(), d.getY(), Plot.LINE);
		}
		plot.setLimitsToFit(true);
	}

	private void sortDataSeries() {
		Collections.sort(data);
		final Iterator<DataSeries> dataIter = data.iterator();
		int nullIndex = 0;
		xmin = Double.MAX_VALUE;
		xmax = Double.MIN_VALUE;
		ymin = Double.MAX_VALUE;
		ymax = Double.MIN_VALUE;
		while (dataIter.hasNext()) {
			final DataSeries d = dataIter.next();
			if (d == null) break;
			double xmin2, xmax2, ymin2, ymax2;
			final RealVector x = new ArrayRealVector(d.getX());
			final RealVector y = new ArrayRealVector(d.getY());
			xmin2 = x.getMinValue();
			xmax2 = x.getMaxValue();
			ymin2 = y.getMinValue();
			ymax2 = y.getMaxValue();
			if (xmin2 < xmin) xmin = xmin2;
			if (xmax2 > xmax) xmax = xmax2;
			if (ymin2 < ymin) ymin = ymin2;
			if (ymax2 > ymax) ymax = ymax2;
			nullIndex++;
		}
		data = new ArrayList<>(data.subList(0, nullIndex));
	}

	public void addDataSeries(final DataSeries d) {
		data.add(d);
		sortDataSeries();
	}

	public void addDataSeries(final ArrayList<DataSeries> d) {
		data.addAll(d);
		sortDataSeries();
	}

	public void removeDataSeries(final DataSeries d) {
		data.remove(d);
		sortDataSeries();
	}

	public ArrayList<DataSeries> getDataSeries() {
		return data;
	}

	public Plot getPlot() {
		return plot;
	}

	public int getNumber() {
		return number;
	}

	public double[] getDataMinMax() {
		return new double[] { xmin, xmax, ymin, ymax };
	}
	
	public boolean isSelected() {
		return selected;
	}
	public void setSelected(final boolean s) {
		selected = s;
	}

	public void setSelectedBGColor(final Color color) {
		selectedBG = color;
	}

	@Override
	public int compareTo(final MyPlot p) {
		return (number - p.getNumber());
	}
}
