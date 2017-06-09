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

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.GridLayout;
import java.awt.event.WindowEvent;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.WindowConstants;
import javax.swing.border.EmptyBorder;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.Range;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RectangleEdge;
import org.jfree.ui.RectangleInsets;
import org.jfree.ui.RefineryUtilities;


import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.gui.ProfilePlot;
import ij.gui.Roi;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;

class Plotter2 extends JFrame {

	
	private final ImagePlus plotImage = new ImagePlus();
	private final ImagePlus imp;
	private final ArrayList<MyPlot> plots;

	private final int rows = 2; // Number of plot Rows in display
	private final int cols = 2; // Number of plot Rows in display

	// Colors are listed here for consistent, easy modification
	static final Color vLineRegColor = new Color(127, 127, 127);
	static final Color vLineAddPeakColor = new Color(0, 185, 19);
	static final Color vLineRemovePeakColor = new Color(185, 7, 0);
	// Colors for DataSeries in plots
	private static final Color profileColor = Color.BLACK;
	static final Color gaussColor = Color.RED;
	static final Color bgColor = Color.BLUE;
	static final Color fittedColor = new Color(255, 153, 0);
	// Background Color of selected plot
	static final Color plotSelColor = new Color(240, 240, 240);
	static final Color plotAddSelColor = new Color(192, 255, 185);
	static final Color plotRemoveSelColor = new Color(255, 188, 185);

	Plotter2(final ImagePlus imp, final ArrayList<Roi> rois) {
		this.imp = imp;
		this.plots = new ArrayList<>();
		plotImage.setTitle("Profiles of " + imp.getShortTitle());
		for (final Roi r : rois)
			updateProfile(r);
		plotsMontage();
		plotImage.show();
	}

	public ImagePlus getPlotImage() {
		return plotImage;
	}

	/**
	 * Profile data from Roi, ready for fitting (null if not possible)
	 *
	 * @param roi
	 */
	private DataSeries getLaneProfile(final Roi roi) {
		imp.setRoi(roi);
		final String name = roi.getName();
		final int lane = Integer.parseInt(name.substring(5));
		final ProfilePlot profileP = new ProfilePlot(imp, true); // get the profile
		final RealVector profile = new ArrayRealVector(profileP.getProfile());
		if (profile.getDimension() < 2) return null;
		imp.killRoi();
		final double y0 = roi.getBounds().getMinY();
		final double[] y = new double[profile.getDimension()];
		for (int p = 0; p < y.length; p++)
			y[p] = y0 + p;
		return new DataSeries(name, lane, DataSeries.PROFILE, new ArrayRealVector(
			y), profile, Color.BLACK);
	}

	void updateProfile(final Roi roi) {
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

	void setVLine(final int laneNumber, final double x, final Color color) {
		final Iterator<MyPlot> plotIter = plots.iterator();
		while (plotIter.hasNext()) {
			final MyPlot p = plotIter.next();
			boolean vline = false;
			if (p.getNumber() == laneNumber) {
				final double[] minmax = p.getDataMinMax();
				final RealVector xvals = new ArrayRealVector(new double[] { x, x });
				final RealVector yvals = new ArrayRealVector(new double[] { minmax[2],
					minmax[3] });
				final DataSeries l = new DataSeries("", laneNumber, DataSeries.VLINE,
					xvals, yvals, color);
				for (final DataSeries d : p.getDataSeries()) {
					if (d.getType() == DataSeries.VLINE) {
						d.setX(xvals.toArray());
						d.setY(yvals.toArray());
						vline = true;
						break;
					}
				}
				if (!vline) {
					p.setSelected(true);
					p.addDataSeries(l);
				}
				p.updatePlot();
				break;
			}
		}
	}

	void removeVLine(final int laneNumber) {
		for (final MyPlot p : plots) {
			if (p.getNumber() == laneNumber) {
				final Iterator<DataSeries> dataIter = p.getDataSeries().iterator();
				while (dataIter.hasNext()) {
					if (dataIter.next().getType() == DataSeries.VLINE) {
						dataIter.remove();
						break;
					}
				}
				p.setSelected(false);
				p.updatePlot();
				break;
			}
		}
	}

	private void removeCustomPeaks(final MyPlot p) {
		final Iterator<DataSeries> dataIter = p.getDataSeries().iterator();
		while (dataIter.hasNext()) {
			final DataSeries d = dataIter.next();
			if (d.getType() == DataSeries.CUSTOMPEAKS) {
				dataIter.remove();
			}
		}
		p.updatePlot();
	}

	void updateCustomPeaks(final int lane, final ArrayList<Peak> pl) {
		for (final MyPlot p : plots) {
			if (p.getNumber() == lane) {
				removeCustomPeaks(p);
				RealVector xadd = new ArrayRealVector();
				RealVector yadd = new ArrayRealVector();
				for (final Peak fp : pl) {
					xadd = xadd.append(fp.getMean());
					yadd = yadd.append(fp.getNorm());
				}
				final DataSeries customPeaks = new DataSeries("CP ADD", lane,
					DataSeries.CUSTOMPEAKS, xadd, yadd, Plotter.vLineAddPeakColor);
				p.addDataSeries(customPeaks);
				p.updatePlot();
			}
		}

	}

	/**
	 * Method to output plot collage to plotImage
	 */
	void plotsMontage() {
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
						if (subplot.getPlot().getProcessor() == null) System.out.println(
							"Stop!");
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

	public ArrayList<DataSeries> getProfiles() {
		final ArrayList<DataSeries> profiles = new ArrayList<>();
		for (final MyPlot p : plots) {
			final ArrayList<DataSeries> dlist = p.getDataSeries();
			for (final DataSeries d : dlist) {
				if (d.getType() == DataSeries.PROFILE) profiles.add(d);
			}
		}
		return profiles;
	}

	public ArrayList<MyPlot> getPlots() {
		return plots;
	}

	void closePlot() {
		if (plotImage.getWindow() != null) plotImage.getWindow().close();
	}
}

class Plotter extends JFrame {
	private ImagePlus imp;
	private ArrayList<JFreeChart> plots;
	private ArrayList<ArrayList<DataSeries>> plotsData;
  private JPanel plotsPanel;
	
	private final int rows = 2; // Number of plot Rows in display
	private final int cols = 2; // Number of plot Rows in display

	// Colors are listed here for consistent, easy modification
	static final Color vLineRegColor = new Color(127, 127, 127);
	static final Color vLineAddPeakColor = new Color(0, 185, 19);
	static final Color vLineRemovePeakColor = new Color(185, 7, 0);
	// Colors for DataSeries in plots
	private static final Color profileColor = Color.BLACK;
	static final Color gaussColor = Color.RED;
	static final Color bgColor = Color.BLUE;
	static final Color fittedColor = new Color(255, 153, 0);
	// Background Color of selected plot
	static final Color plotSelColor = new Color(240, 240, 240);
	static final Color plotUnSelColor = Color.WHITE;
	static final Color plotAddSelColor = new Color(192, 255, 185);
	static final Color plotRemoveSelColor = new Color(255, 188, 185);


  public Plotter(ImagePlus imp, ArrayList<Roi> rois) {
		this.imp = imp;
		this.plots = new ArrayList<>();
		this.setTitle("Profiles of " + imp.getShortTitle());
    setDefaultCloseOperation(WindowConstants.HIDE_ON_CLOSE);
    setBounds(100, 100, 700, 500);
    plotsPanel = new JPanel();
    plotsPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
    plotsPanel.setLayout(new GridLayout(rows, cols));
		for (final Roi r : rois) {
			updateProfile(r);
			for (JFreeChart p : plots) {
				if (p.getTitle().getText().equals(r.getName()))
					plotsPanel.add(new ChartPanel(p));
			}
		}
    this.getContentPane().add(plotsPanel, BorderLayout.CENTER);
    this.setVisible(true);
  }


	/**
	 * Profile data from Roi, ready for fitting (null if not possible)
	 *
	 * @param roi
	 */
	private DataSeries getLaneProfile(final Roi roi) {
		imp.setRoi(roi);
		final String name = roi.getName();
		final int lane = Integer.parseInt(name.substring(5));
		final ProfilePlot profileP = new ProfilePlot(imp, true); // get the profile
		final RealVector profile = new ArrayRealVector(profileP.getProfile());
		if (profile.getDimension() < 2) return null;
		imp.killRoi();
		final double y0 = roi.getBounds().getMinY();
		final double[] y = new double[profile.getDimension()];
		for (int p = 0; p < y.length; p++)
			y[p] = y0 + p;
		return new DataSeries(name, lane, DataSeries.PROFILE, new ArrayRealVector(
			y), profile, Plotter.profileColor);
	}
  
	void updateProfile(Roi roi) {
		final DataSeries profile = getLaneProfile(roi);
		final String xLabel = "Distance (px)";
		final String yLabel = "Grayscale Value";
		final XYSeriesCollection dataset = new XYSeriesCollection();
		final XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
		
		dataset.addSeries(profile.getXYSeries());
		renderer.setSeriesPaint(0, Plotter.profileColor);
		final JFreeChart newPlot = ChartFactory.createXYLineChart(roi.getName(),
  		xLabel, yLabel, dataset, PlotOrientation.VERTICAL, true, true, false);
		newPlot.setBackgroundPaint(plotUnSelColor);
		newPlot.getLegend().setVisible(false);
		newPlot.getLegend().setPosition(RectangleEdge.RIGHT);
		newPlot.getTitle().setMargin(new RectangleInsets(15,5,15,5));
		newPlot.getTitle().setPaint(Color.BLACK);
		
		final Iterator<JFreeChart> plotIter = plots.iterator();
		while (plotIter.hasNext()) {
			final JFreeChart plot = plotIter.next();
			if (plot.getTitle().equals(roi.getName())) plotIter.remove();
		}
		plots.add(newPlot);
	}
	
//	void setVLine(final int ln, final double x, final Color color) {
//		
//		for (JFreeChart p : plots) {
//			boolean vline = false;
//			if (p.getTitle().equals("Lane " + ln)) {
//				p.getXYPlot().getDataset().removeSeries();
//				final Range minmax = p.getXYPlot().getDataRange(p.getXYPlot().getRangeAxis());
//				final RealVector xvals = new ArrayRealVector(new double[] { x, x });
//				final RealVector yvals = new ArrayRealVector(new double[] { minmax.getLowerBound(), minmax.getUpperBound() });
//				final XYSeries l = new XYSeries(DataSeries.VLINE, xvals, yvals, color);
//				for (final DataSeries d : p.getDataSeries()) {
//					if (d.getType() == DataSeries.VLINE) {
//						d.setX(xvals.toArray());
//						d.setY(yvals.toArray());
//						vline = true;
//						break;
//					}
//				}
//				if (!vline) {
//					p.setSelected(true);
//					p.addDataSeries(l);
//				}
//
//				break;
//			}
//		}
//	}
//
//	void removeVLine(final int laneNumber) {
//		for (final MyPlot p : plots) {
//			if (p.getNumber() == laneNumber) {
//				final Iterator<DataSeries> dataIter = p.getDataSeries().iterator();
//				while (dataIter.hasNext()) {
//					if (dataIter.next().getType() == DataSeries.VLINE) {
//						dataIter.remove();
//						break;
//					}
//				}
//				p.setSelected(false);
//				p.updatePlot();
//				break;
//			}
//		}
//	}

//	private void removeCustomPeaks(final MyPlot p) {
//		final Iterator<DataSeries> dataIter = p.getDataSeries().iterator();
//		while (dataIter.hasNext()) {
//			final DataSeries d = dataIter.next();
//			if (d.getType() == DataSeries.CUSTOMPEAKS) {
//				dataIter.remove();
//			}
//		}
//		p.updatePlot();
//	}

//	void updateCustomPeaks(final int lane, final ArrayList<Peak> pl) {
//		for (final MyPlot p : plots) {
//			if (p.getNumber() == lane) {
//				removeCustomPeaks(p);
//				RealVector xadd = new ArrayRealVector();
//				RealVector yadd = new ArrayRealVector();
//				for (final Peak fp : pl) {
//					xadd = xadd.append(fp.getMean());
//					yadd = yadd.append(fp.getNorm());
//				}
//				final DataSeries customPeaks = new DataSeries("CP ADD", lane,
//					DataSeries.CUSTOMPEAKS, xadd, yadd, Plotter.vLineAddPeakColor);
//				p.addDataSeries(customPeaks);
//				p.updatePlot();
//			}
//		}
//
//	}
	
	public ArrayList<DataSeries> getProfiles() {
		final ArrayList<DataSeries> profiles = new ArrayList<>();
		for (final ArrayList<DataSeries> p : plotsData) {
			for (DataSeries d : p) {
				if (d.getType() == DataSeries.PROFILE) profiles.add(d);
			}
		}
		return profiles;
	}
	
	void closePlot() {
		this.dispatchEvent(new WindowEvent(this, WindowEvent.WINDOW_CLOSING));
	}
}


class MyPlot implements Comparable<MyPlot> {

	private Plot plot;

	private final int number;

	private final String title;
	private final String xLabel;
	private final String yLabel;

	private ArrayList<DataSeries> data = new ArrayList<>();
	private double xmin;
	private double xmax;
	private double ymin;
	private double ymax;

	private boolean selected = false;

	private final Color unselectedBG = Color.white;

	private Color selectedBG = Color.white;

	public MyPlot(final String title, final String xLabel, final String yLabel) {
		this.title = title;
		this.xLabel = xLabel;
		this.yLabel = yLabel;
		this.number = Integer.parseInt(title.substring(5));
	}

	private void sortDataSeries() {
		Collections.sort(data);
		// Every time the DataSeries are sorted, the min/max are aslo updated
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

	public double[] getDataMinMax() {
		return new double[] { xmin, xmax, ymin, ymax };
	}

	public ArrayList<DataSeries> getDataSeries() {
		return data;
	}

	public int getNumber() {
		return number;
	}

	public Plot getPlot() {
		return plot;
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

	public void updatePlot() {
		final Plot newPlot = new Plot(title, xLabel, yLabel);
		newPlot.useTemplate(plot, Plot.COPY_LABELS + Plot.COPY_LEGEND +
			Plot.COPY_SIZE);
		if (plot != null) plot.dispose();
		plot = newPlot;
		if (selected) plot.setBackgroundColor(selectedBG);
		else plot.setBackgroundColor(unselectedBG);
		sortDataSeries();
		for (final DataSeries d : data) {
			plot.setColor(d.getColor());
			if (d.getType() == DataSeries.CUSTOMPEAKS) {
				plot.addPoints(d.getX().toArray(), d.getY().toArray(), Plot.CIRCLE);
			}
			else {
				plot.addPoints(d.getX().toArray(), d.getY().toArray(), Plot.LINE);
			}
		}
		plot.setLimits(xmin, xmax, ymin, ymax);
	}

	@Override
	public int compareTo(final MyPlot p) {
		return (number - p.getNumber());
	}
}

class DataSeries implements Comparable<DataSeries> {

	private final String name; // Name for Legend
	private final int lane; // Reference
	private final int type; // Type of function

	private RealVector x;
	private RealVector y;
	private Color color; // Plot color

	// Possible Types
	final static int VLINE = 0;

	final static int PROFILE = 1;
	final static int BACKGROUND = 2;
	final static int GAUSS_BG = 3;
	final static int FITTED = 4;
	final static int CUSTOMPEAKS = 5;

	public DataSeries(final String name, final int lane, final int type,
		final RealVector x, final RealVector y, final Color color)
	{
		this.name = name;
		this.lane = lane;
		this.type = type;
		this.x = x;
		this.y = y;
		this.color = color;
	}

	public DataSeries(final String name, final int lane, final int type,
		final RealVector x, final UnivariateFunction[] function,
		final Color color)
	{
		this.name = name;
		this.lane = lane;
		this.type = type;
		this.x = x;
		this.color = color;
		for (int i = 0; i < function.length; i++) {
			this.y = (i == 0) ? new ArrayRealVector(x.map(function[i])) : y.add(x.map(
				function[i]));
		}
	}

	public DataSeries(final String name, final int lane, final int type,
		final RealVector x, final UnivariateFunction function, final Color color)
	{
		this.name = name;
		this.lane = lane;
		this.type = type;
		this.x = x;
		this.color = color;
		this.y = new ArrayRealVector(x.map(function));
	}

	public String getName() {
		return name;
	}

	/**
	 * @return the lane
	 */
	public int getLane() {
		return lane;
	}

	public int getType() {
		return type;
	}

	public RealVector getX() {
		return x;
	}

	public RealVector getY() {
		return y;
	}
	public XYSeries getXYSeries() {
		XYSeries output = new XYSeries(this.getName());
		double[] xval = x.toArray();
		double[] yval = y.toArray();
		for (int i = 0; i < xval.length; i++) {
			output.add(xval[i], yval[i]);
		}
		return output;
	}
	public Color getColor() {
		return color;
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
