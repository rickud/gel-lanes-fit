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

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.GridLayout;
import java.awt.Shape;
import java.awt.event.WindowEvent;
import java.awt.geom.Ellipse2D;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.ScrollPaneConstants;
import javax.swing.WindowConstants;
import javax.swing.border.EmptyBorder;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.SeriesRenderingOrder;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.RectangleInsets;

import ij.ImagePlus;
import ij.gui.Plot;
import ij.gui.ProfilePlot;
import ij.gui.Roi;

class Plotter extends JFrame {

	private final ImagePlus imp;
	// private ArrayList<JFreeChart> plots;
	private ArrayList<ChartPanel> chartPanels;
	private ArrayList<DataSeries> plotsData;
	private final JPanel chartsMainPanel;
	private final JScrollPane chartsScrollPane;

	private int rows = 2; // Number of plot Rows in display
	private final int cols = 2; // Number of plot Rows in display
	private int selected; // Which plot is highlighted
	private int plotMode = regMode;

	// Colors are listed here for consistent, easy modification
	static final Color vLineRegColor = new Color(127, 127, 127);
	static final Color vLineAddPeakColor = new Color(0, 185, 19);
	static final Color vLineRemovePeakColor = new Color(185, 7, 0);
	// Colors for DataSeries in plots
	static final Color profileColor = Color.BLACK;
	static final Color gaussColor = Color.RED;
	static final Color bgColor = Color.BLUE;
	static final Color fittedColor = new Color(255, 153, 0);
	// Background Color of selected plot
	static final int selectedNone = -1;
	static final int regMode = 0;
	static final int addMode = 1;
	static final int remMode = 2;

	static final Color plotSelColor = new Color(240, 240, 240);
	static final Color plotUnselColor = Color.WHITE;
	static final Color plotAddSelColor = new Color(192, 255, 185);
	static final Color plotRemoveSelColor = new Color(255, 188, 185);

	public Plotter(final ImagePlus imp, final ArrayList<Roi> rois,
		final double screenHeight, final double screenWidth)
	{
		setDefaultCloseOperation(WindowConstants.HIDE_ON_CLOSE);
		setBounds((int) screenWidth / 2, 0, (int) screenWidth / 2,
			(int) screenHeight);

		chartsMainPanel = new JPanel();
		chartsMainPanel.setBorder(new EmptyBorder(5, 5, 5, 5));
		chartsMainPanel.setLayout(new GridLayout(rows, cols));
		selected = selectedNone;

		this.imp = imp;
		this.chartPanels = new ArrayList<>();
		this.plotsData = new ArrayList<>();
		this.setTitle("Profiles of " + imp.getShortTitle());

		for (final Roi r : rois) {
			updateProfile(r);
		}
		for (final ChartPanel p : chartPanels) {
			chartsMainPanel.add(p);
		}
		chartsScrollPane = new JScrollPane(chartsMainPanel);
		chartsMainPanel.setBounds(chartsScrollPane.getVisibleRect());
		chartsScrollPane.setHorizontalScrollBarPolicy(
			ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
		this.getContentPane().add(chartsScrollPane, BorderLayout.CENTER);
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

	public void addDataSeries(final ArrayList<DataSeries> d) {
		plotsData.addAll(d);
	}

	void setVLine(final int ln, final double x) {
		for (final ChartPanel p : chartPanels) {
			if (p.getChart().getTitle().getText().equals("Lane " + ln)) {
				boolean vline = false;
				final XYSeriesCollection dataset = (XYSeriesCollection) p.getChart()
					.getXYPlot().getDataset();
				final XYSeries profile = dataset.getSeries(String.format("%d",
					DataSeries.PROFILE));
				final double min = 0.95 * profile.getMinY();
				final double max = 1.05 * profile.getMaxY();
				final RealVector xvals = new ArrayRealVector(new double[] { x, x });
				final RealVector yvals = new ArrayRealVector(new double[] { min, max });
				final DataSeries newVLine = new DataSeries("", ln, DataSeries.VLINE,
					xvals, yvals, vLineRegColor);
				for (final DataSeries d : plotsData) {
					if (d.getLane() == ln && d.getType() == DataSeries.VLINE) {
						vline = true;
						plotsData.set(plotsData.indexOf(d), newVLine);
					}
				}
				if (!vline) {
					plotsData.add(newVLine);
				}
				updatePlot(ln);
			}
		}
	}

	void removeVLine(final int ln) {
		final Iterator<DataSeries> dataIter = plotsData.iterator();
		while (dataIter.hasNext()) {
			final DataSeries d = dataIter.next();
			if (d.getLane() == ln && d.getType() == DataSeries.VLINE) {
				dataIter.remove();
			}
		}
		selected = selectedNone;
		updatePlot(ln);
	}

	public ArrayList<DataSeries> getProfiles() {
		final ArrayList<DataSeries> profiles = new ArrayList<>();
		for (final DataSeries d : plotsData) {
			if (d.getType() == DataSeries.PROFILE) profiles.add(d);
		}
		return profiles;
	}

	int getPlotMode() {
		return this.plotMode;
	}

	void setPlotMode(final int plotMode) {
		this.plotMode = plotMode;
	}

	public void setSelected(final int selected) {
		this.selected = selected;
	}

	public void setMode(final int plotMode) {
		this.plotMode = plotMode;
	}

	private void sortDataSeries() {
		Collections.sort(plotsData);
		int nullIndex = 0;
		final Iterator<DataSeries> dataIter = plotsData.iterator();
		while (dataIter.hasNext()) {
			final DataSeries d = dataIter.next();
			if (d == null) break;
			nullIndex++;
		}
		plotsData = new ArrayList<>(plotsData.subList(0, nullIndex));
	}

	void updateCustomPeaks(final int lane, final ArrayList<Peak> pl) {
		for (final ChartPanel p : chartPanels) {
			final int plotNumber = Integer.parseInt(p.getChart().getTitle().getText()
				.substring(5));
			if (plotNumber == lane) {
				final Iterator<DataSeries> dataIter = plotsData.iterator();
				while (dataIter.hasNext()) {
					final DataSeries d = dataIter.next();
					if (d.getLane() == lane && d.getType() == DataSeries.CUSTOMPEAKS)
						dataIter.remove();
				}
				RealVector xadd = new ArrayRealVector();
				RealVector yadd = new ArrayRealVector();
				for (final Peak fp : pl) {
					xadd = xadd.append(fp.getMean());
					yadd = yadd.append(fp.getNorm());
				}
				final DataSeries customPeaks = new DataSeries("Custom Peaks", lane,
					DataSeries.CUSTOMPEAKS, xadd, yadd, Plotter.vLineAddPeakColor);
				plotsData.add(customPeaks);
				updatePlot(lane);
			}
		}

	}

	void updateProfile(final Roi roi) {
		final DataSeries profile = getLaneProfile(roi);
		// Assume plotsData, chartsMainPanel was reset
		plotsData.add(profile);
		final String xLabel = "Distance (px)";
		final String yLabel = "Grayscale Value";
		final XYSeriesCollection dataset = new XYSeriesCollection();
		XYPlot thePlot = new XYPlot();
		final XYSeries profileSeries = profile.getXYSeries();
		profileSeries.setKey(String.format("%d", DataSeries.PROFILE));
		dataset.addSeries(profile.getXYSeries());
		boolean found = false;
		for (final ChartPanel c : chartPanels) {
			thePlot = c.getChart().getXYPlot();
			if (c.getChart().getTitle().getText().equals(roi.getName())) {
				found = true;
				thePlot.setDataset(dataset);
			}
		}
		if (!found) {
			final JFreeChart newChart = ChartFactory.createXYLineChart(roi.getName(),
				xLabel, yLabel, dataset, PlotOrientation.VERTICAL, false, true, false);
			final XYPlot newPlot = newChart.getXYPlot();
			newChart.getTitle().setMargin(new RectangleInsets(15, 5, 15, 5));
			newChart.getTitle().setPaint(Color.BLACK);

			newPlot.setDomainGridlinePaint(Color.DARK_GRAY);
			newPlot.setRangeGridlinePaint(Color.DARK_GRAY);
			// newPlot.setAxisOffset(new RectangleInsets(0,0,0,0));
			thePlot = newPlot;
			final ChartPanel chartPanel = new ChartPanel(newChart);
			chartPanels.add(chartPanel);
			chartsMainPanel.add(chartPanel);
			if (chartPanels.size() % cols == 0) rows = chartPanels.size() / cols;
			else rows = chartPanels.size() / cols + 1;
			chartsMainPanel.setLayout(new GridLayout(rows, cols));
			chartsMainPanel.validate();
		}
		selected = selectedNone;

		thePlot.getDomainAxis().setLowerMargin(0);
		thePlot.getDomainAxis().setUpperMargin(0);
		thePlot.getRangeAxis().setLowerMargin(0);
		thePlot.getRangeAxis().setUpperMargin(0);
		final double min = 0.95 * profile.getY().getMinValue();
		final double max = 1.05 * profile.getY().getMaxValue();
//		Range r = new Range(min, max);
		thePlot.getRangeAxis().setLowerBound(min);
		thePlot.getRangeAxis().setUpperBound(max);

		updatePlot(Integer.parseInt(roi.getName().substring(5)));
	}

	public void updatePlot(final int ln) {

		Color plotBGColor = plotSelColor;
		final Color pkColor = vLineAddPeakColor;
		Color vLineColor;

		if (plotMode == Plotter.addMode) {
			plotBGColor = plotAddSelColor;
			vLineColor = vLineAddPeakColor;
		}
		else if (plotMode == Plotter.remMode) {
			plotBGColor = plotRemoveSelColor;
			vLineColor = vLineRemovePeakColor;
		}
		else { // (plotMode == Plotter.regMode)
			plotBGColor = plotSelColor;
			vLineColor = vLineRegColor;
		}

		final XYSeriesCollection dataset = new XYSeriesCollection();
		// sortDataSeries();
		int PKCount = 0;
		final ArrayList<String> keys = new ArrayList<>();
		Collections.sort(plotsData);
		for (final DataSeries d : plotsData) {
			if (d.getLane() == ln) {
				final XYSeries dataSeries = d.getXYSeries();
				int k = d.getType();
				if (k == DataSeries.GAUSS_BG) {
					k += PKCount;
					PKCount++;
				}
				keys.add(String.format("%d", k));
				dataSeries.setKey(String.format("%d", k));
				dataset.addSeries(dataSeries);
			}
		}
		for (final ChartPanel p : chartPanels) {
			final int plotNumber = Integer.parseInt(p.getChart().getTitle().getText()
				.substring(5));
			if (plotNumber == ln) {
				final JFreeChart c = p.getChart();
				c.getXYPlot().setDataset(dataset);
				if (plotNumber == selected) {
					c.getXYPlot().setBackgroundPaint(plotBGColor);
					c.setBackgroundPaint(plotBGColor);
				}
				else {
					c.getXYPlot().setBackgroundPaint(plotUnselColor);
					c.setBackgroundPaint(plotUnselColor);
				}

				p.getChart().getXYPlot().setSeriesRenderingOrder(
					SeriesRenderingOrder.FORWARD);
				final XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) p
					.getChart().getXYPlot().getRenderer();
				PKCount = 0;

				for (final String k : keys) {
					final int seriesIdx = dataset.getSeriesIndex(k);
					final int seriesType = Integer.parseInt(k);
					if (seriesType == DataSeries.VLINE) {
						renderer.setSeriesPaint(seriesIdx, vLineColor);
					}
					if (seriesType == DataSeries.PROFILE) {
						renderer.setSeriesPaint(seriesIdx, profileColor);
					}
					if (seriesType == DataSeries.BACKGROUND) {
						renderer.setSeriesPaint(seriesIdx, bgColor);
					}
					if (seriesType == DataSeries.GAUSS_BG + PKCount) {
						renderer.setSeriesPaint(seriesIdx, gaussColor);
						PKCount++;
					}
					if (seriesType == DataSeries.FITTED) {
						renderer.setSeriesPaint(seriesIdx, fittedColor);
					}
					if (seriesType == DataSeries.CUSTOMPEAKS) {
						renderer.setSeriesPaint(seriesIdx, pkColor);
						final Shape dot = new Ellipse2D.Double(0, 0, 6, 6);
						renderer.setSeriesShape(seriesIdx, dot);
						renderer.setSeriesShapesVisible(seriesIdx, true);
						renderer.setSeriesLinesVisible(seriesIdx, false);
					}

					// Update the Legend
//					LegendItemCollection legendItemsOld = c.getXYPlot().getLegendItems();
//					final LegendItemCollection legendItemsNew = new LegendItemCollection();
//
//					for(int i = 0; i< legendItemsOld.getItemCount(); i++){
//						LegendItem li = legendItemsOld.get(i);
//						if (li.getLabel().equals(String.format("%d", DataSeries.PROFILE))){
//							LegendItem newLi = new LegendItem("Lane " + ln + " Profile");
//							newLi.setSeriesIndex(li.getSeriesIndex());
//							li = newLi;
//						} else if (li.getSeriesKey().equals(String.format("%d", DataSeries.BACKGROUND))) {
//							li.setDescription("Background Polynomial");
//							legendItemsNew.add(li);
//						} else if (li.getSeriesKey().equals(String.format("%d", DataSeries.GAUSS_BG))) {
//							li.setDescription("Gaussian Peaks");
//							legendItemsNew.add(li);
//						} else if (li.getSeriesKey().equals(String.format("%d", DataSeries.FITTED))) {
//							li.setDescription("Fitted Curve");
//							legendItemsNew.add(li);
//						}
//					}
//					LegendItemSource source = new LegendItemSource() {
//						LegendItemCollection lic = new LegendItemCollection();
//						{lic.addAll(legendItemsNew);}
//						@Override
//						public LegendItemCollection getLegendItems() {
//							return lic;
//						}
//					};
//					c.addLegend(new LegendTitle(source));
//					c.getLegend().setPosition(RectangleEdge.RIGHT);
//					c.getLegend().setVisible(true);
				}
			}
		}
	}

	void closePlot() {
		this.dispatchEvent(new WindowEvent(this, WindowEvent.WINDOW_CLOSING));
	}

	public void removeFit() {
		final Iterator<DataSeries> dataIter = plotsData.iterator();
		while (dataIter.hasNext()) {
			final DataSeries d = dataIter.next();
			if (d.getType() != DataSeries.PROFILE && d
				.getType() != DataSeries.CUSTOMPEAKS) dataIter.remove();
		}
	}

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

	public void resetData() {
		plotsData = new ArrayList<>();
		chartPanels = new ArrayList<>();
		chartsMainPanel.removeAll();
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
	final static int FITTED = 400;
	final static int CUSTOMPEAKS = 401;

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
		final XYSeries output = new XYSeries(this.getType());
		final double[] xval = x.toArray();
		final double[] yval = y.toArray();
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
		final int t = type - d.getType();
		final int l = lane - d.getLane();

		if (l == 0) {
			return t;
		}
		return l;
	}
}
