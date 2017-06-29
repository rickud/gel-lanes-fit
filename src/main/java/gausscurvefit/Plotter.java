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
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.Shape;
import java.awt.Stroke;
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
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.LegendItem;
import org.jfree.chart.LegendItemCollection;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.SeriesRenderingOrder;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.Range;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.Layer;
import org.jfree.ui.RectangleEdge;
import org.jfree.ui.RectangleInsets;
import org.jfree.ui.TextAnchor;
import org.scijava.Context;

import ij.ImagePlus;
import ij.gui.ProfilePlot;
import ij.gui.Roi;

class Plotter extends JFrame {
	
	private final ImagePlus imp;
	// private ArrayList<JFreeChart> plots;
	private ArrayList<ChartPanel> chartPanels;
	private ArrayList<DataSeries> plotsData;
	private ArrayList<VerticalMarker> verticalMarkers;
	
	private final JPanel chartsMainPanel;
	private final JScrollPane chartsScrollPane;

	private int rows = 2; // Number of plot Rows in display
	private final int cols = 2; // Number of plot Rows in display
	private int selected; // Which plot is highlighted
	private int plotMode = regMode;
	
	private static final int noRefPlot = -1; // no reference plot
	private int referencePlot = noRefPlot;
	
	// Colors are listed here for consistent, easy modification
	static final Color vMarkerRegColor = new Color(127, 127, 127);
	static final Color vMarkerAddPeakColor = new Color(0, 185, 19);
	static final Color vMarkerRemovePeakColor = new Color(185, 7, 0);
	static final Stroke vMarkerStroke = new BasicStroke();
	static final Color bMarkerColor = new Color(255, 128, 128);
	static final Stroke bMarkerStroke = new BasicStroke(1.25f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, new float[] {15.0f, 5.0f, 5.0f, 5.0f}, 0.0f);
	
	// Colors for DataSeries in plots
	static final Color profileColor = Color.BLACK;
	static final Color gaussColor = Color.RED;
	static final Color bgColor = Color.BLUE;
	static final Color fittedColor = new Color(255, 153, 0);
	static final Stroke dataStroke = new BasicStroke(2.0f);
	// Background Color of selected plot
	static final int selectedNone = -1;
	static final int regMode = 0;
	static final int addMode = 1;
	static final int remMode = 2;
	static final int refMode = 2;

	static final Color plotSelColor = new Color(240, 240, 240);
	static final Color plotUnselColor = Color.WHITE;
	static final Color plotAddSelColor = new Color(192, 255, 185);
	static final Color plotRemoveSelColor = new Color(255, 188, 185);

	public Plotter(final Context context, final ImagePlus imp, final ArrayList<Roi> rois,
		final double screenHeight, final double screenWidth)
	{
		context.inject(this);
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
		this.verticalMarkers = new ArrayList<>();
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

	public void addVerticalMarkers(ArrayList<VerticalMarker> vMarkers) {
		verticalMarkers.addAll(vMarkers);
	}

	void setVLine(final int ln, final double x) {
		Color vMarkerColor;
		if (plotMode == Plotter.addMode) {
			vMarkerColor = vMarkerAddPeakColor;
		} else if (plotMode == Plotter.remMode) {
			vMarkerColor = vMarkerRemovePeakColor;
		} else { // (plotMode == Plotter.regMode)
			vMarkerColor = vMarkerRegColor;
		} 
		for (final ChartPanel p : chartPanels) {
			if (p.getChart().getTitle().getText().equals("Lane " + ln)) {
				boolean found = false;
				for (VerticalMarker m : verticalMarkers) {
					if (m.getLane() == ln && m.getType() == VerticalMarker.VMARK) {
						found = true;
						m.setName(String.format("%1$.1f", x));
						m.setValue(x);
						m.setPaint(vMarkerColor);
					}
				}
				if (!found) {					
					verticalMarkers.add(new VerticalMarker(String.format("%1$.1f", x), ln, VerticalMarker.VMARK, x, vMarkerColor, vMarkerStroke));
				}
				updatePlot(ln);
			}
		}
	}

	void removeVLine(final int ln) {
		final Iterator<VerticalMarker> markIter = verticalMarkers.iterator();
		while (markIter.hasNext()) {
			final VerticalMarker m = markIter.next();
			if (m.getLane() == ln && m.getType() == VerticalMarker.VMARK) {
				markIter.remove();
				for (ChartPanel c : chartPanels) {
					if (Integer.parseInt(c.getChart().getTitle().getText().substring(5)) == ln) {
						Iterator<?> aIter = c.getChart().getXYPlot().getAnnotations().iterator();
						while (aIter.hasNext()) {
							XYTextAnnotation a = (XYTextAnnotation) aIter.next();
							if (a.getText().equals(String.format("%1$.1f", m.getValue())))
								aIter.remove();
						}
					}
				}
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
	
	public ArrayList<DataSeries> getPlotsData() {
		return plotsData;
	}
	
	public ArrayList<VerticalMarker> getPlotsVerticalMarkers() {
		return verticalMarkers;
	}

	double[] getPlotLimits(int ln) {
		double[] limits = new double[4];
		for (ChartPanel p : chartPanels) {
			if (p.getChart().getTitle().getText().equals("Lane " + ln)) {
				XYPlot pl = p.getChart().getXYPlot();
				Range xrange = pl.getDataRange(pl.getDomainAxis());
				Range yrange = pl.getDataRange(pl.getRangeAxis());
				limits[0] = xrange.getLowerBound();
				limits[1] = xrange.getUpperBound();
				limits[2] = yrange.getLowerBound();
				limits[3] = yrange.getUpperBound();
			}
		}
		return limits;
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

	public void setReferencePlot(final int ref, ArrayList<Peak> peaks) {
		referencePlot  = ref;
		if (ref != noRefPlot) {
			for (DataSeries d : plotsData) {
				int ln = d.getLane();
				if (ln != ref) {
					ArrayList<VerticalMarker> bands = new ArrayList<>();
					int pk = 0;
					for (Peak p : peaks) {
						double x = p.getMean();
						VerticalMarker m = new VerticalMarker("Band " + ++pk, ln, VerticalMarker.BMARK, x, Plotter.bMarkerColor, Plotter.bMarkerStroke);
						bands.add(m);
					}
					verticalMarkers.addAll(bands);
				}
			}
		}
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
					DataSeries.CUSTOMPEAKS, xadd, yadd, Plotter.vMarkerAddPeakColor);
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
		profile.setKey(String.format("%d", DataSeries.PROFILE));
		dataset.addSeries(profile);
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
				xLabel, yLabel, dataset, PlotOrientation.VERTICAL, true, true, false);
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
		final double min = 0.95 * profile.getMinY();
		final double max = 1.05 * profile.getMaxY();
		thePlot.getRangeAxis().setLowerBound(min);
		thePlot.getRangeAxis().setUpperBound(max);

		updatePlot(Integer.parseInt(roi.getName().substring(5)));
	}

	public void updatePlot(final int ln) {
		Color plotBGColor  = plotSelColor;
		Color vMarkerColor = vMarkerRegColor;
		if (plotMode == Plotter.addMode) {
			plotBGColor = plotAddSelColor;
			vMarkerColor = vMarkerAddPeakColor;
		}
		else if (plotMode == Plotter.remMode) {
			plotBGColor = plotRemoveSelColor;
			vMarkerColor = vMarkerRemovePeakColor;
		}
		else { // (plotMode == Plotter.regMode)
			plotBGColor = plotSelColor;
			vMarkerColor = vMarkerRegColor;
		}

		final XYSeriesCollection dataset = new XYSeriesCollection();
		Collections.sort(plotsData);

		for (final ChartPanel p : chartPanels) {
			final int plotNumber = Integer.parseInt(p.getChart().getTitle().getText()
				.substring(5));
			if (plotNumber == ln) {
				int gcount = 0;				
				final JFreeChart c = p.getChart();
				final XYPlot pl = c.getXYPlot(); 
				
				if (plotNumber == selected) {
					pl.setBackgroundPaint(plotBGColor);
					c.setBackgroundPaint(plotBGColor);
				}
				else {
					pl.setBackgroundPaint(plotUnselColor);
					c.setBackgroundPaint(plotUnselColor);
				}

				// Clear Markers and Annotations
				if (c.getXYPlot().getDomainMarkers(Layer.BACKGROUND) != null)
					c.getXYPlot().clearDomainMarkers();
				if (c.getXYPlot().getAnnotations() != null) 
					c.getXYPlot().clearAnnotations();
				
				// Plot vertical markers
				for (VerticalMarker m : verticalMarkers) {
					if (m.getLane() == ln) {
						double height = c.getXYPlot().getRangeAxis().getUpperBound();
						double offset = 0;
						XYTextAnnotation label = new XYTextAnnotation(m.getName(), m.getValue() - offset, height);
						if (m.getType() == VerticalMarker.VMARK) {
							label.setPaint(vMarkerColor);
						} else if (m.getType() == VerticalMarker.BMARK) {
							label.setPaint(bMarkerColor);
						}
						label.setFont(new Font("Sans Serif", Font.PLAIN, 12));
						label.setRotationAnchor(TextAnchor.BOTTOM_RIGHT);
						label.setTextAnchor(TextAnchor.TOP_RIGHT);
						label.setRotationAngle(-Math.PI / 2);
						
						c.getXYPlot().addAnnotation(label);
						c.getXYPlot().addDomainMarker(m, Layer.BACKGROUND);
					}
				}
				
				// Plot the data series
				LegendItems legendItems = new LegendItems();
				for (final DataSeries d : plotsData) {
					if (d.getLane() == ln) {
						int k = d.getType();
						String key = "";
						if (k == DataSeries.GAUSS_BG) {
							key = String.format("%d", k + gcount);
							gcount++;
						} else {
							key = String.format("%d", k);
						}
						d.setKey(key);
						dataset.addSeries(d);
						pl.setSeriesRenderingOrder(SeriesRenderingOrder.FORWARD);
						final int seriesIdx = dataset.getSeriesIndex(key);
						final XYLineAndShapeRenderer renderer = (XYLineAndShapeRenderer) pl.getRenderer();
						renderer.setSeriesShapesVisible(seriesIdx, false);
						renderer.setSeriesLinesVisible(seriesIdx, true);
						
						if (k == DataSeries.PROFILE) {
							renderer.setSeriesPaint(seriesIdx, profileColor);
							renderer.setSeriesStroke(seriesIdx, dataStroke);
							LegendItem li = new LegendItem("Lane " + ln + " Profile");
							li.setFillPaint(d.getColor());
							legendItems.add(li);
						}
						if (referencePlot == noRefPlot || (referencePlot != noRefPlot && ln == referencePlot)) {
							if (k == DataSeries.BACKGROUND) {
								renderer.setSeriesPaint(seriesIdx, bgColor);
								renderer.setSeriesStroke(seriesIdx, dataStroke);
								LegendItem li = new LegendItem("Background Polynomial");
								li.setFillPaint(d.getColor());
								legendItems.add(li);
							}
							if (k == DataSeries.GAUSS_BG) {
								renderer.setSeriesPaint(seriesIdx, gaussColor);
								renderer.setSeriesStroke(seriesIdx, dataStroke);
								LegendItem li = new LegendItem("Gaussian Peaks");
								if (!legendItems.contains(li)) {
									li.setFillPaint(d.getColor());
									legendItems.add(li);
								}
							}
							if (k == DataSeries.FITTED) {
								renderer.setSeriesPaint(seriesIdx, fittedColor);
								renderer.setSeriesStroke(seriesIdx, dataStroke);
								LegendItem li = new LegendItem("Fitted Curve");
								li.setFillPaint(d.getColor());
								legendItems.add(li);
							}
							if (k == DataSeries.CUSTOMPEAKS) {
								renderer.setSeriesPaint(seriesIdx, vMarkerAddPeakColor);
								final Shape dot = new Ellipse2D.Double(0, 0, 6, 6);
								renderer.setSeriesShape(seriesIdx, dot);
								renderer.setSeriesShapesVisible(seriesIdx, true);
								renderer.setSeriesLinesVisible(seriesIdx, false);
								LegendItem li = new LegendItem("Custom Peaks");
								li.setFillPaint(d.getColor());
								legendItems.add(li);
							}
						}
					}
				}
				c.getXYPlot().setDataset(dataset);
				pl.setFixedLegendItems(legendItems);
				c.getLegend().setPosition(RectangleEdge.RIGHT);
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

	public void resetData() {
		plotsData = new ArrayList<>();
		chartPanels = new ArrayList<>();
		chartsMainPanel.removeAll();
	}
	
	class LegendItems extends LegendItemCollection {
		public LegendItems() {					
			super();
		}
		public boolean contains(LegendItem li) {
			for (int i = 0; i < this.getItemCount(); i++) {
				if (this.get(i).getLabel().equals(li.getLabel()))
					return true;
			}
			return false;
		}
	}
}


class VerticalMarker extends ValueMarker {
	// Possible types
	final static int VMARK = 0;	// Vertical Position
	final static int BMARK = 1;	// Molecular Weight bands
	
	
	private final int lane; 		// Reference GEL LANE
	private final int type;
	private String name;				// Name for reference

	public VerticalMarker(final String name, final int lane, final int type,
		final double x, final Color color, final Stroke stroke)
	{
		super(x, color, stroke);
		this.name = name;
		this.lane = lane;
		this.type = type;
	}
	
	void setName(String name) {
		this.name = name;
	}

	int getLane() {
		return lane;
	}
	
	String getName() {
		return name;
	}
	
	int getType() {
		return type;
	}
}


class DataSeries extends XYSeries implements Comparable<DataSeries> {

	private final String name; // Name for Legend
	private final int lane; // Reference
	private final int type; // Type of function

	private RealVector x;
	private RealVector y;
	private Color color; // Plot color

	// Possible Types
	final static int PROFILE = 0;
	final static int BACKGROUND = 1;
	final static int GAUSS_BG = 2;
	final static int FITTED = 400;
	final static int CUSTOMPEAKS = 401;

	public DataSeries(final String name, final int lane, final int type,
		final RealVector x, final RealVector y, final Color color)
	{
		super(type);
		this.name = name;
		this.lane = lane;
		this.type = type;
		this.color = color;
		this.x = x;
		this.y = y;
		addData();
	}

	public DataSeries(final String name, final int lane, final int type,
		final RealVector x, final UnivariateFunction[] function,
		final Color color)
	{
		super(type);
		this.name = name;
		this.lane = lane;
		this.type = type;
		this.x = x;
		this.color = color;
		for (int i = 0; i < function.length; i++) {
			this.y = (i == 0) ? new ArrayRealVector(x.map(function[i])) : y.add(x.map(
				function[i]));
		}
		addData();
	}

	public DataSeries(final String name, final int lane, final int type,
		final RealVector x, final UnivariateFunction function, final Color color)
	{
		super(type);
		this.name = name;
		this.lane = lane;
		this.type = type;
		this.x = x;
		this.color = color;
		this.y = new ArrayRealVector(x.map(function));
		addData();
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

	void addData() {
		RealMatrix series = new Array2DRowRealMatrix(new double[][] {x.toArray(), y.toArray()});
		for (int r = 0; r < series.getColumnDimension(); r++) {
			add(series.getEntry(0, r), series.getEntry(1, r));
		}
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

	public void setX(final RealVector x) {
		this.x = x;
	}

	public void setY(final RealVector y) {
		this.y = y;
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
