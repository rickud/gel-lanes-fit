/*
 * Gel Lanes Fit - Plotter.java
 * Author: Rick Ziraldo, 2017
 * The University of Texas at Dallas, Richardson, TX
 *
 * Licensed under the GNU Affero General Public License v3.0; see LICENSE.
 * Source: https://github.com/rickud/gel-lanes-fit
 */

package gellanesfit;

import com.itextpdf.awt.PdfGraphics2D;
import com.itextpdf.text.DocumentException;
import com.itextpdf.text.PageSize;
import com.itextpdf.text.Rectangle;
import com.itextpdf.text.pdf.PdfContentByte;
import com.itextpdf.text.pdf.PdfTemplate;
import com.itextpdf.text.pdf.PdfWriter;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.Paint;
import java.awt.RenderingHints;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.event.WindowEvent;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.atomic.AtomicLong;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.WindowConstants;
import javax.swing.border.EmptyBorder;

import org.apache.commons.io.FileUtils;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartMouseEvent;
import org.jfree.chart.ChartMouseListener;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtils;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.LegendItem;
import org.jfree.chart.LegendItemCollection;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.event.ChartProgressEvent;
import org.jfree.chart.labels.StandardXYToolTipGenerator;
import org.jfree.chart.labels.XYToolTipGenerator;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.SeriesRenderingOrder;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.text.TextUtils;
import org.jfree.chart.ui.Layer;
import org.jfree.chart.ui.RectangleEdge;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.jfree.data.general.DatasetUtils;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.graphics2d.svg.SVGGraphics2D;
import org.jfree.graphics2d.svg.SVGUtils;
import org.scijava.Context;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.ProfilePlot;
import ij.gui.Roi;

/**
 * The Profiles window: one chart per lane, in tabs of four, with the lane's
 * profile, the fitted background, peaks and fit, the custom peaks, and the
 * ladder bands as labelled vertical markers. Clicks in a chart add or remove
 * custom peaks in Edit Custom Peaks mode. Also saves the charts as PNG, PDF and
 * SVG.
 */
class Plotter extends JFrame implements ChartMouseListener {

	private static final long serialVersionUID = 1L;

	@Parameter
	private LogService log;

	private final double SW = IJ.getScreenSize().getWidth();
	private final double SH = IJ.getScreenSize().getHeight();

	private int referencePlot = MainDialog.noLadderLane;

	// Colors are listed here for consistent, easy modification
	private static final Color vMarkerRegColor = new Color(127, 127, 127);
	static final Color vMarkerEditPeakColor = new Color(0, 185, 19);
	private static final Stroke vMarkerStroke = new BasicStroke();
	private static final Color bMarkerColor = new Color(255, 128, 128);
	private static final Stroke bMarkerStroke = new BasicStroke(1.25f,
		BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER, 10.0f, new float[] { 15.0f,
			5.0f, 5.0f, 5.0f }, 0.0f);

	// Colors for DataSeries in plots
	private static final Color profileColor = Color.BLACK;
	static final Color gaussColor = Color.RED;
	static final Color bgColor = Color.BLUE;
	static final Color fittedColor = new Color(255, 153, 0);
	private static final Stroke dataStroke = new BasicStroke(2.0f);

	// Background Color of selected plot
	static final int regMode = 0;
	static final int editPeaksMode = 1;

	private static final Color plotSelColor = new Color(240, 240, 240);
	private static final Color plotUnselColor = Color.WHITE;
	private static final Color plotAddSelColor = new Color(192, 255, 185);
	private static final Color plotRefColor = new Color(206, 217, 255);

	private final ImagePlus imp;
	private final List<ChartPanel> chartPanels;
	private final List<DataSeries> plotsData;
	private final List<Integer> plotNumbers;
	// Lanes whose axis ranges are set; kept on redraws until the profiles change
	private final Set<Integer> rangeSet = new HashSet<>();
	// Range each chart was fitted to, to tell a reset view from a user zoom
	private final Map<JFreeChart, Range> fittedRanges = new HashMap<>();
	// Plot height each fit used, and charts whose data changed since their fit
	private final Map<JFreeChart, Double> fittedHeights = new HashMap<>();
	private final Set<JFreeChart> labelFitPending = new HashSet<>();
	private static final Font labelFont = new Font("Sans Serif", Font.PLAIN, 14);
	private static final int labelGap = 4; // px between a label and the curves
	private List<VerticalMarker> verticalMarkers;

	private List<JPanel> chartTabs;
	private final JTabbedPane chartsTabbedPane;

	private final int rows = 2; // Number of plot Rows in display
	private final int cols = 2; // Number of plot Rows in display
	private int selected = MainDialog.noLaneSelected; // Which plot is highlighted
	private int plotMode = regMode;

	public Plotter(final Context context, final ImagePlus imp) {
		context.inject(this);
		setDefaultCloseOperation(WindowConstants.HIDE_ON_CLOSE);
		setBounds((int) (SW * 0.3), 0, (int) (SW * 0.7), (int) (SH * 0.9));
		this.setTitle("Profiles of " + imp.getShortTitle());
		this.imp = imp;
		chartPanels = new ArrayList<>();
		chartTabs = new ArrayList<>();
		plotsData = new ArrayList<>();
		plotNumbers = new ArrayList<>();
		verticalMarkers = new ArrayList<>();
		chartsTabbedPane = new JTabbedPane();
		this.getContentPane().add(chartsTabbedPane, BorderLayout.CENTER);
		this.setVisible(true);
	}

	/**
	 * Profile data from Roi, ready for fitting (null if not possible)
	 *
	 * @param roi
	 */
	/**
	 * A lane's profile: the intensity averaged across the lane's width, row by
	 * row. x is the row's distance from the top of the image, in pixels.
	 *
	 * @param imp the gel image
	 * @param roi the lane, named "Lane n"
	 * @return the profile, or null if the lane is less than 2 pixels tall
	 */
	static DataSeries laneProfile(final ImagePlus imp, final Roi roi) {
		imp.setRoi(roi);
		final RealVector profile = new ArrayRealVector(new ProfilePlot(imp, true)
			.getProfile());
		imp.killRoi();
		if (profile.getDimension() < 2) return null;
		final String name = roi.getName();
		final int lane = Integer.parseInt(name.substring(5));
		final double y0 = roi.getBounds().getMinY();
		final double[] y = new double[profile.getDimension()];
		for (int p = 0; p < y.length; p++)
			y[p] = y0 + p;
		return new DataSeries(name, lane, DataSeries.PROFILE, new ArrayRealVector(
			y), profile, Plotter.profileColor);
	}

	public void addDataSeries(final DataSeries data) {
		plotsData.add(data);
	}

	public void addDataSeries(final List<DataSeries> data) {
		plotsData.addAll(data);
	}

	public void addVerticalMarkers(final List<Peak> peaks) {
		for (final DataSeries d : plotsData) {
			final int ln = d.getLane();
			final ArrayList<VerticalMarker> bands = new ArrayList<>();
			for (final Peak p : peaks) {
				final double x = p.getMean();
				final VerticalMarker m = new VerticalMarker(p.getName(), ln,
					VerticalMarker.BMARK, x, Plotter.bMarkerColor, Plotter.bMarkerStroke);
				bands.add(m);
			}
			verticalMarkers.addAll(bands);
		}
	}

	public void removeVerticalMarkers() {
		verticalMarkers = new ArrayList<>();
		for (final int i : plotNumbers)
			updatePlot(i);
	}

	void setVLine(final int ln, final double x) {
		Color vMarkerColor = vMarkerRegColor;
		if (plotMode == Plotter.editPeaksMode) {
			vMarkerColor = vMarkerEditPeakColor;
		}
		for (final ChartPanel p : chartPanels) {
			if (p.getChart().getTitle().getText().equals("Lane " + ln)) {
				boolean found = false;
				for (final VerticalMarker m : verticalMarkers) {
					if (m.getLane() == ln && m.getType() == VerticalMarker.VMARK) {
						found = true;
						m.setName(String.format("%1$.1f", x));
						m.setValue(x);
						m.setPaint(vMarkerColor);
					}
				}
				if (!found) {
					verticalMarkers.add(new VerticalMarker(String.format("%1$.1f", x), ln,
						VerticalMarker.VMARK, x, vMarkerColor, vMarkerStroke));
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
				for (final ChartPanel c : chartPanels) {
					if (Integer.parseInt(c.getChart().getTitle().getText().substring(
						5)) == ln)
					{
						final Iterator<?> aIter = c.getChart().getXYPlot().getAnnotations()
							.iterator();
						while (aIter.hasNext()) {
							final XYTextAnnotation a = (XYTextAnnotation) aIter.next();
							if (a.getText().equals(String.format("%1$.1f", m.getValue())))
								aIter.remove();
						}
					}
				}
			}
		}
		selected = MainDialog.noLaneSelected;
		updatePlot(ln);
	}

	public ArrayList<DataSeries> getProfiles() {
		final ArrayList<DataSeries> profiles = new ArrayList<>();
		for (final DataSeries d : plotsData) {
			if (d.getType() == DataSeries.PROFILE) profiles.add(d);
		}
		return profiles;
	}

	public DataSeries getPlotsCustomPeaks(final int lane) {
		for (final DataSeries d : plotsData) {
			if (d.getLane() == lane && d.getType() == DataSeries.CUSTOMPEAKS)
				return d;
		}
		return null;
	}

	public List<DataSeries> getPlotsData() {
		return plotsData;
	}

	public List<VerticalMarker> getPlotsVerticalMarkers() {
		return verticalMarkers;
	}

	void setPlotMode(final int plotMode) {
		this.plotMode = plotMode;
	}

	public void setSelected(final int selected) {
		this.selected = selected;
		if (selected != MainDialog.noLaneSelected) {
			final int tab = (selected - 1) / (rows * cols);
			if (chartsTabbedPane.getSelectedIndex() != tab) chartsTabbedPane
				.setSelectedIndex(tab);
		}
	}

	public void setReferencePlot(final int ref) {
		referencePlot = ref;
		removeVerticalMarkers();
		for (final int i : plotNumbers)
			updatePlot(i);
	}

	/**
	 * While the view is reset, raises the top of the range just enough for the
	 * ladder labels, drawn down from the top, to clear the curves under them
	 */
	private void fitRangeToLabels(final ChartPanel p) {
		final JFreeChart c = p.getChart();
		final XYPlot pl = c.getXYPlot();
		final ValueAxis ra = pl.getRangeAxis();
		final Range fitted = fittedRanges.get(c);
		final Rectangle2D area = p.getChartRenderingInfo().getPlotInfo()
			.getDataArea();
		final XYDataset ds = pl.getDataset();
		// Fit once per change (new data, reset view, new plot height), never on
		// every paint: each fit repaints, and a fit per paint never settles
		final boolean atReset = ra.isAutoRange() || fitted == null || ra.getRange()
			.equals(fitted);
		final Double usedHeight = fittedHeights.get(c);
		final boolean changed = labelFitPending.contains(c) || ra.isAutoRange() ||
			usedHeight == null || usedHeight != area.getHeight();
		if (atReset && changed && ds != null && area.getHeight() > 0) {
			labelFitPending.remove(c);
			fittedHeights.put(c, area.getHeight());
			final Range yr = DatasetUtils.findRangeBounds(ds, false);
			if (yr != null) {
				final double h = area.getHeight();
				final double lo = yr.getLowerBound();
				double hi = yr.getUpperBound();
				final ValueAxis da = pl.getDomainAxis();
				final Graphics2D g2 = labelGraphics();
				for (final Object a : pl.getAnnotations()) {
					if (!(a instanceof XYTextAnnotation)) continue;
					final XYTextAnnotation t = (XYTextAnnotation) a;
					// Box of the label anchored at the top of the plot, in px
					g2.setFont(t.getFont());
					final Shape box = TextUtils.calculateRotatedStringBounds(t.getText(),
						g2, (float) da.valueToJava2D(t.getX(), area, pl.getDomainAxisEdge()),
						(float) area.getMinY(), t.getTextAnchor(), t.getRotationAngle(), t
							.getRotationAnchor());
					if (box == null) continue; // Empty label, e.g. an unnamed band
					final Rectangle2D b = box.getBounds2D();
					final double depth = b.getMaxY() - area.getMinY() + labelGap;
					if (depth >= h) continue; // Cannot fit anyway
					// Curves under the label, one data point wider on each side
					final double x0 = da.java2DToValue(b.getMinX(), area, pl
						.getDomainAxisEdge());
					final double x1 = da.java2DToValue(b.getMaxX(), area, pl
						.getDomainAxisEdge());
					final double top = maxInStrip(ds, FastMath.min(x0, x1) - 1, FastMath
						.max(x0, x1) + 1);
					// The label ends depth px below hi: hi - depth * (hi - lo) / h >= top
					if (!Double.isNaN(top)) hi = FastMath.max(hi, (top * h - lo *
						depth) / (h - depth));
				}
				g2.dispose();
				final Range target = new Range(lo, hi);
				fittedRanges.put(c, target);
				if (ra.isAutoRange() || !target.equals(ra.getRange())) ra.setRange(
					target);
			}
		}
		// Keep the labels at the top of the view
		final double upper = ra.getUpperBound();
		for (final Object a : pl.getAnnotations()) {
			if (a instanceof XYTextAnnotation && ((XYTextAnnotation) a)
				.getY() != upper) ((XYTextAnnotation) a).setY(upper);
		}
	}

	private static Graphics2D labelGraphics() {
		final Graphics2D g2 = new BufferedImage(1, 1, BufferedImage.TYPE_INT_RGB)
			.createGraphics();
		g2.setRenderingHint(RenderingHints.KEY_FRACTIONALMETRICS,
			RenderingHints.VALUE_FRACTIONALMETRICS_ON);
		return g2;
	}

	/** Largest y of any series in [x0, x1], NaN if there is none */
	private static double maxInStrip(final XYDataset ds, final double x0,
		final double x1)
	{
		double max = Double.NaN;
		for (int s = 0; s < ds.getSeriesCount(); s++) {
			for (int i = 0; i < ds.getItemCount(s); i++) {
				final double x = ds.getXValue(s, i);
				final double y = ds.getYValue(s, i);
				if (x >= x0 && x <= x1 && !Double.isNaN(y) && !(y <= max)) max = y;
			}
		}
		return max;
	}

	void updateProfile(final Roi roi) {
		final DataSeries profile = laneProfile(imp, roi);
		// Assume plotsData, chartsMainPanel was reset
		plotNumbers.add(profile.getLane());
		plotsData.add(profile);
		final String xLabel = "Distance (px)";
		final String yLabel = "Grayscale Value";
		final XYSeriesCollection dataset = new XYSeriesCollection();
		XYPlot thePlot = new XYPlot();
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
			final NumberFormat format = NumberFormat.getNumberInstance();
			format.setMaximumFractionDigits(1);
			final XYToolTipGenerator generator = new StandardXYToolTipGenerator(
				"({1} {2})", format, format);
			newPlot.getRenderer().setDefaultToolTipGenerator(generator);

			newChart.getTitle().setMargin(new org.jfree.chart.ui.RectangleInsets(15,
				5, 15, 5));
			newChart.getTitle().setPaint(Color.BLACK);
			newPlot.setDomainGridlinePaint(Color.DARK_GRAY);
			newPlot.setRangeGridlinePaint(Color.DARK_GRAY);
			newPlot.setDomainCrosshairVisible(true);
			newPlot.setRangeCrosshairVisible(true);
			thePlot = newPlot;
			final ChartPanel chartPanel = new ChartPanel(newChart);
			chartPanel.addChartMouseListener(this);
			newChart.addProgressListener(e -> {
				if (e.getType() != ChartProgressEvent.DRAWING_FINISHED) return;
				try {
					fitRangeToLabels(chartPanel);
				}
				catch (final RuntimeException ex) {
					// Never let the label fit stop the chart from being painted
					log.error("Could not fit the plot range to the labels", ex);
				}
			});
			chartPanels.add(chartPanel);
		}

		// Auto range (also used by the zoom reset) spans exactly the profile's
		// domain and the smallest to largest value of the curves displayed
		for (final ValueAxis a : new ValueAxis[] { thePlot.getDomainAxis(), thePlot
			.getRangeAxis() })
		{
			a.setLowerMargin(0);
			a.setUpperMargin(0);
			if (a instanceof NumberAxis) ((NumberAxis) a).setAutoRangeIncludesZero(
				false);
			a.setAutoRange(true);
		}
	}

	public void updatePlot(final Roi r) {
		final int n = Integer.parseInt(r.getName().substring((5)));
		updatePlot(n);
	}

	public void updatePlot(final int ln) {
		Color plotBGColor = plotSelColor;
		Color vMarkerColor = vMarkerRegColor;

		if (plotMode == Plotter.editPeaksMode) {
			plotBGColor = plotAddSelColor;
			vMarkerColor = vMarkerEditPeakColor;
		}

		final XYSeriesCollection dataset = new XYSeriesCollection();
		Collections.sort(plotsData);

		for (final ChartPanel p : chartPanels) {
			final int plotNumber = Integer.parseInt(p.getChart().getTitle().getText()
				.substring(5));
			if (plotNumber == ln) {
				final JFreeChart c = p.getChart();
				final XYPlot pl = c.getXYPlot();

				if (plotNumber == referencePlot) {
					pl.setBackgroundPaint(plotRefColor);
					c.setBackgroundPaint(plotRefColor);
				}
				else if (plotNumber == selected) {
					pl.setBackgroundPaint(plotBGColor);
					c.setBackgroundPaint(plotBGColor);
				}
				else {
					pl.setBackgroundPaint(plotUnselColor);
					c.setBackgroundPaint(plotUnselColor);
				}

				// Clear Markers and Annotations
				if (c.getXYPlot().getDomainMarkers(Layer.BACKGROUND) != null) c
					.getXYPlot().clearDomainMarkers();
				if (c.getXYPlot().getAnnotations() != null) c.getXYPlot()
					.clearAnnotations();

				// Plot the data series
				final LegendItems legendItems = new LegendItems();
				for (final DataSeries d : plotsData) {
					if (d.getLane() == ln) {
						if (d.getItemCount() > 0) {
							final int k = d.getType();
							dataset.addSeries(d);
							pl.setSeriesRenderingOrder(SeriesRenderingOrder.FORWARD);
							final int seriesIdx = dataset.getSeriesIndex(d.getKey());
							final XYLineAndShapeRenderer renderer =
								(XYLineAndShapeRenderer) pl.getRenderer();
							renderer.setSeriesShapesVisible(seriesIdx, false);
							renderer.setSeriesLinesVisible(seriesIdx, true);

							if (k == DataSeries.PROFILE) {
								renderer.setSeriesPaint(seriesIdx, profileColor);
								renderer.setSeriesStroke(seriesIdx, dataStroke);
								final LegendItem li = new LegendItem("Profile");
								li.setFillPaint(d.getColor());
								legendItems.add(li);
							}
							if (k == DataSeries.BACKGROUND) {
								renderer.setSeriesPaint(seriesIdx, bgColor);
								renderer.setSeriesStroke(seriesIdx, dataStroke);
								final LegendItem li = new LegendItem("Background");
								li.setFillPaint(d.getColor());
								legendItems.add(li);
							}
							if (k == DataSeries.GAUSS_BG) {
								renderer.setSeriesPaint(seriesIdx, gaussColor);
								renderer.setSeriesStroke(seriesIdx, dataStroke);
								final LegendItem li = new LegendItem("Peaks");
								if (!legendItems.contains(li)) {
									li.setFillPaint(d.getColor());
									legendItems.add(li);
								}
							}
							if (k == DataSeries.FITTED) {
								renderer.setSeriesPaint(seriesIdx, fittedColor);
								renderer.setSeriesStroke(seriesIdx, dataStroke);
								final LegendItem li = new LegendItem("Fit");
								li.setFillPaint(d.getColor());
								legendItems.add(li);
							}
							if (k == DataSeries.CUSTOMPEAKS) {
								renderer.setSeriesPaint(seriesIdx, vMarkerEditPeakColor);
								final Shape dot = new Ellipse2D.Double(0, 0, 6, 6);
								renderer.setSeriesShape(seriesIdx, dot);
								renderer.setSeriesShapesVisible(seriesIdx, true);
								renderer.setSeriesLinesVisible(seriesIdx, false);
								final LegendItem li = new LegendItem("Custom Peaks");
								li.setFillPaint(d.getColor());
								legendItems.add(li);
							}
						}
					}
				}
				final Range domain = pl.getDomainAxis().getRange();
				final Range range = pl.getRangeAxis().getRange();
				c.getXYPlot().setDataset(dataset);
				labelFitPending.add(c); // The curves under the labels may have changed
				if (rangeSet.contains(ln)) { // Keep the current view
					pl.getDomainAxis().setRange(domain);
					pl.getRangeAxis().setRange(range);
				}
				else {
					pl.getDomainAxis().setAutoRange(true);
					pl.getRangeAxis().setAutoRange(true);
					rangeSet.add(ln);
				}
				pl.setFixedLegendItems(legendItems);
				c.getLegend().setPosition(RectangleEdge.RIGHT);

				// Plot vertical markers
				for (final VerticalMarker m : verticalMarkers) {
					if (m.getLane() == ln) {
						final double height = c.getXYPlot().getRangeAxis().getUpperBound();
						final double offset = 0;
						final XYTextAnnotation label = new XYTextAnnotation(m.getName(), m
							.getValue() - offset, height);
						if (m.getType() == VerticalMarker.VMARK) {
							label.setPaint(vMarkerColor);
						}
						else if (m.getType() == VerticalMarker.BMARK) {
							label.setPaint(bMarkerColor);
						}
						label.setFont(labelFont);
						label.setRotationAnchor(TextAnchor.BOTTOM_RIGHT);
						label.setTextAnchor(TextAnchor.TOP_RIGHT);
						label.setRotationAngle(-Math.PI / 2);

						c.getXYPlot().addAnnotation(label);
						c.getXYPlot().addDomainMarker(m, Layer.BACKGROUND);
					}
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

	public void reloadTabs() {
		int t = chartsTabbedPane.getSelectedIndex();
		chartsTabbedPane.removeAll();
		chartTabs = new ArrayList<>();
		if (chartPanels.size() == 0) return;
		final Iterator<ChartPanel> chartIter = chartPanels.iterator();
		int i = 0;
		while (chartIter.hasNext()) {
			JPanel p = new JPanel();
			if (!chartTabs.isEmpty() && i < (rows * cols)) {
				p = chartTabs.get(chartTabs.size() - 1);
			} else {
				p.setBorder(new EmptyBorder(5, 5, 5, 5));
				p.setLayout(new GridLayout(rows, cols));
				chartTabs.add(p);
				final int tab = chartTabs.size() - 1;
				chartsTabbedPane.addTab("Lanes " + (tab * (rows * cols) + 1) + " - " +
					(tab * (rows * cols) + 4), p);
				i = 0;
			}
			p.add(chartIter.next());
			i++;
		}
		while (i < (rows * cols)) {
			// Add placeholders for empty plots
			chartTabs.get(chartTabs.size() - 1).add(new JPanel());
			i++;
		}
		int tabCount = chartsTabbedPane.getTabCount();
		if (t != -1 && tabCount != 0) {
			if (t < tabCount) chartsTabbedPane.setSelectedIndex(t);
			else chartsTabbedPane.setSelectedIndex(tabCount - 1);
		}
		else if (tabCount > 0) chartsTabbedPane.setSelectedIndex(0);

	}

	public void resetData() {
		plotsData.clear();
		plotNumbers.clear();
		rangeSet.clear();
		fittedRanges.clear();
		fittedHeights.clear();
		labelFitPending.clear();
		chartPanels.clear();
		removeVerticalMarkers();
	}

	public void savePlots(final String savePath) {
		final Iterator<File> it = FileUtils.iterateFiles(new File(savePath),
			new String[] { "png", "pdf", "svg" }, false);
		while (it.hasNext()) {
			it.next().delete();
		}

		for (final ChartPanel p : chartPanels) {
			final String plotfile = savePath + p.getChart().getTitle().getText();
			final double x = PageSize.LETTER.getWidth() * 0.8;
			final double y = x / 1.6;
			final Rectangle2D r = new Rectangle2D.Double(0, 0, x, y);

			try { // Save PNG
				ChartUtils.saveChartAsPNG(new File(plotfile + ".png"), p.getChart(),
					(int) x, (int) y);
			}
			catch (final IOException e) {
				log.error("Could not save " + plotfile + ".png", e);
			}

			// Save PDF
			try (FileOutputStream pdf = new FileOutputStream(plotfile + ".pdf")) {
				final Rectangle ps = new Rectangle((float) x, (float) y);
				final com.itextpdf.text.Document doc = new com.itextpdf.text.Document(
					ps, 20, 20, 20, 20);
				final PdfWriter writer = PdfWriter.getInstance(doc, pdf);
				doc.open();
				final PdfContentByte cb = writer.getDirectContent();
				final PdfTemplate t = cb.createTemplate((float) x, (float) y);
				final Graphics2D g = new PdfGraphics2D(t, (float) x, (float) y);

				p.getChart().draw(g, r);
				g.dispose();
				cb.addTemplate(t, 0, 0);
				doc.close();
			}
			catch (DocumentException | IOException e) {
				log.error("Could not save " + plotfile + ".pdf", e);
			}

			try { // Save SVG
				final SVGGraphics2D svgGen = new SVGGraphics2D((int) x, (int) y);
				p.getChart().draw(svgGen, r);
				File out = new File(plotfile + ".svg");
				SVGUtils.writeToSVG(out, svgGen.getSVGElement());
			}
			catch (final IOException e) {
				log.error("Could not save " + plotfile + ".svg", e);
			}
		}
	}

	@Override
	public void chartMouseClicked(final ChartMouseEvent e) {
		if (plotMode != Plotter.regMode) {
			final JFreeChart c = e.getChart();
			for (final ChartPanel p : chartPanels) {
				if (p.getChart().equals(c)) {
					final int ln = Integer.parseInt(p.getChart().getTitle().getText()
						.substring(5));
					final XYPlot plot = (XYPlot) p.getChart().getPlot(); // your plot
					final Point2D pt = p.translateScreenToJava2D(e.getTrigger()
						.getPoint());
					final Rectangle2D plotArea = p.getScreenDataArea();
					final double xi = plot.getDomainAxis().java2DToValue(pt.getX(),
						plotArea, plot.getDomainAxisEdge());
					if (plotArea.contains(pt)) {
						for (final DataSeries d : getProfiles()) {
							if (d.getLane() == ln) {
								final double[] x = d.getX().toArray();
								final double[] y = d.getY().toArray();
								final PolynomialSplineFunction f = new LinearInterpolator()
									.interpolate(x, y);
								// Clicks in the axis margin fall outside the profile
								if (!f.isValidPoint(xi)) continue;
								final double yi = f.value(xi);
								final DataSeries cp = getPlotsCustomPeaks(ln);
								boolean found = false;
								for (int i = 0; i < cp.getItemCount(); i++) {
									if (FastMath.abs((double) cp.getX(i) -
										xi) <= Fitter.peakDistanceTol)
									{
										found = true;
										cp.remove(i); // REMOVE point
									}
								}
								if (!found) {
									cp.addOrUpdate(xi, yi); // ADD Point
								}
								updatePlot(ln);
							}
						}
					}
				}
			}
		}
		else {
			e.getChart().getXYPlot().setDomainCrosshairVisible(false);
			e.getChart().getXYPlot().setRangeCrosshairVisible(false);
		}
	}

	@Override
	public void chartMouseMoved(final ChartMouseEvent e) {
		if (plotMode != Plotter.regMode) {
			final JFreeChart c = e.getChart();
			for (final ChartPanel p : chartPanels) {
				final XYPlot plot = (XYPlot) p.getChart().getPlot();
				if (p.getChart().equals(c)) {
					final Point2D pt = p.translateScreenToJava2D(e.getTrigger()
						.getPoint());
					final Rectangle2D plotArea = p.getScreenDataArea();

					final double x = plot.getDomainAxis().java2DToValue(pt.getX(),
						plotArea, plot.getDomainAxisEdge());
					final double y = plot.getRangeAxis().java2DToValue(pt.getY(),
						plotArea, plot.getRangeAxisEdge());

					Paint crossHairColor = plot.getDomainGridlinePaint();
					if (plotMode == Plotter.editPeaksMode) crossHairColor =
						vMarkerEditPeakColor;

					plot.setDomainCrosshairPaint(crossHairColor);
					plot.setRangeCrosshairPaint(crossHairColor);
					plot.setDomainCrosshairValue(x);
					plot.setRangeCrosshairValue(y);
					plot.setDomainCrosshairVisible(true);
					plot.setRangeCrosshairVisible(true);
				}
				else {
					plot.setDomainCrosshairVisible(false);
					plot.setRangeCrosshairVisible(false);
				}
			}
		}
	}

	/** A chart's legend, with one item per kind of curve. */
	private class LegendItems extends LegendItemCollection {

		private static final long serialVersionUID = 1L;

		public LegendItems() {
			super();
		}

		private boolean contains(final LegendItem li) {
			for (int i = 0; i < this.getItemCount(); i++) {
				if (this.get(i).getLabel().equals(li.getLabel())) return true;
			}
			return false;
		}
	}

}

/** A labelled vertical line on a lane's chart, such as a ladder band. */
class VerticalMarker extends ValueMarker {

	private static final long serialVersionUID = 1L;

	// Possible types
	final static int VMARK = 0; // Vertical Position
	final static int BMARK = 1; // Molecular Weight bands

	private final int lane; // Reference GEL LANE
	private final int type;
	private String name; // Name for reference

	public VerticalMarker(final String name, final int lane, final int type,
		final double x, final Color color, final Stroke stroke)
	{
		super(x, color, stroke);
		this.name = name;
		this.lane = lane;
		this.type = type;
	}

	void setName(final String name) {
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

/**
 * A curve on a lane's chart: the profile, the background, a peak, the fit, or
 * the custom peaks. x is the position along the lane in pixels, y the
 * intensity in gray values.
 */
class DataSeries extends XYSeries implements Comparable<DataSeries> {

	private static final long serialVersionUID = 1L;

	// Series in a chart's dataset need unique keys; they aren't displayed
	private static final AtomicLong keys = new AtomicLong();

	private static String uniqueKey(final int type) {
		return type + "-" + keys.incrementAndGet();
	}

	private final String name; // Name for Legend
	private final int lane; // Reference
	private final int type; // Type of function
	private Color color; // Plot color

	// Possible Types
	final static int PROFILE = 0;
	final static int GAUSS_BG = 2;
	final static int BACKGROUND = 2000;
	final static int FITTED = 2001;
	final static int CUSTOMPEAKS = 2002;

	public DataSeries(final String name, final int lane, final int type,
		final RealVector x, final RealVector y, final Color color)
	{
		super(uniqueKey(type));
		this.name = name;
		this.lane = lane;
		this.type = type;
		this.color = color;
		if (x.getDimension() == y.getDimension()) {
			if (x.getDimension() != 0) {
				for (int r = 0; r < x.getDimension(); r++) {
					add(x.getEntry(r), y.getEntry(r));
				}
			}
		}
		else throw new DimensionMismatchException(x.getDimension(), y
			.getDimension());
	}

	public DataSeries(final String name, final int lane, final int type,
		final RealVector x, final UnivariateFunction[] function,
		final Color color)
	{
		super(uniqueKey(type));
		this.name = name;
		this.lane = lane;
		this.type = type;
		this.color = color;

		if (x.getDimension() != 0) {
			RealVector y = new ArrayRealVector();
			for (int i = 0; i < function.length; i++) {
				y = (i == 0) ? new ArrayRealVector(x.map(function[i])) : y.add(x.map(
					function[i]));
			}

			for (int r = 0; r < x.getDimension(); r++) {
				add(x.getEntry(r), y.getEntry(r));
			}
		}
	}

	public DataSeries(final String name, final int lane, final int type,
		final RealVector x, final UnivariateFunction function, final Color color)
	{
		super(uniqueKey(type));
		this.name = name;
		this.lane = lane;
		this.type = type;
		this.color = color;
		final RealVector y = new ArrayRealVector(x.map(function));
		for (int r = 0; r < x.getDimension(); r++) {
			add(x.getEntry(r), y.getEntry(r));
		}
	}

	public String getName() {
		return name;
	}

	public int getLane() {
		return lane;
	}

	public int getType() {
		return type;
	}

	public RealVector getX() {
		RealVector x = new ArrayRealVector();
		for (int i = 0; i < getItemCount(); i++) {
			x = x.append((double) getX(i));
		}
		return x;
	}

	public RealVector getY() {
		RealVector y = new ArrayRealVector();
		for (int i = 0; i < getItemCount(); i++) {
			y = y.append((double) getY(i));
		}
		return y;
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
