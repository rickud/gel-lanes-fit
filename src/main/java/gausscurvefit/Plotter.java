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
import java.awt.Paint;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.event.WindowEvent;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.WindowConstants;
import javax.swing.border.EmptyBorder;

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
import org.jfree.chart.JFreeChart;
import org.jfree.chart.LegendItem;
import org.jfree.chart.LegendItemCollection;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.labels.StandardXYToolTipGenerator;
import org.jfree.chart.labels.XYToolTipGenerator;
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

import ij.IJ;
import ij.ImagePlus;
import ij.gui.ProfilePlot;
import ij.gui.Roi;

class Plotter extends JFrame implements ChartMouseListener {
	
	final double SW = IJ.getScreenSize().getWidth();
	final double SH = IJ.getScreenSize().getHeight();
	
	static final int noRefPlot = -1; // no reference plot
	private int referencePlot = noRefPlot;

	// Colors are listed here for consistent, easy modification
	static final Color vMarkerRegColor = new Color(127, 127, 127);
	static final Color vMarkerEditPeakColor = new Color(0, 185, 19);
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
	static final int editPeaksMode = 1;
	static final int refMode = 2;

	static final Color plotSelColor = new Color(240, 240, 240);
	static final Color plotUnselColor = Color.WHITE;
	static final Color plotAddSelColor = new Color(192, 255, 185);
	
	private final ImagePlus imp;
	private List<ChartPanel> chartPanels;
	private List<DataSeries> plotsData;
	private List<VerticalMarker> verticalMarkers;
	
	private List<JPanel> chartTabs;
	private final JTabbedPane chartsTabbedPane;

	private final int rows = 2; // Number of plot Rows in display
	private final int cols = 2; // Number of plot Rows in display
	private int selected; // Which plot is highlighted
	private int plotMode = regMode;

	public Plotter(final Context context, final ImagePlus imp, final List<Roi> rois)
	{
		context.inject(this);
		setDefaultCloseOperation(WindowConstants.HIDE_ON_CLOSE);
		setBounds((int) (SW*0.3), 0, (int) (SW*0.7), (int) (SH*0.9));
		this.setTitle("Profiles of " + imp.getShortTitle());
		this.imp = imp;
		chartPanels = new ArrayList<>();
		chartTabs = new ArrayList<>();
		plotsData = new ArrayList<>();
		verticalMarkers = new ArrayList<>();
		
		chartsTabbedPane = new JTabbedPane();
		selected = selectedNone;
		
		for (final Roi r : rois) {
			updateProfile(r);
		}
		reloadTabs();

		this.getContentPane().add(chartsTabbedPane, BorderLayout.CENTER);
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

	public void addDataSeries(DataSeries data) {
		plotsData.add(data);
	}
	
	public void addDataSeries(final List<DataSeries> data) {
		plotsData.addAll(data);
	}

	public void removeVerticalMarkers() {
		verticalMarkers = new ArrayList<>();
		for (ChartPanel c : chartPanels)
			updatePlot(Integer.parseInt(c.getChart().getTitle().getText().substring(5)));
	}

	void setVLine(final int ln, final double x) {
		Color vMarkerColor = vMarkerRegColor;
		if (plotMode == Plotter.editPeaksMode) {
			vMarkerColor = vMarkerEditPeakColor;
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
	
	public ArrayList<DataSeries> getPlotsCustomPeaks() {
		final ArrayList<DataSeries> custom = new ArrayList<>();
		for (final DataSeries d : plotsData) {
			if (d.getType() == DataSeries.CUSTOMPEAKS) custom.add(d);
		}
		return custom;
	}
	
	public List<DataSeries> getPlotsData() {
		return plotsData;
	}
	
	public List<VerticalMarker> getPlotsVerticalMarkers() {
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
		if (selected != selectedNone) {
			int tab = (selected - 1) / (rows * cols);
			if (chartsTabbedPane.getSelectedIndex() != tab)
				chartsTabbedPane.setSelectedIndex(tab);
		}
	}

	public void setMode(final int plotMode) {
		this.plotMode = plotMode;
	}

	public void setReferencePlot(final int ref, List<Peak> peaks) {
		referencePlot  = ref;
		if (ref != noRefPlot) {
			for (DataSeries d : plotsData) {
				int ln = d.getLane();
				if (ln != ref) {
					ArrayList<VerticalMarker> bands = new ArrayList<>();
					for (Peak p : peaks) {
						double x = p.getMean();
						VerticalMarker m = new VerticalMarker(p.getName(), ln, VerticalMarker.BMARK, x, Plotter.bMarkerColor, Plotter.bMarkerStroke);
						bands.add(m);
					}
					verticalMarkers.addAll(bands);
				}
			}
		} else {
			removeVerticalMarkers();
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
			NumberFormat format = NumberFormat.getNumberInstance();
			format.setMaximumFractionDigits(1);
			XYToolTipGenerator generator =
			    new StandardXYToolTipGenerator("({1} {2})", format, format);
			newPlot.getRenderer().setBaseToolTipGenerator(generator);
			newChart.getTitle().setMargin(new RectangleInsets(15, 5, 15, 5));
			newChart.getTitle().setPaint(Color.BLACK);

			newPlot.setDomainGridlinePaint(Color.DARK_GRAY);
			newPlot.setRangeGridlinePaint(Color.DARK_GRAY);
			newPlot.setDomainCrosshairVisible(true);
			newPlot.setRangeCrosshairVisible(true);
			thePlot = newPlot;
			final ChartPanel chartPanel = new ChartPanel(newChart);
			chartPanel.addChartMouseListener(this);
			chartPanels.add(chartPanel);
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
						if (d.getItemCount() > 0) {
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
								renderer.setSeriesPaint(seriesIdx, vMarkerEditPeakColor);
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
			if (d.getType() != DataSeries.PROFILE && d.getType() != DataSeries.CUSTOMPEAKS) 
				dataIter.remove();
		}
	}

	public void reloadTabs() {
		chartsTabbedPane.removeAll();
		chartTabs = new ArrayList<>();
		Iterator<ChartPanel> chartIter = chartPanels.iterator();
		int i = 0;
		while (chartIter.hasNext()) {
			JPanel p = new JPanel();
			if (!chartTabs.isEmpty() && i < (rows * cols)) {
				p = chartTabs.get(chartTabs.size() - 1);
			}	else {
				p.setBorder(new EmptyBorder(5, 5, 5, 5));
				p.setLayout(new GridLayout(rows, cols));
				chartTabs.add(p);
				int tab = chartTabs.size() - 1;
				chartsTabbedPane.addTab("Lanes " + (tab*(rows*cols) + 1) +" - " + (tab*(rows*cols) + 4), p);
				i = 0;
			}
			p.add(chartIter.next());
			i++;
		}
		if (selected == selectedNone)
			chartsTabbedPane.setSelectedIndex(0);
		else 
			chartsTabbedPane.setSelectedIndex(selected/(rows*cols)-1);
	}
	
	public void resetData() {
		plotsData = new ArrayList<>();
		chartPanels = new ArrayList<>();
	}

	
	@Override
	public void chartMouseClicked(ChartMouseEvent e) {
		if (plotMode != Plotter.regMode) {
			JFreeChart c = e.getChart();
			for (ChartPanel p : chartPanels) {
				if (p.getChart().equals(c)) {
					int ln = Integer.parseInt(p.getChart().getTitle().getText().substring(5));
					XYPlot plot = (XYPlot) p.getChart().getPlot(); // your plot
					Point2D pt = p.translateScreenToJava2D(e.getTrigger().getPoint());
					Rectangle2D plotArea = p.getScreenDataArea();
					double xi = plot.getDomainAxis().java2DToValue(pt.getX(), plotArea, plot.getDomainAxisEdge());
					if (plotArea.contains(pt)) {
						for (DataSeries d : getProfiles()) {
							if (d.getLane() == ln) {
								double[] x = d.getX().toArray();
								double[] y = d.getY().toArray();
								final PolynomialSplineFunction f = new LinearInterpolator().interpolate(x, y);
								double yi = f.value(xi);
								
								for (DataSeries cp : getPlotsCustomPeaks()) { 
									if (cp.getLane() == ln) {
										boolean found = false;
										for (int i = 0; i < cp.getItemCount(); i++) {
											if (FastMath.abs((double) cp.getX(i) - xi) <= Fitter.peakDistanceTol) {
												found = true;
												cp.remove(i); // REMOVE point
											}
										} 
										if (!found) {
											cp.addOrUpdate(xi, yi); // ADD Point
										}
									}
								}
								updatePlot(ln);
							}
						}
					}
				}
			}
		} else {
			e.getChart().getXYPlot().setDomainCrosshairVisible(false);
			e.getChart().getXYPlot().setRangeCrosshairVisible(false);
		}
	}

	@Override
	public void chartMouseMoved(ChartMouseEvent e) {
		if (plotMode != Plotter.regMode) {
			JFreeChart c = e.getChart();
			for (ChartPanel p : chartPanels) {
				XYPlot plot = (XYPlot) p.getChart().getPlot(); 
				if (p.getChart().equals(c)) {
					Point2D pt = p.translateScreenToJava2D(e.getTrigger().getPoint());
					Rectangle2D plotArea = p.getScreenDataArea();
					
					double x = plot.getDomainAxis().java2DToValue(pt.getX(), plotArea, plot.getDomainAxisEdge());
					double y = plot.getRangeAxis().java2DToValue(pt.getY(), plotArea, plot.getRangeAxisEdge());

					Paint crossHairColor = plot.getDomainGridlinePaint();
					if (plotMode == Plotter.editPeaksMode)
						crossHairColor = vMarkerEditPeakColor;
				
					plot.setDomainCrosshairPaint(crossHairColor);
					plot.setRangeCrosshairPaint(crossHairColor);
					plot.setDomainCrosshairValue(x);
					plot.setRangeCrosshairValue(y);
					plot.setDomainCrosshairVisible(true);					
					plot.setRangeCrosshairVisible(true);
				} else {
					plot.setDomainCrosshairVisible(false);
					plot.setRangeCrosshairVisible(false);
				}
			}
		}
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
		if (x.getDimension() == y.getDimension()) { 
			if (x.getDimension() != 0) {
				for (int r = 0; r < x.getDimension(); r++) {
					add(x.getEntry(r), y.getEntry(r));
				}
			}
		} else
			throw new DimensionMismatchException(x.getDimension(), y.getDimension()); 
	}

	public DataSeries(final String name, final int lane, final int type,
		final RealVector x, final UnivariateFunction[] function,
		final Color color)
	{
		super(type);
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
		super(type);
		this.name = name;
		this.lane = lane;
		this.type = type;
		this.color = color;
		RealVector y = new ArrayRealVector(x.map(function));
		for (int r = 0; r < x.getDimension(); r++) {
			add(x.getEntry(r), y.getEntry(r));
		}
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

	public void addAll(final RealVector x, final RealVector y) {
		if (x.getDimension() == y.getDimension()) {
			for (int i = 0; i < x.getDimension(); i++) {
				add(x.getEntry(i), y.getEntry(i));
			}
		} else
			throw new DimensionMismatchException(x.getDimension(), y.getDimension()); 
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
