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
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Paint;
import java.awt.Shape;
import java.awt.Stroke;
import java.awt.event.WindowEvent;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;
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

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
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
import org.jfree.chart.labels.StandardXYToolTipGenerator;
import org.jfree.chart.labels.XYToolTipGenerator;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.SeriesRenderingOrder;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.ui.Layer;
import org.jfree.chart.ui.RectangleEdge;
import org.jfree.chart.ui.TextAnchor;
import org.jfree.data.Range;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.scijava.Context;
import org.w3c.dom.DOMImplementation;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.ProfilePlot;
import ij.gui.Roi;

class Plotter extends JFrame implements ChartMouseListener {

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
	private List<ChartPanel> chartPanels;
	private List<DataSeries> plotsData;
	private List<Integer> plotNumbers;
	private List<VerticalMarker> verticalMarkers;

	private List<JPanel> chartTabs;
	private final JTabbedPane chartsTabbedPane;

	private final int rows = 2; // Number of plot Rows in display
	private final int cols = 2; // Number of plot Rows in display
	private int selected = MainDialog.noLaneSelected; // Which plot is highlighted
	private int plotMode = regMode;

	public Plotter(final Context context, final ImagePlus imp)
	{
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

	public void addDataSeries(final DataSeries data) {
		plotsData.add(data);
	}

	public void addDataSeries(final List<DataSeries> data) {
		plotsData.addAll(data);
	}

	public void addVerticalMarkers(List<Peak> peaks) {
		for (final DataSeries d : plotsData) {
			final int ln = d.getLane();
			final ArrayList<VerticalMarker> bands = new ArrayList<>();
			for (final Peak p : peaks) {
				final double x = p.getMean();
				final VerticalMarker m = new VerticalMarker(p.getName(), ln,
						VerticalMarker.BMARK, x, Plotter.bMarkerColor,
						Plotter.bMarkerStroke);
				bands.add(m);
			}
			verticalMarkers.addAll(bands);
		}
	}
	
	public void removeVerticalMarkers() {
		verticalMarkers = new ArrayList<>();
		for (int i : plotNumbers)
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

	public DataSeries getPlotsCustomPeaks(int lane) {
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

	void setPlotMode(int plotMode) {
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
		for (int i : plotNumbers)
			updatePlot(i);
	}

	void updateProfile(final Roi roi) {
		final DataSeries profile = getLaneProfile(roi);
		// Assume plotsData, chartsMainPanel was reset
		plotNumbers.add(profile.getLane());
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
			final NumberFormat format = NumberFormat.getNumberInstance();
			format.setMaximumFractionDigits(1);
			final XYToolTipGenerator generator = new StandardXYToolTipGenerator(
				"({1} {2})", format, format);
			newPlot.getRenderer().setDefaultToolTipGenerator(generator);

			newChart.getTitle().setMargin(new org.jfree.chart.ui.RectangleInsets(15, 5, 15, 5));
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

		thePlot.getDomainAxis().setLowerMargin(0);
		thePlot.getDomainAxis().setUpperMargin(0);
		thePlot.getRangeAxis().setLowerMargin(0);
		thePlot.getRangeAxis().setUpperMargin(0);
		final double min = 0.95 * profile.getMinY();
		final double max = 1.2 * profile.getMaxY();
		thePlot.getRangeAxis().setLowerBound(min);
		thePlot.getRangeAxis().setUpperBound(max);
	}

	public void updatePlot(Roi r) {
		int n = Integer.parseInt(r.getName().substring((5)));
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
				int gcount = 0;
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
//				double min = Integer.MAX_VALUE;
//				double max = Integer.MIN_VALUE;
				for (final DataSeries d : plotsData) {
					if (d.getLane() == ln) {
						if (d.getItemCount() > 0) {
							final int k = d.getType();
							String key = "";
							if (k == DataSeries.GAUSS_BG) {
								key = String.format("%d", k + gcount);
								gcount++;
							}
							else {
								key = String.format("%d", k);
							}
							d.setKey(key);
							dataset.addSeries(d);
							pl.setSeriesRenderingOrder(SeriesRenderingOrder.FORWARD);
							final int seriesIdx = dataset.getSeriesIndex(key);
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
//							if (d.getMaxY() > max) max = d.getMaxY();
//							if (d.getMaxY() < min) min = d.getMinY();
						}
					}
				}
				c.getXYPlot().setDataset(dataset);
				Range rng = dataset.getRangeBounds(true);
				double below = rng.getLength() * 0.05;
				double above = rng.getLength() * 0.3;
				c.getXYPlot().getRangeAxis().setLowerBound(rng.getLowerBound()-below);
				c.getXYPlot().getRangeAxis().setUpperBound(rng.getUpperBound()+above);
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
						label.setFont(new Font("Sans Serif", Font.PLAIN, 14));
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
		chartsTabbedPane.removeAll();
		chartTabs = new ArrayList<>();
		if (chartPanels.size() == 0) return;
		final Iterator<ChartPanel> chartIter = chartPanels.iterator();
		GridBagConstraints c = new GridBagConstraints();
		int i = 0;
		while (chartIter.hasNext()) {
			JPanel p = new JPanel();
			if (!chartTabs.isEmpty() && i < (rows * cols)) {
				p = chartTabs.get(chartTabs.size() - 1);
			}
			else {
				p.setBorder(new EmptyBorder(5, 5, 5, 5));
//				p.setLayout(new GridLayout(rows, cols));
				p.setLayout(new GridBagLayout());
				c.weightx = 0.5; c.weighty = 0.5;
				c.fill = GridBagConstraints.BOTH;
				c.ipady = 2; c.ipadx = 2;
				chartTabs.add(p);
				final int tab = chartTabs.size() - 1;
				chartsTabbedPane.addTab("Lanes " + (tab * (rows * cols) + 1) + " - " +
					(tab * (rows * cols) + 4), p);
				i = 0;
			}

			c.gridx = i % cols; c.gridy = i/cols;
			p.add(chartIter.next(), c);
			i++;
		}
		while (i < (rows*cols)) {
			// Add placeholders for empty plots
			c.gridx = i % cols; c.gridy = i/cols;
			chartTabs.get(chartTabs.size() - 1).add(new JPanel(), c);
			i++;
		}
		if (selected != MainDialog.noLaneSelected)
			chartsTabbedPane.setSelectedIndex(selected / (rows * cols) - 1);
	}

	public void resetData() {
		plotsData.clear();
		plotNumbers.clear();
		chartPanels.clear();
		removeVerticalMarkers();
		selected = MainDialog.noLaneSelected;
	}

	public void savePlots(String savePath) {
		Iterator<File> it = FileUtils.iterateFiles(new File(savePath), 
				new String[] {"png", "pdf"}, false);
		while (it.hasNext()) {
			it.next().delete();
		}
		
		for (ChartPanel p : chartPanels) {
			String plotfile = savePath + p.getChart().getTitle().getText();
			float x = PageSize.A4.getWidth();
			float y = PageSize.A4.getHeight() / 2;
			Rectangle2D r = new Rectangle2D.Double(0, 0, x, y);
			
			try { // Save PNG
				ChartUtils.saveChartAsPNG(new File(plotfile + ".png"),
																	p.getChart(), 800, 600);
			}
			catch (IOException e) {
				e.printStackTrace();
			}
			
			try { // Save PDF
				Rectangle ps = new Rectangle(x, y);
				com.itextpdf.text.Document doc 
					= new com.itextpdf.text.Document(ps, 20, 20, 20, 20);
				PdfWriter writer = PdfWriter.getInstance(doc, new FileOutputStream(plotfile + ".pdf"));
				doc.open();
				PdfContentByte cb = writer.getDirectContent();
				PdfTemplate t = cb.createTemplate(x, y);
				Graphics2D  g = new PdfGraphics2D(t, x, y);
				
				p.getChart().draw(g, r);
				g.dispose();
				cb.addTemplate(t, 0, 0);
				doc.close();
			} catch (DocumentException | IOException e) {
				e.printStackTrace();
			}
			
			try (Writer out = new OutputStreamWriter(
				new FileOutputStream(plotfile + ".svg")))
			{ // Save SVG
				DOMImplementation domImpl =
					GenericDOMImplementation.getDOMImplementation();
				org.w3c.dom.Document document = domImpl.createDocument(null, "svg", null);
				SVGGraphics2D svgGen = new SVGGraphics2D(document);
				p.getChart().draw(svgGen, r);
				svgGen.stream(out, true /* use css */);
				out.flush();
				out.close();
			} catch (IOException e) {
				e.printStackTrace();
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

	private class LegendItems extends LegendItemCollection {

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

class VerticalMarker extends ValueMarker {

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

class DataSeries extends XYSeries implements Comparable<DataSeries> {

	private final String name; // Name for Legend
	private final int lane; // Reference
	private final int type; // Type of function
	private Color color; // Plot color

	// Possible Types
	final static int PROFILE = 0;
	final static int GAUSS_BG = 2;
	final static int BACKGROUND = 400;
	final static int FITTED = 401;
	final static int CUSTOMPEAKS = 402;

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
		}
		else throw new DimensionMismatchException(x.getDimension(), y
			.getDimension());
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
		final RealVector y = new ArrayRealVector(x.map(function));
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

//	public void addAll(final RealVector x, final RealVector y) {
//		if (x.getDimension() == y.getDimension()) {
//			for (int i = 0; i < x.getDimension(); i++) {
//				add(x.getEntry(i), y.getEntry(i));
//			}
//		}
//		else throw new DimensionMismatchException(x.getDimension(), y
//			.getDimension());
//	}

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
