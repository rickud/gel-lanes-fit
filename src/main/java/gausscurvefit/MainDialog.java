
package gausscurvefit;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.event.WindowEvent;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSlider;
import javax.swing.JTextField;
import javax.swing.JToggleButton;
import javax.swing.SwingConstants;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

import org.apache.commons.lang3.ArrayUtils;
import org.scijava.Context;
import org.scijava.app.StatusService;
import org.scijava.display.DisplayService;
import org.scijava.plugin.Parameter;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.gui.TextRoi;

@SuppressWarnings("serial")
class MainDialog extends JFrame implements ActionListener, ChangeListener,
	DocumentListener, ItemListener, MouseMotionListener, MouseListener,
	MouseWheelListener
{

	@Parameter
	private DisplayService displayServ;

	@Parameter
	private StatusService statusServ;

	private boolean auto = true; // AUTO Lane size mode is ON
	private boolean selectionUpdate = false; // active updating is off
	private boolean fitDone = false; // Keep track of whether fit data exists
	private boolean addPoint = false;
	private boolean removePoint = false;

	private String roiSelected = "none";
	private String roiPreviouslySelected = "none";
	
	final String warningFit =
			"The current plots will be reset and the current fitting data will be lost.";

	private Plotter plotter;
	private Fitter fitter;

	private final ImagePlus imp;
	private ArrayList<Roi> rois;

	private int IW, IH, LW, LH, LSp, LHOff, LVOff;

	private int nLanes = 6;
	private int degBG = 3; // Order of Background Polynomial
	private double tolPK = 0.05; // Peak detection tolerance as % of range

	private JPanel roiButtonsPanel;
	private JRadioButton buttonAuto;
	private JRadioButton buttonManual;
	private ButtonGroup roiButtons;

	private JPanel textPanel;
	private JLabel labelNLanes;
	private JTextField textNLanes;

	private JPanel sliderPanel;
	private JSlider sliderW;
	private JSlider sliderH;
	private JSlider sliderSp;
	private JSlider sliderHOff;
	private JSlider sliderVOff;

	private JPanel buttonPanel;
	private JButton buttonMeasure;
	private JButton buttonClose;

	private JPanel settingsPanel;
	private JPanel degPanel;
	private JPanel tolPanel;
	private JLabel labelDegBG;
	private JLabel labelTolPK;
	private JTextField textDegBG;
	private JTextField textTolPK;
	private JCheckBox chkBoxBands;
	private JToggleButton buttonAddPoint;
	private JToggleButton buttonRemovePoint;

	private JPanel dialogPanel;
	private final JFrame frame;

	public MainDialog(final Context context, final String string,
		final ImagePlus imp)
	{
		context.inject(this);
		frame = new JFrame(string);
		this.imp = imp;
		setupMainDialog();
	}

	private void setupMainDialog() {
		imp.getCanvas().addMouseMotionListener(this);
		imp.getCanvas().addMouseListener(this);
		imp.getCanvas().addMouseWheelListener(this);
		imp.getWindow().addMouseListener(this);
		imp.getWindow().addMouseMotionListener(this);

		roiButtonsPanel = new JPanel();
		buttonAuto = new JRadioButton("Automatic Rectangle Selection");
		buttonAuto.setActionCommand("Auto");
		buttonAuto.setSelected(true);
		buttonAuto.addActionListener(this);
		buttonManual = new JRadioButton("Manual Rectangle Selection");
		buttonManual.setActionCommand("Manual");
		buttonManual.addActionListener(this);
		roiButtons = new ButtonGroup();
		roiButtonsPanel.setLayout(new BoxLayout(roiButtonsPanel, BoxLayout.Y_AXIS));
		roiButtons.add(buttonAuto);
		roiButtons.add(buttonManual);

		roiButtonsPanel.add(buttonAuto);
		roiButtonsPanel.add(buttonManual);
		roiButtonsPanel.setBorder(new TitledBorder(BorderFactory
			.createEtchedBorder(), "Lane Selection", TitledBorder.LEADING,
			TitledBorder.BELOW_TOP, new Font("Sans", Font.PLAIN, 11)));

		sliderPanel = new JPanel();
		textPanel = new JPanel();
		textPanel.setLayout(new FlowLayout(FlowLayout.LEADING));
		textNLanes = new JTextField(3);
		textNLanes.setText("" + nLanes);
		textNLanes.setVisible(true);
		textNLanes.getDocument().addDocumentListener(this);
		textPanel.add(textNLanes);
		labelNLanes = new JLabel("Number of Lanes");
		textPanel.add(labelNLanes);
		sliderPanel.add(textPanel);

		// Default lane size/offset (Just center 4 lanes in the image)
		IW = imp.getWidth();
		IH = imp.getHeight();
		LW = (int) Math.round(0.8 * IW / nLanes);
		LH = (int) Math.round(IH * 0.8);
		LSp = Math.round((IW - LW * nLanes) / (nLanes + 1));
		LHOff = LSp / 2;
		LVOff = (IH - LH) / 2;

		sliderPanel.setLayout(new BoxLayout(sliderPanel, BoxLayout.Y_AXIS));
		sliderW = makeTitledSlider("Width ( " + LW + " px )", Color.black, 1, IW /
			4, LW);
		sliderPanel.add(sliderW);
		sliderH = makeTitledSlider("Height ( " + LH + " px )", Color.black, IH / 10,
			IH, LH);
		sliderPanel.add(sliderH);
		sliderSp = makeTitledSlider("Space ( " + LSp + " px )", Color.black, 1, IW /
			nLanes, LSp);
		sliderPanel.add(sliderSp);
		sliderHOff = makeTitledSlider("Horizontal Offset ( " + LHOff + " px )",
			Color.black, 0, (int) Math.round(IW * 0.9), LHOff);
		sliderPanel.add(sliderHOff);
		sliderVOff = makeTitledSlider("Vertical Offset ( " + LVOff + " px )",
			Color.black, 0, (int) Math.round(IH * 0.9), LVOff);
		sliderPanel.add(sliderVOff);

		buttonPanel = new JPanel();
		buttonPanel.setLayout(new FlowLayout(FlowLayout.TRAILING));
		buttonMeasure = new JButton("Measure");
		buttonMeasure.addActionListener(this);
		buttonPanel.add(buttonMeasure);
		buttonClose = new JButton("Close");
		buttonClose.addActionListener(this);
		buttonPanel.add(buttonClose);

		settingsPanel = new JPanel();
		degPanel = new JPanel();
		tolPanel = new JPanel();
		labelDegBG = new JLabel("Deg BG");
		labelTolPK = new JLabel("Tol PK");
		textDegBG = new JTextField(3);
		textDegBG.setText("" + degBG);
		textTolPK = new JTextField(3);
		textTolPK.setText("" + tolPK);
		chkBoxBands = new JCheckBox("Show Bands");
		buttonAddPoint = new JToggleButton("Add Point");
		buttonRemovePoint = new JToggleButton("Remove Point");

		degPanel.setLayout(new FlowLayout(FlowLayout.TRAILING));
		degPanel.add(labelDegBG);
		degPanel.add(textDegBG);
		tolPanel.setLayout(new FlowLayout(FlowLayout.TRAILING));
		tolPanel.add(labelTolPK);
		tolPanel.add(textTolPK);
		chkBoxBands.addItemListener(this);
		chkBoxBands.setSelected(false);
		chkBoxBands.setEnabled(false);
		buttonAddPoint.addActionListener(this);
		buttonAddPoint.setEnabled(false);
		buttonRemovePoint.addActionListener(this);
		buttonRemovePoint.setEnabled(false);

		settingsPanel.setLayout(new GridLayout(10, 1));
		settingsPanel.add(degPanel);
		settingsPanel.add(tolPanel);
		settingsPanel.add(chkBoxBands);
		settingsPanel.add(buttonAddPoint);
		settingsPanel.add(buttonRemovePoint);

		dialogPanel = new JPanel();
		dialogPanel.setBackground(Color.darkGray);
		dialogPanel.setLayout(new BorderLayout());
		dialogPanel.add(roiButtonsPanel, BorderLayout.NORTH);
		dialogPanel.add(sliderPanel, BorderLayout.CENTER);
		dialogPanel.add(settingsPanel, BorderLayout.EAST);
		dialogPanel.add(buttonPanel, BorderLayout.SOUTH);

		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		frame.getContentPane().add(dialogPanel);
		frame.setResizable(true);
		frame.validate();
		frame.pack();
		frame.setVisible(true);

		resetAutoROIs(LW, LH, LSp, LHOff, LVOff, nLanes);
		while (rois.size() == 0) {
			System.out.println("Loading...");
			IJ.wait(500); // delay to make sure ROIs have updated
		}
		reDrawROIs(imp, "none");
		imp.killRoi();
	}

	private boolean askUser(final String question) {
		final GenericDialog gd = new GenericDialog("WARNING!");
		gd.addMessage(question);
		gd.showDialog();
		if (gd.wasOKed()) return true;
		return false;
	}

	private double askPoint(final boolean add, final int lane, final double y,
		final double v)
	{
		String title;
		String action;
		if (add) {
			title = "ADD POINT";
			action = "add";
		}
		else {
			title = "REMOVE POINT";
			action = "remove";
		}
		final GenericDialog gd = new GenericDialog(title);
		gd.addMessage("You are about to " + action + " this peak");
		gd.addMessage(String.format("Lane: \t%2d", lane));
		gd.addMessage(String.format("Distance: \t%10.1f", y));
		gd.addMessage(String.format("Intensity: \t%.0f", v));
		if (add) {
			gd.addMessage("Estimate a Standard deviation, \u03C3");
			gd.addNumericField("\u03C3", 1.0, 1);
		}
		gd.showDialog();
		if (gd.wasOKed()) return gd.getNextNumber();
		return 0.0;
	}

	private JSlider makeTitledSlider(final String string, final Color color,
		final int minVal, final int maxVal, final int val)
	{
		final JSlider slider = new JSlider(SwingConstants.HORIZONTAL, minVal,
			maxVal, val);
		final TitledBorder tb = new TitledBorder(BorderFactory.createEtchedBorder(),
			// empty,
			"", TitledBorder.CENTER, TitledBorder.BELOW_TOP, new Font("Sans",
				Font.PLAIN, 11));
		tb.setTitle(string);
		tb.setTitleJustification(TitledBorder.LEFT);
		tb.setTitleColor(color);
		slider.setBorder(tb);
		slider.setMajorTickSpacing((maxVal - minVal) / 10);
		slider.setPaintTicks(true);
		slider.addChangeListener(this);
		return slider;
	}

	private void setSliderTitle(final JSlider slider, final String str) {
		final TitledBorder tb = (TitledBorder) slider.getBorder();
		tb.setTitle(str);
		slider.setBorder(tb);
	}

	private void setSliderPanelEnabled(final boolean enabled) {
		if (enabled) {
			textNLanes.setEnabled(true);
			sliderW.setEnabled(true);
			sliderH.setEnabled(true);
			sliderSp.setEnabled(true);
			sliderHOff.setEnabled(true);
			sliderVOff.setEnabled(true);
		}
		else {
			textNLanes.setEnabled(false);
			sliderW.setEnabled(false);
			sliderH.setEnabled(false);
			sliderSp.setEnabled(false);
			sliderHOff.setEnabled(false);
			sliderVOff.setEnabled(false);
		}
	}

	private void reDrawROIs(final ImagePlus imgPlus, final String roiName) {
		if (!auto) {
			// HouseKeeping: Sort the rois based on name and remove null elements
			final Comparator<Roi> nameComparator = new Comparator<Roi>() {

				@Override
				public int compare(final Roi r1, final Roi r2) {
					if (r1 == null && r2 == null) {
						return 0;
					}
					if (r1 == null) {
						return -1;
					}
					if (r2 == null) {
						return 1;
					}
					final int n1 = Integer.parseInt(r1.getName().substring(5));
					final int n2 = Integer.parseInt(r2.getName().substring(5));
					final int comp = Integer.compare(n1, n2);
					if (comp != 0) {
						return comp;
					}
					return 0;
				}
			};
			Collections.sort(rois, nameComparator);
			final Iterator<Roi> roisIter = rois.iterator();
			int nullIndex = 0;
			while (roisIter.hasNext()) {
				if (roisIter.next() == null) break;
				nullIndex++;
			}
			rois = new ArrayList<>(rois.subList(0, nullIndex));
		}
		final Overlay overlay = new Overlay();
		final Iterator<Roi> roiIter = rois.iterator();
		while (roiIter.hasNext()) {
			final Roi roi = roiIter.next();
			final String label = roi.getName().substring(5);
			final Font font = new Font("SansSerif", Font.BOLD, 20);
			final TextRoi labelRoi = new TextRoi(0, 0, label, font);
			labelRoi.setAntialiased(true);
			final double lh = font.getSize() * 1.4;
			final double lw = font.getSize() * 1.4;
			final double rh = roi.getBounds().getHeight();
			final double rw = roi.getBounds().getWidth();
			final double x0 = roi.getBounds().getX();
			final double y0 = roi.getBounds().getY();
			if (lh > rh || lw > rw) labelRoi.setLocation(x0 + 0.1 * lw, y0 + rh +
				1.1 * lh);
			else labelRoi.setLocation(x0 + 0.2 * lw, y0 + rh - 1.1 * lh);
			if (roi.getName().equals(roiName)) {
				roi.setStrokeColor(Color.YELLOW);
				roi.setStrokeWidth(3);
				if (buttonManual.isSelected() && !roiName.equals("none")) {
					imgPlus.setRoi(roi);
				}
			}
			else {
				roi.setStrokeWidth(1);
				roi.setStrokeColor(Color.RED);
				labelRoi.setFillColor(Color.RED);
			}
			overlay.add(roi);
			overlay.add(labelRoi);
			if (buttonAuto.isSelected()) imgPlus.killRoi();
		}
		imgPlus.setOverlay(overlay);
	}

	private void resetAutoROIs(final int lw, final int lh, final int lsp,
		final int lhoff, final int lvoff, final int nlanes)
	{
		rois = new ArrayList<>();
		for (int i = 0; i < nlanes; i++) {
			final Roi roi = new Roi(lhoff + lw * i + lsp * i, lvoff, lw, lh);
			roi.setName("Lane " + (i + 1));
			rois.add(roi);
		}
	}

	private void changeLaneNumber() {
		if (fitDone && !askUser(warningFit)) {
			textNLanes.setText(Integer.toString(nLanes));
			return;
		}
		nLanes = getNLanes();
		if (nLanes == -1) return;
		resetAutoROIs(LW, LH, LSp, LHOff, LVOff, nLanes);
		reDrawROIs(imp, "none");
		redoProfilePlots();
	}

	private void redoProfilePlots() {
		fitDone = false;
		plotter.resetPlots();
		for (final Roi r : rois)
			plotter.updateProfile(r);
		for (final int i : getAllLaneNumbers())
			fitter.resetCustomPoints(i);
		plotter.plotsMontage();
	}

	public void cleanupAndClose() {
		// Remove Listeners on imp
		if (imp != null) {
			imp.getCanvas().removeMouseListener(this);
			imp.getCanvas().removeMouseMotionListener(this);
			imp.getCanvas().removeMouseWheelListener(this);
			imp.getWindow().removeMouseListener(this);
			imp.getWindow().removeMouseMotionListener(this);
			imp.setOverlay(null);
		}
		plotter.closePlot();
	}

	public int getNLanes() {
		try {
			return Integer.parseInt(textNLanes.getText().trim());
		}
		catch (final NumberFormatException e1) {
			return -1;
		}
	}


	public int getDegBG() {
		try {
			degBG = Integer.parseInt(textDegBG.getText().trim());
			return degBG;
		}
		catch (final NumberFormatException e1) {
			return -1;
		}
	}

	/**
	 * @return the currently selected @param tolpk, the tolerance on the peak
	 *         selection
	 */
	public double getTolPK() {
		try {
			tolPK = Double.parseDouble(textTolPK.getText().trim());
			return tolPK;
		}
		catch (final NumberFormatException e1) {
			return 0.0;
		}
	}

	public Roi getRoi(final String title) {
		final Iterator<Roi> roiIter = rois.iterator();
		while (roiIter.hasNext()) {
			final Roi roi = roiIter.next();
			if (roi.getName().equals(title)) return roi;
		}
		return null;
	}

	/**
	 * @return @param(rois), the current ROI set
	 */
	public ArrayList<Roi> getRois() {
		return rois;
	}

	public String getRoiSelected() {
		return roiSelected;
	}

	public int[] getAllLaneNumbers() {
		final Iterator<Roi> roisIter = rois.iterator();
		final int[] laneNumbers = new int[rois.size()];
		int i = 0;
		while (roisIter.hasNext()) {
			final Roi roi = roisIter.next();
			if (roi != null) {
				final String name = roi.getName();
				laneNumbers[i] = Integer.parseInt(name.substring(5));
				i++;
			}
		}
		return ArrayUtils.subarray(laneNumbers, 0, i);
	}

	public String getRoiPreviouslySelected() {
		return roiPreviouslySelected;
	}

	public void setPlotter(final Plotter plotter) {
		this.plotter = plotter;
	}
	

	@Override
	public synchronized void stateChanged(final ChangeEvent e) {
		if (fitDone && !askUser(warningFit))	return;
		final JSlider slider = (JSlider) e.getSource();
		
		if (slider == sliderW) {
			LW = sliderW.getValue();
			final String str = "Width ( " + LW + " px )";
			setSliderTitle(sliderW, str);
		}
		else if (slider == sliderH) {
			LH = sliderH.getValue();
			final String str = "Height ( " + LH + " px )";
			setSliderTitle(sliderH, str);
		}
		else if (slider == sliderSp) {
			LSp = sliderSp.getValue();
			final String str = "Spacing ( " + LSp + " px )";
			setSliderTitle(sliderSp, str);
		}
		else if (slider == sliderHOff) {
			LHOff = sliderHOff.getValue();
			final String str = "Horizontal Offset ( " + LHOff + " px )";
			setSliderTitle(sliderHOff, str);
		}
		else if (slider == sliderVOff) {
			LVOff = sliderVOff.getValue();
			final String str = "Vertical Offset ( " + LVOff + " px )";
			setSliderTitle(sliderVOff, str);
		}
		resetAutoROIs(LW, LH, LSp, LHOff, LVOff, nLanes);
		reDrawROIs(imp, "none");
		redoProfilePlots();
	}

	// TextField Listeners
	@Override
	public void changedUpdate(final DocumentEvent e) {
		changeLaneNumber();
	}

	@Override
	public void removeUpdate(final DocumentEvent e) {
		changeLaneNumber();
	}

	@Override
	public void insertUpdate(final DocumentEvent e) {
		changeLaneNumber();
	}

	@Override
	public void actionPerformed(final ActionEvent e) {
		if (e.getSource().equals(buttonMeasure)) {
			if (fitDone && !askUser(warningFit)) return;
			if (displayServ.getDisplay("Results Display") != null) displayServ
				.getDisplay("Results Display").close();
			addPoint = false;
			removePoint = false;
			
			buttonAddPoint.setSelected(false);
			buttonRemovePoint.setSelected(false);
			
			for (final MyPlot p : plotter.getPlots())
				p.setSelectedBGColor(Plotter.plotSelColor);
			
			fitter.setInputData(plotter.getProfiles());
			final ArrayList<ArrayList<DataSeries>> fitted = fitter.doFit();
			fitDone = true;
			for (final ArrayList<DataSeries> f : fitted) {
				final MyPlot plot = plotter.getPlots().get(fitted.indexOf(f));
				plot.addDataSeries(f);
				plot.updatePlot();
			}
			plotter.plotsMontage();
			chkBoxBands.setEnabled(true);
			buttonAddPoint.setEnabled(true);
			buttonRemovePoint.setEnabled(true);
		}

		if (e.getSource().equals(buttonAuto)) {
			if (!auto) { // NOT already AUTO
				if (askUser(
					"When switching to AUTO mode, the current lane selections will be reset!"))
				{
					auto = true;
					setSliderPanelEnabled(true);
				}
				else buttonManual.setSelected(true);
			}
		}

		if (e.getSource().equals(buttonManual)) {
			auto = false;
			setSliderPanelEnabled(false);
		}

		if (e.getSource().equals(buttonAddPoint)) {
			addPoint = !addPoint;
			removePoint = false;
			buttonRemovePoint.setSelected(false);
			Color plotBGColor =  addPoint ? Plotter.plotAddSelColor : Plotter.plotSelColor;
			for (final MyPlot p : plotter.getPlots())
				p.setSelectedBGColor(plotBGColor);
		}

		if (e.getSource().equals(buttonRemovePoint)) {
			addPoint = false;
			removePoint = !removePoint;
			buttonAddPoint.setSelected(false);
			Color plotBGColor =  removePoint ? Plotter.plotRemoveSelColor : Plotter.plotSelColor;
			for (final MyPlot p : plotter.getPlots())
				p.setSelectedBGColor(plotBGColor);

		}

		if (e.getSource().equals(buttonClose)) {
			if (askUser("Would you like to exit?")) {
				frame.dispatchEvent(new WindowEvent(frame, WindowEvent.WINDOW_CLOSING));
			}
		}
	}

	@Override
	public void itemStateChanged(final ItemEvent e) {
		if (e.getItemSelectable() == chkBoxBands) {
			// TODO Add feature to show/hide horizontal ticks
			// for detected bands in original image
		}
	}

	@Override
	public void mouseDragged(final MouseEvent e) {
		// Nothing to do
	}

	@Override
	public void mouseMoved(final MouseEvent e) {
		if (e.getSource() == imp.getCanvas()) {
			final int x = ((ImageCanvas) e.getSource()).offScreenX(e.getX());
			final int y = ((ImageCanvas) e.getSource()).offScreenX(e.getY());
			statusServ.showStatus("[" + x + ":" + y + "]");
			if (!selectionUpdate) { // Not dragging an ROI
				String roiCurrent = "none"; // None selected
				int roiN = 0;
				final Iterator<Roi> roisIter = rois.iterator();
				while (roisIter.hasNext()) {
					final Roi roi = roisIter.next();
					if (roi.contains(x, y)) {
						roiCurrent = roi.getName(); // This selected
						roiN = Integer.parseInt(roiCurrent.substring(5));
					}
				}

				if (!(roiCurrent.equals(roiSelected))) {
					roiPreviouslySelected = roiSelected;
					roiSelected = roiCurrent;
					reDrawROIs(imp, roiCurrent);
					if (!roiPreviouslySelected.equals("none")) {
						final int roiPN = Integer.parseInt(roiPreviouslySelected.substring(
							5));
						plotter.removeVLine(roiPN);
					}
					plotter.plotsMontage();
					if (roiCurrent.equals("none") && imp.getRoi() != null) imp.killRoi();
				}
				else {
					if (!roiSelected.equals("none")) { // Moving inside same ROI
						Color lineColor;
						if (removePoint) lineColor = Plotter.vLineRemovePointColor;
						else if (addPoint) lineColor = Plotter.vLineAddPointColor;
						else if (!addPoint && !removePoint) lineColor =
							Plotter.vLineRegColor;
						else {
							lineColor = Color.white;
							System.out.println("Something wrong with the booleans");
						}
						plotter.setVLine(roiN, y, lineColor);
						plotter.plotsMontage();
					}
				}
			}
		}
	}

	@Override
	public void mouseClicked(final MouseEvent e) {
		final int x = ((ImageCanvas) e.getSource()).offScreenX(e.getX());
		final int y = ((ImageCanvas) e.getSource()).offScreenX(e.getY());
		if (!(addPoint || removePoint) && e.getClickCount() == 2) {
			if (!auto && !roiSelected.equals("none") && !selectionUpdate) {
				// Remove Roi
				final Roi roi = getRoi(roiSelected);
				if (roi == null || roi.isHandle(e.getX(), e.getY()) != -1) {
					// Might be dragging existing roi
					return;
				}

				if (askUser("Delete " + roiSelected + "?")) {
					if (fitDone && !askUser(warningFit))	return;
					final Iterator<Roi> roiIter = rois.iterator();
					while (roiIter.hasNext()) {
						final Roi roiRemove = roiIter.next();
						if (roiRemove.getName().equals(roiSelected)) {
							roiIter.remove();
							// System.out.println(imp.getTitle() + ": " + roiSelected +
							// " - REMOVED");
						}
					}
					reDrawROIs(imp, "none");
					redoProfilePlots();
				}
			}
		}
		else {
			final int lane = Integer.parseInt(roiSelected.substring(5));
			final int intensity = imp.getPixel(x, y)[0];
			boolean found = false;
			if (addPoint && !removePoint) {
				final double sd = askPoint(true, lane, y, intensity);
				for (final CustomPoints p : fitter.getCustomPoints()) {
					if (p.getLane() == lane) {
						p.addToAddList(y, intensity, sd);
						found = true; break;
					}
				}
				if (!found) {
					CustomPoints cp = new CustomPoints(lane);
					cp.addToAddList(y, intensity, sd);
					fitter.getCustomPoints().add(cp);
				}
			}
			else if (!addPoint && removePoint) {
				askPoint(false, lane, y, intensity);
				for (final CustomPoints p : fitter.getCustomPoints()) {
					if (p.getLane() == lane) {
						p.addToRemoveList(y, intensity);
						found = true; break;
					}
				}
				if (!found) {
					CustomPoints cp = new CustomPoints(lane);
					cp.addToRemoveList(y, intensity);
					fitter.getCustomPoints().add(cp);
				}
			}
		}
	}

	@Override
	public void mousePressed(final MouseEvent e) {
		if (!auto) selectionUpdate = true;
	}

	@Override
	public void mouseReleased(final MouseEvent e) {
		// Add or modify Roi
		if (buttonManual.isSelected()) {
			if (selectionUpdate) {
				if (fitDone && !askUser(warningFit))	return;
				final Roi roiNew = imp.getRoi();
				if (roiNew != null && roiNew.getBounds().getWidth() * roiNew.getBounds()
					.getHeight() > 100)
				{
					if (roiSelected.equals("none")) {
						// Creating New ROI
						// Find an available integer for lane name
						int i = 0;
						while (++i <= rois.size() + 1) {
							final Iterator<Roi> roisIter = rois.iterator();
							boolean contains = false;
							while (roisIter.hasNext()) {
								if (i == Integer.parseInt(roisIter.next().getName().substring(
									5)))
								{
									contains = true;
									break;
								}
							}
							if (!contains) {
								roiNew.setName("Lane " + i);
								roiSelected = roiNew.getName();
								break;
							}
						}
						rois.add(roiNew);
						// System.out.println(imp.getTitle() + ": " + roiSelected +
						// " - ADDED");
					}
					else {
						// System.out.println("Modify");
						// Moving/Resizing a specific ROI
						final Iterator<Roi> roisIter = rois.iterator();
						while (roisIter.hasNext()) {
							final Roi roi = roisIter.next();
							if (roi != null) {
								if (roi.getName().equals(roiSelected)) {
									roiNew.setName(roi.getName());
									rois.set(rois.indexOf(roi), roiNew);
									break;
								}
							}
						}
					}
					reDrawROIs(imp, roiSelected);
					redoProfilePlots();
				}
			}
		}
		selectionUpdate = false;
	}

	@Override
	public void mouseEntered(final MouseEvent e) {
		// Nothing to do
		return;
	}

	@Override
	public void mouseExited(final MouseEvent e) {
		// Nothing to do
		return;
	}

	@Override
	public void mouseWheelMoved(final MouseWheelEvent e) {
		if (e.getSource() == imp.getCanvas() && imp.getImageStackSize() > 1) {
			if (fitDone && !askUser(warningFit)) return;
			redoProfilePlots();
		}
	}

	public void setFitter(final Fitter fitter) {
		this.fitter = fitter;
	}
}
