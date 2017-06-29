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
import java.awt.Checkbox;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
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
import javax.swing.text.Document;

import org.apache.commons.math3.analysis.function.Abs;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;
import org.scijava.Context;
import org.scijava.app.StatusService;
import org.scijava.display.Display;
import org.scijava.display.DisplayService;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.Line;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.gui.TextRoi;

@SuppressWarnings("serial") 
class MainDialog extends JFrame implements ActionListener,
	ChangeListener, DocumentListener, MouseMotionListener,
	MouseListener, MouseWheelListener, WindowListener
{

	@Parameter
	private LogService log;

	@Parameter
	private DisplayService displayServ;

	@Parameter
	private StatusService statusServ;

	private static final Color guess = Color.BLUE;
	private static final Color custom = Color.GREEN;
	private static final Color fit = Color.MAGENTA;
	private static final Color selected = Color.YELLOW;
	private static final Color unselected = Color.RED;
	private static final Color reference = Color.BLUE;

	private boolean auto = true; // AUTO Lane size mode is ON
	private boolean selectionUpdate = false; // active updating is off
	private boolean fitDone = false; // Keep track of whether fit data exists
	private boolean addPeak = false;
	private boolean removePeak = false;

	private String roiSelected = "none";
	private String roiPreviouslySelected = "none";
	private String roiReference = "none";
	private String refLabel = "Use Reference Ladder";

	private final String warningFit =
		"The current plots will be reset and the current fitting data will be lost.";

	private Plotter plotter;
	private Fitter fitter;

	private final ImagePlus imp;
	private ArrayList<Roi> rois;

	private int IW, IH, LW, LH, LSp, LHOff, LVOff;

	private int nLanes = 16;
	private int degBG = 2; // Order of Background Polynomial
	private double tolPK = 0.2; // Peak detection tolerance as % of range

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
	private JButton buttonFit;
	private JButton buttonClose;

	private JPanel settingsPanel;
	private JPanel degPanel;
	private JPanel tolPanel;
	private JLabel labelDegBG;
	private JLabel labelTolPK;
	private JTextField textDegBG;
	private JTextField textTolPK;
	private JCheckBox chkBoxBands;
	private JCheckBox chkBoxRef;
	private JToggleButton buttonAddPeak;
	private JToggleButton buttonRemovePeak;
	private JButton buttonResetCustomPeaks;

	private JPanel dialogPanel;
	private final JFrame frame;

	MainDialog(final Context context, final String string, final ImagePlus imp, double screenHeight, double screenWidth) {
		context.inject(this);
		frame = new JFrame(string);
		rois = new ArrayList<>();
		this.imp = imp;
		setupMainDialog();
		frame.setLocation(0, (int) screenHeight/8);
		frame.setSize((int) screenWidth/3, (int) screenHeight/2);
	}

	private void setupMainDialog() {
		// Default lane size/offset (Just center 4 lanes in the image)
		IW = imp.getWidth();
		IH = imp.getHeight();
		 LW = 36;
		 LH = 866;
		 LSp = 4;
		 LHOff = 132;
		 LVOff = 40;
//		LW = (int) Math.round(0.8 * IW / nLanes);
//		LH = (int) Math.round(IH * 0.8);
//		LSp = Math.round((IW - LW * nLanes) / (nLanes + 1));
//		LHOff = LSp / 2;
//		LVOff = (IH - LH) / 2;

		imp.getCanvas().addMouseMotionListener(this);
		imp.getCanvas().addMouseListener(this);
		imp.getCanvas().addMouseWheelListener(this);
		imp.getWindow().addMouseListener(this);
		imp.getWindow().addMouseMotionListener(this);
		imp.getWindow().addWindowListener(this);

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

		sliderPanel.setLayout(new BoxLayout(sliderPanel, BoxLayout.Y_AXIS));
		sliderW = makeTitledSlider("Width ( " + LW + " px )", Color.black, 1, IW /
			nLanes, LW);
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
		buttonFit = new JButton("Fit");
		buttonFit.addActionListener(this);
		buttonPanel.add(buttonFit);
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
		textDegBG.getDocument().addDocumentListener(this);
		textTolPK = new JTextField(3);
		textTolPK.setText("" + tolPK);
		textTolPK.getDocument().addDocumentListener(this);
		chkBoxBands = new JCheckBox("Show Bands");
		chkBoxRef = new JCheckBox(refLabel);
		buttonAddPeak = new JToggleButton("Add Peak");
		buttonRemovePeak = new JToggleButton("Remove Peak");
		buttonResetCustomPeaks = new JButton("Reset Custom Peaks");
		
		degPanel.setLayout(new FlowLayout(FlowLayout.TRAILING));
		degPanel.add(labelDegBG);
		degPanel.add(textDegBG);
		tolPanel.setLayout(new FlowLayout(FlowLayout.TRAILING));
		tolPanel.add(labelTolPK);
		tolPanel.add(textTolPK);
		chkBoxBands.addActionListener(this);
		chkBoxBands.setSelected(false);
		chkBoxBands.setEnabled(false);
		chkBoxRef.addActionListener(this);
		chkBoxRef.setSelected(false);
		chkBoxRef.setEnabled(false);
		buttonAddPeak.addActionListener(this);
		buttonAddPeak.setEnabled(false);
		buttonRemovePeak.addActionListener(this);
		buttonRemovePeak.setEnabled(false);
		buttonResetCustomPeaks.addActionListener(this);
		buttonResetCustomPeaks.setEnabled(false);

		settingsPanel.setLayout(new GridLayout(10, 1));
		settingsPanel.add(degPanel);
		settingsPanel.add(tolPanel);
		settingsPanel.add(chkBoxBands);
		settingsPanel.add(chkBoxRef);
		settingsPanel.add(buttonAddPeak);
		settingsPanel.add(buttonRemovePeak);
		settingsPanel.add(buttonResetCustomPeaks);

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
		frame.addWindowListener(new WindowAdapter() {

			@Override
			public void windowClosing(final WindowEvent e) {
				cleanupAndClose();
			}
		});

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

	private double askPeak(final boolean add, final int lane, double y,
		double a)
	{
		String title, action, message, fromTo;
		final String fwhmString = "FWHM = \u03C3 * (2 * \u221A (2 * ln(2))";
		double fwhmValue = 0.0;

		if (add) {
			title = "ADD PEAK";
			action = "add";
			message = "If the new peak is located less than " +
				(int) Fitter.peakDistanceTol +
				" px away from an exisitng peak,\nthe existing peak will be replaced with the new one.";
			fromTo = "to the custom list";
		}
		else {
			title = "REMOVE PEAK";
			action = "remove";
			message =
				"The custom peak that is closest\n to this peak will be removed";
			fromTo = "from the custom list.";
		}

		boolean foundGuess = false;
		boolean foundCustom = false;
		for (final Peak p : fitter.getGuessPeaks(lane)) {
			if (FastMath.abs(y - p.getMean()) <= Fitter.peakDistanceTol) {
				y = p.getMean();
				a = p.getNorm();
				fwhmValue = p.getSigma();
				foundGuess = true;
				break;
			}
		}
		for (final Peak p : fitter.getCustomPeaks(lane)) {
			if (FastMath.abs(y - p.getMean()) <= Fitter.peakDistanceTol) {
				y = p.getMean();
				a = p.getNorm();
				fwhmValue = p.getSigma();
				foundCustom = true;
				break;
			}
		}

		if (add && (foundCustom || foundGuess)) {
			action = "replace";
			if (foundCustom && foundGuess) {
				fromTo = " in the custom list and the guesslist.";
			}
			else if (foundCustom && !foundGuess) {
				fromTo = " in the custom list.";
			}
			else if (!foundCustom && foundGuess) {
				fromTo = " in the guesslist.";
			}
		}
		if (!foundCustom) fwhmValue = 5.0;
		final GenericDialog gd = new GenericDialog(title);
		if (add || (!add && foundCustom)) {
			gd.addMessage("You are about to " + action + " this peak\n" + fromTo);
			gd.addMessage(message);
			gd.addMessage(String.format("Lane: \t%2d", lane));
			gd.addMessage(String.format("Distance: \t%10.1f", y));
			gd.addMessage(String.format("Intensity: \t%.0f", a));
		}
		else if (!add && !foundCustom) {
			gd.addMessage("No custom peak found nearby. Try again!");
		}
		if (add) {
			gd.addMessage("Estimate the full width at half maximum,\n" + fwhmString);
			gd.addNumericField("FWHM", fwhmValue, 1);
		}
		gd.showDialog();
		if (gd.wasOKed()) {
			if (add) return gd.getNextNumber();
			if (!add && !foundCustom) return 0.0;
			return -1.0; // Removing peak
		}
		return 0.0; // If Dialog Cancelled
	}

	@SuppressWarnings("unchecked")
	private int[] askResetCustomPoints(final boolean[] defaultValues) {
		String title = "RESET CUSTOM PEAKS";
		String message = "Select the plots from which the custom peaks must be removed";
		final GenericDialog gd = new GenericDialog(title);
		gd.addMessage(message);
		final int[] lanes = getAllLaneNumbers();
		final int rows = lanes.length + 1;
		ArrayList<Integer> selectedLanes = new ArrayList<>();
		final String[] labels = new String[rows];
		for (final int i : lanes) {
			labels[i - 1] = "Lane " + i;
			defaultValues[i - 1] = false;
		}
		labels[labels.length - 1] = "All Lanes";
		defaultValues[defaultValues.length - 1] = false;
		gd.addCheckboxGroup(rows, 1, labels, defaultValues);
		
		gd.showDialog();
		if (gd.wasOKed()) {
			Vector<Checkbox> v = gd.getCheckboxes();
			for (Checkbox cb : v) {
				if (cb.getState()) {
					selectedLanes.add(v.indexOf(cb) + 1);
				}
			}
			if (selectedLanes.contains(rows)){
				// 'All Lanes' is selected
				if (selectedLanes.size() != 1) { // Not just 'All lanes'
					boolean[] b = new boolean[rows];
					b[b.length - 1] = true;
					askResetCustomPoints(b);
				}
				return lanes; 
			} 
			int[] intArray = new int[selectedLanes.size()];
			for (int i=0; i<intArray.length; i++) {
				intArray[i] = selectedLanes.get(i);
			}
			return intArray;
		}
		if(gd.wasCanceled()) {
			return new int[0];
		}
		return new int[0];
	}
	
	private void askReferenceLane() {
		String title = "REFERENCE LANE";
		String message = "Select the lane to be used as reference.";
		final GenericDialog gd = new GenericDialog(title);
		gd.addMessage(message);
		final int[] lanes = getAllLaneNumbers();
		int rows = lanes.length;
		final boolean[] defaultValues = new boolean[rows];
		final String[] labels = new String[rows];
		for (final int i : lanes) {
			labels[i - 1] = "Lane " + i;
			defaultValues[i - 1] = false;
		}
		
		gd.addRadioButtonGroup(null, labels, rows, 1, labels[0]);
		gd.showDialog();
		if (gd.wasOKed()) {
			roiReference = gd.getNextRadioButton();
		}
		if(gd.wasCanceled()) {
			roiReference = "none";
		}
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
			// Housekeeping: Sort the rois based on name and remove null elements
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
			final int ln = Integer.parseInt(roi.getName().substring(5));
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
				roi.setStrokeColor(MainDialog.selected);
				roi.setStrokeWidth(3);
				if (buttonManual.isSelected() && !roiName.equals("none")) {
					imgPlus.setRoi(roi);
				}
			} else if (roi.getName().equals(roiReference)) {
					if (!roiReference.equals("none") && !roiReference.equals(roiName)) {
						roi.setStrokeWidth(3);
						roi.setStrokeColor(MainDialog.reference);
						labelRoi.setFillColor(MainDialog.reference);
					}
			} else {
				roi.setStrokeWidth(1);
				roi.setStrokeColor(MainDialog.unselected);
				labelRoi.setFillColor(MainDialog.unselected);
			}
			overlay.add(roi);
			overlay.add(labelRoi);
			if (buttonAuto.isSelected()) imgPlus.killRoi();
			if (chkBoxBands.isSelected() && fitDone) {
				
				if (!chkBoxRef.isSelected() || (chkBoxRef.isSelected() && roi.getName().equals(roiReference))) {
					// Draw a tick where the bands are
					for (final Peak p : fitter.getFittedPeaks(ln)) {
						final double y = p.getMean();
						final Roi band1 = new Line(x0, y, x0 + 10, y);
						band1.setStrokeColor(MainDialog.fit);
						band1.setStrokeWidth(1);
						overlay.add(band1);
						final Roi band2 = new Line(x0 + rw - 10, y, x0 + rw, y);
						band2.setStrokeColor(MainDialog.fit);
						band2.setStrokeWidth(1);
						overlay.add(band2);
					}
					for (final Peak p : fitter.getGuessPeaks(ln)) {
						final double y = p.getMean();
						final Roi band1 = new Line(x0, y, x0 + 5, y);
						band1.setStrokeColor(MainDialog.guess);
						band1.setStrokeWidth(1);
						overlay.add(band1);
						final Roi band2 = new Line(x0 + rw - 5, y, x0 + rw, y);
						band2.setStrokeColor(MainDialog.guess);
						band2.setStrokeWidth(1);
						overlay.add(band2);
					}
					for (final Peak p : fitter.getCustomPeaks(ln)) {
						final double y = p.getMean();
						final Roi band1 = new Line(x0, y, x0 + 5, y);
						band1.setStrokeColor(MainDialog.custom);
						band1.setStrokeWidth(1);
						overlay.add(band1);
						final Roi band2 = new Line(x0 + rw - 5, y, x0 + rw, y);
						band2.setStrokeColor(MainDialog.custom);
						band2.setStrokeWidth(1);
						overlay.add(band2);
					}
				}
			}
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

	private void changeTextFieldVariable(final DocumentEvent e) {
		final Document textBox = e.getDocument();
		if (textBox.equals(textNLanes.getDocument())) {
			if (fitDone && !askUser(warningFit)) {
				textNLanes.setText(Integer.toString(nLanes));
				return;
			}
			nLanes = getNLanes();
			if (nLanes == -1) return;
			resetAutoROIs(LW, LH, LSp, LHOff, LVOff, nLanes);
			reDrawROIs(imp, "none");
			redoProfilePlots();
			fitter.resetAllFitter();
			fitter.setInputData(plotter.getProfiles());
		}
		else if (textBox.equals(textDegBG.getDocument())) {
			fitter.setDegBG(getDegBG());
		}
		else if (textBox.equals(textTolPK.getDocument())) {
			fitter.setTolPK(getTolPK());
		}
	}

	private void redoProfilePlots() {
		fitDone = false;
		chkBoxBands.setSelected(false);
		chkBoxBands.setEnabled(false);
		buttonAddPeak.setEnabled(false);
		buttonRemovePeak.setEnabled(false);
		buttonResetCustomPeaks.setEnabled(false);
		plotter.resetData();
		plotter.setMode(Plotter.regMode);
		for (final Roi r : rois) {
			final Rectangle rect = r.getBounds();
			if (rect.getMinX() < 0.95 * IW && rect.getMinY() < 0.95 * IH) plotter
				.updateProfile(r);
		}
	}

	private void cleanupAndClose() {
		// Remove Listeners on gel image and close plugin-associated windows in
		// preparation for exit
		if (imp != null) {
			imp.getCanvas().removeMouseListener(this);
			imp.getCanvas().removeMouseMotionListener(this);
			imp.getCanvas().removeMouseWheelListener(this);
			imp.getWindow().removeMouseListener(this);
			imp.getWindow().removeMouseMotionListener(this);
			imp.setOverlay(null);
		}
		plotter.closePlot();
		for (final Display<?> d : displayServ.getDisplays())
			d.close();
		log.info("\nGel Lanes Fit terminated.");
	}

	/**
	 * @return @param nLanes, the number of ROIs currently present in the gel
	 *         image
	 */
	public int getNLanes() {
		try {
			return Integer.parseInt(textNLanes.getText().trim());
		}
		catch (final NumberFormatException e1) {
			return -1;
		}
	}

	/**
	 * @return the currently selected @param degBG, the degree of the polynomial
	 *         representing the background signal
	 */
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
	 * @return the currently selected @param tolPK, the tolerance on the peak
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

	/**
	 * @return ROI identified by title. Typically "Lane n", where int n > 0
	 */
	private Roi getRoi(final String title) {
		final Iterator<Roi> roiIter = rois.iterator();
		while (roiIter.hasNext()) {
			final Roi roi = roiIter.next();
			if (roi.getName().equals(title)) return roi;
		}
		return null;
	}

	/**
	 * @return @param rois, the current ROI set
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
		return Arrays.copyOfRange(laneNumbers, 0, i);
	}

	public String getRoiPreviouslySelected() {
		return roiPreviouslySelected;
	}

	public void setFitter(final Fitter fitter) {
		this.fitter = fitter;
	}

	public void setPlotter(final Plotter plotter) {
		this.plotter = plotter;
	}

	@Override
	public synchronized void stateChanged(final ChangeEvent e) {
		if (fitDone && !askUser(warningFit)) return;
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
		chkBoxRef.setSelected(false);
		redoProfilePlots();
		fitter.resetAllFitter();
		fitter.setInputData(plotter.getProfiles());

		// Close results table
		for (final Display<?> d : displayServ.getDisplays()) {
			if (d.getName().equals("Results Display")) d.close();
		}
	}

	// TextField Document Listeners (change, remove, insert)
	@Override
	public void changedUpdate(final DocumentEvent e) {
		changeTextFieldVariable(e);
	}

	@Override
	public void removeUpdate(final DocumentEvent e) {
		changeTextFieldVariable(e);
	}

	@Override
	public void insertUpdate(final DocumentEvent e) {
		changeTextFieldVariable(e);
	}

	@Override
	public void actionPerformed(final ActionEvent e) {
		if (e.getSource().equals(buttonFit)) {
			if (fitDone && !askUser(warningFit)) return;
			
			// Close results table
			// TODO: Chose between these 2 methods
			if (displayServ.getDisplay("Results Display") != null) displayServ
				.getDisplay("Results Display").close();
			
			for (final Display<?> d : displayServ.getDisplays()) {
				if (d.getName().equals("Results Display")) d.close();
			}
			
			addPeak = false;
			removePeak = false;
			buttonAddPeak.setSelected(false);
			buttonRemovePeak.setSelected(false);
			buttonResetCustomPeaks.setSelected(false);

			// Remove everything in the plots, but the profile
			plotter.setSelected(Plotter.selectedNone);
			plotter.removeFit();
			plotter.setMode(Plotter.regMode);

			// Reset fit
			if (fitDone) {
				fitter.resetFit();
			}
			else {
				fitter.setInputData(plotter.getProfiles());
			}

			final ArrayList<DataSeries> fitted = fitter.doFit();
			fitDone = true;
			plotter.addDataSeries(fitted);
			for (int i :getAllLaneNumbers())
				plotter.updatePlot(i);
			
			chkBoxBands.setEnabled(true);
			chkBoxBands.setSelected(true);
			chkBoxRef.setEnabled(true);
			buttonAddPeak.setEnabled(true);
			buttonRemovePeak.setEnabled(true);
			buttonResetCustomPeaks.setEnabled(true);
			reDrawROIs(imp, "none"); // adds the bands
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

		if (e.getSource().equals(buttonAddPeak)) {
			addPeak = !addPeak;
			removePeak = false;
			buttonRemovePeak.setSelected(false);
			if (addPeak) plotter.setMode(Plotter.addMode);
			else plotter.setMode(Plotter.regMode);
			
			}

		if (e.getSource().equals(buttonRemovePeak)) {
			addPeak = false;
			removePeak = !removePeak;
			buttonAddPeak.setSelected(false);
			if (removePeak) plotter.setMode(Plotter.remMode);
			else plotter.setMode(Plotter.regMode);
		}
		
		
		if (e.getSource().equals(buttonResetCustomPeaks)) {
			int[] lanes = askResetCustomPoints(new boolean[rois.size() +1]);
			for (final int i : lanes) {
				fitter.resetCustomPeaks(i);
				plotter.updateCustomPeaks(i, fitter.getCustomPeaks(i));
				plotter.updatePlot(i);
			}
			reDrawROIs(imp, "none");
		}
		
		// CheckBoxes ------------------------------------------------------------------------------------
		if (e.getSource().equals(chkBoxBands) && chkBoxBands.isSelected() && fitDone) {
			reDrawROIs(imp, "none");
		}
		
		if (e.getSource().equals(chkBoxRef) && chkBoxRef.isSelected() && fitDone) {
			//chkBoxRef.setSelected(chkBoxRef.isSelected());
			askReferenceLane();
			if (!roiReference.equals("none")) {
				int refLn = Integer.parseInt(roiReference.substring(5));
				roiReference = "Lane " + refLn;
				String label = refLabel + ": " + roiReference;
				chkBoxRef.setText(label);
				plotter.setReferencePlot(refLn, fitter.getFittedPeaks(refLn));
				reDrawROIs(imp, "none");
				for (int i : getAllLaneNumbers()){
					plotter.updatePlot(i);
				}
			} else { // Cancelled
				chkBoxRef.setSelected(false);
				chkBoxRef.setText(refLabel);
			}	
		} else {
			roiReference = "none";
			chkBoxRef.setText(refLabel);
			ArrayList<VerticalMarker> plotsMarkers = plotter.getPlotsVerticalMarkers();
			Iterator<VerticalMarker> markIter = plotsMarkers.iterator();
			while (markIter.hasNext()) {
				if (markIter.next().getType() == VerticalMarker.BMARK)
					markIter.remove();
			}
			for (int i : getAllLaneNumbers())
				plotter.updatePlot(i);
		}

		
		if (e.getSource().equals(buttonClose)) {
			if (askUser("Would you like to quit Gel Lanes Fit?")) {
				frame.dispatchEvent(new WindowEvent(frame, WindowEvent.WINDOW_CLOSING));
			}
		}
	}

	@Override
	public void mouseDragged(final MouseEvent e) {
		// Nothing to do
	}

	@Override
	public void mouseMoved(final MouseEvent e) {
		if (e.getSource() == imp.getCanvas()) {
			if (rois.size() == 0 || plotter == null) return;
			final int x = ((ImageCanvas) e.getSource()).offScreenX(e.getX());
			final int y = ((ImageCanvas) e.getSource()).offScreenX(e.getY());
			statusServ.showStatus("[" + x + ":" + y + "]");
			if (!selectionUpdate) { // Not dragging an ROI
				String roiCurrent = "none"; // None selected
				for (Roi r : rois) {
					if (r.contains(x, y)) {
						roiCurrent = r.getName(); // This selected
					}
				}

				reDrawROIs(imp, roiCurrent);
				if (roiCurrent.equals("none") && imp.getRoi() != null)
					imp.killRoi();

				if (!(roiCurrent.equals(roiSelected))) {
					roiPreviouslySelected = roiSelected;
					roiSelected = roiCurrent;
					if (!roiCurrent.equals("none")) {
						int ln = Integer.parseInt(roiCurrent.substring(5));
						plotter.setSelected(ln);
						plotter.setVLine(ln, y);
					} else {
						plotter.setSelected(Plotter.selectedNone);
					}

					if (!roiPreviouslySelected.equals("none")) {
						final int roiPS = Integer.parseInt(roiPreviouslySelected.substring(
							5));
						plotter.removeVLine(roiPS);
					}
				} else { // Moving inside same ROI
					if (!roiCurrent.equals("none")) {
						int ln = Integer.parseInt(roiCurrent.substring(5));
						plotter.setVLine(ln, y);
					}
				}
			}
		}
	}

	@Override
	public void mouseClicked(final MouseEvent e) {
		if (e.getSource().equals(imp.getCanvas())) {
			// final int x = imp.getCanvas().offScreenX(e.getX());
			final int y = imp.getCanvas().offScreenX(e.getY());
			if (!(addPeak || removePeak) && e.getClickCount() == 2) {
				if (!auto && !roiSelected.equals("none") && !selectionUpdate) {
					// Remove Roi
					final Roi roi = getRoi(roiSelected);
					if (roi == null || roi.isHandle(e.getX(), e.getY()) != -1) {
						// Might be dragging existing roi
						return;
					}

					if (askUser("Delete " + roiSelected + "?")) {
						if (fitDone && !askUser(warningFit)) return;
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
			else { // ADD/REMOVE custom peak is on
				if (!roiSelected.equals("none")) {
					final int lane = Integer.parseInt(roiSelected.substring(5));
					double intensity = 0;
					for (final DataSeries d : plotter.getProfiles()) {
						if (d.getLane() == lane) {
							final RealVector xvals = new ArrayRealVector(d.getX());
							final RealVector yvals = new ArrayRealVector(d.getY());
							// Value on profile closest to the clicked y value in imp
							intensity = yvals.getEntry(xvals.mapSubtract(y).mapToSelf(
								new Abs()).getMinIndex());
						}
					}
					// Fitter has 1 CustomPeaks object for each plot/lane
					if (addPeak && !removePeak) { // Add Peak
						final double sd = askPeak(true, lane, y, intensity) / (2 * FastMath
							.sqrt(2 * FastMath.log(2)));
						if (sd != 0.0) { // not cancelled
							fitter.addCustomPeak(lane, new Peak(lane, intensity, y, sd));
						}
					}
					else if (!addPeak && removePeak) { // Remove Peak
						// Remove points previously added MANUALLY, not implemented
						// Does not make sense to remove points that where detected
						// automatically.
						// Better to increase threshold to avoid detecting noise.
						// Use remove feature to remove a previously added custom peak.
						final double sd = askPeak(false, lane, y, intensity);
						if (sd != 0.0) {
							fitter.removeCustomPeak(lane, new Peak(lane, intensity, y));
						}
					}
					reDrawROIs(imp, "none");
					plotter.updateCustomPeaks(lane, fitter.getCustomPeaks(lane));
					plotter.updatePlot(lane);
				}
			}
		}
	}

	@Override
	public void mousePressed(final MouseEvent e) {
		if (!(auto || addPeak || removePeak)) selectionUpdate = true;
	}

	@Override
	public void mouseReleased(final MouseEvent e) {
		// Modify Roi, mouseReleased not triggered when creating Roi
		if (buttonManual.isSelected()) {
			if (selectionUpdate) {
				final Roi roiNew = imp.getRoi();
				if (roiNew != null) {
					// If ROI exists and larger that 100 px
					final Rectangle rN = roiNew.getBounds();
					if (rN.getWidth() * rN.getHeight() > 100) {
						if (roiNew.getName() == null) { // It's a new ROI
							// log.info("Add Roi");
							// Creating New ROI
							// Find an available integer for lane name
							boolean inRoiList = false;
							int i = 0;
							while (++i <= rois.size() + 1) {
								for (final Roi r : rois) {
									if (i == Integer.parseInt(r.getName().substring(5))) {
										inRoiList = true;
										break;
									}
								}
								if (!inRoiList) {
									roiNew.setName("Lane " + i);
									break;
								}
								inRoiList = false;
							}
						}
						else {
							// log.info("Modify Roi");
							// Moving/Resizing a specific ROI
							final Iterator<Roi> roisIter = rois.iterator();
							while (roisIter.hasNext()) {
								final Roi r = roisIter.next();
								if (r.getName().equals(roiNew.getName())) {
									roisIter.remove();
									break;
								}
							}
						}
						rois.add(roiNew);
						roiSelected = roiNew.getName();
						reDrawROIs(imp, roiSelected);
						redoProfilePlots();
					}
				}
				selectionUpdate = false;
			}
		}
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

	@Override
	public void windowOpened(final WindowEvent e) {
		// Nothing to do
	}

	@Override
	public void windowClosing(final WindowEvent e) {
		if (e.getSource() == imp.getWindow()) {
			frame.dispatchEvent(new WindowEvent(frame, WindowEvent.WINDOW_CLOSING));
		}
	}

	@Override
	public void windowClosed(final WindowEvent e) {
		// Nothing to do
	}

	@Override
	public void windowIconified(final WindowEvent e) {
		// Nothing to do
	}

	@Override
	public void windowDeiconified(final WindowEvent e) {
		// Nothing to do
	}

	@Override
	public void windowActivated(final WindowEvent e) {
		// Nothing to do
	}

	@Override
	public void windowDeactivated(final WindowEvent e) {
		// Nothing to do
	}
}
