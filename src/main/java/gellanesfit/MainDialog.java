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

import java.awt.BorderLayout;
import java.awt.Checkbox;
import java.awt.Color;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
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
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;
import java.util.prefs.Preferences;

import javax.swing.AbstractButton;
import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSlider;
import javax.swing.JSpinner;
import javax.swing.JToggleButton;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingConstants;
import javax.swing.border.EmptyBorder;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.text.Document;

import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.jfree.data.general.SeriesChangeEvent;
import org.jfree.data.general.SeriesChangeListener;
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

class MainDialog extends JFrame implements ActionListener, ChangeListener, SeriesChangeListener,
        MouseMotionListener, MouseListener, MouseWheelListener, WindowListener {

	@Parameter
	private LogService log;

	@Parameter
	private DisplayService displayServ;

	@Parameter
	private StatusService statusServ;

	// Preference keys for this class
	private static final String DEGBG = "degBG";
	private static final String POLYCONCAVITY = "polyConcavity";
	private static final String TOLPK = "tolPK";
	private static final String NORMDRIFT = "normDrift";
	private static final String SDDRIFT = "sdDrift";

	private static final String NLANES = "nLanes";
	private static final String LW = "lw";
	private static final String LH = "lh";
	private static final String LSP = "lsp";
	private static final String LHOFF = "lhoff";
	private static final String LVOFF = "lvoff";
	private static final String AUTO = "auto";

	private final double SW = IJ.getScreenSize().getWidth();
	private final double SH = IJ.getScreenSize().getHeight();

	private static final Color guess = Color.BLUE;
	private static final Color custom = Color.GREEN;
	private static final Color fit = Color.MAGENTA;
	private static final Color selected = Color.YELLOW;
	private static final Color unselected = Color.RED;
	private static final Color reference = Color.BLUE;
	
	public static final int noLadderLane = -1;
	public static final int noLaneSelected = -1;

	private boolean auto; // AUTO ROI mode; all lanes same size
	private boolean selectionUpdate = false; // active updating is off
	private boolean fitDone = false; // Keep track of whether fit data exists

	private String roiSelected = "none";
	private String roiPreviouslySelected = "none";
	private String roiReferenceName = "none";
	private int roiReferenceNumber = noLadderLane;
	private final String[] ladderStr = { "Select Ladder Type", "Hi-Lo", "100bp" };
	private final String[] distStr = { "Select Fragment Distribution", "AciI-Lambda", "Uniform" };
	
	// Ladders
	private static String[] hilo = { "10 kbp", "8 kbp", "6 kbp", "4 kbp", "3 kbp", "2 kbp",
	        "1.55 kbp", "1.4 kbp", "1 kbp", "750 bp", "500 bp", "400 bp", "300 bp", "200 bp",
	        "100 bp", "50 bp" };
	private static String[] bp100 = { "1 kbp", "900 bp", "800 bp", "700 bp", "600 bp", "500 bp",
	        "400 bp", "300 bp", "200 bp", "100 bp" };
	private final int[] ladderRange = new int[2];
	private final String warningFit = "The current plots will be reset and the current fitting data will be lost.";

	private Plotter plotter;
	private Fitter fitter;
	private final Preferences prefs;

	private final ImagePlus imp;
	private ArrayList<Roi> rois;

	private int iw, ih, lw, lh, lsp, lhoff, lvoff;

	private int nLanes;
	private int degBG; // Order of Background Polynomial
	private double polyConcavity;
	private double tolPK; // Peak detection tolerance as % of range
	private double normDrift;
	private double sdDrift;
	
	private JPanel buttonPanelAutoManual;
	private JRadioButton buttonAuto;
	private JRadioButton buttonManual;
	private ButtonGroup buttonGroupAutoManual;

	private JPanel lanesPanel;
	private JLabel labelNLanes;
	private JSpinner textNLanes;

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
	private JLabel labelDegBG;
	private JLabel labelPolyConcavity;
	private JLabel labelTolPK;
	private JLabel labelNormDrift;
	private JLabel labelSDDrift;

	private JSpinner textDegBG;
	private JSpinner textPolyConcavity;
	private JSpinner textTolPK;
	private JSpinner textNormDrift;
	private JSpinner textSDDrift;

	private JPanel buttonPanelFitType;
	private JRadioButton buttonBands;
	private JRadioButton buttonContinuum;
	private ButtonGroup buttonGroupFitType;

	private JCheckBox chkBoxBands;
	private JComboBox<String> cmbBoxLadderLane;
	private JComboBox<String> cmbBoxLadderType;
	private JComboBox<String> cmbBoxDist;
	private JToggleButton buttonEditPeaks;
	private JButton buttonResetCustomPeaks;

	private JPanel dialogPanel;
	private final JFrame frame;

	MainDialog(final Context context, final String string, final ImagePlus imp,
	        final Preferences prefs) {
		context.inject(this);
		this.prefs = prefs;
		frame = new JFrame(string);
		rois = new ArrayList<>();
		this.imp = imp;
		setupMainDialog();
		frame.setLocation((int) SW / 16, (int) SH / 8);
	}

	private void setupMainDialog() {
		int textWidth = 4;
		
		// Default lane size/offset
		auto = prefs.getBoolean(AUTO, true);
		degBG = prefs.getInt(DEGBG, 2);
		polyConcavity = prefs.getDouble(POLYCONCAVITY, 10.0);
		tolPK = prefs.getDouble(TOLPK, 0.1);
		normDrift = prefs.getDouble(NORMDRIFT, 0.1);
		sdDrift = prefs.getDouble(SDDRIFT, 0.1);
		nLanes = prefs.getInt(NLANES, 4);
		
		iw = imp.getWidth();
		ih = imp.getHeight();
		lw = prefs.getInt(LW, (int) Math.round(0.8 * iw / nLanes));
		lh = prefs.getInt(LH, (int) Math.round(ih * 0.8));
		lsp = prefs.getInt(LSP, Math.round((iw - lw * nLanes) / (nLanes + 1)));
		lhoff = prefs.getInt(LHOFF, lsp / 2);
		lvoff = prefs.getInt(LVOFF, (ih - lh) / 2);
		

		
		buttonPanelAutoManual = new JPanel();
		buttonAuto = new JRadioButton("Automatic Rectangle Selection");
		buttonAuto.setActionCommand("Auto");
		buttonAuto.addActionListener(this);
		buttonManual = new JRadioButton("Manual Rectangle Selection");
		buttonManual.setActionCommand("Manual");
		buttonManual.addActionListener(this);
		buttonGroupAutoManual = new ButtonGroup();
		buttonPanelAutoManual.setLayout(new BoxLayout(buttonPanelAutoManual, BoxLayout.Y_AXIS));
		buttonGroupAutoManual.add(buttonAuto);
		buttonGroupAutoManual.add(buttonManual);
		
		buttonPanelAutoManual.add(buttonAuto);
		buttonPanelAutoManual.add(buttonManual);
		buttonPanelAutoManual.setBorder(new TitledBorder(BorderFactory.createEtchedBorder(),
		        "Lane Selection", TitledBorder.LEADING, TitledBorder.BELOW_TOP,
		        new Font("Sans", Font.PLAIN, 11)));
		
		sliderPanel = new JPanel();
		lanesPanel = new JPanel();
		lanesPanel.setLayout(new GridBagLayout());
		GridBagConstraints c0 = new GridBagConstraints();
		lanesPanel.setBorder(new EmptyBorder(3, 5, 3, 3));
		labelNLanes = new JLabel("Number of Lanes");
		textNLanes = new JSpinner(new SpinnerNumberModel(nLanes, 1, 100, 1));
		textNLanes.setBorder(BorderFactory.createCompoundBorder(
			textNLanes.getBorder(), BorderFactory.createEmptyBorder(0, 2, 0, 2)));
		((JSpinner.DefaultEditor) textNLanes.getEditor()).getTextField().setColumns(textWidth);
		
		c0.gridx = 0; c0.gridy = 0; c0.weightx = 0.2;
		c0.fill = GridBagConstraints.HORIZONTAL;
		lanesPanel.add(labelNLanes, c0);
		
		c0.gridx = 1; c0.gridy = 0; c0.weightx = 0.0;
		c0.fill = GridBagConstraints.NONE;
		c0.ipady = 2; c0.ipadx = 2;		
		lanesPanel.add(textNLanes, c0);
		sliderPanel.add(lanesPanel);
		
		sliderPanel.setLayout(new BoxLayout(sliderPanel, BoxLayout.Y_AXIS));
		sliderW = makeTitledSlider("Width ( " + lw + " px )", Color.black, 1, iw / nLanes, lw);
		sliderPanel.add(sliderW);
		sliderH = makeTitledSlider("Height ( " + lh + " px )", Color.black, ih / 10, ih, lh);
		sliderPanel.add(sliderH);
		sliderSp = makeTitledSlider("Space ( " + lsp + " px )", Color.black, 1, iw / nLanes, lsp);
		sliderPanel.add(sliderSp);
		sliderHOff = makeTitledSlider("Horizontal Offset ( " + lhoff + " px )", Color.black, 0,
		        (int) Math.round(iw * 0.9), lhoff);
		sliderPanel.add(sliderHOff);
		sliderVOff = makeTitledSlider("Vertical Offset ( " + lvoff + " px )", Color.black, 0,
		        (int) Math.round(ih * 0.9), lvoff);
		sliderPanel.add(sliderVOff);
		
		buttonPanel = new JPanel();
		buttonPanel.setLayout(new FlowLayout(FlowLayout.TRAILING));
		buttonFit = new JButton("Fit");
		
		buttonPanel.add(buttonFit);
		buttonClose = new JButton("Close");
		buttonPanel.add(buttonClose);
		
		
		// Settings Panel
		settingsPanel = new JPanel();
		settingsPanel.setLayout(new GridBagLayout());
		settingsPanel.setBorder(new EmptyBorder(2, 2, 2, 2));
		GridBagConstraints c1 = new GridBagConstraints();
		
		labelDegBG = new JLabel("Polynomial Degree");
		labelPolyConcavity = new JLabel("Max Polynomial Concavity");
		labelTolPK = new JLabel("Peak Tolerance");
		labelNormDrift = new JLabel("Norm Drift");
		labelSDDrift = new JLabel("SD Drift");
		
		textDegBG = new JSpinner(new SpinnerNumberModel(degBG, -1, 15, 1));
		textDegBG.setBorder(BorderFactory.createCompoundBorder(
			textDegBG.getBorder(), BorderFactory.createEmptyBorder(0, 2, 0, 2)));
		((JSpinner.DefaultEditor) textDegBG.getEditor()).getTextField().setColumns(textWidth);
		
		textPolyConcavity = new JSpinner(new SpinnerNumberModel(polyConcavity, 0.1, 20.0, 0.1));
		textPolyConcavity.setBorder(BorderFactory.createCompoundBorder(
			textPolyConcavity.getBorder(), BorderFactory.createEmptyBorder(0, 2, 0, 2)));
		((JSpinner.DefaultEditor) textPolyConcavity.getEditor()).getTextField().setColumns(textWidth);
		
		textTolPK = new JSpinner(new SpinnerNumberModel(tolPK, 0.01, 1.00, 0.005));
		textTolPK.setBorder(BorderFactory.createCompoundBorder(
			textTolPK.getBorder(), BorderFactory.createEmptyBorder(0, 2, 0, 2)));
		((JSpinner.DefaultEditor) textTolPK.getEditor()).getTextField().setColumns(textWidth);
		
		textNormDrift = new JSpinner(new SpinnerNumberModel(normDrift, 0.0, 20.0, 0.01));
		textNormDrift.setBorder(BorderFactory.createCompoundBorder(
			textNormDrift.getBorder(), BorderFactory.createEmptyBorder(0, 2, 0, 2)));
		((JSpinner.DefaultEditor) textNormDrift.getEditor()).getTextField().setColumns(textWidth);
		
		textSDDrift = new JSpinner(new SpinnerNumberModel(sdDrift, 0.0, 20.0, 0.01));
		textSDDrift.setBorder(BorderFactory.createCompoundBorder(
			textSDDrift.getBorder(), BorderFactory.createEmptyBorder(0, 2, 0, 2)));
		((JSpinner.DefaultEditor) textSDDrift.getEditor()).getTextField().setColumns(textWidth);
		
		chkBoxBands = new JCheckBox("Show Bands");
		
		buttonPanelFitType = new JPanel();
		buttonPanelFitType.setLayout(new BoxLayout(buttonPanelFitType, BoxLayout.Y_AXIS));
		buttonBands = new JRadioButton("Banded");
		buttonBands.setActionCommand("Banded");
		buttonContinuum = new JRadioButton("Continuum");
		buttonContinuum.setActionCommand("Continuum");
		
		buttonGroupFitType = new ButtonGroup();
		buttonGroupFitType.add(buttonBands);
		buttonGroupFitType.add(buttonContinuum);
		
		buttonPanelFitType.add(buttonBands);
		buttonPanelFitType.add(buttonContinuum);
		buttonPanelFitType.setBorder(new TitledBorder(BorderFactory.createEtchedBorder(),
		        "Lane Contents Pattern", TitledBorder.LEADING, TitledBorder.BELOW_TOP,
		        new Font("Sans", Font.PLAIN, 11)));
		
		final String[] lanes = new String[1 + nLanes];
		lanes[0] = "Select Ladder Lane";
		for (int i = 1; i < 1 + nLanes; i++)
			lanes[i] = "Lane " + (i);
		cmbBoxLadderLane = new JComboBox<>(lanes);
		cmbBoxDist = new JComboBox<>(distStr);
		cmbBoxLadderType = new JComboBox<>(ladderStr);
		
		buttonEditPeaks = new JToggleButton("Edit Custom Peaks");
		buttonResetCustomPeaks = new JButton("Reset Custom Peaks");
		
		c1.fill = GridBagConstraints.HORIZONTAL;
		c1.ipadx = 2; c1.ipady = 2;
		c1.gridx = 0; c1.gridy = 0; c1.weightx = 0.1;		
		settingsPanel.add(labelDegBG, c1);
		c1.gridx = 0; c1.gridy = 1; c1.weightx = 0.1;
		settingsPanel.add(labelPolyConcavity, c1);
		c1.gridx = 0; c1.gridy = 2; c1.weightx = 0.1;
		settingsPanel.add(labelTolPK, c1);
		c1.gridx = 0; c1.gridy = 3; c1.weightx = 0.1;
		settingsPanel.add(labelNormDrift, c1);
		c1.gridx = 0; c1.gridy = 4; c1.weightx = 0.1;
		settingsPanel.add(labelSDDrift, c1);
		
		c1.fill = GridBagConstraints.NONE;
		c1.gridx = 1; c1.gridy = 0; c1.weightx = 0.0;
		settingsPanel.add(textDegBG, c1);
		c1.gridx = 1; c1.gridy = 1; c1.weightx = 0.0;
		settingsPanel.add(textPolyConcavity, c1);
		c1.gridx = 1; c1.gridy = 2; c1.weightx = 0.0;
		settingsPanel.add(textTolPK, c1);
		c1.gridx = 1; c1.gridy = 3; c1.weightx = 0.0;
		settingsPanel.add(textNormDrift, c1);
		c1.gridx = 1; c1.gridy = 4; c1.weightx = 0.0;
		settingsPanel.add(textSDDrift, c1);
			
		c1.gridx = 2; c1.gridy = 0; 
		c1.gridheight = 4;
		c1.fill = GridBagConstraints.VERTICAL;
		c1.weightx = 0.0;
		settingsPanel.add(buttonPanelFitType, c1);
		
		c1.gridx = 0; c1.gridy = 5;
		c1.gridwidth = 1; c1.gridheight = 1; 
		c1.fill = GridBagConstraints.HORIZONTAL;
		settingsPanel.add(chkBoxBands, c1); 
		
		c1.gridx = 0; c1.gridy = 6; c1.gridwidth = 3; c1.gridheight = 1;
		settingsPanel.add(cmbBoxLadderLane, c1);
		
		c1.gridx = 0; c1.gridy = 7;
		settingsPanel.add(cmbBoxLadderType, c1);
		
		c1.gridx = 0; c1.gridy = 8; 
		settingsPanel.add(cmbBoxDist, c1);
		
		c1.gridx = 0; c1.gridy = 9; 
		c1.gridwidth = 1;
		settingsPanel.add(buttonEditPeaks, c1);
		
		c1.gridx = 0; c1.gridy = 10; 
		settingsPanel.add(buttonResetCustomPeaks, c1);
		
		c1.gridx = 0; c1.gridy = 11;
		c1.gridwidth = GridBagConstraints.REMAINDER; 
		c1.gridheight = 1; c1.weighty = 1.0;
		settingsPanel.add(new JLabel(), c1);
		
		dialogPanel = new JPanel();
		dialogPanel.setBackground(Color.darkGray);
		dialogPanel.setLayout(new BorderLayout());
		dialogPanel.add(buttonPanelAutoManual, BorderLayout.NORTH);
		dialogPanel.add(sliderPanel, BorderLayout.CENTER);
		dialogPanel.add(settingsPanel, BorderLayout.EAST);
		dialogPanel.add(buttonPanel, BorderLayout.SOUTH);
		
		// Initial status of all dialog components
		buttonBands.setSelected(true);
		chkBoxBands.setSelected(false);
		chkBoxBands.setEnabled(false);
		cmbBoxLadderLane.setSelectedIndex(0);
		cmbBoxLadderType.setSelectedIndex(0);
		cmbBoxLadderType.setEnabled(false);
		cmbBoxDist.setSelectedIndex(0);		
		cmbBoxDist.setEnabled(false);
		
		buttonEditPeaks.setEnabled(false);
		buttonResetCustomPeaks.setEnabled(false);
		
		// Add this class to all components as listener here for easy reference
		imp.getCanvas().addMouseMotionListener(this);
		imp.getCanvas().addMouseListener(this);
		imp.getCanvas().addMouseWheelListener(this);
		imp.getWindow().addMouseListener(this);
		imp.getWindow().addMouseMotionListener(this);
		imp.getWindow().addWindowListener(this);
		
		textNLanes.addChangeListener(this);
		textDegBG.addChangeListener(this);
		textPolyConcavity.addChangeListener(this);
		textTolPK.addChangeListener(this);
		textSDDrift.addChangeListener(this);
		textNormDrift.addChangeListener(this);

		chkBoxBands.addActionListener(this);
		buttonBands.addActionListener(this);
		buttonContinuum.addActionListener(this);

		cmbBoxDist.addActionListener(this);
		cmbBoxLadderLane.addActionListener(this);
		cmbBoxLadderType.addActionListener(this);
		
		buttonEditPeaks.addActionListener(this);
		buttonResetCustomPeaks.addActionListener(this);
		buttonFit.addActionListener(this);
		buttonClose.addActionListener(this);
		
		if (auto) {
			buttonAuto.doClick();
		} else {
			if (!readRois() || rois.size() == 0) { // go back to AUTO
				auto = true;
				prefs.putBoolean(AUTO, auto);
				buttonAuto.doClick();
			} else { // Go to MANUAL
				buttonManual.doClick();
			}
		}
		buttonAuto.setSelected(auto);
		buttonManual.setSelected(!auto);
		
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
		
		while (rois.size() == 0) {
			log.info("Loading...");
			IJ.wait(250); // delay to make sure ROIs have updated
		}
		
		reDrawROIs(imp, "none");
		imp.killRoi();
	}

	private boolean askUser(final String question) {
		final GenericDialog gd = new GenericDialog("WARNING!");
		gd.addMessage(question);
		gd.showDialog();
		if (gd.wasOKed())
			return true;
		return false;
	}

	private Peak askPeak(final int lane, final double y, final double a) {
		String title, message;
		final String fwhmString = "FWHM = \u03C3 * (2 * \u221A (2 * ln(2))";
		final double fwhmValue = 5.0;

		title = "ADD PEAK";
		message = "If the new peak is located less than " + (int) Fitter.peakDistanceTol
		        + " px away from an exisitng peak,\nthe existing peak will be replaced with the new one.";

		final GenericDialog gd = new GenericDialog(title);

		gd.addMessage("You are about to add this peak\n");
		gd.addMessage(message);
		gd.addMessage(String.format("Lane: \t%2d", lane));
		gd.addMessage(String.format("Distance: \t%10.1f", y));
		gd.addMessage(String.format("Intensity: \t%.0f", a));
		gd.addMessage("Estimate the full width at half maximum,\n" + fwhmString);
		gd.addNumericField("FWHM", fwhmValue, 1);

		gd.showDialog();
		if (gd.wasOKed()) {
			final double sd = gd.getNextNumber() / Fitter.sd2FWHM;
			return new Peak(lane, a, y, sd);
		}
		return null;
	}

	@SuppressWarnings("unchecked")
	private int[] askResetCustomPoints(final boolean[] defaultValues) {
		final String title = "RESET CUSTOM PEAKS";
		final String message = "Select the plots from which the custom peaks must be removed";
		final GenericDialog gd = new GenericDialog(title);
		gd.addMessage(message);
		final int[] lanes = getAllLaneNumbers();
		final int rows = lanes.length + 1;
		final List<Integer> selectedLanes = new ArrayList<>();
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
			final Vector<Checkbox> v = gd.getCheckboxes();
			for (final Checkbox cb : v) {
				if (cb.getState()) {
					selectedLanes.add(v.indexOf(cb) + 1);
				}
			}
			if (selectedLanes.contains(rows)) {
				// 'All Lanes' is selected
				if (selectedLanes.size() != 1) { // Not just 'All lanes'
					final boolean[] b = new boolean[rows];
					b[b.length - 1] = true;
					askResetCustomPoints(b);
				}
				return lanes;
			}
			final int[] intArray = new int[selectedLanes.size()];
			for (int i = 0; i < intArray.length; i++) {
				intArray[i] = selectedLanes.get(i);
			}
			return intArray;
		}
		if (gd.wasCanceled()) {
			return new int[0];
		}
		return new int[0];
	}

	private String[] askLadderRange(final String[] ladder) {
		final String title = "LADDER RANGE";
		final String message = "Select the first and last ladder bands to consider";
		final GenericDialog gd = new GenericDialog(title);
		String ladderFirst = ladder[ladderRange[0]];
		String ladderLast = ladder[ladderRange[1]];
		if (ladderRange[1] == 0) ladderLast = ladder[ladder.length - 1];
		gd.addMessage(message);
		gd.addChoice("First Band", ladder, ladderFirst);
		gd.addChoice("Last Band", ladder, ladderLast);
		gd.showDialog();
		if (gd.wasOKed()) {
			final String[] range = { gd.getNextChoice(), gd.getNextChoice() };
			return range;
		}
		return null;
	}

	private JSlider makeTitledSlider(final String string, final Color color, final int minVal,
	        final int maxVal, int val) {
		if (val < minVal) val = minVal;
		if (val > maxVal) val = maxVal;
		final JSlider slider = new JSlider(SwingConstants.HORIZONTAL, minVal, maxVal, val);
		final TitledBorder tb = new TitledBorder(BorderFactory.createEtchedBorder(),
		        // empty,
		        "", TitledBorder.CENTER, TitledBorder.BELOW_TOP, new Font("Sans", Font.PLAIN, 11));
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
		} else {
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
			// Housekeeping: Sort the rois based on POSITION and remove null
			// elements
			final Comparator<Roi> roiNameComparator = new Comparator<Roi>() {

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
					final double xmin1 = r1.getXBase();
					final double ymin1 = r1.getYBase();
					final double xmin2 = r2.getXBase();
					final double ymin2 = r2.getYBase();
					final int compx = Double.compare(xmin1, xmin2);
					if (compx != 0) {
						return compx;
					}
					return Double.compare(ymin1, ymin2);
				}
			};

			Collections.sort(rois, roiNameComparator);
			// Rename Rois in order and remove null
			final Iterator<Roi> roisIter = rois.iterator();
			int nullIndex = 0;
			while (roisIter.hasNext()) {
				final Roi r = roisIter.next();
				if (r == null)
					break;
				r.setName("Lane " + (nullIndex + 1));
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
			final double lh1 = font.getSize() * 1.4;
			final double lw1 = font.getSize() * 1.4;
			final double rh = roi.getBounds().getHeight();
			final double rw = roi.getBounds().getWidth();
			final double x0 = roi.getBounds().getX();
			final double y0 = roi.getBounds().getY();
			if (lh1 > rh || lw1 > rw)
				labelRoi.setLocation(x0 + 0.1 * lw1, y0 + rh + 1.1 * lh1);
			else
				labelRoi.setLocation(x0 + 0.2 * lw1, y0 + rh - 1.1 * lh1);
			if (roi.getName().equals(roiName)) {
				roi.setStrokeColor(MainDialog.selected);
				roi.setStrokeWidth(3);
				if (buttonManual.isSelected() && !roiName.equals("none")) {
					imgPlus.setRoi(roi);
				}
			} else if (roi.getName().equals(roiReferenceName)) {
				if (!roiReferenceName.equals("none") && !roiReferenceName.equals(roiName)) {
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
			if (buttonAuto.isSelected())
				imgPlus.killRoi();
			if (chkBoxBands.isSelected()) {
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
		imgPlus.setOverlay(overlay);
	}

	private void resetAutoROIs() {
		rois = new ArrayList<>();
		for (int i = 0; i < nLanes; i++) {
			final Roi roi = new Roi(lhoff + lw * i + lsp * i, lvoff, lw, lh);
			roi.setName("Lane " + (i + 1));
			rois.add(roi);
		}
	}

	private void redoProfilePlots() {
		fitDone = false;
		chkBoxBands.setSelected(false);
		chkBoxBands.setEnabled(false);
		buttonEditPeaks.setEnabled(false);
		buttonResetCustomPeaks.setEnabled(false);
		plotter.resetData();
		plotter.setPlotMode(Plotter.regMode);
		//rois.parallelStream().forEach((r) -> {		
		for (final Roi r : rois) {
			final int ln = Integer.parseInt(r.getName().substring(5));
			final Rectangle rect = r.getBounds();
			// Do not ROIs that are outside the image
			if (rect.getMinX() < 0.95 * iw && rect.getMinY() < 0.95 * ih) {
				plotter.updateProfile(r);
				final RealVector empty = new ArrayRealVector();
				final DataSeries d = new DataSeries("Custom Points", ln, DataSeries.CUSTOMPEAKS,
				        empty, empty, Plotter.vMarkerEditPeakColor);
				d.addChangeListener(this);
				plotter.addDataSeries(d);
			}
		//});
		}
		plotter.reloadTabs();
	}

	private double[][] readDistFile(final String filename) {
		String line = null;
		final List<Integer> fragmentLength = new ArrayList<>();
		final List<Integer> fragmentFrequency = new ArrayList<>();

		final URL url = MainDialog.class.getClassLoader().getResource(filename);
		log.info("Loading " + url.getPath() + " ...");
		try (BufferedReader buffer = new BufferedReader(new InputStreamReader(url.openStream()))) {
			while (true) {
				line = buffer.readLine();
				if (line == null)
					break;

				final String[] words = line.split("\t");
				try {
					fragmentLength.add(Integer.parseInt(words[0].trim()));
					fragmentFrequency.add(Integer.parseInt(words[1].trim()));
				} catch (final NumberFormatException e1) {
					log.info("Invalid token: " + words[0].trim() + " " + words[1].trim());
				}
			}
			buffer.close();
		} catch (final Exception e) {
			e.printStackTrace();
		}
		if (fragmentLength.size() != fragmentFrequency.size())
			return null;
		final double[][] out = new double[fragmentLength.size()][3];
		int count = 0;
		for (int i = 0; i < fragmentLength.size(); i++) {
			out[i][0] = fragmentFrequency.get(i);
			out[i][1] = fragmentLength.get(i);
			out[i][2] = fragmentLength.get(i) * 607.4 + 157.9; //MW
			count = count + fragmentFrequency.get(i);
		}
		for (int i = 0; i < fragmentLength.size(); i++) {
			out[i][0] = out[i][0] / count;
		}
		// 3 columns: Relative Frequency, Length, MW
		return out;
	}

	private boolean readRois() {
		rois = new ArrayList<>();
		final String path = "gel-lanes-fit/data/savedrois.bak";
		log.info("Loading " + path + " ...");
		try (ObjectInputStream ois = new ObjectInputStream(new FileInputStream(path))) {
			try {
				int i = 1;
				while (true) {
					final Rectangle r = (Rectangle) ois.readObject();
					final Roi roi = new Roi(r);
					roi.setName("Lane " + i++);
					rois.add(roi);
				}
			} catch (final Exception e) {
				/* Exit */ }
			ois.close();
			return true;
		} catch (final IOException e) {
			log.info("Could not find previously saved ROIs");
			return false;
		}
	}

	private void resetCustomPeaks(final int lane) {
		fitter.resetCustomPeaks(lane);
		DataSeries d = plotter.getPlotsCustomPeaks(lane);
		if (d == null) {
			final RealVector empty = new ArrayRealVector();
			d = new DataSeries("Custom Points", lane, DataSeries.CUSTOMPEAKS, empty, empty, Plotter.vMarkerEditPeakColor);
			plotter.addDataSeries(d);
		} else {
			d.removeChangeListener(this);
			d.clear();
		}
		d.addChangeListener(this);
	}

	private boolean saveRois() {
		final String path = "gel-lanes-fit/data/";
		final String file = "savedrois.bak";
		new File(path).mkdirs();
		log.info("Saving to " + path + file + " ...");
		try (ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(path + file))) {
			for (final Roi r : rois) {
				oos.writeObject(r.getBounds());
			}
			return true;
		} catch (final IOException e) {
			log.error("ROI file not Found. Creating a new one.");
			e.printStackTrace();
			return false;
		}
	}

	private void updateLaneComboBox() {
		final int idx = cmbBoxLadderLane.getSelectedIndex();
		boolean found = false;
		int ref = 0;
		cmbBoxLadderLane.removeActionListener(this);
		if (!roiReferenceName.equals("none")) {
			ref = Integer.parseInt(roiReferenceName.substring(5));
		}
		final String str = cmbBoxLadderLane.getItemAt(0);
		cmbBoxLadderLane.removeAllItems();
		cmbBoxLadderLane.addItem(str);
		for (final Roi r : rois) {
			final String name = r.getName();
			final int ln = Integer.parseInt(name.substring(5));
			cmbBoxLadderLane.addItem(name);
			if (ln == ref) {
				found = true;
			}
		}
		cmbBoxLadderLane.addActionListener(this);
		if (found)
			cmbBoxLadderLane.setSelectedIndex(idx);
		else
			cmbBoxLadderLane.setSelectedIndex(0);
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

		saveRois();
		plotter.closePlot();
		for (final Display<?> d : displayServ.getDisplays()) {
			log.info(d.getName() + " is closing...");
			d.close();
		}
		log.info("\nGel Lanes Fit terminated.");
	}

	
	/**
	 * @return ROI identified by title. Typically "Lane n", where int n > 0
	 */
	private Roi getRoi(final String title) {
		final Iterator<Roi> roiIter = rois.iterator();
		while (roiIter.hasNext()) {
			final Roi roi = roiIter.next();
			if (roi.getName().equals(title))
				return roi;
		}
		return null;
	}

	/**
	 * @return @param rois, the current ROI set
	 */
	public List<Roi> getRois() {
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
		fitter.setDegBG(degBG);
		fitter.setPolyConcavity(polyConcavity);
		fitter.setTolPK(tolPK);
		fitter.setNormDrift(normDrift);
		fitter.setSDDrift(sdDrift);
	}

	public void setPlotter(final Plotter plotter) {
		plotter.initialize(rois);
		// TODO: Think of a better place for this:
		// Add empty Custom Peak DataSeries with listeners
		for (final int i : getAllLaneNumbers()) {
			final RealVector empty = new ArrayRealVector();
			final DataSeries d = new DataSeries("Custom Points", i, DataSeries.CUSTOMPEAKS, empty,
			        empty, Plotter.vMarkerEditPeakColor);
			d.addChangeListener(this);
			plotter.addDataSeries(d);
		}
		this.plotter = plotter;
	}

	@Override
	public synchronized void stateChanged(final ChangeEvent e) {
		final Object o = e.getSource();
		final List<Object> fitDisruptors = new ArrayList<>(
				Arrays.asList(sliderW, sliderH, sliderSp, 
										sliderHOff, sliderVOff, textNLanes));
		
		if (fitDisruptors.contains(o)) {
			if (fitDone && !askUser(warningFit))
				return;
		}
		// ROI SLIDERS
		if (o == sliderW) {
			lw = sliderW.getValue(); prefs.putInt(LW, lw);
			String str = "Width ( " + lw + " px )";
			setSliderTitle(sliderW, str);
			sliderUpdate();
		} else if (o == sliderH) {
			lh = sliderH.getValue(); prefs.putInt(LH, lh);
			String str = "Height ( " + lh + " px )";
			setSliderTitle(sliderH, str);
			sliderUpdate();
		} else if (o == sliderSp) {
			lsp = sliderSp.getValue(); prefs.putInt(LSP, lsp);
			String str = "Spacing ( " + lsp + " px )";
			setSliderTitle(sliderSp, str);
			sliderUpdate();
		} else if (o == sliderHOff) {
			lhoff = sliderHOff.getValue(); prefs.putInt(LHOFF, lhoff);
			String str = "Horizontal Offset ( " + lhoff + " px )";
			setSliderTitle(sliderHOff, str);
			sliderUpdate();
		} else if (o == sliderVOff) {
			lvoff = sliderVOff.getValue(); prefs.putInt(LVOFF, lvoff);
			String str = "Vertical Offset ( " + lvoff + " px )";
			setSliderTitle(sliderVOff, str);
			sliderUpdate();
		} 
		// SPINNERS
		if (o == textNLanes) {
			if (fitDone && !askUser(warningFit)) {
				textNLanes.setValue(nLanes);
				return;
			}
			nLanes = (int) textNLanes.getModel().getValue();
			prefs.putInt(NLANES, nLanes);
			resetAutoROIs();
			reDrawROIs(imp, "none");
			updateLaneComboBox();
			redoProfilePlots();
			fitter.resetAllFitter();
			fitter.setInputData(plotter.getProfiles());
		} else if (o == textDegBG) {
			degBG = (int) textDegBG.getModel().getValue();
			fitter.setDegBG(degBG);
			prefs.putInt(DEGBG, degBG);
		} else if (o == textPolyConcavity) {
			polyConcavity = (double) textPolyConcavity.getModel().getValue();
			fitter.setPolyConcavity(polyConcavity);
			prefs.putDouble(POLYCONCAVITY, polyConcavity);
		} else if (o == textTolPK) {
			tolPK = (double) textTolPK.getModel().getValue();
			fitter.setTolPK(tolPK);
			prefs.putDouble(TOLPK, tolPK);
		} else if (o == textNormDrift) {
			normDrift = (double) textNormDrift.getModel().getValue();
			fitter.setNormDrift(normDrift);
			prefs.putDouble(NORMDRIFT, normDrift);
		} else if (o == textSDDrift) {
			sdDrift = (double) textSDDrift.getModel().getValue();
			fitter.setSDDrift(sdDrift);
			prefs.putDouble(SDDRIFT, sdDrift);
		}
	}
	
	void sliderUpdate() {
		resetAutoROIs();
		reDrawROIs(imp, "none");
		redoProfilePlots();
		plotter.setReferencePlot(MainDialog.noLadderLane);
		fitter.resetAllFitter();
		fitter.setInputData(plotter.getProfiles());
	}

	@Override
	public void seriesChanged(final SeriesChangeEvent e) {
		final DataSeries d = (DataSeries) e.getSource();
		if (d.getType() == DataSeries.CUSTOMPEAKS) {
			final int lane = d.getLane();
			final List<Peak> peaks = fitter.getCustomPeaks(lane);
			final List<Double> m1 = new ArrayList<>();
			final List<Double> m2 = new ArrayList<>();
			final List<Double> n1 = new ArrayList<>();
			final List<Double> n2 = new ArrayList<>();

			for (int i = 0; i < d.getItemCount(); i++) {
				m1.add(i, d.getX(i).doubleValue());
				n1.add(i, d.getY(i).doubleValue());
			}

			for (final Peak p : peaks) {
				m2.add(peaks.indexOf(p), p.getMean());
				n2.add(peaks.indexOf(p), p.getNorm());
			}
			if (m1.size() > m2.size()) { // ADD
				for (final double m : m1) {
					if (!m2.contains(m)) {
						final Peak peak = askPeak(d.getLane(), m, n1.get(m1.indexOf(m)));
						if (peak == null) {
							d.removeChangeListener(this);
							d.remove(m);
							d.addChangeListener(this);
							log.info("Custom Peak not added: " +d.getLane() + ", " + n1.get(m1.indexOf(m)) + ", " + m);
						} else {
							fitter.addCustomPeak(peak);
							log.info("Custom Peak added: " + d.getLane() + ", " + n1.get(m1.indexOf(m)) + ", " + m);
						}
					}
				}
			} else if (m1.size() < m2.size()) { // REMOVE
				for (final double m : m2) {
					if (!m1.contains(m)) {
						if (fitter.removeCustomPeak(new Peak(lane, n2.get(m2.indexOf(m)), m))) {
							log.info("Custom Peak Removed: " + n2.get(m2.indexOf(m)) + ", " + m);
						}
					}
				}
			}
			final int peakNumber = fitter.getCustomPeaks(lane).size();
			if (d.getItemCount() != peakNumber) {
				log.error("Wrong custom peaks count: " + d.getItemCount() + ", " + peakNumber);
			}
		}
	}

	@Override
	public void actionPerformed(final ActionEvent e) {
		// Buttons
		// ----------------------------------------------------------------
		if (e.getSource().equals(buttonFit)) {
			if (fitDone && !askUser(warningFit))
				return;
			
			if (roiReferenceName.equals("none")) {
				askUser("Reference ladder not selected");
				return;
			}
			
			if (cmbBoxLadderType.getSelectedItem().equals("Select Ladder Type")) {
				askUser("Ladder type not selected");
				return;
			}
			
			if (buttonBands.isSelected())
				fitter.setFitMode(Fitter.bandMode);
			else if (buttonContinuum.isSelected()) {
				fitter.setFitMode(Fitter.fragmentMode);
				if (cmbBoxDist.getSelectedItem().equals("Select Fragment Distribution")) {
					askUser("Fragment distribution not selected");
					buttonBands.setSelected(true);
					fitter.setFitMode(Fitter.bandMode);
					return;
				}
			}
			
			buttonEditPeaks.setSelected(false);
			buttonResetCustomPeaks.setSelected(false);
			
			// Remove everything in the plots, but the profile
			plotter.setSelected(MainDialog.noLaneSelected);
			plotter.setPlotMode(Plotter.regMode);
			plotter.removeFit();
			plotter.removeVerticalMarkers();
			
			// Reset fit
			for (final int l : getAllLaneNumbers())
				fitter.resetFit(l);
			fitter.setInputData(plotter.getProfiles());
			
			chkBoxBands.setEnabled(true);
			chkBoxBands.setSelected(true);
			buttonEditPeaks.setEnabled(true);
			buttonResetCustomPeaks.setEnabled(true);
			
			List<DataSeries> fitted = new ArrayList<>();
			fitted = fitter.doFit(roiReferenceNumber);
			fitter.updateResultsTable();
			plotter.addDataSeries(fitted);
			plotter.addVerticalMarkers(fitter.getFittedPeaks(roiReferenceNumber));
			for (final int i : getAllLaneNumbers())
				plotter.updatePlot(i);
			final List<Peak> peaks = fitter.getFittedPeaks(roiReferenceNumber);
			final int l1 = peaks.size();
			final int l2 = ladderRange[1] - ladderRange[0] + 1;
			if (l1 != l2) {
				String message = "The number of peaks detected does not match the ladder range";
				message = message + "\n Peaks detected: " + l1;
				message = message + "\n Bands in range: " + l2;
				askUser(message);
				return;
			}
			for (final Peak p : peaks) { // Rename ladder peaks with band names
				p.setName(hilo[ladderRange[0] + peaks.indexOf(p)]);
			}
			plotter.addVerticalMarkers(peaks);
			
			fitted.addAll(fitter.doFit());
			fitter.updateResultsTable();
			fitDone = true;
			plotter.removeFit();
			plotter.addDataSeries(fitted);
			for (final int i : getAllLaneNumbers()) {
				if (i != roiReferenceNumber) plotter.updatePlot(i);
			}
			reDrawROIs(imp, "none"); // adds the bands to the ROIs
		}
		
		if (e.getSource().equals(buttonAuto)) {
			if (!auto) { // NOT already AUTO
				if (!saveRois()) {
					String warn = "When switching to AUTO mode, the current lane selections will be reset!";
					if (!askUser(warn)) {
						buttonManual.setSelected(true);
						return;
					}
				}
				auto = true;
				prefs.putBoolean(AUTO, auto);
			}
			setSliderPanelEnabled(true);
			resetAutoROIs();
			reDrawROIs(imp, "none");
			if (plotter == null || fitter == null) return;
			redoProfilePlots();
			fitter.resetAllFitter();
			fitter.setInputData(plotter.getProfiles());
		}
		
		if (e.getSource().equals(buttonManual)) {
			auto = false;
			prefs.putBoolean(AUTO, auto);
			setSliderPanelEnabled(false);
			if (!readRois() || rois.size() == 0) { // Use the AUTO rois as a start
				resetAutoROIs();
				reDrawROIs(imp, "none");
			}
			if (plotter == null || fitter == null) return;
			reDrawROIs(imp, "none");
			redoProfilePlots();
			fitter.resetAllFitter();
			fitter.setInputData(plotter.getProfiles());
		}
		
		if (e.getSource().equals(buttonEditPeaks)) {
			if (buttonEditPeaks.isSelected())
				plotter.setPlotMode(Plotter.editPeaksMode);
			else
				plotter.setPlotMode(Plotter.regMode);
		}
		
		if (e.getSource().equals(buttonResetCustomPeaks)) {
			final int[] lanes = askResetCustomPoints(new boolean[rois.size() + 1]);
			for (final int i : lanes) {
				resetCustomPeaks(i);
			}
			reDrawROIs(imp, "none");
		}
		
		if (e.getSource().equals(buttonClose)) {
			if (askUser("Would you like to quit Gel Lanes Fit?")) {
				frame.dispatchEvent(new WindowEvent(frame, WindowEvent.WINDOW_CLOSING));
			}
		}
		
		// CheckBoxes
		// -------------------------------------------------------------
		if (e.getSource().equals(chkBoxBands)) {
			reDrawROIs(imp, "none");
		}
		
		// ComboBoxes
		// -------------------------------------------------------------
		if (e.getSource().equals(cmbBoxLadderLane)) {
			final String ladderLane = (String) cmbBoxLadderLane.getSelectedItem();
			if (ladderLane.equals("Select Ladder Lane")) {
				roiReferenceName   = "none";
				roiReferenceNumber = MainDialog.noLadderLane;
				cmbBoxLadderType.setEnabled(false);
				cmbBoxDist.setEnabled(false);
			} else {
				roiReferenceName = ladderLane;
				roiReferenceNumber = Integer.parseInt(ladderLane.substring(5));
				cmbBoxLadderType.setEnabled(true);
				cmbBoxDist.setEnabled(true);
			}
			fitter.setReferenceLane(roiReferenceNumber);
			plotter.setReferencePlot(roiReferenceNumber);
			reDrawROIs(imp, "none");
		}
		
		if (e.getSource().equals(cmbBoxLadderType)) {
			String[] ladderBandsStrings = null;
			String ladderType = ladderStr[cmbBoxLadderType.getSelectedIndex()];
			
			if (ladderType.equals("Hi-Lo")) {
				ladderBandsStrings = hilo;
			} else if (ladderType.equals("100bp")) {
				ladderBandsStrings = bp100;
			}
			
			if (!ladderType.equals(ladderStr[0]) && ladderBandsStrings != null) {
				// Ask user on limits of ladder
				final String[] bands = askLadderRange(ladderBandsStrings);
				if (bands != null) {
					ladderRange[0] = Arrays.asList(ladderBandsStrings).indexOf(bands[0]);
					ladderRange[1] = Arrays.asList(ladderBandsStrings).indexOf(bands[1]);
					
					RealVector ladder = new ArrayRealVector();
					for (int i = ladderRange[0]; i <= ladderRange[1]; i++) {
						int unit = 0;
						final String band = ladderBandsStrings[i];
						if (band.contains("k"))
							unit = 1000;
						else
							unit = 1;
						final String[] bandWords = StringUtils.split(band);
						ladder = ladder.append(Double.parseDouble(bandWords[0]) * unit);
					}
					ladder.mapMultiplyToSelf(607.4).mapAddToSelf(157.9);
					fitter.setLadder(ladder);
				}
			}
			cmbBoxLadderType.removeAllItems();
			for (final String s : ladderStr) {
				if (s.equals(ladderType) && ladderBandsStrings != null) {
					ladderType = s + " [" + ladderBandsStrings[ladderRange[0]] + " - "
					        + ladderBandsStrings[ladderRange[1]] + "]";
					cmbBoxLadderType.addItem(ladderType);
				} else
					cmbBoxLadderType.addItem(s);
			}
			cmbBoxLadderType.setSelectedItem(ladderType);
		}
		
		if (e.getSource().equals(cmbBoxDist)) {
			if (cmbBoxDist.getSelectedIndex() == 0) return;
			final String filename = "data/" + cmbBoxDist.getSelectedItem() + ".txt";
			final double[][] dist = readDistFile(filename);
			fitter.setFragmentDistribution(dist);
		}
	}

	@Override
	public void mouseDragged(final MouseEvent e) {
		// nothing to do
	}

	@Override
	public void mouseMoved(final MouseEvent e) {
		if (e.getSource() == imp.getCanvas()) {
			if (rois.size() == 0 || plotter == null)
				return;
			final int x = ((ImageCanvas) e.getSource()).offScreenX(e.getX());
			final int y = ((ImageCanvas) e.getSource()).offScreenX(e.getY());
			statusServ.showStatus("[" + x + ":" + y + "]");
			if (!selectionUpdate) { // Not dragging an ROI
				String roiCurrent = "none"; // None selected
				for (final Roi r : rois) {
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
						final int ln = Integer.parseInt(roiCurrent.substring(5));
						plotter.setSelected(ln);
						plotter.setVLine(ln, y);
					} else {
						plotter.setSelected(MainDialog.noLadderLane);
					}

					if (!roiPreviouslySelected.equals("none")) {
						final int roiPS = Integer.parseInt(roiPreviouslySelected.substring(5));
						plotter.removeVLine(roiPS);
					}
				} else { // Moving inside same ROI
					if (!roiCurrent.equals("none")) {
						final int ln = Integer.parseInt(roiCurrent.substring(5));
						plotter.setVLine(ln, y);
					}
				}
			}
		}
	}

	@Override
	public void mouseClicked(final MouseEvent e) {
		if (e.getSource().equals(imp.getCanvas())) {
			if (!auto && !roiSelected.equals("none") && !selectionUpdate) {
				// Remove Roi
				final Roi roi = getRoi(roiSelected);
				if (roi == null || roi.isHandle(e.getX(), e.getY()) != -1) {
					// Might be dragging existing roi
					return;
				}
				
				if (askUser("Delete " + roiSelected + "?")) {
					if (fitDone && !askUser(warningFit))
						return;
					final Iterator<Roi> roiIter = rois.iterator();
					while (roiIter.hasNext()) {
						final Roi roiRemove = roiIter.next();
						if (roiRemove.getName().equals(roiSelected)) {
							roiIter.remove();
						}
					}
					saveRois();
					reDrawROIs(imp, "none");
					updateLaneComboBox();
					redoProfilePlots();
				}
			}
		}
	}

	@Override
	public void mousePressed(final MouseEvent e) {
		if (!auto)
			selectionUpdate = true;
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
						if (roiNew.getName() != null) { // Existing ROI
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
						reDrawROIs(imp, roiSelected);
						redoProfilePlots();
						updateLaneComboBox();
					}
				}
				selectionUpdate = false;
				saveRois();
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
			if (fitDone && !askUser(warningFit))
				return;
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
