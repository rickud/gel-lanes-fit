/*
 * Gel Lanes Fit - MainDialog.java
 * Author: Rick Ziraldo, 2017
 * The University of Texas at Dallas, Richardson, TX
 *
 * Licensed under the GNU Affero General Public License v3.0; see LICENSE.
 * Source: https://github.com/rickud/gel-lanes-fit
 */

package gellanesfit;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Checkbox;
import java.awt.Color;
import java.awt.Component;
import java.awt.Desktop;
import java.awt.Dimension;
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
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.ConcurrentModificationException;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;
import java.util.prefs.Preferences;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSlider;
import javax.swing.JSpinner;
import javax.swing.JSplitPane;
import javax.swing.JTabbedPane;
import javax.swing.JToggleButton;
import javax.swing.SpinnerNumberModel;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.border.EmptyBorder;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.general.SeriesChangeEvent;
import org.jfree.data.general.SeriesChangeListener;
import org.jfree.data.xy.XYSeriesCollection;
import org.scijava.Context;
import org.scijava.app.StatusService;
import org.scijava.display.Display;
import org.scijava.display.DisplayService;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.GUI;
import ij.gui.GenericDialog;
import ij.gui.HTMLDialog;
import ij.gui.ImageCanvas;
import ij.gui.Line;
import ij.gui.MessageDialog;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.gui.TextRoi;
import ij.io.FileInfo;
import ij.io.FileSaver;

/**
 * The plugin's main window, and the controller of an analysis: the lanes drawn
 * on the image, the ladder and fit settings, the fit itself (run through
 * {@link Fitter} and shown in the {@link Plotter}), custom peaks, and the saved
 * state of the image (saved-state.bak in the image's data folder).
 */
class MainDialog extends JFrame implements ActionListener, ChangeListener,
	SeriesChangeListener, MouseMotionListener, MouseListener, MouseWheelListener,
	WindowListener
{
	private static final long serialVersionUID = 1L;

	@Parameter
	private LogService log;

	@Parameter
	private DisplayService displayServ;

	@Parameter
	private StatusService statusServ;

	// Preference keys for this class
	private static final String DEGBG = "degBG";
	private static final String POLYDERIVATIVE = "polyDerivative";
	private static final String TOLPK = "tolPK";
	private static final String AREADRIFT = "areaDrift";
	private static final String SDDRIFT = "sdSDrift";
	private static final String DLO = "dlo"; // Uniform distribution range
	private static final String DHI = "dhi";
	private static final String EVERY = "every";
	
	private static final String NLANES = "nLanes";
	private static final String LADDERLANEINT = "ladderLaneInt";
	private static final String LW = "lw";
	private static final String LH = "lh";
	private static final String LSP = "lsp";
	private static final String LHOFF = "lhoff";
	private static final String LVOFF = "lvoff";
	private static final String AUTO = "auto";
	private static final String OLDIMPTITLE = "oldImpTitle";
	private static final String SKIPFITWARNING = "skipFitWarning";

	private final double SW = IJ.getScreenSize().getWidth();
	private final double SH = IJ.getScreenSize().getHeight();

	private static final Color guess = Color.BLUE;
	private static final Color custom = Color.GREEN;
	private static final Color fit = Color.MAGENTA;
	private static final Color selected = Color.YELLOW;
	private static final Color unselected = Color.RED;
	private static final Color reference = Color.BLUE;

	public static final int noLadderLane = 0;
	public static final int noLaneSelected = -1;

	private boolean auto; // AUTO ROI mode; all lanes same size
	private boolean selectionUpdate = false; // active updating is off
	private boolean fitDone = false; // Keep track of whether fit data exists

	private final String impTitle;
	private String oldImpTitle; // Stored in Prefs
	private final String savePath;
	private String roiSelected = "none";
	private String roiPreviouslySelected = "none";
	private String ladderLaneStr = "none";
	private int ladderLaneInt = noLadderLane;
	private final String[] ladderStr = { "Select Ladder Type", "Hi-Lo", "100bp", "Quick-Load", "Tapestation","Custom Ladder" };
	private final String[] distStr = { "Select Fragment Distribution",
		"AciI-Lambda", "AciI-Lambda4", "AciI-Lambda3", "AciI-Lambda2",
		"Uniform", "Ladder" };
	private Ladder ladder;
	private final String warningFit =
		"The current plots will be reset and the current fitting data will be lost.";

	private final Plotter plotter;
	private final Fitter fitter;
	private final Preferences prefs;

	private final ImagePlus imp;
	private ArrayList<Roi> rois;
	private List<Peak> savedCustomPeaks = new ArrayList<>();
	private List<Rectangle> savedRects = new ArrayList<>(); // Manual ROIs
	// Custom peaks outside the current lanes, kept for when they fit again
	private List<Peak> parkedPeaks = new ArrayList<>();
	private FitState savedState; // Settings of the last session, if any
	private boolean restoring = true; // Do not save while restoring the state

	private int iw, ih, lw, lh, lsp, lhoff, lvoff, dlo, dhi, every;
	private int nLanes;
	private int degBG; // Order of Background Polynomial
	private double polyDerivative;
	private double tolPK; // Peak detection tolerance as % of range
	private double areaDrift;
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
	private JButton buttonOpenFolder;
	private JButton buttonClose;

	private JPanel settingsPanel;
	private JLabel labelDegBG;
	private JLabel labelPolyDerivative;
	private JLabel labelTolPK;
	private JLabel labelAreaDrift;
	private JLabel labelSDDrift;

	private JSpinner textDegBG;
	private JSpinner textPolyDerivative;
	private JSpinner textTolPK;
	private JSpinner textAreaDrift;
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

	private JPanel topButtonsPanel;
	private JPanel dialogPanel;
	private final JFrame frame;

	MainDialog(final Context context, final String string, final ImagePlus imp,
		final Preferences prefs, final Plotter plotter, final Fitter fitter)
	{
		context.inject(this);
		this.prefs = prefs;
		this.fitter = fitter;
		this.plotter = plotter;
		this.imp = imp;
		this.impTitle = imp.getTitle().substring(0, imp.getTitle().indexOf("."));
		final String sep = File.separator;
		// Absolute, next to the image, so it does not depend on how Fiji started
		final FileInfo fi = imp.getOriginalFileInfo();
		final String imageDir = fi != null && fi.directory != null && !fi.directory
			.isEmpty() ? fi.directory : IJ.getDirectory("imagej");
		final String title = imp.getTitle();
		final String imageName = title.lastIndexOf('.') > 0 ? title.substring(0,
			title.lastIndexOf('.')) : title;
		final File dataDir = new File(imageDir, imageName + " - Gel Lanes Fit");
		moveOldDataFolder(new File(new File(imageDir, "gel-lanes-fit"), imp
			.getShortTitle()), dataDir);
		savePath = dataDir.getAbsolutePath() + sep;
		if (!dataDir.isDirectory() && !dataDir.mkdirs())
			log.error("Cannot create the data folder " + savePath);
		log.info("Data folder: " + savePath);
		frame = new JFrame(string);
		rois = new ArrayList<>();
		ladder = null;
		setupMainDialog();
		frame.setLocation((int) SW / 16, (int) SH / 8);
	}

	/**
	 * Data used to go to gel-lanes-fit/(short title)/ next to the image;
	 * moves it to the new folder so the saved state is found again
	 */
	private void moveOldDataFolder(final File oldDir, final File newDir) {
		if (!oldDir.isDirectory() || newDir.exists()) return;
		if (oldDir.renameTo(newDir)) {
			log.info("Moved " + oldDir + " to " + newDir);
			final File parent = oldDir.getParentFile();
			final String[] left = parent.list();
			if (left != null && left.length == 0) parent.delete();
		}
		else log.error("Could not move " + oldDir + " to " + newDir);
	}

	private void setupMainDialog() {
		final int textWidth = 5;
		iw = imp.getWidth();
		ih = imp.getHeight();
		// Default lane size/offset
		oldImpTitle = prefs.get(OLDIMPTITLE, "None");
		degBG = prefs.getInt(DEGBG, 2);
		polyDerivative = prefs.getDouble(POLYDERIVATIVE, 10.0);
		tolPK = prefs.getDouble(TOLPK, 0.1);
		areaDrift = prefs.getDouble(AREADRIFT, 0.1);
		sdDrift = prefs.getDouble(SDDRIFT, 1.0);

		if (oldImpTitle.equals(impTitle)) {
			auto = prefs.getBoolean(AUTO, true);
			nLanes = prefs.getInt(NLANES, 4);
			ladderLaneInt = prefs.getInt(LADDERLANEINT, noLadderLane);
			
			lw = prefs.getInt(LW, (int) Math.round(0.8 * iw / nLanes));
			lh = prefs.getInt(LH, (int) Math.round(ih * 0.8));
			lsp = prefs.getInt(LSP, Math.round((iw - lw * nLanes) / (nLanes + 1)));
			lhoff = prefs.getInt(LHOFF, lsp / 2);
			lvoff = prefs.getInt(LVOFF, (ih - lh) / 2);
			dlo = prefs.getInt(DLO, 100);
			dhi = prefs.getInt(DHI, 500);
			every = prefs.getInt(EVERY, 1);
		}
		else {
			auto = true;
			nLanes = 4;
			ladderLaneInt = noLadderLane;

			lw = (int) Math.round(0.8 * iw / nLanes);
			lh = (int) Math.round(ih * 0.8);
			lsp = Math.round((iw - lw * nLanes) / (nLanes + 1));
			lhoff = lsp / 2;
			lvoff = (ih - lh) / 2;
			dlo = 100;
			dhi = 500;
			every = 1;
			
			// Create/Update prefs node keys
			prefs.putBoolean(AUTO, auto);
			prefs.putInt(NLANES, nLanes);
			prefs.putInt(LADDERLANEINT, ladderLaneInt);
			prefs.put(OLDIMPTITLE, impTitle);
			prefs.putInt(LW, lw);
			prefs.putInt(LH, lh);
			prefs.putInt(LSP, lsp);
			prefs.putInt(LHOFF, lhoff);
			prefs.putInt(LVOFF, lvoff);
			prefs.putInt(DLO, dlo);
			prefs.putInt(DHI, dhi);
			prefs.putInt(EVERY, every);
		}

		// The saved state of this image overrides the preferences
		loadState();
		if (savedState != null) {
			auto = savedState.auto;
			nLanes = savedState.nLanes;
			lw = savedState.lw;
			lh = savedState.lh;
			lsp = savedState.lsp;
			lhoff = savedState.lhoff;
			lvoff = savedState.lvoff;
			ladderLaneInt = savedState.ladderLane;
			degBG = savedState.degBG;
			polyDerivative = savedState.polyDerivative;
			tolPK = savedState.tolPK;
			areaDrift = savedState.areaDrift;
			sdDrift = savedState.sdDrift;
			dlo = savedState.dlo;
			dhi = savedState.dhi;
			every = savedState.every;
		}

		ladderLaneStr = ladderLaneInt == 0 ? "none" : "Lane " + ladderLaneInt;
		buttonPanelAutoManual = new JPanel();
		buttonAuto = new JRadioButton("Automatic Rectangle Selection");
		buttonAuto.setActionCommand("Auto");
		buttonAuto.addActionListener(this);
		buttonManual = new JRadioButton("Manual Rectangle Selection");
		buttonManual.setActionCommand("Manual");
		buttonManual.addActionListener(this);
		buttonGroupAutoManual = new ButtonGroup();
		buttonPanelAutoManual.setLayout(new BoxLayout(buttonPanelAutoManual,
			BoxLayout.Y_AXIS));
		buttonGroupAutoManual.add(buttonAuto);
		buttonGroupAutoManual.add(buttonManual);

		buttonPanelAutoManual.add(buttonAuto);
		buttonPanelAutoManual.add(buttonManual);
		buttonPanelAutoManual.setBorder(new TitledBorder(BorderFactory
			.createEtchedBorder(), "Lane Selection", TitledBorder.LEADING,
			TitledBorder.BELOW_TOP, new Font("Sans", Font.PLAIN, 11)));

		buttonPanelFitType = new JPanel();
		buttonPanelFitType.setLayout(new BoxLayout(buttonPanelFitType,
			BoxLayout.Y_AXIS));
		buttonBands = new JRadioButton("Banded");
		buttonBands.setActionCommand("Banded");
		buttonContinuum = new JRadioButton("Continuum");
		buttonContinuum.setActionCommand("Continuum");

		buttonGroupFitType = new ButtonGroup();
		buttonGroupFitType.add(buttonBands);
		buttonGroupFitType.add(buttonContinuum);

		buttonPanelFitType.add(buttonBands);
		buttonPanelFitType.add(buttonContinuum);
		buttonPanelFitType.setBorder(new TitledBorder(BorderFactory
			.createEtchedBorder(), "Fit Type", TitledBorder.LEADING,
			TitledBorder.BELOW_TOP, new Font("Sans", Font.PLAIN, 11)));

		topButtonsPanel = new JPanel();
		topButtonsPanel.setLayout(new GridBagLayout());
		final GridBagConstraints c0 = new GridBagConstraints();
		c0.gridx = 0;
		c0.gridy = 0;
		c0.weightx = 0.1;
		c0.fill = GridBagConstraints.HORIZONTAL;

		topButtonsPanel.add(buttonPanelAutoManual, c0);
		c0.gridx = 1;
		c0.gridy = 0;
		c0.weightx = 0.0;
		c0.anchor = GridBagConstraints.EAST;
		topButtonsPanel.add(buttonPanelFitType, c0);

		sliderPanel = new JPanel();
		lanesPanel = new JPanel();
		lanesPanel.setLayout(new GridBagLayout());
		final GridBagConstraints c1 = new GridBagConstraints();
		lanesPanel.setBorder(new EmptyBorder(3, 5, 3, 3));
		labelNLanes = new JLabel("Number of Lanes");
		textNLanes = new JSpinner(new SpinnerNumberModel(nLanes, 1, 100, 1));
		textNLanes.setBorder(BorderFactory.createCompoundBorder(textNLanes
			.getBorder(), BorderFactory.createEmptyBorder(0, 2, 0, 2)));
		((JSpinner.DefaultEditor) textNLanes.getEditor()).getTextField().setColumns(
			textWidth);

		c1.gridx = 0;
		c1.gridy = 0;
		c1.weightx = 0.2;
		c1.fill = GridBagConstraints.HORIZONTAL;
		lanesPanel.add(labelNLanes, c1);

		c1.gridx = 1;
		c1.gridy = 0;
		c1.weightx = 0.0;
		c1.fill = GridBagConstraints.NONE;
		c1.ipady = 2;
		c1.ipadx = 2;
		lanesPanel.add(textNLanes, c1);
		sliderPanel.add(lanesPanel);

		sliderPanel.setLayout(new BoxLayout(sliderPanel, BoxLayout.Y_AXIS));
		sliderW = makeTitledSlider("Width ( " + lw + " px )", Color.black, 10, iw,
			lw);
		sliderPanel.add(sliderW);
		sliderH = makeTitledSlider("Height ( " + lh + " px )", Color.black, 10, ih,
			lh);
		sliderPanel.add(sliderH);
		sliderSp = makeTitledSlider("Space ( " + lsp + " px )", Color.black, 1, iw,
			lsp);
		sliderPanel.add(sliderSp);
		sliderHOff = makeTitledSlider("Horizontal Offset ( " + lhoff + " px )",
			Color.black, 0, iw - 1, lhoff);
		sliderPanel.add(sliderHOff);
		sliderVOff = makeTitledSlider("Vertical Offset ( " + lvoff + " px )",
			Color.black, 0, ih - 1, lvoff);
		sliderPanel.add(sliderVOff);

		buttonPanel = new JPanel();
		buttonPanel.setLayout(new FlowLayout(FlowLayout.TRAILING));
		buttonFit = new JButton("Fit");

		buttonPanel.add(buttonFit);
		buttonOpenFolder = new JButton("Open Data Folder");
		buttonOpenFolder.setToolTipText(savePath);
		buttonPanel.add(buttonOpenFolder);
		buttonClose = new JButton("Close");
		buttonPanel.add(buttonClose);

		// Settings Panel
		settingsPanel = new JPanel();
		settingsPanel.setLayout(new GridBagLayout());
		settingsPanel.setBorder(new EmptyBorder(2, 2, 2, 2));
		final GridBagConstraints c2 = new GridBagConstraints();
		
		// JSpinners
		labelDegBG = new JLabel("Polynomial Degree");
		labelPolyDerivative = new JLabel(
			"Max Polynomial Derivative (grayvalue/px)");
		labelTolPK = new JLabel("Peak Tolerance (fraction)");
		labelAreaDrift = new JLabel("Area Drift ");
		labelSDDrift = new JLabel("SD Drift ");

		int minDegBG = -1; int maxDegBG = 15;
		degBG = degBG < minDegBG ? minDegBG : degBG;
		degBG = degBG > maxDegBG ? maxDegBG : degBG;
		textDegBG = new JSpinner(new SpinnerNumberModel(degBG, minDegBG, maxDegBG, 1));
		textDegBG.setBorder(BorderFactory.createCompoundBorder(textDegBG
			.getBorder(), BorderFactory.createEmptyBorder(0, 2, 0, 2)));
		((JSpinner.DefaultEditor) textDegBG.getEditor()).getTextField().setColumns(
			textWidth);

		double minPolyDerivative = 0.00; double maxPolyDerivative = 10.0;
		polyDerivative = polyDerivative < minPolyDerivative ? minPolyDerivative : polyDerivative;
		polyDerivative = polyDerivative > maxPolyDerivative ? maxPolyDerivative : polyDerivative;
		textPolyDerivative = new JSpinner(new SpinnerNumberModel(polyDerivative,
			minPolyDerivative, maxPolyDerivative, 0.01));
		textPolyDerivative.setBorder(BorderFactory.createCompoundBorder(
			textPolyDerivative.getBorder(), BorderFactory.createEmptyBorder(0, 2, 0,
				2)));
		((JSpinner.DefaultEditor) textPolyDerivative.getEditor()).getTextField()
			.setColumns(textWidth);
		
		double minTolPK = 0.01; double maxTolPK = 1.00;
		tolPK = tolPK < minTolPK ? minTolPK : tolPK;
		tolPK = tolPK > maxTolPK ? maxTolPK : tolPK;
		textTolPK = new JSpinner(new SpinnerNumberModel(tolPK, minTolPK, maxTolPK, 0.005));
		textTolPK.setBorder(BorderFactory.createCompoundBorder(textTolPK
			.getBorder(), BorderFactory.createEmptyBorder(0, 2, 0, 2)));
		((JSpinner.DefaultEditor) textTolPK.getEditor()).getTextField().setColumns(
			textWidth);

		double minAreaDrift = 0.001; double maxAreaDrift = 1.00;
		areaDrift = areaDrift < minAreaDrift ? minAreaDrift : areaDrift;
		areaDrift = areaDrift > maxAreaDrift ? maxAreaDrift : areaDrift;
		textAreaDrift = new JSpinner(
				new SpinnerNumberModel(areaDrift, minAreaDrift, maxAreaDrift, 0.001));
		textAreaDrift.setEditor(new JSpinner.NumberEditor(textAreaDrift, "###.###"));
		textAreaDrift.setBorder(BorderFactory.createCompoundBorder(textAreaDrift
			.getBorder(), BorderFactory.createEmptyBorder(0, 2, 0, 2)));
		((JSpinner.DefaultEditor) textAreaDrift.getEditor()).getTextField()
			.setColumns(textWidth);
		textSDDrift = new JSpinner(new SpinnerNumberModel(sdDrift, 1.0, 5.0,
			0.1));
		textSDDrift.setEditor(new JSpinner.NumberEditor(textSDDrift, "###.#"));
		textSDDrift.setBorder(BorderFactory.createCompoundBorder(textSDDrift
			.getBorder(), BorderFactory.createEmptyBorder(0, 2, 0, 2)));
		((JSpinner.DefaultEditor) textSDDrift.getEditor()).getTextField()
			.setColumns(textWidth);

		chkBoxBands = new JCheckBox("Show Bands");

		final String[] lanes = new String[1 + nLanes];
		lanes[0] = "Select Ladder Lane";
		for (int i = 1; i < 1 + nLanes; i++)
			lanes[i] = "Lane " + (i);
		cmbBoxLadderLane = new JComboBox<>(lanes);

		cmbBoxDist = new JComboBox<>(distStr);
		cmbBoxLadderType = new JComboBox<>(ladderStr);

		buttonEditPeaks = new JToggleButton("Edit Custom Peaks");
		buttonResetCustomPeaks = new JButton("Reset Custom Peaks");

		c2.fill = GridBagConstraints.HORIZONTAL;
		c2.ipadx = 2; c2.ipady = 2;
		c2.gridx = 0; c2.gridy = 0;
		c2.weightx = 0.1;
		settingsPanel.add(labelDegBG, c2);
		c2.gridx = 0; c2.gridy = 1;
		c2.weightx = 0.1;
		settingsPanel.add(labelPolyDerivative, c2);
		c2.gridx = 0; c2.gridy = 2;
		c2.weightx = 0.1;
		settingsPanel.add(labelTolPK, c2);
		c2.gridx = 0; c2.gridy = 3;
		c2.weightx = 0.1;
		settingsPanel.add(labelAreaDrift, c2);
		c2.gridx = 0; c2.gridy = 4;
		c2.weightx = 0.1;
		settingsPanel.add(labelSDDrift, c2);

		c2.fill = GridBagConstraints.NONE;
		c2.gridx = 1; c2.gridy = 0;
		c2.weightx = 0.0;
		settingsPanel.add(textDegBG, c2);
		c2.gridx = 1; c2.gridy = 1;
		c2.weightx = 0.0;
		settingsPanel.add(textPolyDerivative, c2);
		c2.gridx = 1; c2.gridy = 2;
		c2.weightx = 0.0;
		settingsPanel.add(textTolPK, c2);
		c2.gridx = 1; c2.gridy = 3;
		c2.weightx = 0.0;
		settingsPanel.add(textAreaDrift, c2);
		c2.gridx = 1; c2.gridy = 4;
		c2.weightx = 0.0;
		settingsPanel.add(textSDDrift, c2);

		c2.gridx = 0; c2.gridy = 5;
		c2.gridwidth = 1; c2.gridheight = 1;
		c2.fill = GridBagConstraints.HORIZONTAL;
		settingsPanel.add(chkBoxBands, c2);

		c2.gridx = 0; c2.gridy = 6;
		c2.gridwidth = 3; c2.gridheight = 1;
		settingsPanel.add(cmbBoxLadderLane, c2);

		c2.gridx = 0; c2.gridy = 7;
		settingsPanel.add(cmbBoxLadderType, c2);

		c2.gridx = 0; c2.gridy = 8;
		settingsPanel.add(cmbBoxDist, c2);

		c2.gridx = 0; c2.gridy = 9;
		c2.gridwidth = 1;
		settingsPanel.add(buttonEditPeaks, c2);

		c2.gridx = 0; c2.gridy = 10;
		settingsPanel.add(buttonResetCustomPeaks, c2);

		c2.gridx = 0; c2.gridy = 11;
		c2.gridwidth = GridBagConstraints.REMAINDER;
		c2.gridheight = 1;
		c2.weighty = 1.0;
		settingsPanel.add(new JLabel(), c2);

		dialogPanel = new JPanel();
		dialogPanel.setBackground(Color.darkGray);
		dialogPanel.setLayout(new BorderLayout());
		dialogPanel.add(topButtonsPanel, BorderLayout.NORTH);
		dialogPanel.add(sliderPanel, BorderLayout.CENTER);
		dialogPanel.add(settingsPanel, BorderLayout.EAST);
		dialogPanel.add(buttonPanel, BorderLayout.SOUTH);

		// Initial status of all dialog components
		final boolean stateLoaded = !savedRects.isEmpty() || ladder != null;
		int ladderType = 0;
		if (stateLoaded && ladder != null) ladderType = ladder.getType();
		if (auto) {
			buttonAuto.doClick();
		}
		else {
			if (!stateLoaded || rois.size() == 0) { // go back to AUTO
				auto = true;
				prefs.putBoolean(AUTO, auto);
				buttonAuto.doClick();
			}
			else { // Go to MANUAL
				buttonManual.doClick();
			}
		}
		buttonAuto.setSelected(auto);
		buttonManual.setSelected(!auto);
		final boolean continuum = savedState != null && savedState.continuum;
		buttonBands.setSelected(!continuum);
		buttonContinuum.setSelected(continuum);
		fitter.setFitMode(continuum ? Fitter.continuumMode : Fitter.bandMode);
		// Only the user changes Show Bands; it is kept with the saved state
		chkBoxBands.setSelected(savedState != null && savedState.showBands);
		chkBoxBands.setEnabled(false);
		cmbBoxLadderType.setSelectedIndex(ladderType);
		cmbBoxLadderType.setEnabled(ladder != null);
		updateLadderLane();
		updateLadderType();
		cmbBoxDist.setSelectedIndex(savedState == null ? 0 : savedState.dist);
		cmbBoxDist.setEnabled(continuum);
		loadDistribution(false);

		buttonEditPeaks.setEnabled(false);
		buttonResetCustomPeaks.setEnabled(false);

		fitter.setDegBG(degBG);
		fitter.setPolyDerivative(polyDerivative);
		fitter.setTolPK(tolPK);
		fitter.setAreaDrift(areaDrift);
		fitter.setSDDrift(sdDrift);

		// Add this class to all components as listener here for easy reference
		imp.getCanvas().addMouseMotionListener(this);
		imp.getCanvas().addMouseListener(this);
		imp.getCanvas().addMouseWheelListener(this);
		imp.getWindow().addMouseListener(this);
		imp.getWindow().addMouseMotionListener(this);
		imp.getWindow().addWindowListener(this);

		textNLanes.addChangeListener(this);
		textDegBG.addChangeListener(this);
		textPolyDerivative.addChangeListener(this);
		textTolPK.addChangeListener(this);
		textAreaDrift.addChangeListener(this);
		textSDDrift.addChangeListener(this);
		
		chkBoxBands.addActionListener(this);
		buttonBands.addActionListener(this);
		buttonContinuum.addActionListener(this);

		cmbBoxDist.addActionListener(this);
		cmbBoxLadderLane.addActionListener(this);
		cmbBoxLadderType.addActionListener(this);

		buttonEditPeaks.addActionListener(this);
		buttonResetCustomPeaks.addActionListener(this);
		buttonFit.addActionListener(this);
		buttonOpenFolder.addActionListener(this);
		buttonClose.addActionListener(this);

		frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
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
		reDrawROIs(imp, "none");
		imp.killRoi();
		redoProfilePlots();
		restoring = false;
		if (savedState != null && savedState.fitDone) {
			log.info("Repeating the last fit ...");
			SwingUtilities.invokeLater(this::runFit);
		}
	}

	private boolean askUser(final String question) {
		final GenericDialog gd = new GenericDialog("WARNING!");
		gd.addMessage(question);
		gd.showDialog();
		if (gd.wasOKed()) return true;
		return false;
	}

	/** Keeps the Edit Custom Peaks button and the plots' mode in step */
	private void setEditPeaks(final boolean on) {
		buttonEditPeaks.setSelected(on);
		plotter.setPlotMode(on ? Plotter.editPeaksMode : Plotter.regMode);
	}

	/** Warns that the fit will be lost, unless the user opted out */
	private boolean askFitReset() {
		if (prefs.getBoolean(SKIPFITWARNING, false)) return true;
		final GenericDialog gd = new GenericDialog("WARNING!");
		gd.addMessage(warningFit);
		gd.addCheckbox("Don't show this again", false);
		gd.showDialog();
		if (!gd.wasOKed()) return false;
		if (gd.getNextBoolean()) prefs.putBoolean(SKIPFITWARNING, true);
		return true;
	}

	private Peak askPeak(final int lane, final double y, final double a) {
		String title, message;
		final String fwhmString = "FWHM = \u03C3 * (2 * \u221A (2 * ln(2))";
		final double fwhmValue = 5.0;

		title = "ADD PEAK";
		message = "If the new peak is located less than " +
			(int) Fitter.peakDistanceTol +
			" px away from an exisitng peak,\nthe existing peak will be replaced with the new one.";

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
		final String message =
			"Select the plots from which the custom peaks must be removed";
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

	private int[] askLadderRange() {
		final String title = "LADDER RANGE";
		final String message = "Select the first and last ladder bands to consider";
		final GenericDialog gd = new GenericDialog(title);
		final String ladderFirst = ladder.getStrings()[ladder.getRange()[0]];
		final String ladderLast = ladder.getStrings()[ladder.getRange()[1]];
		gd.addMessage(message);
		gd.addChoice("First Band", ladder.getStrings(), ladderFirst);
		gd.addChoice("Last Band", ladder.getStrings(), ladderLast);
		gd.showDialog();
		if (gd.wasOKed()) {
			final String[] bands = { gd.getNextChoice(), gd.getNextChoice() };
			final int[] range = new int[2];
			final String[] ladderBands = ladder.getStrings();
			for (final int i : new int[] { 0, 1 })
				for (int j = 0; j < ladderBands.length; j++) {
					if (bands[i].equals(ladderBands[j])) range[i] = j;
				}
			// Bands can be picked in either order
			if (range[0] > range[1]) {
				final int tmp = range[0];
				range[0] = range[1];
				range[1] = tmp;
			}
			return range;
		}
		return ladder.getRange();
	}

	private JSlider makeTitledSlider(final String string, final Color color,
		final int minVal, final int maxVal, final int val)
	{
		final JSlider slider = new JSlider(SwingConstants.HORIZONTAL, minVal,
			maxVal, Math.min(Math.max(val, minVal), maxVal));
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
			// Housekeeping: Sort the rois based on POSITION and remove null
			// elements
			final Comparator<Roi> roiNameComparator = new Comparator<>() {

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
				if (r == null) break;
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
			if (lh1 > rh || lw1 > rw) labelRoi.setLocation(x0 + 0.1 * lw1, y0 + rh +
				1.1 * lh1);
			else labelRoi.setLocation(x0 + 0.2 * lw1, y0 + rh - 1.1 * lh1);
			if (roi.getName().equals(roiName)) {
				roi.setStrokeColor(MainDialog.selected);
				roi.setStrokeWidth(3);
				if (buttonManual.isSelected() && !roiName.equals("none")) {
					imgPlus.setRoi(roi);
				}
			}
			else if (roi.getName().equals(ladderLaneStr)) {
				roi.setStrokeWidth(3);
				roi.setStrokeColor(MainDialog.reference);
				labelRoi.setFillColor(MainDialog.reference);
			}
			else {
				roi.setStrokeWidth(1);
				roi.setStrokeColor(MainDialog.unselected);
				labelRoi.setFillColor(MainDialog.unselected);
			}
			overlay.add(roi);
			overlay.add(labelRoi);
			if (buttonAuto.isSelected()) imgPlus.killRoi();
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

	/**
	 * Custom peaks keep their position on the gel when a lane changes: those
	 * now outside the lane are parked, the others take the new profile height
	 */
	private List<Peak> snapCustomPeaks(final int ln) {
		DataSeries profile = null;
		for (final DataSeries d : plotter.getProfiles()) {
			if (d.getLane() == ln) profile = d;
		}
		if (profile == null) return fitter.getCustomPeaks(ln);
		final PolynomialSplineFunction f = new LinearInterpolator().interpolate(
			profile.getX().toArray(), profile.getY().toArray());

		// Parked peaks inside the lane again go back to the fitter
		final Iterator<Peak> itParked = parkedPeaks.iterator();
		while (itParked.hasNext()) {
			final Peak p = itParked.next();
			if (p.getLane() == ln && f.isValidPoint(p.getMean())) {
				fitter.addCustomPeak(p);
				itParked.remove();
			}
		}

		final List<Peak> peaks = fitter.getCustomPeaks(ln);
		final Iterator<Peak> it = peaks.iterator();
		while (it.hasNext()) {
			final Peak p = it.next();
			if (f.isValidPoint(p.getMean())) p.setNorm(f.value(p.getMean()));
			else {
				fitter.removeCustomPeak(p);
				parkedPeaks.add(p);
				it.remove();
			}
		}
		return peaks;
	}

	/** Custom peaks in the fitter and those parked outside the lanes */
	private List<Peak> allCustomPeaks() {
		final List<Peak> all = new ArrayList<>(fitter.getAllCustomPeaks());
		all.addAll(parkedPeaks);
		return all;
	}

	/** Keeps the custom peaks while the fitter is reset for new lanes */
	private void resetFitterKeepCustomPeaks() {
		final List<Peak> keep = restoring ? savedCustomPeaks : allCustomPeaks();
		fitter.resetAllFitter();
		parkedPeaks = new ArrayList<>();
		for (final Peak p : keep)
			fitter.addCustomPeak(p);
	}

	private void manualROIsFromSavedRects() {
		rois = new ArrayList<>();
		int i = 1;
		for (final Rectangle r : savedRects) {
			final Roi roi = new Roi(r);
			roi.setName("Lane " + i++);
			rois.add(roi);
		}
	}

	private void redoProfilePlots() {
		fitDone = false;
		chkBoxBands.setEnabled(false);
		plotter.resetData();

		try {
			for (final Roi r : rois) {
				final int ln = Integer.parseInt(r.getName().substring(5));
				final Rectangle rect = r.getBounds();
				// Do not ROIs that are outside the image
				if (rect.getMinX() < 0.95 * iw && rect.getMinY() < 0.95 * ih) {
					plotter.updateProfile(r);
					final RealVector empty = new ArrayRealVector();
					final DataSeries d = new DataSeries("Custom Points", ln,
						DataSeries.CUSTOMPEAKS, empty, empty, Plotter.vMarkerEditPeakColor);
					for (final Peak p : snapCustomPeaks(ln)) {
						d.addOrUpdate(p.getMean(), p.getNorm());
					}
					d.addChangeListener(this);
					plotter.addDataSeries(d);
					plotter.updatePlot(r);
				}
			}
			plotter.reloadTabs();
		}
		catch (final ConcurrentModificationException e) {
			log.info(rois.size());
		}
	}

	private double[][] readDistFile(final String filename) {
		String line = null;
		final List<Integer> fragmentLength = new ArrayList<>();
		final List<Integer> fragmentFrequency = new ArrayList<>();

		final URL url = MainDialog.class.getClassLoader().getResource(filename);
		if (url == null) {
			IJ.error("Fragment Distribution", "The plugin does not contain " +
				filename + ".");
			return null;
		}
		log.info("Loading " + url.getPath() + " ...");
		try (BufferedReader buffer = new BufferedReader(new InputStreamReader(url
			.openStream())))
		{
			int lineNumber = 0;
			while (true) {
				line = buffer.readLine();
				if (line == null) break;
				lineNumber++;

				final String[] words = line.split("\t");
				try {
					final int length = Integer.parseInt(words[0].trim());
					final int frequency = Integer.parseInt(words[1].trim());
					fragmentLength.add(length);
					fragmentFrequency.add(frequency);
				}
				catch (final NumberFormatException e1) {
					if (lineNumber == 1) continue; // the header line
					log.info("Invalid token: " + words[0].trim() + " " + words[1].trim());
				}
			}
			buffer.close();
		}
		catch (final Exception e) {
			log.error("Could not read " + filename, e);
		}
		if (fragmentLength.size() != fragmentFrequency.size()) return null;
		final double[][] out = new double[fragmentLength.size()][3];
		int count = 0;
		for (int i = 0; i < fragmentLength.size(); i++) {
			out[i][0] = fragmentFrequency.get(i);
			out[i][1] = fragmentLength.get(i);
			out[i][2] = Ladder.molecularWeight(fragmentLength.get(i));
			count = count + fragmentFrequency.get(i);
		}
		for (int i = 0; i < fragmentLength.size(); i++) {
			out[i][0] = out[i][0] / count;
		}
		// 3 columns: Relative Frequency, Length, MW
		return out;
	}

	private void resetCustomPeaks(final int lane) {
		fitter.resetCustomPeaks(lane);
		parkedPeaks.removeIf(p -> p.getLane() == lane);
		DataSeries d = plotter.getPlotsCustomPeaks(lane);
		if (d == null) {
			final RealVector empty = new ArrayRealVector();
			d = new DataSeries("Custom Peaks", lane, DataSeries.CUSTOMPEAKS, empty,
				empty, Plotter.vMarkerEditPeakColor);
			plotter.addDataSeries(d);
		}
		else {
			d.removeChangeListener(this);
			d.clear();
		}
		d.addChangeListener(this);
	}

	/** Shows the folder with the saved state, results and plots */
	private void openDataFolder() {
		final File dir = new File(savePath);
		dir.mkdirs();
		try {
			if (!Desktop.isDesktopSupported() || !Desktop.getDesktop().isSupported(
				Desktop.Action.OPEN)) throw new IOException("not supported");
			Desktop.getDesktop().open(dir);
		}
		catch (final IOException | RuntimeException e) {
			new MessageDialog(frame, "Data Folder", "The data is saved in\n" +
				savePath);
		}
	}

	private void displayLog() {
		final String logOutput = fitter.getSummary() + "<p>Data folder: " + savePath
			.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;") + "</p>";
		final HTMLDialog logWindow = new HTMLDialog("LOG", logOutput, false);

		if (buttonContinuum.isSelected()) {
			final JTabbedPane distributionsPane = new JTabbedPane();
			final Component[] c = logWindow.getContentPane().getComponents();

			for (final int i : getAllLaneNumbers()) {
				if (i != ladderLaneInt) {
					final String name = String.format("Lane %1$d", i);
					final XYSeriesCollection dataset = new XYSeriesCollection();
					dataset.addSeries(fitter.getFittedDistribution(i));
					final JFreeChart chart = ChartFactory.createXYLineChart("",
						"Length (Base Pairs)", "Normalized, Length-weighed, Intensity",
						dataset, PlotOrientation.VERTICAL, false, true, false);
					final XYPlot plot = chart.getXYPlot();
					plot.getDomainAxis().setInverted(true);
					plot.setBackgroundPaint(Color.white);
					plot.setDomainGridlinePaint(Color.lightGray);
					plot.setRangeGridlinePaint(Color.lightGray);
					plot.getRenderer().setSeriesPaint(0, Color.darkGray);
					plot.getRenderer().setSeriesStroke(0, new BasicStroke(2.0f));
					final ChartPanel p = new ChartPanel(chart);
					distributionsPane.addTab(name, p);
				}
			}
			final JSplitPane sp = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,
				distributionsPane, c[0]);
			logWindow.getContentPane().add(sp);
			final Dimension screenD = IJ.getScreenSize();
			final Dimension dialogD = logWindow.getSize();
			dialogD.width = (int) FastMath.min(0.80 * screenD.width, 1500);
			logWindow.setSize(dialogD);
			GUI.center(logWindow);
			sp.setDividerLocation(0.9);
		}
		final String file = impTitle + "_log.html";
		try (BufferedWriter out = new BufferedWriter(new FileWriter(savePath +
			file)))
		{
			out.write(logOutput);
			out.close();
		}
		catch (final IOException e) {
			log.info("Exception", e);
		}
	}

	private FitState currentState() {
		final FitState st = new FitState();
		st.auto = auto;
		st.nLanes = nLanes;
		st.lw = lw;
		st.lh = lh;
		st.lsp = lsp;
		st.lhoff = lhoff;
		st.lvoff = lvoff;
		st.ladderLane = ladderLaneInt;
		st.degBG = degBG;
		st.polyDerivative = polyDerivative;
		st.tolPK = tolPK;
		st.areaDrift = areaDrift;
		st.sdDrift = sdDrift;
		st.continuum = buttonContinuum.isSelected();
		st.dist = cmbBoxDist.getSelectedIndex();
		st.dlo = dlo;
		st.dhi = dhi;
		st.every = every;
		st.customPeaks = allCustomPeaks();
		st.fitDone = fitDone;
		st.showBands = chkBoxBands.isSelected();
		return st;
	}

	private boolean loadState() {
		rois = new ArrayList<>();
		savedRects = new ArrayList<>();
		savedCustomPeaks = new ArrayList<>();
		savedState = null;
		final String file = "saved-state.bak";
		final String fullPath = savePath + file;
		log.info("Loading " + fullPath + " ...");
		try (ObjectInputStream ois = new ObjectInputStream(new FileInputStream(
			fullPath)))
		{
			try {
				while (true) {
					final Object o = ois.readObject();
					if (o instanceof Rectangle) {
						savedRects.add((Rectangle) o);
					}
					else if (o instanceof Ladder) {
						ladder = (Ladder) o;
						fitter.setLadder(ladder.getMolecularWeights());
						if (cmbBoxLadderType != null) cmbBoxLadderType.setEnabled(true);
					}
					else if (o instanceof FitState) {
						savedState = (FitState) o;
						savedCustomPeaks = savedState.customPeaks;
					}
				}
			}
			catch (final Exception e) {
				/* Exit */ }
			ois.close();
			// Old files have no FitState; the mode comes from the preferences
			final boolean manual = savedState == null ? !auto : !savedState.auto;
			if (manual) manualROIsFromSavedRects();
			return true;
		}
		catch (final IOException e) {
			log.info("Could not find previously saved state");
			return false;
		}
	}

	private boolean saveState() {
		if (restoring) return true;
		final String file = "saved-state.bak";
		final String fullPath = savePath + file;
		new File(savePath).mkdirs();
		log.info("Saving to " + fullPath + " ...");
		try (ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(
			fullPath)))
		{
			// Keep the manual ROIs while in AUTO mode
			if (!auto) {
				savedRects = new ArrayList<>();
				for (final Roi r : rois)
					savedRects.add(r.getBounds());
			}
			for (final Rectangle r : savedRects)
				oos.writeObject(r);
			oos.writeObject(currentState());
			if (ladder != null) {
				oos.writeObject(ladder);
			}
			return true;
		}
		catch (final IOException e) {
			log.error("Could not save " + fullPath, e);
			return false;
		}
	}

	private void updateLadderLane() {
		if (cmbBoxLadderLane.getItemCount() != rois.size() + 1) {
			if (ladderLaneInt > rois.size()) 
				ladderLaneInt = 0;
			final String str = cmbBoxLadderLane.getItemAt(0);
			cmbBoxLadderLane.setEnabled(false);
			cmbBoxLadderLane.removeAllItems();
			cmbBoxLadderLane.addItem(str);
			for (final Roi r : rois) {
				final String name = r.getName();
				cmbBoxLadderLane.addItem(name);
			}
			cmbBoxLadderLane.setEnabled(true);
		}
		
		if (cmbBoxLadderLane.getSelectedIndex() != ladderLaneInt)
			cmbBoxLadderLane.setSelectedIndex(ladderLaneInt);
		
		fitter.setReferenceLane(ladderLaneInt);
		plotter.setReferencePlot(ladderLaneInt);
		prefs.putInt(LADDERLANEINT, ladderLaneInt);
	}

	private void updateLadderType() {
		if (ladder == null) return;
		cmbBoxLadderType.removeAllItems();
		for (int i = 0; i < ladderStr.length; i++) {
			String s = "";
			if (i == ladder.getType() && i != Ladder.CUSTOM) {
				s = ladderStr[i] + " [" + ladder.getStrings()[ladder.getRange()[0]] +
					" - " + ladder.getStrings()[ladder.getRange()[1]] + "]";
			}
			else s = ladderStr[i];

			cmbBoxLadderType.addItem(s);
		}

		cmbBoxLadderType.setSelectedIndex(ladder.getType());
		final RealVector ladderweights = ladder.getMolecularWeights();
		fitter.setLadder(ladderweights);
		saveState();
	}

	private void cleanupAndClose() {
		saveState();
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

		saveState();
		plotter.savePlots(savePath);
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
			if (roi.getName().equals(title)) return roi;
		}
		return null;
	}

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

	@Override
	public synchronized void stateChanged(final ChangeEvent e) {
		final Object o = e.getSource();
		final List<Object> fitDisruptors = new ArrayList<>(Arrays.asList(sliderW,
			sliderH, sliderSp, sliderHOff, sliderVOff, textNLanes));

		if (fitDisruptors.contains(o)) {
			if (fitDone && !askFitReset()) {
				textNLanes.setValue(nLanes);
				return;
			}
		}
		// ROI SLIDERS
		if (o == sliderW) {
			lw = sliderW.getValue();
			prefs.putInt(LW, lw);
			final String str = "Width ( " + lw + " px )";
			setSliderTitle(sliderW, str);
			sliderUpdate();
		}
		else if (o == sliderH) {
			lh = sliderH.getValue();
			prefs.putInt(LH, lh);
			final String str = "Height ( " + lh + " px )";
			setSliderTitle(sliderH, str);
			sliderUpdate();
		}
		else if (o == sliderSp) {
			lsp = sliderSp.getValue();
			prefs.putInt(LSP, lsp);
			final String str = "Spacing ( " + lsp + " px )";
			setSliderTitle(sliderSp, str);
			sliderUpdate();
		}
		else if (o == sliderHOff) {
			lhoff = sliderHOff.getValue();
			prefs.putInt(LHOFF, lhoff);
			final String str = "Horizontal Offset ( " + lhoff + " px )";
			setSliderTitle(sliderHOff, str);
			sliderUpdate();
		}
		else if (o == sliderVOff) {
			lvoff = sliderVOff.getValue();
			prefs.putInt(LVOFF, lvoff);
			final String str = "Vertical Offset ( " + lvoff + " px )";
			setSliderTitle(sliderVOff, str);
			sliderUpdate();
		}
		// SPINNERS
		if (o == textNLanes) {
			nLanes = (int) textNLanes.getModel().getValue();
			prefs.putInt(NLANES, nLanes);
			ladderLaneInt = ladderLaneInt <= nLanes ? ladderLaneInt : noLadderLane;
			ladderLaneStr = ladderLaneInt == 0 ? "none" : "Lane " + ladderLaneInt;
			sliderUpdate();
			updateLadderLane();
		}
		else if (o == textDegBG) {
			degBG = (int) textDegBG.getModel().getValue();
			fitter.setDegBG(degBG);
			prefs.putInt(DEGBG, degBG);
		}
		else if (o == textPolyDerivative) {
			polyDerivative = (double) textPolyDerivative.getModel().getValue();
			fitter.setPolyDerivative(polyDerivative);
			prefs.putDouble(POLYDERIVATIVE, polyDerivative);
		}
		else if (o == textTolPK) {
			tolPK = (double) textTolPK.getModel().getValue();
			fitter.setTolPK(tolPK);
			prefs.putDouble(TOLPK, tolPK);
		}
		else if (o == textAreaDrift) {
			areaDrift = (double) textAreaDrift.getModel().getValue();
			fitter.setAreaDrift(areaDrift);
			prefs.putDouble(AREADRIFT, areaDrift);
		}
		else if (o == textSDDrift) {
			sdDrift = (double) textSDDrift.getModel().getValue();
			fitter.setSDDrift(sdDrift);
			prefs.putDouble(SDDRIFT, sdDrift);
		}
	}

	private void sliderUpdate() {
		resetAutoROIs();
		redoProfilePlots();
		reDrawROIs(imp, "none");
		if (plotter.getProfiles().size() < 1) return;
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
							log.info("Custom Peak not added: " + d.getLane() + ", " + n1.get(
								m1.indexOf(m)) + ", " + m);
						}
						else {
							fitter.addCustomPeak(peak);
							log.info("Custom Peak added: " + d.getLane() + ", " + n1.get(m1
								.indexOf(m)) + ", " + m);
						}
					}
				}
			}
			else if (m1.size() < m2.size()) { // REMOVE
				for (final double m : m2) {
					if (!m1.contains(m)) {
						if (fitter.removeCustomPeak(new Peak(lane, n2.get(m2.indexOf(m)),
							m)))
						{
							log.info("Custom Peak Removed: " + n2.get(m2.indexOf(m)) + ", " +
								m);
						}
					}
				}
			}
			final int peakNumber = fitter.getCustomPeaks(lane).size();
			if (d.getItemCount() != peakNumber) {
				log.error("Wrong custom peaks count: " + d.getItemCount() + ", " +
					peakNumber);
			}
			saveState();
		}
	}

	/** Fits the ladder lane, then the other lanes with the current settings */
	private void runFit() {
		if (ladderLaneInt == noLadderLane) {
			askUser("Reference ladder not selected");
			return;
		}

		if (cmbBoxLadderType.getSelectedItem().equals("Select Ladder Type")) {
			askUser("Ladder type not selected");
			return;
		}

		// Remove everything in the plots, but the profile
		plotter.setSelected(MainDialog.noLaneSelected);
		plotter.removeFit();
		plotter.removeVerticalMarkers();

		// Reset fit
		for (final int l : getAllLaneNumbers())
			fitter.resetFit(l);
		fitter.setInputData(plotter.getProfiles());
		fitDone = false;

		chkBoxBands.setEnabled(true);
		buttonEditPeaks.setEnabled(true);
		buttonResetCustomPeaks.setEnabled(true);
		buttonResetCustomPeaks.setSelected(false);

		// Start new fit, Ladder first
		List<DataSeries> fitted = new ArrayList<>();
		fitter.setFitMode(Fitter.bandMode);
		fitted = fitter.doFit(ladderLaneInt);
		fitter.updateResultsTable(savePath);

		// Update Reference plot with peaks and vertical markers
		plotter.addDataSeries(fitted);
		plotter.addVerticalMarkers(fitter.getFittedPeaks(ladderLaneInt));
		for (final int i : getAllLaneNumbers()) {
			plotter.updatePlot(i);
		}
		final List<Peak> peaks = fitter.getFittedPeaks(ladderLaneInt);
		final int l1 = peaks.size();
		final int l2 = ladder.getRange()[1] - ladder.getRange()[0] + 1;

		// If mismatched, return
		if (l1 != l2) {
			String message =
				"The number of peaks detected does not match the ladder range";
			message = message + "\n Peaks detected: " + l1;
			message = message + "\n Bands in range: " + l2;
			new MessageDialog(frame, "WARNING", message);
			plotter.savePlots(savePath);
			new FileSaver(imp).saveAsTiff(savePath + impTitle + ".tif");
			return;
		}

		// If ladder peaks correct, rename vertical markers
		for (final Peak p : peaks) { // Rename ladder peaks with band names
			p.setName(ladder.getStrings()[ladder.getRange()[0] + peaks.indexOf(p)]);
		}
		plotter.addVerticalMarkers(peaks);

		// Set ladder in Fitter
		fitter.setLadder(ladder.getMolecularWeights());
		// Fit the rest based on selected criterion
		if (buttonBands.isSelected()) fitter.setFitMode(Fitter.bandMode);
		else if (buttonContinuum.isSelected()) {
			fitter.setFitMode(Fitter.continuumMode);
			if (cmbBoxDist.getSelectedItem().equals(
				"Select Fragment Distribution"))
			{
				new MessageDialog(frame, "WARNING!",
					"Fragment distribution not selected");
				return;
			}
		}

		final List<Integer> otherLanes = new ArrayList<>();
		for (final int i : getAllLaneNumbers()) {
			if (i != ladderLaneInt) otherLanes.add(i);
		}
		fitted.addAll(fitter.doFit(otherLanes));
		fitter.updateResultsTable(savePath);
		fitDone = true;

		plotter.removeFit();
		plotter.addDataSeries(fitted);
		for (final int i : getAllLaneNumbers()) {
			plotter.updatePlot(i);
		}
		reDrawROIs(imp, "none"); // adds the bands to the ROIs

		plotter.savePlots(savePath);
		if (!plotter.isVisible()){
			plotter.pack();
			plotter.setVisible(true);
		}
		displayLog();
		new FileSaver(imp).saveAsTiff(savePath + impTitle + ".tif");
		saveState();
	}

	@Override
	public void actionPerformed(final ActionEvent e) {
		// Buttons
		// ----------------------------------------------------------------
		if (e.getSource().equals(buttonFit)) {
			if (fitDone && !askFitReset()) return;
			runFit();
		}

		if (e.getSource().equals(buttonAuto)) {
			if (!auto) { // NOT already AUTO
				if (!saveState()) {
					final String warn =
						"When switching to AUTO mode, the current lane selections will be reset!";
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
			resetFitterKeepCustomPeaks();
			redoProfilePlots();
			fitter.setInputData(plotter.getProfiles());
			saveState();
		}

		if (e.getSource().equals(buttonManual)) {
			auto = false;
			prefs.putBoolean(AUTO, auto);
			setSliderPanelEnabled(false);
			manualROIsFromSavedRects();
			if (rois.size() == 0) { // Use the AUTO rois as a start
				resetAutoROIs();
				reDrawROIs(imp, "none");
			}
			if (plotter == null || fitter == null) return;
			resetFitterKeepCustomPeaks();
			reDrawROIs(imp, "none");
			redoProfilePlots();
			fitter.setInputData(plotter.getProfiles());
			saveState();
		}

		if (e.getSource().equals(buttonBands)) {
			fitter.setFitMode(Fitter.bandMode);
			cmbBoxDist.setEnabled(false);
		}
		else if (e.getSource().equals(buttonContinuum)) {
			fitter.setFitMode(Fitter.continuumMode);
			cmbBoxDist.setEnabled(true);
		}

		if (e.getSource().equals(buttonEditPeaks)) {
			setEditPeaks(buttonEditPeaks.isSelected());
			for (final int i : getAllLaneNumbers())
				plotter.updatePlot(i);
		}

		if (e.getSource().equals(buttonResetCustomPeaks)) {
			final int[] lanes = askResetCustomPoints(new boolean[rois.size() + 1]);
			for (final int i : lanes) {
				resetCustomPeaks(i);
			}
			saveState();
			reDrawROIs(imp, "none");
		}

		if (e.getSource().equals(buttonOpenFolder)) {
			openDataFolder();
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
			saveState();
		}

		// ComboBoxes
		// -------------------------------------------------------------
		if (e.getSource().equals(cmbBoxLadderLane)) {
			final int laneInt = cmbBoxLadderLane.getSelectedIndex();
			if (laneInt == -1) return;
			if (laneInt == 0) {
				ladderLaneInt = MainDialog.noLadderLane;
				cmbBoxLadderType.setEnabled(false);
				cmbBoxDist.setEnabled(false);
			}
			else {
				ladderLaneInt = laneInt;
				cmbBoxLadderType.setEnabled(true);
				cmbBoxDist.setEnabled(true);
			}
			ladderLaneStr = ladderLaneInt == 0 ? "none" : "Lane " + ladderLaneInt;
			updateLadderLane();
			saveState();
			reDrawROIs(imp, "none");
		}

		if (e.getSource().equals(cmbBoxLadderType)) {
			final int type = cmbBoxLadderType.getSelectedIndex();
			if (type != 0) {
				final boolean loaded;
				if (ladder == null) {
					ladder = Ladder.create(type);
					loaded = ladder != null;
				}
				else loaded = ladder.getType() == type || ladder.setType(type);

				if (loaded) {
					ladder.setRange(askLadderRange());
					updateLadderType();
				}
				else {
					// No custom ladder file loaded: go back to the previous selection
					cmbBoxLadderType.setSelectedIndex(ladder == null ? 0 : ladder
						.getType());
				}
			}
		}

		if (e.getSource().equals(cmbBoxDist)) {
			loadDistribution(true);
		}
	}

	/** Sets the fitter's fragment distribution from the selection */
	private void loadDistribution(final boolean ask) {
		if (cmbBoxDist.getSelectedIndex() == 0) return;
		if (cmbBoxDist.getSelectedItem().equals("Uniform")) {
			if (ask) {
				final String title = "DISTRIBUTION RANGE";
				final String message =
					"Enter the basepair-length range for the distribution.";
				final GenericDialog gd = new GenericDialog(title);
				gd.addMessage(message);
				gd.addNumericField("Lower", dlo, 0, 5, "bp");
				gd.addNumericField("Upper", dhi, 0, 5, "bp");
				gd.addNumericField("Every", every, 0, 5, "bp");
				gd.showDialog();
				if (!gd.wasOKed()) return;
				dlo = (int) gd.getNextNumber();
				dhi = (int) gd.getNextNumber();
				every = (int) gd.getNextNumber();
				prefs.putInt(DLO, dlo);
				prefs.putInt(DHI, dhi);
				prefs.putInt(EVERY, every);
			}
			final double[][] dist = new double[(dhi - dlo) / every + 1][3];
			final double f = 1.0 / (dhi - dlo + 1);

			for (int i = 0; i < dist.length; i++) {
				dist[i][0] = f;
				dist[i][1] = dhi - i * every;
				dist[i][2] = Ladder.molecularWeight(dist[i][1]);
			}
			fitter.setFragmentDistribution(dist);
		}
		else {
			final String filename = "sample-distributions/" + cmbBoxDist
				.getSelectedItem() + ".txt";
			final double[][] dist = readDistFile(filename);
			if (dist == null) cmbBoxDist.setSelectedIndex(0);
			else fitter.setFragmentDistribution(dist);
		}
	}

	@Override
	public void mouseDragged(final MouseEvent e) {
		// nothing to do
	}

	@Override
	public void mouseMoved(final MouseEvent e) {
		if (e.getSource() == imp.getCanvas()) {
			if (rois.size() == 0 || plotter == null) return;
			final int x = ((ImageCanvas) e.getSource()).offScreenX(e.getX());
			final int y = ((ImageCanvas) e.getSource()).offScreenY(e.getY());
			statusServ.showStatus("[" + x + ":" + y + "]");
			if (!selectionUpdate) { // Not dragging an ROI
				String roiCurrent = "none"; // None selected
				for (final Roi r : rois) {
					if (r.contains(x, y)) {
						roiCurrent = r.getName(); // This selected
					}
				}

				reDrawROIs(imp, roiCurrent);
				if (roiCurrent.equals("none") && imp.getRoi() != null) imp.killRoi();

				if (!(roiCurrent.equals(roiSelected))) {
					roiPreviouslySelected = roiSelected;
					roiSelected = roiCurrent;
					if (!roiCurrent.equals("none")) {
						final int ln = Integer.parseInt(roiCurrent.substring(5));
						plotter.setSelected(ln);
						plotter.setVLine(ln, y);
					}
					else {
						plotter.setSelected(MainDialog.noLaneSelected);
					}

					if (!roiPreviouslySelected.equals("none")) {
						final int roiPS = Integer.parseInt(roiPreviouslySelected.substring(
							5));
						plotter.removeVLine(roiPS);
					}
				}
				else { // Moving inside same ROI
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
					if (fitDone && !askFitReset()) return;
					final Iterator<Roi> roiIter = rois.iterator();
					while (roiIter.hasNext()) {
						final Roi roiRemove = roiIter.next();
						if (roiRemove.getName().equals(roiSelected)) {
							roiIter.remove();
						}
					}
					saveState();
					reDrawROIs(imp, "none");
					updateLadderLane();
					redoProfilePlots();
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
		// Modify Roi, mouseReleased not triggered when creating Roi
		if (buttonManual.isSelected()) {
			if (selectionUpdate) {
				final Roi roiNew = imp.getRoi();
				if (roiNew != null) {
					if (roiNew.getType() != Roi.RECTANGLE) return;
					// If ROI exists and larger that 100 px
					final Rectangle rN = roiNew.getBounds();
					if (rN.getWidth() * rN.getHeight() > 100) {
						if (roiNew.getName() != null) { // Existing ROI
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
						updateLadderLane();
					}
				}
				selectionUpdate = false;
				saveState();
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
			if (fitDone && !askFitReset()) return;
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

/**
 * A size ladder: its band names and sizes, and the range of bands visible in
 * the ladder lane. Built-in types (Hi-Lo, 100bp, Quick-Load, Tapestation) or
 * CUSTOM, loaded from a text file with one size in bp per line. It's saved in
 * saved-state.bak, so its fields must stay compatible (see SavedStateTest).
 */
class Ladder implements Serializable {

	// Pinned to the value computed before Ladder changed, so saved states still load
	private static final long serialVersionUID = -3262235272790955773L;

	private static final String[] hilo      = { "10 kbp", "8 kbp", "6 kbp", "4 kbp",
		"3 kbp", "2 kbp", "1.55 kbp", "1.4 kbp", "1 kbp", "750 bp", "500 bp",
		"400 bp", "300 bp", "200 bp", "100 bp", "50 bp" };
	private static final String[] bp100     = { "1.5 kbp", "1.2 kbp", "1 kbp",
		"900 bp", "800 bp", "700 bp", "600 bp", "500 bp", "400 bp", "300 bp",
		"200 bp", "100 bp" };
	private static final String[] quickload = { "48.5 kbp", "20 kbp", "15 kbp",
			"10 kbp", "8 kbp", "6 kbp", "5 kbp", "4 kbp", "3 kbp", "2 kbp", 
			"1.5 kbp",  "1 kbp",  "500 bp" };
	private static final String[] tapestation = { "1.5 kbp", "1 kbp", "700 bp",
			"500 bp", "400 bp", "300 bp", "200 bp", "100 bp",  "50 bp",  "25 bp" };
	
	private static final RealVector hilo_bp = new ArrayRealVector(new double[] {
		10000, 8000, 6000, 4000, 3000, 2000, 1550, 1400, 1000, 750, 500, 400, 300,
		200, 100, 50 });
	private static final RealVector bp100_bp = new ArrayRealVector(new double[] {
		1517, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100 });
	private static final RealVector quickload_bp = new ArrayRealVector(new double[] {
		48500, 20000, 15000, 10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 500 });
	private static final RealVector tapestation_bp = new ArrayRealVector(new double[] {
			1500, 1000, 700, 500, 400, 300, 200, 100, 50, 25 });
	private RealVector custom_bp = new ArrayRealVector();
	
	/** Average molecular weight of a base pair of double-stranded DNA, Da */
	private static final double DALTONS_PER_BP = 607.4;
	/** Molecular weight of the two ends of a double-stranded fragment, Da */
	private static final double END_DALTONS = 157.9;

	/** Molecular weight of a double-stranded DNA fragment, in Da */
	static double molecularWeight(final double bp) {
		return bp * DALTONS_PER_BP + END_DALTONS;
	}

	static final int HILO        = 1;
	static final int BP100       = 2;
	static final int QUICKLOAD   = 3;
	static final int TAPESTATION = 4;
	static final int CUSTOM      = 5;
	
	private int type;
	private int[] ladderRange;
	private String[] ladderStrings;

	/**
	 * Returns null if the type is CUSTOM and no ladder file could be loaded.
	 */
	static Ladder create(final int type) {
		final Ladder ladder = new Ladder();
		return ladder.setType(type) ? ladder : null;
	}

	/**
	 * Asks for a text file listing one band size (bp) per line. Returns false,
	 * leaving the ladder unchanged, if the user cancels or no sizes are found.
	 */
	private boolean askLadderFile() {
		final JFileChooser fc = new JFileChooser();
		if (fc.showOpenDialog(null) != JFileChooser.APPROVE_OPTION) return false;
		final File file = fc.getSelectedFile();
		RealVector bp = new ArrayRealVector();
		try (BufferedReader buffer = new BufferedReader(new InputStreamReader(
			new FileInputStream(file))))
		{
			String line;
			while ((line = buffer.readLine()) != null) {
				line = line.trim();
				if (line.isEmpty()) continue;
				try {
					bp = bp.append(Integer.parseInt(line));
				}
				catch (final NumberFormatException e1) {
					// Not a band size, e.g. a header line: skip it
				}
			}
		}
		catch (final IOException e) {
			IJ.error("Custom Ladder", "Could not read " + file.getName() + ":\n" +
				e.getMessage());
			return false;
		}
		if (bp.getDimension() == 0) {
			IJ.error("Custom Ladder", "No band sizes found in " + file.getName() +
				".\nThe file should list one whole number of base pairs per line.");
			return false;
		}

		custom_bp = bp;
		ladderStrings = new String[custom_bp.getDimension()];
		for (int w = 0; w < custom_bp.getDimension(); w++) {
			final double bases = custom_bp.getEntry(w);
			if (bases < 1000) {
				ladderStrings[w] = (int) bases + " bp";
			}
			else {
				ladderStrings[w] = bases / 1000 + " kbp";
			}
		}
		return true;
	}
	
	public RealVector getMolecularWeights() {
		RealVector bp = new ArrayRealVector();
		final int nel = ladderRange[1] - ladderRange[0] + 1;
		if (this.type == HILO) {
			bp = hilo_bp.getSubVector(ladderRange[0], nel);
		}
		else if (this.type == BP100) {
			bp = bp100_bp.getSubVector(ladderRange[0], nel);
		}
		else if (this.type == QUICKLOAD) {
			bp = quickload_bp.getSubVector(ladderRange[0], nel);
		}
		else if (this.type == TAPESTATION) {
			bp = tapestation_bp.getSubVector(ladderRange[0], nel);
		}
		else if (this.type == CUSTOM) {
			bp = custom_bp.getSubVector(ladderRange[0], nel);
		}
		final RealVector mw = bp.map(Ladder::molecularWeight);
		return mw;
	}

	public int[] getRange() {
		return this.ladderRange;
	}

	public void setRange(final int[] ladderRange) {
		this.ladderRange = ladderRange;
	}

	public String[] getStrings() {
		return this.ladderStrings;
	}

	public int getType() {
		return this.type;
	}

	/**
	 * Returns false, leaving the ladder unchanged, if the type is CUSTOM and no
	 * ladder file could be loaded.
	 */
	public boolean setType(final int type) {
		if (type == CUSTOM && !askLadderFile()) return false;
		this.type = type;

		if (type == HILO) ladderStrings = hilo;
		else if (type == BP100) ladderStrings = bp100;
		else if (type == QUICKLOAD) ladderStrings = quickload;
		else if (type == TAPESTATION) ladderStrings = tapestation;

		ladderRange = new int[] { 0, ladderStrings.length - 1 };
		return true;
	}
}

/**
 * Fitting settings of an image, saved with its lanes and ladder so that the
 * last fit can be repeated when the image is opened again
 */
class FitState implements Serializable {

	private static final long serialVersionUID = 1L;

	boolean auto;
	int nLanes, lw, lh, lsp, lhoff, lvoff;
	int ladderLane;
	int degBG;
	double polyDerivative, tolPK, areaDrift, sdDrift;
	boolean continuum;
	int dist; // index of the fragment distribution
	int dlo, dhi, every;
	List<Peak> customPeaks = new ArrayList<>();
	boolean fitDone;
	boolean showBands; // false when loading files saved before it was added
}
