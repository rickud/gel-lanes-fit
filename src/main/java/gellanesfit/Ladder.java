/*
 * Gel Lanes Fit - Ladder.java
 * Author: Rick Ziraldo, 2017
 * The University of Texas at Dallas, Richardson, TX
 *
 * Licensed under the GNU Affero General Public License v3.0; see LICENSE.
 * Source: https://github.com/rickud/gel-lanes-fit
 */

package gellanesfit;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Serializable;

import javax.swing.JFileChooser;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

import ij.IJ;

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
