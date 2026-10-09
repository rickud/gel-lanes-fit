/*
 * Gel Lanes Fit - FitState.java
 * Author: Rick Ziraldo, 2017
 * The University of Texas at Dallas, Richardson, TX
 *
 * Licensed under the GNU Affero General Public License v3.0; see LICENSE.
 * Source: https://github.com/rickud/gel-lanes-fit
 */

package gellanesfit;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

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
	/**
	 * Uniform distribution: set by the spacing {@code every} in bp, or else by
	 * the number of lengths {@code uniformCount} (0 in files saved before it
	 * was added: use the default)
	 */
	boolean uniformBySpacing;
	int every, uniformCount;
	List<Peak> customPeaks = new ArrayList<>();
	boolean fitDone;
	/**
	 * A fit had started but not finished when the state was saved; false in
	 * files saved before it was added
	 */
	boolean fitRunning;
	boolean showBands; // false when loading files saved before it was added
}
