/*
 * Gel Lanes Fit - SavedStateFile.java
 * Author: Rick Ziraldo, 2017
 * The University of Texas at Dallas, Richardson, TX
 *
 * Licensed under the GNU Affero General Public License v3.0; see LICENSE.
 * Source: https://github.com/rickud/gel-lanes-fit
 */

package gellanesfit;

import java.awt.Rectangle;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;

/**
 * The contents of an image's saved-state.bak: the manual lanes, the fit
 * settings and the ladder, written with Java serialization in that order.
 * <p>
 * Files from before FitState existed have only the lanes and the ladder.
 * Users keep these files across plugin versions, so the classes they contain
 * must stay compatible (see SavedStateTest).
 * </p>
 */
final class SavedStateFile {

	static final String NAME = "saved-state.bak";

	/** The manual lanes, left to right */
	final List<Rectangle> lanes;
	/** The fit settings; null in files from before they were saved */
	final FitState state;
	/** The ladder; null if none was chosen */
	final Ladder ladder;

	SavedStateFile(final List<Rectangle> lanes, final FitState state,
		final Ladder ladder)
	{
		this.lanes = lanes;
		this.state = state;
		this.ladder = ladder;
	}

	/**
	 * Reads a saved state. Reading stops at the end of the file, or at the
	 * first object that can't be read, keeping what was read until then.
	 */
	static SavedStateFile read(final InputStream in) throws IOException {
		final List<Rectangle> lanes = new ArrayList<>();
		FitState state = null;
		Ladder ladder = null;
		final ObjectInputStream ois = new ObjectInputStream(in);
		try {
			while (true) {
				final Object o = ois.readObject();
				if (o instanceof Rectangle) lanes.add((Rectangle) o);
				else if (o instanceof Ladder) ladder = (Ladder) o;
				else if (o instanceof FitState) state = (FitState) o;
			}
		}
		catch (final Exception e) {
			// End of the file, or an object from an incompatible version
		}
		return new SavedStateFile(lanes, state, ladder);
	}

	/** Writes the lanes, then the settings, then the ladder if there is one */
	void write(final OutputStream out) throws IOException {
		final ObjectOutputStream oos = new ObjectOutputStream(out);
		for (final Rectangle r : lanes)
			oos.writeObject(r);
		oos.writeObject(state);
		if (ladder != null) oos.writeObject(ladder);
		oos.flush();
	}
}
