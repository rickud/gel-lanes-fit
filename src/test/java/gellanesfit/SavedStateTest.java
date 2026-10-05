package gellanesfit;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.awt.Rectangle;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.EOFException;
import java.io.InputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

/**
 * saved-state.bak files must stay readable: users keep them next to their
 * images across plugin versions. The fixtures in
 * src/test/resources/gellanesfit/saved-state were written by the plugin's
 * classes as of October 2026; renaming or removing fields of FitState, Ladder
 * or Peak breaks this test, and with it every saved analysis.
 */
public class SavedStateTest {

	/** Reads every object in a saved-state file, as MainDialog.loadState() */
	private static List<Object> read(final InputStream in) throws Exception {
		final List<Object> objects = new ArrayList<>();
		try (ObjectInputStream ois = new ObjectInputStream(in)) {
			while (true)
				objects.add(ois.readObject());
		}
		catch (final EOFException e) {
			// end of file
		}
		return objects;
	}

	private static List<Object> fixture(final String name) throws Exception {
		final InputStream in = SavedStateTest.class.getResourceAsStream(
			"saved-state/" + name);
		assertNotNull(name, in);
		return read(in);
	}

	@Test
	public void readsTheCurrentFormat() throws Exception {
		final List<Object> objects = fixture("current.bak");
		assertEquals(5, objects.size());
		assertEquals(new Rectangle(20, 50, 40, 400), objects.get(0));
		assertEquals(new Rectangle(120, 55, 40, 395), objects.get(2));

		final FitState st = (FitState) objects.get(3);
		assertEquals(false, st.auto);
		assertEquals(3, st.nLanes);
		assertEquals(1, st.ladderLane);
		assertEquals(5, st.polyDerivative, 0);
		assertEquals(0.05, st.tolPK, 0);
		assertEquals(1.5, st.sdDrift, 0);
		assertTrue(st.continuum);
		assertEquals(2, st.dist);
		assertEquals(5, st.every);
		assertTrue(st.fitDone);
		assertTrue(st.showBands);
		assertEquals(2, st.customPeaks.size());
		final Peak p = st.customPeaks.get(0);
		assertEquals(2, p.getLane());
		assertEquals(250.25, p.getMean(), 0);
		assertEquals(120.5, p.getNorm(), 0);
		assertEquals(3.5, p.getSigma(), 0);

		final Ladder ladder = (Ladder) objects.get(4);
		assertEquals(Ladder.HILO, ladder.getType());
		assertArrayEquals(new int[] { 7, 9 }, ladder.getRange());
		assertEquals(3, ladder.getMolecularWeights().getDimension());
	}

	@Test
	public void readsFilesFromBeforeFitStateExisted() throws Exception {
		final List<Object> objects = fixture("before-fitstate.bak");
		assertEquals(3, objects.size());
		assertTrue(objects.get(0) instanceof Rectangle);
		assertEquals(Ladder.HILO, ((Ladder) objects.get(2)).getType());
	}

	@Test
	public void fitStateRoundTrip() throws Exception {
		final FitState st = new FitState();
		st.degBG = 3;
		st.customPeaks.add(new Peak(1, 10, 200, 2));
		final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
		try (ObjectOutputStream out = new ObjectOutputStream(bytes)) {
			out.writeObject(st);
		}
		final FitState copy = (FitState) read(new ByteArrayInputStream(bytes
			.toByteArray())).get(0);
		assertEquals(3, copy.degBG);
		assertEquals(200, copy.customPeaks.get(0).getMean(), 0);
		assertEquals(false, copy.showBands);
	}
}
