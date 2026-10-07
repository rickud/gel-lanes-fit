package gellanesfit;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.awt.Rectangle;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import java.util.Arrays;

import org.junit.Test;

/**
 * saved-state.bak files must stay readable: users keep them next to their
 * images across plugin versions. The fixtures in
 * src/test/resources/gellanesfit/saved-state were written by the plugin's
 * classes as of October 2026; renaming or removing fields of FitState, Ladder
 * or Peak breaks this test, and with it every saved analysis.
 */
public class SavedStateTest {

	private static SavedStateFile fixture(final String name) throws Exception {
		final InputStream in = SavedStateTest.class.getResourceAsStream(
			"saved-state/" + name);
		assertNotNull(name, in);
		try (InputStream i = in) {
			return SavedStateFile.read(i);
		}
	}

	@Test
	public void readsTheCurrentFormat() throws Exception {
		final SavedStateFile file = fixture("current.bak");
		assertEquals(3, file.lanes.size());
		assertEquals(new Rectangle(20, 50, 40, 400), file.lanes.get(0));
		assertEquals(new Rectangle(120, 55, 40, 395), file.lanes.get(2));

		final FitState st = file.state;
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

		assertEquals(Ladder.HILO, file.ladder.getType());
		assertArrayEquals(new int[] { 7, 9 }, file.ladder.getRange());
		assertEquals(3, file.ladder.getMolecularWeights().getDimension());
	}

	@Test
	public void readsFilesFromBeforeFitStateExisted() throws Exception {
		final SavedStateFile file = fixture("before-fitstate.bak");
		assertEquals(2, file.lanes.size());
		assertNull(file.state);
		assertEquals(Ladder.HILO, file.ladder.getType());
	}

	@Test
	public void roundTrip() throws Exception {
		final FitState st = new FitState();
		st.degBG = 3;
		st.customPeaks.add(new Peak(1, 10, 200, 2));
		final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
		new SavedStateFile(Arrays.asList(new Rectangle(1, 2, 3, 4)), st, null)
			.write(bytes);

		final SavedStateFile copy = SavedStateFile.read(new ByteArrayInputStream(
			bytes.toByteArray()));
		assertEquals(new Rectangle(1, 2, 3, 4), copy.lanes.get(0));
		assertEquals(3, copy.state.degBG);
		assertEquals(200, copy.state.customPeaks.get(0).getMean(), 0);
		assertEquals(false, copy.state.showBands);
		assertNull(copy.ladder);
	}
}
