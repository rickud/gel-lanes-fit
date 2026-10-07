package gellanesfit;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

import org.apache.commons.math3.linear.RealVector;
import org.junit.Test;

public class LadderTest {

	private static final double[] HILO_BP = { 10000, 8000, 6000, 4000, 3000,
		2000, 1550, 1400, 1000, 750, 500, 400, 300, 200, 100, 50 };

	private static double mw(final double bp) {
		return bp * 607.4 + 157.9;
	}

	@Test
	public void builtInLaddersHaveOneNamePerBand() {
		final int[][] typeAndBands = { { Ladder.HILO, 16 }, { Ladder.BP100, 12 },
			{ Ladder.QUICKLOAD, 13 }, { Ladder.TAPESTATION, 10 } };
		for (final int[] t : typeAndBands) {
			final Ladder ladder = Ladder.create(t[0]);
			assertNotNull(ladder);
			assertEquals(t[0], ladder.getType());
			assertEquals(t[1], ladder.getStrings().length);
			assertArrayEquals(new int[] { 0, t[1] - 1 }, ladder.getRange());
			assertEquals(t[1], ladder.getMolecularWeights().getDimension());
		}
	}

	@Test
	public void hiLoNamesAndWeights() {
		final Ladder ladder = Ladder.create(Ladder.HILO);
		assertEquals("10 kbp", ladder.getStrings()[0]);
		assertEquals("1.55 kbp", ladder.getStrings()[6]);
		assertEquals("50 bp", ladder.getStrings()[15]);
		final RealVector mw = ladder.getMolecularWeights();
		for (int i = 0; i < HILO_BP.length; i++)
			assertEquals(mw(HILO_BP[i]), mw.getEntry(i), 1e-9);
	}

	@Test
	public void rangeSelectsTheBandsInIt() {
		final Ladder ladder = Ladder.create(Ladder.HILO);
		ladder.setRange(new int[] { 7, 9 }); // 1.4 kbp to 750 bp
		final RealVector mw = ladder.getMolecularWeights();
		assertEquals(3, mw.getDimension());
		assertEquals(mw(1400), mw.getEntry(0), 1e-9);
		assertEquals(mw(750), mw.getEntry(2), 1e-9);
	}

	@Test
	public void otherLaddersFirstAndLastBands() {
		final Object[][] cases = { { Ladder.BP100, 1517.0, 100.0, "1.5 kbp",
			"100 bp" }, { Ladder.QUICKLOAD, 48500.0, 500.0, "48.5 kbp", "500 bp" },
			{ Ladder.TAPESTATION, 1500.0, 25.0, "1.5 kbp", "25 bp" } };
		for (final Object[] c : cases) {
			final Ladder ladder = Ladder.create((Integer) c[0]);
			final RealVector mw = ladder.getMolecularWeights();
			assertEquals(mw((Double) c[1]), mw.getEntry(0), 1e-9);
			assertEquals(mw((Double) c[2]), mw.getEntry(mw.getDimension() - 1),
				1e-9);
			assertEquals(c[3], ladder.getStrings()[0]);
			assertEquals(c[4], ladder.getStrings()[ladder.getStrings().length - 1]);
		}
	}

	@Test
	public void survivesSerialization() throws Exception {
		final Ladder ladder = Ladder.create(Ladder.HILO);
		ladder.setRange(new int[] { 7, 9 });
		final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
		try (ObjectOutputStream out = new ObjectOutputStream(bytes)) {
			out.writeObject(ladder);
		}
		try (ObjectInputStream in = new ObjectInputStream(
			new ByteArrayInputStream(bytes.toByteArray())))
		{
			final Ladder copy = (Ladder) in.readObject();
			assertEquals(Ladder.HILO, copy.getType());
			assertArrayEquals(new int[] { 7, 9 }, copy.getRange());
			assertEquals(ladder.getMolecularWeights(), copy.getMolecularWeights());
		}
	}
}
