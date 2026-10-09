package gellanesfit;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;

import java.io.BufferedReader;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

public class FragmentDistributionTest {

	private static double[][] read(final String text, final List<String> warnings)
		throws Exception
	{
		return FragmentDistribution.read(new BufferedReader(new StringReader(
			text)), warnings::add);
	}

	@Test
	public void readsLengthsAndNormalizesFrequencies() throws Exception {
		final List<String> warnings = new ArrayList<>();
		final double[][] d = read("Fragment Length\t# of Fragments\tWeight\n" +
			"300\t1\t1.5\n200\t3\n100\t4\n", warnings);
		assertEquals("the header isn't a warning", 0, warnings.size());
		assertEquals(3, d.length);
		assertArrayEquals(new double[] { 0.125, 300, Ladder.molecularWeight(300) },
			d[0], 1e-12);
		assertEquals(0.375, d[1][FragmentDistribution.FREQUENCY], 1e-12);
		assertEquals(100, d[2][FragmentDistribution.LENGTH], 0);
	}

	@Test
	public void reportsAndSkipsLinesThatAreNotFragments() throws Exception {
		final List<String> warnings = new ArrayList<>();
		final double[][] d = read("length\tcount\n300\t1\nnote\n200\tx\n100\t1\n",
			warnings);
		assertEquals(2, d.length);
		assertEquals(2, warnings.size());
		assertEquals("Invalid line 3: note", warnings.get(0));
	}

	@Test
	public void uniformSpacesLengthsEvenlyFromUpperToLower() {
		final double[][] d = FragmentDistribution.uniform(100, 700, 31);
		assertEquals(31, d.length);
		assertEquals(700, d[0][FragmentDistribution.LENGTH], 0);
		assertEquals(680, d[1][FragmentDistribution.LENGTH], 1e-9);
		assertEquals(100, d[30][FragmentDistribution.LENGTH], 1e-9);
		assertEquals(Ladder.molecularWeight(400), d[15][FragmentDistribution.MW],
			1e-9);
		double sum = 0;
		for (final double[] row : d)
			sum += row[FragmentDistribution.FREQUENCY];
		assertEquals("frequencies add up to 1", 1, sum, 1e-12);
	}

	@Test
	public void uniformCountKeepsLengthsAtMostTheSpacingApart() {
		assertEquals(31, FragmentDistribution.uniformCount(100, 700, 20));
		assertEquals(196, FragmentDistribution.uniformCount(100, 4000, 20));
		// 610 bp doesn't divide by 20: 32 lengths, 19.7 bp apart
		assertEquals(32, FragmentDistribution.uniformCount(100, 710, 20));
		assertEquals("at least both ends", 2, FragmentDistribution.uniformCount(
			100, 110, 50));
	}

	@Test
	public void fitTimeEstimateGrowsSteeply() {
		// The two measurements it's based on
		assertEquals(0.3, FragmentDistribution.estimatedFitSeconds(31), 1e-9);
		assertEquals(25, FragmentDistribution.estimatedFitSeconds(196), 3);
		assertEquals(0, FragmentDistribution.estimatedFitSeconds(0), 0);
	}
}
