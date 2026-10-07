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
	public void uniformListsLengthsFromUpperDownInSteps() {
		final double[][] d = FragmentDistribution.uniform(100, 700, 20);
		assertEquals(31, d.length);
		assertEquals(700, d[0][FragmentDistribution.LENGTH], 0);
		assertEquals(100, d[30][FragmentDistribution.LENGTH], 0);
		assertEquals(Ladder.molecularWeight(400), d[15][FragmentDistribution.MW],
			1e-9);
		// Known issue, kept until it's fixed: 1 / (upper - lower + 1) instead of
		// 1 / 31, so the frequencies don't add up to 1 for steps above 1 bp
		for (final double[] row : d)
			assertEquals(1.0 / 601, row[FragmentDistribution.FREQUENCY], 1e-15);
	}

	@Test
	public void uniformCountMatchesTheDistribution() {
		assertEquals(31, FragmentDistribution.uniformCount(100, 700, 20));
		assertEquals(196, FragmentDistribution.uniformCount(100, 4000, 20));
		assertEquals(FragmentDistribution.uniform(100, 4000, 20).length,
			FragmentDistribution.uniformCount(100, 4000, 20));
	}

	@Test
	public void fitTimeEstimateGrowsSteeply() {
		// The two measurements it's based on
		assertEquals(0.3, FragmentDistribution.estimatedFitSeconds(31), 1e-9);
		assertEquals(25, FragmentDistribution.estimatedFitSeconds(196), 3);
		assertEquals(0, FragmentDistribution.estimatedFitSeconds(0), 0);
	}
}
