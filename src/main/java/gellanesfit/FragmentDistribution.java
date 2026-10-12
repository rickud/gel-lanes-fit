/*
 * Gel Lanes Fit - FragmentDistribution.java
 * Author: Rick Ziraldo, 2017
 * The University of Texas at Dallas, Richardson, TX
 *
 * Licensed under the GNU Affero General Public License v3.0; see LICENSE.
 * Source: https://github.com/rickud/gel-lanes-fit
 */

package gellanesfit;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

/**
 * Fragment distributions for Continuum fits: the fragment lengths expected in
 * a sample and how common each one is.
 * <p>
 * A distribution is an array with one row per fragment length and three
 * columns: the relative frequency ({@link #FREQUENCY}), the length in base
 * pairs ({@link #LENGTH}) and the molecular weight in Da ({@link #MW}).
 * </p>
 */
final class FragmentDistribution {

	static final int FREQUENCY = 0;
	static final int LENGTH = 1;
	static final int MW = 2;

	/** Where the bundled distributions are, on the classpath */
	static final String FOLDER = "sample-distributions/";

	private FragmentDistribution() {}

	/**
	 * A rough estimate of how long a Continuum fit of one lane takes, in
	 * seconds, from the number of fragment peaks it has. Measured on the
	 * sample image's smears: about 0.5 s up to 80 fragments, 2 s for 150, 4 s
	 * for 200, 17 s for 300 and 45 to 90 s for 400; beyond about 100 the time
	 * grows with roughly the 3.6th power of the number. This is for the
	 * computer it was measured on; {@link FitTimer} adapts it to the one the
	 * plugin runs on.
	 */
	static double estimatedFitSeconds(final int fragments) {
		if (fragments <= 0) return 0;
		return 0.4 + 3.6 * Math.pow(fragments / 200.0, 3.6);
	}

	/**
	 * Reads a distribution file: a header line, then one line per fragment
	 * length with the length in bp and its number of copies, separated by a
	 * tab. Further columns are ignored. The frequencies are normalized to add up
	 * to 1.
	 *
	 * @param reader the file's contents
	 * @param warnings receives a message for each line that isn't a fragment,
	 *          apart from the header
	 */
	static double[][] read(final BufferedReader reader,
		final Consumer<String> warnings) throws IOException
	{
		final List<int[]> rows = new ArrayList<>();
		int lineNumber = 0;
		String line;
		while ((line = reader.readLine()) != null) {
			lineNumber++;
			final String[] words = line.split("\t");
			try {
				if (words.length < 2) throw new NumberFormatException();
				rows.add(new int[] { Integer.parseInt(words[0].trim()), Integer
					.parseInt(words[1].trim()) });
			}
			catch (final NumberFormatException e) {
				if (lineNumber > 1) warnings.accept("Invalid line " + lineNumber +
					": " + line);
			}
		}
		double count = 0;
		for (final int[] row : rows)
			count += row[1];
		final double[][] out = new double[rows.size()][3];
		for (int i = 0; i < out.length; i++) {
			out[i][FREQUENCY] = rows.get(i)[1] / count;
			out[i][LENGTH] = rows.get(i)[0];
			out[i][MW] = Ladder.molecularWeight(rows.get(i)[0]);
		}
		return out;
	}

	/** The default number of lengths of a Uniform distribution */
	static final int UNIFORM_COUNT = 30;
	/** The default number of lengths of a Log Spaced distribution */
	static final int LOG_COUNT = 150;

	/**
	 * How many lengths a Uniform distribution over this range needs so that
	 * they're at most {@code spacing} bp apart; at least 2.
	 */
	static int uniformCount(final double lower, final double upper,
		final double spacing)
	{
		return Math.max(2, (int) Math.ceil((upper - lower) / spacing - 1e-9) + 1);
	}

	/**
	 * A uniform distribution: {@code count} lengths evenly spaced from
	 * {@code upper} down to {@code lower} bp, both included, all equally
	 * frequent.
	 */
	static double[][] uniform(final double lower, final double upper,
		final int count)
	{
		return spaced(lower, upper, count, false);
	}

	/**
	 * A log-spaced distribution: {@code count} lengths from {@code upper} down
	 * to {@code lower} bp, both included, evenly spaced in log(bp), so that
	 * small fragments get more of them than large ones; all equally frequent.
	 * Since fragments run roughly by the log of their size, they end up
	 * roughly evenly spaced along the lane.
	 */
	static double[][] logSpaced(final double lower, final double upper,
		final int count)
	{
		return spaced(lower, upper, count, true);
	}

	private static double[][] spaced(final double lower, final double upper,
		final int count, final boolean log)
	{
		final int n = Math.max(2, count);
		final double[][] dist = new double[n][3];
		for (int i = 0; i < n; i++) {
			dist[i][FREQUENCY] = 1.0 / n;
			dist[i][LENGTH] = log ? upper * Math.pow(lower / upper, i / (n -
				1.0)) : upper - i * (upper - lower) / (n - 1);
			dist[i][MW] = Ladder.molecularWeight(dist[i][LENGTH]);
		}
		return dist;
	}
}
