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
	 * sample image: 31 fragments took about 0.3 s per lane and 196 about 25 s,
	 * so the time grows roughly with the 2.4th power of the number. Computers
	 * differ, so it's only an indication.
	 */
	static double estimatedFitSeconds(final int fragments) {
		return 0.3 * Math.pow(fragments / 31.0, 2.4);
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

	/** How many lengths {@link #uniform} gives for these values */
	static int uniformCount(final int lower, final int upper, final int every) {
		return (upper - lower) / every + 1;
	}

	/**
	 * A uniform distribution: every length from {@code upper} down to
	 * {@code lower}, in steps of {@code every} bp, all equally frequent.
	 * <p>
	 * Known issue, kept for now so that results don't change: the frequency is
	 * 1 / (upper - lower + 1), which is right only for a 1 bp step, so the
	 * frequencies don't add up to 1 otherwise. Fits aren't affected, since the
	 * guess rescales them; the Frequency column of the results is.
	 * </p>
	 */
	static double[][] uniform(final int lower, final int upper,
		final int every)
	{
		final double[][] dist = new double[uniformCount(lower, upper, every)][3];
		final double f = 1.0 / (upper - lower + 1);
		for (int i = 0; i < dist.length; i++) {
			dist[i][FREQUENCY] = f;
			dist[i][LENGTH] = upper - i * every;
			dist[i][MW] = Ladder.molecularWeight(dist[i][LENGTH]);
		}
		return dist;
	}
}
