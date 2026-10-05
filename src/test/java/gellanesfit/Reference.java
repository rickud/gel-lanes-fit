package gellanesfit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.fail;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * Compares results with a reference file stored with the tests, so that any
 * change in the fitting results shows up.
 * <p>
 * The files are in src/test/resources/gellanesfit/reference. When results
 * change on purpose, run the tests once with {@code -Dglf.record=true} from
 * the project folder to rewrite them, and review the change in git.
 * </p>
 */
final class Reference {

	/** Allowed relative difference, for floating-point noise only */
	private static final double TOLERANCE = 1e-6;

	private static final String DIR = "src/test/resources/gellanesfit/reference/";

	private Reference() {}

	/** One line per fitted peak: lane, position, height and width */
	static List<String> peaks(final List<Peak> peaks) {
		final List<String> lines = new ArrayList<>();
		for (final Peak p : peaks)
			lines.add(String.format(Locale.ROOT, "%d %.9g %.9g %.9g", p.getLane(), p
				.getMean(), p.getNorm(), p.getSigma()));
		return lines;
	}

	static void check(final String name, final List<String> actual) {
		if (Boolean.getBoolean("glf.record")) {
			record(name, actual);
			return;
		}
		final List<String> expected = read(name);
		assertEquals(name + ": number of lines", expected.size(), actual.size());
		for (int l = 0; l < expected.size(); l++) {
			final String[] e = expected.get(l).trim().split("\\s+");
			final String[] a = actual.get(l).trim().split("\\s+");
			assertEquals(name + " line " + (l + 1) + ": number of values", e.length,
				a.length);
			for (int v = 0; v < e.length; v++) {
				final double ev = Double.parseDouble(e[v]);
				final double av = Double.parseDouble(a[v]);
				if (Math.abs(av - ev) > TOLERANCE * Math.max(1, Math.abs(ev))) fail(
					name + " line " + (l + 1) + " value " + (v + 1) + ": expected " +
						ev + " but was " + av);
			}
		}
	}

	private static List<String> read(final String name) {
		final InputStream in = Reference.class.getResourceAsStream("reference/" +
			name + ".txt");
		assertNotNull("No reference file for " + name +
			"; run once with -Dglf.record=true to create it", in);
		final List<String> lines = new ArrayList<>();
		try (BufferedReader r = new BufferedReader(new InputStreamReader(in,
			StandardCharsets.UTF_8)))
		{
			String line;
			while ((line = r.readLine()) != null)
				if (!line.trim().isEmpty() && !line.startsWith("#")) lines.add(line);
		}
		catch (final IOException e) {
			throw new AssertionError(e);
		}
		return lines;
	}

	private static void record(final String name, final List<String> lines) {
		final File file = new File(DIR + name + ".txt");
		file.getParentFile().mkdirs();
		try {
			Files.write(file.toPath(), lines, StandardCharsets.UTF_8);
		}
		catch (final IOException e) {
			throw new AssertionError(e);
		}
		System.out.println("Recorded " + file);
	}
}
