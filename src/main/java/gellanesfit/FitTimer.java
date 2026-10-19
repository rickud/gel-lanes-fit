/*
 * Gel Lanes Fit - FitTimer.java
 * Author: Rick Ziraldo, 2017
 * The University of Texas at Dallas, Richardson, TX
 *
 * Licensed under the GNU Affero General Public License v3.0; see LICENSE.
 * Source: https://github.com/rickud/gel-lanes-fit
 */

package gellanesfit;

import org.apache.commons.math3.analysis.function.Gaussian;

/**
 * Estimates how long a Continuum fit takes on this computer.
 * <p>
 * The estimate per lane, {@link FragmentDistribution#estimatedFitSeconds}, was
 * measured on one computer. To adapt it to the computer the plugin runs on, a
 * small fixed benchmark fit is timed once per session and compared with its
 * time on that computer: the ratio scales every estimate.
 * </p>
 */
final class FitTimer {

	/**
	 * The benchmark's time on the computer the estimate was measured on (an
	 * Apple M-series MacBook, Java 21), in seconds
	 */
	static final double REFERENCE_BENCHMARK_SECONDS = 0.045;

	/** How much slower this computer is than the reference; NaN until measured */
	private static double machineFactor = Double.NaN;

	private FitTimer() {}

	/** A duration for messages: seconds up to a minute and a half, then minutes */
	static String format(final double seconds) {
		if (seconds < 90) return Math.max(1, Math.round(seconds)) + " s";
		return Math.round(seconds / 60) + " min";
	}

	/**
	 * The estimated time of a Continuum fit of one lane with this many fragment
	 * peaks, on this computer, in seconds. The first call times the benchmark
	 * (well under a second).
	 */
	static double estimatedSeconds(final int fragments) {
		return FragmentDistribution.estimatedFitSeconds(fragments) *
			machineFactor();
	}

	/** This computer's speed relative to the reference, measured once */
	static synchronized double machineFactor() {
		if (Double.isNaN(machineFactor)) {
			machineFactor = benchmarkSeconds() / REFERENCE_BENCHMARK_SECONDS;
		}
		return machineFactor;
	}

	/**
	 * Times a fixed fit: a Continuum fit of a synthetic smear of 50 fragment
	 * peaks. One run warms up the code, then runs are repeated for at least a
	 * quarter of a second.
	 *
	 * @return the average time of one run, in seconds
	 */
	static double benchmarkSeconds() {
		runBenchmark();
		final long start = System.nanoTime();
		int runs = 0;
		do {
			runBenchmark();
			runs++;
		}
		while (System.nanoTime() - start < 250_000_000L);
		return (System.nanoTime() - start) / 1e9 / runs;
	}

	private static void runBenchmark() {
		final int n = 50;
		final double[] x = new double[201];
		final double[] y = new double[x.length];
		final double[] means = new double[n];
		final double[] sds = new double[n];
		final double[] expected = new double[n];
		for (int i = 0; i < n; i++) {
			means[i] = 4 * i;
			sds[i] = 4;
			expected[i] = 1;
		}
		for (int j = 0; j < x.length; j++) {
			x[j] = j;
			for (int i = 0; i < n; i++)
				y[j] += new Gaussian(500 * (1.2 + Math.sin(i / 15.0)), means[i], 4)
					.value(j);
		}
		new ContinuumFit(x, y, means, sds, expected, 1, 10, () -> false).fit();
	}
}
