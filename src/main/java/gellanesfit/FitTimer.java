/*
 * Gel Lanes Fit - FitTimer.java
 * Author: Rick Ziraldo, 2017
 * The University of Texas at Dallas, Richardson, TX
 *
 * Licensed under the GNU Affero General Public License v3.0; see LICENSE.
 * Source: https://github.com/rickud/gel-lanes-fit
 */

package gellanesfit;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.fitting.WeightedObservedPoint;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;

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
	static final double REFERENCE_BENCHMARK_SECONDS = 0.195;

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
	 * Times a fixed fit: 31 Gaussian bands on a flat background, from a
	 * slightly offset start. One run warms up the code, then runs are repeated
	 * for at least a quarter of a second.
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
		final int bands = 31;
		final List<WeightedObservedPoint> obs = new ArrayList<>();
		for (int x = 0; x <= 400; x++) {
			double y = 10;
			for (int b = 0; b < bands; b++)
				y += new Gaussian(50 + b % 5 * 10, 20 + b * 12, 4).value(x);
			obs.add(new WeightedObservedPoint(1, x, y));
		}
		RealVector norm = new ArrayRealVector();
		RealVector mean = new ArrayRealVector();
		RealVector sd = new ArrayRealVector();
		for (int b = 0; b < bands; b++) {
			norm = norm.append(0.9 * (50 + b % 5 * 10));
			mean = mean.append(20 + b * 12 + 1.5);
			sd = sd.append(4.5);
		}
		final SortedParameters start = new SortedParameters(new ArrayRealVector(
			new double[] { 9, 0, 0 }), norm, mean, sd);
		final GaussianArrayCurveFitter cf = GaussianArrayCurveFitter.create(
			Fitter.bandMode, 2, 10, 0.98, 5, 0.1, 1.0).withStartPoint(start);
		cf.getOptimizer().optimize(cf.getProblem(obs));
	}
}
