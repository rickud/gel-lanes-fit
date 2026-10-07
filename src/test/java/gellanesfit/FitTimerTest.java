package gellanesfit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class FitTimerTest {

	@Test
	public void formatsDurations() {
		assertEquals("1 s", FitTimer.format(0.2));
		assertEquals("45 s", FitTimer.format(45.2));
		assertEquals("89 s", FitTimer.format(89));
		assertEquals("2 min", FitTimer.format(95));
		assertEquals("12 min", FitTimer.format(12 * 60 + 10));
	}

	@Test
	public void calibratesToThisComputer() {
		final double factor = FitTimer.machineFactor();
		System.out.printf("This computer's speed factor: %.2f%n", factor);
		assertTrue(factor > 0 && !Double.isInfinite(factor));
		assertEquals("measured once per session", factor, FitTimer.machineFactor(),
			0);
		assertEquals(FragmentDistribution.estimatedFitSeconds(150) * factor,
			FitTimer.estimatedSeconds(150), 1e-9);
	}

	@Test
	public void progressMessageEstimatesTheTimeLeft() {
		assertEquals("Fitted 2 of 10 lanes, about 2 min left", Fitter
			.progressMessage(2, 10, 30));
		assertEquals("Fitted 3 of 4 lanes, about 5 s left", Fitter.progressMessage(
			3, 4, 15));
		assertEquals("Fitted 4 of 4 lanes", Fitter.progressMessage(4, 4, 20));
	}
}
