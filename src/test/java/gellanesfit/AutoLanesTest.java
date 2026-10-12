package gellanesfit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

/** The automatic lanes stay on the image whatever the sliders ask for. */
public class AutoLanesTest {

	/** The sample image's size */
	private static final int IW = 1392;
	private static final int IH = 1032;

	private static void assertOnImage(final AutoLanes l) {
		assertTrue("left edge", l.hOffset >= 0);
		assertTrue("right edge: " + (l.hOffset + l.span()), l.hOffset + l
			.span() <= IW);
		assertTrue("top edge", l.vOffset >= 0);
		assertTrue("bottom edge", l.vOffset + l.height <= IH);
	}

	@Test
	public void lanesThatFitAreLeftAlone() {
		final AutoLanes l = new AutoLanes(12, 74, 570, 30, 100, 318);
		l.keepOnImage(IW, IH, AutoLanes.Changed.NONE);
		assertEquals(74, l.width);
		assertEquals(30, l.space);
		assertEquals(100, l.hOffset);
	}

	@Test
	public void restoresASavedSpaceTooLargeForTheImage() {
		// A saved state where Space had been dragged to the image width:
		// lanes 2 to 12 were off the image
		final AutoLanes l = new AutoLanes(12, 74, 570, 1392, 387, 318);
		l.keepOnImage(IW, IH, AutoLanes.Changed.NONE);
		assertOnImage(l);
		assertEquals("all 12 lanes kept", 12, l.count);
	}

	@Test
	public void spaceStopsAtTheLargestThatFits() {
		final AutoLanes l = new AutoLanes(12, 74, 570, 1392, 100, 318);
		l.keepOnImage(IW, IH, AutoLanes.Changed.SPACE);
		assertEquals("width unchanged", 74, l.width);
		assertEquals((IW - 100 - 12 * 74) / 11, l.space);
		assertOnImage(l);
	}

	@Test
	public void widthAndOffsetStopAtTheLargestThatFits() {
		AutoLanes l = new AutoLanes(4, 1000, 570, 20, 100, 318);
		l.keepOnImage(IW, IH, AutoLanes.Changed.WIDTH);
		assertEquals((IW - 100 - 3 * 20) / 4, l.width);
		assertOnImage(l);

		l = new AutoLanes(4, 100, 570, 20, 1300, 318);
		l.keepOnImage(IW, IH, AutoLanes.Changed.H_OFFSET);
		assertEquals(IW - (4 * 100 + 3 * 20), l.hOffset);
		assertEquals("width unchanged", 100, l.width);
	}

	@Test
	public void heightAndVerticalOffsetStopAtTheBottom() {
		AutoLanes l = new AutoLanes(4, 100, 2000, 20, 0, 100);
		l.keepOnImage(IW, IH, AutoLanes.Changed.HEIGHT);
		assertEquals(IH - 100, l.height);

		l = new AutoLanes(4, 100, 500, 20, 0, 900);
		l.keepOnImage(IW, IH, AutoLanes.Changed.V_OFFSET);
		assertEquals(IH - 500, l.vOffset);
		assertEquals("height unchanged", 500, l.height);
	}

	@Test
	public void moreLanesShrinkWidthAndSpaceInProportion() {
		// 4 lanes 200 px wide, 100 apart, then 8 lanes
		final AutoLanes l = new AutoLanes(8, 200, 570, 100, 50, 318);
		l.keepOnImage(IW, IH, AutoLanes.Changed.COUNT);
		assertOnImage(l);
		assertEquals("offset kept", 50, l.hOffset);
		assertEquals("same proportions", 2.0, l.width / (double) l.space, 0.1);
	}

	@Test
	public void manyLanesAtTheMinimumSizeMoveTheOffset() {
		final AutoLanes l = new AutoLanes(100, 50, 570, 10, 800, 318);
		l.keepOnImage(IW, IH, AutoLanes.Changed.COUNT);
		assertOnImage(l);
		assertEquals(AutoLanes.MIN_SIZE, l.width);
	}
}
