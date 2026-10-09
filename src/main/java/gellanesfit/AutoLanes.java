/*
 * Gel Lanes Fit - AutoLanes.java
 * Author: Rick Ziraldo, 2017
 * The University of Texas at Dallas, Richardson, TX
 *
 * Licensed under the GNU Affero General Public License v3.0; see LICENSE.
 * Source: https://github.com/rickud/gel-lanes-fit
 */

package gellanesfit;

/**
 * The geometry of the automatic lanes: their number, size, spacing and
 * offset from the image's top-left corner, in pixels. {@link #keepOnImage}
 * keeps them all on the image.
 */
final class AutoLanes {

	/** Narrowest and shortest lane, as on the sliders, in px */
	static final int MIN_SIZE = 10;

	/** The value just changed, which {@link #keepOnImage} holds at its limit */
	enum Changed {
		NONE, COUNT, WIDTH, HEIGHT, SPACE, H_OFFSET, V_OFFSET
	}

	int count, width, height, space, hOffset, vOffset;

	AutoLanes(final int count, final int width, final int height,
		final int space, final int hOffset, final int vOffset)
	{
		this.count = count;
		this.width = width;
		this.height = height;
		this.space = space;
		this.hOffset = hOffset;
		this.vOffset = vOffset;
	}

	/** From the left edge of the first lane to the right edge of the last */
	int span() {
		final int n = Math.max(1, count);
		return n * width + (n - 1) * space;
	}

	/**
	 * Keeps the lanes on an image of this size. The value just changed is held
	 * at the largest that still fits. Anything still off the image after that,
	 * e.g. with more lanes, makes Width and Space shrink in proportion, then the
	 * offset.
	 */
	void keepOnImage(final int imageWidth, final int imageHeight,
		final Changed changed)
	{
		final int n = Math.max(1, count);
		switch (changed) {
			case WIDTH:
				width = Math.max(MIN_SIZE, Math.min(width, (imageWidth - hOffset -
					(n - 1) * space) / n));
				break;
			case SPACE:
				if (n > 1) space = Math.max(1, Math.min(space, (imageWidth - hOffset -
					n * width) / (n - 1)));
				break;
			case H_OFFSET:
				hOffset = Math.max(0, Math.min(hOffset, imageWidth - span()));
				break;
			case HEIGHT:
				height = Math.max(MIN_SIZE, Math.min(height, imageHeight - vOffset));
				break;
			case V_OFFSET:
				vOffset = Math.max(0, Math.min(vOffset, imageHeight - height));
				break;
			default:
				break;
		}
		if (hOffset + span() > imageWidth) {
			final double f = Math.max(0, imageWidth - hOffset) / (double) span();
			width = Math.max(MIN_SIZE, (int) (width * f));
			space = Math.max(1, (int) (space * f));
			hOffset = Math.max(0, Math.min(hOffset, imageWidth - span()));
		}
		if (vOffset + height > imageHeight) {
			height = Math.max(MIN_SIZE, Math.min(height, imageHeight - vOffset));
			vOffset = Math.max(0, Math.min(vOffset, imageHeight - height));
		}
	}
}
