package gellanesfit;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import org.scijava.Context;
import org.scijava.app.StatusService;
import org.scijava.display.DisplayService;
import org.scijava.log.LogService;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.ProfilePlot;
import ij.gui.Roi;

/**
 * A reference fit of real lanes of the sample image, to catch any change in
 * the profiles or the fitting results.
 */
public class SampleImageFitTest {

	private static final String IMAGE =
		"src/main/resources/sample-images/tagment-test/gel-camera-1/Long_5s.tif";

	/** Lanes as x-centres: a 100 bp-type ladder, one strong band, bands on a
	 * smear, and a smear */
	private static final int[] LANE_CENTRES = { 510, 590, 675, 915 };
	private static final int LANE_WIDTH = 40;
	private static final int TOP = 170;
	private static final int BOTTOM = 990;

	private static Context context;

	@BeforeClass
	public static void createContext() {
		context = new Context(LogService.class, StatusService.class,
			DisplayService.class);
	}

	@AfterClass
	public static void disposeContext() {
		context.dispose();
	}

	/** The lane's profile, built the way Plotter.getLaneProfile() builds it */
	private static DataSeries profile(final ImagePlus imp, final int lane,
		final Roi roi)
	{
		imp.setRoi(roi);
		final double[] y = new ProfilePlot(imp, true).getProfile();
		imp.killRoi();
		final double y0 = roi.getBounds().getMinY();
		final double[] x = new double[y.length];
		for (int p = 0; p < x.length; p++)
			x[p] = y0 + p;
		return new DataSeries("Lane " + lane, lane, DataSeries.PROFILE,
			new ArrayRealVector(x), new ArrayRealVector(y), java.awt.Color.BLACK);
	}

	@Test
	public void fitsTheSampleLanes() {
		final ImagePlus imp = IJ.openImage(IMAGE);
		assertNotNull(IMAGE, imp);

		final List<DataSeries> lanes = new ArrayList<>();
		for (int i = 0; i < LANE_CENTRES.length; i++)
			lanes.add(profile(imp, i + 1, new Roi(LANE_CENTRES[i] - LANE_WIDTH / 2,
				TOP, LANE_WIDTH, BOTTOM - TOP)));

		final List<String> profileSums = new ArrayList<>();
		for (final DataSeries d : lanes)
			profileSums.add(d.getLane() + " " + d.getItemCount() + " " + d.getY()
				.getL1Norm());
		Reference.check("sample-image-profiles", profileSums);

		final Fitter fitter = BandedFitTest.bandedFitter(context, lanes.toArray(
			new DataSeries[0]));
		final List<Integer> numbers = new ArrayList<>();
		for (final DataSeries d : lanes)
			numbers.add(d.getLane());
		fitter.doFit(numbers);

		final List<Peak> all = new ArrayList<>();
		for (final int n : numbers) {
			final List<Peak> peaks = fitter.getFittedPeaks(n);
			assertTrue("lane " + n + " has bands", peaks.size() > 0);
			all.addAll(peaks);
		}
		assertEquals("the ladder's 12 bands", 12, fitter.getFittedPeaks(1)
			.size());
		Reference.check("sample-image-banded", Reference.peaks(all));
	}
}
