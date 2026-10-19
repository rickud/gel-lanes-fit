/*
 * Gel Lanes Fit - GaussianArray.java
 * Author: Rick Ziraldo, 2017
 * The University of Texas at Dallas, Richardson, TX
 *
 * Licensed under the GNU Affero General Public License v3.0; see LICENSE.
 * Source: https://github.com/rickud/gel-lanes-fit
 */

package gellanesfit;

import org.apache.commons.math3.analysis.differentiation.DerivativeStructure;
import org.apache.commons.math3.analysis.differentiation.UnivariateDifferentiableFunction;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.util.FastMath;

/**
 * The fit's parameters, split up: the background polynomial's degree and
 * coefficients (constant first), and each peak's height, position and standard
 * deviation. {@link #getParameters()} gives them as one array: degree,
 * coefficients, then a triplet per peak.
 */
class SortedParameters {
	private double[] parameters;
	private RealVector mean;
	private RealVector norm;
	private RealVector sd;
	private RealVector area;
	private final int deg;
	private final RealVector poly; // Polynomial background curve
	
	public SortedParameters(double[] parameters) throws DimensionMismatchException
	{
		this.parameters = parameters;
		this.deg = (int) parameters[0];
		RealVector par = new ArrayRealVector(parameters);
		
		this.mean = new ArrayRealVector();
		this.norm = new ArrayRealVector();
		this.sd   = new ArrayRealVector();
		this.poly = par.getSubVector(1, deg + 1);
		this.area = new ArrayRealVector();
		
		for (int i = deg + 2; i < par.getDimension(); i += 3) {
			norm = norm.append(par.getEntry(i));
			mean = mean.append(par.getEntry(i + 1));
			sd = sd.append(par.getEntry(i + 2));
		}
		area = norm.ebeMultiply(sd).mapMultiply(FastMath.sqrt(2 *
			FastMath.PI));
	}
	
	public SortedParameters(RealVector poly, RealVector norm,
			RealVector mean, RealVector sd)
	{
		if (norm.getDimension() != mean.getDimension())
			throw new DimensionMismatchException(norm.getDimension(), mean.getDimension());
		if (mean.getDimension() != sd.getDimension())
			throw new DimensionMismatchException(mean.getDimension(), sd.getDimension());
		
		this.norm = norm;
		this.mean = mean;
		this.sd   = sd;
		this.poly = poly;
		deg = poly.getDimension() - 1;
		RealVector parArray = new ArrayRealVector().append(deg).append(poly);
		for (int i = 0; i< norm.getDimension(); i++) {
			parArray = parArray.append(norm.getEntry(i));
			parArray = parArray.append(mean.getEntry(i));
			parArray = parArray.append(sd.getEntry(i));
		}
		this.parameters = parArray.toArray();
	}
	
	public int getDeg() {
		return deg;
	}
	public RealVector getMean() {
		return mean;
	}
	public RealVector getNorm() {
		return norm;
	}
	public RealVector getSD() {
		return sd;
	}
	public RealVector getPoly() {
		return poly;
	}
	public RealVector getArea() {
		return area;
	}
	public double[] getParameters() {
		return parameters;
	}
}

/**
 * The fitted curve: the sum of the Gaussian peaks and the polynomial
 * background. An empty background (degree -1) counts as zero.
 */
class GaussianArray implements UnivariateDifferentiableFunction {
	// Implements a train of Gaussian peaks
	// with an optional polynomial background function

	private final SortedParameters sp;

	public GaussianArray(final SortedParameters sp)
		throws NotStrictlyPositiveException {
		this.sp = sp;
	}

	@Override
	public double value(final double x) {
		double output = 0;
		for (int i = 0; i < sp.getNorm().getDimension(); i++) {
			output += new Gaussian( sp.getNorm().getEntry(i),
															sp.getMean().getEntry(i), 
															sp.getSD().getEntry(i)).value(x);
		}
		if (sp.getPoly().getDimension() > 0) output += new PolynomialFunction(sp
			.getPoly().toArray()).value(x);
		return output;
	}

	@Override
	public DerivativeStructure value(final DerivativeStructure t)
		throws DimensionMismatchException
	{
		DerivativeStructure output = t.getField().getZero();
		for (int i = 0; i < sp.getNorm().getDimension(); i++) {
			output = output.add(new Gaussian(sp.getNorm().getEntry(i), sp.getMean()
				.getEntry(i), sp.getSD().getEntry(i)).value(t));
		}
		if (sp.getPoly().getDimension() > 0) output = output.add(
			new PolynomialFunction(sp.getPoly().toArray()).value(t));
		return output;
	}
}
