// Example of using libfityk from Java.
// To run this example:
// - compile dynamic module for Java: after compiling fityk go to src/,
//   run "make java" with proper paths. On Fedora Linux it means:
//   $ JAVAINC=/usr/lib/jvm/java/include
//   $ make java CPPFLAGS="-I$JAVAINC -I/$JAVAINC/linux"
//   Put fitykJ.so into a directory in dynamic library path.
//   (libfityk.so should also be installed.)
// - java files in swig/java are in "package fityk", so compile the files:
//   $ cd ..../src/swig/java
//   $ javac *.java
//   and put them in a directory named fityk in CLASSPATH.
//   Link will also do:
//   $ cd ..../samples
//   $ ln -s .../src/swig/java ./fityk; export CLASSPATH=.
//   Now you can try this sample:
//   $ java hello

import fityk.Fityk;
import fityk.Point;
import fityk.PointVector;

public class FitFityk extends Fityk {
    double[] xout;
    double[] yout;
    static {
        System.load("/Users/Rick/fityk-1.3.1/fityk/swig/java/libfitykJ.so");
    }

    public FitFityk(float[] x, float[] y) {
        Fityk f = new Fityk();
        for (int i = 0; i<x.length; i++){
        	f.add_point(x[i], y[i], 0.0);
        }
        
	}

	public void runLaneFit() {
		System.out.println("Fitting ...");
		execute("guess %q0 = Quadratic");
		execute("fit");
		for (int i = 0; i < 10; i++) {
			execute("guess %g" +i+ " = Gaussian");
            execute("fit");
            System.out.println("WSSR ["+ i +"] =" + get_wssr());
            System.out.println("Center ["+ i +"]: " +
            				calculate_expr("%g" + i + ".center"));
		}
    }
	
	public double[] getXdata() { 
		PointVector points = get_data();
		xout = new double[(int) points.size()];
		for (int p=0; p<xout.length; p++)
			xout[p]=points.get(p).getX();
		return xout;
	}
	
	public double[] getYdata() { 
		PointVector points = get_data();
		yout = new double[(int) points.size()];
		for (int p=0; p<yout.length; p++)
			yout[p]=points.get(p).getX();
		return yout;
	}
	
    public void save_session(String sessionFilename) {
        execute(String.format("info state >'%s'", sessionFilename));
    }
}

