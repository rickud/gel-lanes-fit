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
import fityk.RealVector;

public class FitFityk extends Fityk {
    static {
        System.load("/Users/Rick/fityk-1.3.1/fityk/swig/java/libfitykJ.so");
    }

    public FitFityk(float x, float y) {
        Fityk f = new Fityk();
        f.load_data(0, new RealVector((long) x), new RealVector((long) y), null);
        PointVector points = f.get_data();
        for (int i = 0; i<points.size(); i++){
        	Point p = points.get(i);
        	if (p.getIs_active()) {
        		System.out.println(""+i+": [" + p.getX() + "" + p.getY() + "]");
        	}
        }
    }

    public FitFityk(float[] xValues, float[] yValues) {
        Fityk f = new Fityk();
        
	}

	public void run() {
        execute("guess %gauss = Gaussian");
        System.out.println("Fitting ...");
        execute("fit");
        System.out.println("WSSR=" + get_wssr());
        System.out.println("Gaussian center: " +
                           calculate_expr("%gauss.center"));
    }

    public void save_session(String sessionFilename) {
        execute(String.format("info state >'%s'", sessionFilename));
    }
}

