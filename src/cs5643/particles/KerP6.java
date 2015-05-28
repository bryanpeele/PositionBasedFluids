package cs5643.particles;

/**
 * Class used to calculate poly 6 kernel, used for all non-gradient kernel funtions in this
 * implementation
 * @author Peele
 *
 */

public class KerP6 {
	
	// Method for using two particles as input, distance between xNew for each particle is passed to getKernel()
	public double getKernel(Particle p_i, Particle p_j){
		double r = p_i.xNew.distance(p_j.xNew);
		return getKernel(r);
	}
	
	// Method using scaler input as input
	// Returns 0 if distance between particles is greater than kernel radius, h
	public double getKernel(double r) {
		assert r >0.0;
		if(r <= Constants.h) {
			return (315.0/(64.0*Math.PI*Math.pow(Constants.h, 9))) * Math.pow( (Math.pow(Constants.h,2)-Math.pow(r,2)), 3);
		} else {
			return 0.0;
		}
	}
	
	
}
