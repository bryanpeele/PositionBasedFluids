package cs5643.particles;
import javax.vecmath.*;

/*
 * Method used to calculate kernel gradients
 */

public class KerSpiky {

	/**
	 * 
	 * @param p_i: main particle
	 * @param p_j: neighbor of main particle
	 * @return: vector gradient based on particle positions
	 */
	public Vector3d getGradient(Particle p_i, Particle p_j){
		double r = p_i.xNew.distance(p_j.xNew);
		assert r >0.0;
		
		Vector3d gradient = new Vector3d();
		gradient.set(0.0,0.0,0.0);
		
		// Return zero vector if r > kernel radius, h or r = 0
		if (r > Constants.h) {
			return gradient;
		} else if (r < Constants.tiny) {
			return gradient;
		} else {
			double coeff = -(45.0/(Math.PI*Math.pow(Constants.h, 6))) * (Math.pow(Constants.h-r, 2)/r);
			gradient.sub(p_i.xNew, p_j.xNew);
			gradient.scale(coeff);
			return gradient;
		}
	}

	
	
	/**
	 * 
	 * @param p_i main particle
	 * @param p_k neighbor of main particle
	 * @return scalar quantity accumulated to solve denominator of Eq 9, Macklin
	 */
	public double get2Norm(Particle p_i, Particle p_k) {
		Vector3d  gradAccum = new Vector3d();
		gradAccum.set(0.0,0.0,0.0);
		
		//If particle i is particle k, iterate over neighbors of k
		if(p_i == p_k) {
			for (Particle p_j : p_i.neighbors) {
				gradAccum.scaleAdd(1.0/Constants.rho, getGradient(p_i,p_j),gradAccum);
			}
		} else {
			gradAccum.scale(-1.0/Constants.rho,getGradient(p_i,p_k));
		}
		//return square of norm of accumulated gradient vector
		return gradAccum.lengthSquared();
	}
	
}
