package cs5643.particles;

import java.util.*;
import javax.vecmath.*;
import javax.media.opengl.*;
import com.jogamp.opengl.util.glsl.*;


/**
 * Maintains dynamic lists of Particle and Force objects, and provides
 * access to their state for numerical integration of dynamics.
 *
 * @author Doug James, January 2007
 * @author Eston Schweickart, February 2014
 */
public class ParticleSystem //implements Serializable
{
	
	public KerP6 kp6 = new KerP6();
	public KerSpiky ks = new KerSpiky();
	
	
    /** Current simulation time. */
    public double time = 0;

    /** List of Particle objects. */
    public ArrayList<Particle>   P = new ArrayList<Particle>();

    /** List of Force objects. */
    public ArrayList<Force>      F = new ArrayList<Force>();

    /**
     * true iff prog has been initialized. This cannot be done in the
     * constructor because it requires a GL2 reference.
     */
    private boolean init = false;

    /** Filename of vertex shader source. */
    public static final String[] VERT_SOURCE = {"vert.glsl"};

    /** Filename of fragment shader source. */
    public static final String[] FRAG_SOURCE = {"frag.glsl"};

    /** The shader program used by the particles. */
    ShaderProgram prog;


    /** Basic constructor. */
    public ParticleSystem() {}

    /**
     * Set up the GLSL program. This requires that the current directory (i.e. the package in which
     * this class resides) has a vertex and fragment shader.
     */
    public synchronized void init(GL2 gl) {
        if (init) return;

        prog = new ShaderProgram();
        ShaderCode vert_code = ShaderCode.create(gl, GL2ES2.GL_VERTEX_SHADER, 1, this.getClass(), VERT_SOURCE, false);
        ShaderCode frag_code = ShaderCode.create(gl, GL2ES2.GL_FRAGMENT_SHADER, 1, this.getClass(), FRAG_SOURCE, false);
        if (!prog.add(gl, vert_code, System.err) || !prog.add(gl, frag_code, System.err)) {
            System.err.println("WARNING: shader did not compile");
            prog.init(gl); // Initialize empty program
        } else {
            prog.link(gl, System.err);
        }

        init = true;
    }

    /** Adds a force object (until removed) */
    public synchronized void addForce(Force f) {
        F.add(f);
    }

    /** Useful for removing temporary forces, such as user-interaction
     * spring forces. */
    public synchronized void removeForce(Force f) {
        F.remove(f);
    }

    /** Creates particle and adds it to the particle system.
     * @param p0 Undeformed/material position.
     * @return Reference to new Particle.
     */
    public synchronized Particle createParticle(Point3d p0)
    {
        Particle newP = new Particle(p0);
        P.add(newP);
        return newP;
    }

    /**
     * Helper-function that computes nearest particle to the specified
     * (deformed) position.
     * @return Nearest particle, or null if no particles.
     */
    public synchronized Particle getNearestParticle(Point3d x)
    {
        Particle minP      = null;
        double   minDistSq = Double.MAX_VALUE;
        for(Particle particle : P) {
            double distSq = x.distanceSquared(particle.x);
            if(distSq < minDistSq) {
                minDistSq = distSq;
                minP = particle;
            }
        }
        return minP;
    }

    /** Moves all particles to undeformed/materials positions, and
     * sets all velocities to zero. Synchronized to avoid problems
     * with simultaneous calls to advanceTime(). */
    public synchronized void reset()
    {
        for(Particle p : P)  {
            p.x.set(p.x0);
            p.v.set(0,0,0);
            p.f.set(0,0,0);
            p.setHighlight(false);
        }
        time = 0;
    }

    /**
     * Incomplete/Debugging implementation of Forward-Euler
     * step. WARNING: Contains buggy debugging forces.
     */
    public synchronized void advanceTime(double dt)
    {
    	
       //Apply gravity as external force
       for(Particle p : P)   {
    	   p.f.x = 0;
    	   p.f.z = 0;
    	   p.f.y = -p.m * Constants.gravity;
           p.v.scaleAdd(dt, p.f, p.v); //p.v += dt * p.f;
           p.xNew.scaleAdd(dt, p.v, p.x); //p.x += dt * p.v;
        }
       
       // Find neighbors for all particles
       for(Particle p_i : P) {
    	   p_i.neighbors.clear();
    	   for (Particle p_j : P) {
    		   // Add particles that are close enough to be neighbors (e.g. r < h)
    		   double distanceTmp = p_i.xNew.distance(p_j.xNew);
    		   if (distanceTmp <= 2*Constants.h) p_i.neighbors.add(p_j);/// CONSERVATIVE using 2h rather than h to account for particle movement
    	   }
       }
       
       
       // While loop used to calculate deltaP (change in position) based on density constraint
       int iter = 0;
       while (iter < Constants.solverIterations){
       
    	   //Calculate rho, lambda, norm2 for each particle
    	   for (Particle p_i : P) {
    		   p_i.rho = 0;
    		   p_i.lambda = 0;
    		   p_i.norm2 = 0;
    		   
    		   //Accumulate rho, norm2 over neighbors of p_i
    		   for (Particle p_j : p_i.neighbors) {
    			   p_i.rho += p_i.m*kp6.getKernel(p_i,p_j);
    			   p_i.norm2 += ks.get2Norm(p_i, p_j);
    		   }
    		   p_i.C = (p_i.rho/Constants.rho) - 1.0;
    		   p_i.lambda = -p_i.C/(p_i.norm2 +Constants.epsilon);
    	   }

    	   //Calculate deltaP
		   for (Particle p_i : P){
			   p_i.deltaP.set(0.0,0.0,0.0);
			   // Accumulate deltaP over neighbors of p_i
			   for (Particle p_j : p_i.neighbors) {
				   double s_corr = 0.0;
				   if (Constants.s_corrEnable) {
					   s_corr = -Constants.s_corr_k*Math.pow(kp6.getKernel(p_i,p_j)/kp6.getKernel(Constants.s_corr_deltaQ), Constants.s_corr_n);
					   //s_corr = 0.0001;
				   }
				   Utils.acc(p_i.deltaP,(p_i.lambda+p_j.lambda+s_corr)/Constants.rho, ks.getGradient(p_i,p_j));
			   }
		   }
    	   
    	   //Perform update of position with deltaP, then perform collision detection with box for all particles
    	   for (Particle p : P) {  
    		   // Update position based on deltaP (density constraint)
    		   p.xNew.add(p.deltaP);
    		   
    		   // Collision detection with bounding box, if component of position is outside bounds, set it equal to boundary
    		   if (p.xNew.x < 0.0) p.xNew.x = 0.0;
    		   if (p.xNew.x > 1.0) p.xNew.x = 1.0;
    		   if (p.xNew.y < 0.0) p.xNew.y = 0.0;
    		   if (p.xNew.y > 1.0) p.xNew.y = 1.0;
    		   if (p.xNew.z < 0.0) p.xNew.z = 0.0;
    		   if (p.xNew.z > 1.0) p.xNew.z = 1.0;
    		   
    	   }
   	      	   
    	   iter++;
       }
        // Calculate deltaX = xNew - x and use to set updated velocity   
       	for (Particle p_i : P) {
       		p_i.deltaX.sub(p_i.xNew,p_i.x);
       		p_i.v.scaleAdd(1.0/dt, p_i.deltaX, p_i.v);
       		p_i.v.scale(Constants.vd); //Scale factor used to dissipate excessive energy (1: no effect, 0: kills velocity)
       	}	
       	
       	//Apply vorticity confinement (optional)
       	if (Constants.vorticityEnable) {
       		for (Particle p_i : P) {
       			p_i.omega.set(0.0,0.0,0.0);
       			for(Particle p_j : P) {
       				p_i.TEMP.sub(p_j.v,p_i.v);
       				p_i.TEMP2.cross(p_i.TEMP,ks.getGradient(p_i,p_j));
       				p_i.omega.add(p_i.TEMP2);
       			}
       		}
       		for (Particle p_i : P) {
       			p_i.eta.set(0.0,0.0,0.0);
       			for(Particle p_j : P) {
       				p_i.eta.scaleAdd(p_j.omega.length(),ks.getGradient(p_i, p_j),p_i.eta);
       			}
       			p_i.N.set(p_i.eta);
       			p_i.N.scale(1.0/p_i.eta.length());
       			p_i.f_vort.cross(p_i.N,p_i.omega);
       			p_i.f_vort.scale(Constants.vorticityEpsilon);
       			p_i.v.scaleAdd(dt, p_i.f_vort, p_i.v);
       		}
       	}
       	
       	//Apply viscosity smoothing  (optional)
       	if (Constants.viscosityEnable){
       		for (Particle p_i : P) {
       		   	p_i.vNew.set(p_i.v);
       			for (Particle p_j : p_i.neighbors) {
       				p_i.TEMP.sub(p_j.v,p_i.v);
       				p_i.vNew.scaleAdd(Constants.c*kp6.getKernel(p_i,p_j),p_i.TEMP,p_i.vNew);
       			}
       		}
       		for (Particle p_i : P) {
       			p_i.v.set(p_i.vNew);
       		}
       	}
	
       	// Update particle positions
        for (Particle p_i : P) {
        	p_i.x.set(p_i.xNew);
        
       	}
        time += dt;
    }

    /**
     * Displays Particle and Force objects.
     */
    public synchronized void display(GL2 gl)
    {
        for(Force force : F) {
            force.display(gl);
        }

        if(!init) init(gl);

        prog.useProgram(gl, true);

        for(Particle particle : P) {
            particle.display(gl);
        }

        prog.useProgram(gl, false);
    }
}
