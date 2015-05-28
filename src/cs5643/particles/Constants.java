package cs5643.particles;

/**
 * Default constants. Add your own as necessary.
 *
 * @author Doug James, January 2007
 * @author Eston Schweickart, February 2014
 */
public interface Constants
{
    /** Mass of a particle. */
    public static final double PARTICLE_MASS     = 1.0;

    /** Camera rotation speed constants. */
    public static final double CAM_SIN_THETA     = Math.sin(0.2);
    public static final double CAM_COS_THETA     = Math.cos(0.2);
    
    public static final double tiny = 0.00000001;
    
    // Gravitational acceleration constant (point in -y direction)
    public static final double gravity = 10.0;
    
    // Number of iterations for main while loop where deltaP is calculated each time step
    public static final int solverIterations = 10;
    
    //Global density
    public static final double rho = 6378;
    
    //CFM relaxation parameter
    public static final double epsilon = 600.0;
    
    // Kernel radius
    public static final double h = 0.1;
    
    // Constants for surface tension artificial pressure, use enable to toggle effect
    public static final boolean s_corrEnable = true;
    public static final double s_corr_k = 0.0001;
    public static final double s_corr_n = 4;
    public static final double s_corr_deltaQ = 0.3*h;
    
    // Constants for viscous dissipation, use enable to toggle effect
    public static final boolean viscosityEnable = true;
    public static final double c = 0.00001;
    public static final double vd = 0.49;
    
    // Constants for adding kinetic energy with vorticity confinement, use enable to toggle effect
    public static final boolean vorticityEnable = true;
    public static final double vorticityEpsilon = 0.0005;
    
    
}
