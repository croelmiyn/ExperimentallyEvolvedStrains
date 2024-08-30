package ChemSimSimpleGC_2D;

import mpi.rc.IJ.IJutilities.MersenneTwister;

public class CellPathwayR12D {

    double K; // binding affinity
    double N; // cooperativity
    double M0; // Prefactor

    double[] props; // for convenience and speed

    CellPathwayR12D(MersenneTwister rd){

        double[] p = CellPropertiesR12D.instantiatePathwayParameters(rd);

        N = p[0];
        K = p[1];
        M0 = p[2];

        props = new double[3];
        props[0] = N;
        props[1] = K;
        props[2] = M0;

    }

    CellPathwayR12D(double n, double k, double m){
        N = n;
        K = k;
        M0 = m;
    }

    double biasR(double cc){

        return M0 * Math.pow((1+cc/K),N);

    }

    double[] getProps(){

        return props;

    }




}
