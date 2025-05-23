package ChemSimAdaptGC_2D;

import mpi.rc.IJ.IJutilities.MersenneTwister;

public class CellPathwayR12DA {

    double K; // binding affinity
    double N; // cooperativity
    double M0; // Prefactor

    double Mt; // current bias

    double[] props; // for convenience and speed

    CellPathwayR12DA(MersenneTwister rd){

        double[] p = CellPropertiesR12DA.instantiatePathwayParameters(rd);

        N = p[0];
        K = p[1];
        M0 = p[2];

        Mt = M0;

        props = new double[3];
        props[0] = N;
        props[1] = K;
        props[2] = M0;

    }

    CellPathwayR12DA(double n, double k, double m){
        N = n;
        K = k;
        M0 = m;
    }

    double biasR(double cc){

        return Mt * Math.pow((1+cc/K),N);

    }

    double[] getProps(){

        return props;

    }

    void updateBiasR(double cc){

        Mt += CellPropertiesR12DA.adaptRateDt * (M0/ Math.pow((1+cc/K),N) - Mt);

    }




}
