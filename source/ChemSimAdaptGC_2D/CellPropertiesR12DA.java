package ChemSimAdaptGC_2D;

import ij.gui.GenericDialog;
import mpi.rc.IJ.IJutilities.MersenneTwister;
import mpi.rc.IJ.IJutilities.SimulationTime;

public class CellPropertiesR12DA {

    public static double dr; 					//sqrt(2 Dr)
    public static double dv0;
    public static double v0;

    public static double rtR;
    public static double trR;
    public static double dTheta;

    public static double adaptRateDt;

    public static double sdTh; // surface reorientation rate

    //public static int simuType;
    //public static final int UNIFORM_NO_ADAPT=0, VARIABLE_NO_ADAPT=1;

    public static double N0,M0,K0,sigmaM,sigmaK;


    static void updateRunState(CellR12DA c, MersenneTwister rd){

        double  ptr,prt;

        if(c.runState){

            prt = Math.exp(-rtR/c.getBiasR());

            // do we stay in the runs phase in the next step ?
            synchronized(SimulationTime.holder){c.runState=(rd.nextDouble()<prt);}
        }
        else{

            ptr = Math.exp(-trR);

            // do we stay tumble in the next step ?

            synchronized(SimulationTime.holder){c.runState=!(rd.nextDouble()<ptr);}  // since rd2.nextDouble()<ptr : proba to stay tumble.
        }

        //return runState;

    }

    static double returnCurrentBias(double bR){

        return rtR/(rtR + trR*bR);

    }

     static double[] instantiatePathwayParameters(MersenneTwister rd){

        double[] out = new double[3];

        out[0] = N0;
        out[1] = K0 * Math.exp( sigmaK * rd.nextGaussian() );

        double m = -1;
        while(m<=0){
            m = M0 + sigmaM * rd.nextGaussian() ;
        }
        out[2] = 1/m;

        return out;

    }


    static boolean setPropertiesGUI(){

        boolean canceled;


        double v0 = 0.2;
        double dv0 = 0.02;


        double dr = 0.001;

        double rtR = 0.06; // uniform cells value
        double trR = 0.06;  // fixed
        double dTheta = 0.12;

        double sdth = 0.5;

        double n0 = 10;
        double k0 = 750; // MeAsp
        double m0 = 0.20;
        double sigmak = 1.5;
        double sigmam = 0.5;

        double adaptRate = 1.0/(90.0*60.0)*0.01; // 1/90 mins

        GenericDialog gd = new GenericDialog("Bacterial Properties");
        gd.addNumericField("Norm_of_the_velocity (px/step)", v0, 2);
        gd.addNumericField("Std_Dev_of_the_norm_of_the_velocity (px/step)", dv0, 2);
        gd.addNumericField("Rot_diff_coeff (step-1): ", dr, 3);

        gd.addMessage("basal run tumble Parameters");
        gd.addNumericField("Run_to_Tumble_rate (step^-1)", rtR, 2);
        gd.addNumericField("Tumble_to_Run_rate (step^-1)", trR, 2);
        gd.addNumericField("Tumbling_Rotational_rate (rad2/step)", dTheta, 2);
        gd.addNumericField("surface_reorientation_rate (rad/step/px)", sdth, 1);

        gd.addMessage("Pathway Parameters");
        gd.addNumericField("N", n0, 1,6,"0.8 10");
        gd.addNumericField("K", k0, 1,6,"33 750");
        gd.addNumericField("M", m0, 2, 6, "0.37 0.20");
        gd.addNumericField("sigmaK", sigmak, 2,6,"0 1.5");
        gd.addNumericField("sigmaM", sigmam, 2, 6,"0 0.5"); // add adpat rate gui
        gd.addNumericField("adaptRate", adaptRate, 6, 6,"step^-1"); // add adpat rate gui
        //gd.addNumericField("N ", n0, 1);
        //gd.addNumericField("K ", k0, 1);
        //gd.addNumericField("M ", m0, 2);
        //gd.addNumericField("sigmaK ", sigmak, 2);
        //gd.addNumericField("sigmaM ", sigmam, 2);

        gd.showDialog();

        canceled = gd.wasCanceled();

        v0 = (double) gd.getNextNumber();
        dv0 = (double) gd.getNextNumber();
        dr = (double) gd.getNextNumber();

        rtR = (double) gd.getNextNumber();
        trR = (double) gd.getNextNumber();
        dTheta = (double) gd.getNextNumber();

        sdTh = (double) gd.getNextNumber();

        n0 = gd.getNextNumber();
        k0 = gd.getNextNumber();
        m0 = gd.getNextNumber();
        sigmak = gd.getNextNumber();
        sigmam = gd.getNextNumber();
        adaptRate = gd.getNextNumber();


        // assign values
        dr = Math.sqrt(2*dr); // sqrt(2Dr)
        dTheta = Math.sqrt(2*dTheta); // sqrt(2dTheta)

        CellPropertiesR12DA.setVelocity(v0,dv0);
        CellPropertiesR12DA.setDiffcoeff(dr);
        CellPropertiesR12DA.setRTprops(rtR,trR,dTheta);
        CellPropertiesR12DA.setPathwayProps(n0,k0,m0,sigmak,sigmam,adaptRate);


        return canceled;

    }

    static String getPropertiesLog(){
        String buffer = "";
        String lineSep = System.lineSeparator();

        buffer +="v0\t"+v0+lineSep;
        buffer +="dv0\t"+dv0+lineSep;
        buffer +="dr\t"+dr+lineSep;

        buffer +="rtR\t"+rtR+lineSep;
        buffer +="trR\t"+trR+lineSep;
        buffer +="dTheta\t"+dTheta+lineSep;

        buffer += "N0\t"+N0+lineSep;
        buffer += "K0\t"+K0+lineSep;
        buffer += "M0\t"+M0+lineSep;
        buffer += "sigmaK\t"+sigmaK+lineSep;
        buffer += "sigmaM\t"+sigmaM+lineSep;

        buffer += "adaptRate\t"+adaptRateDt+lineSep;

        return buffer;
    }


    static void setVelocity(double v, double dv){
        v0 = v;
        dv0 = dv;
    }
    static void setDiffcoeff(double D){
        dr=D;
    }
    static void setRTprops(double r2t,double t2r,double th){
        rtR = r2t;
        trR = t2r;
        dTheta = th;
    }
    static void setPathwayProps(double n0,double k0,double m0,double sigmak,double sigmam, double adaptRate){
        N0 =n0;
        K0 = k0;
        M0 = m0;
        sigmaK = sigmak;
        sigmaM = sigmam;
        adaptRateDt = adaptRate;
    }



}
