package ChemSimAdaptGC_2D;

import ij.gui.GenericDialog;

public class SimulationBoxPropertiesR12DA {

    // time
    static int nFrames; // nb frame saved
    static double dt; // elementary step duration
    static int nSubSteps; // 1 frame saved every X steps
    // space
    static int frameSizeX;
    static int frameSizeY;
    static double density;

    // attractant concentration
    static double grad;
    static double avgCC;

    SimulationBoxPropertiesR12DA(){
    }

    public static double conc(double x){

        double ratio = (1 + grad * (2*x - frameSizeX));
        if (ratio>2) ratio = 2;
        else if (ratio<0) ratio = 0;

        return avgCC * ratio; // cmax/2 + cmax*slope*(x-xm) = avgc (1+2*slope*(x-xm))

    }

    public static int getN0(){
        return (int) (density*frameSizeY*frameSizeX);
    }

    public static boolean runInitGUI(){

        nFrames = 5000;
        dt = 0.01;
        nSubSteps = 10;

        frameSizeX = 2048;
        frameSizeY = 512;
        density = 0.005;

        grad = 0.0005;
        avgCC = 50;

        GenericDialog gd = new GenericDialog("Simulation Box Properties");
        gd.addNumericField("Time_step_duration (in s) : ", dt, 2);
        gd.addNumericField("Save_every_X_steps", nSubSteps, 0);
        gd.addNumericField("Length_of_the_Simu (in saved frames) : ", nFrames, 0);
        gd.addNumericField("BoxSize_X (px)", frameSizeX, 0);
        gd.addNumericField("BoxSize_Y (px)", frameSizeY, 0);
        gd.addNumericField("Density_of_Objects (px-2)", density, 2);
        gd.addNumericField("gradSlope_X (px-1):", grad, 4);
        gd.addNumericField("avgC (uM):", avgCC, 0);

        gd.showDialog();

        dt = gd.getNextNumber();
        nSubSteps = (int) gd.getNextNumber();
        nFrames = (int) gd.getNextNumber();
        frameSizeX = (int) gd.getNextNumber();
        frameSizeY = (int) gd.getNextNumber();
        density = gd.getNextNumber();
        grad = gd.getNextNumber();
        avgCC = gd.getNextNumber();

        return gd.wasCanceled();

    }

    static String getParameterLog(){
        String buffer="";
        String lineSep = System.lineSeparator();

        buffer +="dt\t"+dt+lineSep;
        buffer +="stepsPerFrame\t"+nSubSteps+lineSep;
        buffer +="nFrames\t"+nFrames+lineSep;

        buffer +="frameSizeX\t"+frameSizeX+lineSep;
        buffer +="frameSizeY\t"+frameSizeY+lineSep;
        buffer +="density\t"+density+lineSep;

        buffer +="gradC\t"+grad+lineSep;
        buffer +="avgCC\t"+avgCC+lineSep;

        return buffer;
    }

    static double normalizedApparitionRate2D(){
        return 2*density* nSubSteps *frameSizeY/(Math.PI);
    }



}
