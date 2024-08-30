package ChemSimSimpleGC_2D;

import ij.IJ;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.io.SaveDialog;
import ij.plugin.PlugIn;

import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;

public class Analysis_SimuGC_2D_runDurHisto implements PlugIn {


    private SimuIOR12D sIO;
    private String dir,fileName,sFileNameH;

    private double xmin,xmax;
    private int tmin,tmax;

    private ArrayList<TrajectoryR12D> trajectories;

    private boolean canRun;

    private ArrayList<Integer> runUp,runDown;

    private double hist_run_Up[],hist_run_Down[];

    @Override
    public void run(String arg) {

        loadData();
        setParams();

        if(canRun){
            compute();
            saveH();
        }

    }

    // execution functions
    private void loadData(){

        sIO = new SimuIOR12D();

        OpenDialog od = new OpenDialog("dataFile");
        dir = od.getDirectory();
        fileName = od.getFileName();

        SaveDialog sd = new SaveDialog("saveFile","runDurHistograms_"+fileName,".txt");
        sFileNameH = sd.getDirectory()+sd.getFileName();

        trajectories = sIO.readDATA(dir+fileName);

    }

    private void setParams(){

        xmin = 0;
        xmax = getXmax(trajectories);
        tmin = 0;
        tmax = getTmax(trajectories);

        GenericDialog gd = new GenericDialog("Parameters");
        gd.addNumericField("xmin",xmin,0);
        gd.addNumericField("xmax",xmax,0);
        gd.addNumericField("tmin",tmin,0);
        gd.addNumericField("tmax",tmax,0);

        gd.showDialog();

        if(gd.wasCanceled()){return;}

        xmin = gd.getNextNumber();
        xmax = gd.getNextNumber();
        tmin = (int) gd.getNextNumber();
        tmax = (int) gd.getNextNumber();

        canRun = true;

    }

    private void compute(){


        hist_run_Up = new double[50];
        hist_run_Down = new double[50];


        runUp = new ArrayList<>();
        runDown = new ArrayList<>();

        int tc;
        double dx;
        int runDur;

        boolean isIn, wasIn, wasRun;

        for(TrajectoryR12D traj:trajectories){

            dx = 0;
            runDur = 0;

            wasIn = false;
            wasRun = false;

            for(int i=0;i< traj.getSize();i++) {

                tc = traj.t.get(i);

                isIn = traj.x.get(i) >= xmin && traj.x.get(i) <= xmax && tc >= tmin && tc <= tmax;

                if (isIn && wasIn) {

                    if (traj.runState.get(i)) { // is run

                        if (wasRun) {  // run continues
                            runDur++;
                        } else {  // new run
                            runDur = 1;
                            wasRun = true;
                            dx = 0;
                        }
                        dx += traj.vX.get(i);
                    } else {  // is tumble

                        if (wasRun) { // run interupts
                            if (dx > 0) {
                                runUp.add(runDur);
                            } else {
                                runDown.add(runDur);
                            }
                            runDur = 0;
                            dx =0;
                        }
                        wasRun = false;
                    }
                } else if (isIn && !wasIn) { // enters

                    dx = 0;
                    runDur = 0;
                    wasIn = true;

                    if (traj.runState.get(i)) { // runs

                        runDur = 1;
                        wasRun = true;

                        dx += traj.vX.get(i);
                    }
                } else if (wasIn && !isIn) { //leaves
                    if (wasRun) {
                        if (dx > 0) {
                            runUp.add(runDur);
                        } else {
                            runDown.add(runDur);
                        }
                        runDur = 0;
                    }
                    wasIn=false;
                    dx=0;
                    wasRun = false;

                }

            } // end of the trajectory

            if (wasRun && wasIn) {
                if (dx > 0) {
                    runUp.add(runDur);
                } else {
                    runDown.add(runDur);
                }
            }

        }

        // make histograms

        hist_run_Up = makeHistogram(runUp);
        hist_run_Down = makeHistogram(runDown);

    }


    private void saveH(){
        String lineSep = System.lineSeparator();
        String buffer;
        try {
            File f = new File(sFileNameH);
            FileWriter file = new FileWriter(f); // makes new file

            buffer = "runDur(fr)\tHistogramUp\tHistogramDown"+lineSep;
            file.write(buffer);

            int max = hist_run_Down.length>hist_run_Up.length?hist_run_Down.length:hist_run_Up.length;

            for(int i=0;i<max;i++) {
                buffer = "" + i ;

                if(i<hist_run_Up.length){
                buffer+= "\t" + hist_run_Up[i];
                }
                else{buffer+= "\t";}

                if(i<hist_run_Down.length){buffer+= "\t" + hist_run_Down[i];}
                else{buffer+= "\t";}

                buffer+= lineSep;
                file.write(buffer);
            }

            file.close();


        }catch(Exception e){
            IJ.log("Error savePopulationStats --> "+e.getMessage());
            IJ.log("Error savePopulationStats --> "+e.getCause());
            IJ.log("Error savePopulationStats --> "+e.getLocalizedMessage());
        }


    }


    // private utilities functions
    private double getXmax(ArrayList<TrajectoryR12D> ts){
        double x=0;

        for(TrajectoryR12D t:ts){
            for(double xx:t.x){
                if(xx>x){
                    x=xx;
                }
            }
        }

        return x;
    }

    private int getTmax(ArrayList<TrajectoryR12D> ts){
        int x=0;

        for(TrajectoryR12D t:ts){
            for(int xx:t.t){
                if(xx>x){
                    x=xx;
                }
            }
        }

        return x;
    }

    private double[] makeHistogram(ArrayList<Integer> a){ // quick histogram making

        int max = getMax(a);

        double[] h = new double[max+1];

        for(int i:a){
            h[i]++;
        }

        return h;

    }

    private int getMax(ArrayList<Integer> a){
        int m=a.size()==0?0:a.get(0);
        for(int i:a){
            if(i>m){m=i;}

        }
        return m;
    }

}
