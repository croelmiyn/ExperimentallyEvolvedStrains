package ChemSimSimpleGC_2D;

import ij.IJ;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.io.SaveDialog;
import ij.plugin.PlugIn;

import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;

public class Analysis_SimuGC_2D implements PlugIn {


    private SimuIOR12D sIO;
    private String dir,fileName,sFileName,sFileNameT,sFileNameH;

    private double xmin,xmax;
    private int tmin,tmax;

    private double biasMax;
    private int nFrames;

    private ArrayList<TrajectoryR12D> trajectories;

    private boolean canRun;

    private double vX,biasUp,biasDown,fractionUp,fractionDown;
    private double count, countOut;

    private double[] vX_t,biasUp_t,biasDown_t,fractionUp_t,n_t;

    private double hist_b[];

    @Override
    public void run(String arg) {

        loadData();
        setParams();

        if(canRun){
            compute();
            save();
            saveT();
            saveH();
        }

    }

    // execution functions
    private void loadData(){

        sIO = new SimuIOR12D();

        OpenDialog od = new OpenDialog("dataFile");
        dir = od.getDirectory();
        fileName = od.getFileName();

        SaveDialog sd = new SaveDialog("saveFile","Analysis_SimuGC",".txt");
        sFileName = sd.getDirectory()+sd.getFileName();

        sd = new SaveDialog("saveFileT","TimeDepData_"+fileName,".txt");
        sFileNameT = sd.getDirectory()+sd.getFileName();

        sd = new SaveDialog("saveFileH","BiasHistogram_"+fileName,".txt");
        sFileNameH = sd.getDirectory()+sd.getFileName();

        trajectories = sIO.readDATA(dir+fileName);

    }

    private void setParams(){
        nFrames = getTmax(trajectories)+1;
        xmin = 0;
        xmax = getXmax(trajectories);
        tmin = 0;
        tmax = getTmax(trajectories);
        biasMax = 1;

        GenericDialog gd = new GenericDialog("Parameters");
        gd.addNumericField("xmin",xmin,0);
        gd.addNumericField("xmax",xmax,0);
        gd.addNumericField("tmin",tmin,0);
        gd.addNumericField("tmax",tmax,0);
        gd.addNumericField("maxBias",biasMax,2);

        gd.showDialog();

        if(gd.wasCanceled()){return;}

        xmin = gd.getNextNumber();
        xmax = gd.getNextNumber();
        tmin = (int) gd.getNextNumber();
        tmax = (int) gd.getNextNumber();
        biasMax = gd.getNextNumber();

        canRun = true;

    }

    private void compute(){
        vX=0;
        biasUp=0;
        biasDown=0;
        int cnt = 0, cntU=0,cntD=0;
        countOut=0;
        double b;

        vX_t = new double[nFrames];
        biasUp_t = new double[nFrames];
        biasDown_t = new double[nFrames];
        fractionUp_t = new double[nFrames];
        n_t = new double[nFrames];

        hist_b = new double[50];

        int tc;
        double trajRun,trajTumble;

        for(TrajectoryR12D traj:trajectories){

            trajRun=0;
            trajTumble=0;

            for(int i=0;i< traj.getSize();i++){

                if(traj.x.get(i)>=xmin && traj.x.get(i)<=xmax ){

                    b = traj.bias.get(i);

                    // time resolved
                    tc = traj.t.get(i);

                    n_t[tc]++;
                    vX_t[tc] += traj.vX.get(i);
                    if (traj.vX.get(i) > 0) {
                        biasUp_t[tc] += b;
                        fractionUp_t[tc]++;
                    } else {
                        biasDown_t[tc] += b;
                    }

                    if(tc>=tmin && tc<=tmax ) {

                        // bias histogram
                        if(traj.runState.get(i)){trajRun++;}
                        else{trajTumble++;}

                        // mean measurements
                        if (b < biasMax) {

                            vX += traj.vX.get(i);
                            cnt++;


                            if (traj.vX.get(i) > 0) {
                                biasUp += b;
                                cntU++;
                            } else {
                                biasDown += b;
                                cntD++;
                            }
                        }
                        else{
                            countOut++;
                        }
                    }


                }

            }

            // bias histogram
            if(trajRun+trajTumble!=0){
                b = trajTumble/(trajRun+trajTumble);
                tc = (int)Math.floor(b*50);
                if(tc==50){tc=49;}
                hist_b[tc]+=(trajRun+trajTumble);
            }

        }

        // normalizations

        // mean
        vX/=(cnt==0)?1:cnt;
        biasUp/=(cntU==0)?1:cntU;
        biasDown/=(cntD==0)?1:cntD;
        fractionUp = (cnt==0)?cntU:((double)cntU/cnt);
        fractionDown = (cnt==0)?cntD:((double)cntD/cnt);
        count = cnt;

        // time dependent
        for(int t=0; t<nFrames; t++){
            vX_t[t]/=n_t[t];
            biasDown_t[t]/=n_t[t]-fractionUp_t[t];
            biasUp_t[t]/=fractionUp_t[t];
            fractionUp_t[t]/=n_t[t];
        }

        // histogram
        cnt=0;
        for(double d:hist_b){
            cnt+=d;
        }
        for(int i=0;i<50;i++){
            hist_b[i]/=cnt*0.02; // probability density
        }
    }


    private void save(){
        String lineSep = System.lineSeparator();
        String buffer;
        try {
            File f = new File(sFileName);
            if (f.isFile()) {
                FileWriter file = new FileWriter(f,true); // appends existing file (for batch process)

                buffer = ""+fileName+"\t"+vX+"\t"+biasUp+"\t"+biasDown+"\t"+fractionUp+"\t"+fractionDown+"\t"+count+"\t"+countOut+"\t"+biasMax+lineSep;
                file.write(buffer);
                file.close();
            }
            else{
                FileWriter file = new FileWriter(f); // makes new file

                buffer = "FileName\tDrift\ttumbleBiasUpGrad\ttumbleBiasDownGrad\tfractionUp\tfractionDown\tNdetectedCells\tNexcludedCells\tBiasThreshold"+lineSep;
                file.write(buffer);

                buffer = ""+fileName+"\t"+vX+"\t"+biasUp+"\t"+biasDown+"\t"+fractionUp+"\t"+fractionDown+"\t"+count+"\t"+countOut+"\t"+biasMax+lineSep;
                file.write(buffer);
                file.close();
            }

        }catch(Exception e){
            IJ.log("Error savePopulationStats --> "+e.getMessage());
            IJ.log("Error savePopulationStats --> "+e.getCause());
            IJ.log("Error savePopulationStats --> "+e.getLocalizedMessage());
        }


    }


    private void saveT(){
        String lineSep = System.lineSeparator();
        String buffer;
        try {
            File f = new File(sFileNameT);
            FileWriter file = new FileWriter(f); // makes new file

            buffer = "t\tDrift\ttumbleBiasUpGrad\ttumbleBiasDownGrad\tfractionUp\tNdetectedCells"+lineSep;
            file.write(buffer);

            for(int t=0;t<nFrames;t++) {
                buffer = "" + t + "\t" + vX_t[t] + "\t" + biasUp_t[t] + "\t" + biasDown_t[t] + "\t" + fractionUp_t[t] + "\t" + n_t[t] + lineSep;
                file.write(buffer);
            }

            file.close();


        }catch(Exception e){
            IJ.log("Error savePopulationStats --> "+e.getMessage());
            IJ.log("Error savePopulationStats --> "+e.getCause());
            IJ.log("Error savePopulationStats --> "+e.getLocalizedMessage());
        }


    }

    private void saveH(){
        String lineSep = System.lineSeparator();
        String buffer;
        try {
            File f = new File(sFileNameH);
            FileWriter file = new FileWriter(f); // makes new file

            buffer = "bias\tpdf"+lineSep;
            file.write(buffer);

            for(int i=0;i<50;i++) {
                buffer = "" + (0.01+0.02*i) + "\t" + hist_b[i] + lineSep;
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

}
