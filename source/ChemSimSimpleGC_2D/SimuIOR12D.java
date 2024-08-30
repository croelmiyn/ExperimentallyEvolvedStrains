package ChemSimSimpleGC_2D;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.io.SaveDialog;
import mpi.rc.IJ.IJutilities.ParallelUtil;
import mpi.rc.IJ.IJutilities.SimulationTime;

import java.io.*;
import java.nio.ByteBuffer;
import java.util.ArrayList;

public class SimuIOR12D {

    String saveNameDATA, saveNameLOG, movieName;

    boolean canceled;
    String lineSep = System.lineSeparator();

    // for movie
    ImagePlus imp;
    String movieType;
    double L;		// length of the particles
    double e;        // half-thickness of the particles

    //

    SimuIOR12D(){
        canceled=false;
        L=2;
        e=1;
    }

    void runInitGUI(){
        canceled = SimulationBoxPropertiesR12D.runInitGUI();
        canceled = CellPropertiesR12D.setPropertiesGUI();
        if(!canceled) {
            setSaveNames();
        }
    }

    void setSaveNames(){

        String sFileName = "Trajectories_SimuChemotaxis2_";;

        SaveDialog sd = new SaveDialog("Output_File_data",sFileName,".data");
        saveNameDATA = sd.getDirectory()+sd.getFileName();
        if(sd.getFileName()==null){saveNameDATA=null;}

        /*
        sd = new SaveDialog("Output_File_cdf",sFileName,".cdf");
        saveNameCDF = sd.getDirectory()+sd.getFileName();
        if(sd.getFileName()==null){saveNameCDF=null;}
        // */

        sFileName = "Log_SimuChemotaxis2_";

        sd = new SaveDialog("Log_File",sFileName,".txt");
        saveNameLOG = sd.getDirectory() + sd.getFileName();

        GenericDialog gd = new GenericDialog("Movie Params");
        gd.addStringField("MovieName","SimuChemotaxis2_");
        gd.addNumericField("L",L,0);
        gd.addNumericField("e",e,0);
        String[] type = new String[]{"8-bit","32-bit","none"};
        gd.addChoice("MovieType",type,type[0]);
        gd.showDialog();
        if(gd.wasCanceled()){
            movieName=null;
        }
        else{
            movieName = gd.getNextString();
            L = gd.getNextNumber();
            e = gd.getNextNumber();
            movieType = gd.getNextChoice();
        }


    }



    void saveDATA(ArrayList<TrajectoryR12D> trajs){

        if(saveNameDATA==null){
            return;
        }

        TrajectoryR12D traj;

        int nTraj = trajs.size();
        int nData = 0;
        for(int i=0; i<nTraj; i++){
            nData+=trajs.get(i).getSize();
        }

        try {
            FileOutputStream file = new FileOutputStream(saveNameDATA);

            // standard format = nb of data arrays, nFrames, ntrajs, t[], data1[][], data2[][], ...

            DataOutputStream dos = new DataOutputStream(file);
            int byteCount = 2 *4+10 *8+1;

            ByteBuffer bb = ByteBuffer.allocate(byteCount);
            //bb.order(ByteOrder.nativeOrder());

            dos.writeInt(nData);
            dos.writeInt(SimulationBoxPropertiesR12D.nFrames);
            dos.writeInt(nTraj);
            dos.writeInt(byteCount);

            for(int i=0; i<nTraj; i++){
                IJ.showProgress((double)i/nTraj);
                traj = trajs.get(i);
                for(int j=0; j<traj.getSize(); j++){

                    bb.putInt(traj.label);
                    bb.putInt(traj.t.get(j));
                    bb.putDouble(traj.x.get(j));
                    bb.putDouble(traj.y.get(j));
                    bb.putDouble(traj.phi.get(j));
                    bb.putDouble(traj.vX.get(j));
                    bb.putDouble(traj.vY.get(j));
                    bb.putDouble(traj.cc.get(j));
                    bb.put((byte)(traj.runState.get(j)?1:0));
                    bb.putDouble(traj.bias.get(j));
                    bb.putDouble(traj.chemProp.get(j)[0]);
                    bb.putDouble(traj.chemProp.get(j)[1]);
                    bb.putDouble(traj.chemProp.get(j)[2]);

                    dos.write(bb.array());

                    bb.clear();
                }
            }
            file.close();
        } catch (Exception e){
            IJ.log("Erreur doSave --> "+e.getMessage());
            IJ.log("Erreur doSave --> "+e.getCause());
            IJ.log("Erreur doSave --> "+e.getLocalizedMessage());
        }

        IJ.showStatus("Done");


    }

    ArrayList<TrajectoryR12D> readDATA(String fileName){

        ArrayList<TrajectoryR12D> trajs = new ArrayList<>();

        try{
            FileInputStream file = new FileInputStream(fileName);
            DataInputStream dis = new DataInputStream(file);

            int nData = dis.readInt();
            int nFr = dis.readInt();
            int nTraj = dis.readInt();
            int byteCount = dis.readInt();

            IJ.log("nb of Traj: "+nTraj);
            IJ.log("nb of Frames: "+nFr);

            int curT = -1;

            int tn,fn;
            double x,y,phi,vx,vy,bias,n0,k0,m0,cc;
            byte rs;

            ByteBuffer bb = ByteBuffer.allocate(byteCount);

            for(int i=0; i<nData; i++){
                IJ.showProgress((double)i/nData);

                dis.read(bb.array());

                tn = bb.getInt();
                fn = bb.getInt();
                x = bb.getDouble();
                y = bb.getDouble();
                phi = bb.getDouble();
                vx = bb.getDouble();
                vy = bb.getDouble();
                cc = bb.getDouble();
                rs = bb.get();
                bias = bb.getDouble();
                n0 = bb.getDouble();
                k0 = bb.getDouble();
                m0 = bb.getDouble();

                bb.clear();

                if(tn!=curT){
                    curT = tn;
                    TrajectoryR12D traj = new TrajectoryR12D(tn, new CellR12D(fn,x,y,phi,Math.sqrt(vx*vx+vy*vy),new CellPathwayR12D(n0,k0,m0)));
                    trajs.add(traj);
                }
                else {
                    trajs.get(tn).t.add(fn);
                    trajs.get(tn).x.add(x);
                    trajs.get(tn).y.add(y);
                    trajs.get(tn).phi.add(phi);
                    trajs.get(tn).vX.add(vx);
                    trajs.get(tn).vY.add(vy);
                    trajs.get(tn).cc.add(cc);
                    trajs.get(tn).runState.add(rs == 1);
                    trajs.get(tn).bias.add(bias);
                    trajs.get(tn).chemProp.add(new double[]{n0, k0, m0});
                }

            }


        } catch (Exception e){
            IJ.log("Erreur doLoad 1 --> "+e.getMessage());
            IJ.log("Erreur doLoad 2 --> "+e.getCause());
            IJ.log("Erreur doLoad 3 --> "+e.getLocalizedMessage());
        }
        IJ.showStatus("Done");

        return trajs;

    }

    void saveLOG(){

        String buffer = "Parameters Simulations 2"+lineSep;

        buffer += SimulationBoxPropertiesR12D.getParameterLog();

        buffer += CellPropertiesR12D.getPropertiesLog();

        buffer += "L_plot"+L+lineSep;
        buffer += "e_plot"+e+lineSep;

        try {
            FileWriter file = new FileWriter(saveNameLOG);
            file.write(buffer);

            file.close();
        } catch (Exception e){
            IJ.log("Erreur doSave --> "+e.getMessage());
            IJ.log("Erreur doSave --> "+e.getCause());
            IJ.log("Erreur doSave --> "+e.getLocalizedMessage());
        }

        IJ.showStatus("Done");

    }

    void generateMovie(ArrayList<TrajectoryR12D> trajs){

        if(movieName==null || movieType.equals("none")){
            return;
        }

        imp = IJ.createImage(movieName,movieType, SimulationBoxPropertiesR12D.frameSizeX, SimulationBoxPropertiesR12D.frameSizeY, SimulationBoxPropertiesR12D.nFrames);
        float[][] zeros = new float[SimulationBoxPropertiesR12D.frameSizeX][SimulationBoxPropertiesR12D.frameSizeX];
        for(int i = 0; i< SimulationBoxPropertiesR12D.nFrames; i++){
            imp.setSlice(i+1);
            imp.getProcessor().setFloatArray(zeros);
        }

        //ImageProcessor ip = imp.getProcessor();
        ImageStack stack = imp.getStack();

        ParallelUtil.parfor(trajs.size(),"printTraj",(k)->{

            TrajectoryR12D t = trajs.get(k);

            double x_per,y_per;
            double dr2, cz,X,Y, ck,sk;
            int x0,y0,xl,yl;
            double a = 1.0;

            for(int i=0;i<t.getSize();i++){

                x_per = t.x.get(i) - SimulationBoxPropertiesR12D.frameSizeX * Math.floor(t.x.get(i)/ SimulationBoxPropertiesR12D.frameSizeX);
                y_per = t.y.get(i) - SimulationBoxPropertiesR12D.frameSizeY * Math.floor(t.y.get(i)/ SimulationBoxPropertiesR12D.frameSizeY);

                x0 = (int) x_per;
                y0 = (int) y_per;

                cz = 0.25;

                ck = Math.cos(t.phi.get(i));
                sk = Math.sin(t.phi.get(i));

                for(int xx=-(int)(L+e); xx<=(L+e); xx++){
                    xl = (x0+xx)% SimulationBoxPropertiesR12D.frameSizeX;
                    for(int yy=-(int)(L+e); yy<=(L+e); yy++){
                        yl = (yy+y0)% SimulationBoxPropertiesR12D.frameSizeY;

                        X =  (xx+x0-x_per) * ck + (yy+y0-y_per) * sk;
                        Y = -(xx+x0-x_per) * sk + (yy+y0-y_per) * ck;

                        dr2 = 4*X*X/L/L + 4*Y*Y/e/e;

                        double I =  (255 * Math.exp(-a*dr2) * cz);

                        synchronized (SimulationTime.holder) {
                            I += stack.getVoxel(xl, yl, t.t.get(i));
                            stack.setVoxel(xl, yl,t.t.get(i), I);
                        }
                    }

                }

            }


        });

        imp.setSlice(1);

        imp.show();

        IJ.log("done");

    }
}
