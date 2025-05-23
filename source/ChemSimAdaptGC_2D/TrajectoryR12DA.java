package ChemSimAdaptGC_2D;

import mpi.rc.IJ.IJutilities.MersenneTwister;

import java.util.ArrayList;

public class TrajectoryR12DA {

    static byte byteTrue = 1;
    static byte byteFalse = 0;


    CellR12DA cell;

    int label;

    ArrayList<Integer> t;
    ArrayList<Double> x,y,phi,cc,vX,vY,bias;
    ArrayList<Boolean> runState;
    ArrayList<double[]> chemProp;

    TrajectoryR12DA(CellR12DA c){
        label = (int)System.currentTimeMillis();
        t = new ArrayList<>();
        x = new ArrayList<>();
        y = new ArrayList<>();
        phi = new ArrayList<>();
        vX = new ArrayList<>();
        vY = new ArrayList<>();
        cc = new ArrayList<>();
        bias = new ArrayList<>();
        runState = new ArrayList<>();

        chemProp = new ArrayList<>();

        cell = c;

        add(cell);

    }

    TrajectoryR12DA(int l, CellR12DA c){
        label = l;
        t = new ArrayList<>();
        x = new ArrayList<>();
        y = new ArrayList<>();
        phi = new ArrayList<>();
        vX = new ArrayList<>();
        vY = new ArrayList<>();
        cc = new ArrayList<>();
        bias = new ArrayList<>();
        runState = new ArrayList<>();

        chemProp = new ArrayList<>();

        cell = c;

        add(cell);

    }

    void add(CellR12DA c){

        t.add(c.t);
        x.add(c.x);
        y.add(c.y);
        phi.add(c.phi);
        cc.add(SimulationBoxPropertiesR12DA.conc(c.x));

        vX.add(c.getVX());
        vY.add(c.getVY());

        bias.add(c.getBias());
        runState.add(c.runState);

        chemProp.add(c.cp.getProps());

    }


    // for simu

    void advanceCell(int nSteps, MersenneTwister rd){
        for(int i=0; i<nSteps;i++){
            cell.advanceStep(rd);
        }
        cell.updateTime();
        add(cell);
    }
    boolean isFinished(){
        return cell.finished();
    }


    int getSize(){
        return t.size();
    }






}
