package ChemSimAdaptGC_2D;

import mpi.rc.IJ.IJutilities.MersenneTwister;
import mpi.rc.IJ.IJutilities.SimulationTime;

public class CellR12DA {

    CellPathwayR12DA cp;

    // current variables
    int t;
    double x,y,phi,v;
    boolean runState;

    CellR12DA(int t0, double x0, double y0, double phi0, double v, CellPathwayR12DA c0){

        t = t0;

        x=x0;
        y=y0;
        phi=phi0;
        this.v=v;

        cp = c0;

    }


    // get methods
    double getBias() {return CellPropertiesR12DA.returnCurrentBias(getBiasR()); }
    double getVX() { return v*Math.cos(phi); }
    double getVY(){
        return v*Math.sin(phi);
    }

    // update methods
    // to reach steady state
    void equilibrateNetwork(int nSteps, MersenneTwister rd){

        CellPropertiesR12DA.updateRunState(this,rd);

    }

    double getBiasR(){return cp.biasR(SimulationBoxPropertiesR12DA.conc(x));}

    // one simulation step
    void advanceStep(MersenneTwister rd){

        // update Network
        cp.updateBiasR(SimulationBoxPropertiesR12DA.conc(x));
        // not necessary if no motor adapt

        //prepare displacement computation
        double PhiM=0;

        double ux = Math.cos(phi);
        double uy = Math.sin(phi);

        synchronized(SimulationTime.holder){
            // rotational diffusion
            PhiM = CellPropertiesR12DA.dr * rd.nextGaussian();

            if(!runState){
                // reorientation in the tumble phase
                PhiM += CellPropertiesR12DA.dTheta *rd.nextGaussian();// dTheta is scaled by dt

            }
        }

        // interaction with surface z in ~[0 frameSizeZ]
        //*
        if(y<0 || y> SimulationBoxPropertiesR12DA.frameSizeY){
            double dy= Math.max(y- SimulationBoxPropertiesR12DA.frameSizeY,-y); //dz always +
            PhiM -= dy * CellPropertiesR12DA.sdTh * uy * ux;  // swimming in circle
            y -= v * uy + Math.signum(y) * CellPropertiesR12DA.sdTh * dy;
        }
        // */

        // straight movement in both phases, but motion reorients much faster during tumbles
        x += v * ux; // v is scaled by dt
        y += v * uy;

        CellPropertiesR12DA.updateRunState(this,rd);

        // position and orientation iteration
        phi += PhiM;

    }

    // updateTimeStep (useful to have separate if no save every time step)
    void updateTime(){
        t++;
    }

    // check if cell fell off simulation box or simulation is finished
    boolean finished(){
        return x<0 || x> SimulationBoxPropertiesR12DA.frameSizeX || t== SimulationBoxPropertiesR12DA.nFrames-1;
    }



}
