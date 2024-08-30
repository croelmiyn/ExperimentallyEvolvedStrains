package ChemSimSimpleGC_2D;

import mpi.rc.IJ.IJutilities.MersenneTwister;
import mpi.rc.IJ.IJutilities.SimulationTime;

public class CellR12D {

    CellPathwayR12D cp;

    // current variables
    int t;
    double x,y,phi,v;
    boolean runState;

    CellR12D(int t0, double x0, double y0, double phi0, double v, CellPathwayR12D c0){

        t = t0;

        x=x0;
        y=y0;
        phi=phi0;
        this.v=v;

        cp = c0;

    }


    // get methods
    double getBias() {return CellPropertiesR12D.returnCurrentBias(getBiasR()); }
    double getVX() { return v*Math.cos(phi); }
    double getVY(){
        return v*Math.sin(phi);
    }

    // update methods
    // to reach steady state
    void equilibrateNetwork(int nSteps, MersenneTwister rd){

        CellPropertiesR12D.updateRunState(this,rd);

    }

    double getBiasR(){
        return cp.biasR(SimulationBoxPropertiesR12D.conc(x));
    }

    // one simulation step
    void advanceStep(MersenneTwister rd){

        // update Network
        // not necessary if no motor adapt

        //prepare displacement computation
        double PhiM=0;

        double ux = Math.cos(phi);
        double uy = Math.sin(phi);

        synchronized(SimulationTime.holder){
            // rotational diffusion
            PhiM = CellPropertiesR12D.dr * rd.nextGaussian();

            if(!runState){
                // reorientation in the tumble phase
                PhiM += CellPropertiesR12D.dTheta *rd.nextGaussian();// dTheta is scaled by dt

            }
        }

        // interaction with surface z in ~[0 frameSizeZ]
        //*
        if(y<0 || y>SimulationBoxPropertiesR12D.frameSizeY){
            double dy= Math.max(y-SimulationBoxPropertiesR12D.frameSizeY,-y); //dz always +
            PhiM -= dy * CellPropertiesR12D.sdTh * uy * ux;  // swimming in circle
            y -= v * uy + Math.signum(y) * CellPropertiesR12D.sdTh * dy;
        }
        // */

        // straight movement in both phases, but motion reorients much faster during tumbles
        x += v * ux; // v is scaled by dt
        y += v * uy;

        CellPropertiesR12D.updateRunState(this,rd);

        // position and orientation iteration
        phi += PhiM;

    }

    // updateTimeStep (useful to have separate if no save every time step)
    void updateTime(){
        t++;
    }

    // check if cell fell off simulation box or simulation is finished
    boolean finished(){
        return x<0 || x> SimulationBoxPropertiesR12D.frameSizeX || t== SimulationBoxPropertiesR12D.nFrames-1;
    }



}
