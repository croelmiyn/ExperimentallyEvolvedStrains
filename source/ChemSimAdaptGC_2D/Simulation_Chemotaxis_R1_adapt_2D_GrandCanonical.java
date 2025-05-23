package ChemSimAdaptGC_2D;

import ij.IJ;
import ij.plugin.PlugIn;
import mpi.rc.IJ.IJutilities.MersenneTwister;
import mpi.rc.IJ.IJutilities.ParallelUtil;

import java.util.ArrayList;

public class Simulation_Chemotaxis_R1_adapt_2D_GrandCanonical implements PlugIn {

    ArrayList<TrajectoryR12DA> trajs;
    int trajCount;

    SimuIOR12DA simu_io;

    MersenneTwister rd;

    @Override
    public void run(String arg) {

        simu_io = new SimuIOR12DA();

        simu_io.runInitGUI();

        if(simu_io.canceled){
            return;
        }

        computeTrajectories();

        simu_io.saveLOG();
        simu_io.saveDATA(trajs);
        simu_io.generateMovie(trajs);

    }

    private void computeTrajectories(){


        int n0 = SimulationBoxPropertiesR12DA.getN0();
        rd = new MersenneTwister();
        trajCount = 0;
        trajs = new ArrayList<>();

        // initial particles

        for(int i=0; i<n0; i++){

            startNewTraj(0);

        }

        // rate of cells appearing between 2 frames at both ends of the channel
        double prodRate = CellPropertiesR12DA.v0* SimulationBoxPropertiesR12DA.normalizedApparitionRate2D();
        // in Cells per frame

        // Gillespie-like Apparition throughout the simulation
        double tC=1;
        while( tC< SimulationBoxPropertiesR12DA.nFrames){

            // pick time step
            double u1 = rd.nextDouble();
            double dt = -(1/prodRate)*Math.log(u1);// pick time step

            tC +=dt;
            startNewTraj((int)Math.floor(tC));

        }

        IJ.log("trajCount="+trajCount);

        // compute the evolution of all the trajectories
        ParallelUtil.parfor(trajCount,"CellTrajectoryCalculation",(k)->{

            TrajectoryR12DA traj = trajs.get(k);

            while(!traj.isFinished()){

                traj.advanceCell(SimulationBoxPropertiesR12DA.nSubSteps,rd);

            }

        });



    }

    private void startNewTraj(int time){

        if(time>= SimulationBoxPropertiesR12DA.nFrames){ // just in case
            return;
        }

        double x0,y0,phi0,v;

        if (time==0) {
            x0 = SimulationBoxPropertiesR12DA.frameSizeX * rd.nextDouble();
            y0 = SimulationBoxPropertiesR12DA.frameSizeY * rd.nextDouble();
            phi0 = 2 * Math.PI * rd.nextDouble();
        }
        else{
            y0 = SimulationBoxPropertiesR12DA.frameSizeY * rd.nextDouble();
            phi0 = 2 * Math.PI * rd.nextDouble();
            x0 = CellPropertiesR12DA.v0* SimulationBoxPropertiesR12DA.nSubSteps * rd.nextDouble();
            if(Math.cos(phi0)<0){ // going towards lower x, start at top
                x0 = SimulationBoxPropertiesR12DA.frameSizeX-x0;
            }
        }

        v = CellPropertiesR12DA.v0+ CellPropertiesR12DA.dv0 * rd.nextGaussian();

        CellPathwayR12DA cp = new CellPathwayR12DA(rd);

        TrajectoryR12DA t = new TrajectoryR12DA(trajCount, new CellR12DA(time,x0,y0,phi0,v,cp)); // ,Network2.getInitMeth(x0)
        trajCount++;

        trajs.add(t);

    }
}
