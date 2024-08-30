dir = "E:\\Remy\\Simu_Manika\\2024-06-12_SimuGC_Pheno2D\\heterogeneousPop\\";


c = newArray(0,100,1000);

for(i=0; i<c.length; i++){

for(k=0; k<9; k++){

name = "R1-c-"+c[i]+"_"+k;


params = "time_step_duration=0.01 save_every_x_steps=10 length_of_the_simu=5000 ";
params+= "boxsize_x=2048 boxsize_y=512 ";
params+= "density_of_objects=0.01 ";
params+= "gradslope_x=0.0005 avgc="+(c[i]/2)+" ";
params+= "norm_of_the_velocity=0.25 std_dev_of_the_norm_of_the_velocity=0.025 ";
params+= "rot_diff_coeff=0.001 run_to_tumble_rate_(step^-1)=0.06 tumble_to_run_rate_(step^-1)=0.06 tumbling=0.12 ";
params+= "surface_reorientation_rate=0.5 ";

params+= "n=10 k=750 m=0.20 sigmak=1.50 sigmam=0.50 ";  // heterogeneous
// params+= "n=0.8 k=33 m=0.37 sigmak=0 sigmam=0 ";  // homogeneous

params+= "output_file_data="+dir+"Trajectories_SimuChemotaxisGC2D_"+name+".data ";
params+= "log_file="+dir+"Log_SimuChemotaxisGC2D_"+name+".txt ";
//params+= "moviename=SimuChemotaxisGC2D_"+name+" l=2 e=1 h=25 dh=15 movietype=8-bit";

run("Simulation Chemotaxis R1 2D GrandCanonical", params);

makeRectangle(768, 0, 512, 512);
run("Crop");
run("Save", "save="+dir+"SimuChemotaxisGC2D_"+name+".tif");
close();

}
}


dir = "E:\\Remy\\Simu_Manika\\2024-06-12_SimuGC_Pheno2D\\homogeneousPop\\";


c = newArray(0,100,1000);

for(i=0; i<c.length; i++){

for(k=0; k<9; k++){

name = "R1-c-"+c[i]+"_"+k;

params = "time_step_duration=0.01 save_every_x_steps=10 length_of_the_simu=5000 ";
params+= "boxsize_x=2048 boxsize_y=512 ";
params+= "density_of_objects=0.01 ";
params+= "gradslope_x=0.0005 avgc="+(c[i]/2)+" ";
params+= "norm_of_the_velocity=0.25 std_dev_of_the_norm_of_the_velocity=0.025 ";
params+= "rot_diff_coeff=0.001 run_to_tumble_rate_(step^-1)=0.06 tumble_to_run_rate_(step^-1)=0.06 tumbling=0.12 ";
params+= "surface_reorientation_rate=0.5 ";

//params+= "n=10 k=750 m=0.20 sigmak=1.50 sigmam=0.50 ";  // heterogeneous
params+= "n=0.8 k=33 m=0.37 sigmak=0 sigmam=0 ";  // homogeneous

params+= "output_file_data="+dir+"Trajectories_SimuChemotaxisGC2D_"+name+".data ";
params+= "log_file="+dir+"Log_SimuChemotaxisGC2D_"+name+".txt ";
//params+= "moviename=SimuChemotaxisGC2D_"+name+" l=2 e=1 h=25 dh=15 movietype=8-bit";

run("Simulation Chemotaxis R1 2D GrandCanonical", params);

makeRectangle(768, 0, 512, 512);
run("Crop");
run("Save", "save="+dir+"SimuChemotaxisGC2D_"+name+".tif");
close();

}
}
