dir = "E:\\Remy\\Simu_Manika\\2024-06-12_SimuGC_Pheno2D\\heterogeneousPop\\";

c = newArray(0,100,1000);


for(i=0; i<c.length; i++){

for(k=0; k<9; k++){

name = "R1-c-"+c[i]+"_"+k;

params ="datafile="+dir+"Trajectories_SimuChemotaxisGC2D_"+name+".data ";
params+="savefile="+dir+"runDurHistogram_Trajectories_SimuChemotaxisGC2D_"+name+".txt";
params+=" xmin=650 xmax=1400 tmin=2000 tmax=4999";
run("Analysis SimuGC 2D runDurHisto", params);

}}

dir = "E:\\Remy\\Simu_Manika\\2024-06-12_SimuGC_Pheno2D\\homogeneousPop\\";

for(i=0; i<c.length; i++){

for(k=0; k<9; k++){

name = "R1-c-"+c[i]+"_"+k;

params ="datafile="+dir+"Trajectories_SimuChemotaxisGC2D_"+name+".data ";
params+="savefile="+dir+"runDurHistogram_Trajectories_SimuChemotaxisGC2D_"+name+".txt";
params+=" xmin=650 xmax=1400 tmin=2000 tmax=4999";
run("Analysis SimuGC 2D runDurHisto", params);

}}
