/*
 * main.cpp
 *
 *  Created on: 17-Dec-2008
 *      Author: hbenamor
 */

#include "WilsonCowan.h"
#include "Pendule2.h"
#include "EnzymePFK.h"
#include "EnzymePFK2.h"
#include "Lorenz.h"
#include "VanDerPol2.h"
#include <pthread.h>
#include "VariationPendule.h"
#include "Rabinovich.h"
#include "Pendule2Oscillators.h"
#include "Pendule2Test.h"
#include "TestSampleDS.h"
#include "YeastGlycosis.h"
#include "DSGrid.h"
#include <iostream>

#include <iostream>
#include <fstream>
#include <istream>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <cstring>
#include <cmath>
#include <ctime>

#include <sys/types.h>
#include <sys/stat.h>



using namespace std;

CyclicAttractor * cyclicAttractor;
DSGrid * dsGrid;
int nbPoints, numberOfIterations, gridSimulationSelector, dynamicalSystem, paramIso[5], simPortraitParam, gridWidth, gridHeight, gridNbCompute, gridSteps, simGridNbIterations, indexOfParameterToIterate;
char name[255], iterationName[255];
double radius, deltaT, paramPortrait, simGridStimulationStrength, parameterMinValue, parameterMaxValue, parameterStepValue;
double myParam[255], simPPhiParam[3], simPXYParam[4], simPParam[2], coorOneTrajectory[255], window[255];
bool dsParameters, simIteration, simGrid, simIso, simTrajectory, simAttractor, simPortrait, simPXY, simPPhi, simP, simGridPerturbation, windowPulling, exitForced;
bool simPerturbationResponse;

void initDefault(){
	deltaT = 0.001; // Default Value
	numberOfIterations = 10000;
	dynamicalSystem = -1;
	//name = "CyclicAttractor";
	paramIso[0] = 30;		// Nb Isochrons
	paramIso[1] = 10000;	// Nb points per Isochrons
	paramIso[2] = 10;		// Nb IsoK
	paramIso[3] = 10000;	// Nb points per Isoks
	paramIso[4] = 0;		// Boolean generateIsoks?
	paramIso[5] = 10;		// Precision parameter

	window[0] = -10;		// X1
	window[1] = 10;			// X2
	window[2] = -10;		// Y1
	window[3] = 10;			// Y2
	window[4] = -10;		// Z1
	window[5] = 10;			// Z2

	gridSimulationSelector = 1;  // Id of the simulation to perform
	gridWidth = 100;		// Width of the grid of oscillators
	gridHeight = 100;		// Height of the grid of oscillators
	gridNbCompute = 1; 		// Number Of Computations to perform on a grid
	gridSteps = 10;			// Step of capturing images

	simGridNbIterations = 1000;		// NumberOfIterations for the Grid Max phase shift simulation
	simGridStimulationStrength = 10;	// Strength of the perturbation to apply to the Grid oscillators

	simPXYParam[0] = 20;	// Nb of points to extract and to test for the sensitivity
	simPXYParam[1] = 3;		// Rx max
	simPXYParam[2] = 100;	// Ry max
	simPXYParam[3] = 100;	// Nb of points per dimension (determin the step within consecutive perturbations)

	coorOneTrajectory[0] = 10;	// Sim Trajectory X
	coorOneTrajectory[0] = 10;	// Sim Trajectory Y
	coorOneTrajectory[0] = 0;	// Sim Trajectory Z

	simPParam[0] = 20;		// Nb of points to extract and to test for the sensitivity
	simPParam[1] = 10;		// RMax (radius of the circle)

	simPPhiParam[0] = 20;		// Nb of points to extract and to test for the sensitivity
	simPPhiParam[1] = 45;		// Phi (the direction in radian of the perturbations)
	simPPhiParam[2] = 20;		// Rmax (max of the perturbation strength)

	simPortraitParam = 10;		// Rayon du cercle pour tracer le portrait

	nbPoints = 100;
	radius = 7;

	//lorenz.tracerIsochrons(30, 1000000, 30, 1000000, false);
	simIso = false;			// Compute isochrons of the system
	windowPulling = false;	// Pull in a window space
	simPXY = false;		// Compute sensitivity to phase perturbation following a grid mapping for 2D systems
	simPPhi = false;		// Compute sensitivity to phase perturbation through a given direction for 2D systems
	simP = false;		// Compute sensitivity to phase perturbation through a circular are for 2D Systems
	simAttractor = false;	// Compute Attractor of the system
	simPortrait = false;	// Compute the phase portrait of the system
	simTrajectory = false;	// Compute one trajectory
	dsParameters = false;	// DS Parameters are by default not furnished otherwise true
	simPortrait = false;	// Compute the profile of the system
	simGrid = false;		// Simulate a grid of oscillators;
	simGridPerturbation = false; // Simulation of the evolution of the maximum phase shift of a population of oscillators
	simIteration = false;		// Iteration simulation
	exitForced = false;
	simPerturbationResponse = false; // Simulate the response to instant perturbations



}

/*
 * Help for the usage of this program.
 *
 */
void Help(bool ExitForced)
{
	cout << endl << "'-ds type name' --> specify the type of the dynamical system to study and the directory name of the output data files\n";
	cout << endl << "'-g simulationID(1) width(100) height(100) gridNbCompute(1) steps(10)' --> Construct a matrix of oscillators\n";
	cout << endl << "type: 1 - Anharmonic Pendulum, 2 - Wilson Cowan, 3 - Van der Pol, 4 - Lorenz, 5 - Rabinovich, 6 - PFK, 7 - Yeast Glyocosis\n";
	cout << endl << "'-dt deltaT(0.001)' -> the time step to use default value is 0.001\n";
	cout << endl << "'-m maxNumberOfIterations' -> Maximum number of iterations, the default value is 10 000\n";
	cout << endl << "'-p N parameter1 parameter2 ... paramterN' -> N is the numberOfParameters including the dimensions, firsts parameters are the initial conditions\n";
	cout << endl << "'-sIso NumberOfIsochrons(30) NumberOfPointsPerIsochron(10000) NumberOfIsoKs(10) NumberOfPointsPerIsoKs(10000) GenerateIsoKs(0:false) Precision(10)' -> compute the isochrons of the system\n";
	cout << endl << "'-sPXY NumberOfPointsToTest(20) RXMax(3) RYMax(100) NumberOfPointsPerDimension(100)' -> Compute the maximal phase shift following a grid for 2D systems\n";
	cout << endl << "'-sPPhi NumberOfPointsToTest(20) Phi(45°) RMax(100)' -> Compute the maximal phase shift following a direction (in degree) for 2D systems\n";
	cout << endl << "'-sP NumberOfPointsToTest(20) RMax(10)' -> Compute the maximal phase shift following a circular area for 2D systems\n";
	cout << endl << "'-sT N(3) coor1(10) coor2(10) coor3(0) ... coorN' -> compute one trajectory N is the dimension\n";
	cout << endl << "'-sPor Radius(10)' -> compute the profile of the dynamical system following a circle\n";
	cout << endl << "'-w NbParam(4) X1min(-10) X1max(10) X2min(-10) X2max(10)' -> restrict random shots on a window in the phase space\n";
	cout << endl << "'-gP simGridNbIterations(1000) simGridStimulationStrength(10) ' -> Trace the evolution of the maximum phase shift after a perturbation\n";
	cout << endl << "'-simIt nameOfSimulation indexOfParameterToIterate parameterMinValue parameterMaxValue parameterStepValue  ' -> Iterate the command line considering all the values of parameter giving a step\n";
	cout << endl << "'-simPR Radius NbPoints ' -> Return phase response of the oscillator\n";
	cout << endl << "'-h' or '--help' --> Print this help menu.\n";
	cout << endl << endl << "usage e.g.1 ./DSAnalyser -ds 1 Pendule -dt 0.01 -m 10000 -p 2 5 5 -sT 2 3 2 -sP 20 10\n";
	cout << endl << "usage e.g.2 ./DSAnalyser -ds 1 WilsonCowan -dt 0.01 -m 10000 -p 2 5 5 -sT 2 3 2 -sP 20 10 -sIso 30 100 10 100 0\n";
	cout << endl;
	if (ExitForced) exit(0);
}

//void CyclicAttractor::tracerInstantPerturbationResponse(double deltaT, double R, double NbPoints){
/*-------------------------------------------------------------------------*/
/*------------------- Lecture des arguments d'entree  ---------------------*/
/*-------------------------------------------------------------------------*/

void ReadArgs(int argc, char *argv[]){
for (int j = 1; j < argc; j++){
    if (argv[j][0] == '-')
        {
        if(strcmp(argv[j],"-ds") == 0 || strcmp(argv[j],"--dynamicalsystem") == 0)   {
        	dynamicalSystem = atoi(argv[j+1]);
        	sprintf(name,"%s", argv[j+2]);
        	j=j+2;
       	}
        else if(strcmp(argv[j],"-simPR") == 0 || strcmp(argv[j],"--PerturbationResponse") == 0)   {
           	radius = atof(argv[j+1]);
           	nbPoints = atoi(argv[j+2]);
           	simPerturbationResponse = true;
           	simAttractor = true;
           	j=j+2;
        }
        else if(strcmp(argv[j],"-simIt") == 0 || strcmp(argv[j],"--simIteration") == 0)   {
			sprintf(iterationName, "%s",argv[j+1]);
			indexOfParameterToIterate = atoi(argv[j+2]);
			parameterMinValue = atof(argv[j+3]);
			parameterMaxValue = atof(argv[j+4]);
			parameterStepValue = atof(argv[j+5]);
			simAttractor = true;
			simIteration = true;
			j=j+5;
        }
        else if(strcmp(argv[j],"-g") == 0 || strcmp(argv[j],"--grid") == 0)   {
        	gridSimulationSelector = atoi(argv[j+1]);
        	gridWidth = atoi(argv[j+2]);
			gridHeight = atoi(argv[j+3]);
			gridNbCompute = atoi(argv[j+4]);
			gridSteps = atoi(argv[j+5]);
			simGrid = true;
			simAttractor = true;
			j=j+5;
        }
        else if(strcmp(argv[j],"-gP") == 0 || strcmp(argv[j],"--gridPerturbation") == 0)   {
        	simGridNbIterations = atoi(argv[j+1]);
        	simGridStimulationStrength = atof(argv[j+2]);
			simGridPerturbation = true;
			simAttractor = true;
			j=j+2;
        }
        else if(strcmp(argv[j],"-dt") == 0 || strcmp(argv[j],"--deltaT") == 0)   {
        	deltaT = atof(argv[j+1]);
        	j++;
        }
        else if(strcmp(argv[j],"-m") == 0 || strcmp(argv[j],"--maxNumberOfIterations") == 0)   {
			numberOfIterations = atoi(argv[j+1]);
			j++;
		}
        else if(strcmp(argv[j],"-p") == 0 || strcmp(argv[j],"--parameters") == 0)   {
        	int nbParam = atoi(argv[j+1]);
        	dsParameters = true;
        	// L'element qui suit est le nombre de parametres dans la liste
        	for(int i = 0; i < nbParam;i++) myParam[i] = atof(argv[j+2+i]);
        	j = j+nbParam+1;
        }
        else if(strcmp(argv[j],"-sIso") == 0 || strcmp(argv[j],"--simIsochrons") == 0)   {
        	simIso = true;
        	simAttractor = true;
        	for(int i=0; i < 6;i++) paramIso[i] = atoi(argv[j+i+1]);
        	j=j+6;
        }
        else if(strcmp(argv[j],"-sPXY") == 0 || strcmp(argv[j],"--simPXY") == 0)   {
               	simPXY = true;
               	simAttractor = true;
               	for(int i=0; i < 4;i++) simPXYParam[i] = atoi(argv[j+i+1]);
            	j=j+4;
        }
        else if(strcmp(argv[j],"-sPPhi") == 0 || strcmp(argv[j],"--simPPhi") == 0)   {
				simPPhi = true;
				simAttractor = true;
				for(int i=0; i < 3;i++) simPPhiParam[i] = atof(argv[j+i+1]);
				j=j+3;
         }
        else if(strcmp(argv[j],"-sP") == 0 || strcmp(argv[j],"--simP") == 0)   {
				simP = true;
				simAttractor = true;
				for(int i=0; i < 2;i++) simPParam[i] = atoi(argv[j+i+1]);
				j=j+2;
        }
        else if(strcmp(argv[j],"-w") == 0 || strcmp(argv[j],"--window") == 0)   {
        	int nbParam = atof(argv[j+1]);
        	windowPulling = true;
			// L'element qui suit est le nombre de parametres dans la liste
			for(int i = 0; i < nbParam;i++) window[i] = atof(argv[j+2+i]);
			for(int i=0; i < nbParam/2; i++) if(window[i*2]>window[i*2+1]){
				printf("\n Option window wrong: make sure that the interval are correct\n");
				printf("The interval -w 2 x1 x2 is valid if x1 <= x2\n");
				printf("Usage e.g: -w 4 -10 10 -3 -2 \n");
			}
			j = j+nbParam+1;
        }
        else if(strcmp(argv[j],"-sT") == 0 || strcmp(argv[j],"--simT") == 0)   {
				simTrajectory = true;
				simAttractor = true;
				int nbParam = atoi(argv[j+1]);
				for(int i=0; i < nbParam;i++) coorOneTrajectory[i] = atof(argv[j+i+2]);
				j=j+nbParam+1;
        }
        else if(strcmp(argv[j],"-sPor") == 0 || strcmp(argv[j],"--simPortrait") == 0)   {
				simPortrait = true;
				simPortraitParam = atoi(argv[j+1]);
				j++;
         }
		else if(strcmp(argv[j],"-h") == 0 || strcmp(argv[j],"--help") == 0)   Help(exitForced);
        }
    }
}

//Debug/DSAnalyser -ds 3 VANFirst -p 3 1 1 5 -sIso 30 10 10 10 0 1000  -dt 0.01 -m 50000

void runSimulations(char * name){
	//CyclicAttractor * cyclicAttractor;
	printf("here1\n");
	if(dynamicalSystem == 1) {
		if(dsParameters) cyclicAttractor = new Pendule2(myParam,name);
		else cyclicAttractor = new Pendule2(name);
	}
	else if(dynamicalSystem == 2) {
			if(dsParameters) cyclicAttractor = new WilsonCowan(myParam,name);
			else cyclicAttractor = new WilsonCowan(name);
	}
	else if(dynamicalSystem == 3){
		if(dsParameters) cyclicAttractor = new VanDerPol2(myParam,name);
		else cyclicAttractor = new VanDerPol2(name);
	}
	else if(dynamicalSystem == 4) {
		if(dsParameters)	cyclicAttractor = new Lorenz(myParam,name);
		else 	cyclicAttractor = new Lorenz(name);
	}
	else if(dynamicalSystem == 5) {
		if(dsParameters) cyclicAttractor = new Rabinovich(myParam,name);
		else cyclicAttractor = new Rabinovich(name);
	}
	else if(dynamicalSystem == 6) {
		if(dsParameters) cyclicAttractor = new EnzymePFK2(myParam,name);
		else cyclicAttractor = new EnzymePFK2(name);
	}
	else if(dynamicalSystem == 7) {
			if(dsParameters) cyclicAttractor = new YeastGlycosis(myParam,name);
			else cyclicAttractor = new YeastGlycosis(name);
	}
	else Help(true);

	printf("here 2 - Instanced\n");

	if(simAttractor) cyclicAttractor->tracerAttractor(deltaT, numberOfIterations);
	if(simGrid and simAttractor){
		dsGrid = new DSGrid(gridHeight, gridWidth, cyclicAttractor);
		if(simGridPerturbation)
			dsGrid->traceMaximumPhaseShift(simGridStimulationStrength,simGridNbIterations);
		else{
			if(gridSimulationSelector==1) dsGrid->computeGrid(gridNbCompute, gridSteps);
			else dsGrid->computeAnimation(gridNbCompute, gridSteps);
		}
		//dsGrid->compute(gridNbCompute);
		//dsGrid->saveStates();
	}
	if(simPortrait and not simAttractor)
		cyclicAttractor->tracerPortraitWithoutCalculusOfAttractor(deltaT, numberOfIterations, simPortraitParam);
	else if(simAttractor and simPortrait)
		cyclicAttractor->tracerPortrait(deltaT, numberOfIterations, simPortraitParam);
	if(simTrajectory) cyclicAttractor->traceOneTrajectory(coorOneTrajectory, deltaT, numberOfIterations);
	if(simPPhi)		cyclicAttractor->tracerPhaseSensitivityToPerturbationAccordingToPhi((int)simPPhiParam[0],(double)M_PI*(simPPhiParam[1]/180),(double)simPPhiParam[2], false);
	if(simP)		cyclicAttractor->tracerPhaseSensitivityToPerturbation((int)simPParam[0],(int)simPParam[1]);
	if(simPXY)		cyclicAttractor->tracerPhaseSensitivityToPerturbationXY((int)simPXYParam[0],(double)simPXYParam[1],(double)simPXYParam[2], (int)simPXYParam[3]);
	if(simPerturbationResponse) cyclicAttractor->tracerInstantPerturbationResponse(deltaT, radius,nbPoints);
	if(simIso)		{
		if(windowPulling)
			cyclicAttractor->tracerIsochronsInAWindow(paramIso[0],paramIso[1],paramIso[2],paramIso[3],(bool)paramIso[4], paramIso[5],window);
		else
			cyclicAttractor->tracerIsochrons(paramIso[0],paramIso[1],paramIso[2],paramIso[3],(bool)paramIso[4], paramIso[5]);
	}
	printf("here 3 - Simulation Finished\n");
	delete cyclicAttractor;
}


int main(int argc, char *argv[]){
	try{
		char simulationFile[255];
		time_t rawtime;
		time(&rawtime);

		initDefault();
		ReadArgs(argc, argv);

		char rootPathName[255];
		char tmpSimName[255];

		//if(!simAttractor) Help(true);
		if(simIteration and dsParameters){
			sprintf(rootPathName,"./%s", iterationName);
			int status;
			status = mkdir(rootPathName, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
			if(status!=0) printf("Directory %s already exists\n", rootPathName);
			sprintf(rootPathName,"%s", name);
			// Generate Simulation Infor Header
			sprintf(simulationFile, "./%s/SimulationInfo.txt", iterationName);
			ofstream out(simulationFile, ios::app);
			out << "###DSAnalyser simulations information file###\n";
			out << "Simulation started by: \n";
			for(int i=0; i < argc; i++) out << argv[i] << " ";
			out <<"\n";
			out << "Begin: " << ctime(&rawtime);
			out.close();
			// End of simulation info file header generation
			for(double selectedParameter=parameterMinValue; selectedParameter <= parameterMaxValue; selectedParameter+=parameterStepValue){
				sprintf(tmpSimName,"%s/%s_%f", iterationName, name,selectedParameter);
				myParam[indexOfParameterToIterate] = selectedParameter;
				runSimulations(tmpSimName);
			}
		}
		else{
			sprintf(rootPathName,"%s", name);
			sprintf(simulationFile, "./%s/SimulationInfo.txt", rootPathName);
			ofstream out(simulationFile, ios::app);
			out << "###DSAnalyser simulations information file###\n";
			out << "Simulation started by: \n";
			for(int i=0; i < argc; i++) out << argv[i] << " ";
			out <<"\n";
			out << "Begin: " << ctime(&rawtime);
			out.close();
			runSimulations(name);
		}


		ofstream out1(simulationFile, ios::app);
		time(&rawtime);
		out1 << "Finished: " << ctime(&rawtime) << endl << endl;
		out1.close();
		exit(0);

	} catch(...){
		Help(true);
	}
//	for (int j = 0; j<1000; j++)
	//	cout << "mmmmmmh I'm a so dirty boy !" << endl;

	/*Rabinovich rabinovich;
	rabinovich.tracerPortrait(0.001,10000,50);
*/
/*
	TestSampleDS testDS;
	testDS.tracerPortrait(0.001,10000,10);
	*/

	//double coor[3];
	/*coor[0] = -40.005;
	coor[1] = -107.719;
	coor[2] = -234.033;*/

	/*coor[0] = -35.0576;
	coor[1] = -829.588;
	coor[2] = 770.799;*/
/*
	coor[0] = -64.7064;
	coor[1] = 4.07577;
	coor[2] = -200.277;
*/
/*
	Lorenz lorenz;
	lorenz.tracerAttractor(0.001, 10000000);
	*/
	//lorenz.tracerPortrait(0.001,100000,10);
	//lorenz.traceOneTrajectory(coor,0.001,100000);//2238);
	//lorenz.tracerAttractor(0.001, 1000000);
	//lorenz.tracerIsochrons(30, 1000000, 30, 1000000, false);


/*
	VanDerPol2 vanderpol2;
	vanderpol2.tracerAttractor(0.001, 1000000);
	vanderpol2.tracerPortrait(0.01, 10000000, 10);
	vanderpol2.tracerIsochrons(30, 10000, 30, 10000, false);
*/
/*
		VanDerPol2 vanderpol2;
		vanderpol2.tracerAttractor(0.001, 1000000);
		//vanderpol2.deltaT = 0.001;
		//vanderpol2.tracerPortrait(0.01, 10000000, 10);
		vanderpol2.tracerIsochrons(30, 10000, 30, 10000, false);
*/
//	vanderpol2.gridExploration(0.001,100,100);

/*
	WilsonCowan wilsonCowan;
	wilsonCowan.tracerAttractor(0.01, 100000);
	//wilsonCowan.tracerPhaseSensitivityToPerturbationAccordingToPhi(0.001, 100000, 100, 0, 30); // Duree prevu 12h de calcul commence a 21h13

	wilsonCowan.tracerIsochrons(30, 10000,10,10000,false);
*/
/*
	Pendule2Test pendule2Test;
	pendule2Test.tracerPortraitPhase(0.01,10000,4);
	*/
/*
	Pendule2 pendule2;
	pendule2.tracerPortrait(0.001,2000000,5);
	pendule2.tracerAttractor(0.001, 1000000);
	pendule2.tracerIsochrons(30, 10000, 10, 10000, false);
	pendule2.tracerPortrait(0.01,20000,1);
*/

/*	Pendule2Oscillators pendule2Oscillators;
	//pendule2Oscillators.tracerPhase(0.01,0,10000);
	//pendule2Oscillators.tracerPhaseLag(0.01,0,10000);
	pendule2Oscillators.tracerPortrait(0.01,7000,1);*/
/*
	EnzymePFK enzymePFK;
	enzymePFK.tracer(0.001, 1000000);
	enzymePFK.tracerIsochrons(30, 10000);
*/

/*
	EnzymePFK2 enzymePFK2;
	enzymePFK2.tracerAttractor(0.001, 10000); //1000000 Delta T =0.001 est la meilleur valeur pour la PFK
	//enzymePFK2.tracerIsochrons(30, 10000, 10, 10000,false);
	//enzymePFK2.numberOfIterations = 10000000;
	//enzymePFK2.tracerPhaseSensitivityToPerturbation(20, 5);
	enzymePFK2.tracerPhaseSensitivityToPerturbationXY(20, 3, 100, 100);
*/
/*
	VariationPendule variationP;
	//variationP.tracerAttractor(0.00005, 100000;
	variationP.tracerPortrait(0.001,100,25);
	//lorenz.tracerAttractor(0.001, 1000000);
	//variationP.tracerIsochrons(30, 1000000);
*/
	return 0;
}
