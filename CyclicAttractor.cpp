/*
 * CyclicAttractor.cpp
 *
 *  Created on: 19-Jan-2009
 *      Author: hbenamor
 */

//#include "CyclicAttractor.h"
#include "CyclicAttractorSmartPencil.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
#include "random-singleton.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define PI 3.14159265

using namespace std;

CyclicAttractor::CyclicAttractor() {

}

CyclicAttractor::CyclicAttractor(int nbDimensions, int nbParam, double * myParam, char * name){
	initializePath(name);
	this->nbDimension = nbDimensions;
	this->nbParam = nbParam;
	for(int i=0;i < nbDimensions; i++) {
		coordinates[i] = myParam[i];
		minFirstCoor[i] = 0;
		maxFirstCoor[i] = 0;
	}
	for(int i=0;i < nbParam; i++) param[i] = myParam[i+nbDimensions];
}

CyclicAttractor::CyclicAttractor(char * name) {
	// TODO Auto-generated constructor stub
	initializePath(name);
}

CyclicAttractor::CyclicAttractor(const CyclicAttractor & cyclicAttractor){
	sprintf(this->name, "%s", cyclicAttractor.name);
	this->nbDimension = cyclicAttractor.nbDimension;
	this->nbParam= cyclicAttractor.nbParam;
	this->periode = cyclicAttractor.periode;
	this->deltaT = cyclicAttractor.deltaT;
	this->hypothenus = cyclicAttractor.hypothenus;
	this->nbPointsPerPeriode = cyclicAttractor.nbPointsPerPeriode;
	this->nbPointsPerIsoKs = cyclicAttractor.nbPointsPerIsoKs;
	this->nbPointsPerIsochrons = cyclicAttractor.nbPointsPerIsochrons;
	this->isochronsFitness = cyclicAttractor.isochronsFitness;
	this->nbIsoKs = cyclicAttractor.nbIsoKs;
	this->nbIsochrons = cyclicAttractor.nbIsochrons;
	this->numberOfIterations = cyclicAttractor.numberOfIterations;
	this->resolution = cyclicAttractor.resolution;
	for(int i=0;i < nbDimension; i++) {
		this->coordinates[i] = cyclicAttractor.coordinates[i];
		this->radius[i] = cyclicAttractor.radius[i];
	}
	for(int i=0;i < nbParam; i++) this->param[i] = cyclicAttractor.param[i];
	for(int i=0; i < (int)cyclicAttractor.attractor.size();i++){
		this->attractor.push_back(new point3D(*(cyclicAttractor.attractor[i])));
	}

}

void CyclicAttractor::initializePath(char * name){
	sprintf(this->name, "%s", name);
	sprintf(rootPath, "./%s", name);
	sprintf(isochronsPath, "./%s/Isochrons", name);
	sprintf(trajectoriesPath, "./%s/Trajectories", name);
	int status;
	status = mkdir(this->rootPath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if(status!=0) printf("Directory %s already exists\n", this->rootPath);
	status = mkdir(this->isochronsPath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if(status!=0) printf("Directory %s already exists\n", this->isochronsPath);
	status = mkdir(this->trajectoriesPath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if(status!=0) printf("Directory %s already exists\n", this->trajectoriesPath);
	srand(time(NULL));
}

CyclicAttractor::~CyclicAttractor() {
	// TODO Auto-generated destructor stub
	/*
	while((int)attractor.size()>0){
		delete attractor[attractor.size()-1];
		attractor.pop_back();
	}
	while((int)isochrons.size()>0){
			delete isochrons[isochrons.size()-1];
			isochrons.pop_back();
	}
	while((int)isoKs.size()>0){
			delete isoKs[isoKs.size()-1];
			isoKs.pop_back();
	}
	while((int)smartPencils.size()>0){
			delete smartPencils[smartPencils.size()-1];
			smartPencils.pop_back();
	}
	while((int)myThreads.size()>0){
			pthread_cancel(*myThreads[myThreads.size()-1]);
			delete myThreads[myThreads.size()-1];
			myThreads.pop_back();
	}
	while((int)mutexs.size()>0){
			delete mutexs[mutexs.size()-1];
			mutexs.pop_back();
	}
	while((int)isoKsMutexs.size()>0){
			delete isoKsMutexs[isoKsMutexs.size()-1];
			isoKsMutexs.pop_back();
	}
	*/
}

/*
 * Trace a Trajectory begining from the point coor and stocke it
 * to the file DSTest_TrajectoryOne.dat"
 */
void CyclicAttractor::traceOneTrajectory(double * coor, double deltaT, int numberOfIteration){
	point3D point(0,0,0,0,0,0,0);
	char outPutTrajectoryFile[255];
	sprintf(outPutTrajectoryFile, "./%s/DSTest_TrajectoryOne.dat", rootPath);
	ofstream out(outPutTrajectoryFile);
	int j=0;
	double phase = 0;
	MTRand myRandom;

	point.x = coor[0];
	point.y = coor[1];
	point.z = coor[2];
	phase = phaseOf(&point,0);
	out << point.x << "\t" << point.y << "\t" << point.z << "\t" << point.angle << "\t" << point.norme << "\t" << point.phase << "\t" << point.k << "\n";
	while(j < numberOfIteration){
		computeMethod(coor, 1);
		point.x = coor[0];
		point.y = coor[1];
		point.z = coor[2];
		phase = phaseOf(&point,0);
		//if(myRandom.rand() > 0.92)
		out << point.x << "\t" << point.y << "\t" << point.z << "\t" << point.angle << "\t" << point.norme << "\t" << point.phase << "\t" << point.k << "\n";
		j++;
	}
	out.close();
}

/*
 * Trace le portrait de phase du système dynamique selon le shéma suivant
 * en 2D tracer les trajectoires à partir de points situés sur le cercle
 * de rayon l'hypothenus (selon le tableau radius)
 * Puis faire de même pour le cercle de rayon 1000 fois plus petit situé
 * à l'intérieur du cycle. Les informations relatifs aux points des trajectoire
 * tel que la norme, l'angle et la vitesse, la phase et le temps de convergence
 * vers l'attracteur sont assurées par l'appel à la méthode phaseOf(&point, 0);
 *
 */
void CyclicAttractor::tracerPortrait(double deltaT, int numberOfIteration, int nbShots){

	char trajecFilePath[255];
	double hypothenus = 0.0, rho = 0.0, phi=0.0, theta=0.0;
	this->deltaT = deltaT;
	double angularStep = PI / nbShots;
	double phase;
	int j=0;
	point3D point(0,0,0,0,0,0,0);

	for(int i=0; i < nbDimension; i++){
		hypothenus += radius[i] * radius[i];
	}
	hypothenus = sqrt(hypothenus);
	rho = hypothenus;

	for(int m=0; m < 2;m++){
		if(m==1) rho = rho/1000;
		for(int i=0; i < nbShots; i++){
			if(m==1)
				sprintf(trajecFilePath, "%s/DSTest_Trajectory_%d.dat", trajectoriesPath, i+nbShots);
			else
				sprintf(trajecFilePath, "%s/DSTest_Trajectory_%d.dat", trajectoriesPath, i);
			ofstream out(trajecFilePath);
			coordinates[0] = rho * cos(phi);
			coordinates[1] = rho * sin(phi) * cos(theta);
			coordinates[2] = rho * sin(phi) * sin(theta);
			computeMethod(coordinates,0); // Pour le pendule anharmonique
			point.x = coordinates[0];
			point.y = coordinates[1];
			point.z = coordinates[2];
			phase = phaseOf(&point,0);
			out << point.x << "\t" << point.y << "\t" << point.z << "\t" << point.angle << "\t" << point.norme << "\t" << point.phase << "\t" << point.k << "\n";
			while(j < numberOfIteration and point.k!=1){
				computeMethod(coordinates, 1);
				point.x = coordinates[0];
				point.y = coordinates[1];
				point.z = coordinates[2];
				phase = phaseOf(&point,0);
				out << point.x << "\t" << point.y << "\t" << point.z << "\t" << point.angle << "\t" << point.norme << "\t" << point.phase << "\t" << point.k << "\n";
				j++;
			}
			phi += angularStep;
			theta += angularStep;
			out.close();
		}
	}
}

/*
 * Without Extracting the phase on the attractor
 */
void CyclicAttractor::tracerPortraitWithoutCalculusOfAttractor(double deltaT, int numberOfIteration, int nbShots){

	char trajecFilePath[255];
	double hypothenus = 0.0, rho = 100;
	this->deltaT = deltaT;
	double angularStep = 2*PI / nbShots;
	double phase;
	int j=0;
	point3D point(0,0,0,0,0,0,0);

	rho = 10;

	int count=0;
	for(int m=0; m < 2;m++){
		if(m==1) rho = rho/1000;
		for(double phi=0; phi < 2*M_PI; phi+=angularStep){
			for(double theta=0; theta < M_PI; theta+=(angularStep/2)){
				sprintf(trajecFilePath, "%s/DSTest_Trajectories_%d.dat", trajectoriesPath,count);
				count++;
				ofstream out(trajecFilePath);
				coordinates[0] = rho * cos(phi);
				coordinates[1] = rho * sin(phi) * cos(theta);
				coordinates[2] = rho * sin(phi) * sin(theta);
				computeMethod(coordinates,0); // Pour le pendule anharmonique
				point.x = coordinates[0];
				point.y = coordinates[1];
				point.z = coordinates[2];
				if(not finite(point.x) or not finite(point.y) or not finite(point.z)){
					out.close();
					exit;
				}
				out << point.x << "\t" << point.y << "\t" << point.z << "\n";
				while(j < numberOfIteration and point.k!=1){
					computeMethod(coordinates, 1);
					point.x = coordinates[0];
					point.y = coordinates[1];
					point.z = coordinates[2];
					if(not finite(point.x) or not finite(point.y) or not finite(point.z)){
						out.close();
						exit;
					}
					out << point.x << "\t" << point.y << "\t" << point.z << "\n";
					j++;
				}
				j=0;
				out.close();
			}
		}
	}
	char tmpString[255];
	sprintf(trajecFilePath, "%s/DSTest_Trajectories.gp", trajectoriesPath);
	ofstream out(trajecFilePath);
	out << "set key off\n";
	out << "splot \"DSTest_Trajectories_0.dat\" w l";
	for(int i=0;i<count;i++){
		sprintf(tmpString, "DSTest_Trajectories_%d.dat", i);
		out << ", \"" << tmpString << "\" w l";
	}
	out << "\n pause -1\n";

	out.close();
}


void CyclicAttractor::tracerPortraitPhase(double deltaT, int numberOfIteration, int nbShots){

	char trajecFilePath[255];
	double hypothenus = 0.0, rho = 0.0, phi=0.0;
	this->deltaT = deltaT;
	double angularStep = 2*PI / nbShots;
	double phase;
	int j=0;
	point3D point(0,0,0,0,0,0,0);

	hypothenus = 20;
	hypothenus = sqrt(hypothenus);
	rho = hypothenus;
	sprintf(trajecFilePath, "%s/DSTest_Portrait.dat", trajectoriesPath);
	ofstream out(trajecFilePath);
	for(int m=0; m < 2;m++){
		if(m==1) {
			rho = rho/1000;
			phi = 0;
		}
		for(int i=0; i < nbShots; i++){
			coordinates[0] = rho * cos(phi);
			coordinates[1] = rho * sin(phi);
			coordinates[2] = 0;
			//computeMethod(coordinates,0); // Pour le pendule anharmonique
			point.x = coordinates[0];
			point.y = coordinates[1];
			point.z = coordinates[2];
			//phase = phaseOf(&point,0);
			if(not finite(point.x) or not finite(point.y) or not finite(point.z)){
				out.close();
				exit;
			}
			out << point.x << " " << point.y << " " << point.z << "\n";
			//out << point.x << "\t" << point.y << "\t" << point.z << "\t" << point.angle << "\t" << point.norme << "\t" << point.phase << "\t" << point.k << "\n";
			j = 0;
			while(j < numberOfIteration){
				computeMethod(coordinates, 1);
				point.x = coordinates[0];
				point.y = coordinates[1];
				point.z = coordinates[2];
				if(not finite(point.x) or not finite(point.y) or not finite(point.z)){
					out.close();
					exit;
				}
				//phase = phaseOf(&point,0);
				out << point.x << " " << point.y << " " << point.z << "\n";
				//out << point.x << "\t" << point.y << "\t" << point.z << "\t" << point.angle << "\t" << point.norme << "\t" << point.phase << "\t" << point.k << "\n";
				j++;
			}
			out << "\n\n";
			phi += angularStep;
		}
	}
	out.close();
}

/*
 * Analyse the response
 */
void CyclicAttractor::tracerInstantPerturbationResponse(double deltaT, double R, double NbPoints){
	char perturbationResponse[255];
	sprintf(perturbationResponse, "%s/perturbationResponse.dat", rootPath);
	ofstream out(perturbationResponse);
	this->deltaT = deltaT;
	this->numberOfIterations = numberOfIterations;
	double tmpState[3];
	double phase = 0;

	for(double angle =0; angle <= 2 *M_PI; angle += (2*PI)/NbPoints){
		tmpState[0]= R*cos(angle);
		tmpState[1] = R* sin(angle);
		tmpState[2] = 0;
		phase = (this->phaseOf(tmpState)/periode)*(2*M_PI);
		out << angle << " " << phase << "\n";
	}
}

void CyclicAttractor::tracerPhaseLag(double deltaT, double initial, int numberOfIterations){
	char phaseLagPath[255];
	sprintf(phaseLagPath, "%s/phaseLag.dat", rootPath);
	ofstream out(phaseLagPath);
	this->deltaT = deltaT;
	this->numberOfIterations = numberOfIterations;
	double phaseLag=0;
	double oldPhaseLag=0;

	for(int i=0; i < numberOfIterations; i++){
		oldPhaseLag = -atan2(coordinates[1], coordinates[0]) +atan2(coordinates[3], coordinates[2]);
		computeMethod(coordinates,1);
		phaseLag = -atan2(coordinates[1], coordinates[0]) +atan2(coordinates[3], coordinates[2]);
		if(fabs(oldPhaseLag-phaseLag)>3)	out << "\n";
		//else out << oldPhaseLag << " " << phaseLag << "\n";
		else{
			if(fabs(phaseLag) >5) phaseLag += 2*PI;
			out << i*deltaT << " " << phaseLag << "\n";
		}
	}
}
void CyclicAttractor::tracerPhaseSensitivityToPerturbation(int nbPointsToTest, double Rmax){
	vector<point3D> initialPoints;
	extractPoints(nbPointsToTest, &initialPoints);
	for(double phi=0; phi <2*PI; phi+=0.1)
		tracerPhaseSensitivityToPerturbationAccordingToPhi(&initialPoints, phi, Rmax, false);

	char sensitivityGnuplotPath[255];
	sprintf(sensitivityGnuplotPath, "%s/DSTest_Perturbation_Phi.gp", rootPath);
	ofstream outMe(sensitivityGnuplotPath);

	outMe << "set title \"Maximum Phase Shift of the "<< name <<" Attractor\"\n";
	outMe << "set xlabel \"X\"\n";
	outMe << "set xrange["<<-Rmax<<":"<<Rmax<<"]\n";
	outMe << "set ylabel \"Y\"\n";
	outMe << "set yrange["<<-Rmax<<":"<<Rmax<<"]\n";
	outMe << "set zrange[0:"<< 2*M_PI << "]\n";
	outMe << "set zlabel \"Maximum Phase Shift\"\n";
	outMe << "unset key\n";
	outMe << "splot \"DSTest_Perturbation_Phi.dat\" title 'Maximum Phase Shift of the "<< name << " Attractor' with pm3d\n";
	outMe << "pause -1\n";
	outMe.close();
}

/*
 * Choose the extracting methods
 */
void CyclicAttractor::extractPoints(int numberOfPoints, vector<point3D>* initialPoints){
	extractPointsFromTheAttractor(numberOfPoints, initialPoints);
	//generateRandomPointsOnTheAttractor(numberOfPoints, initialPoints);
}

/*
 * Extract equally spaced point from the already found attractor
 */
void CyclicAttractor::extractPointsFromTheAttractor(int numberOfPoints, vector<point3D>* initialPoints){
	int step= (int)attractor.size()/numberOfPoints;
	(*initialPoints).clear();
	for(int j=0; j< (int)attractor.size();j+=step)
			(*initialPoints).push_back(point3D(attractor[j]->x,attractor[j]->y,0));

}

/*
 * Generate random points on the attractor
 */
void CyclicAttractor::generateRandomPointsOnTheAttractor(int numberOfPoints, vector<point3D>* initialPoints){
	MTRand myRandom;
	double tmpPoint[255];
	(*initialPoints).clear();
	for(int i=0; i< numberOfPoints;i++){
			tmpPoint[0] = myRandom.rand() * 2 * 10 - 10;
			tmpPoint[1] = myRandom.rand() * 2 * 10 - 10;
			computeMethod(tmpPoint,numberOfIterations);
			(*initialPoints).push_back(point3D(tmpPoint[0],tmpPoint[1],0));
	}
}


/*
 *
 */
void CyclicAttractor::tracerPhaseSensitivityToPerturbationXY(int nbPointsToTest, double Rxmax, double Rymax, int nbPointsPerDimension){
	char perturbationPath[255];
	vector<point3D> initialPoints;
	extractPoints(20, &initialPoints);
	sprintf(perturbationPath, "%s/DSTest_PerturbationXY.dat", rootPath);

	double phase;
	double numPoints=nbPointsPerDimension; //100
	double phasePI;
	double maxDistances = 0;
	double sumDistances = 0;
	bool linePrintedBefore = false;
	int numberOfPoints = (int)initialPoints.size();
	vector<point3D> vectPoint3D;
	vector<double> distances;
	MTRand myRandom;

	int n=10000;
	int partitionPhase[10000];
	for(int i=0;i<n;i++) partitionPhase[i] = 0;
	int tmpIndices=0;

	for(int i=0; i< numberOfPoints;i++){
		distances.push_back(0);
		vectPoint3D.push_back(point3D(0,0,0));
	}

	for(double Rx = -Rxmax; Rx  <= Rxmax; Rx +=(Rxmax/numPoints)){
		for(double Ry = -Rymax; Ry <= Rymax; Ry += (Rymax/numPoints)){
			maxDistances = -1;
			sumDistances = 0;
			ofstream out(perturbationPath, ios::app);
			for(int i=0; i < numberOfPoints; i++){
				vectPoint3D[i].x = initialPoints[i].x + Rx;
				vectPoint3D[i].y = initialPoints[i].y + Ry;
				phase = phaseOf(&vectPoint3D[i],0);
				phasePI = (2*PI*phase/periode);
				while(phasePI >= 2*PI) phasePI -= 2*PI;
				if(phase==-1)
					vectPoint3D[i].phase= -1;
				else
					vectPoint3D[i].phase= phasePI;
			}
			sort(vectPoint3D.begin(),vectPoint3D.end());
			for(int i=0; i< numberOfPoints-1;i++){
				if((vectPoint3D[i].phase == -1) or (vectPoint3D[i+1].phase==-1))
					distances[i] = -1;
				else
					distances[i] = vectPoint3D[i+1].phase - vectPoint3D[i].phase;
				if(distances[i]>maxDistances) maxDistances = distances[i];
				//sumDistances += distances[i];
			}
				if(vectPoint3D[numberOfPoints-1].phase == -1 or vectPoint3D[0].phase ==-1)
					distances[numberOfPoints-1] = -1;
				else
					distances[numberOfPoints-1] = 2*PI - (vectPoint3D[numberOfPoints-1].phase - vectPoint3D[0].phase);
				if(distances[numberOfPoints-1] > maxDistances) maxDistances = distances[numberOfPoints-1];
				//sumDistances += distances[numberOfPoints-1];
				if(maxDistances != -1){
					linePrintedBefore = true;
					sumDistances = 2*PI - maxDistances;
					printf("Rx = %f, Ry = %f, D = %f\n",Rx, Ry, sumDistances);
					out << Rx << " " << Ry << " " << sumDistances <<"\n";

					tmpIndices = (int)(sumDistances / (2*M_PI/n));
					if(tmpIndices >=0 and tmpIndices<=n)
						partitionPhase[tmpIndices]++;

				}
				else if(linePrintedBefore){
					linePrintedBefore = false;
					printf("Rx = %f, Ry = %f, D = -1\n",Rx, Ry);
					out <<"\n";

				}
				else printf("Rx = %f, Ry = %f, D = -1\n",Rx, Ry);
				out.close();
		}
		if(linePrintedBefore){
			ofstream out1(perturbationPath, ios::app);
			out1 << "\n";
			out1.close();
		}
	}

	char partitionPhaseFile[255];
	sprintf(partitionPhaseFile, "%s/DSTest_partitionPhase.dat", rootPath);
	ofstream outP(partitionPhaseFile);
	for(int i=0; i < n;i++)
		outP << (i * 2*M_PI/(double)n) << "\t" << partitionPhase[i]/(4*numPoints*numPoints)<<"\n";
	outP.close();

	char sensitivityXYGnuplotPath[255];
	sprintf(sensitivityXYGnuplotPath, "%s/DSTest_PerturbationXY.gp", rootPath);
	ofstream outMe(sensitivityXYGnuplotPath);

	outMe << "set title \"Maximum Phase Shift of the "<< name <<" Attractor\"\n";
	outMe << "set xlabel \"X\"\n";
	outMe << "set xrange["<<-Rxmax<<":"<<Rxmax<<"]\n";
	outMe << "set ylabel \"Y\"\n";
	outMe << "set yrange["<<-Rymax<<":"<<Rymax<<"]\n";
	outMe << "set zrange[0:"<< 2*M_PI << "]\n";
	outMe << "set zlabel \"Maximum Phase Shift\"\n";
	outMe << "unset key\n";
	outMe << "splot \"DSTest_PerturbationXY.dat\" title 'Maximum Phase Shift of the "<< name << " Attractor' with pm3d\n";
	outMe << "pause -1\n";
	outMe.close();
}

void CyclicAttractor::tracerPhaseSensitivityToPerturbationAccordingToPhi(std::vector<point3D>* initialPoints, double phi, double Rmax, bool gnuplotFile){
	char perturbationPath[255];
		sprintf(perturbationPath, "%s/DSTest_Perturbation_Phi.dat", rootPath,phi);
		double phase;
		double phasePI;
		double maxDistances = 0;
		double sumDistances = 0;
		int numberOfPoints = initialPoints->size();
		vector<point3D> vectPoint3D;
		vector<double> distances;
		MTRand myRandom;
		for(int i=0; i< numberOfPoints;i++){
			distances.push_back(0);
			vectPoint3D.push_back(point3D(0,0,0));
		}

		for(double R = 0; R  <= Rmax; R+=0.1){
			maxDistances = 0;
			sumDistances = 0;
			ofstream out(perturbationPath, ios::app);
			for(int i=0; i < numberOfPoints; i++){
				vectPoint3D[i].x = (*initialPoints)[i].x + R*cos(phi);
				vectPoint3D[i].y = (*initialPoints)[i].y + R*sin(phi);
				phase = phaseOf(&vectPoint3D[i],0);
				phasePI = (2*PI*phase/periode);
				while(phasePI >= 2*PI) phasePI -= 2*PI;
				vectPoint3D[i].phase= phasePI;
			}
			sort(vectPoint3D.begin(),vectPoint3D.end());
			for(int i=0; i< numberOfPoints-1;i++){
				distances[i] = vectPoint3D[i+1].phase - vectPoint3D[i].phase;
				if(distances[i]>maxDistances) maxDistances = distances[i];
				sumDistances += distances[i];
			}
				distances[numberOfPoints-1] = 2*PI - (vectPoint3D[numberOfPoints-1].phase - vectPoint3D[0].phase);
				if(distances[numberOfPoints-1] > maxDistances) maxDistances = distances[numberOfPoints-1];
				sumDistances += distances[numberOfPoints-1];
				sumDistances -= maxDistances;
				printf("R = %f, Phi = %f, D = %f\n",R, phi, sumDistances);
				out << R*cos(phi) << " " << R*sin(phi) << " " << sumDistances <<"\n";
				out.close();
		}
		ofstream out1(perturbationPath, ios::app);
		out1 << "\n";
		out1.close();
		if(gnuplotFile){
			char sensitivityGnuplotPath[255];
			sprintf(sensitivityGnuplotPath, "%s/DSTest_Perturbation_Phi.gp", rootPath);
			ofstream outMe(sensitivityGnuplotPath);

			outMe << "set title \"Maximum Phase Shift of the "<< name <<" Attractor\"\n";
			outMe << "set xlabel \"X\"\n";
			outMe << "set xrange["<<-Rmax<<":"<<Rmax<<"]\n";
			outMe << "set ylabel \"Y\"\n";
			outMe << "set yrange["<<-Rmax<<":"<<Rmax<<"]\n";
			outMe << "set zrange[0:"<< 2*M_PI << "]\n";
			outMe << "set zlabel \"Maximum Phase Shift\"\n";
			outMe << "unset key\n";
			outMe << "splot \"DSTest_Perturbation_Phi.dat\" title 'Maximum Phase Shift of the "<< name << " Attractor' with pm3d\n";
			outMe << "pause -1\n";
			outMe.close();
		}
}

void CyclicAttractor::tracerPhaseSensitivityToPerturbationAccordingToPhi(int nbPointsToTest, double phi, double Rmax, bool gnuplotFile){
	vector<point3D> initialPoints;
	extractPoints(nbPointsToTest, &initialPoints);
	tracerPhaseSensitivityToPerturbationAccordingToPhi(&initialPoints,phi,Rmax, gnuplotFile);
}

void CyclicAttractor::tracerPhase(double deltaT, double initial, int numberOfIterations){
	char phasePath[255];
	sprintf(phasePath, "%s/phase.dat", rootPath);
	ofstream out(phasePath);
	this->deltaT = deltaT;
	this->numberOfIterations = numberOfIterations;
	double phase;

	phase = atan2(coordinates[1], coordinates[0]);
	out << 0 << " " << phase << "\n";
	for(int i=0; i < numberOfIterations; i++){
		computeMethod(coordinates,1);
		phase = atan2(coordinates[1], coordinates[0]);
		out << (i+1)*deltaT << " " << phase << "\n";
	}
}

void CyclicAttractor::tracerAttractor(double deltaT, int numberOfIterations){

	char attractorPath[255];
	sprintf(attractorPath, "%s/DSTest_Attractor.dat", rootPath);
	ofstream out(attractorPath);

	this->deltaT = deltaT;
	this->numberOfIterations = numberOfIterations;

	// Boucle d'initialisation afin de rapprocher le point initial de l'attracteur.
	computeMethod(coordinates,0); // Pour reinitializer le calcul, exemple dans le cas du pendule couplé
	computeMethod(coordinates, 9999);

	attractor.push_back(new point3D(coordinates[0],coordinates[1],coordinates[2]));

	// Tracer le portrait de phase du système dynamique
	for(int i=0; i < numberOfIterations; i++){
		computeMethod(coordinates,1);
		attractor.push_back(new point3D(coordinates[0],coordinates[1],coordinates[2]));
	}

	// Trouve les points appartenant à l'attracteur
	bool found = false;
	bool inSegment = false;
	double maxDistance = distancePointSegment(attractor[numberOfIterations-2],attractor[numberOfIterations-3], attractor[numberOfIterations-1], &inSegment);
	double distance;
	int i=numberOfIterations - 4;
	attractor[i+1]->inTheAttractor = true;
	attractor[i+2]->inTheAttractor = true;
	attractor[i+3]->inTheAttractor = true;

	/*
	 * Boucle de détéction de l'attracteur
	 */
	double minDistance = -1;
	double tmpDistance;
	int minIndex;
	int maxVerification = 5;
	int count = 0;
	while((i >= 0 and not found) or (found and count < maxVerification)){
		if(not found) attractor[i]->inTheAttractor = true;
		distance = distancePointSegment(attractor[i],attractor[numberOfIterations-3], attractor[numberOfIterations-1], &inSegment);
		if(distance <= maxDistance){
			tmpDistance = distanceP(attractor[i], attractor[numberOfIterations-2]);
			printf("i = %d and distance %f\n", i, tmpDistance);
			if(tmpDistance < minDistance or minDistance < 0){
				minDistance = tmpDistance;
				minIndex = i;
			}
			found=true;
		}
		if(found) count++;
		i--;
	}
	i = minIndex;
	int j = i;
	while(attractor[j]->inTheAttractor == false){
		attractor[j]->inTheAttractor = true;
		j++;
	}

	/*
	 * En partant du dernier point ajouter chercher le premier point
	 * qui pourrait être confondu avec.
	 *
	 */


	// Imprime que l'attracteur
	printf("i = %d\n", i);

	radius[0] = abs(attractor[i+1]->x);
	radius[1] = abs(attractor[i+1]->y);
	radius[2] = abs(attractor[i+1]->z);

	double angle, norme, phase;
	for(int j=i; j< numberOfIterations-1; j++){//j=i+1
		if(abs(attractor[j]->x) > radius[0]) radius[0] = abs(attractor[j]->x);
		if(abs(attractor[j]->y) > radius[1]) radius[1] = abs(attractor[j]->y);
		if(abs(attractor[j]->z) > radius[2]) radius[2] = abs(attractor[j]->z);
		angle = atan2(attractor[j+1]->y -attractor[j]->y, attractor[j+1]->x-attractor[j]->x) + 2 * PI;
		norme = sqrt(pow(attractor[j+1]->y -attractor[j]->y,2.0) + pow(attractor[j+1]->x-attractor[j]->x,2.0));
		phase = ((j-i) * 2 * PI) / (numberOfIterations-i);
		out << attractor[j]->x << "\t" << attractor[j]->y << "\t" << attractor[j]->z << "\t" << angle << "\t" << norme << "\t" << phase << "\t" << "\t0\n";
	}

	angle = atan2(attractor[numberOfIterations-1]->y -attractor[i]->y, attractor[numberOfIterations-1]->x-attractor[i]->x) + 2 * PI;
	norme = sqrt(pow(attractor[i]->y -attractor[numberOfIterations-1]->y,2.0) + pow(attractor[i]->x-attractor[numberOfIterations-1]->x,2.0));
	phase = ((numberOfIterations-1-i) * 2 * PI) / (numberOfIterations-i);
	out << attractor[numberOfIterations-1]->x << "\t" << attractor[numberOfIterations-1]->y << "\t" << attractor[numberOfIterations-1]->z << "\t" << angle << "\t" << norme << "\t" << phase << "\t" << "0\n";
	out.close();
	// ça pourrait etre utile de délimiter l'attracteur pas une région dans l'espace
	hypothenus = sqrt(radius[0]*radius[0] + radius[1]*radius[1] + radius[2]*radius[2]);
	spatialStep = hypothenus / 100; //100
	printf("hypothenus = %f, spatial Step=%f\n", hypothenus, spatialStep);
	radius[0] *= 1.5;
	radius[1] *= 1.5;
	radius[2] *= 1.5;

	nbPointsPerPeriode = numberOfIterations-i;
	periode  = (nbPointsPerPeriode-1)*deltaT;
	printf("attracteur de période = %f \n deltaT=%f\n", periode, this->deltaT);

	/*
	 * supprimer tout les points qui sont entre [0..i-1]
	 */
	for(int l=0; l<i; l++){
		delete attractor[l];
	}
	if(i > 0) attractor.erase(attractor.begin(), attractor.begin()+i-1);

	/*
	 * Les points contenu dans attractor sont les points de l'attracteur + un points supplémentaire
	 *
	 */
/*
	int nbPointsToRemove = ((int)attractor.size()-1) - nbPointsPerPeriode * (int)((attractor.size()-1) / nbPointsPerPeriode);
	printf("NbAttractor=%d, NBPeriod=%d, Nb point to remove=%d\n", attractor.size(), nbPointsPerPeriode, nbPointsToRemove);
	for(int l=attractor.size()-nbPointsToRemove; l < (int)attractor.size(); l++){
		delete attractor[l];
	}
	for(int l=0; l < nbPointsToRemove; l++){
			attractor.pop_back();
	}*/

	printf("Attractor = %d, pointPerPeriode=%d\n", attractor.size(), nbPointsPerPeriode);
	printf("X firsPoint-lastpoint=%f ", attractor[0]->x -attractor[attractor.size()-1]->x);
	printf("Y firsPoint-lastpoint=%f ", attractor[0]->y -attractor[attractor.size()-1]->y);
	printf("Z firsPoint-lastpoint=%f\n", attractor[0]->z -attractor[attractor.size()-1]->z);

	/*
	 * We compute 2 aditional Points to compare
	 */
	/*double firstPoint[3];
	firstPoint[0] = attractor[attractor.size()-1]->x;
	firstPoint[1] = attractor[attractor.size()-1]->y;
	firstPoint[2] = attractor[attractor.size()-1]->z;
	computeMethod(firstPoint,1);
	attractor.push_back(new point3D(firstPoint[0], firstPoint[1], firstPoint[2]));
	computeMethod(firstPoint,1);
	attractor.push_back(new point3D(firstPoint[0], firstPoint[1], firstPoint[2]));

	printf("/n the first x = %f  y = %f z = %f\n", firstPoint[0], firstPoint[1], firstPoint[2]);
	double secondPoint[3];
	secondPoint[0] = firstPoint[0];
	secondPoint[1] = firstPoint[1];
	secondPoint[2] = firstPoint[2];
	computeMethod(secondPoint, nbPointsPerPeriode - 3);
	printf("/n the first x = %f  y = %f z = %f\n", secondPoint[0], secondPoint[1], secondPoint[2]);
	printf(" differen is %f \n", distanceP(new point3D(firstPoint[0], firstPoint[1], firstPoint[2]), new point3D(secondPoint[0], secondPoint[1], secondPoint[2])));*/
	//this->numberOfIterations = nbPointsPerPeriode * 4;
	char attractorGnuplotPath[255];
	sprintf(attractorGnuplotPath, "%s/DSTest_Attractor.gp", rootPath);
	ofstream outMe(attractorGnuplotPath);

	outMe << "set title \""<< name <<" Attractor\"\n";
	outMe << "set xlabel \"X\"\n";
	outMe << "set xrange["<<-(this->radius[0])<<":"<<this->radius[0]<<"]\n";
	outMe << "set ylabel \"Y\"\n";
	outMe << "set yrange["<<-(this->radius[1])<<":"<<this->radius[1]<<"]\n";
	outMe << "unset key\n";
	if(this->nbDimension > 2){
		outMe << "set zlabel \"Z\"\n";
		outMe << "set zrange["<<-(this->radius[2])<<":"<<this->radius[2]<<"]\n";
		outMe << "splot \"DSTest_Attractor.dat\" title '"<< name << " Attractor' with lines\n";
	}
	else{
		outMe << "plot \"DSTest_Attractor.dat\" title '"<< name << " Attractor' with lines\n";
	}
	outMe << "pause -1\n";
}

/*
 * Return the index of the first point in the attractor that respect
 */
int CyclicAttractor::seekTheIndexOfTheFirstPointIn(double * min, double * max){
		for(int i=0; i < (int)attractor.size();i++){
				if((min[0] <= attractor[i]->x) and (max[0] > attractor[i]->x) and (min[1] <= attractor[i]->y) and (max[1] > attractor[i]->y)) return i;
		}
		return 0;
}

void CyclicAttractor::tracerIsochrons(int nbIsochrons, int nbPointsPerIsochrons, int nbIsoKs, int nbPointsPerIsoKs, bool generateIsoks, double F){
/*
 * Cet algo est inspiré de la méthode de monte-carlo: on fait des tirs dans l'espace
 * des phases et pour chaque tirs (s'il est admissible: contenu dans le bassin d'attraction)
 * on détérmine sa phase (là ou il rejoint l'attracteur) si on a suffisament de points
 * pour une phase donnée on les classes pour tracer la fibre isochrone correspondante.
 * Cet algo est barbare mais est parallélisable BINGO!! ALler il faut l'ecrire assez
 * de description, envisageons alors d'extraire un ensemble de points de l'attracteur
 * et donnons à chacun d'entre eux un score qui correspond au nombre de point trouvés
 * qui ont cette phase (A une formulation prés! :p ).
 *
 */

	double point[255];
	this->nbIsochrons = nbIsochrons;
	this->nbPointsPerIsochrons = nbPointsPerIsochrons;
	this->nbIsoKs = nbIsoKs;
	this->nbPointsPerIsoKs = nbPointsPerIsoKs;
	int indexOfIsochronsInTheAttractor;

	indexOfFirstPhase = seekTheIndexOfTheFirstPointIn(minFirstCoor,maxFirstCoor);
	//double F = 10;//pow(10.0, 10.0);
	isochronsFitness = this->periode / (nbIsochrons*2*F); //0.00000000001;//0.00000000001;//
	printf("isochronsFitness = %f of periode %f\n", isochronsFitness, periode);
	point3D p3d(0,0,0);
	char isofilePath[255];
	sprintf(isofilePath, "%s/Isochrons.dat", isochronsPath);
	ofstream out(isofilePath);
	char fpaths[255];
	char fKpaths[255];
	ofstream iso;
	int maxInMemory = 3;
	this->generateIsoks = generateIsoks;


	// Generate the gnuplot file to charge isochrons
	char isochronsGnuplotPath[255];
	sprintf(isochronsGnuplotPath, "%s/DSTest_Isochrons.gp", rootPath);
	ofstream outMe(isochronsGnuplotPath);

	outMe << "#Isochrons of "<< name << " Attractor\n";
	if(this->nbParam > 0){
		outMe << "#Attractor having the following parameters: \n";
		for(int i=0;i<nbParam;i++){
			outMe << "# \t Param " << i << ": " << param[i] << "\n";
		}
	}
	outMe << "# Attractor period is " << periode << "\n";
	outMe << "# Time step deltaT = " << deltaT << "\n";
	outMe << "# Number of Isochrons: " << nbIsochrons << "\n";
	outMe << "# The precision parameter is " << F << "\n";
	outMe << "# The fitness parameter is " << isochronsFitness << "\n\n";
	outMe << "set title \"Isochrons of the "<< name <<" Attractor\"\n";
	outMe << "set xlabel \"X\"\n";
	outMe << "set ylabel \"Y\"\n";
	//outMe << "set autoscale\n";
	outMe << "unset key\n";
	if(this->nbDimension > 2){
		outMe << "set zlabel \"Z\"\n";
		outMe << "splot \"DSTest_Attractor.dat\" with lines";// title '"<< name << " Attractor' with lines\n";
		for(int i=0; i < nbIsochrons; i++) outMe << ",\"Isochrons/DSTest_Isochron_"<<i<<".dat\" with points ";
		outMe << "title 'Isochrons of the "<< name << " Attractor'\n";
	}
	else{
		outMe << "plot \"DSTest_Attractor.dat\" with lines"; // title '"<< name << " Attractor' with lines\n";
		for(int i=0; i < nbIsochrons; i++) outMe << ",\"Isochrons/DSTest_Isochron_"<<i<<".dat\" with points ";
		outMe << "title 'Isochrons of the "<< name << " Attractor'\n";
	}
	outMe << "pause -1\n";
	outMe.close();

	printf("Thread shitty parts is over\n");
	//Random::Uniform<double>(0.0,radius[0]);
	Random::Randomize(5);
	if(nbIsochrons < (1+(int)(periode/deltaT)))  deltaIsochrons = nbPointsPerPeriode / nbIsochrons;
	else deltaIsochrons = 1;
	 if(nbPointsPerPeriode > 2){
		 for(int i=0; i <= nbIsoKs; i++){
			 sprintf(fKpaths, "%s/DSTest_IsoK_%d.dat",isochronsPath, i+1);
			 isoKs.push_back(new Isochrons(i, fKpaths, maxInMemory, i, indexOfIsochronsInTheAttractor, attractor[indexOfIsochronsInTheAttractor]));
			 isoKsMutexs.push_back(new pthread_mutex_t());
			 pthread_mutex_init(isoKsMutexs[i], NULL);
		 }
		 //pthread_mutex_init(randoM, NULL);
		 for(int i=0; i < nbIsochrons; i++){
				sprintf(fpaths, "%s/DSTest_Isochron_%d.dat",isochronsPath, i);
				saved[i] = false;
				// Considering the first point
				indexOfIsochronsInTheAttractor = (i * nbPointsPerPeriode / nbIsochrons)+indexOfFirstPhase;
				while(indexOfIsochronsInTheAttractor>=(int)attractor.size()) indexOfIsochronsInTheAttractor-=attractor.size();
				// End Considering the first point
//				isochrons.push_back(new Isochrons(indexOfIsochronsInTheAttractor*deltaT, fpaths, maxInMemory, i, indexOfIsochronsInTheAttractor, attractor[indexOfIsochronsInTheAttractor]));
				isochrons.push_back(new Isochrons((i * nbPointsPerPeriode / nbIsochrons)*deltaT, fpaths, maxInMemory, i, indexOfIsochronsInTheAttractor, attractor[indexOfIsochronsInTheAttractor]));
				point[0] = attractor[indexOfIsochronsInTheAttractor]->x;
				point[1] = attractor[indexOfIsochronsInTheAttractor]->y;
				point[2] = attractor[indexOfIsochronsInTheAttractor]->z;
				printf("iso %d index=%d,  x=%f, y=%f, z=%f, phase=%f\n",i, indexOfIsochronsInTheAttractor, point[0], point[1], point[2], indexOfIsochronsInTheAttractor*deltaT);
				isochrons[(int)isochrons.size() - 1]->points.push_back(new point3D(point[0],point[1],point[2]));
				smartPencils.push_back(new SmartPencil(this, isochrons[(int)isochrons.size() - 1]));
				myThreads.push_back(new pthread_t());
				mutexs.push_back(new pthread_mutex_t());
				pthread_mutex_init(mutexs[i], NULL);
		 }
		/* myThreads.push_back(new pthread_t());
		 pthread_create(myThreads[0], NULL, runMe, (void*) smartPencils[0]);*/
		 for(int i=0; i < nbIsochrons; i++){
			 if(pthread_create(myThreads[i], NULL, runMe, (void*) smartPencils[i]) != 0){
				 printf("Thread %d creation failed: %d\n", i, 0);
			 }
		 }
/*
		 while (k<totalPoints){
			 //smartPencils[10]->runMe();

			 if(selectedIsochrons != -1){
				 k=0;
				 for(int l=0; l < isochrons.size(); l++){
					 if(isochrons[l]->nbPoints > nbPointsPerIsochrons)	 k+=nbPointsPerIsochrons;
					 else k+=isochrons[l]->nbPoints;
				 }
				 //printf("k = %d of %d\n", k, totalPoints);
			 }
		 }*/
		 // Solution provisoire normalement attendre aussi le dernier sauf que le dernier n'existe pas du coup on s'en fout pour l'instant 27/01/2010
		 for(int i=0; i < nbIsochrons-1; i++){
			 pthread_join(*myThreads[i], NULL);
			 printf("Pencil %d finished\n", i);
			 delete smartPencils[i];
			 smartPencils[i] = NULL;
		 }
		 // Provisoire
		 pthread_cancel(*myThreads[nbIsochrons-1]);
		 delete smartPencils[nbIsochrons-1];
		 smartPencils[nbIsochrons-1] = NULL;

		 // Provisoire
		 while(mutexs.size()>0){
			 delete mutexs[mutexs.size()-1];
			 mutexs.pop_back();
			 delete isochrons[isochrons.size()-1];
			 isochrons.pop_back();
		 }
		 while(isoKsMutexs.size()>0){
		 	 delete isoKsMutexs[isoKsMutexs.size()-1];
		 	 isoKsMutexs.pop_back();
		 	 delete isoKs[isoKs.size()-1];
		 	 isoKs.pop_back();
		 }
		 printf("done\n");
		 printf("saved\n");
	 }

}


void CyclicAttractor::tracerIsochronsInAWindow(int nbIsochrons, int nbPointsPerIsochrons, int nbIsoKs, int nbPointsPerIsoKs, bool generateIsoks, double F, double window[]){
/*
 * On fait des tirs dans l'espace
 * des phases et pour chaque tirs (s'il est admissible: contenu dans le bassin d'attraction)
 * on détérmine sa phase (là ou il rejoint l'attracteur) si on a suffisament de points
 * pour une phase donnée on les classes pour tracer la fibre isochrone correspondante.
 * Cet algo est barbare mais est parallélisable BINGO!! ALler il faut l'ecrire assez
 * de description, envisageons alors d'extraire un ensemble de points de l'attracteur
 * et donnons à chacun d'entre eux un score qui correspond au nombre de point trouvés
 * qui ont cette phase (A une formulation prés! :p ).
 *
 */

	double point[255];
	this->nbIsochrons = nbIsochrons;
	this->nbPointsPerIsochrons = nbPointsPerIsochrons;
	this->nbIsoKs = nbIsoKs;
	this->nbPointsPerIsoKs = nbPointsPerIsoKs;
	int indexOfIsochronsInTheAttractor;
	indexOfFirstPhase = seekTheIndexOfTheFirstPointIn(minFirstCoor,maxFirstCoor);
	//double F = 10;//pow(10.0, 10.0);
	isochronsFitness = this->periode / (nbIsochrons*2*F); //0.00000000001;//0.00000000001;//
	printf("isochronsFitness = %f of periode %f\n", isochronsFitness, periode);
	point3D p3d(0,0,0);
	char isofilePath[255];
	sprintf(isofilePath, "%s/Isochrons.dat", isochronsPath);
	ofstream out(isofilePath);
	char fpaths[255];
	char fKpaths[255];
	ofstream iso;
	int maxInMemory = 10;
	this->generateIsoks = generateIsoks;


	// Generate the gnuplot file to charge isochrons
	char isochronsGnuplotPath[255];
	sprintf(isochronsGnuplotPath, "%s/DSTest_Isochrons.gp", rootPath);
	ofstream outMe(isochronsGnuplotPath);

	outMe << "#Isochrons of "<< name << " Attractor\n";
	if(this->nbParam > 0){
		outMe << "#Attractor having the following parameters: \n";
		for(int i=0;i<nbParam;i++){
			outMe << "# \t Param " << i << ": " << param[i] << "\n";
		}
	}
	outMe << "# Attractor period is " << periode << "\n";
	outMe << "# Time step deltaT = " << deltaT << "\n";
	outMe << "# Number of Isochrons: " << nbIsochrons << "\n";
	outMe << "# The precision parameter is " << F << "\n";
	outMe << "# The fitness parameter is " << isochronsFitness << "\n\n";
	outMe << "set title \"Isochrons of the "<< name <<" Attractor\"\n";
	outMe << "set xlabel \"X\"\n";
	outMe << "set ylabel \"Y\"\n";
	//outMe << "set autoscale\n";
	outMe << "unset key\n";
	if(this->nbDimension > 2){
		outMe << "set zlabel \"Z\"\n";
		outMe << "splot \"DSTest_Attractor.dat\" with lines";// title '"<< name << " Attractor' with lines\n";
		for(int i=0; i < nbIsochrons; i++) outMe << ",\"Isochrons/DSTest_Isochron_"<<i<<".dat\" with points ";
		outMe << "title 'Isochrons of the "<< name << " Attractor'\n";
	}
	else{
		outMe << "plot \"DSTest_Attractor.dat\" with lines"; // title '"<< name << " Attractor' with lines\n";
		for(int i=0; i < nbIsochrons; i++) outMe << ",\"Isochrons/DSTest_Isochron_"<<i<<".dat\" with points ";
		outMe << "title 'Isochrons of the "<< name << " Attractor'\n";
	}
	outMe << "pause -1\n";
	outMe.close();


	//Random::Uniform<double>(0.0,radius[0]);
	Random::Randomize(5);
	if(nbIsochrons < (1+(int)(periode/deltaT)))  deltaIsochrons = nbPointsPerPeriode / nbIsochrons;
	else deltaIsochrons = 1;
	 if(nbPointsPerPeriode > 2){
		 for(int i=0; i <= nbIsoKs; i++){
			 sprintf(fKpaths, "%s/DSTest_IsoK_%d.dat",isochronsPath, i+1);
			 isoKs.push_back(new Isochrons(i, fKpaths, maxInMemory, i, indexOfIsochronsInTheAttractor, attractor[indexOfIsochronsInTheAttractor]));
			 isoKsMutexs.push_back(new pthread_mutex_t());
			 pthread_mutex_init(isoKsMutexs[i], NULL);
		 }
		 //pthread_mutex_init(randoM, NULL);
		 for(int i=0; i < nbIsochrons; i++){
				sprintf(fpaths, "%s/DSTest_Isochron_%d.dat",isochronsPath, i);
				saved[i] = false;
				// Considering the first point
				indexOfIsochronsInTheAttractor = (i * nbPointsPerPeriode / nbIsochrons)+indexOfFirstPhase;
				while(indexOfIsochronsInTheAttractor>=(int)attractor.size()) indexOfIsochronsInTheAttractor-=attractor.size();
				// End Considering the first point
				//isochrons.push_back(new Isochrons(indexOfIsochronsInTheAttractor*deltaT, fpaths, maxInMemory, i, indexOfIsochronsInTheAttractor, attractor[indexOfIsochronsInTheAttractor]));
				isochrons.push_back(new Isochrons((i * nbPointsPerPeriode / nbIsochrons)*deltaT, fpaths, maxInMemory, i, indexOfIsochronsInTheAttractor, attractor[indexOfIsochronsInTheAttractor]));
				point[0] = attractor[indexOfIsochronsInTheAttractor]->x;
				point[1] = attractor[indexOfIsochronsInTheAttractor]->y;
				point[2] = attractor[indexOfIsochronsInTheAttractor]->z;
				printf("iso %d index=%d,  x=%f, y=%f, z=%f, phase=%f\n",i, indexOfIsochronsInTheAttractor, point[0], point[1], point[2], indexOfIsochronsInTheAttractor*deltaT);
				isochrons[(int)isochrons.size() - 1]->points.push_back(new point3D(point[0],point[1],point[2]));
				smartPencils.push_back(new SmartPencil(this, isochrons[(int)isochrons.size() - 1]));
				smartPencils[(int)smartPencils.size()-1]->setWindow(window);
				myThreads.push_back(new pthread_t());
				mutexs.push_back(new pthread_mutex_t());
				pthread_mutex_init(mutexs[i], NULL);
		 }
		/* myThreads.push_back(new pthread_t());
		 pthread_create(myThreads[0], NULL, runMe, (void*) smartPencils[0]);*/
		 for(int i=0; i < nbIsochrons; i++){
			 if(pthread_create(myThreads[i], NULL, runMe, (void*) smartPencils[i]) != 0){
				 printf("Thread %d creation failed: %d\n", i, 0);
			 }
		 }
/*
		 while (k<totalPoints){
			 //smartPencils[10]->runMe();

			 if(selectedIsochrons != -1){
				 k=0;
				 for(int l=0; l < isochrons.size(); l++){
					 if(isochrons[l]->nbPoints > nbPointsPerIsochrons)	 k+=nbPointsPerIsochrons;
					 else k+=isochrons[l]->nbPoints;
				 }
				 //printf("k = %d of %d\n", k, totalPoints);
			 }
		 }*/
		 for(int i=0; i < nbIsochrons; i++){
			 pthread_join(*myThreads[i], NULL);
			 printf("Pencil %d finished\n", i);
			 delete smartPencils[i];
			 smartPencils[i] = NULL;
		 }
		 while(mutexs.size()>0){
			 delete mutexs[mutexs.size()-1];
			 mutexs.pop_back();
			 delete isochrons[isochrons.size()-1];
			 isochrons.pop_back();
		 }
		 while(isoKsMutexs.size()>0){
		 	 delete isoKsMutexs[isoKsMutexs.size()-1];
		 	 isoKsMutexs.pop_back();
		 	 delete isoKs[isoKs.size()-1];
		 	 isoKs.pop_back();
		 }
		 printf("done\n");
		 printf("saved\n");
	 }

}


/*
 * Fonction du choix de la méthode d'intégration
 */
inline void CyclicAttractor::computeMethod(double * coordinates, int nbCompute){
	//compute(coordinates, nbCompute);
	compute_RK4(coordinates, nbCompute);
}

void * runMe(void * arg){
	double epsilon[3];
	epsilon[0]= 5;
	epsilon[1]= 5;
	epsilon[2]= 7;
	SmartPencil * smartPencil = (SmartPencil*) arg;
	//smartPencil->initializeOldPoint();
	//smartPencil->runMe();
	//smartPencil->runMeOnNeighbourhood(epsilon);
	smartPencil->runMeNoDep();
	return NULL;
}

int CyclicAttractor::seekIsochronsAndAdd3(point3D* p3d, double* myPhase, bool addPoint){
	double tmpPhase = phaseOf(p3d, nbPointsPerIsochrons);
	*myPhase = tmpPhase;
	int inIsochrons = -1;
	bool added = false;
	int indexOfIsochronsInTheAttractor;
	double minDistance, tmpDistance;
	int minIndex;
	double localRadius, tmpLocalRadius;
	 /*
	  * Chercher l'isochrons le plus proche et voir si le seuil de tolérance
	  * est accépté.
	  *
	  */
		if(tmpPhase < 0) return -1;
		int indexPointInTheAttractor = (int)(tmpPhase / deltaT);
		 int indexOfIsochrons = (int)(indexPointInTheAttractor / (double)deltaIsochrons);
		 if(indexOfIsochrons >= (int)isochrons.size()){
			 indexOfIsochrons = 0;
			 //indexPointInTheAttractor = 0;
		 }
		 indexOfIsochronsInTheAttractor = indexOfIsochrons * deltaIsochrons;

		 minDistance = distanceP(attractor[indexPointInTheAttractor], attractor[indexOfIsochronsInTheAttractor]);
		 minIndex = indexOfIsochrons;
		 if(indexOfIsochrons == (int)isochrons.size() - 1){
			 indexOfIsochrons = 0;
			 indexOfIsochronsInTheAttractor = 0;
		 }
		 else{
			 indexOfIsochrons++;
			 indexOfIsochronsInTheAttractor = deltaIsochrons * indexOfIsochrons;
		 }

		 tmpDistance = distanceP(attractor[indexPointInTheAttractor], attractor[indexOfIsochronsInTheAttractor]);
		 if(tmpDistance < minDistance){
			 minDistance = tmpDistance;
			 minIndex = indexOfIsochrons;
		 }
		 indexOfIsochrons = minIndex;
		 indexOfIsochronsInTheAttractor = indexOfIsochrons * deltaIsochrons;

		 /*
		  * Détermine la résolution locale
		  */
		 //printf("%d and %d for a total of %d\n", indexPointInTheAttractor, indexOfIsochronsInTheAttractor, attractor.size());
	/*	 localRadius = distanceP(attractor[indexOfIsochronsInTheAttractor], attractor[indexOfIsochronsInTheAttractor+1]);
		 if(indexOfIsochronsInTheAttractor == 0){
			 tmpLocalRadius = distanceP(attractor[attractor.size()-3], attractor[attractor.size()-1]);
		 }
		 else{
			 tmpLocalRadius = distanceP(attractor[indexOfIsochronsInTheAttractor], attractor[indexOfIsochronsInTheAttractor-1]);
		 }
		 if(localRadius > tmpLocalRadius) localRadius = tmpLocalRadius;*/
/*
		 if(indexOfIsochronsInTheAttractor == (int)attractor.size()-2){
			maxDistance = distancePointSegment(attractor[indexOfIsochronsInTheAttractor+1], attractor[indexOfIsochronsInTheAttractor], attractor[1], &inSegment);
			distance = distancePointSegment(p3d, attractor[indexOfIsochronsInTheAttractor], attractor[1], &inSegment);
		}
		else if(indexOfIsochronsInTheAttractor == (int)attractor.size()-1){
			maxDistance = distancePointSegment(attractor[1], attractor[indexOfIsochronsInTheAttractor], attractor[2], &inSegment);
			distance = distancePointSegment(p3d, attractor[indexOfIsochronsInTheAttractor], attractor[2], &inSegment);
		}
		else {
			maxDistance = distancePointSegment(attractor[indexOfIsochronsInTheAttractor+1], attractor[indexOfIsochronsInTheAttractor], attractor[indexOfIsochronsInTheAttractor+2], &inSegment);
			distance = distancePointSegment(p3d, attractor[indexOfIsochronsInTheAttractor], attractor[indexOfIsochronsInTheAttractor+2], &inSegment);
		}
		if(distance <= maxDistance) {*/
		//if(minDistance <= (localRadius * isochronsFitness)){
		//if(minDistance <= (0.001 * isochronsFitness)){
		 if(indexOfIsochrons == 0){
			 tmpPhase -= periode;
		 }
		 if(fabs(tmpPhase - indexOfIsochrons*((double)deltaT*deltaIsochrons))<= isochronsFitness){
			//printf("my phase is %f and isochron phase is %f\n", tmpPhase, isochrons[indexOfIsochrons]->phase);
			pthread_mutex_lock(mutexs[indexOfIsochrons]);
			//pthread_mutex_lock(mutexs[0]);
			//if(isochrons[indexOfIsochrons]->nbPoints < nbPointsPerIsochrons){
				 if(addPoint){
					 added = isochrons[indexOfIsochrons]->addPoint(p3d);
					 for(int l=0; l < (int)(p3d->isoPoints.size()-1); l++){
						 added = isochrons[indexOfIsochrons]->addPoint(p3d->isoPoints[l]);
					 }
				 }
			 /*}
			 else if(not saved[indexOfIsochrons]){
				 isochrons[indexOfIsochrons]->saveAll();
				 saved[indexOfIsochrons] = true;
			 }*/
			 inIsochrons = indexOfIsochrons;
			 //pthread_mutex_lock(mutexs[0]);
			pthread_mutex_unlock(mutexs[indexOfIsochrons]);
		 }

		//printf("init done %s\n", fpaths);
		if(generateIsoks){
			if(p3d->k <= nbIsoKs){
				pthread_mutex_lock(isoKsMutexs[p3d->k-1]);
				added = isoKs[p3d->k-1]->addPoint(p3d);
				pthread_mutex_unlock(isoKsMutexs[p3d->k-1]);
			}
			else{
				pthread_mutex_lock(isoKsMutexs[nbIsoKs-1]);
				added = isoKs[nbIsoKs-1]->addPoint(p3d);
				pthread_mutex_unlock(isoKsMutexs[nbIsoKs-1]);
			}
			 for(int l=0; l < (int)(p3d->isoPoints.size()-1); l++){
				 if(p3d->isoPoints[l]->k <= nbIsoKs){
					 pthread_mutex_lock(isoKsMutexs[p3d->isoPoints[l]->k-1]);
					 added = isoKs[p3d->isoPoints[l]->k-1]->addPoint(p3d->isoPoints[l]);
					 pthread_mutex_unlock(isoKsMutexs[p3d->isoPoints[l]->k-1]);
				 }
				 else{
					 pthread_mutex_lock(isoKsMutexs[nbIsoKs-1]);
					 added = isoKs[nbIsoKs-1]->addPoint(p3d->isoPoints[l]);
					 pthread_mutex_unlock(isoKsMutexs[nbIsoKs-1]);
				 }
			 }
		}
		return inIsochrons;
}

int CyclicAttractor::seekIsochronsAndAdd(point3D* p3d, double* myPhase, bool addPoint){
	double tmpPhase = phaseOf(p3d, nbPointsPerIsochrons);
	*myPhase = tmpPhase;
	int inIsochrons = -1;
	bool added = false;
	double minDistance, tmpDistance;
	 /*
	  * Chercher l'isochrons le plus proche et voir si le seuil de tolérance
	  * est accépté.
	  *
	  */
		if(tmpPhase < 0) return -1;
		/*if(tmpPhase < isochrons[1]->phase and tmpPhase <= isochronsFitness) printf("OK\n");
		if(tmpPhase > isochrons[29]->phase and (periode-tmpPhase) <= isochronsFitness) printf("OK\n");*/
		int indexOfIsochron = tmpPhase / (periode/nbIsochrons);
		if(indexOfIsochron >= nbIsochrons) indexOfIsochron = nbIsochrons - 1;

		minDistance = fabs(tmpPhase - isochrons[indexOfIsochron]->phase);

		if(indexOfIsochron == (int)nbIsochrons - 1){
			 tmpDistance = fabs(periode - tmpPhase);
			 if(tmpDistance < minDistance) {
				 indexOfIsochron = 0;
				 tmpPhase = tmpDistance;
				 minDistance = tmpDistance;
			 }
		 }
		 else{
			tmpDistance = fabs(tmpPhase - isochrons[indexOfIsochron+1]->phase);
			if(tmpDistance < minDistance) {
				indexOfIsochron++;
				minDistance = tmpDistance;
			}
		 }
		if(minDistance<= isochronsFitness){
			pthread_mutex_lock(mutexs[indexOfIsochron]);
				 if(addPoint){
					 added = isochrons[indexOfIsochron]->addPoint(p3d);
					 for(int l=0; l < (int)(p3d->isoPoints.size()-1); l++){
						 added = isochrons[indexOfIsochron]->addPoint(p3d->isoPoints[l]);
					 }
				 }
			 inIsochrons = indexOfIsochron;
			pthread_mutex_unlock(mutexs[indexOfIsochron]);
		 }

		if(generateIsoks){
			if(p3d->k <= nbIsoKs){
				pthread_mutex_lock(isoKsMutexs[p3d->k-1]);
				added = isoKs[p3d->k-1]->addPoint(p3d);
				pthread_mutex_unlock(isoKsMutexs[p3d->k-1]);
			}
			else{
				pthread_mutex_lock(isoKsMutexs[nbIsoKs-1]);
				added = isoKs[nbIsoKs-1]->addPoint(p3d);
				pthread_mutex_unlock(isoKsMutexs[nbIsoKs-1]);
			}
			 for(int l=0; l < (int)(p3d->isoPoints.size()-1); l++){
				 if(p3d->isoPoints[l]->k <= nbIsoKs){
					 pthread_mutex_lock(isoKsMutexs[p3d->isoPoints[l]->k-1]);
					 added = isoKs[p3d->isoPoints[l]->k-1]->addPoint(p3d->isoPoints[l]);
					 pthread_mutex_unlock(isoKsMutexs[p3d->isoPoints[l]->k-1]);
				 }
				 else{
					 pthread_mutex_lock(isoKsMutexs[nbIsoKs-1]);
					 added = isoKs[nbIsoKs-1]->addPoint(p3d->isoPoints[l]);
					 pthread_mutex_unlock(isoKsMutexs[nbIsoKs-1]);
				 }
			 }
		}
		return inIsochrons;
}

int CyclicAttractor::seekIsochronsAndAdd2(point3D* p3d, bool addPoint){
	int temps = -1;
	int indexPointInTheAttractor = inPhaseWith(p3d, nbPointsPerIsochrons, &temps);
	int inIsochrons = -1;
	bool added = false;
	double minDistance, tmpDistance;
	double localRadius, tmpLocalRadius;
	int indexOfIsochrons, indexOfIsochronsInTheAttractor;
	 /*
	  * Chercher l'isochrons le plus proche et voir si le seuil de tolérance
	  * est accépté.
	  *
	  */
		if(indexPointInTheAttractor < 0) return -1;

		// Cherche l'isochron le plus proche
		minDistance = -1;
		indexOfIsochrons = -1;
		for(int i=0; i < (int)isochrons.size(); i++){
			tmpDistance = distanceP(attractor[indexPointInTheAttractor], isochrons[i]->phasePoint);
			if(tmpDistance < minDistance or minDistance < 0){
				minDistance = tmpDistance;
				indexOfIsochrons = i;
			}
		}
		indexOfIsochronsInTheAttractor = isochrons[indexOfIsochrons]->myIndexInTheAttractor;

		 minDistance = distanceP(attractor[indexPointInTheAttractor], isochrons[indexOfIsochrons]->phasePoint);

		 /*
		  * Détermine la résolution locale
		  */
		 localRadius = distanceP(attractor[indexOfIsochronsInTheAttractor], attractor[indexOfIsochronsInTheAttractor+1]);
		 if(indexOfIsochronsInTheAttractor == 0){
			 tmpLocalRadius = distanceP(attractor[attractor.size()-3], attractor[attractor.size()-4]);
		 }
		 else{
			 tmpLocalRadius = distanceP(attractor[indexOfIsochronsInTheAttractor], attractor[indexOfIsochronsInTheAttractor-1]);
		 }
		 if(localRadius > tmpLocalRadius) localRadius = tmpLocalRadius;

		if(minDistance <= (localRadius * isochronsFitness)){
			//printf("my phase is %f and isochron phase is %f\n", tmpPhase, isochrons[indexOfIsochrons]->phase);
			pthread_mutex_lock(mutexs[indexOfIsochrons]);
			//if(isochrons[indexOfIsochrons]->nbPoints < nbPointsPerIsochrons){
				 if(addPoint){
					 added = isochrons[indexOfIsochrons]->addPoint(p3d);
					 for(int l=0; l < p3d->numberOfIsoPoints; l++){
						 added = isochrons[indexOfIsochrons]->addPoint(p3d->isoPoints[l]);
					 }
				 }
			 /*}
			 else if(not saved[indexOfIsochrons]){
				 isochrons[indexOfIsochrons]->saveAll();
				 saved[indexOfIsochrons] = true;
			 }*/
			 inIsochrons = indexOfIsochrons;
		}
			 pthread_mutex_unlock(mutexs[indexOfIsochrons]);

			 if(generateIsoks){
				if(p3d->k <= nbIsoKs){
					pthread_mutex_lock(isoKsMutexs[p3d->k-1]);
					added = isoKs[p3d->k-1]->addPoint(p3d);
					pthread_mutex_unlock(isoKsMutexs[p3d->k-1]);
				}
				else{
					pthread_mutex_lock(isoKsMutexs[nbIsoKs-1]);
					added = isoKs[nbIsoKs-1]->addPoint(p3d);
					pthread_mutex_unlock(isoKsMutexs[nbIsoKs-1]);
				}
				 for(int l=0; l < (int)(p3d->isoPoints.size()-1); l++){
					 if(p3d->isoPoints[l]->k <= nbIsoKs){
						 pthread_mutex_lock(isoKsMutexs[p3d->isoPoints[l]->k-1]);
						 added = isoKs[p3d->isoPoints[l]->k-1]->addPoint(p3d->isoPoints[l]);
						 pthread_mutex_unlock(isoKsMutexs[p3d->isoPoints[l]->k-1]);
					 }
					 else{
						 pthread_mutex_lock(isoKsMutexs[nbIsoKs-1]);
						 added = isoKs[nbIsoKs-1]->addPoint(p3d->isoPoints[l]);
						 pthread_mutex_unlock(isoKsMutexs[nbIsoKs-1]);
					 }
				 }
			}
		 return inIsochrons;
}

/*
 * Return the phase of a state defined by an array
 */
double CyclicAttractor::phaseOf(double* state){
	point3D p3d(state[0],state[1],state[2]);
	return phaseOf(&p3d,0);
}

// La phase est calculée en prenant comme point de départ le premier point dans
// Un phase st égale à 1 correspond à une phase non défini

double CyclicAttractor::phaseOf(point3D * p3d, int numberMaxOfIsoPoints){
	double point[255];
	point[0] = p3d->x;
	point[1] = p3d->y;
	point[2] = p3d->z;
	double result;
	int timeToAttractor = 0;
	bool found = false;
	int i = 0;
	int j = 0;
	double distance;
	double minDistance;
	double localRadius;
	int n=1;
	int computingSteps = (nbPointsPerPeriode-1)/n;
	computeMethod(point,0);
	//computeMethod(point,5*computingSteps);
	bool inSegment;

	point3D tmpPoint3D(point[0],point[1],point[2]);
	p3d->k = 0;
	while(not(found) and (i < numberOfIterations)){
		j=0;
		if(!finite(point[0]) or !finite(point[1]) or !finite(point[2])) return -1;
		tmpPoint3D.x = point[0];
		tmpPoint3D.y = point[1];
		tmpPoint3D.z = point[2];
		minDistance = -1;
		while(j < (int)attractor.size()-2 and not(found)){
			distance = distancePointSegment(&tmpPoint3D, attractor[j], attractor[j+2], &inSegment);
			localRadius = distancePointSegment(attractor[j+1], attractor[j], attractor[j+2], &inSegment);
			//if(distance <= localRadius) {
			if((minDistance == -1) or (distance < minDistance)) minDistance = distance;
			if(distance <= spatialStep or distance <= localRadius) {
				j++;
				found = true;
				break;
			}
			j++;
		}
		computeMethod(point, computingSteps);
		timeToAttractor++;
		i += computingSteps;
		if(p3d->numberOfIsoPoints < numberMaxOfIsoPoints){
				if(((point[0] <= radius[0]*5) or (point[1] <= radius[1]*5) or (point[2] <= radius[2]*5)) and (minDistance >= spatialStep *10)){
				// Modification pour boxon
					//if((fabs(point[0]) <= 600) and (fabs(point[1]) <= 600) and (point[2] <= 800) and (point[2] >= -400)){
						p3d->isoPoints.push_back(new point3D(point[0], point[1], point[2],0,0,0,timeToAttractor));
						p3d->numberOfIsoPoints++;
					//}
				}
		}
	}
	if(found){
		result = (j - indexOfFirstPhase) *deltaT;
		while (result >= periode) result -= periode;
		while (result < 0) result += periode;
		p3d->phase = result;
		p3d->k = timeToAttractor;
		if(p3d->numberOfIsoPoints > 0)
			p3d->isoPoints[p3d->numberOfIsoPoints-1]->k = 0;
		getAngleAndNorme(p3d);
		for(int l=0; l<p3d->numberOfIsoPoints-1;l++){
			p3d->isoPoints[l]->k = timeToAttractor-p3d->isoPoints[l]->k;
		}
		//printf("found\n");
		for(int l=0; l < p3d->numberOfIsoPoints; l++){
			getAngleAndNorme(p3d->isoPoints[l]);
			p3d->isoPoints[l]->phase = result;
		}
		//printf("END found\n");
		//printf("tta=%d,nip=%d\n",timeToAttractor,p3d->numberOfIsoPoints);
		return result;
	}
	else return -1;
}
/*
 * Write the angle and the norme of the flow field
 */
void CyclicAttractor::getAngleAndNorme(point3D * p3d){
	double angle, norme;
	double tmpPoint[3];
	tmpPoint[0] = p3d->x;
	tmpPoint[1] = p3d->y;
	tmpPoint[2] = p3d->z;
	computeMethod(tmpPoint,1);
	angle = atan2(tmpPoint[1] - p3d->y, tmpPoint[0] - p3d->x);
	norme = sqrt(pow(tmpPoint[1] - p3d->y,2.0)+pow(tmpPoint[0] - p3d->x,2.0));
	p3d->angle = angle;
	p3d->norme = norme;
}

/*
 * Routine principale pour détécter la phase d'un point
 * Sachant qu'on a une approximation de la période il faudrait
 * alors faire aux alentours d'une periode/deltaT computations
 * On considere une certaine marge d'erreur qui nous permettera
 * d'évaluer si on est ou pas sur l'attracteur. Le point le plus
 * optimale est retenu.
 */
int CyclicAttractor::inPhaseWith(point3D * p3d, int numberMaxOfIsoPoints, int* temps){
	bool found = false;
	int i = 0, j = 0, n = 1, deltaErrorMax = 12, errorCount = 0, maxVerification = 8,
	count = 0, minIndex;
	double point[255], distance, minDistance, localRadius, tmpLocalRadius;
	int computingSteps = (nbPointsPerPeriode-1)/n;

	point[0] = p3d->x;
	point[1] = p3d->y;
	point[2] = p3d->z;
	p3d->k = 0;
	point3D tmpPoint3D(point[0],point[1],point[2]);

	point3D pointToHold(0.0,0.0,0.0);
	point3D tmpPointToHold(0.0,0.0,0.0);
	double tmpDistance2, minDistance2;
	*temps = 0;
	int timeToAttractor=0;
	while(not(found) and (i < numberOfIterations)){
		j=0;
		if(!finite(point[0]) or !finite(point[1]) or !finite(point[2])) return -1;
		if(point[0]> 5 * radius[0] or point[1]> 5 * radius[1] or point[2] > 5 * radius[2]) return -1;
		tmpPoint3D.x = point[0];
		tmpPoint3D.y = point[1];
		tmpPoint3D.z = point[2];

		count = 0;
		minDistance = -1;
		while((j < (int)attractor.size()-2 and not(found)) or (found and (count < maxVerification))){
			if(found) count++;
			// Calcul la résolution locale (le min de distance avec les 2 points voisins)
			localRadius = distanceP(attractor[j+1], attractor[j]);
			tmpLocalRadius = distanceP(attractor[j+1], attractor[j+2]);
			if(localRadius > tmpLocalRadius) localRadius = tmpLocalRadius;

			distance = distanceP(&tmpPoint3D, attractor[j+1]);
			// Compare la résolution local avec la distance du point
			if(distance <= localRadius) {
				found = true;		// Une solution au moins existe
				if(distance < minDistance or minDistance < 0){	// Retenir la solution optimale
					minDistance = distance;
					minIndex = j+1;
				}
				//break;
			}
			j++;
		}
		if(errorCount == deltaErrorMax and not(found)){
			if(p3d->numberOfIsoPoints < numberMaxOfIsoPoints){
				if((point[0] <= radius[0]*2) or (point[1] <= radius[1]*2) or (point[2] <= radius[2]*2)){
					//p3d->isoPoints.push_back(new point3D(point[0], point[1], point[2]));
					p3d->isoPoints.push_back(new point3D(pointToHold.x, pointToHold.y, pointToHold.z,0,0,0,timeToAttractor));
					p3d->numberOfIsoPoints++;
				}
			}
		}
		if(found) {
			j = minIndex;
			break;
		}
		//computeMethod(point, computingSteps);
		if(errorCount == 0){
			computeMethod(point, computingSteps-deltaErrorMax/2);
			timeToAttractor++;
			*temps++;
			i += computingSteps-deltaErrorMax/2;
			pointToHold.x = point[0];
			pointToHold.y = point[1];
			pointToHold.z = point[2];
			if(p3d->isoPoints.size() == 0){
				minDistance2 = distanceP(&pointToHold, p3d);
			}
			else minDistance2 = distanceP(&pointToHold, p3d->isoPoints[p3d->isoPoints.size()-1]);
		}
		else{
			tmpPointToHold.x = point[0];
			tmpPointToHold.y = point[1];
			tmpPointToHold.z = point[2];
			if(p3d->isoPoints.size() == 0){
				tmpDistance2 = distanceP(&tmpPointToHold, p3d);
			}
			else tmpDistance2 = distanceP(&tmpPointToHold, p3d->isoPoints[p3d->isoPoints.size()-1]);
			if(tmpDistance2 < minDistance2){
				pointToHold.x = tmpPointToHold.x;
				pointToHold.y = tmpPointToHold.y;
				pointToHold.z = tmpPointToHold.z;
				minDistance2 = tmpDistance2;
			}
			computeMethod(point, 1);
			i ++;
		}
		errorCount++;
		if(errorCount > deltaErrorMax) {
			errorCount = 0;
		}
	}
	if(found){
		if(j == (int)attractor.size()-3) j = 0;
		p3d->phase = j*deltaT;
		p3d->k = timeToAttractor;
		getAngleAndNorme(p3d);
		if(p3d->numberOfIsoPoints != 0){
			p3d->isoPoints[p3d->numberOfIsoPoints-1]->k = 0;
			for(int l=0; l<p3d->numberOfIsoPoints-1;l++){
				p3d->isoPoints[l]->k = timeToAttractor-p3d->isoPoints[l]->k;
			}
			//printf("found\n");
			for(int l=0; l < p3d->numberOfIsoPoints; l++){
				getAngleAndNorme(p3d->isoPoints[l]);
				p3d->isoPoints[l]->phase = j*deltaT;
			}
		}
		return j;
	}
	else return -1;
}


/*
 * Renvoie la distance entre un point et un segment si le point se projete en dehors du segment
 * inSegment prend la valeur false, si le point est aligné avec le segment la distance retournée
 * est celle qui le relie avec le plus proche des deux points
 */

double CyclicAttractor::distancePointSegment(point3D * p, point3D * p1, point3D * p2, bool * inSegment){

	double result, a, d, ratio;
	double dx, dy, dz;


	if((p1->x == p2->x) and (p1->y == p2->y)){
		result = sqrt(pow(p->x - p1->x, 2.0) + pow(p->y - p1->y, 2.0) + pow(p->z - p1->z, 2.0));
		*inSegment = true;
	}
	else{
		dx = p2->x - p1->x;
		dy = p2->y - p1->y;
		dz = p2->z - p1->z;
		a = dx * (p->x - p1->x) + dy * (p->y - p1->y) + dz * (p->z - p1->z);
		d = sqrt(dx*dx + dy*dy + dz*dz);
		ratio = a / d;
		if(ratio < 0){
			result = sqrt(pow(p->x - p1->x, 2.0) + pow(p->y - p1->y, 2.0) + pow(p->z - p1->z, 2.0));
			*inSegment = false;
		}
		else if(ratio > 1){
			result = sqrt(pow(p->x - p2->x, 2.0) + pow(p->y - p2->y, 2.0) + pow(p->z - p2->z, 2.0));
			*inSegment = false;
		}
		else{
			result = sqrt(pow(p->x - p1->x - ratio * dx, 2.0) + pow(p->y - p1->y - ratio * dy, 2.0) + pow(p->z - p1->z - ratio * dz, 2.0));
			*inSegment = true;
		}
	}
	return result;
}

/*
 * Calcule la distance entre 2 points
 */
double CyclicAttractor::distanceP(point3D * p, point3D * q){
	double distance = 0.0;
	distance += pow((p->x -q->x), 2.0);
	distance += pow((p->y -q->y), 2.0);
	distance += pow((p->z -q->z), 2.0);
	return sqrt(distance);
}

/*
 * Génnère des points selon un type de fichier bien particulier destiné
 * à être affiché sur un terminal openGL élaboré par Nicolas Glade
 * Le format de fichier est le suivant:
 * X (double) \t   Y (double) \t   Angle (double) \t   Norme (double) \t   Phase (double) \t   k (int)
 * ou X, Y: point initial; Angle (entre point initial et le point suivant); norme;
 * phase (sur l'attracteur); k (le temps qu'il met pour rejoindre l'attracteur)
 */
void CyclicAttractor::gridExploration(double deltaT, int nLines, int nColumns){
	point3D tmpPoint(0,0,0);
	int phasePoint = 0;
	int tempsToAttractor = 0;
	char gridPath[255];
	sprintf(gridPath, "%s/DSTest_Map.dat", rootPath);
	int m=0, n=0;
	ofstream out(gridPath);
	out << "#" << nColumns << "," << nLines << "#\n";
	for(double i=-radius[0]/2.0; i < radius[0]/2.0; i=i+radius[0]/nLines){
		for(double j=-radius[1]/2.0; j < radius[1]/2.0; j = j + radius[1]/nColumns){
			tmpPoint.x = i;
			tmpPoint.y = j;
			tmpPoint.z = 0;
			getAngleAndNorme(&tmpPoint);
			phasePoint = inPhaseWith(&tmpPoint,0,&tempsToAttractor);
			tmpPoint.phase = phasePoint * deltaT/periode;
			tmpPoint.k = tempsToAttractor;
			if(n < (nColumns-1))
				out << tmpPoint.x << "\t" << tmpPoint.y << "\t" << tmpPoint.z << "\t" << tmpPoint.angle << "\t" << tmpPoint.norme << "\t" << tmpPoint.phase << "\t" << tmpPoint.k << "\t";
			else
				out << tmpPoint.x << "\t" << tmpPoint.y << "\t" << tmpPoint.z << "\t" << tmpPoint.angle << "\t" << tmpPoint.norme << "\t" << tmpPoint.phase << "\t" << tmpPoint.k << "\n";
			n++;
		}
		m++;
		printf("%d%% map exploration done\n", m);
	}
	out.close();
}

/*
 * Renvoie un angle compris entre 0..2PI qui relie deux points
 */
double CyclicAttractor::getAngle(double firstPoint[], double nextPoint[]){
	return atan2(nextPoint[1]-firstPoint[1], nextPoint[0]-firstPoint[0]) + 2 * PI;
}

/*
 * L'encapsulation de la méthode Uniform de la classe singleton de Random
 * Assure un accés exclusif à cette méthode. Dans le cas ou plusieurs thread peuvent
 * y accéder en même temps il y a une possiblité d'erreur de type segmentation faults
 * L'utilisation d'un mutex_t permet de remédier à cela
 * temps de trackage de l'erreur (a pris des mois :( ) mais ça a permis de corriger
 * d'autres!!!
 */

double CyclicAttractor::getRandomNum(double minB, double maxB){
	double ret=0;
	double range=(maxB-minB);
	pthread_mutex_lock(&randoM);
	//ret = Random::Uniform(minB, maxB);
	  //ret = fmod(rand(),range)+minB;
	  ret = rand() * ((double)range) / RAND_MAX + minB ;
	pthread_mutex_unlock(&randoM);
	return ret;
}

/*
 * Runge Kutta Integration: order 4
 *
 */
void CyclicAttractor::compute_RK4(double * coordinates, int nbCompute){
	double k1[255],k2[255],k3[255],k4[255];
	double tmpCoor[255];
	int i = 0;

	for(int k=0;k<nbCompute;k++){
		for(i=0; i < nbDimension; i++) tmpCoor[i] = coordinates[i];
		compute(tmpCoor,1);

		for(i=0; i < nbDimension; i++) {
			k1[i] = tmpCoor[i] - coordinates[i];
			tmpCoor[i]=coordinates[i] + (k1[i]/2);
		}
		compute(tmpCoor,1);

		for(i=0; i < nbDimension; i++) {
			k2[i] = (tmpCoor[i] - (coordinates[i] + (k1[i]/2)));
			tmpCoor[i]=coordinates[i] + (k2[i]/2);
		}
		compute(tmpCoor,1);

		for(i=0; i < nbDimension; i++) {
			k3[i] = (tmpCoor[i] - (coordinates[i] + (k2[i]/2)));
			tmpCoor[i]=coordinates[i] + k3[i];
		}
		compute(tmpCoor,1);
		for(i=0; i < nbDimension; i++) {
			k4[i] = tmpCoor[i] - (coordinates[i] + k3[i]);
			coordinates[i] += ((k1[i]/6) + (k2[i]/3) + (k3[i]/3) + (k4[i]/6));
		}
	}
}
