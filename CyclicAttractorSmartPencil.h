/*
 * CyclicAttractorSmartPencil.h
 *
 *  Created on: 20-Jan-2009
 *      Author: hbenamor
 */

#ifndef CYCLICATTRACTORSMARTPENCIL_H_
#define CYCLICATTRACTORSMARTPENCIL_H_
#include <pthread.h>
#include <semaphore.h>
#include "Isochrons.h"
#include <fstream>
#include "point3D.h"
#include <vector>
#include "random-singleton.h"
#include "MersenneTwister.h"


class CyclicAttractor;
void * runMe(void* arg);
class SmartPencil {
public:
	SmartPencil(CyclicAttractor* cyclicAttractor, Isochrons * isochrons);
	SmartPencil(CyclicAttractor* cyclicAttractor, Isochrons * isochrons, double * oldOne);
	virtual ~SmartPencil();
	void * runMe();
	friend void * runMe(void* arg);
	void setOldOne(double * oldOne);
	friend void * runMeToo(void * arg);
	void * runMeNoDep();
	void * runMeOnNeighbourhood(double * epsilon);
	SmartPencil * secondPencil;
	void setWindow(double window[]);

private:
	SmartPencil();
	void initialize(CyclicAttractor* cyclicAttractor, Isochrons * isochrons);
	double center[255];
	double average_vector[255];
	int nb_points_found, nbMaxPointsPerStep, nb_max_dep, nb_dep;
	double cote[255];
	double deltaDep[255];
	int nbDimension;
	Isochrons * isochrons;
	CyclicAttractor * cyclicAttractor;
	int direction;
	double hypothenus;
	int normalizeIt(double * myVector);
	void initializeOldPoint();
	double oldOne[255];
	bool headPencil;
	MTRand * mtRand;
};

class CyclicAttractor {
public:
	std::vector<point3D*> attractor;
	std::vector<Isochrons*> isochrons;
	std::vector<Isochrons*> isoKs;
	CyclicAttractor(char * name);
	CyclicAttractor(int nbDimensions, int nbParam, double * myParam, char * name);
	CyclicAttractor(const CyclicAttractor & cyclicAttractor);
	virtual ~CyclicAttractor();

	virtual void compute(double * coordinates, int nbCompute) = 0; // Euler Computation
	void compute_RK4(double * coordinates, int nbCompute);	// nbCompute steps integration using Runge Kutta Methods
	// Fonction du choix de la méthode d'intégration numérique
	inline void computeMethod(double * coordinates, int nbCompute);
	//////////////// ADDED FOR PENDULE
	virtual void resetFather(double coordinates[]) = 0;
	///////////////////
	void tracerInstantPerturbationResponse(double deltaT, double R, double NbPoints);
	void traceOneTrajectory(double * coor, double deltaT, int numberOfIterations);
	void tracerPhase(double deltaT, double initial, int numberOfIterations);
	void tracerPhaseLag(double deltaT, double initial, int numberOfIterations);
	void tracerAttractor(double deltaT, int numberOfIterations);
	void tracerIsochrons(int nbIsochrons, int nbPointsPerIsochrons, int nbIsoKs, int nbPointsPerIsoKs, bool generateIsoks, double F);
	void tracerIsochronsInAWindow(int nbIsochrons, int nbPointsPerIsochrons, int nbIsoKs, int nbPointsPerIsoKs, bool generateIsoks, double F, double window[]);
	void tracerPortraitPhase(double deltaT, int numberOfIteration, int nbShots);
	void tracerPortrait(double deltaT, int numberOfIteration, int nbShots);
	void tracerPhaseSensitivityToPerturbationAccordingToPhi(int nbPointsToTest, double phi, double Rmax, bool gnuplotFile);
	void tracerPhaseSensitivityToPerturbationAccordingToPhi(std::vector<point3D>* initialPoints, double phi, double Rmax, bool gnuplotFile) ;
	void tracerPhaseSensitivityToPerturbation(int nbPointsToTest, double Rmax);
	void tracerPhaseSensitivityToPerturbationXY(int nbPointsToTest, double Rxmax, double Rymax, int nbPointsPerDimension);
	void tracerPortraitWithoutCalculusOfAttractor(double deltaT, int numberOfIteration, int nbShots);

	void generateRandomPointsOnTheAttractor(int numberOfPoints, std::vector<point3D>* initialPoints);
	void extractPointsFromTheAttractor(int numberOfPoints, std::vector<point3D>* initialPoints);
	void extractPoints(int numberOfPoints, std::vector<point3D>* initialPoints);
	int seekTheIndexOfTheFirstPointIn(double * min, double * max);

	void gridExploration(double deltaT, int nLines, int nColumns);

	double phaseOf(double* state);
	double phaseOf(point3D * p3d, int numberMaxOfIsoPoints);
	int inPhaseWith(point3D * p3d, int numberMaxOfIsoPoints, int * temps);
	double distancePointSegment(point3D * p, point3D * p1, point3D * p2, bool * inSegment);
	double distanceP(point3D * p, point3D * q);
	int seekIsochronsAndAdd(point3D* p3d, double * myPhase, bool addPoint);
	int seekIsochronsAndAdd3(point3D* p3d, double * myPhase, bool addPoint);
	int seekIsochronsAndAdd2(point3D* p3d, bool addPoint);
	double getAngle(double firstPoint[], double nextPoint[]);
	void getAngleAndNorme(point3D * p3d);

	double getRandomNum(double minB, double maxB);

	//bool compare(point3D* p1,point3D* p2);

	double radius[255];
	int numberOfIterations;
	double deltaT;
	double periode;
	double resolution;
	int nbPointsPerPeriode;
	char name[255];
	char rootPath[255];
	char isochronsPath[255];
	char trajectoriesPath[255];
	double coordinates[255];
	double param[255];
	int nbDimension, nbParam;
	double minFirstCoor[255],maxFirstCoor[255];
	int indexOfFirstPhase;

	int deltaIsochrons, nbIsochrons, nbPointsPerIsochrons, nbIsoKs, nbPointsPerIsoKs;
	double isochronsFitness;
	bool saved[255];
	std::vector<SmartPencil*> smartPencils;
	std::vector<pthread_t*> myThreads;
	std::vector<pthread_mutex_t *> mutexs;
	std::vector<pthread_mutex_t *> isoKsMutexs;
	pthread_mutex_t randoM;
	bool generateIsoks;
	double spatialStep, hypothenus;
private:
	CyclicAttractor();
	void initializePath(char * name);
};

#endif /* CYCLICATTRACTORSMARTPENCIL_H_ */
