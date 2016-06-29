/*
 * DSGrid.h
 *
 *  Created on: Jan 12, 2010
 *      Author: hbenamor
 *      Adress: TIMC-IMAG Laboratory
 *      The grid of coupled oscillators
 */

#ifndef DSGRID_H_
#define DSGRID_H_
#define MAX_INT 50
#include "CyclicAttractorSmartPencil.h"
#include "random-singleton.h"

class DSGrid {
private:
	int height, width;
	char path[255];
	double Xmax;
	double PXLearningZone;
	double alpha;
	double initialCoupling;
	double filterCoeff;
	//std::vector<CyclicAttractor> oscillators;
	CyclicAttractor * cyclicAttractor;
	/*
	double *** couplages;//[MAX_INT*MAX_INT][MAX_INT*MAX_INT][4]; // (xi xj) (xi yj) (yi xj) (yi yj) Couplage de j sur i
	double ***states;//[MAX_INT][MAX_INT][3];
	double globalSignal[3];			// The global state is an addition of states of the oscillators forming the grid
	double ***tmpStates;//[MAX_INT][MAX_INT][3];
	double **seuils;//[MAX_INT][MAX_INT];
	*/
	double couplages[MAX_INT*MAX_INT][MAX_INT*MAX_INT][4]; // (xi xj) (xi yj) (yi xj) (yi yj) Couplage de j sur i
	double states[MAX_INT][MAX_INT][3];
	double globalSignal[3];			// The global state is an addition of states of the oscillators forming the grid
	double tmpStates[MAX_INT][MAX_INT][3];
	double seuils[MAX_INT][MAX_INT];
	MTRand * mtRand;
public:
	DSGrid(int height, int width, CyclicAttractor * cyclicAttractor);
	DSGrid();

	void compute(int nbCompute);
	void compute_RK4(int nbCompute);
	inline void computeMethod(int nbCompute);

	int getID(int i, int j);
	void Create_Save_BMP();

	void buildCoupling_all(double epsilon);
	void buildSeuils(double epsilon);
	void buildCoupling_VonNeumann(double epsilon);
	void buildCoupling_Moore(double epsilon);
	void buildCoupling_Chain(double epsilon);
	void buildCoupling_PerceptronOneLayerToOneOutput(double epsilon);
	inline double couplingFunction(double input);

	void computeAnimation(int gridNbCompute, int steps);
	void computeGrid(int gridNbCompute, int steps);
	void saveStates(char outPutFile[255]);
	void saveFilterState(char outPutFile[255]);
	void saveLocalSynchrony(char outPutFile[255]);
	void savePhases(char outPutFile[255]);
	inline void saveMethod(char outPutFile[255]);
	void savePhaseLag(int iteration, int simNumber);
	void savePhaseLagResponse();

	void refreshGlobalSignal();
	void saveGlobalSignal(int iteration);
	void saveSignals(int iteration);

	inline void applyNoise(int iteration);
	void perturbate_burnH(double stimulationStrength, double angle);
	void perturbate_burnC(double stimulationStrength, double angle);
	void perturbate_burn8(double stimulationStrength, double angle);
	void perturbate_burn1(double stimulationStrength, double angle);
	void perturbate_burn2(double stimulationStrength, double angle);
	void perturbate_burnD(double stimulationStrength, double angle);
	void perturbate_all(double stimulationStrength, double angle);
	void perturbate_randomly();
	void perturbate_toEvoke1(double stimulationStrength,double angle);

	void perturbate_burnImage(char * imagePath, double stimulationStrength);
	void mapPerturbation(double * state, int grayValue, double stimulationStrength);
	void perturbate_evokeImage(char * imagePath, double stimulationStrength);

	double getPhaseLag_Radian(double * stateMaster, double * stateSlave);

	void setTheSameStateForAll(double x, double y, double z);

	void learningFromNeighbourhood(int i, int j);
	void performLearning();


	void traceMaximumPhaseShift(double stimulationStrength, int nbIterations);
	virtual ~DSGrid();
};

#endif /* DSGRID_H_ */
