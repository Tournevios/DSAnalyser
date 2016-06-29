/*
 * VariationPendule.cpp
 *
 *  Created on: 24-Feb-2009
 *      Author: hbenamor
 */

#include "VariationPendule.h"

VariationPendule::VariationPendule() : CyclicAttractor("VariationPendule") {
	// TODO Auto-generated constructor stub
	initDefault();
}

VariationPendule::VariationPendule(char * name) : CyclicAttractor(name) {
	// TODO Auto-generated constructor stub
	initDefault();
}

VariationPendule::VariationPendule(double * myParam, char * name) : CyclicAttractor(2, 1, myParam, name) {
	// TODO Auto-generated constructor stub
}

VariationPendule::~VariationPendule() {
	// TODO Auto-generated destructor stub
}

void VariationPendule::initDefault(){
	nbParam = 1;
	//param[0] = 2 * 3.14;
	param[0] = 3.14/2.0;
	nbDimension = 2;
	coordinates[0] = 2.03;
	coordinates[1] = -0.02;
	coordinates[2] = 0;
}
void VariationPendule::resetFather(double coordinates[]){

}


void VariationPendule::compute(double coordinates[], int nbCompute){
	double tmpCoordinates[255];
	double localDeltaX, localDeltaY;
	double theta0;

	tmpCoordinates[0] = coordinates[0];
	tmpCoordinates[1] = coordinates[1];

	for(int i=0; i < nbCompute; i++){
//		localDeltaX = (exp(param[0])/sqrt(1+(pow(coordinates[0],2.0)+pow(coordinates[1],2.0))*(exp(2*param[0]) - 1))) - 1;
	//	localDeltaY = (exp(param[0])/sqrt(1+(pow(coordinates[0],2.0)+pow(coordinates[1],2.0))*(exp(2*param[0]) - 1))) - 1;
		localDeltaX = -tmpCoordinates[1]*(exp(param[0])/sqrt(1+(pow(coordinates[0],2.0)+pow(coordinates[1],2.0))*(exp(2*param[0]) - 1))) - tmpCoordinates[0];
		localDeltaY = tmpCoordinates[0] *(exp(param[0])/sqrt(1+(pow(coordinates[0],2.0)+pow(coordinates[1],2.0))*(exp(2*param[0]) - 1))) - tmpCoordinates[1];

		//tmpCoordinates[0] += tmpCoordinates[0] * localDeltaX;
		//tmpCoordinates[1] += tmpCoordinates[1] * localDeltaX;
		tmpCoordinates[0] += localDeltaX;
		tmpCoordinates[1] += localDeltaY;
		coordinates[0] = tmpCoordinates[0];
		coordinates[1] = tmpCoordinates[1];
	}

}
