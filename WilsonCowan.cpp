/*
 * WilsonCowan.cpp
 *
 *  Created on: 29-Jan-2009
 *      Author: hbenamor
 */

#include "WilsonCowan.h"
#include <cmath>

WilsonCowan::WilsonCowan(): CyclicAttractor("WilsonCowan") {
	// TODO Auto-generated constructor stub
	initDefault();
}

WilsonCowan::WilsonCowan(char * name): CyclicAttractor(name) {
	// TODO Auto-generated constructor stub
	initDefault();
}

WilsonCowan::WilsonCowan(double * myParam, char * name): CyclicAttractor(2, 2, myParam, name) {
	// TODO Auto-generated constructor stub
	// For the system that evolve in the 4th state space area
	minFirstCoor[1] = -0.1;
	maxFirstCoor[1] = 0.1;

	minFirstCoor[0] = 0.01;
	maxFirstCoor[0] = 65000;

}

void WilsonCowan::initDefault(){
	nbParam = 3;
	nbDimension = 2;

	param[0] = 0.6; //1; //0.6;
	param[1] = 2;//1.1; //2;
	coordinates[0] = 2; //0.0001;
	coordinates[1] = 2; //0.0002;

	// For the system that evolve in the 4th state space area
	minFirstCoor[1] = -0.1;
	maxFirstCoor[1] = 0.1;

	minFirstCoor[0] = 0.01;
	maxFirstCoor[0] = 65000;

}

void WilsonCowan::resetFather(double coordinates[]){

}

WilsonCowan::~WilsonCowan() {
	// TODO Auto-generated destructor stub
}

void WilsonCowan::compute(double coordinates[], int nbCompute){
	double tmpCoordinates[255];
	tmpCoordinates[0] = coordinates[0];
	tmpCoordinates[1] = coordinates[1];

	for(int i=0; i < nbCompute; i++){
		//tmpCoordinates[0] = tanh(param[0] * coordinates[0] - param[1] * coordinates[1]);
		//tmpCoordinates[1] = tanh(param[1] * coordinates[0] + param[0] * coordinates[1]);
		tmpCoordinates[0] += deltaT * ((-coordinates[0] / param[0]) + tanh(param[1] * coordinates[0]) - tanh(param[1] * coordinates[1]));
		tmpCoordinates[1] += deltaT * ((-coordinates[1] / param[0]) + tanh(param[1] * coordinates[0]) + tanh(param[1] * coordinates[1]));

		coordinates[0] = tmpCoordinates[0];
		coordinates[1] = tmpCoordinates[1];
		//printf("x = %f y=%f\n", coordinates[0], coordinates[1]);
	}
}
