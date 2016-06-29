/*
 * YeastGlycosis.cpp
 *
 *  Created on: Jan 8, 2010
 *      Author: hbenamor
 */

#include "YeastGlycosis.h"

YeastGlycosis::YeastGlycosis() : CyclicAttractor("YeastGlycosis") {
	// TODO Auto-generated constructor stub
	initDefault();
}

YeastGlycosis::YeastGlycosis(char * name) : CyclicAttractor(name) {
	// TODO Auto-generated constructor stub
	initDefault();
}

YeastGlycosis::YeastGlycosis(double * myParam, char * name) : CyclicAttractor(2,4, myParam, name) {
	// TODO Auto-generated constructor stub
	initDefault();
}

YeastGlycosis::~YeastGlycosis() {
	// TODO Auto-generated destructor stub
}

void YeastGlycosis::initDefault(){
	nbParam = 4;
	nbDimension = 2;

	param[0] = 0.36;
	param[1] = 0.02;
	param[2] = 6.0;
	param[3] = 13.0;
	coordinates[0] = 6.6;
	coordinates[1] = 7.6;
	coordinates[2] = 0;
}

void YeastGlycosis::resetFather(double coordinates[]){

}


void YeastGlycosis::compute(double coordinates[], int nbCompute){
	double tmpCoordinates[255];
	tmpCoordinates[0] = coordinates[0];
	tmpCoordinates[1] = coordinates[1];

	for(int i=0; i < nbCompute; i++){
		tmpCoordinates[0] += deltaT * (param[0] - param[1] * coordinates[0]*coordinates[1]);
		tmpCoordinates[1] += deltaT * (2*param[1] * coordinates[0]*coordinates[1] - (param[2] * coordinates[1]/(coordinates[1] + param[3])));
		coordinates[0] = tmpCoordinates[0];
		coordinates[1] = tmpCoordinates[1];
	}
}
