/*
 * Pendule2Oscillators.cpp
 *
 *  Created on: Apr 6, 2009
 *      Author: hbenamor
 */

#include "Pendule2Oscillators.h"

Pendule2Oscillators::Pendule2Oscillators() : CyclicAttractor("P2Oscillators") {
	// TODO Auto-generated constructor stub
	initDefault();
}

Pendule2Oscillators::Pendule2Oscillators(char * name) : CyclicAttractor(name) {
	// TODO Auto-generated constructor stub
	initDefault();
}

Pendule2Oscillators::Pendule2Oscillators(double * myParam,char * name) : CyclicAttractor(4,4,myParam,name) {
	// TODO Auto-generated constructor stub
}

Pendule2Oscillators::~Pendule2Oscillators() {
	// TODO Auto-generated destructor stub
}

void Pendule2Oscillators::initDefault(){
	nbParam = 4;		// Couplage des Cx, Cy, Cxy, Cyx
	nbDimension = 4;
	coordinates[0] = 0.0001;
	coordinates[1] = 0;
	coordinates[2] = 0.0001;
	coordinates[3] = 0;
	param[0] = 0;
	param[1] = 100;
	param[2] = 0;
	param[3] = 0;
}

void Pendule2Oscillators::resetFather(double coordinates[]){

}


void Pendule2Oscillators::compute(double coordinates[], int nbCompute){
	double tmpCoordinates[255];
	double nPhase;
	tmpCoordinates[0] = coordinates[0];
	tmpCoordinates[1] = coordinates[1];
	tmpCoordinates[2] = coordinates[2];
	tmpCoordinates[3] = coordinates[3];


	for(int i=0; i < nbCompute; i++){
		// Systeme entrainé
		tmpCoordinates[0] += (deltaT * (coordinates[1]+coordinates[0] * (1 - coordinates[0] * coordinates[0] - coordinates[1] * coordinates[1]) + param[0] * coordinates[2] + param[3] * coordinates[3]));
		tmpCoordinates[1] += (deltaT * (-coordinates[0] + coordinates[1] * (1-coordinates[0]*coordinates[0]-coordinates[1]*coordinates[1] + param[1] * coordinates[3] + param[2] * coordinates[2])));
		// SecondOscillator -- Entraineur
		tmpCoordinates[2] += (deltaT * (coordinates[3]+coordinates[2] * (1 - coordinates[2] * coordinates[2] - coordinates[3] * coordinates[3])));
		tmpCoordinates[3] += (deltaT * (-coordinates[2] + coordinates[3] * (1-coordinates[2]*coordinates[2]-coordinates[3]*coordinates[3])));
		// Correction de tir
		// Pour le systeme de base
		nPhase = atan2(tmpCoordinates[1], tmpCoordinates[0]);
		tmpCoordinates[0] = cos(nPhase);
		tmpCoordinates[1] = sin(nPhase);


		coordinates[0] = tmpCoordinates[0];
		coordinates[1] = tmpCoordinates[1];
		coordinates[2] = tmpCoordinates[2];
		coordinates[3] = tmpCoordinates[3];
	}
}

