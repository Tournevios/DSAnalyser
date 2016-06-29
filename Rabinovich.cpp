/*
 * Rabinovich.cpp
 *
 *  Created on: Mar 26, 2009
 *      Author: hbenamor
 */

#include "Rabinovich.h"

Rabinovich::~Rabinovich() {
	// TODO Auto-generated destructor stub
}

Rabinovich::Rabinovich() : CyclicAttractor("Rabinovich"){
	// TODO Auto-generated constructor stub
	initDefault();
}

Rabinovich::Rabinovich(char * name) : CyclicAttractor(name){
	// TODO Auto-generated constructor stub
	initDefault();
}

Rabinovich::Rabinovich(double * myParam, char * name) : CyclicAttractor(3, 4, myParam, name){
	// TODO Auto-generated constructor stub
}

void Rabinovich::initDefault(){
	nbParam = 4;
	nbDimension = 3;

	param[0] = 1;
	param[1] = 1;
	param[2] = 1;
	param[3] = 1;

	coordinates[0] = 1;
	coordinates[1] = 1.2;
	coordinates[2] = 1.1;
}

void Rabinovich::resetFather(double coordinates[]){

}

void Rabinovich::compute(double coordinates[], int nbCompute){
	double tmpCoordinates[255];
	tmpCoordinates[0] = coordinates[0];
	tmpCoordinates[1] = coordinates[1];
	tmpCoordinates[2] = coordinates[2];

	for(int i=0; i < nbCompute; i++){
		tmpCoordinates[0] += (deltaT * (param[0] * coordinates[1]-param[1]*coordinates[0]-coordinates[1]*coordinates[2]));
		tmpCoordinates[1] += (deltaT * (param[0] * coordinates[0] - param[2] * coordinates[1] + param[0] * param[1]));
		tmpCoordinates[2] += (deltaT * (-param[3] * coordinates[2] + coordinates[0] * coordinates[1]));
		//if((tmpCoordinates[0] > radius[0]) or (tmpCoordinates[1] > radius[1]) or (tmpCoordinates[2] > radius[2])) break;
		coordinates[0] = tmpCoordinates[0];
		coordinates[1] = tmpCoordinates[1];
		coordinates[2] = tmpCoordinates[2];
		if(!finite(coordinates[0]) or !finite(coordinates[1]) or !finite(coordinates[2])) break;
	}
}
