/*
 * EnzymePFK2.cpp
 *
 *  Created on: 06-Feb-2009
 *      Author: hbenamor
 */

#include "EnzymePFK2.h"

EnzymePFK2::EnzymePFK2(char * name):CyclicAttractor(name) {
	// TODO Auto-generated constructor stub
	initDefault();
}

void EnzymePFK2::initDefault(){
	nbParam = 9;
	nbDimension = 2;

	param[0] = 3000;  				// L0
	param[1] = 3;					// n
	param[2] = 0.02;				// C
	param[3] = 1.0;					// D
	param[4] = 0.01;				// N
	param[5] = 0.1;					// A
	param[6] = 3.0;					// P
	param[7] = pow(10.0,param[6]);	// R

	coordinates[0] = -0.3561;
	coordinates[1] = -30.9889 ;
	coordinates[2] = 0.0;

	param[8] = coordinates[0] * (pow((1+coordinates[0]),(param[1]-1)) * pow((1+param[3] * coordinates[1]),param[1]) + param[0] * param[2] * pow((1+param[2] * coordinates[0]),(param[1]-1))) / (pow((1+coordinates[0]),param[1]) * pow((1+param[3] * coordinates[1]),param[1]) + param[0] * pow((1+param[2] * coordinates[0]),param[1]));
}

EnzymePFK2::EnzymePFK2():CyclicAttractor("EnzymePFK2") {
	// TODO Auto-generated constructor stub
	initDefault();
}

EnzymePFK2::EnzymePFK2(double * myParam,char * name) : CyclicAttractor(3,9,myParam,name) {
	// TODO Auto-generated constructor stub
	param[8] = coordinates[0] * (pow((1+coordinates[0]),(param[1]-1)) * pow((1+param[3] * coordinates[1]),param[1]) + param[0] * param[2] * pow((1+param[2] * coordinates[0]),(param[1]-1))) / (pow((1+coordinates[0]),param[1]) * pow((1+param[3] * coordinates[1]),param[1]) + param[0] * pow((1+param[2] * coordinates[0]),param[1]));
	coordinates[2] = 0;
}


EnzymePFK2::~EnzymePFK2() {
	// TODO Auto-generated destructor stub
	initDefault();
}

void EnzymePFK2::resetFather(double coordinates[]){

}

void EnzymePFK2::compute(double coordinates[], int nbCompute){
	double tmpCoordinates[255];
	tmpCoordinates[0] = coordinates[0];
	tmpCoordinates[1] = coordinates[1];
	tmpCoordinates[2] = coordinates[2];

	for(int i=0; i < nbCompute; i++){
		tmpCoordinates[0] += (deltaT * (param[5]-param[8]));
		tmpCoordinates[1] += (deltaT * (param[7] *(param[8] - param[4] * coordinates[1])));
		if(not(finite(tmpCoordinates[0])) or not(finite(tmpCoordinates[1])) or not(finite(tmpCoordinates[2]))) break;
		//if((tmpCoordinates[0] > radius[0]) or (tmpCoordinates[1] > radius[1]) or (tmpCoordinates[2] > radius[2])) break;
		coordinates[0] = tmpCoordinates[0];
		coordinates[1] = tmpCoordinates[1];
		param[8] = coordinates[0] * (pow((1+coordinates[0]),(param[1]-1)) * pow((1+param[3] * coordinates[1]),param[1]) + param[0] * param[2] * pow((1+param[2] * coordinates[0]),(param[1]-1))) / (pow((1+coordinates[0]),param[1]) * pow((1+param[3] * coordinates[1]),param[1]) + param[0] * pow((1+param[2] * coordinates[0]),param[1]));
	}
}
