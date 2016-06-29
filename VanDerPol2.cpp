/*
 * VanDerPol2.cpp
 *
 *  Created on: 19-Jan-2009
 *      Author: hbenamor
 */

#include "VanDerPol2.h"

VanDerPol2::VanDerPol2(): CyclicAttractor("VanDerPol2") {
	// TODO Auto-generated constructor stub
	initDefault();
}

VanDerPol2::VanDerPol2(char * name): CyclicAttractor(name) {
	// TODO Auto-generated constructor stub
	initDefault();
}

VanDerPol2::VanDerPol2(double * myParam, char * name): CyclicAttractor(2, 1, myParam, name) {
	// TODO Auto-generated constructor stub
	// For the system that evolve in the 4th state space area
	minFirstCoor[1] = -0.1;
	maxFirstCoor[1] = 0.1;

	minFirstCoor[0] = 0.01;
	maxFirstCoor[0] = 65000;

}

VanDerPol2::~VanDerPol2() {
	// TODO Auto-generated destructor stub
}

void VanDerPol2::initDefault(){
	nbParam = 1;
	nbDimension = 2;
	param[0] = 0.01; //0.01; //2
	coordinates[0] = 2.03;
	coordinates[1] = -0.02;
	coordinates[2] = 0;

	// For the system that evolve in the 4th state space area
	minFirstCoor[1] = -0.1;
	maxFirstCoor[1] = 0.01;

	minFirstCoor[0] = 0.01;
	maxFirstCoor[0] = 65000;

}

void VanDerPol2::resetFather(double coordinates[]){

}

void VanDerPol2::compute(double coordinates[], int nbCompute){
	double tmpCoordinates[255];
	tmpCoordinates[0] = coordinates[0];
	tmpCoordinates[1] = coordinates[1];
	double norme;

	for(int i=0; i < nbCompute; i++){
		// systeme initial
//		tmpCoordinates[0] += (deltaT * (coordinates[1]));
//		tmpCoordinates[1] += (deltaT * (-coordinates[0] + coordinates[1] * param[0] * (1 -  coordinates[0] * coordinates[0])));

		tmpCoordinates[0] += (deltaT * param[0] * (coordinates[0]*(1-(coordinates[0]*coordinates[0]/3))- coordinates[1]));
		tmpCoordinates[1] += (deltaT * (coordinates[0] / param[0]));

		// systeme initial normé
		//norme = sqrt(pow(coordinates[1],2.0)+pow(-coordinates[0] + coordinates[1] * param[0] * (1 -  coordinates[0] * coordinates[0]),2.0));
		/*tmpCoordinates[0] += (deltaT * norme *(coordinates[1]));
		tmpCoordinates[1] += (deltaT * norme * (-coordinates[0] + coordinates[1] * param[0] * (1 -  coordinates[0] * coordinates[0])));*/

		// Systeme duale
		/*tmpCoordinates[0] += (deltaT * (-coordinates[0] + coordinates[1] * param[0] * (1 -  coordinates[0] * coordinates[0])));
		tmpCoordinates[1] += (deltaT * -(coordinates[1]));*/

		// Partie potentielle
		/*tmpCoordinates[0] += (deltaT * (3*coordinates[0]-0.75*pow(coordinates[0],3.0)+0.25 * coordinates[0] * pow(coordinates[1],2.0)));
		tmpCoordinates[1] += (deltaT * (-coordinates[1] + pow(coordinates[0],2.0)*coordinates[1]*0.25 - pow(coordinates[1],3.0)/12.0));*/

		// Partie Hamiltonienne
		/*tmpCoordinates[0] += (deltaT * (coordinates[1]-3*coordinates[0] + 0.75 * pow(coordinates[0],3.0)-0.25 * coordinates[0] * pow(coordinates[1],2.0)));
		tmpCoordinates[1] += (deltaT * -(coordinates[0] -3*coordinates[1] + 2.25 * pow(coordinates[0],2.0)*coordinates[1] - pow(coordinates[1],3.0)/12.0));*/

		//if(!finite(tmpCoordinates[0]) or !finite(tmpCoordinates[1])) return;
		coordinates[0] = tmpCoordinates[0];
		coordinates[1] = tmpCoordinates[1];
		if(!finite(coordinates[0]) or !finite(coordinates[1])) return;
	}

}
