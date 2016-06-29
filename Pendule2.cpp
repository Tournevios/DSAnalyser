/*
 * Pendule2.cpp
 *
 *  Created on: 19-Jan-2009
 *      Author: hbenamor
 */

#include "Pendule2.h"

Pendule2::Pendule2() : CyclicAttractor("Pendule2Couple") {
	// TODO Auto-generated constructor stub
	initDefault();
}

Pendule2::Pendule2(char * name) : CyclicAttractor(name) {
	// TODO Auto-generated constructor stub
	initDefault();
}

Pendule2::Pendule2(double * param,char * name) : CyclicAttractor(2,0,param,name) {
	// TODO Auto-generated constructor stub
	// For the system that evolve in the 4th state space area
	minFirstCoor[1] = -0.1;
	maxFirstCoor[1] = 0.1;

	minFirstCoor[0] = 0.01;
	maxFirstCoor[0] = 65000;

}

void Pendule2::initDefault(){
	nbParam = 0;
	//nbDimension = 2;
	nbDimension = 2;
	coordinates[0] = 1;
	coordinates[1] = 0;
	coordinates[2] = 0;

	// For the system that evolve in the 4th state space area
	minFirstCoor[1] = -0.1;
	maxFirstCoor[1] = 0.1;

	minFirstCoor[0] = 0.01;
	maxFirstCoor[0] = 65000;
}

Pendule2::~Pendule2() {
	// TODO Auto-generated destructor stub
}


void Pendule2::resetFather(double coordinates[]){
/*	double phase = atan2(coordinates[0], coordinates[1]);
	coordinates2[0] = sin(phase);
	coordinates2[1] = cos(phase);

	coordinates3[0] = cos(phase);
	coordinates3[1] = sin(phase/3);*/
}

void Pendule2::compute(double coordinates[], int nbCompute){
	double tmpCoordinates[255];
	/*double tmpCoordinates2[255];
	double tmpCoordinates3[255];
	double nPhase;
	double nPhase2;*/


	//if(nbCompute==0) resetFather(coordinates);
	tmpCoordinates[0] = coordinates[0];
	tmpCoordinates[1] = coordinates[1];
	tmpCoordinates[2] = 0;
	coordinates[2] = 0;
	//tmpCoordinates[2] = coordinates[2];
	/*
	tmpCoordinates2[0] = coordinates2[0];
	tmpCoordinates2[1] = coordinates2[1];
	tmpCoordinates3[0] = coordinates3[0];
	tmpCoordinates3[1] = coordinates3[1];

	double c1x = 0.0001;
	double c1y = 0;
	double c2x = 0;
	double c2y = 0;
	double c3x = 0;
	double c3y = 0;
	double kz = -3;		// relaxation suivant z
*/
	/*
	double myAngle = 0;
	double h;
	double hmin=0.1;
	double hmax = 10;*/

	for(int i=0; i < nbCompute; i++){
		// Systeme de base
		/*h=hmin;
			if(coordinates[1]==0)	{myAngle = 0;h=hmax;}
			else if(coordinates[1] > 0) {
				myAngle = fabs(atan(coordinates[0]/coordinates[1]));
				h=hmin+(hmax-hmin)*((M_PI/2)-myAngle)/(M_PI/2);
			}*/
		tmpCoordinates[0] += deltaT * (coordinates[1]+coordinates[0] * (1 - coordinates[0] * coordinates[0] - coordinates[1] * coordinates[1]));
		tmpCoordinates[1] += deltaT * (-coordinates[0] + coordinates[1] * (1-coordinates[0]*coordinates[0]-coordinates[1]*coordinates[1]));
/*
		tmpCoordinates[0] += deltaT * h*(coordinates[1]+coordinates[0] * (1 - coordinates[0] * coordinates[0] - coordinates[1] * coordinates[1]));
		tmpCoordinates[1] += deltaT * h*(-coordinates[0] + coordinates[1] * (1-coordinates[0]*coordinates[0]-coordinates[1]*coordinates[1]));
*/
/*
		tmpCoordinates[0] += deltaT * (coordinates[1] - coordinates[2]);
		tmpCoordinates[1] += deltaT * (-coordinates[0] + coordinates[2]);
		tmpCoordinates[2] += deltaT * (coordinates[0] - coordinates[1]);
*/
		// Correction de tir
		// Pour le systeme de base
		/*nPhase = atan2(tmpCoordinates[1], tmpCoordinates[0]);
		tmpCoordinates[0] = cos(nPhase);
		tmpCoordinates[1] = sin(nPhase);*/

		// Systeme 2
		/*
		tmpCoordinates2[0] += (deltaT * (coordinates2[1]+coordinates2[0] * (1 - coordinates2[0] * coordinates2[0] - coordinates2[1] * coordinates2[1])) + c2x * coordinates3[0]);
		tmpCoordinates2[1] += (deltaT * (-coordinates2[0] + coordinates2[1] * (1-coordinates2[0]*coordinates2[0]-coordinates2[1]*coordinates2[1]) + c2y * coordinates3[1]));
		 */
		// Correction de tir
		// Systeme 2
		/*nPhase2 = atan2(tmpCoordinates2[1], tmpCoordinates2[0]);
		tmpCoordinates2[0] = cos(nPhase2);
		tmpCoordinates2[1] = sin(nPhase2);*/

		// Systeme meneur
		/*tmpCoordinates3[0] += (deltaT * (coordinates3[1]+coordinates3[0] * (1 - coordinates3[0] * coordinates3[0] - coordinates3[1] * coordinates3[1]) + c3x * coordinates[0]));
		tmpCoordinates3[1] += (deltaT * (-coordinates3[0] + coordinates3[1] * (1-coordinates3[0]*coordinates3[0]-coordinates3[1]*coordinates3[1]))+ c3y * coordinates[1]);
		nPhase = atan2(tmpCoordinates3[1], tmpCoordinates3[0]);
		tmpCoordinates3[0] = cos(nPhase);
		tmpCoordinates3[1] = sin(nPhase);*/

		// systeme duale
		/*tmpCoordinates[0] += (deltaT * (-coordinates[0] + coordinates[1] * (1-coordinates[0]*coordinates[0]-coordinates[1]*coordinates[1])));
		tmpCoordinates[1] += (deltaT * -(coordinates[1]+coordinates[0] * (1 - coordinates[0] * coordinates[0] - coordinates[1] * coordinates[1])));
		*/
		// Systeme mirroir
/*		tmpCoordinates[0] += (deltaT * (-coordinates[1] + coordinates[0] * (1-coordinates[1]*coordinates[1]-coordinates[0]*coordinates[0])));
		tmpCoordinates[1] += (deltaT * (coordinates[0]+coordinates[1] * (1 - coordinates[1] * coordinates[1] - coordinates[0] * coordinates[0])));
		*/
		// Systeme duale mirroir
	/*	tmpCoordinates[0] += (-deltaT * (coordinates[0]+coordinates[1] * (1 - coordinates[1] * coordinates[1] - coordinates[0] * coordinates[0])));
		tmpCoordinates[1] += (deltaT * (-coordinates[1] + coordinates[0] * (1-coordinates[1]*coordinates[1]-coordinates[0]*coordinates[0])));
*/

		coordinates[0] = tmpCoordinates[0];
		coordinates[1] = tmpCoordinates[1];
		coordinates[2] = tmpCoordinates[2];
/*
		coordinates2[0] = tmpCoordinates2[0];
		coordinates2[1] = tmpCoordinates2[1];
*/
		/*coordinates3[0] = tmpCoordinates3[0];
		coordinates3[1] = tmpCoordinates3[1];*/

	}

}
