/*
 * Pendule2Test.cpp
 *
 *  Created on: Oct 30, 2009
 *      Author: hbenamor
 */

#include "Pendule2Test.h"

Pendule2Test::Pendule2Test() : CyclicAttractor("Pendule2Test") {
	// TODO Auto-generated constructor stub
	initDefault();
}

Pendule2Test::Pendule2Test(char * name) : CyclicAttractor(name) {
	// TODO Auto-generated constructor stub
	initDefault();
}

Pendule2Test::Pendule2Test(double * myParam, char * name) : CyclicAttractor(2, 0, myParam, name) {
	// TODO Auto-generated constructor stub
}

Pendule2Test::~Pendule2Test() {
	// TODO Auto-generated destructor stub
}

void Pendule2Test::initDefault() {
	// TODO Auto-generated destructor stub
	nbParam = 0;
	nbDimension = 2;
	coordinates[0] = 1;
	coordinates[1] = 0;
}

void Pendule2Test::resetFather(double coordinates[]){
	double phase = atan2(coordinates[0], coordinates[1]);
}

void Pendule2Test::compute(double coordinates[], int nbCompute){
	double tmpCoordinates[255];

	if(nbCompute==0) resetFather(coordinates);
	tmpCoordinates[0] = coordinates[0];
	tmpCoordinates[1] = coordinates[1];
	tmpCoordinates[2] = 0;
	tmpCoordinates[3] = 0;




	for(int i=0; i < nbCompute; i++){
		// Systeme de base
		// Evaluation de la dérivé seconde en x et en y
		tmpCoordinates[2] = -(tmpCoordinates[1]+tmpCoordinates[0] * (1 - tmpCoordinates[0] * tmpCoordinates[0] - tmpCoordinates[1] * tmpCoordinates[1]));
		tmpCoordinates[3] = -(-tmpCoordinates[0] + tmpCoordinates[1] * (1-tmpCoordinates[0]*tmpCoordinates[0]-tmpCoordinates[1]*tmpCoordinates[1]));

		tmpCoordinates[0] += (deltaT * (coordinates[1]+coordinates[0] * (1 - coordinates[0] * coordinates[0] - coordinates[1] * coordinates[1])));
		tmpCoordinates[1] += (deltaT * (-coordinates[0] + coordinates[1] * (1-coordinates[0]*coordinates[0]-coordinates[1]*coordinates[1])));

		tmpCoordinates[2] += tmpCoordinates[1]+tmpCoordinates[0] * (1 - tmpCoordinates[0] * tmpCoordinates[0] - tmpCoordinates[1] * tmpCoordinates[1]);
		tmpCoordinates[3] += -tmpCoordinates[0] + tmpCoordinates[1] * (1-tmpCoordinates[0]*tmpCoordinates[0]-tmpCoordinates[1]*tmpCoordinates[1]);

		tmpCoordinates[2] = tmpCoordinates[2]/deltaT;
		tmpCoordinates[3] = tmpCoordinates[3]/deltaT;
		tmpCoordinates[4] = tmpCoordinates[2] * tmpCoordinates[2] + tmpCoordinates[3] * tmpCoordinates[3];
		tmpCoordinates[4] = sqrt(tmpCoordinates[4]);

		coordinates[0] = tmpCoordinates[0];
		coordinates[1] = tmpCoordinates[1];
		coordinates[2] = tmpCoordinates[4];

		/*coordinates3[0] = tmpCoordinates3[0];
		coordinates3[1] = tmpCoordinates3[1];*/

	}

}
