/*
 * Lorenz.cpp
 *
 *  Created on: 19-Jan-2009
 *      Author: hbenamor
 */

#include "Lorenz.h"

Lorenz::Lorenz(double * myParam,char * name) : CyclicAttractor(3,3,myParam,name) {
	// TODO Auto-generated constructor stub
}

Lorenz::Lorenz(char * name) : CyclicAttractor(name){
	initDefault();
}

void Lorenz::initDefault(){
	nbParam = 3;
	nbDimension = 3;
	// 8 boucles
/*	param[0] = 10.00;
	param[1] = 194;
	param[2] = 2.666666;
*/

	// 4 boucles
/*
	param[0] = 10.00;
	param[1] = 195;
	param[2] = 2.666666;
*/

	// chiée de boucles
/*
	param[0] = 10.00;
	param[1] = 190;
	param[2] = 2.666666;
*/


	// 2 boucles
	/*
	param[0] = 10.00;
	param[1] = 200;
	param[2] = 2.666666;
*/

  	// 1 boucle

	param[0] = 10.00;
	param[1] = 350;
	param[2] = 2.66666;



 	// chaotique
	/*param[0] = 10.00;
	param[1] = 170;
	param[2] = 2.66666;*/

	coordinates[0] = 1;
	coordinates[1] = 1.2;
	coordinates[2] = 1.1;
}

Lorenz::Lorenz() : CyclicAttractor("Lorenz"){
	// TODO Auto-generated constructor stub
	initDefault();
}

void Lorenz::resetFather(double coordinates[]){

}

void Lorenz::compute(double * coordinates, int nbCompute){
	double tmpCoordinates[255];
	tmpCoordinates[0] = coordinates[0];
	tmpCoordinates[1] = coordinates[1];
	tmpCoordinates[2] = coordinates[2];

	for(int i=0; i < nbCompute; i++){
		tmpCoordinates[0] += (deltaT * param[0] * (coordinates[1]-coordinates[0]));
		tmpCoordinates[1] += (deltaT * (param[1] * coordinates[0] - coordinates[1] - coordinates[0] * coordinates[2]));
		tmpCoordinates[2] += (deltaT * (coordinates[0] * coordinates[1] - param[2] * coordinates[2]));
		//if((tmpCoordinates[0] > radius[0]) or (tmpCoordinates[1] > radius[1]) or (tmpCoordinates[2] > radius[2])) break;
		coordinates[0] = tmpCoordinates[0];
		coordinates[1] = tmpCoordinates[1];
		coordinates[2] = tmpCoordinates[2];
		if(!finite(coordinates[0]) or !finite(coordinates[1]) or !finite(coordinates[2])) break;
	}
}

Lorenz::~Lorenz() {
	// TODO Auto-generated destructor stub
}
