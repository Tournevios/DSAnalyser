/*
 * TestSampleDS.cpp
 *
 *  Created on: May 4, 2009
 *      Author: hbenamor
 */

#include "TestSampleDS.h"

TestSampleDS::~TestSampleDS() {
	// TODO Auto-generated destructor stub
}


TestSampleDS::TestSampleDS() : CyclicAttractor("TestSampleDS") {
	// TODO Auto-generated constructor stub
	nbParam = 0;
	nbDimension = 2;
	coordinates[0] = 10;
	coordinates[1] = 10;
	coordinates[2] = 0;
}


void TestSampleDS::resetFather(double coordinates[]){
}

void TestSampleDS::compute(double coordinates[], int nbCompute){
	double tmpCoordinates[255];
	tmpCoordinates[0] = coordinates[0];
	tmpCoordinates[1] = coordinates[1];

	for(int i=0; i < nbCompute; i++){
		tmpCoordinates[0] += (deltaT * (-coordinates[0]));
		tmpCoordinates[1] += (deltaT * (-coordinates[1]));

		coordinates[0] = tmpCoordinates[0];
		coordinates[1] = tmpCoordinates[1];
	}

}
