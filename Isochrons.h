/*
 * Isochrons.h
 *
 *  Created on: 18-Dec-2008
 *      Author: hbenamor
 */

#ifndef ISOCHRONS_H_
#define ISOCHRONS_H_

#include <fstream>
#include <vector>
#include "point3D.h"

class Isochrons {
public:

	Isochrons();
	Isochrons(double phase, char fpaths[255], int maxInMemory, int myIndex, int myIndexInTheAttractor, point3D * p3D);
	virtual ~Isochrons();
	std::vector<point3D*> points;
	double phase;
	int nbPoints;
	int nbSavedPoints;
	int maxInMemory;
	bool addPoint(point3D* p3D);
	char fpaths[255];
	void saveAll();
	int myIndex;
	int myIndexInTheAttractor;
	point3D * phasePoint;
};

#endif /* ISOCHRONS_H_ */
