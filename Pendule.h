/*
 * Pendule.h
 *
 *  Created on: 16-Jan-2009
 *      Author: hbenamor
 */

#ifndef PENDULE_H_
#define PENDULE_H_

#include <fstream>
#include "point3D.h"
#include <vector>
#include "random-singleton.h"
#include "Isochrons.h"

class Pendule {
public:
	std::vector<point3D*> attractor;
	std::vector<Isochrons*> isochrons;
	Pendule();
	virtual ~Pendule();

	void tracer(double deltaT, int numberOfIterations);
	void tracerIsochrons(int nbIsochrons, int nbPointsPerIsochrons);
	double phaseOf(point3D * p3d, int numberMaxOfIsoPoints);
	double distancePointSegment(point3D * p, point3D * p1, point3D * p2, bool * inSegment);
private:
	double xradius, yradius, zradius;
	double sigma, rho, beta;
	double mu;
	int numberOfIterations;
	double deltaT;
	double periode;
	double resolution;
	int nbPointsPerPeriode;
};

#endif /* PENDULE_H_ */


