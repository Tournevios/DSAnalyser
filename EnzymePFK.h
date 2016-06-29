/*
 * EnzymePFK.h
 *
 *  Created on: 16-Jan-2009
 *      Author: hbenamor
 */

#ifndef ENZYMEPFK_H_
#define ENZYMEPFK_H_

#include <fstream>
#include "point3D.h"
#include <vector>
#include "random-singleton.h"
#include "Isochrons.h"

class EnzymePFK {
public:
	std::vector<point3D*> attractor;
	std::vector<Isochrons*> isochrons;
	EnzymePFK();
	virtual ~EnzymePFK();

	void tracer(double deltaT, int numberOfIterations);
	void tracerIsochrons(int nbIsochrons, int nbPointsPerIsochrons);
	double phaseOf(point3D * p3d, int numberMaxOfIsoPoints);
	double distancePointSegment(point3D * p, point3D * p1, point3D * p2, bool * inSegment);

private:
	double xradius, yradius, zradius;
	int numberOfIterations;
	double deltaT;
	double periode;
	double resolution;
	int nbPointsPerPeriode;

	double L0, n, C, D, N, A, P, R, L;
};
#endif /* ENZYMEPFK_H_ */
