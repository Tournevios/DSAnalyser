/*
 * point3D.h
 *
 *  Created on: 17-Dec-2008
 *      Author: hbenamor
 */

#ifndef POINT3D_H_
#define POINT3D_H_

#include <vector>

class point3D {
public:
	double x,y,z, angle, norme, phase;
	int k;
	int numberOfIsoPoints;
	bool inTheAttractor;
	point3D();
	point3D(double x, double y, double z, double angle, double norme, double phase, int k);
	point3D(double x, double y, double z);
	point3D(const point3D & p3d);
	virtual ~point3D();
	void reInitialize();
	std::vector<point3D*> isoPoints;
	bool operator<(const point3D& p) const;
};
#endif /* POINT3D_H_ */
