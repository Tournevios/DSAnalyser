/*
 * point3D.cpp
 *
 *  Created on: 17-Dec-2008
 *      Author: hbenamor
 */
#include "point3D.h"

point3D::point3D() {
	// TODO Auto-generated constructor stub
	x = 0;
	y = 0;
	z = 0;
	angle = 0;
	phase = 0;
	norme = 0;
	k = 0;
	inTheAttractor = false;
	numberOfIsoPoints = 0;

}

point3D::point3D(double x, double y, double z, double angle, double norme, double phase, int k){
	this->x = x;
	this->y = y;
	this->z = z;
	this->phase = phase;
	this->norme = norme;
	this->angle = angle;
	this->k = k;
	inTheAttractor = false;
	numberOfIsoPoints = 0;
}

point3D::point3D(double x, double y, double z){
	this->x = x;
	this->y = y;
	this->z = z;
	this->phase = 0;
	this->norme = 0;
	this->angle = 0;
	this->k = 0;
	inTheAttractor = false;
	numberOfIsoPoints = 0;
}

point3D::~point3D() {
	// TODO Auto-generated destructor stub
	while(isoPoints.size() > 0) {
		delete isoPoints[isoPoints.size()-1];
		isoPoints.pop_back();
	}
}

point3D::point3D(const point3D & p3d){
	this->x=p3d.x;
	this->y=p3d.y;
	this->z=p3d.z;
	this->norme=p3d.norme;
	this->angle=p3d.angle;
	this->phase=p3d.phase;
	this->k=p3d.k;
	this->inTheAttractor=p3d.inTheAttractor;
	this->numberOfIsoPoints=p3d.numberOfIsoPoints;
}

void point3D::reInitialize(){
	while(isoPoints.size() > 0){
		delete isoPoints[isoPoints.size()-1];
		isoPoints.pop_back();
	}
	x = 0;
	y = 0;
	z = 0;
	phase = 0;
	norme = 0;
	angle = 0;
	k = 0;
	inTheAttractor = false;
	numberOfIsoPoints = 0;
}

bool point3D::operator<(const point3D& p2) const{
	return phase < p2.phase;
}
