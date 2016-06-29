/*
 * Isochrons.cpp
 *
 *  Created on: 18-Dec-2008
 *      Author: hbenamor
 */

#include "Isochrons.h"
#include <string.h>
#include <stdio.h>
#include <fstream>

using namespace std;
Isochrons::Isochrons() {
	// TODO Auto-generated constructor stub
	phase = 0;
	nbPoints = 0;
	maxInMemory = 1000;	//10
	nbSavedPoints = 0;
}


Isochrons::Isochrons(double phase, char fpaths[255], int maxInMemory, int myIndex, int myIndexInTheAttractor, point3D * p3D){
	this->phase = phase;
	nbPoints = 0;
	strcpy(this->fpaths,fpaths);
	this->maxInMemory = maxInMemory;
	nbSavedPoints = 0;
	this->myIndex = myIndex;
	this->myIndexInTheAttractor = myIndexInTheAttractor;
	phasePoint = new point3D(p3D->x, p3D->y, p3D->z);
}

Isochrons::~Isochrons() {
	// TODO Auto-generated destructor stub
	while(points.size() > 0){
		delete points[points.size()-1];
		points.pop_back();
	}
	delete phasePoint;
}

bool Isochrons::addPoint(point3D * p3D){
	ofstream out;
	bool found = false;
	int k=0;
	while(not found and k <(int) points.size()){
		if((points[k]->x == p3D->x) and (points[k]->y == p3D->y) and (points[k]->z == p3D->z)){
			found = true;
			break;
		}
		else k++;
	}
	if(not found) {
		points.push_back(new point3D(p3D->x, p3D->y, p3D->z, p3D->angle, p3D->norme, p3D->phase, p3D->k));
		nbPoints++;
		if((int)points.size()>=maxInMemory){
			ofstream out(fpaths, ios::app);
			while(points.size() > 1){
				out << points[0]->x << "\t" << points[0]->y << "\t" << points[0]->z << "\t" << points[0]->angle << "\t" << points[0]->norme << "\t" << points[0]->phase << "\t" << points[0]->k << "\n";
				delete points[0];
				points.erase(points.begin());
				nbSavedPoints++;
			}
			out.close();
		}
		return true;
	}
	return false;
}


void Isochrons::saveAll(){
	ofstream out(fpaths, ios::app);
	while(points.size() > 0){
		out << points[points.size()-1]->x << "\t" << points[points.size()-1]->y << "\t" << points[points.size()-1]->z << "\t" << points[points.size()-1]->angle << "\t" << points[points.size()-1]->norme << "\t" << "\t" << points[points.size()-1]->phase << "\t" << points[points.size()-1]->k << "\n";
		delete points[points.size()-1];
		points.pop_back();
		nbSavedPoints++;
	}
	out.close();
}
