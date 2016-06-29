/*
 * Rabinovich.h
 *
 *  Created on: Mar 26, 2009
 *      Author: hbenamor
 */

#ifndef RABINOVICH_H_
#define RABINOVICH_H_
#include "CyclicAttractorSmartPencil.h"

class Rabinovich : public CyclicAttractor{
public:
	Rabinovich();
	Rabinovich(char * name);
	Rabinovich(double * myParam, char * name);
	void compute(double coordinates[], int nbCompute);
	virtual ~Rabinovich();
	void resetFather(double coordinates[]);
private:
	void initDefault();
};

#endif /* RABINOVICH_H_ */
