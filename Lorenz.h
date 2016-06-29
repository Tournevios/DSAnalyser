/*
 * Lorenz.h
 *
 *  Created on: 19-Jan-2009
 *      Author: hbenamor
 */

#ifndef LORENZ_H_
#define LORENZ_H_

#include "CyclicAttractorSmartPencil.h"

class Lorenz: public CyclicAttractor {
public:
	Lorenz();
	Lorenz(char * name);
	Lorenz(double * myParam,char * name);
	virtual ~Lorenz();
	void compute(double coordinates[], int nbCompute);
	void resetFather(double coordinates[]);

	/*double param[255];
	double coordinates[255];
	int nbDimension;
	int nbParam;*/
private:
	void initDefault();
};

#endif /* LORENZ_H_ */
