/*
 * Pendule2Oscillators.h
 *
 *  Created on: Apr 6, 2009
 *      Author: hbenamor
 */

#ifndef PENDULE2OSCILLATORS_H_
#define PENDULE2OSCILLATORS_H_

#include "CyclicAttractorSmartPencil.h"

class Pendule2Oscillators: public CyclicAttractor {
public:
	Pendule2Oscillators();
	Pendule2Oscillators(char * name);
	Pendule2Oscillators(double * myParam,char * name);
	virtual ~Pendule2Oscillators();
	void compute(double coordinates[], int nbCompute);
	void resetFather(double coordinates[]);
private:
	void initDefault();
};

#endif /* PENDULE2OSCILLATORS_H_ */
