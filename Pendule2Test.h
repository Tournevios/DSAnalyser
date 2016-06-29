/*
 * Pendule2Test.h
 *
 *  Created on: Oct 30, 2009
 *      Author: hbenamor
 */

#ifndef PENDULE2TEST_H_
#define PENDULE2TEST_H_

#include "CyclicAttractorSmartPencil.h"

class Pendule2Test: public CyclicAttractor {
public:
	Pendule2Test();
	Pendule2Test(char * name);
	Pendule2Test(double * myParam, char * name);
	virtual ~Pendule2Test();
	void compute(double coordinates[], int nbCompute);
	void resetFather(double coordinates[]);
private:
	void initDefault();
};
#endif /* PENDULE2TEST_H_ */
