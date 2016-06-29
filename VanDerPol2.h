/*
 * VanDerPol2.h
 *
 *  Created on: 19-Jan-2009
 *      Author: hbenamor
 */

#ifndef VANDERPOL2_H_
#define VANDERPOL2_H_

#include "CyclicAttractorSmartPencil.h"

class VanDerPol2: public CyclicAttractor {
public:
	VanDerPol2();
	VanDerPol2(char * name);
	VanDerPol2(double * myParam, char * name);
	virtual ~VanDerPol2();
	void compute(double coordinates[], int nbCompute);
	void resetFather(double coordinates[]);

private:
	void initDefault();
};

#endif /* VANDERPOL2_H_ */
