/*
 * VariationPendule.h
 *
 *  Created on: 24-Feb-2009
 *      Author: hbenamor
 */

#ifndef VARIATIONPENDULE_H_
#define VARIATIONPENDULE_H_

#include "CyclicAttractorSmartPencil.h"

class VariationPendule: public CyclicAttractor {
public:
	VariationPendule();
	VariationPendule(char * name);
	VariationPendule(double * myParam, char * name);
	virtual ~VariationPendule();
	void compute(double coordinates[], int nbCompute);
	void resetFather(double coordinates[]);
private:
	void initDefault();
};


#endif /* VARIATIONPENDULE_H_ */
