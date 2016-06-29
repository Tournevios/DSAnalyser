/*
 * YeastGlycosis.h
 *
 *  Created on: Jan 8, 2010
 *      Author: hbenamor
 */

#ifndef YEASTGLYCOSIS_H_
#define YEASTGLYCOSIS_H_
#include "CyclicAttractorSmartPencil.h"

class YeastGlycosis : public CyclicAttractor {
public:
	YeastGlycosis();
	YeastGlycosis(char * name);
	YeastGlycosis(double * myParam, char * name);
	virtual ~YeastGlycosis();
	void compute(double coordinates[], int nbCompute);
	void resetFather(double coordinates[]);
private:
	void initDefault();
};

#endif /* YEASTGLYCOSIS_H_ */
