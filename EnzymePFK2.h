/*
 * EnzymePFK2.h
 *
 *  Created on: 06-Feb-2009
 *      Author: hbenamor
 */

#ifndef ENZYMEPFK2_H_
#define ENZYMEPFK2_H_

#include "CyclicAttractorSmartPencil.h"

class EnzymePFK2: public CyclicAttractor {
public:
	EnzymePFK2();
	EnzymePFK2(char * name);
	EnzymePFK2(double * myParam,char * name);
	virtual ~EnzymePFK2();
	void compute(double coordinates[], int nbCompute);
	void resetFather(double coordinates[]);
private:
	void initDefault();
};

#endif /* ENZYMEPFK2_H_ */
