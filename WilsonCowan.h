/*
 * WilsonCowan.h
 *
 *  Created on: 29-Jan-2009
 *      Author: hbenamor
 */

#ifndef WILSONCOWAN_H_
#define WILSONCOWAN_H_
#include "CyclicAttractorSmartPencil.h"

class WilsonCowan: public CyclicAttractor {
public:
	WilsonCowan();
	WilsonCowan(char * name);
	WilsonCowan(double * myParam, char * name);
	virtual ~WilsonCowan();
	void compute(double coordinates[], int nbCompute);
	void resetFather(double coordinates[]);
private:
	void initDefault();
};

#endif /* WILSONCOWAN_H_ */
