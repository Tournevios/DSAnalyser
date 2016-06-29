/*
 * Pendule2.h
 *
 *  Created on: 19-Jan-2009
 *      Author: hbenamor
 */

#ifndef PENDULE2_H_
#define PENDULE2_H_

#include "CyclicAttractorSmartPencil.h"

class Pendule2: public CyclicAttractor {
public:
	Pendule2();
	Pendule2(char * name);
	Pendule2(double * param,char * name);
	virtual ~Pendule2();
	void compute(double coordinates[], int nbCompute);
	void resetFather(double coordinates[]);
	double coordinates2[255];
	//double coordinates3[255];
private:
	void initDefault();
};

#endif /* PENDULE2_H_ */
