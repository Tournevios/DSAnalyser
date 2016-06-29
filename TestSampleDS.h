/*
 * TestSampleDS.h
 *
 *  Created on: May 4, 2009
 *      Author: hbenamor
 */

#ifndef TESTSAMPLEDS_H_
#define TESTSAMPLEDS_H_
#include "CyclicAttractorSmartPencil.h"
class TestSampleDS : public CyclicAttractor {
public:
	TestSampleDS();
	virtual ~TestSampleDS();
	void compute(double coordinates[], int nbCompute);
	void resetFather(double coordinates[]);

};

#endif /* TESTSAMPLEDS_H_ */
