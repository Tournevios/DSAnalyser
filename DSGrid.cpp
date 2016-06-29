/*
 * DSGrid.cpp
 *  Created on: Jan 12, 2010
 *      Author: Hedi BEN AMOR
 *      Adress: TIMC-IMAG Laboratory
 */

#include "DSGrid.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include "Bitmap.h"
#include "random-singleton.h"
#include <vector>
#include <algorithm>

#define PI 3.14159265
using namespace std;

DSGrid::DSGrid(int height, int width, CyclicAttractor * cyclicAttractor) {
	// TODO Auto-generated constructor stub
	this->height = height;
	this->width = width;
	this->cyclicAttractor = cyclicAttractor;
	char signalsPath[255];
/*
	unsigned int area= width*height;
	couplages = new double[area][area][4];//[MAX_INT*MAX_INT][MAX_INT*MAX_INT][4]; // (xi xj) (xi yj) (yi xj) (yi yj) Couplage de j sur i
	states = new double[height][width][3];//[MAX_INT][MAX_INT][3];
	tmpStates = new double[height][width][3];//[MAX_INT][MAX_INT][3];
	seuils = new double[height][width];//[MAX_INT][MAX_INT];
*/
	sprintf(this->path, "%s/DSGrid", cyclicAttractor->rootPath);
	int status;
	status = mkdir(this->path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if(status!=0) printf("Directory %s already exists\n", this->path);

	sprintf(signalsPath, "%s/DSGrid/Signals", cyclicAttractor->rootPath);
	status = mkdir(signalsPath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if(status!=0) printf("Directory %s already exists\n", signalsPath);


	perturbate_randomly();
	//setTheSameStateForAll(1, 1, 0);

	/*
	states[0][0][0] = 0.6;
	states[0][0][1] = 0;
	states[0][0][2] = 0;
	states[0][1][0] = -0.6;
	states[0][1][1] = 0;
	states[0][1][2] = 0;
*/
/*
	states[10][10][0] = 0.5;
	states[10][10][1] = 0.5;
	states[10][10][2] = 0;
	this->cyclicAttractor->computeMethod(states[10][10],cyclicAttractor->numberOfIterations);
*/
	//buildCoupling_all(0);
	initialCoupling = 0.00000001;
	buildCoupling_Moore(initialCoupling);
	//buildCoupling_VonNeumann(0.1);
	//buildCoupling_Chain(initialCoupling);
	//buildCoupling_PerceptronOneLayerToOneOutput(initialCoupling);
	//buildSeuils(10);
	/*
	for(int j=0; j < width*height;j++){
			couplages[154][j][0] = 0;
			couplages[287][j][0] = 0;
	}
	*/
	//states[7][14][0] += 5;
	//states[14][7][0] += 5;
	Xmax = 0.95;
	PXLearningZone = 4;
	alpha = 10000;
	filterCoeff = 0.7777;
}
/*
 * Perturbate to Evoke the image: only 25% of the image is evoked
 */
void DSGrid::perturbate_evokeImage(char * imagePath, double stimulationStrength){
	// Read picture from a file
	Bitmap * bitmapToBurn = new Bitmap(imagePath);

	// Adjust width and height to the minimal values
	// (dimension of the area of intersection of the image with the grid)
	int iWidth = bitmapToBurn->width;
	int iHeight = bitmapToBurn->height;
	if(iWidth > this->width) iWidth = width;
	if(iHeight > this->height) iHeight = height;

	for(int i=0; i < iHeight/2; i++)
		for(int j=0; j < iWidth/2; j++){
			mapPerturbation(states[i][j], (int)(unsigned char)bitmapToBurn->data[(iHeight-1-i)*iWidth*3+j*3], stimulationStrength);
			//printf("test=%d\n",(int)(unsigned char)bitmapToBurn->data[(iHeight-i)*iWidth*3+j]);
		}
	// Convert color to gray intensity
	// USE: Y = 0.3*R + 0.59*G + 0.11*B


	// Adapt to the dimension of the grid

	// Map the gray value to instant perturbations
}

/*
 * Print a picture on a grid : filter the image to gray scale
 * adapt the picture dimension to the grid dimension
 * map the gray value of the adapted picture to the phase to follow in
 * the phase space.
 * Developpement begin the 5th of July 2010 by Hedi Ben Amor
 */
void DSGrid::perturbate_burnImage(char * imagePath, double stimulationStrength){
	// Read picture from a file
	Bitmap * bitmapToBurn = new Bitmap(imagePath);

	// Adjust width and height to the minimal values
	// (dimension of the area of intersection of the image with the grid)
	int iWidth = bitmapToBurn->width;
	int iHeight = bitmapToBurn->height;
	if(iWidth > this->width) iWidth = width;
	if(iHeight > this->height) iHeight = height;

	for(int i=0; i < iHeight; i++)
		for(int j=0; j < iWidth; j++){
			mapPerturbation(states[i][j], (int)(unsigned char)bitmapToBurn->data[(iHeight-1-i)*iWidth*3+j*3], stimulationStrength);
			//printf("test=%d\n",(int)(unsigned char)bitmapToBurn->data[(iHeight-i)*iWidth*3+j]);
		}
	// Convert color to gray intensity
	// USE: Y = 0.3*R + 0.59*G + 0.11*B


	// Adapt to the dimension of the grid

	// Map the gray value to instant perturbations
}

void DSGrid::mapPerturbation(double * state, int grayValue, double stimulationStrength){
	double directionOfPerturbation = ((double)grayValue / 256 ) * M_PI;
	*state += cos(directionOfPerturbation)*stimulationStrength;
	*(state+1) += sin(directionOfPerturbation)*stimulationStrength;
}

/*
 * Apply a stimulation following an H form
 * So following the first (i=0) and last colum (i=width-1) and the middle (j=height/2) line
 */
void DSGrid::perturbate_burnH(double stimulationStrength, double angle){
	// Left line of H
	for(int i=0; i < height; i++){
		states[i][0][0] += stimulationStrength * cos(angle);
		states[i][1][0] += stimulationStrength * cos(angle);
		states[i][2][0] += stimulationStrength * cos(angle);

		states[i][0][1] += stimulationStrength * sin(angle);
		states[i][1][1] += stimulationStrength * sin(angle);
		states[i][2][1] += stimulationStrength * sin(angle);
	}
	// Right line of H
	for(int i=0; i < height; i++){
			states[i][width-1][0] += stimulationStrength * cos(angle);
			states[i][width-2][0] += stimulationStrength * cos(angle);
			states[i][width-3][0] += stimulationStrength * cos(angle);

			states[i][width-1][1] += stimulationStrength * sin(angle);
			states[i][width-2][1] += stimulationStrength * sin(angle);
			states[i][width-3][1] += stimulationStrength * sin(angle);
	}
	// Middle line of H
	for(int j=3; j < width-3; j++){
			states[(height/2)-1][j][0] += stimulationStrength * cos(angle);
			states[height/2][j][0] += stimulationStrength * cos(angle);
			states[(height/2)+1][j][0] += stimulationStrength * cos(angle);

			states[(height/2)-1][j][1] += stimulationStrength * sin(angle);
			states[height/2][j][1] += stimulationStrength * sin(angle);
			states[(height/2)+1][j][1] += stimulationStrength * sin(angle);
	}
}


/*
 * Apply a stimulation following an 8 form
 * So following the first (i=0,1,2) and last colum (i=width-1) and the up, middle and down (j=0,1,2, height-1,height-2,height-3) line
 */
void DSGrid::perturbate_burn8(double stimulationStrength, double angle){
	// Left line of 8
	for(int i=0; i < height; i++){
		states[i][0][0] += stimulationStrength * cos(angle);
		states[i][1][0] += stimulationStrength * cos(angle);
		states[i][2][0] += stimulationStrength * cos(angle);

		states[i][0][1] += stimulationStrength * sin(angle);
		states[i][1][1] += stimulationStrength * sin(angle);
		states[i][2][1] += stimulationStrength * sin(angle);
	}
	// Right line of 8
	for(int i=0; i < height; i++){
			states[i][width-1][0] += stimulationStrength * cos(angle);
			states[i][width-2][0] += stimulationStrength * cos(angle);
			states[i][width-3][0] += stimulationStrength * cos(angle);

			states[i][width-1][1] += stimulationStrength * sin(angle);
			states[i][width-2][1] += stimulationStrength * sin(angle);
			states[i][width-3][1] += stimulationStrength * sin(angle);
	}

	// Up line of 8
	for(int j=3; j < width-3; j++){
				states[0][j][0] += stimulationStrength * cos(angle);
				states[1][j][0] += stimulationStrength * cos(angle);
				states[2][j][0] += stimulationStrength * cos(angle);

				states[0][j][1] += stimulationStrength * sin(angle);
				states[1][j][1] += stimulationStrength * sin(angle);
				states[2][j][1] += stimulationStrength * sin(angle);
	}

	// Middle line of 8
	for(int j=3; j < width-3; j++){
			states[(height/2)-1][j][0] += stimulationStrength * cos(angle);
			states[height/2][j][0] += stimulationStrength * cos(angle);
			states[(height/2)+1][j][0] += stimulationStrength * cos(angle);

			states[(height/2)-1][j][1] += stimulationStrength * sin(angle);
			states[height/2][j][1] += stimulationStrength * sin(angle);
			states[(height/2)+1][j][1] += stimulationStrength * sin(angle);
	}
	// Down line of 8
	for(int j=3; j < width-3; j++){
			states[height-1][j][0] += stimulationStrength * cos(angle);
			states[height-2][j][0] += stimulationStrength * cos(angle);
			states[height-3][j][0] += stimulationStrength * cos(angle);

			states[height-1][j][1] += stimulationStrength * sin(angle);
			states[height-2][j][1] += stimulationStrength * sin(angle);
			states[height-3][j][1] += stimulationStrength * sin(angle);
	}
}

/*
 * Apply a stimulation following an 1 form
 * So following the first (i=0,1,2) and last colum (i=width-1) and the up, middle and down (j=0,1,2, height-1,height-2,height-3) line
 */
void DSGrid::perturbate_burn1(double stimulationStrength, double angle){
	// Right line of 1
	int position = width/2;
	for(int i=0; i < height; i++){
			states[i][position-1][0] += stimulationStrength * cos(angle);
			states[i][position-2][0] += stimulationStrength * cos(angle);
			states[i][position-3][0] += stimulationStrength * cos(angle);

			states[i][position-1][1] += stimulationStrength * sin(angle);
			states[i][position-2][1] += stimulationStrength * sin(angle);
			states[i][position-3][1] += stimulationStrength * sin(angle);
	}
}


/*
 * Apply a stimulation following an 1 form
 * So following the first (i=0,1,2) and last colum (i=width-1) and the up, middle and down (j=0,1,2, height-1,height-2,height-3) line
 */
/*
void DSGrid::perturbate_FromPicture(double stimulationStrength, double angle){
	// Right line of 1
	int position = width/2;
	for(int i=0; i < height; i++){
			states[i][position-1][0] += stimulationStrength * cos(angle);
			states[i][position-2][0] += stimulationStrength * cos(angle);
			states[i][position-3][0] += stimulationStrength * cos(angle);

			states[i][position-1][1] += stimulationStrength * sin(angle);
			states[i][position-2][1] += stimulationStrength * sin(angle);
			states[i][position-3][1] += stimulationStrength * sin(angle);
	}
}
*/

void DSGrid::performLearning(){
	// perform learning
	for(int i=0; i < height; i++){
		for(int j=0; j < width; j++){
			learningFromNeighbourhood(i, j);
		}
	}
}

void DSGrid::perturbate_toEvoke1(double stimulationStrength,double angle){
	int position = width/2;
	for(int i=0; i < 3; i++){
				states[i][position-1][0] += stimulationStrength * cos(angle);
				states[i][position-2][0] += stimulationStrength * cos(angle);
				states[i][position-3][0] += stimulationStrength * cos(angle);

				states[i][position-1][1] += stimulationStrength * sin(angle);
				states[i][position-2][1] += stimulationStrength * sin(angle);
				states[i][position-3][1] += stimulationStrength * sin(angle);
		}
}

/*
 * Apply a stimulation following an 2 form
 * So following the first (i=0,1,2) and last colum (i=width-1) and the up, middle and down (j=0,1,2, height-1,height-2,height-3) line
 */
void DSGrid::perturbate_burn2(double stimulationStrength, double angle){
	// Left line of 2
	for(int i= (height/2)+2; i < height-3; i++){
		states[i][0][0] += stimulationStrength * cos(angle);
		states[i][1][0] += stimulationStrength * cos(angle);
		states[i][2][0] += stimulationStrength * cos(angle);

		states[i][0][1] += stimulationStrength * sin(angle);
		states[i][1][1] += stimulationStrength * sin(angle);
		states[i][2][1] += stimulationStrength * sin(angle);
	}
	// Right line of 2
	for(int i=3; i < (height/2)-1; i++){
			states[i][width-1][0] += stimulationStrength * cos(angle);
			states[i][width-2][0] += stimulationStrength * cos(angle);
			states[i][width-3][0] += stimulationStrength * cos(angle);

			states[i][width-1][1] += stimulationStrength * sin(angle);
			states[i][width-2][1] += stimulationStrength * sin(angle);
			states[i][width-3][1] += stimulationStrength * sin(angle);
	}

	// Up line of 2
	for(int j=0; j < width; j++){
				states[0][j][0] += stimulationStrength * cos(angle);
				states[1][j][0] += stimulationStrength * cos(angle);
				states[2][j][0] += stimulationStrength * cos(angle);

				states[0][j][1] += stimulationStrength * sin(angle);
				states[1][j][1] += stimulationStrength * sin(angle);
				states[2][j][1] += stimulationStrength * sin(angle);
	}

	// Middle line of 2
	for(int j=0; j < width; j++){
			states[(height/2)-1][j][0] += stimulationStrength * cos(angle);
			states[height/2][j][0] += stimulationStrength * cos(angle);
			states[(height/2)+1][j][0] += stimulationStrength * cos(angle);

			states[(height/2)-1][j][1] += stimulationStrength * sin(angle);
			states[height/2][j][1] += stimulationStrength * sin(angle);
			states[(height/2)+1][j][1] += stimulationStrength * sin(angle);
	}
	// Down line of 2
	for(int j=0; j < width; j++){
			states[height-1][j][0] += stimulationStrength * cos(angle);
			states[height-2][j][0] += stimulationStrength * cos(angle);
			states[height-3][j][0] += stimulationStrength * cos(angle);

			states[height-1][j][1] += stimulationStrength * sin(angle);
			states[height-2][j][1] += stimulationStrength * sin(angle);
			states[height-3][j][1] += stimulationStrength * sin(angle);
	}
}


/*
 * Apply a stimulation following an D (or a square)
 */
void DSGrid::perturbate_burnD(double stimulationStrength, double angle){
	// Left column of D
	for(int i=0; i < height; i++){
		states[i][0][0] += stimulationStrength * cos(angle);
		states[i][1][0] += stimulationStrength * cos(angle);
		states[i][2][0] += stimulationStrength * cos(angle);

		states[i][0][1] += stimulationStrength * sin(angle);
		states[i][1][1] += stimulationStrength * sin(angle);
		states[i][2][1] += stimulationStrength * sin(angle);
	}
	// Right column of D
	for(int i=0; i < height; i++){
			states[i][width-1][0] += stimulationStrength * cos(angle);
			states[i][width-2][0] += stimulationStrength * cos(angle);
			states[i][width-3][0] += stimulationStrength * cos(angle);

			states[i][width-1][1] += stimulationStrength * sin(angle);
			states[i][width-2][1] += stimulationStrength * sin(angle);
			states[i][width-3][1] += stimulationStrength * sin(angle);
	}

	// Up line of D
	for(int j=3; j < width-3; j++){
				states[0][j][0] += stimulationStrength * cos(angle);
				states[1][j][0] += stimulationStrength * cos(angle);
				states[2][j][0] += stimulationStrength * cos(angle);

				states[0][j][1] += stimulationStrength * sin(angle);
				states[1][j][1] += stimulationStrength * sin(angle);
				states[2][j][1] += stimulationStrength * sin(angle);
	}

	// Down line of D
	for(int j=3; j < width-3; j++){
			states[height-1][j][0] += stimulationStrength * cos(angle);
			states[height-2][j][0] += stimulationStrength * cos(angle);
			states[height-3][j][0] += stimulationStrength * cos(angle);

			states[height-1][j][1] += stimulationStrength * sin(angle);
			states[height-2][j][1] += stimulationStrength * sin(angle);
			states[height-3][j][1] += stimulationStrength * sin(angle);
	}
}


/*
 * Apply a stimulation following a C form
 * So following the first column (i=0,1,2) and the up and down line (j=0,1,2, height-1,height-2,height-3) line
 */
void DSGrid::perturbate_burnC(double stimulationStrength, double angle){
	// Left line of C
	for(int i=0; i < height; i++){
		states[i][0][0] += stimulationStrength * cos(angle);
		states[i][1][0] += stimulationStrength * cos(angle);
		states[i][2][0] += stimulationStrength * cos(angle);

		states[i][0][1] += stimulationStrength * sin(angle);
		states[i][1][1] += stimulationStrength * sin(angle);
		states[i][2][1] += stimulationStrength * sin(angle);
	}
	// Up and down line of C
	for(int j=3; j < width; j++){
			states[0][j][0] += stimulationStrength * cos(angle);
			states[1][j][0] += stimulationStrength * cos(angle);
			states[2][j][0] += stimulationStrength * cos(angle);

			states[0][j][1] += stimulationStrength * sin(angle);
			states[1][j][1] += stimulationStrength * sin(angle);
			states[2][j][1] += stimulationStrength * sin(angle);

			states[height-1][j][0] += stimulationStrength * cos(angle);
			states[height-2][j][0] += stimulationStrength * cos(angle);
			states[height-3][j][0] += stimulationStrength * cos(angle);

			states[height-1][j][1] += stimulationStrength * sin(angle);
			states[height-2][j][1] += stimulationStrength * sin(angle);
			states[height-3][j][1] += stimulationStrength * sin(angle);
	}
}


void DSGrid::traceMaximumPhaseShift(double stimulationStrength, int nbIterations){
		char perturbationPath[255];
		sprintf(perturbationPath, "%s/DSTest_PhaseShift.dat", this->path);
		ofstream out(perturbationPath);

		double phase;
		double phasePI;
		double maxDistances = 0;
		double sumDistances = 0;
		int numberOfPoints;
		vector<point3D> vectPoint3D;
		vector<double> distances;


		// First apply the perturbation on all oscillators
		perturbate_all(stimulationStrength,0);
		// Initialize the vector
		for(int i=0; i< height;i++){
			for(int j=0; j< width;j++){
				distances.push_back(0);
				vectPoint3D.push_back(point3D(0,0,0));
			}
		}
		numberOfPoints = (int)vectPoint3D.size();

		for(int k=0; k < nbIterations; k++){
				ofstream out(perturbationPath, ios::app);
				maxDistances = -1;
				sumDistances = 0;

				for(int i=0; i< height;i++){
					for(int j=0; j< width;j++){
						distances[i*width+j] = 0;
						vectPoint3D[i*width+j].x = states[i][j][0];
						vectPoint3D[i*width+j].y = states[i][j][1];
						vectPoint3D[i*width+j].z = states[i][j][2];
						//vectPoint3D[i*width+j].phase = 0;
					}
				}

				for(int i=0; i < numberOfPoints; i++){
					phase = cyclicAttractor->phaseOf(&vectPoint3D[i],0);
					phasePI = (2*PI*phase/cyclicAttractor->periode);
					while(phasePI >= 2*PI) phasePI -= 2*PI;
					if(phase==-1)
						vectPoint3D[i].phase= -1;
					else
						vectPoint3D[i].phase= phasePI;
				}
				sort(vectPoint3D.begin(),vectPoint3D.end());
				// Calculus of the phase difference between consecutive oscillators and retaining the
				// total phase shift
				for(int i=0; i< numberOfPoints-1;i++){
					if((vectPoint3D[i].phase == -1) or (vectPoint3D[i+1].phase==-1))
						distances[i] = -1;
					else
						distances[i] = vectPoint3D[i+1].phase - vectPoint3D[i].phase;
					if(distances[i]>maxDistances) maxDistances = distances[i];
				}
				// The last element is tested alone
				if(vectPoint3D[numberOfPoints-1].phase == -1 or vectPoint3D[0].phase ==-1)
					distances[numberOfPoints-1] = -1;
				else
					distances[numberOfPoints-1] = 2*PI - (vectPoint3D[numberOfPoints-1].phase - vectPoint3D[0].phase);
				if(distances[numberOfPoints-1] > maxDistances) maxDistances = distances[numberOfPoints-1];

				// Calculus of the maximum phase shift and Writing the result
				if(maxDistances != -1){
					sumDistances = 2*PI - maxDistances;
					printf("Iteration %d , phaseShift=%f\n", k, sumDistances);
					out << k* cyclicAttractor->deltaT << " " << sumDistances <<"\n";
				}
				else {
					printf("Iteration = %d, D = -1\n",k);
					out.close();
					break;
				}
			computeMethod(1);
			out.close();
			if(k==nbIterations-1) 	out.close();
		}

		char phaseShiftGnuplotPath[255];
		sprintf(phaseShiftGnuplotPath, "%s/DSTest_PhaseShift.gp", this->path);
		ofstream outMe(phaseShiftGnuplotPath);

		outMe << "set title \"Evolution of the Maximum Phase Shift of a population of "<< cyclicAttractor->name <<" Attractor\"\n";
		outMe << "set xlabel \"t\"\n";
		outMe << "set xrange["<<0<<":"<<nbIterations<<"]\n";
		outMe << "set yrange[0:"<< 2*M_PI << "]\n";
		outMe << "set ylabel \"Maximum Phase Shift\"\n";
		outMe << "unset key\n";
		outMe << "plot \"DSTest_PhaseShift.dat\" title 'Maximum Phase Shift of a population of "<< cyclicAttractor->name << " Attractor' with pm3d\n";
		outMe << "pause -1\n";
		outMe.close();
}

void DSGrid::perturbate_randomly(){
	MTRand myRandom;
	for(int i=0; i < height;i++){
		for(int j=0; j < width;j++){
			states[i][j][0] = myRandom.rand() * 2 * cyclicAttractor->radius[0] - cyclicAttractor->radius[0];
			states[i][j][1] = myRandom.rand() * 2 * cyclicAttractor->radius[1] - cyclicAttractor->radius[1];
			states[i][j][2] = myRandom.rand() * 2 * cyclicAttractor->radius[2] - cyclicAttractor->radius[2];
/*			states[i][j][0] = 0;
			states[i][j][1] = 1;
			states[i][j][2] = 0;*/
			this->cyclicAttractor->computeMethod(states[i][j],cyclicAttractor->numberOfIterations);
		}
	}
}

/*
 * Set All the oscillators at the same phase
 */
void DSGrid::setTheSameStateForAll(double x, double y, double z){
	for(int i=0; i < height;i++){
		for(int j=0; j < width;j++){
			states[i][j][0] = x;
			states[i][j][1] = y;
			states[i][j][2] = z;
/*			states[i][j][0] = 0;
			states[i][j][1] = 1;
			states[i][j][2] = 0;*/
			this->cyclicAttractor->computeMethod(states[i][j],cyclicAttractor->numberOfIterations);
		}
	}
}


/*
 * Perturbate all the oscillators along the X line
 */
void DSGrid::perturbate_all(double stimulationStrength, double angle){
	for(int i=0; i < height; i++){
		for(int j=0; j < width; j++){
			states[i][j][0] += stimulationStrength * cos(angle);
			states[i][j][1] += stimulationStrength * sin(angle);
		}
	}
}

/*
 * All oscillators have a coupling equal to epsilon
 */
void DSGrid::buildCoupling_all(double epsilon){
	for(int i=0; i < width*height; i++){
		for(int j=0; j < width*height; j++){
			couplages[i][j][0] = epsilon;	//xi xj
			couplages[i][j][1] = 0;//epsilon;	//xi yj
			couplages[i][j][2] = 0;//epsilon;	//yi xj
			couplages[i][j][3] = 0;//epsilon;	//yi yj
		}
	}
/*
	for(int j=0; j < width*height;j++){
		couplages[210][j][0] = 0;
	}
	*/
}

/*
 * Left to Right chaine d'oscillateur sur une dimension
 */
void DSGrid::buildCoupling_Chain(double epsilon){
	buildCoupling_all(0);
	for(int i=1; i < width*height; i++){
				couplages[i][i-1][0] = epsilon;	//xi xj
				couplages[i][i-1][1] = 0;//epsilon;	//xi yj
				couplages[i][i-1][2] = 0;//epsilon;	//yi xj
				couplages[i][i-1][3] = 0;//epsilon;	//yi yj

				couplages[i-1][i][0] = 0;	//xi xj
				couplages[i-1][i][1] = 0;//epsilon;	//xi yj
				couplages[i-1][i][2] = 0;//epsilon;	//yi xj
				couplages[i-1][i][3] = 0;//epsilon;	//yi yj
	}
}

/*
 * Perceptron One layer to one output
 */
void DSGrid::buildCoupling_PerceptronOneLayerToOneOutput(double epsilon){
	buildCoupling_all(0);
	/*
	for(int i=0; i < width*height-1; i++){
				couplages[width*height-1][i][0] = epsilon;	//xi xj
				couplages[width*height-1][i][1] = 0;//epsilon;	//xi yj
				couplages[width*height-1][i][2] = 0;//epsilon;	//yi xj
				couplages[width*height-1][i][3] = 0;//epsilon;	//yi yj
	}
	*/
	couplages[width*height-1][0][0] = epsilon*0.1;	//xi xj
	couplages[width*height-1][1][0] = epsilon;	//xi xj

}


/*
 * Coupling exists only between up, up left, up right ,down, down left, down right, left, right
 * Voisinage = 8
 */
void DSGrid::buildCoupling_Moore(double epsilon){
	int upLeftID, downRightID, upRightID, downLeftID;

	buildCoupling_VonNeumann(epsilon);
	for(int i=0; i < height - 1; i++){
		for(int j=0; j < width - 1; j++){
			upLeftID = getID(i,j);
			downRightID = getID(i+1,j+1);
			couplages[upLeftID][downRightID][0] = epsilon;
			couplages[upLeftID][downRightID][1] = -epsilon;
			couplages[upLeftID][downRightID][3] = 0;//epsilon;
			couplages[downRightID][upLeftID][0] = epsilon;
			couplages[downRightID][upLeftID][1] = -epsilon;
			couplages[downRightID][upLeftID][3] = 0;//epsilon;
		}
		// For a toric coupling
		upLeftID = getID(i,width-1);
		downRightID = getID(i+1,0);
		couplages[upLeftID][downRightID][0] = epsilon;
		couplages[upLeftID][downRightID][1] = -epsilon;
		couplages[upLeftID][downRightID][3] = 0;//epsilon;
		couplages[downRightID][upLeftID][0] = epsilon;
		couplages[downRightID][upLeftID][1] = -epsilon;
		couplages[downRightID][upLeftID][3] = 0;//epsilon;
	}
	for(int j=0; j < width-1;j++){
		upLeftID = getID(height-1,j);
		downRightID = getID(0,j+1);
		couplages[upLeftID][downRightID][0] = epsilon;
		couplages[upLeftID][downRightID][1] = -epsilon;
		couplages[upLeftID][downRightID][3] = 0;//epsilon;
		couplages[downRightID][upLeftID][0] = epsilon;
		couplages[downRightID][upLeftID][1] = -epsilon;
		couplages[downRightID][upLeftID][3] = 0;//epsilon;
	}
	upLeftID = getID(height-1,width-1);
	downRightID = getID(0,0);
	couplages[upLeftID][downRightID][0] = epsilon;
	couplages[upLeftID][downRightID][1] = -epsilon;
	couplages[upLeftID][downRightID][3] = 0;//epsilon;
	couplages[downRightID][upLeftID][0] = epsilon;
	couplages[downRightID][upLeftID][1] = -epsilon;
	couplages[downRightID][upLeftID][3] = 0;//epsilon;


	// Couplage gauche droite
	for(int i=0; i < height-1; i++){
		for(int j=1; j < width; j++){
			downLeftID = getID(i,j);
			upRightID = getID (i+1,j-1);
			couplages[upRightID][downLeftID][0] = epsilon;
			couplages[upRightID][downLeftID][1] = -epsilon;
			couplages[upRightID][downLeftID][3] = 0;//epsilon;
			couplages[downLeftID][upRightID][0] = epsilon;
			couplages[downLeftID][upRightID][1] = -epsilon;
			couplages[downLeftID][upRightID][3] = 0;//epsilon;
		}
		// Toric coupling

		downLeftID = getID(i,0);
		upRightID = getID (i+1,width-1);
		couplages[upRightID][downLeftID][0] = epsilon;
		couplages[upRightID][downLeftID][1] = -epsilon;
		couplages[upRightID][downLeftID][3] = 0;//epsilon;
		couplages[downLeftID][upRightID][0] = epsilon;
		couplages[downLeftID][upRightID][1] = -epsilon;
		couplages[downLeftID][upRightID][3] = 0;//epsilon;
	}
	for(int j=1; j < width; j++){
		downLeftID = getID(height-1,j);
		upRightID = getID (0,j-1);
		couplages[upRightID][downLeftID][0] = epsilon;
		couplages[upRightID][downLeftID][1] = -epsilon;
		couplages[upRightID][downLeftID][3] = 0;//epsilon;
		couplages[downLeftID][upRightID][0] = epsilon;
		couplages[downLeftID][upRightID][1] = -epsilon;
		couplages[downLeftID][upRightID][3] = 0;//epsilon;
	}
	downLeftID = getID(height-1,0);
	upRightID = getID (0,width-1);
	couplages[upRightID][downLeftID][0] = epsilon;
	couplages[upRightID][downLeftID][1] = -epsilon;
	couplages[upRightID][downLeftID][3] = 0;//epsilon;
	couplages[downLeftID][upRightID][0] = epsilon;
	couplages[downLeftID][upRightID][1] = -epsilon;
	couplages[downLeftID][upRightID][3] = 0;//epsilon;

	/*
	for(int j=0; j < width*height;j++){
			couplages[210][j][0] = 0;
		}*/
}

/*
 * Coupling exists only between up,down,left and right cells
 * Voisinage = 4
 */
void DSGrid::buildCoupling_VonNeumann(double epsilon){
	// Couplage haut bas droite gauche
		int upID, downID, leftID, rightID;

		buildCoupling_all(0);

		for(int i=0; i < height; i++){
			for(int j=0; j < width - 1; j++){
				upID = getID(i,j);
				downID = getID(i,j+1);
				couplages[upID][downID][0] = epsilon;
				couplages[upID][downID][1] = -epsilon;
				couplages[upID][downID][3] = 0;//epsilon;
				couplages[downID][upID][0] = epsilon;
				couplages[downID][upID][1] = -epsilon;
				couplages[downID][upID][3] = 0;//epsilon;
			}

			// For a toric coupling : The final line with the first one
			upID = getID(i, width-1);
			downID = getID(i, 0);
			couplages[upID][downID][0] = epsilon;
			couplages[upID][downID][1] = -epsilon;
			couplages[upID][downID][3] = 0;//epsilon;
			couplages[downID][upID][0] = epsilon;
			couplages[downID][upID][1] = -epsilon;
			couplages[downID][upID][3] = 0;//epsilon;
		}
		// Couplage gauche droite
		for(int i=0; i < height-1; i++){
			for(int j=0; j < width; j++){
				leftID = getID(i,j);
				rightID = getID (i+1,j);
				couplages[leftID][rightID][0] = epsilon;
				couplages[leftID][rightID][1] = -epsilon;
				couplages[leftID][rightID][3] = 0;//epsilon;
				couplages[rightID][leftID][0] = epsilon;
				couplages[rightID][leftID][1] = -epsilon;
				couplages[rightID][leftID][3] = 0;//epsilon;
			}
		}
		// For a toring coupling: the final column is coupled with the first one
		for(int j=0; j < width; j++){
			leftID = getID(height-1,j);
			rightID = getID (0,j);
			couplages[leftID][rightID][0] = epsilon;
			couplages[leftID][rightID][1] = -epsilon;
			couplages[leftID][rightID][3] = 0;//epsilon;
			couplages[rightID][leftID][0] = epsilon;
			couplages[rightID][leftID][1] = -epsilon;
			couplages[rightID][leftID][3] = 0;//epsilon;
		}

		/*
		for(int j=0; j < width*height;j++){
				couplages[210][j][0] = 0;
			}*/
}

void DSGrid::buildSeuils(double initialThreshold){
	for(int i = 0; i < height; i++){
		for(int j = 0;j < width; j++){
			seuils[i][j] = initialThreshold;
		}
	}
}

/*
 * Return the label of the oscillator from the coordinates on the matrices
 */
int DSGrid::getID(int i, int j){
	return i*width+j;
}

/*
 * Compute the new state of the grid after nbCompute time
 * using a simple euler integration
 */
void DSGrid::compute(int nbCompute){
	int tmpCurrentID;
	double tmpSumVoisinage = 0;
	double tc; // terme de couplage
	for(int k=0; k < nbCompute; k++){
		for(int i=0; i < height;i++){
			for(int j=0; j < width;j++){
				tmpStates[i][j][0] = states[i][j][0];
				tmpStates[i][j][1] = states[i][j][1];
				tmpStates[i][j][2] = states[i][j][2];
				cyclicAttractor->compute(tmpStates[i][j],1);
				// Add the coupling to the new states following the rung kutta 4th order methods
				tmpSumVoisinage = 0;
				for(int l =0; l < height*width;l++) {
					tmpCurrentID = getID(i,j);
					 // Couplage simple
					tmpStates[i][j][0] += cyclicAttractor->deltaT * (couplages[tmpCurrentID][l][0] * couplingFunction(states[l/width][l%width][0])+couplages[tmpCurrentID][l][1] * couplingFunction(states[l/width][l%width][1]));
					tmpStates[i][j][1] += cyclicAttractor->deltaT * (couplages[tmpCurrentID][l][2] * couplingFunction(states[l/width][l%width][0])+couplages[tmpCurrentID][l][3] * couplingFunction(states[l/width][l%width][1]));
					// Couplage fireflies  IIIIIIIIIIIIII
					/*
					if(couplages[tmpCurrentID][l][0] > initialCoupling)
						tmpSumVoisinage += couplages[tmpCurrentID][l][0] * states[l/width][l%width][0];
						*/
					// Couplage par écart sur les variables correspondantes
					//tmpStates[i][j][0] += cyclicAttractor->deltaT * couplingFunction(couplages[tmpCurrentID][l][0] * (states[i][j][0] - states[l/width][l%width][0]) + couplages[tmpCurrentID][l][1] * (states[i][j][0] - states[l/width][l%width][1]));
					//tmpStates[i][j][1] += cyclicAttractor->deltaT * couplingFunction(couplages[tmpCurrentID][l][2] * (states[i][j][1] - states[l/width][l%width][0]) + couplages[tmpCurrentID][l][3] * (states[i][j][1] - states[l/width][l%width][1]));
					// Couplages par pulsation
					//tmpStates[i][j][0] += cyclicAttractor->deltaT * (couplages[tmpCurrentID][l][0] * (states[i][j][0] - states[l/width][l%width][0]) * couplingFunction(states[i][j][0]) * couplingFunction(states[l/width][l%width][0]) + couplages[tmpCurrentID][l][1] * (states[i][j][0] - states[l/width][l%width][1]) * couplingFunction(states[i][j][0])* couplingFunction(states[l/width][l%width][1]));
					//tmpStates[i][j][0] += cyclicAttractor->deltaT * (couplages[tmpCurrentID][l][2] * (states[i][j][1] - states[l/width][l%width][0]) * couplingFunction(states[i][j][1]) * couplingFunction(states[l/width][l%width][0]) + couplages[tmpCurrentID][l][3] * (states[i][j][1] - states[l/width][l%width][1]) * couplingFunction(states[i][j][1])* couplingFunction(states[l/width][l%width][1]));
				}
				/* Simple learning IIIIIIIIIIIIIIII
				 tc = cyclicAttractor->deltaT * alpha * (Xmax-states[i][j][0]) * (1 - couplingFunction(states[i][j][0]-Xmax)) * couplingFunction(tmpSumVoisinage - seuils[i][j]);
				 if(tc>PXLearningZone){
	 					 tc = PXLearningZone;
				 }
				 if(tc > 0){
					 tmpStates[i][j][0] += tc;
					// printf("TC = %f\n",tc);
				 }
				 */
			}
		}

		for(int i=0; i < height;i++){
			for(int j=0; j < width;j++){
				states[i][j][0] = tmpStates[i][j][0];
				states[i][j][1] = tmpStates[i][j][1];
				states[i][j][2] = tmpStates[i][j][2];
			}
		}
	}
}

/*
 * Compute the next state of the grid using a runge kutta 4th order integration
 */
inline void DSGrid::compute_RK4(int nbCompute){

	double k1[MAX_INT][MAX_INT][3], k2[MAX_INT][MAX_INT][3], k3[MAX_INT][MAX_INT][3], k4[MAX_INT][MAX_INT][3];
	double tmpCoor[MAX_INT][MAX_INT][3], originalState[MAX_INT][MAX_INT][3];
	int nbDimension = cyclicAttractor->nbDimension;

	for(int k=0;k < nbCompute;k++){
		for(int i=0; i < nbDimension; i++){
			for(int l=0; l < height; l++){
				for(int m=0; m < width; m++){
					tmpCoor[l][m][i] = states[l][m][i];
					originalState[l][m][i] = states[l][m][i];
				}
			}
		}

		compute(1);
		for(int i=0; i< nbDimension;i++){
			for(int l=0; l < height; l++){
					for(int m=0; m < width; m++){
						tmpCoor[l][m][i] = states[l][m][i];
					}
			}
		}

		for(int i=0; i < nbDimension; i++) {
			for(int l=0; l < height; l++){
					for(int m=0; m < width; m++){
						k1[l][m][i] = tmpCoor[l][m][i] - originalState[l][m][i];
						tmpCoor[l][m][i]=originalState[l][m][i] + (k1[l][m][i]/2);
						states[l][m][i] = tmpCoor[l][m][i];
					}
			}
		}

		compute(1);
		for(int i=0; i< nbDimension;i++){
			for(int l=0; l < height; l++){
				for(int m=0; m < width; m++){
					tmpCoor[l][m][i] = states[l][m][i];
				}
			}
		}

		for(int i=0; i < nbDimension; i++) {
			for(int l=0; l < height; l++){
				for(int m=0; m < width; m++){
					k2[l][m][i] = (tmpCoor[l][m][i] - (originalState[l][m][i] + (k1[l][m][i]/2)));
					tmpCoor[l][m][i]=originalState[l][m][i] + (k2[l][m][i]/2);
					states[l][m][i] = tmpCoor[l][m][i];
				}
			}
		}

		compute(1);
		for(int i=0; i< nbDimension;i++){
			for(int l=0; l < height; l++){
				for(int m=0; m < width; m++){
					tmpCoor[l][m][i] = states[l][m][i];
				}
			}
		}

		for(int i=0; i < nbDimension; i++) {
			for(int l=0; l < height; l++){
				for(int m=0; m < width; m++){
					k3[l][m][i] = (tmpCoor[l][m][i] - (originalState[l][m][i] + (k2[l][m][i]/2)));
					tmpCoor[l][m][i]=originalState[l][m][i] + k3[l][m][i];
					states[l][m][i] = tmpCoor[l][m][i];
				}
			}
		}

		compute(1);
		for(int i=0; i< nbDimension;i++){
			for(int l=0; l < height; l++){
				for(int m=0; m < width; m++){
					tmpCoor[l][m][i] = states[l][m][i];
				}
			}
		}

		for(int i=0; i < nbDimension; i++) {
			for(int l=0; l < height; l++){
				for(int m=0; m < width; m++){
					k4[l][m][i] = tmpCoor[l][m][i] - (originalState[l][m][i] + k3[l][m][i]);
					originalState[l][m][i] += ((k1[l][m][i]/6) + (k2[l][m][i]/3) + (k3[l][m][i]/3) + (k4[l][m][i]/6));
					states[l][m][i] = originalState[l][m][i];
				}
			}
		}
	}

}

/*
 * The global signal is a simple addition of the grid states
 */
void DSGrid::refreshGlobalSignal(){
	globalSignal[0] = 0;
	globalSignal[1] = 0;
	globalSignal[2] = 0;

	for(int i=0; i < height;i++){
		for(int j=0; j < width;j++){
			globalSignal[0] += states[i][j][0];
			globalSignal[1] += states[i][j][1];
			globalSignal[2] += states[i][j][2];
		}
	}
}

void DSGrid::saveStates(char outPutFile[255]){
	//char outPutFile[255];
	char * image = (char*) malloc(width*height);
	Bitmap * myBitMap = new Bitmap();
	int colorIndex;
	//sprintf(outPutFile, "./%s/DSTest_States.bmp", path);
	//ofstream out(outPutFile);
	for(int i=0; i < height;i++){
		for(int j=0; j < width;j++){
			/*colorIndex = ((double)((states[i][j][0] + cyclicAttractor->radius[0])/(2 * (cyclicAttractor->radius[0]))) * 256);
			if(colorIndex > 255 or colorIndex < 0) colorIndex = 0;
			image[i*width+j] = (char) colorIndex;*/
			colorIndex = (double)(((states[i][j][0] + 1)/3) * 255);
			//printf("colorIndex =%d\n", colorIndex);
			if(colorIndex > 255) colorIndex = 255;
			if(colorIndex < 0) colorIndex = 0;
			image[(height-1-i)*width+j] = (unsigned char) (int) colorIndex;
			//printf("relaxed %f\n",states[i][j][0]);
			//out << i << "\t" << j << "\t" << states[i][j][0] << "\n";
			/*states[i][j][0] = tmpStates[i][j][0];
			states[i][j][1] = tmpStates[i][j][1];
			states[i][j][2] = tmpStates[i][j][2];*/
		}
		//out << "\n";
	}
	//out.close();
	myBitMap->save8BMP(outPutFile, image, width, height);
	free(image);
}

void DSGrid::saveFilterState(char outPutFile[255]){
	//char outPutFile[255];
	char * image = (char*) malloc(width*height);
	Bitmap * myBitMap = new Bitmap();
	//sprintf(outPutFile, "./%s/DSTest_States.bmp", path);
	//ofstream out(outPutFile);
	for(int i=0; i < height;i++){
		for(int j=0; j < width;j++){
			if(fabs(states[i][j][0] - 1) <= 0.01)
				image[i*width+j] = (char) (255);
			else
				image[i*width+j] = (char) (0);
			//out << i << "\t" << j << "\t" << states[i][j][0] << "\n";
			/*states[i][j][0] = tmpStates[i][j][0];
			states[i][j][1] = tmpStates[i][j][1];
			states[i][j][2] = tmpStates[i][j][2];*/
		}
		//out << "\n";
	}
	//out.close();
	myBitMap->save8BMP(outPutFile, image, width, height);
	free(image);
}

void DSGrid::saveLocalSynchrony(char outPutFile[255]){
	//char outPutFile[255];
	char * image = (char*) malloc(width*height);
	Bitmap * myBitMap = new Bitmap();
	//sprintf(outPutFile, "./%s/DSTest_States.bmp", path);
	//ofstream out(outPutFile);
	double sumSync=0;
	int nbOfElements;
	int tmpCurrentID;
	for(int i=0; i < height;i++){
		for(int j=0; j < width;j++){
			tmpCurrentID = getID(i,j);
			sumSync = states[i][j][0];
			nbOfElements = 1;
			for(int k=0;k<width*height;k++){
				if(couplages[tmpCurrentID][k][0] > 0){
					sumSync += states[k/height][k%width][0];
					nbOfElements++;
				}
			}
			if(sumSync > (double)nbOfElements*filterCoeff)
				image[i*width+j] = (char) (255);
			else
				image[i*width+j] = (char) (0);
			//out << i << "\t" << j << "\t" << states[i][j][0] << "\n";
			/*states[i][j][0] = tmpStates[i][j][0];
			states[i][j][1] = tmpStates[i][j][1];
			states[i][j][2] = tmpStates[i][j][2];*/
		}
		//out << "\n";
	}
	//out.close();
	myBitMap->save8BMP(outPutFile, image, width, height);
	free(image);
}

inline void DSGrid::saveMethod(char outPutFile[255]){
	//savePhases(outPutFile);
	saveStates(outPutFile);
	//saveFilterState(outPutFile);
	//saveLocalSynchrony(outPutFile);
}

void DSGrid::savePhases(char outPutFile[255]){
	char * image = (char*) malloc(width*height);
	Bitmap * myBitMap = new Bitmap();
	point3D p3d(0,0,0);
	double phase;
	int colorIndex = 0;
	for(int i=0; i < height;i++){
		for(int j=0; j < width;j++){
			p3d.x = states[i][j][0];
			p3d.y = states[i][j][1];
			p3d.z = states[i][j][2];
			phase = cyclicAttractor->phaseOf(&p3d,0);
			phase = (2*PI*phase/cyclicAttractor->periode);
			while(phase >= 2*PI) phase -= 2*PI;
			colorIndex = ((phase/(2*PI)) * 255);
			if(colorIndex < 0 or colorIndex > 255) colorIndex = 0;
			image[(height-1-i)*width+j] = (char) colorIndex;
		}
	}
	myBitMap->save8BMP(outPutFile, image, width, height);
	free(image);

}

/*
 * Save The global signal into appropriate files
 */
void DSGrid::saveGlobalSignal(int iteration){
		char signalOutputFileX[255];
		sprintf(signalOutputFileX, "./%s/DSTest_Global_Signal_X.dat", path);

		char signalOutputFileY[255];
		sprintf(signalOutputFileY, "./%s/DSTest_Global_Signal_Y.dat", path);

		ofstream outX;
		ofstream outY;
		// First time opening?!
		if(iteration == 0){
			outX.open(signalOutputFileX);
			outY.open(signalOutputFileY);
		}
		else{
			outX.open(signalOutputFileX, ios::app);
			outY.open(signalOutputFileY, ios::app);
		}
		refreshGlobalSignal();
		outX << iteration*cyclicAttractor->deltaT << " " << globalSignal[0] << "\n";
		outY << iteration*cyclicAttractor->deltaT << " " << globalSignal[1] << "\n";

		outX.close();
		outY.close();
}

/*
 *	Save The signals of the grid
 */
void DSGrid::saveSignals(int iteration){
	char signalOutputFileX[255];
	for(int i=0; i < height*width; i++){
		sprintf(signalOutputFileX, "./%s/Signals/DSTest_Signal_X_%d.dat", path,i);

		ofstream outX;
		// First time opening?!
		if(iteration == 0){
			outX.open(signalOutputFileX);
		}
		else{
			outX.open(signalOutputFileX, ios::app);
		}
		outX << iteration*cyclicAttractor->deltaT << " " << states[i/width][i%width][0] << "\n";

		outX.close();
	}
}

/*
 *	Save The phaseLag correlation (phaseLagBetweenInputs, phaseLagbetween_the_output_and_one_input
 */
void DSGrid::savePhaseLagResponse(){
	char PhaseLagResponseFile[255];
	double phi1, phi2;
	sprintf(PhaseLagResponseFile, "./%s/Signals/DSTest_PhaseLagResponse.dat", path);
	phi1 = getPhaseLag_Radian(states[0][0],states[0][1]);		// Reference
	phi2 = getPhaseLag_Radian(states[0][0],states[0][2]);		// PhaseLag to evaluate (old value on 2nd index is 1 for chain coupling)
	ofstream outPhaseLagResponse;
	// First time opening?!
	outPhaseLagResponse.open(PhaseLagResponseFile, ios::app);
	outPhaseLagResponse << phi1 << " " << phi2 << "\n";
	outPhaseLagResponse.close();
}


/*
 *	Save The phaseLag between oscillator 0 and 1
 */
void DSGrid::savePhaseLag(int iteration, int simNumber){
	char signalOutputFileX[255];
	double phaseLag = 0;
	double phi1, phi2;
	sprintf(signalOutputFileX, "./%s/Signals/DSTest_PhaseLag%d.dat", path,simNumber);
	phi1 = cyclicAttractor->phaseOf(states[0][0]);		// Reference
	phi2 = cyclicAttractor->phaseOf(states[0][1]);		// PhaseLag to evaluate (old value on 2nd index is 1 for chain coupling)
	phaseLag =  fabs(phi2-phi1);
	if(fabs(cyclicAttractor->periode - phi2 + phi1) < phaseLag)
		phaseLag = fabs(cyclicAttractor->periode - phi2 + phi1);
	if(fabs(cyclicAttractor->periode - phi1 + phi2) < phaseLag)
			phaseLag = fabs(cyclicAttractor->periode - phi1 + phi2);
	ofstream outX;
	// First time opening?!
	if(iteration == 0){
		outX.open(signalOutputFileX);
	}
	else{
		outX.open(signalOutputFileX, ios::app);
	}
	outX << iteration*cyclicAttractor->deltaT << " " << phaseLag << "\n";

	outX.close();
}
/*
 * Return the phase lag of the slave from the master
 * The returned value is between [0 , 2*PI[
 */
double DSGrid::getPhaseLag_Radian(double * stateMaster, double * stateSlave){
	double phi1 = cyclicAttractor->phaseOf(stateMaster);		// Master's absolute phase
	double phi2 = cyclicAttractor->phaseOf(stateSlave);		// Slave's absolute phase
	double phaseLag =  phi2-phi1;
	if(phaseLag < 0) phaseLag += cyclicAttractor->periode;
	return phaseLag * (2*M_PI/cyclicAttractor->periode);
}

/*
 * Specify here the type of the noise to produce on the grid at specific moment
 */
inline void DSGrid::applyNoise(int iteration){
	/*if ((iteration>(1000/3)-30) and (iteration<(1000/3)+30) and (iteration%4==0)) perturbate_burnH(5,0);
	if ((iteration>(2*1000/3)-10) and (iteration<(2*1000/3)+10) and (iteration%4==0)) perturbate_all(1,0);
	if ((iteration>(2100-10)) and (iteration<(2100+10)) and (iteration%4==0)) perturbate_all(1,0);*/
	if(iteration == 2) {
		char imagePath[255];
		sprintf(imagePath, "/home/hbenamor/workspace/DSAnalyser/Debug/imgt.bmp");
		perturbate_burnImage(imagePath, 10);
	}
	//else if(iteration > 2 and iteration <20){
		//performLearning();
	//}
	else if(iteration == 20) {
		perturbate_randomly();
	}
	else if(iteration == 30) {
		char imagePath[255];
		sprintf(imagePath, "/home/hbenamor/workspace/DSAnalyser/Debug/imgt.bmp");
		perturbate_evokeImage(imagePath, 10);
	}
	/*
	if(iteration == 0){
		perturbate_all(5,0);
	}
	else if((iteration - (int)(2*cyclicAttractor->periode/cyclicAttractor->deltaT))==1){
		perturbate_burn1(3,0);
		//perturbate_burnD(5,0);
		//perturbate_burnH(5,PI/4);
		//perturbate_burn2(5,PI/2);
		//perturbate_burnH(4,27*PI/230);
		//perturbate_burnC(4,109*PI/230);
		//perturbate_burn8(7, 8*PI/10);
	}
	*/
}

/*
 * The unit (i,j) adjust its weights corresponding to the neighbourhood relationship
 */
/*
void DSGrid::learningFromNeighbourhood(int i, int j){
	int tmpCurrentID = getID(i,j);
	int tmpSeuils = -1;
	for(int l =0; l < height*width;l++) {
		if((couplages[tmpCurrentID][l][0] == initialCoupling) and (couplingFunction(states[l/width][l%width][0]-PXLearningZone)>0)){
				couplages[tmpCurrentID][l][0] = couplingFunction(states[l/width][l%width][0]-PXLearningZone);
				tmpSeuils++;
		}
		else if(couplages[tmpCurrentID][l][0]>initialCoupling) tmpSeuils++;
	}
	// Le seuils prend toujours le minimum
	if((tmpSeuils < seuils[i][j]) and (tmpSeuils > -1)){
		seuils[i][j] = tmpSeuils;
	}
}
*/

/*
 * Learning From Neighbourhood
 * Dwij = TauLearning * s(xi)*s(xj)*(maxcoupling/numOfNeighboor)
 * Be sure that sum(wij,j) < maxcoupling
 */
void DSGrid::learningFromNeighbourhood(int i, int j){
	int tmpCurrentID = getID(i,j);

	double maxCoupling = 0.1;
	int nbNeighboors = 8;
	double dw=0;
	double weakestValue = 0.00000001;
	double Taulearning = 2;
	for(int l =0; l < height*width;l++) {
		for(int k=0; k<4;k++){
			if(couplages[tmpCurrentID][l][k] != 0){
				if(k==0) dw = Taulearning * couplingFunction(states[i][j][0]) * couplingFunction(states[l/width][l%width][0])*maxCoupling/nbNeighboors;
				else if(k==1) dw = Taulearning * couplingFunction(states[i][j][0]) * couplingFunction(states[l/width][l%width][1])*maxCoupling/nbNeighboors;
				else if(k==2) dw = Taulearning * couplingFunction(states[i][j][1]) * couplingFunction(states[l/width][l%width][0])*maxCoupling/nbNeighboors;
				else if(k==4) dw = Taulearning * couplingFunction(states[i][j][1]) * couplingFunction(states[l/width][l%width][1])*maxCoupling/nbNeighboors;
				if(couplages[tmpCurrentID][l][k] > 0){	// Cas de l'excitation
					if(couplages[tmpCurrentID][l][k]+dw<=0) couplages[tmpCurrentID][l][k] = weakestValue;
					else if(couplages[tmpCurrentID][l][k]+dw>maxCoupling/nbNeighboors) couplages[tmpCurrentID][l][k] = maxCoupling/nbNeighboors;
					else couplages[tmpCurrentID][l][k]+=dw;
				}
				else if(couplages[tmpCurrentID][l][k] < 0) {	// Cas de l'inhibition
					if(couplages[tmpCurrentID][l][k]+dw>=0) couplages[tmpCurrentID][l][k] = -weakestValue;
					else if(couplages[tmpCurrentID][l][k]+dw<-(maxCoupling/nbNeighboors)) couplages[tmpCurrentID][l][k] = -maxCoupling/nbNeighboors;
					else couplages[tmpCurrentID][l][k]+=dw;
				}
				printf("dw = %f\n",dw);
			}
		}
	}
}


void DSGrid::computeAnimation(int gridNbCompute, int steps){
	char outPutFile[255];

	int j=0;
	for(int i=0; i < gridNbCompute;i++){
		applyNoise(i);
		/*if(i==10){
			states[7][14][0] += 10;
			states[14][7][0] -= 10;
		}*/

		//saveGlobalSignal(i);
		if(i%steps==0){
			sprintf(outPutFile, "./%s/DSTest_States_%d.bmp", path, j);
			saveMethod(outPutFile);
			j++;
		}
		computeMethod(1);
	}
	/*
	j=0;

	for(int i=0; i < gridNbCompute;i++){
			//applyNoise(i);
			saveGlobalSignal(i);
			if(i%steps==0){
				sprintf(outPutFile, "./%s/DSTest_States_%d.bmp", path, j);
				saveMethod(outPutFile);
				j++;
			}
			computeMethod(1);
	}*/
	if(gridNbCompute%steps==0){
		sprintf(outPutFile, "./%s/DSTest_States_%d.bmp", path, j+1);
		saveMethod(outPutFile);
	}
}


/*
 * Simple Grid Computation
 */
void DSGrid::computeGrid(int gridNbCompute, int steps){
	//for(int j=0; j < 10; j++){ // Use with savePhaseLag
		//perturbate_randomly();
		for(int i=0; i < gridNbCompute;i++){
			//saveSignals(i);
			//savePhaseLag(i, j);				// Pour avoir l'évolution du phase Lag
			//saveFinalPhaseLag(i);
			//saveGlobalSignal(i);
			computeMethod(1);
		}
		savePhaseLagResponse();
	//}
}


/*
 * Define the coupling function:
 * e.g. tanh for the wilson cowan oscillator
 * or identity, anything else... more ideas?!
 */
inline double DSGrid::couplingFunction(double input){
	/*if(input>0) return 1;
	else return 0;*/
	//return input;
	/*if(input >=0.8) return 1;
	else return 0;*/
	//double seuil = 1.1*0.5;
	return(tanh(1.1*input));
}

inline void DSGrid::computeMethod(int nbCompute){
	compute_RK4(nbCompute);
	//compute(nbCompute);
}

DSGrid::DSGrid() {

}

DSGrid::~DSGrid() {
	// TODO Auto-generated destructor stub
}
