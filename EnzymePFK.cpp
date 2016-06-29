/*
 * EnzymePFK.cpp
 *
 *  Created on: 16-Jan-2009
 *      Author: hbenamor
 */

#include "EnzymePFK.h"
#include <vector>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>
using namespace std;

EnzymePFK::EnzymePFK() {
	// TODO Auto-generated constructor stub
	int status;
	status = mkdir("./EnzymePFK", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if(status!=0) printf("Directory EnzymePFK already exists\n");
	status = mkdir("./EnzymePFK/Isochrons", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if(status!=0) printf("Directory EnzymePFK/Isochrons already exists\n");

}

EnzymePFK::~EnzymePFK() {
	// TODO Auto-generated destructor stub
}

void EnzymePFK::tracer(double deltaT, int numberOfIterations){

	ofstream out("./EnzymePFK/attractor.dat");

	double x = -0.3561;//0.17811457;//= 0.3;
	double y = -30.9889 ;//0.67;
	double z = 0;
	this->deltaT = deltaT;
	this->numberOfIterations = numberOfIterations;
	//double deltaT = 0.001;
	L0 = 3000;
	n = 3;
	C = 0.02;
	D = 1;
	N = 0.01;
	A = 0.1;
	P=3; //(c'est le paramètre sur lequel on peut jouer, avec P disons entre 2 et 5 -> 'durci' le cycle);
	R = 1000.0;//pow(10, P);


	double xTemp = x;
	double yTemp = y;

	// Boucle d'initialisation afin de rapprocher le point initial de l'attracteur.
	/*for(int i=0; i < 999; i++){
			L = x * (pow((1+x),(n-1)) * pow((1+D * y),n) + L0 * C * pow((1+C * x),(n-1))) / (pow((1+x),n) * pow((1+D * y),n) + L0 * pow((1+C * x),n));
			xTemp += (deltaT * (A-L));
			yTemp += (deltaT * (R *(L - N * y)));
			x = xTemp;
			y = yTemp;

	}*/

	attractor.push_back(new point3D(x,y,0));

	// Tracer le portrait de phase du système dynamique
	for(int i=0; i < numberOfIterations; i++){
		L = x * (pow((1+x),(n-1)) * pow((1+D * y),n) + L0 * C * pow((1+C * x),(n-1))) / (pow((1+x),n) * pow((1+D * y),n) + L0 * pow((1+C * x),n));
		xTemp += (deltaT * (A-L));
		yTemp += (deltaT * (R *(L - N * y)));
		x = xTemp;
		y = yTemp;
		//out << x << " " << y << " " << z << "\n";
		attractor.push_back(new point3D(x,y,0));
	}

	// Trouve les points appartenant à l'attracteur
	bool found = false;
	bool inSegment = false;
	double maxDistance = distancePointSegment(attractor[numberOfIterations-2],attractor[numberOfIterations-3], attractor[numberOfIterations-1], &inSegment);
	double distance;
	int i=numberOfIterations - 4;
	attractor[i+1]->inTheAttractor = true;
	attractor[i+2]->inTheAttractor = true;
	attractor[i+3]->inTheAttractor = true;

	while(i >= 0 and not found){
		attractor[i]->inTheAttractor = true;
		distance = distancePointSegment(attractor[i],attractor[numberOfIterations-3], attractor[numberOfIterations-1], &inSegment);
		if(distance <= maxDistance)
			found=true;
		else i--;
	}

	/*
	 * En partant du dernier point ajouter chercher le premier point
	 * qui pourrait être confondu avec.
	 *
	 */


	// Imprime que l'attracteur
	printf("i = %d\n", i);

	xradius = abs(attractor[i+1]->x);
	yradius = abs(attractor[i+1]->y);
	zradius = 0;

	for(int j=i; j< numberOfIterations; j++){//j=i+1
		if(abs(attractor[j]->x) > xradius) xradius = abs(attractor[j]->x);
		if(abs(attractor[j]->y) > yradius) yradius = abs(attractor[j]->y);
		if(abs(attractor[j]->z) > zradius) zradius = abs(attractor[j]->z);
		out << attractor[j]->x << " " << attractor[j]->y << "\n";
	}


	// ça pourrait etre utile de délimiter l'attracteur pas une région dans l'espace
	xradius *= 1.5;
	yradius *= 1.5;

	periode  = (numberOfIterations-i-1)*deltaT;
	nbPointsPerPeriode = numberOfIterations-i-1;
	printf("attracteur de période = %f\n", periode);
	//out << attractor[i+1]->x << " " << attractor[i+1]->y << " " << attractor[i+1]->z << "\n";

	//for(int k=0; k<=i;k++) attractor.pop_back();
	/*
	 *
	 * Ce bout de code permet de garder qu'une période de l'attracteur
	 * Les points obtenus correspondent à une seule phase
	 */
	/*
	for(int l=i+1; l<numberOfIterations;l++){
		delete attractor[l-i-1];
		attractor[l-i-1] = attractor[l];
	}

	for(int l=numberOfIterations-1-i; l<numberOfIterations;l++){
		delete attractor[l];
	}
	for(int l=numberOfIterations-1-i; l<numberOfIterations;l++){
			attractor.pop_back();
	}
	*/
	/*
	 *
	 * Gardez plusieurs points par soucis de précision
	 * Penser à ne garder qu'un multiple de point par période dans le vecteur
	 */
	for(int l=0; l<i; l++){
		delete attractor[l];
	}
	attractor.erase(attractor.begin(), attractor.begin()+i-1);
	int nbPointsToRemove = ((int)attractor.size()-1) - nbPointsPerPeriode * (int)((attractor.size()-1) / nbPointsPerPeriode);
	printf("NbAttractor=%d, NBPeriod=%d, Nb point to remove=%d\n", attractor.size(), nbPointsPerPeriode, nbPointsToRemove);
	for(int l=attractor.size()-nbPointsToRemove; l < attractor.size(); l++){
		delete attractor[l];
	}
	for(int l=0; l < nbPointsToRemove; l++){
			attractor.pop_back();
	}

	printf("firsPoint-lastpoint=%f", attractor[0]->x -attractor[attractor.size()-1]->x);
	printf("firsPoint-lastpoint=%f", attractor[0]->y -attractor[attractor.size()-1]->y);
	printf("firsPoint-lastpoint=%f", attractor[0]->z -attractor[attractor.size()-1]->z);
	ofstream outMe("simulation.gp");
	outMe << "set title \"EnzymePFK strange Attractor" << " \"\n";
	outMe << "set xlabel \"X\"\n";
	outMe << "set ylabel \"Y %\"\n";
	outMe << "set zlabel \"Z %\"\n";
	//outMe << "set yrange [0:101]\n";
	//outMe << "set xtic auto\n";
	//outMe << "set ytic auto\n";
//	outMe << "plot \"simulation.dat\" title 'Lorenz Attractor' with lines\n";
	outMe << "splot \"simulation.dat\" using 1:2:3 title 'Van Der Pol Attractor' with lines\n";
}

void EnzymePFK::tracerIsochrons(int nbIsochrons, int nbPointsPerIsochrons){
/*
 * Cet algo est inspiré de la méthode de monte-carlo: on fait des tirs dans l'espace
 * des phases et pour chaque tirs (s'il est admissible: contenu dans le bassin d'attraction)
 * on détérmine sa phase (là ou il rejoint l'attracteur) si on a suffisament de points
 * pour une phase donnée on les classes pour tracer la fibre isochrone correspondante.
 * Cet algo est barbare mais est parallélisable BINGO!! ALler il faut l'ecrire assez
 * de description, envisageons alors d'extraire un ensemble de points de l'attracteur
 * et donnons à chacun d'entre eux un score qui correspond au nombre de point trouvés
 * qui ont cette phase (A une formulation prés! :p ).
 *
 */

	double x, y, z;
	double phase;
	int totalPoints = nbIsochrons * nbPointsPerIsochrons;
	int k=0;
	double tmpPhase;
	bool added = false;
	int deltaIsochrons;
	int indexPointInTheAttractor;
	int indexOfIsochronsInTheAttractor;
	int indexOfIsochrons;
	double isochronsFitness = deltaT;
	point3D p3d(0,0,0);
	bool saved[nbIsochrons];
	ofstream out("./EnzymePFK/Isochrons/isochrons.dat");
	char fpaths[255];
	char fKpaths[255];
	ofstream iso;
	int maxInMemory = 10;

	Random::Uniform<double>(0.0,xradius);
	if(nbIsochrons < (1+((int)(periode/deltaT))))  deltaIsochrons = nbPointsPerPeriode / nbIsochrons;
	else deltaIsochrons = 1;

	 if(nbPointsPerPeriode > 2){
		 for(int i=0; i < nbIsochrons; i++){
				sprintf(fpaths, "./EnzymePFK/Isochrons/iso%d.dat",i);
				saved[i] = false;
				indexOfIsochronsInTheAttractor = i * nbPointsPerPeriode / nbIsochrons;
				//isochrons.push_back(new Isochrons(i * periode/(double)nbIsochrons, fpaths, maxInMemory));
				isochrons.push_back(new Isochrons(indexOfIsochronsInTheAttractor*deltaT, fpaths, maxInMemory, i, indexOfIsochronsInTheAttractor, attractor[indexOfIsochronsInTheAttractor]));
				x = attractor[indexOfIsochronsInTheAttractor]->x;
				y = attractor[indexOfIsochronsInTheAttractor]->y;
				z = 0;
				printf("iso %d index=%d,  x=%f, y=%f, phase=%f\n",i, indexOfIsochronsInTheAttractor, x, y, indexOfIsochronsInTheAttractor*deltaT);
				isochrons[(int)isochrons.size() - 1]->points.push_back(new point3D(x,y,z));
				//printf("indexOfIsochronsInTheAttractor = %d\n", indexOfIsochronsInTheAttractor);
		 }

		 while (k<totalPoints){
			 x = Random::Uniform(-xradius, xradius);
			 y = Random::Uniform(-yradius, yradius);
			 z = 0;
			 p3d.x = x;
			 p3d.y = y;
			 p3d.z = z;
			 tmpPhase = phaseOf(&p3d, nbPointsPerIsochrons);
			 //printf("tmpPhase = %f, entier = %f\n", tmpPhase, periode/tmpPhase);
			 indexPointInTheAttractor = tmpPhase / deltaT;
			 indexOfIsochrons = indexPointInTheAttractor / deltaIsochrons;
			 if(indexOfIsochrons >= isochrons.size()) indexOfIsochrons = isochrons.size() - 1;
			 if(fabs((double)indexPointInTheAttractor*deltaT - indexOfIsochrons*(periode/(double)nbIsochrons)) < isochronsFitness){
				 if(isochrons[indexOfIsochrons]->nbPoints < nbPointsPerIsochrons){
					 added = isochrons[indexOfIsochrons]->addPoint(&p3d);
					 //isochrons[indexOfIsochrons]->nbPoints++;
					 if(added) k++;
					 for(int l=0; l < p3d.numberOfIsoPoints; l++){
						 added = isochrons[indexOfIsochrons]->addPoint(p3d.isoPoints[l]);
						 if(added) k++;
					 }
					 printf("k = %d phase=%f\n", k, tmpPhase);
					 //delete p3d;
				 }
				 else if(not saved[indexOfIsochrons]){
				 	 isochrons[indexOfIsochrons]->saveAll();
					 saved[indexOfIsochrons] = true;
					 //delete p3d;
				 }
				 //else delete p3d;
			 }
			 //else delete p3d;
			 k=0;
			 for(int l=0; l < isochrons.size(); l++){
				 if(isochrons[l]->nbPoints > nbPointsPerIsochrons)	 k+=nbPointsPerIsochrons;
				 else k+=isochrons[l]->nbPoints;
			 }
		 }
		 printf("done\n");
		 printf("saved\n");
	 }

}
// La phase est calculée en prenant comme point de départ le premier point dans

double EnzymePFK::phaseOf(point3D * p3d, int numberMaxOfIsoPoints){
	double x = p3d->x;
	double y = p3d->y;
	double z = p3d->z;
	double xTemp = x;
	double yTemp = y;
	double zTemp = z;
	double result;
	bool found = false;
	int i = 0;
	int j = 0;
	double distance;
	double maxDistance;
	bool inSegment = false;

	point3D tmpPoint3D(x,y,z);
	while(not(found) and (i < numberOfIterations)){
		j=0;
		tmpPoint3D.x = x;
		tmpPoint3D.y = y;
		while(not(found) and j < attractor.size()){
			if(j == attractor.size()-2){
				maxDistance = distancePointSegment(attractor[j+1], attractor[j], attractor[1], &inSegment);
				distance = distancePointSegment(&tmpPoint3D, attractor[j], attractor[1], &inSegment);
			}
			else if(j == attractor.size()-1){
				maxDistance = distancePointSegment(attractor[1], attractor[j], attractor[2], &inSegment);
				distance = distancePointSegment(&tmpPoint3D, attractor[j], attractor[2], &inSegment);
			}
			else {
				maxDistance = distancePointSegment(attractor[j+1], attractor[j], attractor[j+2], &inSegment);
				distance = distancePointSegment(&tmpPoint3D, attractor[j], attractor[j+2], &inSegment);
			}
			//diffNorm = pow(attractor[j]->x -x,2.0) + pow(attractor[j]->y - y,2.0) + pow(attractor[j]->z - z,2.0);
			if(distance <= (maxDistance)) {
				found = true;
				break;
			}
			else{
				if(p3d->numberOfIsoPoints < numberMaxOfIsoPoints){
					p3d->isoPoints.push_back(new point3D(x, y, z));
					p3d->numberOfIsoPoints++;
				//printf("numberOfIsoPoints = %d\n", p3d->numberOfIsoPoints);
				}

			}
			j++;
		}
		for(int k=0; k < nbPointsPerPeriode -1; k++){
			L = x * (pow((1+x),(n-1)) * pow((1+D * y),n) + L0 * C * pow((1+C * x),(n-1))) / (pow((1+x),n) * pow((1+D * y),n) + L0 * pow((1+C * x),n));
			xTemp += (deltaT * (A-L));
			yTemp += (deltaT * (R *(L - N * y)));
			x = xTemp;
			y = yTemp;

		}
		z=0;
		i = i + nbPointsPerPeriode - 1;
		//printf("x = %f, y=%f, z = %f\n",x,y,z);
	}
	if(found){
		result = j*deltaT;
		while (result >= periode) result -= periode;
		//result = periode - result;
		return result;
	}
	else return -1;
}

double EnzymePFK::distancePointSegment(point3D * p, point3D * p1, point3D * p2, bool * inSegment){

	double result, a, d, ratio;
	double dx, dy, dz;


	if((p1->x == p2->x) and (p1->y == p2->y)){
		result = sqrt(pow(p->x - p1->x, 2.0) + pow(p->y - p1->y, 2.0) + pow(p->z - p1->z, 2.0));
		*inSegment = true;
	}
	else{
		dx = p2->x - p1->x;
		dy = p2->y - p1->y;
		dz = p2->z - p1->z;
		a = dx * (p->x - p1->x) + dy * (p->y - p1->y) + dz * (p->z - p1->z);
		d = sqrt(dx*dx + dy*dy + dz*dz);
		ratio = a / d;
		if(ratio < 0){
			result = sqrt(pow(p->x - p1->x, 2.0) + pow(p->y - p1->y, 2.0) + pow(p->z - p1->z, 2.0));
			*inSegment = false;
		}
		else if(ratio > 1){
			result = sqrt(pow(p->x - p2->x, 2.0) + pow(p->y - p2->y, 2.0) + pow(p->z - p2->z, 2.0));
			*inSegment = false;
		}
		else{
			result = sqrt(pow(p->x - p1->x - ratio * dx, 2.0) + pow(p->y - p1->y - ratio * dy, 2.0) + pow(p->z - p1->z - ratio * dz, 2.0));
			*inSegment = true;
		}
	}
	return result;
}
