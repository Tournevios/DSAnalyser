
/*
 * SmartPencil.cpp
 *
 *  Created on: 19-Jan-2009
 *      Author: hbenamor
 */

#include "CyclicAttractorSmartPencil.h"
#include "random-singleton.h"
#include "MersenneTwister.h"
#include "point3D.h"
#include <time.h>

SmartPencil::SmartPencil() {
	// TODO Auto-generated constructor stub

}

SmartPencil::~SmartPencil() {
	// TODO Auto-generated destructor stub
}

void SmartPencil::initialize(CyclicAttractor* cyclicAttractor, Isochrons * isochrons){
	this->isochrons = isochrons;
		this->cyclicAttractor = cyclicAttractor;
		nbDimension = cyclicAttractor->nbDimension;
		nb_points_found = 0;
		nb_max_dep = 10; //2
		nb_dep = 0;
		/*center[0] = isochrons->points[0]->x;
		center[1] = isochrons->points[0]->y;
		center[2] = isochrons->points[0]->z;*/
		center[0] = 0;
		center[1] = 0;
		center[2] = 0;
		hypothenus = 0;
		//deltaDep = 0;
		//Modification pour boxon
		//cote[0] = 600;
		//cote[1] = 600;
		//cote[2] = 600;

		for(int i=0; i < nbDimension; i++){
			cote[i] = cyclicAttractor->radius[i];
			deltaDep[i] = cote[i]/2.0;
			hypothenus += deltaDep[i] * deltaDep[i];
			//deltaDep += ((cote[i]/2.0) * (cote[i]/2.0));
		}
		hypothenus = sqrt(hypothenus);
		//deltaDep = sqrt(deltaDep);
		//cote = 0.01;
		//deltaDep = 0.005;
		nbMaxPointsPerStep = 50; //50
		headPencil = true;
		secondPencil = NULL;
		mtRand = new MTRand(isochrons->myIndex +time(NULL));
}

/*
 * To pull points in a state window
 */
void SmartPencil::setWindow(double window[]){
	center[0] = (window[0]+window[1])/2;
	center[1] = (window[2]+window[3])/2;
	center[2] = (window[4]+window[5])/2;
	for(int i=0; i < nbDimension; i++){
		cote[i] = fabs(window[i*2]-window[i*2+1]);
		deltaDep[i] = cote[i]/2.0;
		hypothenus += deltaDep[i] * deltaDep[i];
	}
	hypothenus = sqrt(hypothenus);

}

SmartPencil::SmartPencil(CyclicAttractor* cyclicAttractor, Isochrons * isochrons){
	initialize(cyclicAttractor, isochrons);
}

SmartPencil::SmartPencil(CyclicAttractor* cyclicAttractor, Isochrons * isochrons, double * oldOne){
	initialize(cyclicAttractor, isochrons);
	setOldOne(oldOne);
}

void * runMeToo(void * arg){
	SmartPencil * pencil2 = (SmartPencil*) arg;
	if(pencil2 != NULL)	pencil2->runMe();
	else printf("no object\n");
	return NULL;
}

/*void SmartPencil::neighbourhoodShots(double * center, double epsilon, double * point){
// Pour tirer aux alentours de center selon un epsilon pret
	point[0] = mtRand->rand() * 2 * epsilon - epsilon;
	point[1] = mtRand->rand() * 2 * epsilon - epsilon;
	point[2] = 0;
	point[0] += center[0];
	point[1] += center[1];
}*/
void * SmartPencil::runMeOnNeighbourhood(double * epsilon){
	double point[255];
	point3D p3d(0,0,0,0,0,0,0);
	int indexOfIsochrons;
	int myIndex;
	int coeffCentre[3];
	double tmpPhase = -1;
	double phasePoint[3];

	phasePoint[0] = isochrons->phasePoint->x;
	phasePoint[1] = isochrons->phasePoint->y;
	phasePoint[2] = isochrons->phasePoint->z;

		// Boucle de recherche des points appartenant au même isochrons
		while(nb_dep < nb_max_dep){
			p3d.reInitialize();
			// Générer un point en faisant un tire de Monte-Carlo aux alentours du centre à un epsilon pret
			point[0] = mtRand->rand() * 2 * epsilon[0] - epsilon[0];
			point[1] = mtRand->rand() * 2 * epsilon[1] - epsilon[1];
			point[2] = mtRand->rand() * 2 * epsilon[2] - epsilon[2];
			point[0] += phasePoint[0];
			point[1] += phasePoint[1];
			point[2] += phasePoint[2];

			// Récupérer les 3 coordonnées (dimension 3)
			p3d.x = point[0];
			p3d.y = point[1];
			p3d.z = point[2];

			// Identifier l'isochrone d'appartenance du point
			indexOfIsochrons = cyclicAttractor->seekIsochronsAndAdd(&p3d, &tmpPhase, true);
			//indexOfIsochrons = cyclicAttractor->seekIsochronsAndAdd2(&p3d, true);
			if(indexOfIsochrons == isochrons->myIndex) {nb_dep++; printf("found\n");}
			else if(indexOfIsochrons != -1) printf("Isochron %d\n", indexOfIsochrons);
		}
		return NULL;
}

void * SmartPencil::runMeNoDep(){
	double point[255];
	point3D p3d(0,0,0,0,0,0,0);
	int indexOfIsochrons;
	int myIndex;
	int coeffCentre[3];
	double tmpPhase = -1;
		/*
		 * Pour tirer dans la 27 régions autour de l'attracteur
		 */
		myIndex = isochrons->myIndex % 27;
		if(myIndex == 13) myIndex++; // testme

		for(int i=0; i < nbDimension; i++){
			point[i] = 0;
			coeffCentre[i] = myIndex / (int)pow(3,2-i);
			myIndex -= coeffCentre[i] * pow(3,2-i);
			coeffCentre[i] -= 1;
		}
		point[2] = 0;

		// Boucle de recherche des points appartenant au même isochrons
		while(nb_dep < nb_max_dep){
			p3d.reInitialize();
			// Générer un point en faisant un tire de Monte-Carlo
			for(int i=0; i < nbDimension; i++){
				//point[i] = Random::Uniform(-cote[i], cote[i]);
				//point[i] = cyclicAttractor->getRandomNum(-cote[i], cote[i]);
				// Plus rapide que Lucky luke
				//point[i] = mtRand->rand() * 2 * cote[i] - cote[i] + coeffCentre[i] * cote[i];
			/*	if(cote[i] == 0)
					point[i] = mtRand->rand() * 2 * 10 - 10;
				else*/
					point[i] = mtRand->rand() * 2 * cote[i] - cote[i] + center[i];
			}

			// Modification pour boxon
			//point[2] += 200;

			// Récupérer les 3 coordonnées (dimension 3)
			p3d.x = point[0];
			p3d.y = point[1];
			p3d.z = point[2];

			// Identifier l'isochrone d'appartenance du point
			indexOfIsochrons = cyclicAttractor->seekIsochronsAndAdd(&p3d, &tmpPhase, true);
			//indexOfIsochrons = cyclicAttractor->seekIsochronsAndAdd2(&p3d, true);
			if(indexOfIsochrons == isochrons->myIndex) {
				nb_dep++;
				//printf("found %f%%\n",((double)nb_dep*100/nb_max_dep));
			}
			//else if(indexOfIsochrons != -1) printf("Isochron %d\n", indexOfIsochrons);
		}
		return NULL;
}

void * SmartPencil::runMe(){
	double point[255];
	point3D p3d(0,0,0);
	double tmpPhase;
	int indexOfIsochrons;
	double projection;
	double oldOne[3];
	bool nullVector;
	pthread_t * secondThread; // = new pthread_t();


	if(headPencil){
		secondThread = new pthread_t();
		if(pthread_create(secondThread, NULL, runMeToo, (void*) secondPencil) != 0){
			 printf("Thread %d-%d creation failed: %d\n", isochrons->myIndex, 2, 0);
		}
	}

	for(int i=0; i < nbDimension; i++){
		average_vector[i] = 0;
		point[i] = 0;
		oldOne[i] = 0;
		//oldOne[i] = direction * deltaDep[i];
		//norme += oldOne[i] * oldOne[i];
	}


	// Boucle de recherche des points appartenant au même isochrons
	int count = 0;
	while(nb_dep < nb_max_dep){
		p3d.reInitialize();
		// Générer un point en faisant un tire de Monte-Carlo
		for(int i=0; i < nbDimension; i++){
			point[i] = cyclicAttractor->getRandomNum(center[i] - (cote[i]/2.0), center[i] + (cote[i]/2.0));
		}

		// Récupérer les 3 coordonnées (dimension 3)
		p3d.x = point[0];
		p3d.y = point[1];
		p3d.z = point[2];

		count++;
		// Identifier l'isochrone d'appartenance du point
		indexOfIsochrons = cyclicAttractor->seekIsochronsAndAdd(&p3d, &tmpPhase, true);
		// Traitement dans le cas ou l'isochrone d'appartenance est celui du pinceau
		if(indexOfIsochrons == isochrons->myIndex){
			nb_points_found++;
			// Si le point se situe sur la direction de déplacement
			// alors calculer la prochaine trajectoire
			double proj = 0.0;
			for(int i=0; i < nbDimension; i++){
				proj += (point[i]-center[i])*oldOne[i];
			}
			if(proj >= 0){
				for(int i=0; i < nbDimension; i++){
					average_vector[i] += (point[i] - center[i])/nbMaxPointsPerStep;
				}
			}
		}
		/*
		 * Si le nombre de points trouvés dépasse le nombre de points
		 * autorisés par petit cube alors
		 * bouger le cube d'une unité suivant le vecteur moyen
		 *
		 */
		if((nb_points_found >= nbMaxPointsPerStep) or (count > 10 * nbMaxPointsPerStep)){
			// Calcul de la norme du vecteur moyen
			nullVector = normalizeIt(average_vector);
			for(int i=0; i < nbDimension; i++){
					average_vector[i] += oldOne[i]*0.25;
			}

			// Calculer la nouvelle norme et renormaliser le vecteur déplacement résultant
			normalizeIt(average_vector);
			// Tester si le pinceau tente de repartir dans le sens contraire
			// Et eviter qu'il s'arrete

			projection = 0;
			for(int i=0; i < nbDimension; i++){
				projection += average_vector[i] * oldOne[i];
			}
			//printf("sens %d = %f\n", isochrons->myIndex, projection);
			/*
			if(projection < 0){ // Empeche le pinceau de repartir dans le sens contraire
				for(int i=0; i < nbDimension; i++){
					average_vector[i] = - average_vector[i];
				}
				printf("ca bouge!!!\n");
			}*/
			/*
			else if(projection ==0){ // Empeche le pinceau de s'arrêter
				for(int i=0; i < nbDimension; i++){
					average_vector[i] = oldOne[i]+ pow(deltaDep, 1.0/nbDimension) * pow(-1.0, i);
				}
				printf("glandage\n");
			}*/
			double ratio = (double)(nb_points_found/(double)nbMaxPointsPerStep) * 100;
			if(ratio < 8){
				printf("pencil %d %d dep number %d ratio=%f%% x=%f y=%f \n", isochrons->myIndex , headPencil, nb_dep, ratio, average_vector[0], average_vector[1]);
				for(int i=0; i < nbDimension; i++){
					center[i] -= (oldOne[i] * deltaDep[i]);
					average_vector[i] = 0;
					oldOne[i] = 0;
					cote[i] = cote[i] * (1 - cyclicAttractor->getRandomNum(0.05,0.1));
					deltaDep[i] = cote[i] / 2.0;
				}
				initializeOldPoint();
				normalizeIt(oldOne);
			}
			else if(ratio > 50){
				for(int i=0; i < nbDimension; i++){
					cote[i] = cote[i] * (1 + cyclicAttractor->getRandomNum(0.05,0.1));
					deltaDep[i] = cote[i] / 2.0;
					center[i] += (average_vector[i] * deltaDep[i]);
					oldOne[i] = average_vector[i];
					average_vector[i] = 0;
				}
				nb_dep++;
			}
			else{
				for(int i=0; i < nbDimension; i++){
					//Calcul du nouveau centre
					//cote[i] = cote[i] * (1 - cyclicAttractor->getRandomNum(0.01,0.05));
					//deltaDep[i] = cote[i] / 2.0;
					center[i] += (average_vector[i] * deltaDep[i]);
					oldOne[i] = average_vector[i];
					average_vector[i] = 0;
				}
				nb_dep++;
			}
			nb_points_found = 0;
			count = 0;
		}
	}
	if(headPencil)	pthread_join(*secondThread, NULL);
	return NULL;
}

void SmartPencil::initializeOldPoint(){

	/*
	 * Methode d'exploration du voisinage
	 *
	 */
	double point[255];
	double ratios[255];
	double last_max = 0;
	int secondMax = -1;
	int maxIndex = -1;
	bool nullVector;
	int initialMaxShots = 100;
	int nbInitialVectors = (int)(pow(3.0, nbDimension));
	double initialVector[255][255];
	point3D p3d(0,0,0);
	double tmpPhase;
	int indexOfIsochrons;
	double norme = 0;

	// Algorithme exhaustif de construction de tout les vecteurs radiaux

	for(int j=0; j < nbDimension; j++){
		initialVector[0][j] = -1;
	}
	int tmpNumberOfVectors = 1;
	int progIndex = 0;
	int count = 0;
	while (progIndex < nbDimension){
		count = 0;
		for(int i=0; i < tmpNumberOfVectors; i++){
			for(int j=0; j<nbDimension;j++){
				if(j==progIndex){
					initialVector[tmpNumberOfVectors+count][progIndex] = 0;
					initialVector[tmpNumberOfVectors+count+1][progIndex] = 1;
				}
				else{
					initialVector[tmpNumberOfVectors+count][j] = initialVector[i][j];
					initialVector[tmpNumberOfVectors+count+1][j] = initialVector[i][j];
				}
			}
			count += 2;
		}
		tmpNumberOfVectors += tmpNumberOfVectors*2;
		progIndex++;
	}

	// normalisation des vecteurs radiaux
	for(int i=0; i < nbInitialVectors; i++){
		normalizeIt(initialVector[i]);
	}

	if(isochrons->myIndex ==0 and direction==1){
		for(int i=0; i < nbInitialVectors; i++){
			printf("Vector %d %d x=%f y=%f z=%f:\n", isochrons->myIndex, i, initialVector[i][0], initialVector[i][1], initialVector[i][2]);
		}
	}

	for(int k=0; k < nbInitialVectors; k++){
		ratios[k] = 0;
		nullVector = normalizeIt(initialVector[k]);;

		point[2] = 0;
		if(not nullVector){
			// Séléction du vecteur candidats
			int count = 0;
			nb_points_found = 0;
			while(count < initialMaxShots){
				for(int i=0; i < nbDimension; i++){
					point[i] = cyclicAttractor->getRandomNum(center[i] - initialVector[k][i] * (cote[i]/2.0), center[i] + initialVector[k][i] * (cote[i]/2.0));
				}

				p3d.x = point[0];
				p3d.y = point[1];
				p3d.z = point[2];
				indexOfIsochrons = cyclicAttractor->seekIsochronsAndAdd(&p3d, &tmpPhase, false);
				count++;
				if(indexOfIsochrons == isochrons->myIndex){
					nb_points_found++;
				}
				// Si on trouve le nombre de points demandés pour cet isochrons
				// alors retenir le vecteur le plus adéquat
				if(count == initialMaxShots){
					ratios[k] = (double)nb_points_found / (double)initialMaxShots;
					if(ratios[k] > last_max){
						last_max = ratios[k];
						maxIndex = k;
					}
				}
			}
		}
	}
	// Chercher le 2ème plus grand vecteur et crée le 2ème pinceau
	last_max = 0;
	secondMax = -1;
	for(int i = 0; i < nbInitialVectors; i++){
		// Retenir la plus grande projection
		if(i != maxIndex){
			if(ratios[i] > last_max){
				last_max = ratios[i];
				secondMax = i;
			}

		}
	}
	if((secondMax!=-1) and (secondPencil == NULL)) secondPencil = new SmartPencil(cyclicAttractor, isochrons, initialVector[secondMax]);

	norme = 0;
	if(headPencil){
		for(int j=0; j < nbDimension; j++){
			oldOne[j] = initialVector[maxIndex][j];
		}
	}
	else{
		for(int j=0; j < nbDimension; j++){
			oldOne[j] = initialVector[secondMax][j];
		}
	}

	printf("#################################\n");
	printf("vector %d:ratio = %f%%  %f  %f\n", isochrons->myIndex, ratios[maxIndex] * 100, oldOne[0], oldOne[1]);
	printf("second vector %d:ratio = %f%%  %f  %f\n", isochrons->myIndex, ratios[secondMax] * 100, initialVector[secondMax][0], initialVector[secondMax][1]);
	printf("#################################\n");
	/*
	 * Fin de la méthode d'exploration du voisinage
	 */

/*
	while(score < 0.9 * hypothenus){

		// Générer les candidats et initialiser les variables associées
		for(int i=0; i < nbInitialVectors; i++){
			norme = 0;
			initialProjection[i] = 0;
			// N'est admis qu'un vecteur non nul!!
			while(norme == 0){
				for(int j=0; j < nbDimension; j++){
						initialVector[i][j] = cyclicAttractor->getRandomNum(-1.0,1.0);
						norme += initialVector[i][j] * initialVector[i][j];
				}
				norme = sqrt(norme);
				if(norme != 0){
					for(int j=0; j < nbDimension; j++){
						initialVector[i][j] = initialVector[i][j]/norme;
					}
				}
			}
		}
		// Fin d'initialisation
		norme = 0;
		point[2] = 0;

		// Séléction du vecteur candidats
		while(nb_points_found < nbMaxPointsPerStep){
			for(int i=0; i < nbDimension; i++){
				point[i] = Random::Uniform(center[i] - (cote[i]/2.0), center[i] + (cote[i]/2.0));
			}

			p3d.x = point[0];
			p3d.y = point[1];
			p3d.z = point[2];
			indexOfIsochrons = cyclicAttractor->seekIsochronsAndAdd(&p3d, &tmpPhase);

			if(indexOfIsochrons == isochrons->myIndex){
				nb_points_found++;
				for(int i=0; i < nbInitialVectors; i++){
					// Calcul de la projection du point sur le vecteur
					// et retenir la projection maximale
					double tmpProj = 0;
					for(int j=0; j < nbDimension; j++){
						tmpProj += fabs((point[j]-center[j])*initialVector[i][j]);
					}
					if(tmpProj > initialProjection[i]) initialProjection[i] = tmpProj;
				}
			}
			// Si on trouve le nombre de points demandés pour cet isochrons
			// alors retenir le vecteur le plus adéquat
			if(nb_points_found == nbMaxPointsPerStep){
				int indexOfMax = 0;
				for(int i = 0; i < nbInitialVectors; i++){
					// Retenir la plus grande projection
					if(initialProjection[i] > maxProjection){
						maxProjection = initialProjection[i];
						indexOfMax = i;
						score = maxProjection;
						theMax = i;
						for(int j=0; j < nbDimension; j++){
							oldOne[j] = direction * initialVector[theMax][j];
						}
					}
				}
				// Fin de sélection pour cette génération
				nb_points_found = 0;
				//nb_dep++;
				break;
			}
		}
	}
	printf("Final score %d = %f %%\n", this->isochrons->myIndex, score*100/(hypothenus));
	for(int j=0; j < nbDimension; j++){
		oldOne[j] = direction * initialVector[theMax][j];
	}
	printf("vector %d: %f  %f\n", isochrons->myIndex, oldOne[0], oldOne[1]);
*/
}

/*
 * Set the initial vector
 *
 */
void SmartPencil::setOldOne(double * oldOne){
	for(int i=0; i < nbDimension; i++){
		this->oldOne[i] = oldOne[i];
	}
}

/*
 * Take a vector and normalize it
 * Return false if the norm is equal to zero
 * else return true
 */
int SmartPencil::normalizeIt(double * myVector){
	double norme = 0;
	for(int i=0; i < nbDimension; i++){
		norme += myVector[i] * myVector[i];
	}
	norme = sqrt(norme);
	if(norme != 0){
		for(int i=0; i < nbDimension; i++){
				myVector[i] = myVector[i] / norme;
		}
		return 0;
	}
	else return 1;
}
