#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int nbre_al(int borne_inf, int borne_sup) {

	int randomNumber = 0;
	randomNumber = rand()%(borne_sup-borne_inf) + borne_inf;
	return randomNumber;
}

float nbre_al_flottant(int borne_inf, int borne_sup) {

	float randomNumber = 0;
	randomNumber = rand()/(double)RAND_MAX*(borne_sup-borne_inf) + borne_inf;
	return randomNumber;
}


int choix_Ax(int compteur, char seq_exclue[], int taille_segment,float tab_Ax[strlen(seq_exclue)-taille_segment])
{
	float choix = 0.0;
	int ind_choisi = 0;
	float cumsum_Ax = tab_Ax[ind_choisi];

	choix = nbre_al_flottant(0,1);
	while(cumsum_Ax < choix){
	ind_choisi++;
	cumsum_Ax += tab_Ax[ind_choisi];}

	if(compteur==0){
	printf("\nNombre_Aleatoire : %.4f\nSomme_Cumulee : %.4f\nIndice_selectionnee : %i\n",choix,cumsum_Ax,ind_choisi);}
	return ind_choisi;
}

int choix_Ax_tab(int shifts, char seq_exclue[], int taille_segment,float **tab_Ax)
{
	int i,j;
	float choix = 0.0;
	int ind_choisi = 0;
	float cumsum_Ax = 0.0;
	int stop = 0;

	choix = nbre_al_flottant(0,1);
	//printf("Choix = %f\n",choix);

	while(cumsum_Ax < choix){
	for (i = 0; i<shifts && !stop; i++){
		for (j = 0; j<strlen(seq_exclue)-taille_segment+1 && !stop; j++){
			if(cumsum_Ax > choix){stop=1;}
			cumsum_Ax += tab_Ax[i][j];//printf("Somme cumu : %i %f\n",i,cumsum_Ax);
			 }}}
	ind_choisi = i-1;
	//printf("ind_choisi = %i\n",ind_choisi);
	return ind_choisi;
}
