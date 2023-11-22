// Interdiction de renvoyer un tableau dans une fonction !
// Le mettre en arguments.


/*-----------------------------------------------------------------------------------------------*/
/*											GIBBS												*/
/*-----------------------------------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


#include "Gibbs_fct_collect.h"
#include "Gibbs_fct_fill.h"

// pb entre lorsque nbr shif < W ?
// pb lors du shift gauche droite
// pb dans Ax des fois apparitions de Nan a cause de 0/0 normalement
#define NOMBRE_MAX_SEQUENCE 200
#define W 16 //Taille des segments a extraire
#define nbr_Aas 20
#define nbr_iterations 6000
#define nbr_shifts 10


int Gibbs_Sampling(char *Fichier_Fasta)
{	printf("%s\n","");
	/*Initialise la fonction srand sur le temps actuel
	Nb : Il faut l'appeler qu'une seule fois dans la fonction !*/
	srand(time(NULL));

	/*Init du nombre de sequence*/
	int N = 0;

	/*Tableau des 20 Acides amines*/
	char lesAas[nbr_Aas] = {'E','D','A','R','N','C','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V'};

	/*Initialisation du tableau contenant les sequences*/
	char **Liste_Sequences = (char**) malloc(NOMBRE_MAX_SEQUENCE*sizeof(char*));
	if (Liste_Sequences == NULL) //Verification de lallocation memoire
	    {printf("Probleme d'allocation\n"); //essayer memset(&Liste_Sequence) pour afficher le msg d'erreur
	    exit(0);} // Arret immediat du programme

	/*Initialisation du tableau contenant les titires des sequences*/
	char **Titres = (char**) malloc(NOMBRE_MAX_SEQUENCE*sizeof(char*));
	if (Titres == NULL)
	    {printf("Probleme d'allocation\n");
	    exit(0);}

	/*Nombre de sequences
	1 Sequence = la sequence + '\n' + '\0' */
	N = read_fasta(Liste_Sequences,Fichier_Fasta,Titres);
	printf("%s\n","");

	// Pour parcourir le tableau !
	// Deuxieme boucle pour verifier que l'on a bien les elements terminaux
	int i=0,j=0;
	while (Liste_Sequences[i] != NULL) {
	    printf("%s %i\n","Sequence",i);
	    for (j=0 ; j<strlen(Liste_Sequences[i]) ; j++) {
	        printf("%c",Liste_Sequences[i][j]);}printf("%s\n","");
	    i++;
	}
	/*Calcul de beta*/
	float Beta = sqrt(N);
	printf("\n%s %f\n","Beta est egale a : ",Beta);
	printf("%s\n","");

	/*Les segments de chaque sequence*/
	int *Liste_ind_segments = (int*)malloc(50*sizeof(int*));
	if (Liste_ind_segments == NULL) //Verification de lallocation memoire
	    {printf("Probleme d'allocation\n");
	    exit(0);} // Arret immediat du programme ;

	//Remplissage du tableau contenant ma liste de segments et celui contenant les sequences sans motifs
	remplissage_ind_segments(Liste_ind_segments,Liste_Sequences,W);
	// Pour parcourir le tableau !
	i=0;
	j=0;
	printf("%s%i %s\n","Il y a ",N,"segments : ");
	printf("%s\n","");
	//if (Liste_ind_segments[2] == NULL) { printf("%s","coucou");}
	//printf("%s",Liste_ind_segments[2]);
	while (Liste_Sequences[i] != NULL) {
	    printf("%s %i\n","Segment",i);
	    for (j = Liste_ind_segments[i] ; j < Liste_ind_segments[i] + W ; j++) {
	        printf("%c",Liste_Sequences[i][j]);}printf("%s\n","");
	    i++;
	}

/* A REPETER Jusqu'a nbr_iterations */
	int compteur = 0;
	int shift = 0; /*variable qui va compter le nombre de tour et des que shift == nbr_shifts => phases de decalage*/
	int taille_mat = 0;

	for(compteur=0;compteur<nbr_iterations;compteur++)
{
	//printf("iteration : %i\n",compteur);
	//printf("shift : %i\n",shift);
	int index_seq_ex = 0;
	index_seq_ex = nbre_al(0,N);

	if(compteur == 0){
	printf("\nIndice sequence exclue : %i\n",index_seq_ex);
	printf("%s %s","Sequence exclue :",Liste_Sequences[index_seq_ex]);
	printf("\nSegment de la sequence exclue : %.*s\n\n",W,Liste_Sequences[index_seq_ex]+Liste_ind_segments[index_seq_ex]);}

	/* Creation d'un tableau contenant les indices deplacees a gauche et a droite => pour eviter les min. locaux */

	int nshift = W;
	int *tab_index = malloc((nshift*2+1)*sizeof(int*));

	if (shift == nbr_shifts) // Toutes les nbr_shifts iterations : decalage gauche/droite
	{taille_mat = 0;
	if (tab_index == NULL) {printf("Probleme d'allocation\n"); exit(0);}
	taille_mat = remplissage_tab_index(tab_index,W,index_seq_ex,Liste_Sequences,Liste_ind_segments);
	}else
	{taille_mat = 1;
	tab_index[taille_mat-1] = index_seq_ex;} // Entre ces compteur iterations, un seul indice
	/*int m=0;
	for(m=0;m<taille_mat;m++){printf("Indice %i: %i\n",m,tab_index[m]);} */

	int iter = 0;
	float **tab_Ax = (float**)malloc(taille_mat*(strlen(Liste_Sequences[index_seq_ex])-W+1)*sizeof(float**));

	for(iter = 0; iter < taille_mat; iter++) /* On recalcule toutes les variables pour retrouver les Ax de chaque decalage sauf si taille_mat = 1 dans ce cas on calcul un seul set de Ax */
{
	/*
	printf("Iter %i\n",iter);
	printf("tab[Iter] %i\n",tab_index[iter]);*/

	//On remplace le nouvel indice dans la liste contenant les
	Liste_ind_segments[index_seq_ex] = tab_index[iter];

	/*Initialisation du tableau contenant les sequences sans leur motif*/
	char **Liste_Seq_sans_motifs = (char**) malloc(NOMBRE_MAX_SEQUENCE*sizeof(char*));
	if (Liste_Seq_sans_motifs == NULL)
	    {printf("Probleme d'allocation\n");
	    exit(0);}
	remplissage_seg_sans_motifs(Liste_Sequences,Liste_Seq_sans_motifs,Liste_ind_segments,W);
	// Pour parcourir le tableau !
	i=0;
	if(compteur == 0){
	while (Liste_Seq_sans_motifs[i] != NULL) {
	    printf("%s %i\n","Sequence sans motif",i);
	    for (j=0 ; j<strlen(Liste_Seq_sans_motifs[i]); j++) {
	        printf("%c",Liste_Seq_sans_motifs[i][j]);}printf("%s\n","");
	    i++;
	}
	printf("%s\n"," ");}
	//Initialisation du tableau Cij pour les motifs
	int Cij_seg[nbr_Aas][W];
	for (i=0; i<nbr_Aas; i++) {
	    for (j=0; j<W; j++) {
		Cij_seg[i][j] = 0;
	}}

	//Remplissage de Cij pour les motifs sauf de la sequence exclue
	remplissage_Cij_seg(nbr_Aas, W, Cij_seg, Liste_ind_segments, lesAas, index_seq_ex,Liste_Sequences, N);
	// Pour parcourir le tableau !

	if(compteur == 0){
	printf("%s\n","Le tableau Cij pour les segments (sauf de la sequence exclue) : ");
	//printf("%s %i","par ici camarade",Cij_seg[10][9]);
	printf("%s\n","");
	for (i=0; i<nbr_Aas; i++) {
	    printf("%c ",lesAas[i]);
	    for (j=0; j<W; j++) {
		printf("%i",Cij_seg[i][j]);
	}
	printf("%s\n"," ");}}

	//Initialisation du tableau Cj hors les motifs
	int Cj_hors_motifs[nbr_Aas];
	for (i=0; i<nbr_Aas; i++) {
	    Cj_hors_motifs[i] = 0.0;
	}

	//Remplissage de Cij pour les motifs
	remplissage_Cj_hors_motifs(nbr_Aas,Cj_hors_motifs, Liste_Seq_sans_motifs, lesAas);

	if(compteur == 0){
	printf("%s\n","");
	printf("%s\n","Le tableau Cj pour toutes les sequences hors blocs : ");
	printf("%s\n","");
	for (i=0; i<nbr_Aas; i++) {
	    printf("%c %i",lesAas[i],Cj_hors_motifs[i]);
	printf("%s\n"," ");}
	}
	//Initialisation du tableau bj
	float bj[nbr_Aas];
	for (i=0; i<nbr_Aas; i++) {
	    bj[i] = 0.0;
	}
	remplissage_bj(nbr_Aas, N, bj, Liste_Sequences, lesAas);

	// Pour parcourir le tableau !
	if(compteur == 0){
	printf("%s\n"," ");
	i=0;
	printf("%s\n","Le tableau bj pour toutes les donnees : ");
	printf("%s\n","");
	for (i=0; i<nbr_Aas; i++) {
	    printf("%c %.2f\n",lesAas[i],bj[i]);
	}
	printf("%s\n"," ");}

	/*Initialisation du tableau des Qij*/
	float Qij[nbr_Aas][W];
	for (i=0; i<nbr_Aas; i++) {
	    //printf("%i ",i);
	    for (j=0; j<W; j++) {
		Qij[i][j] = 0;
	}}

	//Remplissage de Qij pour les motifs
	remplissage_Qij(nbr_Aas, W, N, Qij, Cij_seg, bj);

	// Pour parcourir le tableau !
	if(compteur == 0){
	printf("%s\n","");
	printf("%s\n","Le tableau Qij : ");
	printf("%s\n","");
	printf("%s%4s","Position ",":");
	for (j=0; j<W; j++) { printf("%11i", j);}
	printf("%s\n","");
	for (i=0; i<nbr_Aas; i++) {
	    printf("%s %c ","Acide amine :",lesAas[i]);
	    for (j=0; j<W; j++) {
		printf("%10.3f",Qij[i][j]);
	}
	printf("%s\n"," ");}
	printf("%s\n"," ");
}
	//Calculs des Pj dans les sequences hors motif
	//Initialisation des Pj hors motifs
	float Pj[nbr_Aas];
	for (i=0; i<nbr_Aas; i++)
	{
	    Pj[i] = 0.0;
	}
	//Remplissage de Pj
	remplissage_Pj(nbr_Aas,N,Pj,bj, Cj_hors_motifs,lesAas);

	//Affichage de Pj
	if(compteur < 10){
	printf("%s\n","Le tableau Pj pour les sequences hors motfis : ");
	printf("%s\n","");
	for (i=0; i<nbr_Aas; i++) {
	    printf("%c %f\n",lesAas[i],Pj[i]);
	}printf("\n");}

	//Initialisation des Ax sur la sequence exclue
	float Ax[strlen(Liste_Sequences[index_seq_ex])-W+1];
	for (i=0; i<strlen(Liste_Sequences[index_seq_ex])-W+1; i++)
	{Ax[i] = 1.0;}

	//Remplissage des Ax
	remplissage_Ax(Liste_Sequences[index_seq_ex], nbr_Aas, W, Qij, Pj, Ax, lesAas);

	//Affichage des Ax
	if(compteur < 10){
	printf("%s\n","Le tableau des Ax : \n");
	for (i=0; i<strlen(Liste_Sequences[index_seq_ex])-W+1; i++) {
	    printf("indice %i : %5.5f\n",i,Ax[i]);}}

	remplissage_Ax_normalise(Liste_Sequences[index_seq_ex], W, Ax);
	if(compteur == 0){
	printf("%s\n","");
	printf("%s\n","Le tableau des Ax normalises : \n");
	for (i=0; i<strlen(Liste_Sequences[index_seq_ex])-W+1; i++) {
	    printf("indice %i : %5.5f\n",i,Ax[i]);
	}}
	//Choix du Ax
	int ind_choisi = 0;
	ind_choisi = choix_Ax(compteur, Liste_Sequences[index_seq_ex], W, Ax);

	//On remplace l ancien indice par le nouveau choisi
	Liste_ind_segments[index_seq_ex] = ind_choisi;

	if(compteur == 0 ){
	printf("\nSequence exclue (indice) : %s (%i)\n",Liste_Sequences[index_seq_ex],index_seq_ex);
	printf("Le nouvel indice du motif de la sequence exclue est : %i\n",ind_choisi);
	printf("Le nouveau motif de la seq exclue est : %.*s\n\n",W,Liste_Sequences[index_seq_ex]+Liste_ind_segments[index_seq_ex]);

	i = 0;
	while (Liste_Sequences[i] != NULL) {
		if (i != index_seq_ex) { printf("%s %i\n","Ancien segment",i);}
		else { printf("%s %i\n","Nouveau segment",i);}
	    for (j = Liste_ind_segments[i] ; j < Liste_ind_segments[i] + W ; j++) {
	        printf("%c",Liste_Sequences[i][j]);}printf("%s\n","");
	    i++;
	}	printf("\n");}

	i = 0;

	if(compteur == nbr_iterations-1){
	printf("\n");
	printf("%s\n","Les motifs convergents : ");
	while (Liste_Sequences[i] != NULL) {
	printf("de la sequence %i : ",i);
	    for (j = Liste_ind_segments[i] ; j < Liste_ind_segments[i] + W ; j++) {
	        printf("%c",Liste_Sequences[i][j]);}printf("%s\n","");
	    i++;
	}}

 if(shift == nbr_shifts)
	{
	remplissage_tab_Ax(iter,W,tab_Ax,Liste_Sequences[index_seq_ex],Ax);}

 if(shift == nbr_shifts && iter == taille_mat-1){ /* == A la derniere iteration */
	remplissage_tab_Ax_normalise(iter,Liste_Sequences[index_seq_ex], W, tab_Ax);
	//Choix du Ax
	int ind_choisi = 0;
	ind_choisi = tab_index[choix_Ax_tab(iter, Liste_Sequences[index_seq_ex], W, tab_Ax)];
	Liste_ind_segments[index_seq_ex] = ind_choisi;
	free(tab_Ax);/* On libere le tableau des Ax */
	shift =0 ;
	}
}//boucle iter
	shift++;

} // boucle compteur

	// A faire a la fin
	/* Création d'un dossier qui va contenir le fichier de sortie */
    make_dir_file(N,W,Liste_Sequences,Liste_ind_segments,Titres);
	/* On libèere l'espace allouee aux sequences et titres */
    for(i=0; i<N; i++){free(Liste_Sequences[i]);} free(Liste_Sequences);
    for(i=0; i<N; i++){free(Titres[i]);} free(Titres);
    free(Liste_ind_segments);
    return 0;
}
int main()
{
	/* Le nom du fichier format fasta qui contient les sequences */
	char *namefile = "lipocalin.fst";
	//char *namefile = "ExDuProjet.fst"; //Dans lexemple il y a des B (pas un Aa => remplace ici par R
	//char *namefile = "TestSeq2.fasta";

	/* Renvoie le fichier motif contenu dans le dossier ./Gibbs/motifs_convergents/, contenant les motifs convergeants de ces sequences */

	Gibbs_Sampling(namefile);
	//UPGMA();

	return 0;
}
