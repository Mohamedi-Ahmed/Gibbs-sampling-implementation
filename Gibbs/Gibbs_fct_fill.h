#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/* Pour la creation du dossier */
#include <sys/types.h>
#include <sys/stat.h>

void remplissage_ind_segments(int *tab, char **Liste_Sequences,int taille_seg)
{
	int i=0,index_deb=0;
	//char dest[TAILLE_MAX_SEGMENT] = "";

	while (Liste_Sequences[i] != NULL)
	{
		index_deb = 0;
	        index_deb = nbre_al(0,strlen(Liste_Sequences[i])-taille_seg);
		//printf("ici %c\n",Liste_Sequences[i][strlen(Liste_Sequences[i])-2]);
		//tab[i] = malloc(sizeof(int));
		//memset(tab[i],'0',sizeof(tab));
		tab[i] = index_deb;
		//memset(dest,'\0',sizeof(dest));
		//strncpy(dest,Liste_Sequences[i],index_deb);
		//strncat(dest,Liste_Sequences[i]+index_deb+taille_seg,strlen(Liste_Sequences[i])-index_deb-taille_seg);
		//strncpy(tab[i],(Liste_Sequences[i]+index_deb),taille_seg);
		//tab[i][taille_seg] = '\0';
	        i++;
	        }


}

void remplissage_seg_sans_motifs(char **Liste_Sequences,char **Seq_sans_motif, int *Liste_ind_seg, int taille_seg)
{
	int i = 0, index_deb = 0;


	while (Liste_Sequences[i] != NULL)
	{
		char dest[TAILLE_MAX_SEGMENT] = "";
		index_deb = Liste_ind_seg[i];
		Seq_sans_motif[i] = malloc(strlen(Liste_Sequences[i])-taille_seg);
		memset(Seq_sans_motif[i],'a',sizeof(Seq_sans_motif));

		strncpy(dest,Liste_Sequences[i],index_deb);
		strncat(dest,Liste_Sequences[i]+index_deb+taille_seg,strlen(Liste_Sequences[i])-index_deb-taille_seg);
		strcpy(Seq_sans_motif[i],dest);
		//printf("\ncoucou : %s",dest);
		i++;
	}

}


void remplissage_Cij_seg(int nb_Aas,int taille_segment,int tab_Cij[nb_Aas][taille_segment], int *tab_seg,char lesAas[], int id_exclue, char **LesSeq, int nbr_seg)
{
	int i = 0, compt=0, loc_Aa = 0, j = 0, pos = 0;
	char Aa;
	while (i != nbr_seg)
	{
	    if(i != id_exclue) {
		pos = 0;
	        for (j = tab_seg[i]; j < tab_seg[i]+taille_segment; j++)
		    {
	            Aa = LesSeq[i][j];
		    for(compt = 0; compt < nb_Aas; compt++)
		    	if (Aa == lesAas[compt]) {
			    loc_Aa = compt;}
		    //printf("\n ici %c %i %i",Aa,loc_Aa,pos);
	            tab_Cij[loc_Aa][pos] += 1;
		    pos++;
	            }

	    }i++;
	}
}

void remplissage_Cj_hors_motifs(int nb_Aas,int tab_Cj[nb_Aas], char **tab_hors_seg,char lesAas[])
{
	int i = 0, j = 0;
	int compt=0, loc_Aa = 0, Aa = 0;
	while (tab_hors_seg[i] != NULL)
	{
	    for (j=0; j<strlen(tab_hors_seg[i]); j++)
		{

	        Aa = tab_hors_seg[i][j];
		//printf("Salut %c\n",tab_hors_seg[i][j]);

		for(compt = 0; compt < nb_Aas; compt++)
		    if (Aa == lesAas[compt]) {
		        loc_Aa = compt;}

		//printf("coucou %i\n",compt);
	        tab_Cj[loc_Aa] += 1;
	        }

	i++;
	}
}
void remplissage_bj(int nbr_Aas, int nbr_seq, float tab_bj[nbr_Aas], char **tab_seq, char lesAas[])
{
	int i = 0, j = 0, Aa = 0, compt=0, nb_Tot = 0, loc_Aa = 0;

	while (tab_seq[i] != NULL)
	{
		nb_Tot += strlen(tab_seq[i]);
		//printf("par ici %li\n",strlen(tab_seq[i]));
	        for (j=0; j<strlen(tab_seq[i]); j++)
		    {
	            Aa = tab_seq[i][j];

		for(compt = 0; compt < nbr_Aas; compt++)
		    if (Aa == lesAas[compt]) {
		        loc_Aa = compt;}

	            tab_bj[loc_Aa] += 1.0;
	           ;}//printf("\n");

	    i++;
	    }

	for(i=0; i<nbr_Aas; i++)
	{
	    //printf("COUCOU %f %i %f\n",tab_bj[i],nb_Tot,sqrt(nbr_seq));
	    tab_bj[i] = sqrt(nbr_seq)*tab_bj[i]/(nb_Tot);
	}
}

void remplissage_Qij(int nb_Aas, int taille_segment, int nbr_seq, float Qij[nb_Aas][taille_segment], int tab_Cij[nb_Aas][taille_segment], float tab_bj[nb_Aas])
{
	float num = 0, denom = nbr_seq-1+sqrt(nbr_seq);
	int i=0,j=0;
	//printf("%s %i\n","Look here",tab_Cij[10][9]);
	for(j = 0; j < nb_Aas; j++)
	{
	    for (i=0; i < taille_segment; i++)
	    {
	        //Qij[i]) = malloc(sizeof(int));
	        //memset(Qij[i],'\0',sizeof(Qij[i]));
		num = tab_Cij[j][i]+tab_bj[j];
		//printf("PAR LA ! j = %i i = %i : %i %f %f\n",j,i,tab_Cij[j][i],tab_bj[j],num);
	        Qij[j][i] = num/denom;
	    }
	}
}

void remplissage_Pj(int nb_Aas,int N, float tab_Pj[nb_Aas],float tab_bj[nb_Aas], int tab_Cj_hors[nb_Aas], char lesAas[])
{
	int j, k, nb_Res_Hors;
	float sum = 0.0;
	nb_Res_Hors = 0; for(k=0;k<nb_Aas;k++){nb_Res_Hors += tab_Cj_hors[k];}
	//printf("Coucou %i\n",nb_Res_Hors);
	for (j = 0; j<nb_Aas; j++){
	    //printf("Coucou %i %f %i %f\n",tab_Cj_hors[j],tab_bj[j],nb_Res_Hors,sqrt(N));
	    sum += (tab_Cj_hors[j]+tab_bj[j])/(nb_Res_Hors+sqrt(N));
	    //printf("Int Aa =%c : %f\n",Aa,(tab_Cj_hors[j]+tab_bj[j])/(nb_Res_Hors+sqrt(N)));
	    tab_Pj[j] = (tab_Cj_hors[j]+tab_bj[j])/(nb_Res_Hors+sqrt(N)); }
	    //printf("Somme des Pj : %f\n",sum);
}


void remplissage_Ax(char seq_exclue[], int nb_Aas, int taille_segment, float Qij[nb_Aas][taille_segment], float tab_Pj[nb_Aas], float tab_Ax[strlen(seq_exclue)-taille_segment],char lesAas[])
{	//printf("%s\n","");

	//printf("%s %c\n","coucou",seq_exclue[2]);
	int i, j, compt = 0, loc_Aa = 0, Aa = 0, res=0;

	for (i = 0; i<strlen(seq_exclue)-taille_segment+1; i++)
	{
	    for (j = 0; j<taille_segment; j++)
	    {
        loc_Aa = 0;
		res = i + j;
		Aa = seq_exclue[res];
		for (compt = 0; compt < nb_Aas; compt++){
		    if (Aa == lesAas[compt]) {
		        loc_Aa = compt;}}
	        tab_Ax[i] *= (fabsf(Qij[loc_Aa][j])/fabsf(tab_Pj[loc_Aa]));
	    }
	}
}
void remplissage_Ax_normalise(char seq_exclue[], int taille_segment,float tab_Ax[strlen(seq_exclue)-taille_segment]){
	int i = 0;
	float sum_Ax = 0.0;

	for (i = 0; i<strlen(seq_exclue)-taille_segment+1; i++){sum_Ax += tab_Ax[i];}
	//printf("Here %f\n",sum_Ax);
	for (i = 0; i<strlen(seq_exclue)-taille_segment+1; i++){tab_Ax[i] /= sum_Ax;}


}

int remplissage_tab_index(int *tab_index, int W, int index_seq_ex, char **Liste_Sequences, int *Liste_Seg)
{
	int nshift = W/2, i=0, j=0, left = 0, right =0;

	/*
	printf("Ind ex %i\n",Liste_Seg[index_seq_ex]);
	printf("Nshift : %i\n",nshift);
	printf("Borne gauche : %i\n",Liste_Seg[index_seq_ex]);
	printf("Borne droite : %li\n",strlen(Liste_Sequences[index_seq_ex])-Liste_Seg[index_seq_ex]-W);
	*/

	/* Si deplacement gauche et droite : OK alors */
	if ( (Liste_Seg[index_seq_ex] >= nshift) && (Liste_Seg[index_seq_ex] <= strlen(Liste_Sequences[index_seq_ex])-W-nshift) )
	{   /*printf("Probleme 1 : 2 bornes OK\n");*/
		for (i=0; i<= nshift; i++) {
	        if(i==0){ tab_index[i] = Liste_Seg[index_seq_ex];}
	        else { tab_index[i] = Liste_Seg[index_seq_ex]+i; tab_index[i+nshift] = Liste_Seg[index_seq_ex]-i;}}
	return nshift*2+1;}

	/* Si seulement deplacement gauche : OK */
	else if (Liste_Seg[index_seq_ex] >= nshift)
	{/*printf("Probleme 2 : borne droite\n");*/
	if (strlen(Liste_Sequences[index_seq_ex])-Liste_Seg[index_seq_ex]-W == 0) //Pas de deplacement possible
		{for (j = 0; j <= nshift; j++)
		{	if(j==0){tab_index[j] = Liste_Seg[index_seq_ex];}
			else {tab_index[j] = Liste_Seg[index_seq_ex] -j;}
		}return j;}
	else {
	/* right compte le nombre de deplacement possible vers la droite*/
	right = strlen(Liste_Sequences[index_seq_ex])-Liste_Seg[index_seq_ex]-W;
	for (j = 0; j <= nshift; j++) {
		if(j==0){tab_index[j] = Liste_Seg[index_seq_ex];}
		else { tab_index[j] = Liste_Seg[index_seq_ex] -j;}}
	for (i = 1; i <=right; i++) {
		tab_index[j] =  Liste_Seg[index_seq_ex] + i;
		j++;}
	}return j;}

	/* Si seulement deplacement droite : OK */
	else if (Liste_Seg[index_seq_ex] <= (strlen(Liste_Sequences[index_seq_ex])-W-nshift))
	{/*printf("Probleme 3 : borne gauche\n");*/
	if (Liste_Seg[index_seq_ex] == 0) //Pas de deplacement possible
		{for (j = 0; j <= nshift; j++)
		{
			if(j==0){tab_index[j] = Liste_Seg[index_seq_ex];}
			else {tab_index[j] = Liste_Seg[index_seq_ex] + j;}
		}return j;}
	else {
	/* left compte le nombre de deplacement possible vers la gauche*/
	left = Liste_Seg[index_seq_ex];
	for (j = 0; j <= nshift; j++) {
		if(j==0){tab_index[j] = Liste_Seg[index_seq_ex];}
		else { tab_index[j] = Liste_Seg[index_seq_ex] + j;}}
	for (i = 1; i <=left; i++) {
		tab_index[j] =  Liste_Seg[index_seq_ex] - i;
		j++;}
	}return j;}
return 0;}

void remplissage_tab_Ax (int n, int taille_segment, float **tab_Ax, char seq_exclue[],float Ax[strlen(seq_exclue)-taille_segment+1])
{ 	int i =0;
	tab_Ax[n] = malloc((strlen(seq_exclue)-taille_segment+1)*sizeof(float));
	tab_Ax[n] = memset(tab_Ax[n],0,sizeof(tab_Ax));
	for(i=0;i<strlen(seq_exclue)-taille_segment+1;i++){
	tab_Ax[n][i] = Ax[i];}
 }

void remplissage_tab_Ax_normalise(int n,char seq_exclue[], int taille_segment,float *tab_Ax[strlen(seq_exclue)-taille_segment]){
	int i = 0,j;
	float sum_Ax = 0.0;
	for (i = 0; i<n; i++){
	for (j = 0; j<strlen(seq_exclue)-taille_segment+1; j++){sum_Ax += tab_Ax[i][j];}}
	//printf("Here %f\n",sum_Ax);
	for (i = 0; i<n; i++){
	for (j = 0; j<strlen(seq_exclue)-taille_segment+1; j++){tab_Ax[i][j] /= sum_Ax;}}

	/*for (i = 0; i<n; i++){
	for (j = 0; j<strlen(seq_exclue)-taille_segment+1; j++){printf("i = %i j = %i %f\n",i,j,tab_Ax[i][j]);}}*/
}

void make_dir_file(int nbr_seq, int W, char **tab_sequences,int *vect_indice,char **titres)
{
	int i,j,ind;
	/* Creation des repertoires qui vont contenir le fichier les motifs convergents */
	mkdir("./Gibbs_/",0775);
	mkdir("./Gibbs_/motifs_convergents/",0775);

	/* Initialisation du fichier qui va contenir les sequences */
    FILE* fichier = NULL;
	/* Creation et ouverture du fichier qui va contenir les motifs */
    fichier = fopen("./Gibbs_/motifs_convergents/motifs.fst", "w");

	if (fichier != NULL){
	for(i=0;i<nbr_seq;i++){
		fputs(titres[i],fichier);
		ind = vect_indice[i];
	for(j=0;j<W;j++){
			fputc((tab_sequences[i]+ind)[j],fichier);}
	fprintf(fichier,"\n");}}
	else { printf("Probleme avec l'ouverture et/ou la creation du fichier contenant les motifs convergents : PAS DE SAUVEGARDE\n");}
    fclose(fichier);
}
