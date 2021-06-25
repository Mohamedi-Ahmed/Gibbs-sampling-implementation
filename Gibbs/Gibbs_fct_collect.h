#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Gibbs_fct_calculus.h"

#define TAILLE_MAX_SEGMENT 1000
#define TAILLE_MAX_SEQUENCE 500

int read_fasta(char **tab, char *Nom_Fichier, char **titres)
{
    // PARAMETRES
    FILE* fichier = NULL;
    int i = 0;
    char chaine[TAILLE_MAX_SEQUENCE] = "";
    const char *retour = "\n";

    //Lecture du fichier
    fichier = fopen(Nom_Fichier, "r");
    if (fichier != NULL)
    {
        while (fgets(chaine,TAILLE_MAX_SEQUENCE, fichier) != NULL) // Lecture chaine a chaine
        {
	    if (chaine[0] != '>' ) //62 == '>' en ASCII
	       {tab[i] = malloc(strlen(chaine));
		memset(tab[i],'\0',sizeof(tab));
		strncpy(tab[i],chaine,strlen(chaine)-1);
		/* strlen(chaine) renvoie la taille de la chaine en prenant en compte \n mais pas \0
		pour ne pas prendre en compte ces deux caract. speciaux, on rajoute le -1
		On copie la chaine de sa taille - le caractere special \n */
	        i++;}
		else{
			titres[i] = malloc(strlen(chaine));
			memset(titres[i],'\0',sizeof(titres));
			strncpy(titres[i],chaine,strlen(chaine)-2);
            strcat(titres[i],retour);
        }
        }

        fclose(fichier);
    }
	return i;
}


