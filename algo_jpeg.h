#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct arbre_code arbre_code;
struct arbre_code
{
   int nbr;
   int freq;
   struct arbre_code *branche0;
   struct arbre_code *branche1;
};

typedef struct dico_huff dico_huff;
struct dico_huff
{
   int nbr;
   int freq;
};

typedef struct dico_encodage dico_encodage;
struct dico_encodage
{
   int symb;
   char *code;
   short taille;
};

long int min_8(long int dim);
void discret_cos_8(long int i,long int j,double* inter,long int dimx_8);
void quantif(long int i, long int j,double* inter,double fq,long int dimx_8,int *quantificateur);
void zigzag(double* inter,int i,int j,long int dimx_8,int* code);
void DCT_8(double bloc_in[64], double bloc_out[64]);
int rlc(int *code,int *Vlrc,long int dimx_8,long int dimy_8);
int nb_dif_VLRC(int* Vlrc,dico_huff *verif,int taille_code);
void creer_arbre_huff(arbre_code *tab_arbre[],int nbr_dif);
int creer_code_huff(arbre_code *noeud,char *code,int nbr_dif,unsigned short pos,dico_encodage *dico_huff,int place, arbre_code *origine);
int decode_huff(dico_encodage *dico,char *sequ,int *Vlrc,int nbr_dif);
void Vlrc_zigzag(int taille_code,int *Vlrc,int *zigzag);
void zigzag_inv(int *zigzag, long int dimx, long int dimy,int* iu,long int i,long int j);
void quantification_inv(double fq,int *iu,long int i,long int j,int *quantificateur,long int dimx_8);
void discret_cos_8_inv(long int i,long int j,long int dimx_8,int* iu);
void DCT_8_inv(double block_in[64], int block_out[64]);