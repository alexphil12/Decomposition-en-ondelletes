#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct complexe complexe;
struct complexe
{
   double re;
   double im;
};

int inv( unsigned int n, unsigned int max);
long int min_fft(long int dim);
void Fourier_2D(int *tab_fft,long int dimx_ff,long int dimy_ff,complexe* ir_fft,complexe* inter_fft);
void inversion_tab_bit_ligne(int *tab_fft,long int dimx_ff,long int num_ligne, complexe *inter_fft);
void inversion_tab_bit_colonne(complexe *inter_fft,long int dimy_ff,long int dimx_ff,long int num_colonne, complexe *ir_fft);
void fft_1D_colonne(complexe *inter_fft,long int dimy_ff,long int dimx_ff,long int num_colonne,complexe *ir_fft);
void fft_1D_ligne(int* tab_fft,long int dimx_ff,long int num_ligne,complexe *inter_fft);
void Conv(double *A, double *B,long int lenA,long int lenB,int *C);
void Extract_ligne(int*tab,long int dimx,long int nu_ligne,double *ligne,int u);
void Extract_colonne(int*tab,long int dimy,long int nu_colonne,long int dimx,double *colonne,int u);
complexe Produit_comp_re(complexe z1,double r1);
complexe Somme_complexe(complexe z1,complexe z2);
complexe Produit_complexe(complexe z1,complexe z2);
double val_abs(complexe z1);
complexe expo_complexe(double theta);
double argument(complexe z1);
int power(int a,int b);
void max_tab(int* ir,long int dimx,long int dimy);
void fft_shift_ligne(int* ir, long int dimx_ff,long int num_ligne);
void fft_shift_colonne(int* ir, long int dimx_ff,long int dimy_ff,long int num_colonne);
void shift(int *ir,long int dimx_ff,long int dimy_ff);
void fft_1D_ligne_recur(int *tab_fft,long int dimx_ff,long int num_ligne,complexe *inter_fft);
void fft_1D_colonne_recur(complexe *inter_fft,long int dimy_ff,long int dimx_ff,long int num_colonne,complexe *ir_fft);
void fft(complexe lig_col[], long int dim);
void _fft(complexe buf[], complexe out[], long int dim, int step);
void rotation (int *tab,long int dimx_ff,long int dimy_ff,double theta);