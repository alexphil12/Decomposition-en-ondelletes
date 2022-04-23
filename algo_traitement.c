
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "algo_traitement.h"
#include <fftw3.h>


long int min_fft(long int dim){
int k;
k=1;
while(k<dim){ //Cet algorithme calcule et renvoie l première puissance de 2 supérieur à dim
k*=2;}
return(k);}

int inv( unsigned int n, unsigned int max)
{
 unsigned int b=0;
 int u;
 while(max)      //Cet algorithme calcule le nombre obtenu par inversion de lecture binaire du nombre n avec comme pivot d'inversion max
 {               //Note:max est une puissance de 2 et n est inférieur à max
  if(n&1) b |= max ;
  max >>=1;
  n >>=1;
 }
 u=b;
 return u;
}

void inversion_tab_bit_ligne(int *tab_fft,long int dimx_ff,long int num_ligne, complexe *inter_fft){
long int k,u;
for(k=0;k<dimx_ff;k++){
  u=inv((unsigned int) k,(unsigned int)dimx_ff/2); //exécute l'inversion de lecture binnaire et change la place des coef[k] avec coef[inv(k)]
  inter_fft[u+dimx_ff*num_ligne].re=(double)tab_fft[k+dimx_ff*num_ligne];
  inter_fft[u+dimx_ff*num_ligne].im=0;}// pour chaque ligne
}

void inversion_tab_bit_colonne(complexe *inter_fft,long int dimy_ff,long int dimx_ff,long int num_colonne, complexe *ir_fft){
long int k,u;
for(k=0;k<dimy_ff;k++){
  u=inv((unsigned int)k,(unsigned int)dimy_ff/2);//idem avec les colonnes
  ir_fft[num_colonne+dimx_ff*u]=inter_fft[num_colonne+dimx_ff*k];}//idem que la fonction précédente pour les colonnes
}

void fft_1D_ligne(int* tab_fft,long int dimx_ff,long int num_ligne,complexe *inter_fft){
inversion_tab_bit_ligne(tab_fft,dimx_ff,num_ligne,inter_fft);
long int n,s,m; //applique une fft itérative sur la ligne num_ligne et écris le résultat dans inter fft
double pi;
pi=4*atan(1.0);
long int k,j;
complexe w,wm,t,u;
n=log2(dimx_ff);
for(s=1;s<n+1;s++){
   m=power(2,s);
   wm=expo_complexe(-2*pi/m);
   for(k=0;k<dimx_ff;k+=m){
       w.re=1;
       w.im=0;
       for(j=0;j<m/2;j++){
          t=Produit_complexe(w,inter_fft[dimx_ff*num_ligne+k+j+m/2]);
          u=inter_fft[dimx_ff*num_ligne+k+j];
          inter_fft[dimx_ff*num_ligne+k+j]=Somme_complexe(u,t);
          inter_fft[dimx_ff*num_ligne+k+j+m/2]=Somme_complexe(u,Produit_comp_re(t,-1));
          w=Produit_complexe(w,wm);}
}
}
}

void fft_1D_colonne(complexe *inter_fft,long int dimy_ff,long int dimx_ff,long int num_colonne,complexe *ir_fft){
inversion_tab_bit_colonne(inter_fft,dimy_ff,dimx_ff,num_colonne,ir_fft);
int n,s,m; //reprend le résultat de inter_fft et écris une fft itérative pour chaque colonne dans le tableau ir_fft
double pi;
pi=4*atan(1.0);
long int k,j;
complexe w,wm,t,u;
n=log2(dimy_ff);
for(s=1;s<n+1;s++){
   m=power(2,s);
   wm=expo_complexe(-2*pi/m);
   for(k=0;k<dimy_ff;k+=m){
       w.re=1;
       w.im=0;
       for(j=0;j<m/2;j++){
          t=Produit_complexe(w,ir_fft[dimx_ff*(k+j+m/2)+num_colonne]);
          u=ir_fft[dimx_ff*(k+j)+num_colonne];
          ir_fft[dimx_ff*(k+j)+num_colonne]=Somme_complexe(u,t);
          ir_fft[dimx_ff*(k+j+m/2)+num_colonne]=Somme_complexe(u,Produit_comp_re(t,-1));
          w=Produit_complexe(w,wm);}
}
}
}

void fft_1D_ligne_recur(int *tab_fft,long int dimx_ff,long int num_ligne,complexe *inter_fft){
complexe *ligne;
ligne=(complexe *)malloc(dimx_ff*sizeof(complexe));//applique une fft recusive su les ligne inter_fft
long int k;
for(k=0;k<dimx_ff;k++){
   ligne[k].re=tab_fft[num_ligne*dimx_ff+k];
   ligne[k].im=0;}

fft(ligne,dimx_ff);

for(k=0;k<dimx_ff;k++){
   inter_fft[num_ligne*dimx_ff+k].re=ligne[k].re;
   inter_fft[num_ligne*dimx_ff+k].im=ligne[k].im;
}
free(ligne);
}

void fft_1D_colonne_recur(complexe *inter_fft,long int dimy_ff,long int dimx_ff,long int num_colonne,complexe *ir_fft){
complexe *colonne;
colonne=(complexe *)malloc(dimy_ff*sizeof(complexe));//idem sur les colonnes
long int k;
for(k=0;k<dimy_ff;k++){
   colonne[k].re=inter_fft[k*dimx_ff+num_colonne].re;
   colonne[k].im=inter_fft[k*dimx_ff+num_colonne].im;}

fft(colonne,dimy_ff);

for(k=0;k<dimy_ff;k++){
   ir_fft[k*dimx_ff+num_colonne].re=colonne[k].re;
   ir_fft[k*dimx_ff+num_colonne].im=colonne[k].im;
}
free(colonne);
}






void Fourier_2D(int *tab_fft,long int dimx_ff,long int dimy_ff,complexe* ir_fft,complexe* inter_fft){
long int k;
for(k=0;k<dimy_ff;k++){
   fft_1D_ligne(tab_fft,dimx_ff,k,inter_fft);}//apllique les fonctions de fft 1D sur les lignes puis les colonnes
for(k=0;k<dimx_ff;k++){
   fft_1D_colonne(inter_fft,dimy_ff,dimx_ff,k,ir_fft);}
}



void Conv(double *A, double *B,long int lenA,long int lenB,int *C){
	int nconv;
	int i, j, i1;
	double tmp;

	//mémoire pour le tableau résultat
	nconv = lenA+lenB-1;


	//La convolution
	for (i=0; i<nconv; i++)
	{
		i1 = i;
		tmp = 0.0;
		for (j=0; j<lenB; j++)
		{
			if(i1>=0 && i1<lenA)
				tmp = tmp + (A[i1]*B[j]);

			i1 = i1-1;
			C[i] = tmp;
		}
	}





}
void Extract_ligne(int*tab,long int dimx,long int nu_ligne,double *ligne,int u){
long int k;
for(k=0;k<dimx/power(2,u);k++){//extrait la ligne necessaire pour l'analyse ou la reconstruction au niveau u
   ligne[k]=tab[nu_ligne*dimx+k];}}

void Extract_colonne(int*tab,long int dimy,long int nu_colonne,long int dimx,double *colonne,int u){
long int k;
for(k=0;k<dimy/power(2,u);k++){//même principe pour les colonnes
   colonne[k]=tab[k*dimx+nu_colonne];}
}
complexe Somme_complexe(complexe z1,complexe z2){//ensemble d'opération élémentaire pour le nouveau type complexe
complexe z3;
z3.re=z1.re+z2.re;
z3.im=z1.im+z2.im;
return(z3);}

complexe Produit_complexe(complexe z1,complexe z2){
complexe z3;
z3.re=(z1.re*z2.re)-(z1.im-z2.im);
z3.im=(z1.re*z2.im)+(z1.im*z2.re);
return(z3);}

complexe Produit_comp_re(complexe z1,double r1){
complexe z3;
z3.re=z1.re*r1;
z3.im=z1.im*r1;
return(z3);}

double val_abs(complexe z1){
double abs;
abs=sqrt(z1.re*z1.re+z1.im*z1.im);
return(abs);}

complexe expo_complexe(double theta){
complexe z1;
z1.re=cos(theta);
z1.im=sin(theta);
return(z1);}

double argument(complexe z1){
double theta;
if(z1.im==0 && z1.re<0){
theta=4*atan(1.0);}
else{
theta=2*atan(z1.im/(z1.re+val_abs(z1)));
}
return(theta);}

int power(int a,int b){
int k,resultat;
resultat=a;
if(b==0){
   resultat=1;}
else{
for(k=1;k<b;k++){
resultat *=a;}}
return(resultat);}

void max_tab(int* ir,long int dimx,long int dimy){//normalisation pour mieux voir en paçant le max à 255
long int max,k,l;
max=0;
for(k=0;k<dimx;k++){
   for(l=0;l<dimy;l++){
      if(ir[dimx*l+k]>=max){
         max=ir[dimx*l+k];}
   }
}
for(k=0;k<dimx;k++){
   for(l=0;l<dimy;l++){
      ir[dimx*l+k]=ir[dimx*l+k]*(255/max);}
}
}

void fft_shift_ligne(int* ir, long int dimx_ff,long int num_ligne){
long int k;
int* ligne_ff;
ligne_ff=(int *)malloc(dimx_ff*sizeof(int));
for(k=0;k<dimx_ff/2;k++){
   ligne_ff[k]=ir[num_ligne*dimx_ff+(dimx_ff/2)+k];//remet les fréquences nulles au centre en partenariat avec fft_shift_colonne
   ligne_ff[(dimx_ff/2)+k]=ir[num_ligne*dimx_ff+k];
}
for(k=0;k<dimx_ff;k++){
   ir[num_ligne*dimx_ff+k]=ligne_ff[k];
}
}


void fft_shift_colonne(int* ir, long int dimx_ff,long int dimy_ff,long int num_colonne){
long int k;
int* colonne_ff;
colonne_ff=(int *)malloc(dimy_ff*sizeof(int));
for(k=0;k<dimy_ff/2;k++){
   colonne_ff[k]=ir[(k+(dimy_ff/2))*dimx_ff+num_colonne];
   colonne_ff[k+(dimy_ff/2)]=ir[k*dimx_ff+num_colonne];
}
for(k=0;k<dimy_ff;k++){
   ir[k*dimx_ff+num_colonne]=colonne_ff[k];
}
}



void fft(complexe lig_col[], long int dim)
{
	complexe out[dim];
	for (int i = 0; i < dim; i++){
           out[i].re=lig_col[i].re;
           out[i].im=lig_col[i].im;}//calcul à propremenent parler la fft recursive.
        _fft(lig_col, out, dim, 1);
}

void _fft(complexe buf[], complexe out[], long int dim, int step){
double pi;
pi=4*atan(1.0);
	if(step<dim){
		_fft(out,buf,dim,step*2);
		_fft(out+step,buf+step,dim,step*2);

		for (int i = 0; i < dim; i += 2 * step) {
			complexe t=Produit_complexe(expo_complexe(-1*pi*i/dim),out[i+step]);
			buf[i/2]=Somme_complexe(out[i],t);
			buf[(i+dim)/2]=Somme_complexe(out[i],Produit_comp_re(t,-1));
		}
	}
}


void rotation (int *tab,long int dimx_ff,long int dimy_ff,double theta){
int *out;
out=(int *)malloc(dimx_ff*dimy_ff*sizeof(int)); permet d'effectuer une rotation d'ange theta pour tout image
int x,y;
int tempx,tempy;
for(x=0;x<dimx_ff;x++){
   for(y=0;y<dimy_ff;y++){
      tempx=floor((cos(theta)*(x-dimx_ff/2))+(-1*sin(theta)*(y-dimy_ff/2))+dimx_ff/2);
      tempy=floor((sin(theta)*(x-dimx_ff/2))+(cos(theta)*(y-dimy_ff/2))+dimy_ff/2);
      if(tempx<dimx_ff && tempx>0 && tempy<dimy_ff && tempy>0){
         out[tempy*dimx_ff+tempx]=tab[y*dimx_ff+x];}
      }}
for(x=0;x<dimx_ff;x++){
   for(y=0;y<dimy_ff;y++){
      tab[y*dimx_ff+x]=out[y*dimx_ff+x];}}


}

void DFT_2D(int *tab_fft,long int dimx_ff,long int dimy_ff,complexe* out){
long int x,y,u,v;
double pi;
pi=4*atan(1.0);
complexe *in;
in=(complexe *)malloc(dimx_ff*dimy_ff*sizeof(complexe));// transformée Discrète naïve de Fourier 
for(x=0;x<dimx_ff;x++){
   for(y=0;y<dimy_ff;y++){
      in[y*dimx_ff+x].re=tab_fft[y*dimx_ff+x];
      in[y*dimx_ff+x].im=0;
      out[y*dimx_ff+x].re=0;
      out[y*dimx_ff+x].im=0;
}}
for(u=0;u<dimx_ff;u++){
   printf("%d\n",u);
   for(v=0;v<dimy_ff;v++){
      for(x=0;x<dimx_ff;x++){
         for(y=0;y<dimy_ff;y++){
            out[v*dimx_ff+u]=Somme_complexe(out[v*dimx_ff+u],Produit_complexe(in[y*dimx_ff+x],expo_complexe(-2*pi*((u*x/dimx_ff)+(v*y/dimy_ff)))));
      }}
   out[v*dimx_ff+u]=Produit_comp_re(out[v*dimx_ff+u],(1/(dimx_ff*dimy_ff)));
}}
}




