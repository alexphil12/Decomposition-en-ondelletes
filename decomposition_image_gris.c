#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "algo_traitement.c"
FILE *fio,*fir;


int main()
{
unsigned char chi;

int *io,*ir,*e_li,*w_li,*e_col,*w_col,*iu;
int *inter2,*inter;
long int tour,nbg;
register long int i,j;
long int dimx,dimy,taille;
long int dimx_ff,dimy_ff;
double angle;
char entree[35],sortie[35],chaine[10],Traitement[6];


printf("donnez le fichier entree sans .pgm\n");
scanf("%s",entree);
strcat(entree,".pgm");

printf("donnez le fichier sortie sans .pgm\n");
scanf("%s",sortie);
strcat(sortie,".pgm");

printf("Donnez le traitement souhaite parmis les suivants:\n 'Har' pour Haar,'Dbn' pour Daubechies 4,'5/3' pour biortho 5/3, '9/7' pour biortho 9/7,'Fou' pour Fourier\n");
scanf("%s",Traitement);

if((strcmp(Traitement,"Fou")!=0)){
printf("donner le nombre de tour (pas superieur a log2(min(dimx;dimy))\n");
scanf("%d",&tour);
}
else{
printf("Entrez l'angle de rotation souhaitee en radian\n");
scanf("%lf",&angle);}

fio=fopen(entree,"rb");
fir=fopen(sortie,"wb+");

if(fio !=NULL){printf("ouverture du fichier in\n");}
if(fir !=NULL){printf("ouverture du fichier out\n");}

//if(im!=0){printf("allocation de im\n")};

//lecture entete pgm
fscanf(fio,"%s\n",chaine);
fprintf(fir,"%s\n",chaine);

i=0;
while(i<1)
{chi=fgetc(fio);
// if((int)chi==6)  fputc((unsigned char)5,fr);
 fputc(chi,fir);
 if(chi=='\n') i++;
}

fscanf(fio,"%ld %ld\n",&dimx,&dimy);
dimx_ff=min_fft(dimx);
dimy_ff=min_fft(dimy);

if(strcmp(Traitement,"Fou")==0){
   fprintf(fir,"%ld %ld\n",dimx_ff,dimy_ff);}
else{
   fprintf(fir,"%ld %ld\n",dimx,dimy);}



fscanf(fio,"%d\n",&nbg);
fprintf(fir,"%d\n",nbg);

printf("dimx=%d dimy=%d nbg=%d\n",dimx,dimy,nbg);

io=(int *)malloc(dimx*dimy*sizeof(int));

iu=(int *)malloc(dimx*dimy*sizeof(int));

if(strcmp(Traitement,"Fou")==0){
   ir=(int *)malloc(dimx_ff*dimy_ff*sizeof(int));}
else{
   ir=(int *)malloc(dimx*dimy*sizeof(int));}

inter=(int *)malloc(dimx*dimy*sizeof(int));

inter2=(int *)malloc(dimx*dimy*sizeof(int));

if(ir==NULL) printf("allocation ir impossible\n");

taille=dimx*dimy;
printf("taille=%ld\n",taille);

//lecture des pixels
for(i=0;i<dimy;i++)
 for(j=0;j<dimx;j++)
 {
  chi=(unsigned char)fgetc(fio);
  io[i*dimx+j]=(int)chi;
 }
//Traitement

if(strcmp(Traitement,"Fou")==0){
   int *tab_fft,*ligne_ff,*colonne_ff;
   complexe *inter_fft,*ir_fft;
   tab_fft=(int *)malloc(dimy_ff*dimx_ff*sizeof(int));
   ir_fft=(complexe *)malloc(dimy_ff*dimx_ff*sizeof(complexe));
   inter_fft=(complexe *)malloc(dimy_ff*dimx_ff*sizeof(complexe));
   ligne_ff=(int*)malloc(dimx_ff*sizeof(int));
   colonne_ff=(int*)malloc(dimy_ff*sizeof(double));
   long int k,l,u,v,p,m;
   for(k=0;k<dimx_ff;k++){
      for(l=0;l<dimy_ff;l++){
         tab_fft[dimx_ff*l+k]=0;}}
   for(u=0;u<dimx;u++){
      for(v=0;v<dimy;v++){
         tab_fft[dimx_ff*v+u]=io[dimx*v+u]*power(-1,u+v);}}
   DFT_2D(tab_fft,dimx_ff,dimy_ff,ir_fft);
   for(p=0;p<dimx_ff;p++){
      for(m=0;m<dimy_ff;m++){
         ir[dimx_ff*m+p]=(int)val_abs(ir_fft[dimx_ff*m+p]);}}
   //shift(ir,dimx_ff,dimy_ff);
   for(p=0;p<dimx_ff;p++){
      fft_shift_colonne(ir,dimx_ff,dimy_ff,p);
   }
   for(p=0;p<dimy_ff;p++){
      fft_shift_ligne(ir,dimx_ff,p);
   }
   rotation (ir,dimx_ff,dimy_ff,angle);
   free(colonne_ff);
   free(ligne_ff);
   free(ir_fft);
   free(tab_fft);
   free(inter_fft);
   goto fin_fft;
   }

if(strcmp(Traitement,"Har")==0){
   int k;
   long int dimx_t=dimx,dimy_t=dimy;
   for(i=0;i<dimy;i++){
         for(j=0;j<dimx;j++){
            iu[i*dimx+j]=io[i*dimx+j];}}
   for(k=0;k<tour;k++){
      for(i=0;i<dimy_t/2;i++){
         for(j=0;j<dimx_t;j++){
	   inter[i*dimx+j]=(iu[(2*i)*dimx+j]+iu[2*(i+1)*dimx+j])/2;}}
      for(i=0;i<dimy_t/2;i++){
         for(j=0;j<dimx_t;j++){
	   inter[(i+dimy_t/2)*dimx+j]=(iu[(2*i)*dimx+j]-iu[2*(i+1)*dimx+j])/2;}}
      for(i=0;i<dimy_t;i++){
         for(j=0;j<dimx_t/2;j++){
	   ir[i*dimx+j]=(inter[i*dimx+2*j]+inter[i*dimx+2*j+1])/2;}}
      for(i=0;i<dimy_t;i++){
         for(j=0;j<dimx_t/2;j++){
	   ir[i*dimx+(j+dimx_t/2)]=(inter[i*dimx+2*j]-inter[i*dimx+2*j+1])/2;}}
      for(i=0;i<dimy;i++){
         for(j=0;j<dimx;j++){
            iu[i*dimx+j]=ir[i*dimx+j];}}
      dimx_t=dimx_t/2;
      dimy_t=dimy_t/2;
}
free(inter);
goto fin_norm;
}

if(strcmp(Traitement,"Dbn")==0){
int k;
long int dimx_t=dimx,dimy_t=dimy;
for(i=0;i<dimy;i++){
      for(j=0;j<dimx;j++){
         iu[i*dimx+j]=io[i*dimx+j];}}

double rg_daub4[4]={-0.4830,0.8365,-0.2241,-0.1294};
double rh_daub4[4]={-0.1294,0.2241,0.8365,0.4830};

for(k=0;k<tour;k++){
e_li=(int *)malloc((dimx_t+3)*sizeof(int));
w_li=(int *)malloc((dimx_t+3)*sizeof(int));
e_col=(int *)malloc((dimy_t+3)*sizeof(int));
w_col=(int *)malloc((dimy_t+3)*sizeof(int));
double ligne[dimx_t];
double colonne[dimy_t];
  for(i=0;i<dimy_t;i++){
   Extract_ligne(iu,dimx,i,ligne,k);
   Conv(rh_daub4,ligne,4,dimx_t,e_li);
   Conv(rg_daub4,ligne,4,dimx_t,w_li);
   for(j=0;j<dimx_t;j++){
      if(j<dimx_t/2){
         inter2[i*dimx+j]=e_li[j*2];}
      if(j>=dimx_t/2){
         inter2[i*dimx+j]=w_li[j*2-dimx_t];}
         }
}

for(j=0;j<dimx_t;j++){
   Extract_colonne(inter2,dimy,j,dimx,colonne,k);
   Conv(rh_daub4,colonne,4,dimy_t,e_col);
   Conv(rg_daub4,colonne,4,dimy_t,w_col);
   for(i=0;i<dimy_t;i++){
      if(i<dimy_t/2){
         ir[i*dimx+j]=e_col[i*2];}
      if(i>=dimy_t/2){
         ir[i*dimx+j]=w_col[i*2-dimy_t];}
         }
}
for(i=0;i<dimy;i++){
         for(j=0;j<dimx;j++){
            iu[i*dimx+j]=ir[i*dimx+j];}}
dimx_t=dimx_t/2;
dimy_t=dimy_t/2;
free(e_li);
free(w_li);
free(e_col);
free(w_col);
}
goto fin_norm;
}

if(strcmp(Traitement,"9/7")==0){
int k;
long int dimx_t=dimx,dimy_t=dimy;
for(i=0;i<dimy;i++){
      for(j=0;j<dimx;j++){
         iu[i*dimx+j]=io[i*dimx+j];}}
double rg_9_7[7]={0.0913,-0.0575,-0.5913,1.1151,-0.5913,-0.0575,0.0913};
double rh_9_7[9]={0.0267,-0.0169,-0.0782,0.2669,0.6029,0.2669,-0.0782,-0.0169,0.0267};
for(k=0;k<tour;k++){
e_li=(int *)malloc((dimx_t+8)*sizeof(int));
w_li=(int *)malloc((dimx_t+6)*sizeof(int));
e_col=(int *)malloc((dimy_t+8)*sizeof(int));
w_col=(int *)malloc((dimy_t+6)*sizeof(int));
double ligne[dimx_t];
double colonne[dimy_t];
  for(i=0;i<dimy_t;i++){
   Extract_ligne(iu,dimx,i,ligne,k);
   Conv(rh_9_7,ligne,9,dimx_t,e_li);
   Conv(rg_9_7,ligne,7,dimx_t,w_li);
   for(j=0;j<dimx_t;j++){
      if(j<dimx_t/2){
         inter2[i*dimx+j]=e_li[j*2];}
      if(j>=dimx_t/2){
         inter2[i*dimx+j]=w_li[j*2-dimx_t];}
         }
}

for(j=0;j<dimx_t;j++){
   Extract_colonne(inter2,dimy,j,dimx,colonne,k);
   Conv(rh_9_7,colonne,9,dimy_t,e_col);
   Conv(rg_9_7,colonne,7,dimy_t,w_col);
   for(i=0;i<dimy_t;i++){
      if(i<dimy_t/2){
         ir[i*dimx+j]=e_col[i*2];}
      if(i>=dimy_t/2){
         ir[i*dimx+j]=w_col[i*2-dimy_t];}
         }
}
for(i=0;i<dimy;i++){
         for(j=0;j<dimx;j++){
            iu[i*dimx+j]=ir[i*dimx+j];}}
dimx_t=dimx_t/2;
dimy_t=dimy_t/2;
free(e_li);
free(w_li);
free(e_col);
free(w_col);
}
goto fin_norm;
}

if(strcmp(Traitement,"5/3")==0){
int k;
long int dimx_t=dimx,dimy_t=dimy;
for(i=0;i<dimy;i++){
      for(j=0;j<dimx;j++){
         iu[i*dimx+j]=io[i*dimx+j];}}
double rg_5_3[3]={-0.5,1,-0.5};
double rh_5_3[5]={-0.1250,0.25,0.75,0.25,-0.1250};
for(k=0;k<tour;k++){
e_li=(int *)malloc((dimx_t+4)*sizeof(int));
w_li=(int *)malloc((dimx_t+2)*sizeof(int));
e_col=(int *)malloc((dimy_t+4)*sizeof(int));
w_col=(int *)malloc((dimy_t+2)*sizeof(int));
double ligne[dimx_t];
double colonne[dimy_t];
for(i=0;i<dimy_t;i++){
   Extract_ligne(iu,dimx,i,ligne,k);
   Conv(rh_5_3,ligne,5,dimx_t,e_li);
   Conv(rg_5_3,ligne,3,dimx_t,w_li);
   for(j=0;j<dimx_t;j++){
      if(j<dimx_t/2){
         inter2[i*dimx+j]=e_li[j*2];}
      if(j>=dimx_t/2){
         inter2[i*dimx+j]=w_li[j*2-dimx_t];}
         }
}

for(j=0;j<dimx_t;j++){
   Extract_colonne(inter2,dimy,j,dimx,colonne,k);
   Conv(rh_5_3,colonne,5,dimy_t,e_col);
   Conv(rg_5_3,colonne,3,dimy_t,w_col);
   for(i=0;i<dimy_t;i++){
      if(i<dimy_t/2){
         ir[i*dimx+j]=e_col[i*2];}
      if(i>=dimy_t/2){
         ir[i*dimx+j]=w_col[i*2-dimy_t];}
         }
}
for(i=0;i<dimy;i++){
         for(j=0;j<dimx;j++){
            iu[i*dimx+j]=ir[i*dimx+j];}}
dimx_t=dimx_t/2;
dimy_t=dimy_t/2;
free(e_li);
free(w_li);
free(e_col);
free(w_col);
}
goto fin_norm;
}
if(strcmp(Traitement,"Fou")==0){
fin_fft:
for(i=0;i<dimy_ff;i++){
  for(j=0;j<dimx_ff;j++){
     fputc((unsigned char)ir[i*dimx_ff+j],fir);}}}
else{
fin_norm:
   for(i=0;i<dimy;i++){
     for(j=0;j<dimx;j++){
        fputc((unsigned char)ir[i*dimx+j],fir);}}
}
free(io); //...

free(ir);

fclose(fio);fclose(fir);
return 0;
}
