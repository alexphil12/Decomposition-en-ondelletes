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
int *inter;
int *inter2;
long int nbg;
register long int i,j;
long int dimx,dimy,taille;
long int dimx_ff,dimy_ff;
int tour,niveau;


char entree[35],sortie[35],chaine[10],Traitement[6]; 

printf("donner le fichier entree sans .pgm\n");
scanf("%s",entree);
strcat(entree,".pgm");

printf("donner le fichier sortie sans .pgm\n");
scanf("%s",sortie);
strcat(sortie,".pgm");

printf("Donner le traitement souhaite parmis les suivants:\n 'Har' pour Haar,'Dbn' pour Daubechies 4,'5/3' pour biortho 5/3, '9/7' pour biortho 9/7,'Fou' pour Fourrier\n");
scanf("%s",Traitement);

if(strcmp(Traitement,"Fou")!=0){
   printf("Donner le niveau d'analyse effectue sur l'image a reconstruire\n");
   scanf("%d",&niveau);

   printf("Donner le nombre de remonter souhaite\n");
   scanf("%d",&tour);
}

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


fscanf(fio,"%ld\n",&nbg);
fprintf(fir,"%ld\n",nbg);

printf("dimx=%ld dimy=%ld nbg=%ld\n",dimx,dimy,nbg);

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
   Fourier_2D(tab_fft,dimx_ff,dimy_ff,ir_fft,inter_fft);
   for(p=0;p<dimx_ff;p++){
      for(m=0;m<dimy_ff;m++){
         ir[dimx_ff*m+p]=(int)val_abs(ir_fft[dimx_ff*m+p]);}}
   shift(ir,dimx_ff,dimy_ff);
   for(p=0;p<dimx_ff;p++){
      fft_shift_colonne(ir,dimx_ff,dimy_ff,p);
   }
   for(p=0;p<dimy_ff;p++){
      fft_shift_ligne(ir,dimx_ff,p);
   }
   free(colonne_ff);
   free(ligne_ff);
   free(ir_fft);
   free(tab_fft);
   free(inter_fft);
   goto fin_fft;
   }

if(strcmp(Traitement,"Har")==0){
long int k;
int u;
niveau=niveau-1;
long int dimx_t=dimx/(power(2,niveau)),dimy_t=dimy/power(2,niveau);
for(i=0;i<dimy;i++){
        for(j=0;j<dimx;j++){
           iu[i*dimx+j]=io[i*dimx+j];}}
for(u=tour;u>=1;u--){
e_li=(int *)malloc((2*dimx_t+2)*sizeof(int));
w_li=(int *)malloc((2*dimx_t+2)*sizeof(int));
e_col=(int *)malloc((2*dimy_t+2)*sizeof(int));
w_col=(int *)malloc((2*dimy_t+2)*sizeof(int));
double g_Har[3]={0.5,-0.5,0};
double h_Har[3]={0.5,0.5,0};
double ligne[dimx_t];
double colonne[dimy_t];
double tab_ligne[2*dimx_t];
double tab_colonne[2*dimy_t];

for(i=0;i<dimy_t;i++){
   Extract_ligne(iu,dimx,i,ligne,niveau);
   for(k=0;k<dimx_t;k++){
      tab_ligne[2*k]=0;
      tab_ligne[2*k+1]=ligne[k];}
   Conv(h_Har,tab_ligne,3,2*dimx_t,e_li);
   Conv(g_Har,tab_ligne,3,2*dimx_t,w_li);
   for(j=0;j<dimx_t;j++){
   inter2[i*dimx+j]=e_li[j+2]+w_li[j+2];
}
}
for(i=0;i<dimy;i++){
   for(j=0;j<dimx;j++){
      if((i>dimy_t)||(j>dimx_t)){
        inter2[i*dimx+j]=io[i*dimx+j];}}}

for(j=0;j<dimx_t;j++){
   Extract_colonne(inter2,dimy,j,dimx,colonne,niveau);
   for(k=0;k<dimy_t;k++){
      tab_colonne[2*k]=0;
      tab_colonne[2*k+1]=colonne[k];}
   Conv(h_Har,tab_colonne,3,2*dimy_t,e_col);
   Conv(g_Har,tab_colonne,3,2*dimy_t,w_col);
   for(i=0;i<dimy_t;i++){
   ir[i*dimx+j]=e_col[i+2]+w_col[i+2];
}
}

for(i=0;i<dimy;i++){
   for(j=0;j<dimx;j++){
      if((i>dimy_t)||(j>dimx_t)){
        ir[i*dimx+j]=io[i*dimx+j];}}}
for(i=0;i<dimy;i++){
   for(j=0;j<dimx;j++){
        iu[i*dimx+j]=ir[i*dimx+j];}}
niveau=niveau-1;
dimx_t=dimx_t*2;
dimy_t=dimy_t*2;
free(e_li);
free(w_li);
free(e_col);
free(w_col);
} 

goto fin_norm;
}

if(strcmp(Traitement,"Dbn")==0){
long int k;
int u;
niveau=niveau-1;
long int dimx_t=dimx/(power(2,niveau)),dimy_t=dimy/power(2,niveau);
for(i=0;i<dimy;i++){
        for(j=0;j<dimx;j++){
           iu[i*dimx+j]=io[i*dimx+j];}}
for(u=tour;u>=1;u--){
e_li=(int *)malloc((2*dimx_t+4)*sizeof(int));
w_li=(int *)malloc((2*dimx_t+4)*sizeof(int));
e_col=(int *)malloc((2*dimy_t+4)*sizeof(int));
w_col=(int *)malloc((2*dimy_t+4)*sizeof(int));
double g_daub4[5]={-0.1294,-0.2241,0.8365,-0.4830,0};
double h_daub4[5]={-0.4830,0.8365,0.2241,-0.1294,0};
double ligne[dimx_t];
double colonne[dimy_t];
double tab_ligne[2*dimx_t];
double tab_colonne[2*dimy_t];

for(i=0;i<dimy_t;i++){
   Extract_ligne(iu,dimx,i,ligne,u);
   for(k=0;k<dimx_t;k++){
      tab_ligne[2*k]=0;
      tab_ligne[2*k+1]=ligne[k];}
   Conv(h_daub4,tab_ligne,5,2*dimx_t,e_li);
   Conv(g_daub4,tab_ligne,5,2*dimx_t,w_li);
   for(j=0;j<dimx_t;j++){
   inter2[i*dimx+j]=e_li[j]+w_li[j];
}
}

for(j=0;j<dimx_t;j++){
   Extract_colonne(inter2,dimy,j,dimx,colonne,u);
   for(k=0;k<dimy_t;k++){
      tab_colonne[2*k]=0;
      tab_colonne[2*k+1]=colonne[k];}
   Conv(h_daub4,tab_colonne,5,2*dimy_t,e_col);
   Conv(g_daub4,tab_colonne,5,2*dimy_t,w_col);
   for(i=0;i<dimy_t;i++){
   ir[i*dimx+j]=e_col[i]+w_col[i];
}
     
}
for(i=0;i<dimy;i++){
   for(j=0;j<dimx;j++){
      if((i>dimy_t)||(j>dimx_t)){
        ir[i*dimx+j]=io[i*dimx+j];}}}
for(i=0;i<dimy;i++){
   for(j=0;j<dimx;j++){
        iu[i*dimx+j]=ir[i*dimx+j];}}
niveau=niveau-1;
dimx_t=dimx_t*2;
dimy_t=dimy_t*2;
free(e_li);
free(w_li);
free(e_col);
free(w_col);
}
goto fin_norm;
}

if(strcmp(Traitement,"9/7")==0){ 
long int k;
int u;
niveau=niveau-1;
long int dimx_t=dimx/(power(2,niveau)),dimy_t=dimy/power(2,niveau);
for(i=0;i<dimy;i++){
        for(j=0;j<dimx;j++){
           iu[i*dimx+j]=io[i*dimx+j];}}
for(u=tour;u>=1;u--){
e_li=(int *)malloc((2*dimx_t+6)*sizeof(int));
w_li=(int *)malloc((2*dimx_t+8)*sizeof(int));
e_col=(int *)malloc((2*dimy_t+6)*sizeof(int));
w_col=(int *)malloc((2*dimy_t+8)*sizeof(int));
double g_9_7[9]={0.267,0.0169,-0.0782,-0.2669,0.6029,-0.2669,-0.0782,0.0169,0.0267};
double h_9_7[7]={-0.0913,-0.0575,0.5913,1.1151,0.5913,-0.0575,-0.0913};
double ligne[dimx_t];
double colonne[dimy_t];
double tab_ligne[2*dimx_t];
double tab_colonne[2*dimy_t];

for(i=0;i<dimy_t;i++){
   Extract_ligne(iu,dimx,i,ligne,u);
   for(k=0;k<dimx_t;k++){
      tab_ligne[2*k]=0;
      tab_ligne[2*k+1]=ligne[k];}
   Conv(h_9_7,tab_ligne,7,2*dimx_t,e_li);
   Conv(g_9_7,tab_ligne,9,2*dimx_t,w_li);
   for(j=0;j<dimx_t;j++){
   inter2[i*dimx+j]=e_li[j]+w_li[j];
}
}

for(j=0;j<dimx_t;j++){
   Extract_colonne(inter2,dimy,j,dimx,colonne,u);
   for(k=0;k<dimy_t;k++){
      tab_colonne[2*k]=0;
      tab_colonne[2*k+1]=colonne[k];}
   Conv(h_9_7,tab_colonne,7,2*dimy_t,e_col);
   Conv(g_9_7,tab_colonne,9,2*dimy_t,w_col);
   for(i=0;i<dimy_t;i++){
   ir[i*dimx+j]=e_col[i]+w_col[i];
}
     
}
for(i=0;i<dimy;i++){
   for(j=0;j<dimx;j++){
      if((i>dimy_t)||(j>dimx_t)){
        ir[i*dimx+j]=io[i*dimx+j];}}}
for(i=0;i<dimy;i++){
   for(j=0;j<dimx;j++){
        iu[i*dimx+j]=ir[i*dimx+j];}}
niveau=niveau-1;
dimx_t=dimx_t*2;
dimy_t=dimy_t*2;
free(e_li);
free(w_li);
free(e_col);
free(w_col);
}
goto fin_norm;
}

if(strcmp(Traitement,"5/3")==0){
long int k;
int u;
niveau=niveau-1;
long int dimx_t=dimx/(power(2,niveau)),dimy_t=dimy/power(2,niveau);
for(i=0;i<dimy;i++){
        for(j=0;j<dimx;j++){
           iu[i*dimx+j]=io[i*dimx+j];}}
for(u=tour;u>=1;u--){
e_li=(int *)malloc((2*dimx_t+2)*sizeof(int));
w_li=(int *)malloc((2*dimx_t+4)*sizeof(int));
e_col=(int *)malloc((2*dimy_t+2)*sizeof(int));
w_col=(int *)malloc((2*dimy_t+4)*sizeof(int));
double g_5_3[5]={-0.125,-0.25,0.75,-0.25,-0.125};
double h_5_3[3]={0.5,1,0.5};
double ligne[dimx_t];
double colonne[dimy_t];
double tab_ligne[2*dimx_t];
double tab_colonne[2*dimy_t];

for(i=0;i<dimy;i++){
   Extract_ligne(iu,dimx,i,ligne,u);
   for(k=0;k<dimx_t;k++){
      tab_ligne[2*k]=0;
      tab_ligne[2*k+1]=ligne[k];}
   Conv(h_5_3,tab_ligne,3,2*dimx_t,e_li);
   Conv(g_5_3,tab_ligne,5,2*dimx_t,w_li);
   for(j=0;j<dimx_t;j++){
   inter2[i*dimx+j]=e_li[j]+w_li[j];
}
}

for(j=0;j<dimx_t;j++){
   Extract_colonne(inter2,dimy,j,dimx,colonne,u);
   for(k=0;k<dimy_t;k++){
      tab_colonne[2*k]=0;
      tab_colonne[2*k+1]=colonne[k];}
   Conv(h_5_3,tab_colonne,3,2*dimy_t,e_col);
   Conv(g_5_3,tab_colonne,5,2*dimy_t,w_col);
   for(i=0;i<dimy_t;i++){
   ir[i*dimx+j]=e_col[i]+w_col[i];
}
     
}
for(i=0;i<dimy;i++){
   for(j=0;j<dimx;j++){
      if((i>dimy_t)||(j>dimx_t)){
        ir[i*dimx+j]=io[i*dimx+j];}}}
for(i=0;i<dimy;i++){
   for(j=0;j<dimx;j++){
        iu[i*dimx+j]=ir[i*dimx+j];}}
niveau=niveau-1;
dimx_t=dimx_t*2;
dimy_t=dimy_t*2;
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
free(inter);

fclose(fio);fclose(fir);
return 0;
}