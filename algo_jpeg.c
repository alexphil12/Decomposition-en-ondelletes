#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "algo_jpeg.h"

long int min_8(long int dim){
double k;
long int u;
k=dim/8;
if(dim%8==0){
   u=dim;}
else{
u=floor(k)+1;
u=u*8;
return(u);}
}

void discret_cos_8(long int i,long int j,double* inter,long int dimx_8){
int k,l;
double *bloc_in,*bloc_out;
bloc_in=(double *)malloc(sizeof(double)*64);
bloc_out=(double *)malloc(sizeof(double)*64);// implémente la transformée en cos dans les bloc,c'est DCT 8 qui fait la transformée à proprement parler
int u,v;
int c;
double s,pi;
pi=4*atan(1.0);
for(k=0;k<8;k++){
   for(l=0;l<8;l++){
      bloc_in[k*8+l]=inter[(j+k)*dimx_8+i+l];
      bloc_out[k*8+l]=0;
}}
DCT_8(bloc_in,bloc_out);
for(k=0;k<8;k++){
   for(l=0;l<8;l++){
      inter[(j+k)*dimx_8+i+l]=bloc_out[k*8+l];//recopie de la sortie de la transformée
}}
free(bloc_in);
free(bloc_out);
}

void quantif(long int i, long int j,double* inter,double fq,long int dimx_8, int *quantificateur){
int u,p;
double s;
for(p=0;p<8;p++){
   for(u=0;u<8;u++){
      inter[(j+p)*dimx_8+(u+i)]=floor(inter[(j+p)*dimx_8+(u+i)]/quantificateur[p*8+u]);
}}
}
void zigzag(double* inter,int i,int j, long int dimx_8,int* code){

int col[64] ={1,2,1,1,2,3,4,3,2,1,1,2,3,4,5,6,5,4,3,2,1,1,2,3,4,5,6,7,8,7,6,5,4,3,2,1,2,3,4,5,6,7,8,8,7,6,5,4,3,4,5,6,7,8,8,7,6,5,6,7,8,8,7,8};
int lig[64] ={1,1,2,3,2,1,1,2,3,4,5,4,3,2,1,1,2,3,4,5,6,7,6,5,4,3,2,1,1,2,3,4,5,6,7,8,8,7,6,5,4,3,2,3,4,5,6,7,8,8,7,6,5,4,5,6,7,8,8,7,6,7,8,8};
int *block;
block=(int *)malloc(sizeof(int)*64);
int k,l;
for(k=0;k<8;k++){
   for(l=0;l<8;l++){
      block[k*8+l]=inter[(j+k)*dimx_8+i+l]; //les tableaux de zigzag sont issue du TP de majeur.
}}
for(k=0;k<64;k++){
   code[k+i*8+j*dimx_8]=block[(col[k]-1)*8+lig[k]-1];
}
free(block);
}

int rlc(int* code,int *Vlrc,long int dimx_8,long int dimy_8){
int val,count,place;
int compte;
compte=0;
count = 0;
place=0;
int i;
for(i=0;i<dimx_8*dimy_8-1;i++){// le code effectue un compte de zéros successifs puis ajoute 257 et ce nombre,dans tout les autres cas il recopie la valeur lu
   val=code[i];
   if(val==0){
      count=count+1;}
   else{
      if(count>0){
         Vlrc[place]=257;
         Vlrc[place+1]=count;
         compte=compte+count;
         count=0;
         place=place+2;
         }
      Vlrc[place]=val;
      place=place+1;
      }
}
Vlrc=(int *)realloc(Vlrc,(dimx_8*dimy_8-compte)*sizeof(int));
return(dimx_8*dimy_8-compte);
}

void DCT_8(double bloc_in[64], double bloc_out[64]){
    int u=0, v=0;
    int x=0, y=0;
    double au=0, av=0,PI;
    PI=4*atan(1.0);
    for(v=0;v<=7;v++)
    {
        for(u=0;u<=7;u++)
        {
            //Calcule des doubles sommes
            for(y=0;y<=7;y++)
            {
                for(x=0;x<=7;x++)
                {
                    bloc_out[v*8+u]+=bloc_in[y*8+x]*(cos((PI*(2*x+1)*u)/16)*cos((PI*(2*y+1)*v)/16));
                }
            }

            //Calcule des facteurs d'ortogonalité
            if(u==0){au=(1/sqrt(2));}
            else {au=1;}

            if(v==0){av=(1/sqrt(2));}
            else{ av=1;}

            //Calcule du coefficient complet
            bloc_out[v*8+u]=bloc_out[v*8+u]*0.25*au*av;
        }
    }
}

int nb_dif_VLRC(int* Vlrc,dico_huff *verif,int taille_code){
   int k,l,check;
   verif[0].nbr=Vlrc[0];//crée un dictionnaire contenant tout les symbole de VLrc avec leur fréquences classés par fréquence croissante
   int avance;
   avance=1;
   for(k=1;k<taille_code;k++){//Cette partie calcul quels sont les différents symbole de VLRC
      check=0;
      for(l=0;l<avance;l++){
         if(Vlrc[k]==verif[l].nbr){
            check=1;}
      }
      if(check==0){
         verif[avance].nbr=Vlrc[k];
         avance=avance+1;}
   }
   for(k=0;k<avance;k++){//à partir de ces deux boucles, on calcule les fréquences associéà chaque symbole
      verif[k].freq=0;
      for(l=0;l<taille_code;l++){
         if(verif[k].nbr==Vlrc[l]){
             verif[k].freq=verif[k].freq+1;
         }
      }
   }
   dico_huff c;
   int i,j;
   for(i=0;i<avance-1;i++){  //On tri les fréquences par ordres croissant avec un tri par sélection
      for(j=i+1;j<avance;j++){
         if ( verif[i].freq > verif[j].freq) {
            c = verif[i];
            verif[i] = verif[j];
            verif[j] = c;
        }
      }
   }
   return(avance);
   }


void creer_arbre_huff(arbre_code *tab_arbre[],int nbr_dif){
int k,new_freq;//On commence avec le nombre de symbole différent et un tableau contenant les futures feuilles.
int compte,rang;
compte=nbr_dif;
arbre_code *tmp;
long int *adresse;
while(compte!=1){//ici on aurait put faire une boucle for car il y aur toujours exactement nbr_dif-1 noeud
   tmp=malloc(sizeof(arbre_code));
   adresse=malloc(sizeof(long int));
   rang=0;
   new_freq=(tab_arbre[0]->freq)+(tab_arbre[1]->freq);//On créée un nouveau noeud donc la fréquence est la somme de ses branches
   tmp->freq=new_freq;
   tmp->nbr=-1000;//Il permettra la reconnaisance de fin d'algorithme de creer_code_huff
   tmp->branche0=tab_arbre[0];
   tmp->branche1=tab_arbre[1];
   for(k=1;k<compte;k++){
      tab_arbre[k-1]=tab_arbre[k];}
   for(k=1;k<compte;k++){
      tab_arbre[k-1]=tab_arbre[k];}//on décale de 2 pas vers la gauche pour supprimer les noeuds 0 et 1 qui sont les branches du noeud tmp
   for(k=0;k<compte;k++){
      rang+=1;
      if(new_freq<tab_arbre[rang]->freq){//On cherche à insérer tmp dans la liste triée
         break;}
   }
   *adresse=tab_arbre[rang+1];
   for(k=compte-1;k>rang;k--){
      tab_arbre[k]=tab_arbre[k-1];//On replace les autres arbres après l'insertion de tmp
   }
   *tab_arbre[rang]=*tmp;
   tab_arbre[rang+1]=NULL;
   tab_arbre[rang+1]=*adresse;
   tab_arbre[compte-1]=NULL;
   compte=compte-1;
   free(tmp);
   free(adresse);
}
}

int creer_code_huff(arbre_code *noeud,char *code,int nbr_dif,unsigned short pos,dico_encodage *dico_huff,int place, arbre_code *origine){

if((noeud->branche0==NULL)&&(noeud->branche1==NULL)){//Cette algorithme recursif parcourt tout l'arbre jusqu'à ce qu'il tombe sur une feuille 
   if(noeud->nbr==-1000){
      return(place);}//alors sauvegarde le chemin parcouru ainsi que le symbole et la longueur du chemin puis supprime la feuille et repart
   else{
   char *tmp;
   tmp=(char *)malloc(pos*sizeof(char));
   strncpy(tmp,code,pos*sizeof(char));
   dico_huff[place].symb=noeud->nbr;
   dico_huff[place].code=tmp;
   dico_huff[place].taille=pos;
   free(tmp);
   place=place+1;
   pos=0;
   noeud=NULL;
   code="";
   creer_code_huff(origine,code,nbr_dif,(short)pos,dico_huff,place,origine);}
}

else{
   if(noeud->branche0!=NULL){
      code[pos]=48;
      creer_code_huff(noeud->branche0,code,nbr_dif,(short)(pos+1),dico_huff,place,origine);}//selon la valeur des branches on suit l'une ou l'autre
   if(noeud->branche1!=NULL){
      code[pos]=49;
      creer_code_huff(noeud->branche1,code,nbr_dif,(short)(pos+1),dico_huff,place,origine);}
}
}

int decode_huff(dico_encodage *dico,char *sequ,int *Vlrc,int nbr_dif){
int pos,k;
pos=0;
int place;
place=0;
while(sequ[pos]!='\0'){
   for(k=0;k<nbr_dif;k++){
       if(strncmp(dico[k].code,sequ,dico[k].taille)==0){//si la séquence vaut une de taille dico.taille de sequ vaut un code dans le dico,comme le code de huffman est préfixe c'est qu'il code effecivement la séquence donné.
          break;}
    }
    pos+=dico[k].taille;
    Vlrc[place]=dico[k].symb;
    place+=1;
}
return(place);
}

void Vlrc_zigzag(int taille_code,int *Vlrc,int *zigzag){
int k;
int i=0;
int j;
int l;
for(k=0;k<taille_code-1;k++){
   if(Vlrc[k]!=257){
      zigzag[k+i]=Vlrc[k];}//l'inverce de Vlrc
   else{
      j=Vlrc[k+1];
      for(l=0;l<j;l++){
         zigzag[k+l]=0;}
      i+=j;
     }
}
}
void zigzag_inv(int *zigzag, long int dimx, long int dimy,int* iu,long int i,long int j){
int k,l;
int bloc[64];
int col[64] ={1,2,1,1,2,3,4,3,2,1,1,2,3,4,5,6,5,4,3,2,1,1,2,3,4,5,6,7,8,7,6,5,4,3,2,1,2,3,4,5,6,7,8,8,7,6,5,4,3,4,5,6,7,8,8,7,6,5,6,7,8,8,7,8};
int lig[64] ={1,1,2,3,2,1,1,2,3,4,5,4,3,2,1,1,2,3,4,5,6,7,6,5,4,3,2,1,1,2,3,4,5,6,7,8,8,7,6,5,4,3,2,3,4,5,6,7,8,8,7,6,5,4,5,6,7,8,8,7,6,7,8,8};
for(k=0;k<64;k++){
   bloc[lig[k]-1+(col[k]-1)*8]=zigzag[k+i*8+j*dimx];}
for(k=0;k<8;k++){
   for(l=0;l<8;l++){
      iu[(j+k)*dimx+i+l]=bloc[k*8+l];
}}
}
void quantification_inv(double fq,int *iu,long int i,long int j,int *quantificateur,long int dimx_8){
double S_inv;
if(fq<50){
   S_inv=5000/fq;}
else{
   S_inv=200-2*fq;
}
for(int p=0;p<8;p++){
   for(int u=0;u<8;u++){
      iu[(j+p)*dimx_8+(u+i)]=floor(quantificateur[p*8+u]*iu[(j+p)*dimx_8+(u+i)]);
}}
}
void discret_cos_8_inv(long int i,long int j,long int dimx_8,int* iu){
int k,l;
int *bloc_out;
double *bloc_in;
bloc_in=(double *)malloc(sizeof(double)*64);//l'inverse de la dct plus haut coder selon le même canevas
bloc_out=(int *)malloc(sizeof(int)*64);
int u,v;
int c;
double s,pi;
pi=4*atan(1.0);
for(k=0;k<8;k++){
   for(l=0;l<8;l++){
      bloc_in[k*8+l]=iu[(j+k)*dimx_8+i+l];
      bloc_out[k*8+l]=0;
}}
DCT_8_inv(bloc_in,bloc_out);
for(k=0;k<8;k++){
   for(l=0;l<8;l++){
      iu[(j+k)*dimx_8+i+l]=bloc_out[k*8+l];
}}
free(bloc_in);
free(bloc_out);
}

void DCT_8_inv(double block_in[64], int block_out[64])
{
    int u=0, v=0;
    int x=0, y=0;
    double PI;
    PI=4*atan(1.0);
    double au=0, av=0;
    //Matrice intermédiaire pour travailler en nombre réel
    double f[64]={0};

    for(y=0;y<=7;y++)
    {
        for(x=0;x<=7;x++)
        {
            for(v=0;v<=7;v++)
            {
                for(u=0;u<=7;u++)
                {
                    if(u==0){au=1/sqrt(2);}
                    else {au=1;}

                    if(v==0){av=1/sqrt(2);}
                    else{ av=1;}
                    f[y*8+x]+=block_in[v*8+u]*au*av*cos(((2*x+1)*PI*u)/16)*cos(((2*y+1)*PI*v)/16);
                }
            }
            f[y*8+x]*=0.25;
            //Conversion réel entier avec arrondi
            block_out[y*8+x]=floor(f[y*8+x]+0.5);
        }
    }
}


void ecr_huffman(dico_encodage *dico,int *Vlrc,char *huff,int taille_code, int nbr_dif){
int k,p;
for(k=0;K<taille_code;k++){
   p=0;
   while(Vlrc[k]!=dico[p].symb){//choisie la suite du code de huffman associé au symbole Vlrc[k]
      p+=1;}
   strcat(huff,dico[p].code);
}
}
   

















