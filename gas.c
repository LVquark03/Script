#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L 80
int n=0,T;
float ro;

void movimento(int start[],int pos[],int pos_true[],int reticolo[]);

void inizializzazione(int reticolo[]);
void inizio(int start[],int pos[],int pos_true[],int reticolo[]);

void stampa(int reticolo[]);
void print_array(int v[],int dim);

///funzioni di prova
int count(int array[]);

int main(int argc,char *argv[2])
{
ro = (float)atof(argv[1]);
T = (int)atoi(argv[2]);
int seed=10000001;
srand48(seed);
    
int count=0;
FILE *fp;
fp=fopen("DR2_L80.dat","w");

int *start,*pos,*pos_true;
long int DR2=0;
float coeff_diff;
int reticolo[L*L]={0};
inizializzazione(reticolo);

start = (int*)malloc(2*n*sizeof(int));
pos = (int*)calloc(2*n,sizeof(int));
pos_true = (int *)calloc(2*n,sizeof(int));
inizio(start,pos,pos_true,reticolo);

stampa(reticolo);

printf("\n######\n");
for(int i=1;i<T+1;i++) 
{
  DR2=0;
  movimento(start,pos,pos_true,reticolo);
  for(int j=0;j<n;j++){DR2+=  pow(pos_true[j] -start[j],2) + pow(pos_true[j+n] -start[j+n],2);}
  coeff_diff = (float) DR2/(4*n*i);
  fprintf(fp,"%i \t%f\n",i, coeff_diff);
}
printf("\n######\n");
  
stampa(reticolo);
fclose(fp);
for(int i=0;i<L*L;i++)
{
  count+=reticolo[i];
}
printf("%i \t %i\n",n,count);
}



void inizializzazione(int reticolo[])
{
  for(int i=0; i<L*L;i++)
  {
    if( (double)lrand48()/RAND_MAX <ro)
    {
    reticolo[i]=1;
    n++;
    }	 
  }
}  

void inizio(int start[],int pos[],int pos_true[],int reticolo[])
{
  int count=-1;
  for(int i=0;i<L;i++)
  {
    for(int j=0;j<L;j++)
    {
      if (reticolo[i*L+j]==1) 
      {
      count++;
      start[count]=j;
      start[count+n]=i;
      pos[count]=j;
      pos[count+n]=i;
      pos_true[count]=j;
      pos_true[count+n]=i;
      }
    }

  }
}


void stampa(int reticolo[])
{
  for(int i=0;i<L;i++)
    {
      for(int j=0;j<L;j++)
	{
	  printf("%i ",reticolo[L*i+j]);
	}
      printf("\n");

    }
}

void print_array(int v[],int dim)
{
  printf("\n");
  for(int i=0;i<dim;i++) {printf("%i ",v[i]);}
  printf("\n");
}


void movimento(int start[],int pos[],int pos_true[],int reticolo[])
{
  double r;
  int new_x,new_y,true_x,true_y;
  for(int i=0;i<n;i++)
    {
    r=(double)lrand48()/RAND_MAX;
    new_x=pos[i];
    new_y=pos[i+n];
    true_x = 0;
    true_y = 0;
    if (r<0.25)
    { 
    new_y=(L+pos[i+n] +1)%L; 
    true_y = 1;
    }
    if (r>0.25 && r<0.5) 
    {
    new_x=(L+pos[i] -1)%L;
    true_x = -1;
    }
    if (r>0.5 && r<0.75)
    {
    new_y=(L+pos[i+n] -1)%L;
    true_y = -1;
    }
    if (r>0.75)
    {
    new_x=(L+pos[i] +1)%L;
    true_x = 1;
    }

    if(reticolo[new_x+new_y*L]==0)
    {
    reticolo[pos[i]+pos[i+n]*L]=0;
    pos[i] = new_x;
    pos[i+n] = new_y;
    pos_true[i] += true_x;
    pos_true[i+n] += true_y;
    reticolo[pos[i]+pos[i+n]*L]=1;
    }
    }  
} 


int count(int array[])
{
   int sum=0;
  for(int i=0;i<L*L;i++){sum+=array[i];}
  return sum;
}

