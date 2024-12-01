//ABBIAMO USATO QUESTO SCRIPT PER I METODI EULERO, EULERO-CROMER, PUNTO CENTRALE, MEZZO PASSO

#include <stdio.h>
#include <math.h>

#define k 10
#define m 1
#define dt 0.01
#define T 100

void eulero(float *x, float *v,float *E);
void eu_cro(float *x, float *v,float *E);
void pc(float *x, float *v,float *E);
    
int main()
{
    float v,x,E;
    int i,algoritmo;
    
    //APERTURA FILE
    FILE *fp;
    fp = fopen("alg.dat","w");
    fprintf(fp,"#t \t\t x \t\t v \n");
    
    //Scelta algoritmo
    printf("Inserire algoritmo 1:e, 2:ec, 3:pc, 4:mp \n");
    scanf("%i",&algoritmo);
    
    v = 0;
    x = 1;
    E = 5.;
    fprintf(fp,"%f \t %f \t %f \t %f\n",0.,x,v,E);
   
   if (algoritmo == 1)
   {
       for(i=0;i<T/dt;i++)
       {
           eulero(&x,&v,&E);
           fprintf(fp,"%f \t %f \t %f \t %f \n",(i+1)*dt,x,v,E);
       }
   }
   if (algoritmo == 2)
   {
       for(i=0;i<T/dt;i++)
       {
           eu_cro(&x,&v,&E);
           fprintf(fp,"%f \t %f \t %f \t %f \n",(i+1)*dt,x,v,E);
       }
   }
    if (algoritmo == 3)
    {
        for(i=0;i<T/dt;i++)
        {
            pc(&x,&v,&E);
            fprintf(fp,"%f \t %f \t %f \t %f \n",(i+1)*dt,x,v,E);
        }
    }
    if (algoritmo == 4) ///mezzo passo
    {
        v = v -(k/m)*x*dt/2;
        for(i=1;i<T/dt;i++)
        {
            eu_cro(&x,&v,&E);
            fprintf(fp,"%f \t %f \t %f \t %f \n",(i+1)*dt,x,v,E);
        }
    }
    
    fclose(fp);    
}


void eulero(float *x, float *v,float *E)
{ 
    float fi = -k/m *(*x);
    *x = *x + (*v)*dt;
    *v = *v + fi*dt;
    *E = (k*(*x)*(*x) + m*(*v)*(*v))/2;
}

void eu_cro(float *x, float *v,float *E)
{ 
    float fi = -k/m *(*x);
    *v = *v + fi*dt;
    *x = *x + (*v)*dt;
    *E = (k*(*x)*(*x) + m*(*v)*(*v))/2;
}

void pc(float *x, float *v,float *E)   //passo centrale
{
    float v1 = (*v);
    float fi = -k/m *(*x);
    *v = *v + fi*dt;
    *x = *x + ((*v)+v1)*dt/2;
    *E = (k*(*x)*(*x) + m*(*v)*(*v))/2;
    
}

