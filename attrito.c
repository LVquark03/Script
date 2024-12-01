//Script usato per il calcolo del coefficiente beta critico oscillatore armonico
//IMPIEGA CIRCA 10 SECONDI PER L'ESECUZIONE
//VALORI USATI PER IL CALCOLO X0 = 1, V0=0 da passare con riga di comando
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define k 10
#define m 1
#define dt 0.001
#define T 30

double beta = 0.;
int status = 0;


struct stato 
{
    double x,v,E;
};
struct delta
{
    double x[4],v[4];
};


void rk4(struct stato s);

double f(double x, double v);

int main(int argc,char *argv[2])
{
    if (argc!=3) 
    {
    fprintf(stderr,"\nExiting... tree parameters are needed.\n ");
    fprintf(stderr,"%i parameters have been passed.\n\n",argc -1);
    exit(EXIT_FAILURE); 
    }
    int alg;
    double x0,v0,E0;
    x0=(double)atof(argv[1]);
    v0=(double)atof(argv[2]);
    E0=(k*pow((x0),2)+m*pow((v0),2))/2;
    struct stato p = {x0,v0,E0};

    while(status == 0) 
    {
        rk4(p);
        beta += 0.0001;
    }

    printf("beta critico circa = %f \n", (beta-0.00005));
        
        
}

//FUNZIONE 
double f(double x, double y)
{
        return -(k/m)*x -beta*y;
}

//ALGORITMI

void rk4(struct stato s)
{
        int i=0;
        struct delta d;
        // per l'attrito
        double x1;
        int q=1;

        for(i=0;i<T/dt;i++)
        {
            x1 = s.x;
            d.v[0] = f(s.x,s.v) *dt;
            d.x[0] = s.v *dt; 
            //
            d.v[1] = f(s.x + d.x[0]/2,s.v + d.v[0]/2)*dt; 
            d.x[1] = (s.v + d.v[0]/2)*dt;
            //
            d.v[2] = f(s.x + d.x[1]/2,s.v + d.v[1]/2)*dt; 
            d.x[2] = (s.v + d.v[1]/2)*dt;
            //
            d.v[3] = f(s.x + d.x[2],s.v + d.v[2])*dt; 
            d.x[3] = (s.v + d.v[2])*dt;
            //
            s.x = s.x + (d.x[0]+ 2*d.x[1]+ 2*d.x[2]+ d.x[3]) * 1/6;
            s.v = s.v + (d.v[0]+ 2*d.v[1]+ 2*d.v[2]+ d.v[3]) * 1/6;
            s.E = (k*pow((s.x),2)+m*pow((s.v),2))/2; 

            
            if(s.x<0.)    
                {
                    i=T/dt;
                    q=0;    
                }  
        } 
        if(q==1){status = 1;}
}


