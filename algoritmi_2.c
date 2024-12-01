//ABBIAMO USATO QUESTO SCRIPT PER I METODI EULERO, EULERO-CROMER, PUNTO CENTRALE, MEZZO PASSO
//SI DEVE PASSARE IL COEFFICIENTE D'ATTRITO BETA DA RIGA DI COMANDO

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define k 10
#define m 1
#define dt 0.01
#define T 20
double beta;

struct stato 
{
    double x,v,E;
};
struct delta
{
    double x[4],v[4];
};

void rk2(struct stato s, FILE *fp);
void verlet_auto(struct stato s,FILE *fp);
void rk4(struct stato s, FILE *fp);
void predizione(struct stato s, FILE *fp);

double f(double x, double v);

int main(int argc,char *argv[1]) 
{
        if (argc!=2) 
        {
        fprintf(stderr,"\nExiting... one parameter is needed.\n ");
        fprintf(stderr,"%i parameters have been passed.\n\n",argc -1);
        exit(EXIT_FAILURE); 
        }
        int alg;
        FILE *fp;
        fp = fopen("alg.dat","w");
        struct stato p = {1.,0.,5.};
        beta = (double)atof(argv[1]);

        printf("Inserire algoritmo 1:verl_auto, 2:rk2, 3:rk4, 4:predizione \n");
        scanf("%i",&alg);
        if(alg==1) {verlet_auto(p,fp);}
        if(alg==2) {rk2(p,fp);}
        if(alg==3) {rk4(p,fp);}
        if(alg==4) {predizione(p,fp);}
        
        fclose(fp);
}

//FUNZIONE 
double f(double x, double y)
{
        return -(k/m)*x -beta*y;
}

//ALGORITMI
void rk2(struct stato s, FILE *fp)
{
        int i=0;
        double  xt,vt;
        struct stato new;
        
        fprintf(fp,"%lf %lf %lf %lf\n",0.,s.x,s.v,s.E);
        for(i=0;i<T/dt;i++)
        {
                xt = s.v * dt;
                vt = f(s.x,s.v) * dt;
    
                new.x = s.x + (s.v +vt/2)*dt;
                new.v = s.v +  f(s.x+xt/2,s.v+vt/2) * dt;
        
                s.x = new.x;
                s.v = new.v;
                s.E = (k*s.x*s.x + m*s.v*s.v)/2;

                fprintf(fp,"%lf %lf %lf %lf\n",(i+1)*dt,s.x,s.v,s.E);
        }
}

void verlet_auto(struct stato s,FILE *fp)
{
        int i=0;
        double  xt,vt;
        
        fprintf(fp,"%lf %lf %lf %lf\n",0.,s.x,s.v,s.E);
        for(i=0;i<T/dt;i++)
        {
                xt= s.x + s.v*dt + f(s.x,s.v)*dt*dt/2;
                vt = s.v + (f(s.x,s.v) + f(xt,s.v))*dt/2;
                
                s.x = xt;
                s.v = vt;
                s.E = (k*pow((s.x),2)+m*pow((s.v),2))/2;
                fprintf(fp,"%lf %lf %lf %lf\n",(i+1)*dt,s.x,s.v,s.E);
        }
}

void rk4(struct stato s, FILE *fp)
{
        int i=0;
        struct delta d;

        fprintf(fp,"%lf %lf %lf %lf\n",0.,s.x,s.v,s.E);
        for(i=0;i<T/dt;i++)
        {
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
                fprintf(fp,"%lf %lf %lf %lf\n",(i+1)*dt,s.x,s.v,s.E);
        }
        
}

void predizione(struct stato s, FILE *fp)
{
        int i=0;
        double pre_x,xt,vt;

        fprintf(fp,"%lf %lf %lf %lf\n",0.,s.x,s.v,s.E);
        pre_x = s.x - s.v*dt -f(s.x,s.v)*dt*dt/2 ;
        
        for(i=0;i<T/dt;i++)
        {
                vt = s.v +(f(pre_x + 2*s.v*dt,s.v)+ f(s.x,s.v))*dt/2;
                xt = s.x +(s.v + vt)*dt/2;
                pre_x = s.x;
                s.x = xt;
                s.v = vt;
                s.E = (k*pow((s.x),2)+m*pow((s.v),2))/2; 
                fprintf(fp,"%lf %lf %lf %lf\n",(i+1)*dt,s.x,s.v,s.E);    
        }

}
