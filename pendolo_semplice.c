///ABBIAMO USATO QUESTO SCRIPT PER VERLET AUTOSUFFICIENTE E RUNGE-KUTTA AL II ORDINE
#include <stdio.h>
#include <math.h>

#define g 9.81
#define l 1
#define m 1
#define dt 0.01
#define T 10

struct stato 
{
    double x,v,E; //dove x = angolo e v = velocita ang
};

struct delta
{
    double x[4],v[4];
};

struct stato rk2(struct stato s,double(*f)(double));
struct stato rk4(struct stato s,double(*f)(double));
struct stato predizione(struct stato s,double *pre_x,double(*f)(double));
struct stato verlet_auto(struct stato s,double(*f)(double));
double f(double x);

int max(double v[]);

int main()
{
    FILE *fp;
    fp = fopen("ps.dat","w");
    
    int i,algoritmo;
    //VALORI INIZIALI
    struct stato p = {-M_PI/2,0.,9.81};
    //fprintf(fp,"%lf %lf %lf %lf\n",0.,p.x,p.v,p.E);
    fprintf(fp,"%lf %lf %lf %lf %lf\n",0.,p.x,l*sin(p.x),-l*cos(p.x),p.E);  
    ///////
    printf("Inserire algoritmo  1:rk2, 1:rk4, 3:predizione, 4:verl_auto \n");
    scanf("%i",&algoritmo);
    
   

   if (algoritmo == 1)
   {
       for(i=0;i<T/dt;i++)
        {
            p = rk2(p,f);
            fprintf(fp,"%lf %lf %lf %lf %lf\n",(i+1)*dt,p.x,l*sin(p.x),-l*cos(p.x),p.E);
        }
   }
   
    if (algoritmo == 2)
    {
        for(i=0;i<T/dt;i++)
        {
            p = rk4(p,f); 
            fprintf(fp,"%lf %lf %lf %lf %lf\n",(i+1)*dt,p.x,l*sin(p.x),-l*cos(p.x),p.E);
        }
    }
    if (algoritmo == 3)
    {
        double pre_x = p.x - p.v*dt -f(p.x)*dt*dt/2 ;
        for(i=1;i<T/dt;i++)
        {
            p = predizione(p,&pre_x,f); 
            fprintf(fp,"%lf %lf %lf %lf %lf\n",(i+1)*dt,p.x,l*sin(p.x),-l*cos(p.x),p.E);
        }
    }

     if (algoritmo == 4)
   {
       for(i=0;i<T/dt;i++)
        {
            p = verlet_auto(p,f);
            fprintf(fp,"%lf %lf %lf %lf %lf\n",(i+1)*dt,p.x,l*sin(p.x),-l*cos(p.x),p.E);
        }
   }
    
    
    
    fclose(fp);
}



////FUNZIONE
double f(double x)
{
    return -g/l * sin(x);
}

///ALGORITMI
struct stato rk2(struct stato s,double(*f)(double))
{
    double fi = f(s.x);
    struct stato temp, new;
    
    temp.x = s.v*dt;
    temp.v = fi*dt;
    
    new.x = s.x + (s.v +temp.v/2)*dt;
    fi = f(s.x+temp.x/2); //Fi che tiene conto delle dei nuovi t,v
    new.v = s.v + fi*dt;
    new.E = m*pow(l*s.x,2)/2 + m*g*l*(1-cos(s.x));
    
    return new;
}


struct stato rk4(struct stato s,double(*f)(double))
{
    struct delta d;
    struct stato new;
    //
    d.v[0] = f(s.x) *dt;
    d.x[0] = s.v *dt; 
    //
    d.v[1] = f(s.x + d.x[0]/2)*dt; 
    d.x[1] = (s.v + d.v[0]/2)*dt;
    //
    d.v[2] = f(s.x + d.x[1]/2)*dt; 
    d.x[2] = (s.v + d.v[1]/2)*dt;
    //
    d.v[3] = f(s.x + d.x[2])*dt; 
    d.x[3] = (s.v + d.v[2])*dt;
    //
    new.x = s.x + (d.x[0]+ 2*d.x[1]+ 2*d.x[2]+ d.x[3]) * 1/6;
    new.v = s.v + (d.v[0]+ 2*d.v[1]+ 2*d.v[2]+ d.v[3]) * 1/6;
    new.E = m*pow(l*s.x,2)/2 + m*g*l*(1-cos(s.x));
    return new;
}

struct stato predizione(struct stato s,double *pre_x,double(*f)(double))
{
    struct stato new;
    double fi;
    fi = f(*pre_x + 2*s.v*dt);
    
    new.v = s.v +(fi+ f(s.x))*dt/2;
    new.x = s.x +(s.v + new.v)*dt/2;
    new.E = m*pow(l*s.x,2)/2 + m*g*l*(1-cos(s.x));
    *pre_x = s.x;
    return new;
}

struct stato verlet_auto(struct stato s,double(*f)(double))
{
  struct stato new;
  
  new.x= s.x + s.v*dt +f(s.x) *dt*dt/2;

  new.v = s.v + (f(s.x) +f(new.x))*dt/2;
  new.E = m*pow(l*s.x,2)/2 + m*g*l*(1-cos(s.x));

  return new;
}




