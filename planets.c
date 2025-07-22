#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define A 1 //coefficiente di interazione, di default A=1
#define gamma 4*pow(M_PI,2)
#define h 0.0001
#define T 100.
#define muA 3E-6          
#define muB 3.2E-7   //mu Marte 3.2E-7   luna 0.24E-6

struct pianeti
{
    double x,y;
};

struct pianeti f(struct pianeti,struct pianeti,int);
void rk2(struct pianeti,struct pianeti,struct pianeti rb,struct pianeti vb,FILE *fp);


int main()
{
    struct pianeti p0a,p0b,v0a,v0b;
    p0a.x=1.016;
    p0a.y=0.;
    v0a.x=0.;
    v0a.y=6.067;

    p0b.x=1.666;   //marte 1.666 luna 1.02
    p0b.y=0.;
    v0b.x=0.;
    v0b.y=4.64; // marte 4.64 luna 6.1



    FILE *fp;
    fp=fopen("pianeti.dat","w");
    rk2(p0a,v0a,p0b,v0b,fp);
    fclose(fp);
    
}

struct pianeti f(struct pianeti p1,struct pianeti p2,int i)
{
    struct pianeti new;
    double mu;
    long double r = pow(p1.x,2)+pow(p1.y,2);
    long double r12 = pow(p1.x-p2.x,2)+pow(p1.y-p2.y,2);

    if(i==1){mu  = muB;}
    else{mu = muA;}
    new.x = -gamma*p1.x/pow(r,1.5) + A*mu*gamma*(p2.x-p1.x)/pow(r12,1.5);
    new.y = -gamma*p1.y/pow(r,1.5) + A*mu*gamma*(p2.y-p1.y)/pow(r12,1.5);

    return new;
}

void rk2(struct pianeti ra,struct pianeti va,struct pianeti rb,struct pianeti vb,FILE *fp)
{
    struct pianeti rat, vat, sumA;
    struct pianeti rbt, vbt, sumB;

    int i;
    double L = muA*(ra.x*va.y-ra.y*va.x)+muB*(rb.x*vb.y-rb.y*vb.x);
    
 for ( i = 0; i < T/h; i++)
    {
        if (i%100==0){ fprintf(fp,"%.7lf \t %.7lf \t%.7lf \t %.7lf \t %.15lf \t %i \n",ra.x,ra.y,rb.x,rb.y,L,i);}
        
        //Pianeta A
        rat.x=va.x*h;
        rat.y=va.y*h;
        vat.x=(f(ra,rb,1)).x*h;
        vat.y=(f(ra,rb,1)).y*h;
        //Pianeta B
        rbt.x=vb.x*h;
        rbt.y=vb.y*h;
        vbt.x=(f(rb,ra,-1)).x*h;
        vbt.y=(f(rb,ra,-1)).y*h;
        //A
        sumA.x=ra.x+rat.x*0.5;
        sumA.y=ra.y+rat.y*0.5;
        //B
        sumB.x=rb.x+rbt.x*0.5;
        sumB.y=rb.y+rbt.y*0.5;       
        //A
        ra.x=ra.x+(va.x+vat.x*0.5)*h;
        ra.y=ra.y+(va.y+vat.y*0.5)*h;
        //B
        rb.x=rb.x+(vb.x+vbt.x*0.5)*h;
        rb.y=rb.y+(vb.y+vbt.y*0.5)*h;     
        //A
        va.x=va.x+(f(sumA,sumB,1)).x*h;
        va.y=va.y+(f(sumA,sumB,1)).y*h;
        //B
        vb.x=vb.x+(f(sumB,sumA,-1)).x*h;
        vb.y=vb.y+(f(sumB,sumA,-1)).y*h;
        L = muA*(ra.x*va.y-ra.y*va.x)+muB*(rb.x*vb.y-rb.y*vb.x);

        //Check
        if ( (fabs(ra.x-rb.x)<0.001) && (fabs(ra.y-rb.y)<0.001) )
        {
            i=T/h;
            printf("PLANETARY CRASH\n");
        }        
    }
    fprintf(fp,"%.7lf \t %.7lf \t%.7lf \t %.7lf \t %.15lf \t %i \n",ra.x,ra.y,rb.x,rb.y,L,i);
}
