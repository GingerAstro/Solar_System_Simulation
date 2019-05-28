// Program 8 N System with Sun wobble and asteroid

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//Runge Kutta function

double accx(double,double,double,double,double);
double accy(double,double,double,double,double);

int main(void) {
    setvbuf(stdout,NULL,_IONBF,0); setvbuf(stderr,NULL,_IONBF,0);
    
    
    double x[10],y[10]={0},vx[10]={0},vy[10],m[10];
    double ax[10]={0},avx[10]={0},ay[10]={0},avy[10]={0};
    double bx[10]={0},bvx[10]={0},by[10]={0},bvy[10]={0};
    double cx[10]={0},cvx[10]={0},cy[10]={0},cvy[10]={0};
    double dx[10]={0},dvx[10]={0},dy[10]={0},dvy[10]={0};
    
    double t=0, dt=1000, day=86400, year=365, mom, start, stop, timer;           // time & orbit settings
    double au=1.495978707e11, g=6.67408e-11;   // constants
    int p,s,planets,bodies,counter=0;
    
    x[0]=0,             m[0]=1.98892e30;   // Sun
    x[1]=0.387098*au,   m[1]=3.3011e23;    // Mercury
    x[2]=0.723332*au,   m[2]=4.8675e24;    // Venus
    x[3]=1.000000*au,   m[3]=5.9742e24;    // Earth
    x[4]=1.523679*au,   m[4]=6.4171e23;    // Mars
    x[5]=5.20260*au,    m[5]=1.8986e27;    // Jupiter
    x[6]=9.55490*au,    m[6]=5.6836e26;    // Saturn
    x[7]=19.2184*au,    m[7]=8.6810e25;    // Uranus
    x[8]=30.1104*au,    m[8]=1.0243e26;    // Neptune
    x[9]=30*au,     m[9]=0.95e21;      // Asteroid
    
    planets=9;//planets
    bodies=9;//bodies
    
    FILE * finout; //opening data file
    finout = fopen("data.txt","w");
    
    //start timer
    start = clock()/CLOCKS_PER_SEC;
    
    // setting initial velocities
    for (p=1; p<=planets-1; p++)
    {
        vy[p]=sqrt((g*m[0])/(fabs(x[p]-x[0]) ));
        printf ("vy[%d]=%g\n",p,vy[p]);
        mom+=(m[p]*vy[p]);
    }
    
    // setting velocity of asteroid
    vy[9]=0.3*sqrt((g*m[0])/(fabs(x[9]-x[0]) ));
    printf ("vy[%d]=%g\n",p,vy[p]);
    mom+=(m[p]*vy[p]);
    
    // setting velocity of Sun
    vy[0]=(-mom/m[0]);
    printf("vy[0]=%g\tmoment=%g\n",  vy[0],mom);
    
    // simulation for n year
    for (t=dt ; t<year*day*10 ; t+=dt)
    {
        
        //************ calculating a
        for (p=0; p<=planets;p++)
        {
            ax[p]=vx[p];
            ay[p]=vy[p];
            avx[p]=0; avy[p]=0;
            bvx[p]=0; bvy[p]=0;
            cvx[p]=0; cvy[p]=0;
            dvx[p]=0; dvy[p]=0;
            
            for (s=0; s<=bodies; s++)
            {
                if(p==s) continue;
                else{
                    avx[p]+=accx(x[p],y[p],x[s],y[s],m[s]);
                    avy[p]+=accy(x[p],y[p],x[s],y[s],m[s]);
                }
            }
        }
        
        //************ calculating b
        for (p=0; p<=planets;p++)
        {
            bx[p]=vx[p]+(dt/2)*avx[p];
            by[p]=vy[p]+(dt/2)*avy[p];
            for (s=0; s<=bodies; s++)
            {
                if(p==s) continue;
                else{
                    bvx[p]+=accx(x[p]+(dt/2)*ax[p], y[p]+(dt/2)*ay[p], x[s]+(dt/2)*ax[s], y[s]+(dt/2)*ay[s] ,m[s]);
                    bvy[p]+=accy(x[p]+(dt/2)*ax[p], y[p]+(dt/2)*ay[p], x[s]+(dt/2)*ax[s], y[s]+(dt/2)*ay[s] ,m[s]);
                }
            }
        }
        
        //************ calculating c
        for (p=0; p<=planets;p++)
        {
            cx[p]=vx[p]+(dt/2)*bvx[p];
            cy[p]=vy[p]+(dt/2)*bvy[p];
            for (s=0; s<=bodies; s++)
            {
                if(p==s) continue;
                else{
                    cvx[p]+=accx(x[p]+(dt/2)*bx[p], y[p]+(dt/2)*by[p], x[s]+(dt/2)*bx[s], y[s]+(dt/2)*by[s] ,m[s]);
                    cvy[p]+=accy(x[p]+(dt/2)*bx[p], y[p]+(dt/2)*by[p], x[s]+(dt/2)*bx[s], y[s]+(dt/2)*by[s] ,m[s]);
                }
            }
        }
        
        //************ calculating d
        for (p=0; p<=planets;p++)
        {
            dx[p]=vx[p]+(dt*cvx[p]);
            dy[p]=vy[p]+(dt*cvy[p]);
            for (s=0; s<=bodies; s++)
            {
                if(p==s) continue;
                else{
                    dvx[p]+=accx(x[p]+(dt)*cx[p], y[p]+(dt)*cy[p], x[s]+(dt)*cx[s], y[s]+(dt)*cy[s] ,m[s]);
                    dvy[p]+=accy(x[p]+(dt)*cx[p], y[p]+(dt)*cy[p], x[s]+(dt)*cx[s], y[s]+(dt)*cy[s] ,m[s]);
                }
            }
        }
        
        for (p=0; p<=planets;p++)
        {
            
            vx[p] +=(dt/6)*(avx[p]+(2*bvx[p])+(2*cvx[p])+dvx[p]);
            x[p] +=(dt/6)*(ax[p]+(2*bx[p])+(2*cx[p])+dx[p]);
            
            vy[p] += (dt/6)*(avy[p]+(2*bvy[p])+(2*cvy[p])+dvy[p]);
            y[p] += (dt/6)*(ay[p]+(2*by[p])+(2*cy[p])+dy[p]);
            
        }
        counter+=1;
        
        // used to print out the data to file for each planet per dt
        if (counter==1000){
            
            for (p=0;p<=9;p++)
            {
                fprintf (finout,"%g\t%g\t",x[p]/au,y[p]/au);
            }
            
            fprintf (finout,"\n");
            // printf ("time=%g\n",t/(86400*365));
            counter=0;
        }
    }
    
    //stoping timer and printing result
    stop = clock()/CLOCKS_PER_SEC;
    
    timer = stop - start;
    printf ("Time for program to run = %g seconds\n",timer);
    return 0;
    
}

double accx(double xplanet, double yplanet, double xb, double yb, double m)
{
    //  function calculates a values for both x and y, then b values etc. r changes each time to
    //  adjust for the change in step using the RK method.
    
    double r,x;
    double g = 6.67408e-11;

    r=sqrt( (xplanet-xb)*(xplanet-xb)+(yplanet-yb)*(yplanet-yb) );
    x=(g*m*(xb-xplanet))/(r*r*r);
    
    return x;
}

double accy(double xplanet, double yplanet, double xb, double yb, double m)
{
    //  function calculates a values for both x and y, then b values etc. r changes each time to
    //  adjust for the change in step using the RK method.
    
    double r,y;
    double g = 6.67408e-11;

    r=sqrt( (xplanet-xb)*(xplanet-xb)+(yplanet-yb)*(yplanet-yb) );
    y=(g*m*(yb-yplanet))/(r*r*r);
    
    return y;
}







