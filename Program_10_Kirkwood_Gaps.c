//  Program 10 Kirkwood Gaps

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

//Runge Kutta function

double accx(double,double,double,double,double);
double accy(double,double,double,double,double);

int main(void) {
    setvbuf(stdout,NULL,_IONBF,0); setvbuf(stderr,NULL,_IONBF,0);
    
    int tb=50;
    double x[50]={0},y[50]={0},vx[50]={0},vy[50]={0},m[50]={0};
    double ax[50]={0},avx[50]={0},ay[50]={0},avy[50]={0};
    double bx[50]={0},bvx[50]={0},by[50]={0},bvy[50]={0};
    double cx[50]={0},cvx[50]={0},cy[50]={0},cvy[50]={0};
    double dx[50]={0},dvx[50]={0},dy[50]={0},dvy[50]={0};
    
    double t=0, dt=3000, day=86400, year=365, mom=0, start, stop, timer;           // time & orbit settings
    double au=1.495978707e11, g=6.67408e-11, dist=0;  // constants
    int p,s,planets,bodies,counter=0, counter3=0;
    double counter2=0;
    
    x[0]=0,             m[0]=1.98892e30;   // Sun
    x[1]=0.387098*au,   m[1]=3.3011e23;    // Mercury
    x[2]=0.723332*au,   m[2]=4.8675e24;    // Venus
    x[3]=1.000000*au,   m[3]=5.9742e24;    // Earth
    x[4]=1.523679*au,   m[4]=6.4171e23;    // Mars
    x[5]=5.20260*au,    m[5]=1.8986e27;    // Jupiter
    x[6]=9.55490*au,    m[6]=5.6836e26;    // Saturn
    x[7]=19.2184*au,    m[7]=8.6810e25;    // Uranus
    x[8]=30.1104*au,    m[8]=1.0243e26;    // Neptune
    x[9]=10.5529*au,        m[9]=7.34767309e22;      // asteroid
    
    //start timer
    start = clock()/CLOCKS_PER_SEC;
    
    // setting asteroid positions and velocities
    dist=2.35;
    for (tb=10; tb<=49; tb++)
        {
            x[tb]=au*dist;
            m[tb]=1000;
            dist+=0.01;
        }
    for(tb=0; tb<=49; tb++)
    {printf ("p=%d\tx=%g\ty=%g\tvx=%g\tvy=%g\tm=%g\n",tb,x[tb]/au,y[tb],vx[tb],vy[tb],m[tb]);
    }
    
    planets=49;//planets
    bodies=9;//bodies
    
    FILE * finout; //opening data file
    finout = fopen("data.txt","w");
    
    // positions of asteroids
    for (p=0;p<=49;p++)
    {
        fprintf (finout,"%g au\t%g au\t",x[p]/au,x[p]/au);
    }
    
    fprintf (finout,"\n");
   
    // Velocity of planets
    for (p=1; p<=planets; p++)
    {
        vy[p]=sqrt((g*m[0])/(fabs(x[p]-x[0]) ));;
        mom+=(m[p]*vy[p]);
    }
    
    // velocity of Sun
    vy[0]=(-mom/m[0]);

    // simulation for n year
    for (t=dt ; t<year*day*50 ; t+=dt)
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
        
        for (p=0; p<=planets;p++){
            
            vx[p] +=(dt/6)*(avx[p]+(2*bvx[p])+(2*cvx[p])+dvx[p]);
            x[p] +=(dt/6)*(ax[p]+(2*bx[p])+(2*cx[p])+dx[p]);
            
            vy[p] += (dt/6)*(avy[p]+(2*bvy[p])+(2*cvy[p])+dvy[p]);
            y[p] += (dt/6)*(ay[p]+(2*by[p])+(2*cy[p])+dy[p]);
            
        }
        
        // counter used to ensure program is running for long periods
        counter3+=1;
        if (counter3==50000){
                printf ("Time=%g years\n",t/(86400.0*365.2));
                counter3=0;
        }

        // counters used to simulate n years before printing out to file
        counter2=4994*year*day;

        if (counter2<t){
            counter+=1;
            
                // used to print out the data to file for each planet per step
            if (counter==100){
            
                for (p=0;p<=49;p++){
                    fprintf (finout,"%g\t%g\t",x[p]/au,y[p]/au);
                }
                
                fprintf (finout,"\n");

                counter=0;
            }
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







