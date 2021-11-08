#include <stdio.h>
#include <math.h>
#include <stdlib.h>

struct Point{
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
};

struct Point rungekutta(double dt, double a, double b, double c, struct Point point);

double fx(double x, double y,double z,double a,double b,double c);
double fy(double x, double y,double z,double a,double b,double c);
double fz(double x, double y,double z,double a,double b,double c);

int main(int argc, char* argv[]){
    if(argc!=9){
        fprintf(stderr,"Per l\'esecuzione del programma è necessario passare come argomenti: a, b, c, x0, y0, z0, dt, T.\nSarà utilizzato l'algoritmo di Runge Kutta.\n");
        exit(1);
    }

    int i,npassi;
    double a,b,c,x0,y0,z0,dt,T,tempo=0;
    struct Point point;

    a=atof(argv[1]);
    b=atof(argv[2]);
    c=atof(argv[3]);
    x0=atof(argv[4]);
    y0=atof(argv[5]);
    z0=atof(argv[6]);
    dt=atof(argv[7]);
    T=atof(argv[8]);
    npassi=T/dt;

    point.x=x0;
    point.y=y0;
    point.z=z0;
    point.vx=0;
    point.vy=0;
    point.vz=0;

    fprintf(stdout,"#Integrazione con metodo di Runge Kutta, a=%g b=%g c=%g x0=%g y0=%g z0=%g dt=%g T=%g\n",a,b,c,x0,y0,z0,dt,T);
    fprintf(stdout,"#    T                    x                    y                    z                   vx                   vy                   vz\n");
    fprintf(stdout,"%6.10g %20.10g %20.10g %20.10g %20.10g %20.10g %20.10g\n",tempo,point.x,point.y,point.z,point.vx,point.vy,point.vz);
    for(i=1;i<=npassi;i++){
        point=rungekutta(dt,a,b,c,point);
        tempo=i*dt;
        fprintf(stdout,"%6.10g %20.10g %20.10g %20.10g %20.10g %20.10g %20.10g\n",tempo,point.x,point.y,point.z,point.vx,point.vy,point.vz);
    }
    

}

struct Point rungekutta(double dt, double a, double b, double c, struct Point point){
    struct Point point_new;
    double deltaX,deltaY,deltaZ;

    deltaX=fx(point.x,point.y,point.z,a,b,c)*dt;
    deltaY=fy(point.x,point.y,point.z,a,b,c)*dt;
    deltaZ=fz(point.x,point.y,point.z,a,b,c)*dt;

    point_new.x=point.x+fx(point.x+deltaX/2.,point.y+deltaY/2.,point.z+deltaZ/2.,a,b,c)*dt;
    point_new.y=point.y+fy(point.x+deltaX/2.,point.y+deltaY/2.,point.z+deltaZ/2.,a,b,c)*dt;
    point_new.z=point.z+fz(point.x+deltaX/2.,point.y+deltaY/2.,point.z+deltaZ/2.,a,b,c)*dt;

    //salvo anche le velocità perchè serviranno per la sezione di poincare
    point_new.vx=fx(point.x,point.y,point.z,a,b,c);
    point_new.vy=fy(point.x,point.y,point.z,a,b,c);
    point_new.vz=fz(point.x,point.y,point.z,a,b,c);

    return point_new;
}

double fx(double x, double y,double z,double a,double b,double c){
    return -y-z;
}

double fy(double x, double y,double z,double a,double b,double c){
    return x+a*y;
}

double fz(double x, double y,double z,double a,double b,double c){
    return b+(x-c)*z;
}