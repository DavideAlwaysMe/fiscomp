#include <stdio.h>
#include <math.h>
#define ALGORITMI_NUM 4

struct valori{
    double x;
    double v;
    double t;
};

struct valori eulero(double dt, double omegaquadro, struct valori valori_n);
struct valori eulerocromer(double dt, double omegaquadro, struct valori valori_n);
struct valori puntocentrale(double dt, double omegaquadro, struct valori valori_n);
struct valori mezzopasso(double dt, double omegaquadro, struct valori valori_n);

int sceglialgoritmo();
double energia(double m,double v,double k,double x);

int main() {
    //dati iniziali da inserire: x0 v0 T dt k m
    struct valori valori_n={.x=5,.v=0,.t=0}; //posizione, velocità e tempo iniziale
    double T=10, dt=0.01, k=30, m=10,npassi=T/dt,energia_n;
    int  i,algoritmo;
    double omegaquadro=k/m;
    FILE *fptr;

    //algoritmo_lista contiene i puntatori alle funzioni dei vari algoritmi
    struct valori (*algoritmo_lista[ALGORITMI_NUM])(double dt, double omegaquadro, struct valori valori_n)={eulero,eulerocromer,puntocentrale,mezzopasso};

    //quale algoritmo verrà utilizzato
    algoritmo=sceglialgoritmo();

    //calcolo i valori di x, v e t con l'algoritmo scelto e salvo in un file di testo
    fptr=fopen("valori.dat","w+");
    fprintf(fptr,"#t         x         v         E\n");
    energia_n=energia(m,valori_n.v,k,valori_n.x);
    fprintf(fptr,"%lf %lf %lf %lf\n",valori_n.t,valori_n.x,valori_n.v,energia_n);
    for(i=0;i<=npassi;i++){
        //calcolo i valori con la funzione prelevata dall'array algoritmi_lista
        valori_n=(*algoritmo_lista[algoritmo])(dt,omegaquadro,valori_n);
        energia_n=energia(m,valori_n.v,k,valori_n.x); //calcolo energia utilizzando la rispettiva energia
        fprintf(fptr,"%lf %lf %lf %lf\n",valori_n.t,valori_n.x,valori_n.v,energia_n); //stampo su un file t, x, v e energia
    }
    fclose(fptr);
}

struct valori eulero(double dt, double omegaquadro, struct valori valori_n){
    struct valori valori_new;
    //la nuova posizione viene calcolata utilizzando la velocità al passo n
    valori_new.x=valori_n.x+valori_n.v*dt;

    //la nuova velocità viene calcolata utilizzando a=-omega^2*x
    valori_new.v=valori_n.v - omegaquadro*valori_n.x*dt;

    valori_new.t=valori_n.t+dt;
    
    return valori_new;
}

struct valori eulerocromer(double dt, double omegaquadro, struct valori valori_n){
    struct valori valori_new;

    //la nuova velocità viene calcolata utilizzando a=-omega^2*x
    valori_new.v=valori_n.v - omegaquadro*valori_n.x*dt;

    //la nuova posizione viene calcolata utilizzando la velocità al passo n+1
    valori_new.x=valori_n.x+valori_new.v*dt;

    valori_new.t=valori_n.t+dt;
    
    return valori_new;
}

struct valori puntocentrale(double dt, double omegaquadro, struct valori valori_n){
    struct valori valori_new;

    //la nuova velocità viene calcolata utilizzando a=-omega^2*x
    valori_new.v=valori_n.v - omegaquadro*valori_n.x*dt;

    //la nuova posizione viene calcolata utilizzando la media tra la velocità al passo n e quella al passo n+1
    valori_new.x=valori_n.x+(valori_new.v+valori_n.v)/2.*dt;

    valori_new.t=valori_n.t+dt;
    
    return valori_new;
}

struct valori mezzopasso(double dt, double omegaquadro, struct valori valori_n){
    struct valori valori_new;

    //la nuova velocità viene calcolata utilizzando a=-omega^2*x
    valori_new.v=valori_n.v - omegaquadro*valori_n.x*dt;

    //la nuova posizione viene calcolata utilizzando la velocità al passo n+1
    valori_new.x=valori_n.x+valori_new.v*dt;

    valori_new.t=valori_n.t+dt;
    
    return valori_new;
}

double energia(double m,double v,double k,double x){
    double e_potenziale,e_cinetica;
    e_potenziale=1./2. * k * pow(x,2);
    e_cinetica=1./2. * m * pow(v,2);
    return e_potenziale+e_cinetica;
}

int sceglialgoritmo(){
    int algoritmo;

    printf("Selezionare l\'algoritmo da utilizzare:\n0 - Eulero\n1 - Eulero-Cromer\n2 - Punto centrale\n3 - Mezzo passo\n");
    scanf(" %d",&algoritmo);
    return algoritmo;
}