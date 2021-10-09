#include <stdio.h>
#include <math.h>

struct valori{
    double x;
    double v;
    double t;
};

struct valori eulero(double dt, double omegaquadro, struct valori valori_n);
struct valori eulerocromer(double dt, double omegaquadro, struct valori valori_n);

int sceglialgoritmo();
double energia(double m,double v,double k,double x);

int main() {
    //dati iniziali da inserire: x0 v0 T dt k m
    struct valori valori_n={.x=5,.v=0,.t=0}; //posizione, velocità e tempo iniziale
    double T=10, dt=0.01, k=30, m=10,npassi=T/dt,energia_n;
    int  i,algoritmo;
    double omegaquadro=k/m;
    FILE *fptr;


    //quale algoritmo verrà utilizzato
    algoritmo=sceglialgoritmo();

    //seleziona algoritmo e itera algoritmo npassi volte
    switch (algoritmo){
        case 0:
            fptr=fopen("valori.dat","w+");
            fprintf(fptr,"#t  x   v\n");
            fprintf(fptr,"%lf %lf %lf\n",valori_n.t,valori_n.x,valori_n.v);
            for(i=0;i<=npassi;i++){
                //calcolo i valori con la funzione eulero
                valori_n=eulero(dt,omegaquadro,valori_n);
                energia_n=energia(m,valori_n.v,k,valori_n.x);
                fprintf(fptr,"%lf %lf %lf %lf\n",valori_n.t,valori_n.x,valori_n.v,energia_n);
            }
            fclose(fptr);
            break;

        case 1:
            fptr=fopen("valori.dat","w+");
            //aggiungere dettagli all'iniio del file
            fprintf(fptr,"#t  x   v\n");
            fprintf(fptr,"%lf %lf %lf\n",valori_n.t,valori_n.x,valori_n.v);
            for(i=0;i<=npassi;i++){
                //calcolo i valori con la funzione eulerocromer
                valori_n=eulerocromer(dt,omegaquadro,valori_n);
                energia_n=energia(m,valori_n.v,k,valori_n.x);
                fprintf(fptr,"%lf %lf %lf %lf\n",valori_n.t,valori_n.x,valori_n.v,energia_n);
            }
            fclose(fptr);
            break;
    }


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

double energia(double m,double v,double k,double x){
    double e_potenziale,e_cinetica;
    e_potenziale=1/2 * k * pow(x,2);
    e_cinetica=1/2 * m * pow(v,2);
    return e_potenziale+e_cinetica;
}

int sceglialgoritmo(){
    int algoritmo;

    printf("Selezionare l\'algoritmo da utilizzare:\n0 - Eulero\n1 - Eulero-Cromer\n2 - Punto centrale\n3 - Mezzo passo\n");
    scanf(" %d",&algoritmo);
    return algoritmo;
}