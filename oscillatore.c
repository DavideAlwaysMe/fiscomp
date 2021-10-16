#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define ALGORITMI_NUM 4

struct valori{
    double x;
    double v;
};

struct valori eulero(double dt, double omegaquadro, struct valori valori_n);
struct valori eulerocromer(double dt, double omegaquadro, struct valori valori_n);
struct valori puntocentrale(double dt, double omegaquadro, struct valori valori_n);
struct valori mezzopasso(double dt, double omegaquadro, struct valori valori_n);

int sceglialgoritmo();
double energia(double m,double v,double k,double x);

int main(int argc, char* argv[]) {

    if(argc!=7){
        fprintf(stderr,"Per l\'esecuzione del programma è necessario passare come argomenti: x0, v0, dt, T, k, m.\nIl programma eseguira tutti gli algoritmi e salverà i dati in %d file del tipo Eulero.dat.\n",ALGORITMI_NUM);
        exit(1);
    }


    //dati iniziali da inserire: x0 v0 T dt k m
    struct valori valori_n; //posizione, velocità e tempo iniziale
    int  i,algoritmo;
    double tempo=0, T, dt, k, m,npassi,energia_n,energia_0,energia_rapporto=0.,omegaquadro;
    FILE *fptr;
    //t x e v sono array dinamici nei quali vengono salvati i valori di t x e v
    double *t,*x,*v;

    //assegno i parametri di esecuzione alle variabili iniziali
    valori_n.x=atof(argv[1]);
    valori_n.v=atof(argv[2]);
    dt=atof(argv[3]);
    T=atof(argv[4]);
    k=atof(argv[5]);
    m=atof(argv[6]);
    npassi=T/dt;
    omegaquadro=k/m;

    //uso malloc per assegnare le giuste celle di memoria agli array x,v e t
    x=(double*)malloc(sizeof(double)*npassi);
    v=(double*)malloc(sizeof(double)*npassi);
    t=(double*)malloc(sizeof(double)*npassi);

    //algoritmo_lista contiene i puntatori alle funzioni dei vari algoritmi
    struct valori (*algoritmo_lista[ALGORITMI_NUM])(double dt, double omegaquadro, struct valori valori_n)={eulero,eulerocromer,puntocentrale,mezzopasso};

    //calcolo E(0)
    energia_0=energia(m,valori_n.v,k,valori_n.x);

    //quale algoritmo verrà utilizzato
    algoritmo=sceglialgoritmo();

    //se l'algoritmo è mezzopasso devo calcolare v(1/2) con eulero prima di eseguirlo
    if(algoritmo==3){
        valori_n.v=valori_n.v+0.5*(- omegaquadro*valori_n.x*dt);
    }

    //calcolo i valori di x, v e t con l'algoritmo scelto e salvo in un file di testo
    fptr=fopen("valori.dat","w+");
    fprintf(fptr,"#t         x         v         E        delta_E/E(0)\n");
    fprintf(fptr,"%g %g %g %g %g\n",tempo,valori_n.x,valori_n.v,energia_0,energia_rapporto);
    for(i=1;i<=npassi;i++){
        //calcolo i valori con la funzione prelevata dall'array algoritmi_lista
        valori_n=(*algoritmo_lista[algoritmo])(dt,omegaquadro,valori_n);
        energia_n=energia(m,valori_n.v,k,valori_n.x); //calcolo energia utilizzando la rispettiva funzione energia
        energia_rapporto=(energia_n-energia_0)/energia_0;
        tempo=(double)i*dt; //calcolo il tempo così per non perdere precisione per via dell'approssimazione di numeri molto piccoli
        fprintf(fptr,"%g %g %g %g %g\n",tempo,valori_n.x,valori_n.v,energia_n,energia_rapporto); //stampo su un file t, x, v e energia
        x[i-1]=valori_n.x; //i-1 perché i parte da 1 e non da 0
        v[i-1]=valori_n.v;
        t[i-1]=tempo;
    }
    fclose(fptr);
    return(0);
}

struct valori eulero(double dt, double omegaquadro, struct valori valori_n){
    struct valori valori_new;
    //la nuova posizione viene calcolata utilizzando la velocità al passo n
    valori_new.x=valori_n.x+valori_n.v*dt;

    //la nuova velocità viene calcolata utilizzando a=-omega^2*x
    valori_new.v=valori_n.v - omegaquadro*valori_n.x*dt;

    return valori_new;
}

struct valori eulerocromer(double dt, double omegaquadro, struct valori valori_n){
    struct valori valori_new;

    //la nuova velocità viene calcolata utilizzando a=-omega^2*x
    valori_new.v=valori_n.v - omegaquadro*valori_n.x*dt;

    //la nuova posizione viene calcolata utilizzando la velocità al passo n+1
    valori_new.x=valori_n.x+valori_new.v*dt;
    
    return valori_new;
}

struct valori puntocentrale(double dt, double omegaquadro, struct valori valori_n){
    struct valori valori_new;

    //la nuova velocità viene calcolata utilizzando a=-omega^2*x
    valori_new.v=valori_n.v - omegaquadro*valori_n.x*dt;

    //la nuova posizione viene calcolata utilizzando la media tra la velocità al passo n e quella al passo n+1
    valori_new.x=valori_n.x+(valori_new.v+valori_n.v)/2.*dt;
    
    return valori_new;
}

struct valori mezzopasso(double dt, double omegaquadro, struct valori valori_n){
    struct valori valori_new;

    //la nuova velocità viene calcolata utilizzando a=-omega^2*x
    valori_new.v=valori_n.v - omegaquadro*valori_n.x*dt;

    //la nuova posizione viene calcolata utilizzando la velocità al passo n+1
    valori_new.x=valori_n.x+valori_new.v*dt;
    
    return valori_new;
}

double energia(double m,double v,double k,double x){
    double e_potenziale,e_cinetica;
    e_potenziale=1./2. * k * x * x;
    e_cinetica=1./2. * m * v * v;
    return e_potenziale+e_cinetica;
}

int sceglialgoritmo(){
    int algoritmo;

    printf("Selezionare l\'algoritmo da utilizzare:\n0 - Eulero\n1 - Eulero-Cromer\n2 - Punto centrale\n3 - Mezzo passo\n");
    scanf(" %d",&algoritmo);
    return algoritmo;
}