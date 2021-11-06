#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define ALGORITMI_NUM 9

struct valori{
    double x;
    double v;
};

double phi(double,double);

struct valori rungekutta(double dt, double omegaquadro, struct valori valori_n);
//void trovaperiodo()

double energia(double m,double v,double k,double x);

int main(int argc, char* argv[]) {

    if(argc!=7){
        //Il programma eseguira tutti gli algoritmi e salverà i dati in %d file del tipo Eulero.dat.
        fprintf(stderr,"Per l\'esecuzione del programma è necessario passare come argomenti: x0, v0, dt, T, k, m.\nAlgoritmo è Runge Kutta.\n");
        exit(1);
    }

    struct valori valori_n; //posizione, velocità
    int  i;
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
    algoritmo=atof(argv[7]);
    npassi=T/dt;
    omegaquadro=k/m;

    //uso malloc per assegnare le giuste celle di memoria agli array x,v e t
    x=(double*)malloc(sizeof(double)*npassi);
    v=(double*)malloc(sizeof(double)*npassi);
    t=(double*)malloc(sizeof(double)*npassi);
    if(t==NULL || x==NULL || v==NULL){
        printf("Errore nella creazione degli array dinamici.\n");
        exit(EXIT_FAILURE);
    }

    //calcolo E(0)
    energia_0=energia(m,valori_n.v,k,valori_n.x);

    fptr=fopen(nomefile,"w+");
    fprintf(fptr,"#t         x         v         E        delta_E/E(0)\n");
    fprintf(fptr,"%g %g %g %g %g\n",tempo,valori_n.x,valori_n.v,energia_0,energia_rapporto);
    x[0]=valori_n.x;
    v[0]=valori_n.v;
    t[0]=0.;

    for(i=1;i<=npassi;i++){
        //calcolo i valori con la funzione prelevata dall'array algoritmi_lista
        valori_n=rungekutta(dt,omegaquadro,valori_n);
        energia_n=energia(m,valori_n.v,k,valori_n.x); //calcolo energia utilizzando la rispettiva funzione energia
        energia_rapporto=(energia_n-energia_0)/energia_0;
        tempo=(double)i*dt; //calcolo il tempo così per non perdere precisione per via dell'approssimazione di numeri molto piccoli
        fprintf(fptr,"%g %g %g %g %g\n",tempo,valori_n.x,valori_n.v,energia_n,energia_rapporto); //stampo su un file t, x, v e energia
        x[i]=valori_n.x;
        v[i]=valori_n.v;
        t[i]=tempo;
    }
    fclose(fptr);
    return(0);
}

//funzione che calcola l'accelerazione
double phi(double omegaquadro,double x,double gamma, double v){
    double phi;
    phi= -omegaquadro*x-gamma*v;
    return phi;
}

struct valori eulero(double dt, double omegaquadro, struct valori valori_n){
    struct valori valori_new;
    //la nuova posizione viene calcolata utilizzando la velocità al passo n
    valori_new.x=valori_n.x+valori_n.v*dt;

    //la nuova velocità viene calcolata utilizzando a=-omega^2*x
    valori_new.v=valori_n.v + phi(omegaquadro,valori_n.x)*dt;

    return valori_new;
}

struct valori eulerocromer(double dt, double omegaquadro, struct valori valori_n){
    struct valori valori_new;

    //la nuova velocità viene calcolata utilizzando a=-omega^2*x
    valori_new.v=valori_n.v + phi(omegaquadro,valori_n.x)*dt;

    //la nuova posizione viene calcolata utilizzando la velocità al passo n+1
    valori_new.x=valori_n.x+valori_new.v*dt;
    
    return valori_new;
}

struct valori puntocentrale(double dt, double omegaquadro, struct valori valori_n){
    struct valori valori_new;

    //la nuova velocità viene calcolata utilizzando a=-omega^2*x
    valori_new.v=valori_n.v + phi(omegaquadro,valori_n.x)*dt;

    //la nuova posizione viene calcolata utilizzando la media tra la velocità al passo n e quella al passo n+1
    valori_new.x=valori_n.x+(valori_new.v+valori_n.v)/2.*dt;
    
    return valori_new;
}

struct valori mezzopasso(double dt, double omegaquadro, struct valori valori_n){
    struct valori valori_new;

    //la nuova velocità viene calcolata utilizzando a=-omega^2*x
    valori_new.v=valori_n.v + phi(omegaquadro,valori_n.x)*dt;

    //la nuova posizione viene calcolata utilizzando la velocità al passo n+1
    valori_new.x=valori_n.x+valori_new.v*dt;
    
    return valori_new;
}

struct valori verlet(double dt, double omegaquadro, struct valori valori_n, double x_old){
    struct valori valori_new;
    double x_futuro;

    //ho bisogno di x_old perchè nella formula compare x(n-2)
    valori_new.x=2*valori_n.x-x_old-omegaquadro*valori_n.x*dt*dt;
    //per calcolare la velocità mi serve x(n+2)
    x_futuro=2*valori_new.x-valori_n.x-omegaquadro*valori_new.x*dt*dt;
    valori_new.v=(x_futuro-x_old)/(2*dt);
    return valori_new;
}

struct valori verletautosufficiente(double dt, double omegaquadro, struct valori valori_n){
    struct valori valori_new;

    valori_new.x=valori_n.x+valori_n.v*dt-(0.5)*omegaquadro*valori_n.x*dt*dt;
    valori_new.v=valori_n.v+((-omegaquadro*valori_n.x-omegaquadro*valori_new.x)/2.)*dt;
    return valori_new;
}

struct valori predcorr(double dt, double omegaquadro, struct valori valori_n, double x_old){
    struct valori valori_new;
    double x_futuro,phi_new,phi_n;

    //ho bisogno di x_old per calcolare un x(n+1) approssimativo per calcolare poi phi(n+1)
    x_futuro=x_old+2*valori_n.v*dt;
    phi_new=phi(omegaquadro,x_futuro);
    phi_n=phi(omegaquadro,valori_n.x);
    
    valori_new.v=valori_n.v+((phi_new+phi_n)/2.)*dt;
    valori_new.x=valori_n.x+((valori_n.v+valori_new.v)/2.)*dt;
    return valori_new;
}

struct valori rungekutta(double dt, double omegaquadro, struct valori valori_n){
    struct valori valori_new;
    double dx,dv;

    dx=valori_n.v*dt;
    dv=phi(omegaquadro,valori_n.x)*dt;

    valori_new.x=valori_n.x+(valori_n.v+0.5*dv)*dt;

    valori_new.v=valori_n.v + phi(omegaquadro,valori_n.x+0.5*dx)*dt;
    
    return valori_new;
}

struct valori rungekutta4(double dt, double omegaquadro, struct valori valori_n){
    //DA CORREGGERE
    struct valori valori_new;
    double dx1,dx2,dx3,dx4,dv1,dv2,dv3,dv4;

    dx1=valori_n.v*dt;
    dv1=phi(omegaquadro,valori_n.x)*dt;
    dx2=(valori_n.v+0.5*dv1)*dt;
    dv2=phi(omegaquadro,(valori_n.x+0.5*dx1))*dt;
    dx3=(valori_n.v+0.5*dv2)*dt;
    dv3=phi(omegaquadro,(valori_n.x+0.5*dx2))*dt;
    dx4=(valori_n.v+0.5*dv3)*dt;
    dv4=phi(omegaquadro,(valori_n.x+0.5*dx3))*dt;

    valori_new.x=valori_n.x+(valori_n.v+(1/6)*(dv1+2*dv2+2*dv3+dv4))*dt;
    valori_new.v=valori_n.v + phi(omegaquadro,valori_n.x+(1/6)*(dx1+2*dx2+2*dx3+dx4))*dt;
    
    return valori_new;
}

double energia(double m,double v,double k,double x){
    double e_potenziale,e_cinetica;
    e_potenziale=1./2. * k * x * x;
    e_cinetica=1./2. * m * v * v;
    return e_potenziale+e_cinetica;
}