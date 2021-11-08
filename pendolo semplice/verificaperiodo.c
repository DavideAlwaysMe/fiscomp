#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#define g 9.81

struct valori{
    double alfa;
    double omega;
};

double phi(double,double);
double energia(double m,double v,double l,double x);
struct valori rungekutta(double dt, double coefficiente, struct valori valori_n);
double trovaperiodo(double omega[],double t[],int npassi);
double interpolazionelin(double y1,double y2,double x1, double x2,double xstar);

int main(int argc, char* argv[]){

    if(argc!=6){
        fprintf(stderr,"Per l\'esecuzione del programma è necessario passare come argomenti: omega0, dt, T, l, errore_percentuale_massimo .\n Sarà utilizzato l'algoritmo di Runge Kutta.\n");
        exit(1);
    }

    struct valori valori_n; //posizione, velocità
    int  i,count;
    double T, dt, l ,npassi,coefficiente,periodo,periodoteorico,diff,diff_perc,delta_alpha=0.01,err_max;
    double alfa, *omega, *t;
    //assegno i parametri di esecuzione alle variabili iniziali
    //alfa_max=atof(argv[1]);
    valori_n.omega=atof(argv[1]);
    dt=atof(argv[2]);
    T=atof(argv[3]);
    l=atof(argv[4]);
    err_max=atof(argv[5]);
    npassi=T/dt;
    coefficiente=g/l;

    //uso malloc per assegnare le giuste celle di memoria agli array x,v e t
    omega =(double*)malloc(sizeof(double)*npassi);
    t=(double*)malloc(sizeof(double)*npassi);
    if(t==NULL || omega==NULL){
        printf("Errore nella creazione degli array dinamici.\n");
        exit(EXIT_FAILURE);
    }

    omega[0]=valori_n.omega;
    t[0]=0.;

    while(diff_perc<err_max){
        alfa=delta_alpha*count;
        valori_n.alfa=alfa;
        for(i=1;i<=npassi;i++){
            valori_n=rungekutta(dt,coefficiente,valori_n);
            t[i]=(double)i*dt; //calcolo il tempo così per non perdere precisione per via dell'approssimazione di numeri molto piccoli
            omega[i]=valori_n.omega; //mi salvo solo omega e t perchè la funzione trovaperiodo cerca gli zeri di omega per trovare il periodo
        }
        
        periodo=trovaperiodo(omega,t,npassi);
        periodoteorico=2*M_PI*sqrt(l/g);
        diff=periodo-periodoteorico;
        diff_perc=(diff)/periodoteorico*100;
        printf("alfa: %lf differenza percentuale: %lf\n",alfa,diff_perc);
        count++;
    }
}

//funzione che calcola l'accelerazione
double phi(double coefficiente,double alfa){
    double phi;
    phi= -coefficiente*sin(alfa);
    return phi;
}

struct valori rungekutta(double dt, double coefficiente, struct valori valori_n){
    struct valori valori_new;
    double dalfa,domega;

    dalfa=valori_n.omega*dt;
    domega=phi(coefficiente,valori_n.alfa)*dt;

    valori_new.alfa=valori_n.alfa+(valori_n.omega+0.5*domega)*dt;

    valori_new.omega=valori_n.omega + phi(coefficiente,valori_n.alfa+0.5*dalfa)*dt;
    
    return valori_new;
}

double trovaperiodo(double omega[],double t[],int npassi){
    int i,zericounter=-1; //zericounter sarà il numero di zeri (quindi semiperiodi) trovati
    double sumsemiperiodi=0,periodo,ultimo_tempo=0.,tempo;

    for(i=0;i<npassi;i++){
        if(omega[i]*omega[i+1]<0){
            //c'è stato un cambio di segno nella velocità
            //devo interpolare per stimare il tempo in cui la velocità ha cambiato di segno
            tempo=interpolazionelin(t[i],t[i+1],omega[i],omega[i+1],0.);
            //considero il tempo trovato come un semiperiodo trovato solo se non è il primo 
            if(zericounter!=-1){
                sumsemiperiodi+=tempo-ultimo_tempo;
            }
            zericounter++;
            //la variabile ultimo_tempo contiene l'ultimo tempo in cui è stato misurato un cambio di direzione di omega
            ultimo_tempo=tempo;
        }
    }
    //calcolo la media dei semiperiodi e moltiplico per due
    periodo=(sumsemiperiodi/zericounter)*2;
    return periodo;
}

double interpolazionelin(double y1,double y2,double x1, double x2,double xstar){
    double a,b,ystar;
    a = y1-(y1-y2)/(x1-x2)*x1;
    b = (y1-y2)/(x1-x2);
    ystar= a + b * xstar;
    return ystar;
}