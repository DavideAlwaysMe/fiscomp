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
double calcmedia(double array[],double dim);

int main(int argc, char* argv[]){

    if(argc!=7){
        fprintf(stderr,"Per l\'esecuzione del programma è necessario passare come argomenti: aplha0, omega0, dt, T, l, m.\n Sarà utilizzato l'algoritmo di Runge Kutta.\n");
        //\nAlgoritmo è un numero intero, scegliere tra:\n0 Eulero\n1 Eulero-Cromer\n2 Punto centrale\n3 Mezzo passo\n4 Verlet\n5 Verlet autosufficiente\n6 Predizione Correzione\n7 Runge Kutta\n
        exit(1);
    }

    struct valori valori_n; //posizione, velocità
    int  i;
    double tempo=0, T, dt, l, m ,npassi,energia_n,energia_0,energia_rapporto=0.,coefficiente,periodo,periodoteorico;
    FILE *fptr;
    double *alfa, *omega, *t;
    //assegno i parametri di esecuzione alle variabili iniziali
    valori_n.alfa=atof(argv[1]);
    valori_n.omega=atof(argv[2]);
    dt=atof(argv[3]);
    T=atof(argv[4]);
    l=atof(argv[5]);
    m=atof(argv[6]);
    npassi=T/dt;
    coefficiente=g/l;

    //uso malloc per assegnare le giuste celle di memoria agli array x,v e t
    alfa =(double*)malloc(sizeof(double)*npassi);
    omega =(double*)malloc(sizeof(double)*npassi);
    t=(double*)malloc(sizeof(double)*npassi);
    if(t==NULL || alfa==NULL || omega==NULL){
        printf("Errore nella creazione degli array dinamici.\n");
        exit(EXIT_FAILURE);
    }

    //calcolo E(0)
    energia_0=energia(m,valori_n.omega,l,valori_n.alfa);

    //calcolo i valori di x, v e t con l'algoritmo scelto e salvo in un file di testo
    fptr=fopen("pendolosemplice.dat","w+");
    fprintf(fptr,"#t         x         v         E        delta_E/E(0)\n");
    fprintf(fptr,"%g %g %g %g %g\n",tempo,valori_n.alfa,valori_n.omega,energia_0,energia_rapporto);
    alfa[0]=valori_n.alfa;
    omega[0]=valori_n.omega;
    t[0]=0.;

    for(i=1;i<=npassi;i++){
        //calcolo i valori con la funzione prelevata dall'array algoritmi_lista
        valori_n=rungekutta(dt,coefficiente,valori_n);
        energia_n=energia(m,valori_n.omega,l,valori_n.alfa); //calcolo energia utilizzando la rispettiva funzione energia
        energia_rapporto=(energia_n-energia_0)/energia_0;
        tempo=(double)i*dt; //calcolo il tempo così per non perdere precisione per via dell'approssimazione di numeri molto piccoli
        fprintf(fptr,"%g %g %g %g %g\n",tempo,valori_n.alfa,valori_n.omega,energia_n,energia_rapporto); //stampo su un file t, x, v e energia
        alfa[i]=valori_n.alfa;
        omega[i]=valori_n.omega;
        t[i]=tempo;
    }
    
    periodo=trovaperiodo(omega,t,npassi);
    periodoteorico=2*M_PI*sqrt(l/g);
    printf("%lf",periodo-periodoteorico);
}

//funzione che calcola l'accelerazione
double phi(double coefficiente,double alfa){
    double phi;
    phi= -coefficiente*alfa;
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
    int i,zericounter=0; //zericounter sarà il numero di zeri (quindi semiperiodi) trovati
    double sumsemiperiodi,periodo,ultimo_tempo=0.,tempo; //da ritoccare la dimensione di semiperiodo

    for(i=0;i<npassi;i++){
        if(omega[i]*omega[i+1]<0){
            //c'è stato un cambio di segno nella velocità
            //devo interpolare per stimare il tempo in cui la velocità ha cambiato di segno
            tempo=interpolazionelin(t[i],t[i+1],omega[i],omega[i+1],0.);
            //considero il semiperiodo trovato solo se non è il primo 
            if(zericounter!=0){
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

double energia(double m,double omega,double l,double alfa){
    double e_potenziale,e_cinetica;
    e_potenziale=m*g*l*(1-cos(alfa));
    e_cinetica=(m*l*l)/2.*omega*omega;
    return e_potenziale+e_cinetica;
}

double interpolazionelin(double y1,double y2,double x1, double x2,double xstar){
    double a,b,ystar;
    a = y1-(y1-y2)/(x1-x2)*x1;
    b = (y1-y2)/(x1-x2);
    ystar= a + b * xstar;
    return ystar;
}

double calcmedia(double array[],double dim){
    int i;
    double sum;
    for(i=0;i<dim;i++){
        sum+=array[i];
    }
    return sum/dim;
}