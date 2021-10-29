#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define g 9.81

struct valori{
    double alfa;
    double omega;
};

double phi(double,double);
double energia(double m,double v,double l,double x);
struct valori rungekutta(double dt, double coefficiente, struct valori valori_n);

int main(int argc, char* argv[]){

    if(argc!=7){
        fprintf(stderr,"Per l\'esecuzione del programma è necessario passare come argomenti: aplha0, omega0, dt, T, l, m.\n Sarà utilizzato l'algoritmo di Runge Kutta.\n");
        //\nAlgoritmo è un numero intero, scegliere tra:\n0 Eulero\n1 Eulero-Cromer\n2 Punto centrale\n3 Mezzo passo\n4 Verlet\n5 Verlet autosufficiente\n6 Predizione Correzione\n7 Runge Kutta\n
        exit(1);
    }

    struct valori valori_n; //posizione, velocità
    int  i;
    double tempo=0, T, dt, l, m ,npassi,energia_n,energia_0,energia_rapporto=0.,coefficiente;
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

double energia(double m,double omega,double l,double alfa){
    double e_potenziale,e_cinetica;
    e_potenziale=m*g*l*(1-cos(alfa));
    e_cinetica=(m*l*l)/2.*omega*omega;
    return e_potenziale+e_cinetica;
}