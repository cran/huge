#include<R.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void SFGen(int *dd0, int*dd, int *G, int *seed){
    int i,j;
    int d0,d;
    d0 = dd0[0];
    d = dd[0];
    double x;
    int *size_a = (int*) malloc((d)*sizeof(int));
    int tmp;
    int total;

    srand(*seed);

    for(i=0;i<(d0-1);i++){
        G[i*d+i+1] = 1;
        G[(i+1)*d+i] = 1;
    }
    G[d0-1] = 1;
    G[(d0-1)*d] = 1;
        
    for(i=0;i<d0;i++){
        size_a[i] = 2;
    }
    for(i=d0;i<d;i++){
        size_a[i] = 0;
    }
    total = 2*d0;
    
    for(i=d0;i<d;i++){
        x = (double) total*rand()/RAND_MAX;
        tmp = 0;
        j = 0;
        while(tmp<x&&j<i){
            tmp += size_a[j];
            j++;
        }
        j = j-1;
        G[i*d+j] = 1;
        G[j*d+i] = 1;
        total+=2;
        size_a[j]++;
        size_a[i]++;
    }
    free(size_a);
}
