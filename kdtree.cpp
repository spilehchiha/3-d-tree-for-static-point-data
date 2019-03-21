//
//  main.c
//  Nearest Neighbour by Scanning
//
//  Created by Sina Pilehchiha on 2019-03-21.
//  Copyright © 2019 TeamNumber02. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include<set>
// Algortihm: Nearest Neighbour by Scanning (Naïve Nearest Neighbour)

typedef struct{
    double vector[3];
} domainVector;
typedef struct{
    double vector[3];
} rangeVector;

typedef struct {
    domainVector exemplarDomainVector;
    rangeVector exemplarRangeVector;
} exemplar;

/*
 int main(int argc, const char * argv[]) {
 exemplar exlist[4] =  {{-1, -1, -1}, {1, -1, 5}, {6, 7, 8}, {9, 0, -13}};
 domainVector dom   =  {2, 3, 3};
 exemplar nearest;
 double nearestDist = INFINITY;
 exemplar ex;
 
 for (int i = 0; i < sizeof(exlist) / sizeof(exemplar); ++i) {
 ex.exemplarDomainVector.vector[0] = exlist[i].exemplarDomainVector.vector[0];
 ex.exemplarDomainVector.vector[1] = exlist[i].exemplarDomainVector.vector[1];
 ex.exemplarDomainVector.vector[2] = exlist[i].exemplarDomainVector.vector[2];
 
 double dist = sqrt(pow((dom.vector[0] - ex.exemplarDomainVector.vector[0]), 2) +
 pow((dom.vector[1] - ex.exemplarDomainVector.vector[1]), 2) +
 pow((dom.vector[2] - ex.exemplarDomainVector.vector[2]), 2));
 
 if (dist < nearestDist) {
 nearestDist = dist;
 nearest = ex;
 }
 }
 return 0;
 }
 */
// Algorithm: Constructing a k_d-tree
typedef struct {
    std::set<exemplar> exemplar_set;
} exemplarSet;

typedef struct {
    struct node {
        domainVector domElt;
        rangeVector  rangeElt;
        int          split;
        struct node  *left;
        struct node  *right;
    };
} kdtree;

int main(int argc, const char * argv[]) {
    return 0;
}
