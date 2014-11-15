/* 
 * File:   RandomizedSearchEngine.cpp
 * Author: Eric
 * 
 * Created on November 15, 2014, 1:43 PM
 */
#include <stdlib.h>
#include <tgmath.h>
#include <stdexcept>
#include "RandomizedSearchEngine.h"


int RandomizedSearchEngine::randomPosition(int sequence){
    int maxPos = dna->Size(sequence) - 1;
    int i = rand();
    double random = (i / (double)(RAND_MAX)) * maxPos;
    return (int)floor(random);
}

void RandomizedSearchEngine::SetRepo(IDnaRepository* input){
    if(input == 0)
        throw std::invalid_argument("input cannot be 0.");
    dna = input;
    scoreEngine.SetRepo(input);
}

RandomizedSearchEngine::~RandomizedSearchEngine() {
    delete &scoreEngine;
    delete &bestMotif;
    delete &startingLoci;
    delete &profileMatrix;
}

