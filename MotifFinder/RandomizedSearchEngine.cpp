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
    int maxPos = dna.Size(sequence) - 1;
    int i = rand();
    double random = (i / (double)(RAND_MAX)) * maxPos;
    return (int)floor(random);
}

void RandomizedSearchEngine::SetRepo(IDnaRepository& input){
    dna = input;
    scoreEngine.SetRepo(input);
}

vector<int> RandomizedSearchEngine::randomLoci() {
    vector<int> loci (dna.Size());
    for(int i = 0; i < dna.Size(); i++)
    {
        loci.push_back(i);
    }
    return loci;
}

vector<Nucleotide_t> RandomizedSearchEngine::lociToMotif(vector<int> loci, int motifLength, int dontCares) {
    //vector<vector<int> > 
}



void RandomizedSearchEngine::Search(int time, int motifLength, int dontCares) {
    
}


RandomizedSearchEngine::~RandomizedSearchEngine() {
    delete &scoreEngine;
    delete &bestMotif;
    delete &startingLoci;
    delete &profileMatrix;
}

