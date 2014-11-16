/* 
 * File:   RandomizedSearchEngine.cpp
 * Author: Eric
 * 
 * Created on November 15, 2014, 1:43 PM
 */
#include <stdlib.h>
#include <tgmath.h>
#include <time.h>
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
    vector<vector<int> > profileMatrix(4);
    for(int i = 0; i < 4; i++)
    {
        profileMatrix[i](motifLength);
    }
}



void RandomizedSearchEngine::Search(double time, int motifLength, int dontCares) {
    time_t epoch = time(NULL);
    vector<int> currentLoci;
    vector<Nucleotide_t> currentMotif;
    double bestScore = -10, currentScore;
    do{
        currentLoci = randomLoci();
        currentMotif = lociToMotif(currentLoci, motifLength, dontCares);
        currentScore = scoreEngine.Score(&currentMotif, &currentLoci);
        if(currentScore > bestScore) {
            startingLoci = currentLoci;
            motif = currentMotif;
        }
    }while(difftime(epoch, time(NULL) < time));
}

std::vector<Nucleotide_t> RandomizedSearchEngine::GetMotif() {
    return motif;
}

std::vector<int> RandomizedSearchEngine::GetStartingLoci() {
    return startingLoci;
}




RandomizedSearchEngine::~RandomizedSearchEngine() {
    delete &scoreEngine;
    delete &motif;
    delete &startingLoci;
    delete &profileMatrix;
}

