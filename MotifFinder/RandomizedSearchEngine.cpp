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
    int maxPos = dna.Size(sequence) - motifLength;
    int i = rand();
    double random = (i / (double)(RAND_MAX)) * (maxPos);
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
        loci.push_back(randomPosition(i, motifLength));
    }
    return loci;
}

vector<Nucleotide_t> RandomizedSearchEngine::lociToMotif(vector<int> loci) {
    vector<Nucleotide_t> bestMotif = vector<Nucleotide_t>();
    double bestMotifScore = -1000;
    vector<Nucleotide_t> startMotif = vector<Nucleotide_t>();
    vector<Nucleotide_t> currentMotif;
    profileMatrix = createProfileMatrix(loci); //get a motif from matrix with the loci
    startMotif = getStartingMotif(profileMatrix, startMotif, motifLength);
    currentMotif = new(startMotif);         //copies over contents, not just reference
    for (int j = 0; j < dontCares; j++){    //create initial starting motif with dont cares
        currentMotif.erase(j+1);
        currentMotif.insert(j+1, 4);
    }     //try all valid positions of the motif with don't cares in valid locations
    while (canIncrement(currentMotif)){ //change the location of the dont care
        currentMotif = increment(currentMotif, startMotif);
        //if current motif is better score than best motif, then best motif = current
        double currentMotifScore = scoreEngine.Score(currentMotif, loci);
        if (bestMotifScore < currentMotifScore){
            bestMotifScore = currentMotifScore;
            bestMotif = new(currentMotif);
        }
    } //return the best scoring motif w/ don't cares placed
    return bestMotif;
}

vector<Nucleotide_t> RandomizedSearchEngine::increment(vector<Nucleotide_t> motif, vector<Nucleotide_t> startMotif){
    //find dc that isn't all the way right. Move it right one. Then move all 
    //dcs that are to the right, next to what we just incremented
    vector<int> locations = vector<int>();
    bool moved = false;
    int vectorSize = motif.size();
    int location = vectorSize - 3;
    while (!moved){
        if (location != 0 && motif.at(location) == 4){//not at end and is a dc
            if (motif.at(location + 1) != 4){
                moved = true;
                motif.erase(location);
                motif.insert(location, startMotif.at(location));
            } else {
                locations.push_back(location);
            }
        }
        location--;
    } //iterated the dc one spot forward. move all to the right of loc leftmost
    int size = locations.size();
    for (int i = 0; i < size; i++){
        int whereInMotif = locations.pop_back();
        motif.erase(whereInMotif);
        motif.insert(whereInMotif, startMotif.at(whereInMotif));
    }
    return motif;
}
bool RandomizedSearchEngine::canIncrement(vector<Nucleotide_t> motif){
    //start at back of vector and check if all dcs (4s) are at the back 
    int vectorSize = motif.size();
    for (int i = 1; i <= dontCares; i++){
        if (motif.at(vectorSize - (1 + i)) != 4){
            return false;
        }
    }
    return true;
}

void RandomizedSearchEngine::Search(double runTime) {
    time_t epoch;
    epoch = time(NULL);
    vector<int> currentLoci;
    vector<Nucleotide_t> currentMotif;
    double bestScore = -10, currentScore;
    do{
        currentLoci = randomLoci();
        currentMotif = lociToMotif(currentLoci);
        currentScore = scoreEngine.Score(currentMotif, currentLoci);
        if(currentScore > bestScore) {
            startingLoci = currentLoci;
            motif = currentMotif;
        }
    }while(difftime(epoch, time(NULL) < runTime));
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
    delete &dontCares;
    delete &motifLength;
}

RandomizedSearchEngine::RandomizedSearchEngine(IDnaRepository* input, int motiflength, 
        int dontcares){
    dontCares = dontcares;
    motifLength = motiflength;
    scoreEngine = ScoreEngine();
    //create profile matrix and fill with all 0s of 4 columns and motiflength rows
    profileMatrix.resize(4, vector<int>(motifLength, 0));
    dna = input;
}

int RandomizedSearchEngine::getMax(int a, int t, int g, int c){
    max = a;
    if (max > t)   
        max = t;
    if (max < g)
        max = g;
    if (max < c)
        max = c;
    return max;
}
vector<Nucleotide_t> RandomizedSearchEngine::getStartingMotif(vector<vector<int> > profileMatrix, 
        vector<Nucleotide_t> motif, int length){
    for (int k = 0; k < length; k++) { //get most probable (0, 1, 2, 3) at k and add to startMotif
        int a = profileMatrix[0][k];
        int t = profileMatrix[1][k];
        int g = profileMatrix[2][k];
        int c = profileMatrix[3][k];
        int max = getMax(a, t, g, c);
        if (max == a)  
            motif.push_back(0);
        else if (max == t) 
            motif.push_back(1);
        else if (max == g) 
            motif.push_back(2);
        else    
            motif.push_back(3);
    }
    return motif;
}

vector<vector<int> > RandomizedSearchEngine :: createProfileMatrix(vector<int> loci){
    for (int i = 0; i < motifLength; i++){
        int dnaSize = dna.Size();
        vector<int> counts(4);
        counts = {0, 0, 0, 0};
        for (int j = 0; j < dnaSize; j++){
            int nucleotide = dna.Get(j, loci[j] + i);
            counts[nucleotide] += 1;
        }
        profileMatrix[0][i] = counts.at(0);
        profileMatrix[1][i] = counts.at(1);
        profileMatrix[2][i] = counts.at(2);
        profileMatrix[3][i] = counts.at(3);
    }
    return profileMatrix;
}