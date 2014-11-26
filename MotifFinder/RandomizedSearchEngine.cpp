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

RandomizedSearchEngine::RandomizedSearchEngine(IDnaRepository& input, int motiflength, int dontcares) : dna(input), scoreEngine(ScoreEngine(input)){
    dontCares = dontcares;
    motifLength = motiflength;
    score = -10;
    //create profile matrix and fill with all 0s of 4 columns and motiflength rows
    profileMatrix.resize(4, vector<int>(motifLength, 0));
}

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
    vector<int> loci = vector<int>();
    for(int i = 0; i < dna.Size(); i++)
    {
        loci.push_back(randomPosition(i));
    }
    return loci;
}

vector<Nucleotide_t> RandomizedSearchEngine::lociToMotif(vector<int> loci) {
    vector<Nucleotide_t> bestMotif = vector<Nucleotide_t>();
    double bestMotifScore = -1000;
    vector<Nucleotide_t> startMotif;
    vector<Nucleotide_t> currentMotif;
    profileMatrix = createProfileMatrix(loci); //get a motif from matrix with the loci
    startMotif = getStartingMotif(profileMatrix);
    currentMotif = startMotif;
    for (int i = 0; i < dontCares; i++){
        currentMotif[i + 1] = DC;
    }//copies over contents, not just reference
    bestMotif = currentMotif;
    if (dontCares == 0){
        return currentMotif;
    }
    for (int j = 0; j < dontCares; j++){    //create initial starting motif with dont cares
        currentMotif[j+1] = DC;
    }     //try all valid positions of the motif with don't cares in valid locations
    while (canIncrement(currentMotif)){ //change the location of the dont care
        currentMotif = increment(currentMotif, startMotif);
        //if current motif is better score than best motif, then best motif = current
        double currentMotifScore = scoreEngine.Score(currentMotif, loci);
        if (bestMotifScore < currentMotifScore){
            bestMotifScore = currentMotifScore;
            bestMotif = currentMotif;
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
        if (location != 0 && motif[location] == 4){
            if(motif[location + 1] != 4)
            {
                motif[location + 1] = DC;
                motif[location] = startMotif[location];
                moved = true;
            } else {
                locations.push_back(location);
                motif[location] = startMotif[location];
            }
        }
        location--;
    } //iterated the dc one spot forward. move all to the right of loc leftmost
    location += 3;
    int size = locations.size();
    for (int i = 0; i < size; i++){
        motif[location + i] = DC;
    }
    return motif;
}
bool RandomizedSearchEngine::canIncrement(vector<Nucleotide_t> motif){
    //start at back of vector and check if all dcs (4s) are at the back 
    int vectorSize = motif.size();
    for (int i = 0; i < dontCares; i++){
        if (motif[vectorSize - (2 + i)] != 4){
            return true;
        }
    }
    return false;
}

void RandomizedSearchEngine::Search(double runTime) {
    time_t epoch;
    epoch = time(NULL);
    vector<int> currentLoci;
    vector<Nucleotide_t> currentMotif;
    double currentScore;
    do{
        currentLoci = randomLoci();
        currentMotif = lociToMotif(currentLoci);
        currentScore = scoreEngine.Score(currentMotif, currentLoci);
        if(currentScore > score) {
            startingLoci = currentLoci;
            motif = currentMotif;
            score = currentScore;
        }
    }while(difftime(time(NULL), epoch) < runTime);
}

double RandomizedSearchEngine::GetBestScore() {
    return score;
}


std::vector<Nucleotide_t> RandomizedSearchEngine::GetMotif() {
    return motif;
}

std::vector<int> RandomizedSearchEngine::GetStartingLoci() {
    return startingLoci;
}

Nucleotide_t RandomizedSearchEngine::getMax(int a, int t, int g, int c){
    int max = a, index = 0;
    if (max < t){   
        max = t;
        index = 1;
    }
    if (max < g){
        max = g;
        index = 2;
    }
    if (max < c){
        max = c;
        index = 3;
    }
    return (Nucleotide_t)index;
}
vector<Nucleotide_t> RandomizedSearchEngine::getStartingMotif(vector<vector<int> > profileMatrix){
    vector<Nucleotide_t> motif(profileMatrix[0].size(),A);
    for (int k = 0; k < motifLength; k++) { //get most probable (0, 1, 2, 3) at k and add to startMotif
        int a = profileMatrix[0][k];
        int t = profileMatrix[1][k];
        int g = profileMatrix[2][k];
        int c = profileMatrix[3][k];
        Nucleotide_t maxNucl = getMax(a, t, g, c);
        motif[k] = maxNucl;
    }
    return motif;
}

vector<vector<int> > RandomizedSearchEngine :: createProfileMatrix(vector<int> loci){
    for (int i = 0; i < motifLength; i++){
        int dnaSize = dna.Size();
        vector<int> counts(4, 0);
        for (int j = 0; j < dnaSize; j++){
            int nucleotide = dna.Get(j, loci[j] + i);
            counts[nucleotide] += 1;
        }
        profileMatrix[0][i] = counts[0];
        profileMatrix[1][i] = counts[1];
        profileMatrix[2][i] = counts[2];
        profileMatrix[3][i] = counts[3];
    }
    return profileMatrix;
}

int RandomizedSearchEngine::GetDontCareCount() {
    return dontCares;
}

int RandomizedSearchEngine::GetMotifLength() {
    return motifLength;
}

void RandomizedSearchEngine::SetDontCares(int number) {
    dontCares = number;
}

void RandomizedSearchEngine::SetMotifLength(int length) {
    motifLength = length;
}


RandomizedSearchEngine::~RandomizedSearchEngine() {
}