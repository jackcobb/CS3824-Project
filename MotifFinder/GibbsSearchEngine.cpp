/* 
 * File:   GibbsSearchEngine.cpp
 * Author: Eric
 * 
 * Created on November 26, 2014, 2:06 PM
 */
#include <iostream>
#include <vector>
#include <time.h>
#include <stdlib.h>
#include <tgmath.h>
#include "GibbsSearchEngine.h"
#include "RandomEnumerator.h"
#include "ScoreEngine.h"

using namespace std;

GibbsSearchEngine::GibbsSearchEngine(IDnaRepository& repo, int motifLength, int dontCares) : repo(repo), scoreEngine(ScoreEngine(repo)){
    this->motifLength = motifLength;
    this->dontCares = dontCares;
    profileMatrixFull = vector<vector<double> >(4, vector<double>(motifLength, 0));
    profileMatrixGibbs = vector<vector<double> >(4, vector<double>(motifLength, 0));
}

void GibbsSearchEngine::Search(double runTime) {
    time_t epoch = time(NULL);
    SearchResult currentResult, bestResult;
    bestResult = Search(runTime, 140);
    double timeRemaining = difftime(runTime, difftime(time(NULL), epoch));
    while(timeRemaining > 0)
    {
        currentResult = Search(timeRemaining, 140);
        if(currentResult.score > bestResult.score)
        {
            bestResult.loci = currentResult.loci;
            bestResult.motif = currentResult.motif;
            bestResult.score = currentResult.score;
        }
        timeRemaining = difftime(runTime, difftime(time(NULL), epoch));
    }    
    
    loci = bestResult.loci;
    motif = bestResult.motif;
    score = bestResult.score;
}

SearchResult GibbsSearchEngine::Search(double maxRunTime, int maxIterWithoutImprovement) {
    vector<int> startingLoci = getRandomLoci();
    return Search(maxRunTime, maxIterWithoutImprovement, startingLoci);
}

SearchResult GibbsSearchEngine::Search(double maxRunTime, int maxIterWithoutImprovement, vector<int>& startingLoci) {
    time_t epoch = time(NULL);
    vector<vector<double> > profileMatrix(4, vector<double>(motifLength, 0));
    SearchResult bestResult;
    SearchResult currentResult;
    currentResult.loci = startingLoci;
    currentResult.motif = lociToMotif(currentResult.loci, profileMatrix);
    currentResult.score = scoreEngine.Score(currentResult.motif, currentResult.loci);
    bestResult.loci = currentResult.loci;
    bestResult.motif = currentResult.motif;
    bestResult.score = currentResult.score;
    double oldScore = -10;
    int currentOptIndex = 0;
    RandomEnumerator enumerator;
    int iterWithoutImprovement = 0;
    do{
        oldScore = currentResult.score;
        currentOptIndex = getRandomSequence();
        updateProfileMatrixExcluding(currentOptIndex, currentResult.loci, profileMatrix);
        probDistri = probabilityLMer(currentOptIndex, profileMatrix);
        enumerator.InitializeDistribution(probDistri);
        currentResult.loci[currentOptIndex] = enumerator.EnumerateRandomVar();
        currentResult.motif = lociToMotif(currentResult.loci, profileMatrix);
        currentResult.score = scoreEngine.Score(currentResult.motif, currentResult.loci);
        if(currentResult.score > bestResult.score)
        {
            bestResult.loci = currentResult.loci;
            bestResult.motif = currentResult.motif;
            bestResult.score = currentResult.score;
            iterWithoutImprovement = 0;
        }
        else
        {
            iterWithoutImprovement++;
        }
    }while((difftime(time(NULL), epoch) < maxRunTime) && iterWithoutImprovement <= maxIterWithoutImprovement);
    return bestResult;
}


vector<int> GibbsSearchEngine::getRandomLoci() {
    vector<int> loci = vector<int>(repo.Size(), 0);
    for(int i = 0; i < repo.Size(); i++)
    {
        loci[i] = randomPosition(i);
    }
    return loci;
}

int GibbsSearchEngine::randomPosition(int i) {
    int maxPos = repo.Size(i) - motifLength + 1;
    int randI = rand();
    double random = (randI / (double)(RAND_MAX)) * (maxPos);
    return (int)floor(random);
}

int GibbsSearchEngine::getRandomSequence() {
    int randI = rand();
    double random = (randI / (double)(RAND_MAX)) * (repo.Size() - 1);
    return (int)floor(random);    
}

vector<Nucleotide_t> GibbsSearchEngine::lociToMotif(vector<int> loci, vector<vector<double> >& profileMatrix) {
    vector<Nucleotide_t> bestMotif;
    vector<Nucleotide_t> startMotif;
    updateProfileMatrix(loci, profileMatrix); //get a motif from matrix with the loci
    startMotif = getStartingMotif(profileMatrix);
    bestMotif = scoreEngine.optimizeDontCaresInMotif(startMotif, dontCares, loci);
    return bestMotif;
}

Nucleotide_t GibbsSearchEngine::getMax(int a, int t, int g, int c){
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

vector<Nucleotide_t> GibbsSearchEngine::getStartingMotif(vector<vector<double> >& profileMatrix){
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

void GibbsSearchEngine::updateProfileMatrixExcluding(int excludedIndex, vector<int> loci, vector<vector<double> >& profileToUpdate ) {
    int motifIterator;
    int sequenceIterator;
    int nucleotide;
    int dnaSize = repo.Size();
    int counts[4];
    for (motifIterator = 0; motifIterator < motifLength; motifIterator++){
        counts[0] = 0;
        counts[1] = 0;
        counts[2] = 0;
        counts[3] = 0;
        for (sequenceIterator = 0; sequenceIterator < excludedIndex; sequenceIterator++){
            nucleotide = repo.Get(sequenceIterator, loci[sequenceIterator] + motifIterator);
            counts[nucleotide] += 1;
        }
        for (sequenceIterator = excludedIndex + 1; sequenceIterator < dnaSize; sequenceIterator++)
        {
            nucleotide = repo.Get(sequenceIterator, loci[sequenceIterator] + motifIterator);
            counts[nucleotide] += 1;            
        }
        profileToUpdate[0][motifIterator] = counts[0];
        profileToUpdate[1][motifIterator] = counts[1];
        profileToUpdate[2][motifIterator] = counts[2];
        profileToUpdate[3][motifIterator] = counts[3];
    }
}

void GibbsSearchEngine::updateProfileMatrix(vector<int> loci, vector<vector<double> >& profileToUpdate) {
    int motifIterator;
    int sequenceIterator;
    int nucleotide;
    int dnaSize = repo.Size();
    vector<int> counts(4, 0);
    for (motifIterator = 0; motifIterator < motifLength; motifIterator++){
        counts[0] = 0;
        counts[1] = 0;
        counts[2] = 0;
        counts[3] = 0;
        for (sequenceIterator = 0; sequenceIterator < dnaSize; sequenceIterator++)
        {
            nucleotide = repo.Get(sequenceIterator, loci[sequenceIterator] + motifIterator);
            counts[nucleotide] += 1;            
        }
        profileToUpdate[0][motifIterator] = counts[0];
        profileToUpdate[1][motifIterator] = counts[1];
        profileToUpdate[2][motifIterator] = counts[2];
        profileToUpdate[3][motifIterator] = counts[3];
    }
}

vector<double> GibbsSearchEngine::probabilityLMer(int sequence, vector<vector<double> >& profileMatrix) {
    int lastIndex = repo.Size(sequence) - motifLength + 1;
    vector<double> probLmers(lastIndex, 0);
    double runningSum = 0;
    for(int i = 0; i < lastIndex; i++)
    {
        probLmers[i] = probOfLmerAt(sequence, i, profileMatrix);
        runningSum += probLmers[i];
    }
    if(runningSum == 0)
    {
        double equalLikely = 1/(double)lastIndex;
        for(int i = 0; i < lastIndex; i++)
            probLmers[i] = equalLikely;
    }
    else
    {
        for(int i = 0; i < lastIndex; i++)
        {
            probLmers[i] /= runningSum;
        }
    }
    return probLmers;
}

double GibbsSearchEngine::probOfLmerAt(int sequenceID, int sequenceIndex, vector<vector<double> >& profileMatrix) {
    double prob = 1;
    Nucleotide_t lMerIterator;
    double columnCount;
    double columnSize = repo.Size() - 1;
    for(int i = 0; i < motifLength; i++)
    {
        lMerIterator = repo.Get(sequenceID, i + sequenceIndex);
        columnCount = profileMatrix[lMerIterator][i];
        prob *= columnCount / columnSize;
    }
    return prob;
}



double GibbsSearchEngine::GetBestScore() {
    return score;
}

int GibbsSearchEngine::GetDontCareCount() {
    return dontCares;
}

std::vector<Nucleotide_t> GibbsSearchEngine::GetMotif() {
    return motif;
}

int GibbsSearchEngine::GetMotifLength() {
    return motifLength;
}

std::vector<int> GibbsSearchEngine::GetStartingLoci() {
    return loci;
}

void GibbsSearchEngine::SetDontCares(int number) {
    if(number <= motifLength - 2  &&  number >= 0)
        dontCares = number;
    else
        cout << "Dont cares cannot be greater than " << motifLength - 2 << ".\n"; 
}

void GibbsSearchEngine::SetMotifLength(int length) {
    if(length > 0)
        motifLength = length;
    else
        cout << "Motif length must be positive.\n";       
}

void GibbsSearchEngine::SetRepo(IDnaRepository& input) {
    repo = input;
}




