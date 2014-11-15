/* 
 * File:   RandomizedSearchEngine.h
 * Author: Eric
 *
 * Created on November 15, 2014, 1:43 PM
 */

#ifndef RANDOMIZEDSEARCHENGINE_H
#define	RANDOMIZEDSEARCHENGINE_H

#include "ISearchEngine.h"
#include "IScoreEngine.h"

using namespace std;

class RandomizedSearchEngine : ISearchEngine{
public:
    RandomizedSearchEngine();
    RandomizedSearchEngine(IDnaRepository& input);
    virtual void SetRepo(IDnaRepository* input);
    virtual std::vector<Nucleotide_t> GetMotif();
    virtual std::vector<int> GetStartingLoci();
    virtual void Search(int time, int motifLength, int dontCares);



    virtual ~RandomizedSearchEngine();
private:
    int randomPosition(int sequence);
    IScoreEngine& scoreEngine;
    IDnaRepository* dna;
    vector<Nucleotide_t> bestMotif;
    vector<int> startingLoci;
    vector<vector<Nucleotide_t> > profileMatrix;

};

#endif	/* RANDOMIZEDSEARCHENGINE_H */

