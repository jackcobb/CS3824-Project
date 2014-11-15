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
    virtual ~RandomizedSearchEngine();
private:
    int randomPosition(int sequence);
    IScoreEngine& scoreEngine;
    vector<Nucleotide_t> bestMotif;
    vector<int> startingLoci;
    vector<vector<Nucleotide_t> > ProfileMatrix;

};

#endif	/* RANDOMIZEDSEARCHENGINE_H */

