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
    virtual void SetRepo(IDnaRepository& input);
    virtual std::vector<Nucleotide_t> GetMotif();
    virtual std::vector<int> GetStartingLoci();
    virtual void Search(double time, int motifLength, int dontCares);



    virtual ~RandomizedSearchEngine();
private:
    int randomPosition(int sequence, int motifLength);
    vector<int> randomLoci(int motifLength);
    vector<Nucleotide_t> lociToMotif(vector<int> loci, int motifLength, int dontCares);
    IScoreEngine& scoreEngine;
    IDnaRepository& dna;
    vector<Nucleotide_t> motif;
    vector<int> startingLoci;
    vector<vector<Nucleotide_t> > profileMatrix;

};

#endif	/* RANDOMIZEDSEARCHENGINE_H */

