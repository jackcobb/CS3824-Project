/* 
 * File:   RandomizedSearchEngine.h
 * Author: Eric
 *
 * Created on November 15, 2014, 1:43 PM
 */

#ifndef RANDOMIZEDSEARCHENGINE_H
#define	RANDOMIZEDSEARCHENGINE_H

#include "ISearchEngine.h"
#include "ScoreEngine.h"

using namespace std;

class RandomizedSearchEngine : ISearchEngine{
public:
    RandomizedSearchEngine(IDnaRepository* input, int motiflength, int dontcares);
    virtual void SetRepo(IDnaRepository& input);
    virtual std::vector<Nucleotide_t> GetMotif();
    virtual std::vector<int> GetStartingLoci();
    virtual void Search(double time);



    virtual ~RandomizedSearchEngine();
private:
    int randomPosition(int sequence);
    vector<int> randomLoci();
    vector<Nucleotide_t> lociToMotif(vector<int> loci);
    ScoreEngine& scoreEngine;
    IDnaRepository& dna;
    vector<Nucleotide_t> motif;
    int motifLength;
    int dontCares;
    vector<int> startingLoci;
    vector<vector<int> > profileMatrix;
    bool canIncrement(vector<Nucleotide_t> motif);
    vector<Nucleotide_t> increment(vector<Nucleotide_t> motif, vector<Nucleotide_t> startMotif);
    int getMax(int a, int t, int g, int c);
    vector<Nucleotide_t> getStartingMotif(vector<vector<int> > profileMatrix, vector<Nucleotide_t> motif);
    vector<vector<int> > createProfileMatrix(vector<int> loci);
};

#endif	/* RANDOMIZEDSEARCHENGINE_H */

