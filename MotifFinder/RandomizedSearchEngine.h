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

class RandomizedSearchEngine : public ISearchEngine{
public:
    RandomizedSearchEngine(IDnaRepository& input, int motiflength, int dontcares);
    virtual void SetRepo(IDnaRepository& input);
    virtual std::vector<Nucleotide_t> GetMotif();
    virtual std::vector<int> GetStartingLoci();
    virtual void Search(double time);
    virtual int GetDontCareCount();
    virtual int GetMotifLength();
    virtual void SetDontCares(int number);
    virtual void SetMotifLength(int length);
    virtual double GetBestScore();
    virtual ~RandomizedSearchEngine();
private:
    int randomPosition(int sequence);
    vector<int> randomLoci();
    vector<Nucleotide_t> lociToMotif(vector<int> loci);
    ScoreEngine scoreEngine;
    IDnaRepository& dna;
    vector<Nucleotide_t> motif;
    int motifLength;
    int dontCares;
    double score;
    vector<int> startingLoci;
    vector<vector<int> > profileMatrix;
    bool canIncrement(vector<Nucleotide_t> motif);
    vector<Nucleotide_t> increment(vector<Nucleotide_t> motif, vector<Nucleotide_t> startMotif);
    Nucleotide_t getMax(int a, int t, int g, int c);
    vector<Nucleotide_t> getStartingMotif(vector<vector<int> > profileMatrix);
    vector<vector<int> > createProfileMatrix(vector<int> loci);
};

#endif	/* RANDOMIZEDSEARCHENGINE_H */

