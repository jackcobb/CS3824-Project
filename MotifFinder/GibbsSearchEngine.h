/* 
 * File:   GibbsSearchEngine.h
 * Author: Eric
 *
 * Created on November 26, 2014, 2:06 PM
 */

#ifndef GIBBSSEARCHENGINE_H
#define	GIBBSSEARCHENGINE_H

#include <vector>
#include "ISearchEngine.h"
#include "ScoreEngine.h"
#include "SearchResult.h"

using namespace std;

class GibbsSearchEngine : public ISearchEngine {
public:
    GibbsSearchEngine(IDnaRepository& repo, int motifLength, int dontCares);
    virtual double GetBestScore();

    virtual std::vector<Nucleotide_t> GetMotif();

    virtual int GetDontCareCount();

    virtual int GetMotifLength(); 

    virtual std::vector<int> GetStartingLoci();

    virtual void Search(double time);
    
    virtual SearchResult Search(double maxRunTime, int maxIterations);
    
    virtual SearchResult Search(double maxRunTime, int maxIterations, vector<int>& startingLoci);

    virtual void SetDontCares(int number);

    virtual void SetMotifLength(int length);

    virtual void SetRepo(IDnaRepository& input);
private:
    int dontCares;
    int motifLength;
    vector<Nucleotide_t> motif;
    vector<int> loci;
    double score;
    IDnaRepository& repo;
    ScoreEngine scoreEngine;
    vector<vector<double> > profileMatrixGibbs;
    vector<vector<double> > profileMatrixFull;
    vector<double> probDistri;
    
    vector<int> getRandomLoci();
    int randomPosition(int i);
    int getRandomSequence();
    void updateProfileMatrix(vector<int> loci);
    void updateProfileMatrixExcluding(int excludedIndex, vector<int> loci);
    vector<double> probabilityLMer(int sequence);
    double probOfLmerAt(int sequenceID, int sequenceIndex);
    vector<Nucleotide_t> lociToMotif(vector<int> loci);
    Nucleotide_t getMax(int a, int t, int g, int c);
    vector<Nucleotide_t> getStartingMotif();
    
    

};

#endif	/* GIBBSSEARCHENGINE_H */

