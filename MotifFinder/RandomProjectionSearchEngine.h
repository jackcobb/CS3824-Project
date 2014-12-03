/* 
 * File:   RandomProjectionSearchEngine.h
 * Author: Eric
 *
 * Created on December 2, 2014, 3:26 PM
 */

#ifndef RANDOMPROJECTIONSEARCHENGINE_H
#define	RANDOMPROJECTIONSEARCHENGINE_H

#include "ISearchEngine.h"
#include "SearchResult.h"


class RandomProjectionSearchEngine : public ISearchEngine{
public:
    RandomProjectionSearchEngine();

    virtual double GetBestScore();
    
    virtual int GetDontCareCount();

    virtual std::vector<Nucleotide_t> GetMotif();

    virtual int GetMotifLength();

    virtual std::vector<int> GetStartingLoci();

    virtual void Search(double time);

    virtual SearchResult Search(double runTime, int threshold, int numberIterations);
    
    virtual void SetDontCares(int number);
    
    virtual void SetMotifLength(int length);

    virtual void SetRepo(IDnaRepository& input);
    
private:
    int dontCares;
    int motifLength;
    vec

};

#endif	/* RANDOMPROJECTIONSEARCHENGINE_H */

