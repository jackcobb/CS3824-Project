/* 
 * File:   ISearchEngine.h
 * Author: Eric
 *
 * Created on November 15, 2014, 1:31 PM
 */

#ifndef ISEARCHENGINE_H
#define	ISEARCHENGINE_H

#include <vector>
#include "Nucleotide.h"
#include "IDnaRepository.h"

class ISearchEngine
{
    public:
        virtual void SetRepo(IDnaRepository& input) = 0;
        virtual void Search(double time) = 0;
        virtual int GetMotifLength() = 0;
        virtual int GetDontCareCount() = 0;
        virtual void SetMotifLength(int length) = 0;
        virtual void SetDontCares(int number) = 0;
        virtual double GetBestScore() = 0;
        virtual std::vector<Nucleotide_t> GetMotif() = 0;
        virtual std::vector<int> GetStartingLoci() = 0;
};

#endif	/* ISEARCHENGINE_H */

