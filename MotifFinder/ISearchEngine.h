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
        virtual void SetRepo(IDnaRepository* input) = 0;
        virtual void Search(int time) = 0;
        virtual std::vector<Nucleotide_t> GetMotif();
        virtual std::vector<int> GetStartingLoci();
};

#endif	/* ISEARCHENGINE_H */

