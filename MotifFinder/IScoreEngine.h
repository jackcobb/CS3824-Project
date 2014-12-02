/* 
 * File:   IScoreEngine.h
 * Author: Eric
 *
 * Created on November 3, 2014, 6:59 PM
 */

#ifndef ISCOREENGINE_H
#define	ISCOREENGINE_H

#include "Nucleotide.h"
#include "IDnaRepository.h"

using namespace std;

class IScoreEngine {
public:
    virtual void SetRepo(IDnaRepository& input) =0;
    virtual double Score(vector<Nucleotide_t>& Motif, vector<int>& StartingLoci) = 0;
    virtual vector<Nucleotide_t> optimizeDontCaresInMotif(vector<Nucleotide_t>& motif, int dontCares, vector<int>& loci) =0;
    
private:

};

#endif	/* ISCOREENGINE_H */

