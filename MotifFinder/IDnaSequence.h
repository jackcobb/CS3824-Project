/* 
 * File:   IDnaSequence.h
 * Author: Eric
 *
 * Created on November 3, 2014, 1:47 PM
 */

#ifndef IDNASEQUENCE_H
#define	IDNASEQUENCE_H
#include "Nucleotide.h"

class IDnaSequence {
public:
    virtual int Size() = 0;
    virtual Nucleotide_t Get(int i) = 0;
    virtual void Set(int i, Nucleotide_t input) = 0;
private:

};
#endif