/* 
 * File:   IDnaSequence.h
 * Author: Eric
 *
 * Created on November 3, 2014, 1:47 PM
 */

#ifndef IDNASEQUENCE_H
#define	IDNASEQUENCE_H
#include "Neocleotide.h"9

class IDnaSequence {
public:
    virtual int Count() = 0;
    virtual Neocleotide_t Get(int i) = 0;
    virtual void Set(int i, Neocleotide_t input) = 0;
    virtual ~IDnaSequence();
private:

};

#endif	/* IDNASEQUENCE_H */
