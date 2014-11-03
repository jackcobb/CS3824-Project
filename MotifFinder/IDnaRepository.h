/* 
 * File:   IDnaRepository.h
 * Author: Eric
 *
 * Created on November 3, 2014, 2:07 PM
 */

#ifndef IDNAREPOSITORY_H
#define	IDNAREPOSITORY_H

#include "IDnaSequence.h"


class IDnaRepository {
public:
    virtual int Count();
    virtual IDnaSequence& Get(int i) =0;
    Neocleotide_t Get(int sequence, int position);
    void Set(int sequence, int position, Neocleotide_t input);
    virtual void Add(IDnaSequence& input) = 0;
    virtual ~IDnaRepository() = 0;
private:

};

#endif	/* IDNAREPOSITORY_H */

