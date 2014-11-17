/* 
 * File:   DnaRepository.h
 * Author: d
 *
 * Created on November 16, 2014, 6:44 PM
 */

#ifndef DNAREPOSITORY_H
#define	DNAREPOSITORY_H
#include <vector>
#include "DnaSequenceVector.h"
#include "IDnaRepository.h"

class DnaRepository : public IDnaRepository{
public:
    DnaRepository(int max_size);
    int Size();
    int Size(int i);
    int Count(Nucleotide_t input);
    int Count(int sequence, Nucleotide_t input);
    IDnaSequence& Get(int i);
    Nucleotide_t Get(int sequence, int position);
    void Set(int sequence, int position, Nucleotide_t input);
    void Add(DnaSequenceVector& input);

private:
    std::vector<DnaSequenceVector> sequences;
};

#endif	/* DNAREPOSITORY_H */

