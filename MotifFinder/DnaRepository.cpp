/* 
 * File:   DnaRepository.cpp
 * Author: d
 * 
 * Created on November 16, 2014, 6:44 PM
 */
#include <cstdlib>
#include <vector>
#include "DnaRepository.h"
#include "DnaSequenceVector.h"
using namespace std;

DnaRepository::DnaRepository(int max_size) {
            sequences = std::vector<DnaSequenceVector> ();
}

void DnaRepository::Add(DnaSequenceVector& input) {
    sequences.push_back(input);
}

int DnaRepository::Count(Nucleotide_t input) {
    int count = 0;
    for(int i = 0; i < Size(); i++)
    {
        count += Count(i, input);
    }
    return count;
}

int DnaRepository::Count(int sequence, Nucleotide_t input) {
    int count = 0;
    for(int i = 0; i < Size(sequence); i++)
    {
        if(Get(sequence, i) == input)
            count++;
    }
    return count;
}

IDnaSequence& DnaRepository::Get(int i) {
    return sequences[i];
}

Nucleotide_t DnaRepository::Get(int sequence, int position) {
    return Get(sequence).Get(position);
}

void DnaRepository::Set(int sequence, int position, Nucleotide_t input) {
    Get(sequence).Set(position, input);
}

int DnaRepository::Size() {
    return sequences.size();
}

int DnaRepository::Size(int i) {
    return Get(i).Size();
}   


