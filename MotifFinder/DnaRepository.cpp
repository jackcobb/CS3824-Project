/* 
 * File:   DnaSequenceVector.cpp
 * Author: Eric
 * 
 * Created on November 15, 2014, 3:03 PM
 */
#include <cstdlib>
#include <vector>
#include "DnaRepository.h"
#include "DnaSequenceVector.h"
using namespace std;

DnaRepository::DnaRepository(int max_size) {
            sequences = std::vector<DnaSequenceVector> (max_size, DnaSequenceVector(512));
}

void DnaRepository::Add(IDnaSequence& input) {

}

int DnaRepository::Count(Nucleotide_t input) {
    return 0;
}

int DnaRepository::Count(int sequence, Nucleotide_t input) {
    return 0;
}

IDnaSequence DnaRepository::Get(int i) {
    DnaSequenceVector a(123);
    return a;
}

Nucleotide_t DnaRepository::Get(int sequence, int position) {
    return DC;
}

void DnaRepository::Set(int sequence, int position, Nucleotide_t input) {

}

int DnaRepository::Size() {
    return 0;
}

int DnaRepository::Size(int i) {
    return 0;
}   










