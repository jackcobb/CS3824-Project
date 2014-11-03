/* 
 * File:   IDnaRepository.cpp
 * Author: Eric
 * 
 * Created on November 3, 2014, 2:07 PM
 */

#include "IDnaRepository.h"

Neocleotide_t IDnaRepository::Get(int sequence, int position)
{
    IDnaSequence& parentSequence = Get(sequence);
    return parentSequence.Get(position);
}

void IDnaRepository::Set(int sequence, int position, Neocleotide_t input)
{
    IDnaSequence& parentSequence = Get(sequence);
    parentSequence.Set(position, input);
}
