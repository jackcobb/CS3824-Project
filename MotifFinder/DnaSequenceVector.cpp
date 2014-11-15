/* 
 * File:   DnaSequenceVector.cpp
 * Author: Eric
 * 
 * Created on November 15, 2014, 3:03 PM
 */
#include<iostream>
#include<vector>
#include "DnaSequenceVector.h"
using namespace std;

        DnaSequenceVector::DnaSequenceVector(int max_size) {
            members = std::vector<Nucleotide_t>(max_size);
        }
        int DnaSequenceVector::Size()
        {
            return members.size();
        }
        
        Nucleotide_t DnaSequenceVector::Get(int i)
        {
            try{
                if(i >= Size() || i < 0)
                    throw i;
                return members[i];
            }
            catch (int badIterator)
            {
                cout << i << " is a bad iterator. Must be in range [0," << members.size() - 1 << "]\n";
            }
        }
    
        void DnaSequenceVector::Set(int i, Nucleotide_t input)
        {
            try{
                if(i >= DnaSequenceVector::Size() || i < 0)
                    throw i;
                members[i] = input;
            }
            catch (int badIterator)
            {
                cout << i << " is a bad iterator. Must be in range [0," << DnaSequenceVector::members.size() - 1 << "]\n";
            }

        }

