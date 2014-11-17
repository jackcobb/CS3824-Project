/* 
 * File:   DnaSequenceVector.cpp
 * Author: Eric
 * 
 * Created on November 15, 2014, 3:03 PM
 */
#include<iostream>
#include<vector>
#include<stdexcept>
#include "IDnaSequence.h"
#include "DnaSequenceVector.h"

using namespace std;

        DnaSequenceVector::DnaSequenceVector() {
            members = std::vector<Nucleotide_t>();
        }
        int DnaSequenceVector::Size()
        {
            return members.size();
        }
        
        Nucleotide_t DnaSequenceVector::Get(int i)
        {
            Nucleotide_t nucleotide;
            try{
                nucleotide = members[i];
            }
            catch (const std::out_of_range& oor)
            {
                cout << i << " is a bad iterator. Must be in range [0," << members.size() - 1 << "]\n";
                std::cerr << "Out of Range error: " << oor.what() << '\n';
            }

            return nucleotide;
        }
    
        void DnaSequenceVector::Set(int i, Nucleotide_t input)
        {
            try{
                members[i] = input;
            }
            catch (const std::out_of_range& oor)
            {
                cout << i << " is a bad iterator. Must be in range [0," << DnaSequenceVector::members.size() - 1 << "]\n";
                std::cerr << "Out of Range error: " << oor.what() << '\n';
            }

        }
        
        void DnaSequenceVector::Push(Nucleotide_t input)
        {
            members.push_back(input);
        }
        
        DnaSequenceVector::~DnaSequenceVector()
        {
        }
        
