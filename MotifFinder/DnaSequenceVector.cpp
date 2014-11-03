/* 
 * File:   DnaSequenceVector.cpp
 * Author: Eric
 * 
 * Created on November 3, 2014, 2:18 PM
 */
#include<iostream>
#include<vector>
using namespace std;

class DnaSequenceVector : IDnaSequence
{
    public:
        
        DnaSequenceVector(int max_size) {
            members = std::vector<Nucleotide_t>(max_size);

        }
        int Size()
        {
            return members.size();
        }
        
        Nucleotide_t Get(int i)
        {
            try{
                if(i >= Count() || i < 0)
                    throw i;
                return members[i];
            }
            catch (int badIterator)
            {
                cout << i << " is a bad iterator. Must be in range [0," << members.size() - 1 << "]\n";
            }
        }
    
        void Set(int i, Nucleotide_t input)
        {
            try{
                if(i >= DnaSequenceVector::Count() || i < 0)
                    throw i;
                members[i] = input;
            }
            catch (int badIterator)
            {
                cout << i << " is a bad iterator. Must be in range [0," << DnaSequenceVector::members.size() - 1 << "]\n";
            }

        }
        
        ~DnaSequenceVector()
        {
            delete &members;
        }
    private:
        std::vector<Nucleotide_t> members;
};

