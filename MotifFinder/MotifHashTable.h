/* 
 * File:   MotifHashTable.h
 * Author: Eric
 *
 * Created on December 2, 2014, 4:21 PM
 */

#ifndef MOTIFHASHTABLE_H
#define	MOTIFHASHTABLE_H
#include <vector>

#include "Nucleotide.h"
using namespace std;
class MotifHashTable {
public:
    MotifHashTable();
    MotifHashTable(const MotifHashTable& orig);
    int Count(vector<Nucleotide_t> input);
    void Increment(vector<Nucleotide_t> input);
    virtual ~MotifHashTable();
private:
    vector<vector<int> > hashTable;

};

#endif	/* MOTIFHASHTABLE_H */

