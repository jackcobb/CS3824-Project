/* 
 * File:   Runner.cpp
 * Author: Eric
 * 
 * Created on November 16, 2014, 8:10 PM
 */

#include<iostream>
#include "Runner.h"

Runner::Runner(IDnaRepository& dna, ISearchEngine& searchEngine) : repo(dna), searchEngine(searchEngine), scoreEngine(ScoreEngine(dna)) {

}

void Runner::Run(double runTime) {
    searchEngine.Search(runTime);
    int motifLength = searchEngine.GetMotifLength();
    int dontCares = searchEngine.GetDontCareCount();
    vector<Nucleotide_t> motif = searchEngine.GetMotif();
    string motifString = MotifToString(motif);
    double score = searchEngine.GetBestScore();
    cout << "Best motif of length " << motifLength << " and " << dontCares << " don't cares is " <<
            motifString << "\n";
    cout << "Log likelihood is " << searchEngine.GetBestScore() << "\n";
    cout << "Loci of the best motif are here:\n";
    vector<int> loci = searchEngine.GetStartingLoci();
    for(int i = 0; i < loci.size(); i++) {
       cout << loci[i] << "\n";
    }
    
           
}

string Runner::MotifToString(vector<Nucleotide_t>& input) {
    string output;
    for(int i = 0; i < input.size(); i++)
    {
        if(input[i] == A){
            output.append("A");
        }
        else if(input[i] == T){
            output.append("T");
        }
        else if(input[i] == G){
            output.append("G");
        }
        else if(input[i] == C){
            output.append("C");
        }
        else if(input[i] == DC){
            output.append("*");
        }
    }
    return output;
}



