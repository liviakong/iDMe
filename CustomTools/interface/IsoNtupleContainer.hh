#ifndef ISONTUPLECONTAINER_HH
#define ISONTUPLECONTAINER_HH

#include <vector>
using std::vector;
#include <iostream>

#include <TTree.h>
#include "NtupleContainer.hh"

class IsoNtupleContainer : public NtupleContainer {

public:
    IsoNtupleContainer();
    virtual ~IsoNtupleContainer();
    void CreateTreeBranches();
    void ClearTreeBranches();

    // PF Candidates
    int nPFCand_;
    vector<int> pfID_;
    vector<int> pfCharge_;
    vector<float> pfPt_;
    vector<float> pfEta_;
    vector<float> pfPhi_;
    vector<float> pfE_;
};


#endif