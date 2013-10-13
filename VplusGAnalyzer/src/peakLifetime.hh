#ifndef PEAKLIFETIME_HH_
#define PEAKLIFETIME_HH_

#include "fftjet/Peak.hh"
#include "fftjet/ProximityClusteringTree.hh"
#include "fftjet/SparseClusteringTree.hh"

typedef fftjet::SparseClusteringTree<fftjet::Peak,long> SparseTree;

double peakSplitTime(const SparseTree& tree, SparseTree::NodeId id,
                     double minScale);

double peakMergeTime(const SparseTree& tree, SparseTree::NodeId id,
                     double maxScale);

#endif // PEAKLIFETIME_HH_
