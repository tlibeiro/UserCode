#include <cmath>
#include <cassert>

#include "peakLifetime.hh"

double peakSplitTime(const SparseTree& tree, const SparseTree::NodeId id,
                     const double minScale)
{
    assert(id);
    assert(minScale > 0.0);
    const SparseTree::Node& myNode = tree.getNode(id);
    const fftjet::Peak& myPeak = myNode.getCluster();
    const double myScale = myPeak.scale();
    assert(myScale > 0.0);
    const unsigned nDaus = myNode.nDaus();
    const double normalization(log(myScale)-log(minScale));
    
    if (nDaus && myScale > minScale)
    {
        double firstDaugherFraction = 1.0;
        if (nDaus > 1)
        {
            double et0 = 0.0, etSum = 0.0;
            for (unsigned idau=0; idau<nDaus; ++idau)
            {
                const SparseTree::NodeId dauId = myNode.daus()[idau];
                const SparseTree::Node& dauNode = tree.uncheckedNode(dauId);
                const fftjet::Peak& dauPeak = dauNode.getCluster();
                const double s = dauPeak.scale();
                const double et = s*s*dauPeak.magnitude();
                etSum += et;
                if (!idau)
                    et0 = et;
            }
            assert(et0 > 0.0);
            assert(etSum > 0.0);
            firstDaugherFraction = et0/etSum;
        }
        const SparseTree::NodeId dauId = myNode.daus()[0];
        const SparseTree::Node& dauNode = tree.uncheckedNode(dauId);
        const fftjet::Peak& dauPeak = dauNode.getCluster();
        const double dauScale = dauPeak.scale();
        assert(myScale > dauScale);
        double t;
        if (dauScale > minScale)
            t = log(myScale/dauScale) + peakSplitTime(tree, dauId, minScale);
        else
            t = log(myScale/minScale);
        return firstDaugherFraction*t;
    }
    else
        return log(myScale/minScale)/normalization;
}

double peakMergeTime(const SparseTree& tree, const SparseTree::NodeId id,
                     const double maxScale)
{
    assert(id);
    assert(maxScale > 0.0);
    const SparseTree::Node& myNode = tree.getNode(id);
    const fftjet::Peak& myPeak = myNode.getCluster();
    const double myScale = myPeak.scale();
    assert(myScale > 0.0);
    const SparseTree::NodeId parentId = myNode.parent();
    const double normalization(log(maxScale)-log(myScale));

    if (parentId && myScale < maxScale)
    {
        const SparseTree::Node& parentNode = tree.getNode(parentId);
        const fftjet::Peak& parentPeak = parentNode.getCluster();
        const double parentScale = parentPeak.scale();
        assert(parentScale > myScale);
        double t;
        if (parentScale >= maxScale)
            t = log(maxScale/myScale);
        else
        {
            t = log(parentScale/myScale);
            double myFraction = 1.0;
            const unsigned nDaus = parentNode.nDaus();
            if (nDaus > 1)
            {
                double et0 = 0.0, etSum = 0.0;
                for (unsigned idau=0; idau<nDaus; ++idau)
                {
                    const SparseTree::NodeId dauId = parentNode.daus()[idau];
                    const SparseTree::Node& dauNode = tree.uncheckedNode(dauId);
                    const fftjet::Peak& dauPeak = dauNode.getCluster();
                    const double s = dauPeak.scale();
                    const double et = s*s*dauPeak.magnitude();
                    etSum += et;
                    if (dauId == id)
                        et0 = et;
                }
                assert(et0 > 0.0);
                assert(etSum > 0.0);
                myFraction = et0/etSum;
            }
            t += myFraction*peakMergeTime(tree, parentId, maxScale);
        }
        return t;
    }
    else
        return log(maxScale/myScale)/normalization;
}
