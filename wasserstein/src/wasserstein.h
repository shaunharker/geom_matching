/*
 
Copyright (c) 2015, M. Kerber, D. Morozov, A. Nigmetov
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 

You are under no obligation whatsoever to provide any bug fixes, patches, or
upgrades to the features, functionality or performance of the source code
(Enhancements) to anyone; however, if you choose to make your Enhancements
available either publicly, or directly to copyright holder,
without imposing a separate written license agreement for such Enhancements,
then you hereby grant the following license: a  non-exclusive, royalty-free
perpetual license to install, use, modify, prepare derivative works, incorporate
into other computer software, distribute, and sublicense such enhancements or
derivative works thereof, in binary and source code form.

  */

#ifndef WASSERSTEIN_H
#define WASSERSTEIN_H

#include <vector>
#include <queue>
#include <memory>
#include <boost/heap/d_ary_heap.hpp>

#include "basic_defs.h"

#include <dnn/geometry/euclidean-fixed.h>
#include <dnn/local/kd-tree.h>

using IdxType = int;
using IdxValPair = std::pair<IdxType, double>;

struct CompPairsBySecondStruct {
    bool operator()(const IdxValPair& a, const IdxValPair& b) const
    {
        return a.second < b.second;
    }
};

// 
struct CompPairsBySecondGreaterStruct {
    bool operator()(const IdxValPair& a, const IdxValPair& b) const
    {
        return a.second > b.second;
    }
};


using ItemsTimePair = std::pair<IdxType, int>;

using UpdateList = std::list<ItemsTimePair>;
using UpdateListIter = UpdateList::iterator;


using BiddersPriorityQueue = std::priority_queue< IdxValPair,
                                                  std::vector<IdxValPair>,
                                                  CompPairsBySecondStruct >;

struct DebugOptimalBid {
    DebugOptimalBid() : bestItemIdx(-1), bestItemValue(-666.666), secondBestItemIdx(-1), secondBestItemValue(-666.666) {};
    IdxType bestItemIdx;
    double bestItemValue;
    IdxType secondBestItemIdx;
    double secondBestItemValue;
};

struct AuctionOracleAbstract {
    AuctionOracleAbstract(const std::vector<DiagramPoint>& _bidders, const std::vector<DiagramPoint>& _items, const double _wassersteinPower, const double _internal_p = std::numeric_limits<double>::infinity());
    ~AuctionOracleAbstract() {}
    virtual IdxValPair getOptimalBid(const IdxType bidderIdx) = 0;
    virtual void setPrice(const IdxType itemsIdx, const double newPrice) = 0;
    virtual void adjustPrices(void) = 0;
    double getEpsilon() { return epsilon; };
    virtual void setEpsilon(double newEpsilon) { assert(newEpsilon >= 0.0); epsilon = newEpsilon; };
protected:
    const std::vector<DiagramPoint>& bidders;
    const std::vector<DiagramPoint>& items;
    std::vector<double> prices;
    double wassersteinPower;
    double epsilon;
    double internal_p;
    double getValueForBidder(size_t bidderIdx, size_t itemsIdx);
};

struct AuctionOracleLazyHeap final : AuctionOracleAbstract {
    using AuctionHeap = boost::heap::d_ary_heap<IdxValPair, boost::heap::arity<2>, boost::heap::mutable_<true>, boost::heap::compare<CompPairsBySecondStruct>>;
    AuctionOracleLazyHeap(const std::vector<DiagramPoint>& bidders, const std::vector<DiagramPoint>& items, const double wassersteinPower, const double _internal_p = std::numeric_limits<double>::infinity());
    ~AuctionOracleLazyHeap();
    // data members
    // temporarily make everything public
    std::vector<std::vector<double>> weightMatrix;
    double weightAdjConst;
    double maxVal;
    // vector of heaps to find the best items
    std::vector<AuctionHeap*> profitHeap;
    std::vector<std::vector<AuctionHeap::handle_type>> profitHeapHandles;
    // methods
    void fillInProfitHeap(void);
    void setPrice(const IdxType itemsIdx, const double newPrice) override final;
    IdxValPair getOptimalBid(const IdxType bidderIdx) override final;
    double getMatchingWeight(const std::vector<IdxType>& biddersToItems) const;
    void adjustPrices(void) override final;
    // to update the queue in lazy fashion
    std::vector<UpdateListIter> itemsIterators;
    UpdateList updateList;
    std::vector<int> biddersUpdateMoments;
    int updateCounter;
    void updateQueueForBidder(const IdxType bidderIdx);
    // debug
    DebugOptimalBid getOptimalBidDebug(const IdxType bidderIdx);
};

struct AuctionOracleLazyHeapRestricted final : AuctionOracleAbstract {
    using AuctionHeap = boost::heap::d_ary_heap<IdxValPair, boost::heap::arity<2>, boost::heap::mutable_<true>, boost::heap::compare<CompPairsBySecondGreaterStruct>>;
    using DiagPricesHeap = boost::heap::d_ary_heap<IdxValPair,
                                                   boost::heap::arity<2>,
                                                   boost::heap::mutable_<true>,
                                                   boost::heap::compare<CompPairsBySecondGreaterStruct>>;
     AuctionOracleLazyHeapRestricted(const std::vector<DiagramPoint>& bidders, const std::vector<DiagramPoint>& items, const double wassersteinPower, const double _internal_p = std::numeric_limits<double>::infinity());
    ~AuctionOracleLazyHeapRestricted();
    // data members
    // temporarily make everything public
    std::vector<std::vector<double>> weightMatrix;
    double weightAdjConst;
    double maxVal;
    // vector of heaps to find the best items
    std::vector<AuctionHeap*> profitHeap;
    std::vector<std::vector<size_t>> itemsIndicesForHeapHandles;
    std::vector<std::vector<AuctionHeap::handle_type>> profitHeapHandles;
    // methods
    void fillInProfitHeap(void);
    void setPrice(const IdxType itemsIdx, const double newPrice) override final;
    IdxValPair getOptimalBid(const IdxType bidderIdx) override final;
    double getMatchingWeight(const std::vector<IdxType>& biddersToItems) const;
    void adjustPrices(void) override final;
    // to update the queue in lazy fashion
    std::vector<UpdateListIter> itemsIterators;
    UpdateList updateList;
    std::vector<int> biddersUpdateMoments;
    int updateCounter;
    std::vector<size_t> biddersToProjItems;
    void updateQueueForBidder(const IdxType bidderIdx);
    DiagPricesHeap diagItemsHeap;
    std::vector<DiagPricesHeap::handle_type> diagHeapHandles;
    std::vector<size_t> heapHandlesIndices;
    // debug
   
    DebugOptimalBid getOptimalBidDebug(const IdxType bidderIdx);
    
    // for diagonal points
    bool bestDiagonalItemsComputed;
    size_t bestDiagonalItemIdx;
    double bestDiagonalItemValue;
    size_t secondBestDiagonalItemIdx;
    double secondBestDiagonalItemValue;
};

struct AuctionOracleKDTree final : AuctionOracleAbstract {
    typedef dnn::Point<2, double> DnnPoint;
    typedef dnn::PointTraits<DnnPoint> DnnTraits;

    using DiagPricesHeap = boost::heap::d_ary_heap<IdxValPair,
                                                   boost::heap::arity<2>,
                                                   boost::heap::mutable_<true>,
                                                   boost::heap::compare<CompPairsBySecondGreaterStruct>>;
    
    AuctionOracleKDTree(const std::vector<DiagramPoint>& bidders, const std::vector<DiagramPoint>& items, const double wassersteinPower, const double _internal_p = std::numeric_limits<double>::infinity());
    ~AuctionOracleKDTree();
    // data members
    // temporarily make everything public
    double maxVal;
    double weightAdjConst;
    dnn::KDTree<DnnTraits>* kdtree;
    std::vector<DnnPoint> dnnPoints;
    std::vector<DnnPoint*> dnnPointHandles;
    dnn::KDTree<DnnTraits>* kdtreeAll;
    std::vector<DnnPoint> dnnPointsAll;
    std::vector<DnnPoint*> dnnPointHandlesAll;
    DiagPricesHeap diagItemsHeap;
    std::vector<DiagPricesHeap::handle_type> diagHeapHandles;
    std::vector<size_t> heapHandlesIndices;
    std::vector<size_t> kdtreeItems;
    // vector of heaps to find the best items
    void setPrice(const IdxType itemsIdx, const double newPrice) override final;
    IdxValPair getOptimalBid(const IdxType bidderIdx) override final;
    void adjustPrices(void) override final;
    // debug routines
    DebugOptimalBid getOptimalBidDebug(IdxType bidderIdx);
    void setEpsilon(double newVal) override final;
};

struct AuctionOracleKDTreeRestricted final : AuctionOracleAbstract {
    typedef dnn::Point<2, double> DnnPoint;
    typedef dnn::PointTraits<DnnPoint> DnnTraits;

    using DiagPricesHeap = boost::heap::d_ary_heap<IdxValPair,
                                                   boost::heap::arity<2>,
                                                   boost::heap::mutable_<true>,
                                                   boost::heap::compare<CompPairsBySecondGreaterStruct>>;
    
    AuctionOracleKDTreeRestricted(const std::vector<DiagramPoint>& bidders, const std::vector<DiagramPoint>& items, const double wassersteinPower, const double _internal_p = std::numeric_limits<double>::infinity());
    ~AuctionOracleKDTreeRestricted();
    // data members
    // temporarily make everything public
    double maxVal;
    double weightAdjConst;
    dnn::KDTree<DnnTraits>* kdtree;
    std::vector<DnnPoint> dnnPoints;
    std::vector<DnnPoint*> dnnPointHandles;
    std::vector<DnnPoint> dnnPointsAll;
    std::vector<DnnPoint*> dnnPointHandlesAll;
    DiagPricesHeap diagItemsHeap;
    std::vector<DiagPricesHeap::handle_type> diagHeapHandles;
    std::vector<size_t> heapHandlesIndices;
    std::vector<size_t> kdtreeItems;
    std::vector<size_t> biddersToProjItems;
    // vector of heaps to find the best items
    void setPrice(const IdxType itemsIdx, const double newPrice) override final;
    IdxValPair getOptimalBid(const IdxType bidderIdx) override final;
    void adjustPrices(void) override final;
    // debug routines
    DebugOptimalBid getOptimalBidDebug(IdxType bidderIdx);
    void setEpsilon(double newVal) override final;


    bool bestDiagonalItemsComputed;
    size_t bestDiagonalItemIdx;
    double bestDiagonalItemValue;
    size_t secondBestDiagonalItemIdx;
    double secondBestDiagonalItemValue;
};

struct AuctionOracleRestricted final : AuctionOracleAbstract {
    AuctionOracleRestricted(const std::vector<DiagramPoint>& bidders, const std::vector<DiagramPoint>& items, const double wassersteinPower, const double _internal_p = std::numeric_limits<double>::infinity());
    IdxValPair getOptimalBid(const IdxType bidderIdx) override;
    void setPrice(const IdxType itemsIdx, const double newPrice) override;
    void adjustPrices(void) override {};
    void setEpsilon(double newEpsilon) override { assert(newEpsilon >= 0.0); epsilon = newEpsilon; };
    // data 
    std::vector<std::vector<double>> weightMatrix;
    double maxVal;
    constexpr static bool isRestricted = true;
};

 

//using AuctionOracle = AuctionOracleKDTree;
using AuctionOracle = AuctionOracleKDTreeRestricted;
//using AuctionOracle = AuctionOracleLazyHeap;
//using AuctionOracle = AuctionOracleLazyHeapRestricted;
//using AuctionOracle = AuctionOracleRestricted;

class AuctionRunner {
public:
    AuctionRunner(DiagramPointSet& A, DiagramPointSet& B, const double q,  const double _delta, const double _internal_p);
    void setEpsilon(double newVal) { assert(epsilon > 0.0); epsilon = newVal; };
    double getEpsilon(void) const { return epsilon; }
    double getWassersteinDistance(void);
//private:
    // private data
    const size_t numBidders;
    const size_t numItems;
    std::vector<IdxType> allIndices;
    std::vector<IdxType> itemsToBidders;
    std::vector<IdxType> biddersToItems;
    double wassersteinPower;
    std::vector<DiagramPoint> bidders, items;
    double epsilon;
    double delta;
    double internal_p;
    double weightAdjConst;
    double wassersteinDistance;
    std::vector<IdxValPair> bidTable;
    // to get the 2 best items
    std::unique_ptr<AuctionOracle> oracle;
    std::list<size_t> unassignedBidders;
    std::vector< std::list<size_t>::iterator > unassignedBiddersIterators;
    std::vector< short > itemReceivedBidVec;
    std::list<size_t> itemsWithBids;
    // private methods
    void assignGoodToBidder(const IdxType bidderIdx, const IdxType itemsIdx);
    void assignToBestBidder(const IdxType itemsIdx);
    void clearBidTable(void);
    void runAuction(void);
    void runAuctionPhase(void);
    void submitBid(IdxType bidderIdx, const IdxValPair& itemsBidValuePair);
    void flushAssignment(void);

    // for debug only
    void sanityCheck(void);
    void printDebug(void);
    int countUnhappy(void);
    void printMatching(void);
    double getDistanceToQthPowerInternal(void);
    static constexpr double epsilonCommonRatio { 5 }; // next epsilon = current epsilon / epsilonCommonRatio
    static constexpr int maxIterNum { 25 };
};

// get Wasserstein distance between two persistence diagrams
double wassersteinDist(DiagramPointSet& A, DiagramPointSet& B, const double q, const double delta, const double _internal_p = std::numeric_limits<double>::infinity());

#endif
