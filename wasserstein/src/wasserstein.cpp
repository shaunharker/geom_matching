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

#include <math.h>
#include <assert.h>
#include <algorithm>
#include <functional>
#include <iterator>

#include "def_debug.h"
#include "wasserstein.h"

std::ostream& operator<< (std::ostream& output, const DebugOptimalBid& db)
{
    std::cout << "bestItemValue = " << db.bestItemValue;
    std::cout << "; bestItemIdx = " << db.bestItemIdx;
    std::cout << "; secondBestItemValue = " << db.secondBestItemValue;
    std::cout << "; secondBestItemIdx = " << db.secondBestItemIdx;
    return output;
}


double wassersteinDistSlow(DiagramPointSet& A, DiagramPointSet& B, const double q, const double delta)
{
    if(A==B) {
      return 0.0;
    }
    if (q < 1) {
        std::cerr << "Wasserstein distance not defined for q = " << q << ", must be >= 1" << std::endl;
        throw "Bad q in Wasserstein";
    }
    if (delta < 0.0) {
        std::cerr << "Relative error  " << delta << ", must be > 0" << std::endl;
        throw "Bad delta in Wasserstein";
    }
    AuctionRunner auction(A, B, q,  delta);
    return auction.getWassersteinDistance();
}

// *****************************
// AuctionRunner 
// *****************************

AuctionRunner::AuctionRunner(DiagramPointSet& A, DiagramPointSet& B, const double q, const double _delta) :
    numBidders(A.size()),
    numItems(A.size()),
    itemsToBidders(A.size(), -1),
    biddersToItems(A.size(), -1),
    wassersteinPower(q),
    bidTable(A.size(), std::make_pair(-1, std::numeric_limits<double>::lowest()) ),
    itemReceivedBidVec(B.size(), 0 ),
    delta(_delta)
{
    assert(A.size() == B.size());
    items.reserve(numBidders);
    bidders.reserve(numItems);
    IdxType idx { 0 };
    for(const auto& pointA : A) {
        bidders.push_back(pointA);
        allIndices.push_back(idx++);
    }
    for(const auto& pointB : B) {
        items.push_back(pointB);
    }
    oracle = std::unique_ptr<AuctionOracle>(new AuctionOracle(bidders, items, wassersteinPower));
}

void AuctionRunner::assignGoodToBidder(IdxType itemIdx, IdxType bidderIdx)
{
    //sanityCheck();
    IdxType myOldItem = biddersToItems[bidderIdx];
    IdxType currItemOwner = itemsToBidders[itemIdx];

    // set new owner
    biddersToItems[bidderIdx] = itemIdx;
    itemsToBidders[itemIdx] = bidderIdx;


    // remove bidder from the list of unassigned bidders
    unassignedBidders.erase( unassignedBiddersIterators[bidderIdx] );
    assert( 0 <= bidderIdx and bidderIdx < unassignedBiddersIterators.size() );
    unassignedBiddersIterators[bidderIdx] = unassignedBidders.end();

    if (-1 == currItemOwner) {
        // the item we want to assign does not belong to anybody,
        // just free myOldItem, if necessary
        // RE: this cannot be necessary. I submitted the best bid, hence I was
        // an unassigned bidder.
        if (myOldItem != -1) {
            std::cout << "This is not happening" << std::endl;
            assert(false);
            itemsToBidders[myOldItem] = -1;
        }
    } else {
        // the current owner of itemIdx gets my old item (OK if it's -1)
        biddersToItems[currItemOwner] = myOldItem;
        // add the old owner of bids to the list of 
        if ( -1 != myOldItem ) {
            std::cout << "This is not happening" << std::endl;
            assert(false);
            // if I had something, update itemsToBidders, too
            // RE: nonsense: if I had something, I am not unassigned and did not
            // sumbit any bid
            itemsToBidders[myOldItem] = currItemOwner;
        }
        unassignedBidders.push_back(currItemOwner);
        assert( unassignedBiddersIterators[currItemOwner] == unassignedBidders.end() );
        unassignedBiddersIterators[currItemOwner] = std::prev( unassignedBidders.end() );
    }
    //sanityCheck();
}


void AuctionRunner::assignToBestBidder(IdxType itemIdx)
{
    assert( itemIdx >= 0 and itemIdx < static_cast<IdxType>(numItems) );
    assert( bidTable[itemIdx].first != -1);
    IdxValPair bestBid { bidTable[itemIdx] };
    assignGoodToBidder(itemIdx, bestBid.first);
    //std::cout << "About to call setPrice" << std::endl;
    oracle->setPrice(itemIdx,  bestBid.second);
    //dynamic_cast<AuctionOracleKDTree*>(oracle)->setNai
}

void AuctionRunner::clearBidTable(void)
{
    for(auto& itemWithBidIdx : itemsWithBids) {
        itemReceivedBidVec[itemWithBidIdx] = 0;
        bidTable[itemWithBidIdx].first = -1;
        bidTable[itemWithBidIdx].second = std::numeric_limits<double>::lowest();
    }
    itemsWithBids.clear();
}

void AuctionRunner::submitBid(IdxType bidderIdx, const IdxValPair& itemsBidValuePair)
{
    IdxType itemIdx = itemsBidValuePair.first;
    double bidValue = itemsBidValuePair.second;
    assert( itemIdx >= 0 );
    if ( bidTable[itemIdx].second < itemsBidValuePair.second ) {
        bidTable[itemIdx].first = bidderIdx;
        bidTable[itemIdx].second = bidValue;
    }
    if (0 == itemReceivedBidVec[itemIdx]) {
        itemReceivedBidVec[itemIdx] = 1;
        itemsWithBids.push_back(itemIdx);
    }
}

void AuctionRunner::printDebug(void)
{
#ifdef DEBUG_AUCTION
    sanityCheck();
    std::cout << "**********************" << std::endl;
    std::cout << "Current assignment:" << std::endl;
    for(size_t idx = 0; idx < biddersToItems.size(); ++idx) {
        std::cout << idx << " <--> " << biddersToItems[idx] << std::endl;
    }
    std::cout << "Weights: " << std::endl;
    //for(size_t i = 0; i < numBidders; ++i) {
        //for(size_t j = 0; j < numItems; ++j) {
            //std::cout << oracle->weightMatrix[i][j] << " ";
        //}
        //std::cout << std::endl;
    //}
    std::cout << "Prices: " << std::endl;
    for(auto price : oracle->prices) {
        std::cout << price << std::endl;
    }
    //std::cout << "Value matrix: " << std::endl;
    //for(size_t i = 0; i < numBidders; ++i) {
        //for(size_t j = 0; j < numItems; ++j) {
            //std::cout << oracle->weightMatrix[i][j] - oracle->prices[j] << " ";
        //}
        //std::cout << std::endl;
    //}
    std::cout << "**********************" << std::endl;
#endif
}

void AuctionRunner::flushAssignment(void)
{
    for(auto& b2g : biddersToItems) {
        b2g = -1;
    }
    for(auto& g2b : itemsToBidders) {
        g2b = -1;
    }
    //oracle->adjustPrices();
}

void AuctionRunner::runAuction(void)
{
    // relative error 
    // choose some initial epsilon
    oracle->setEpsilon(oracle->maxVal / 4.0);
    assert( oracle->getEpsilon() > 0 );
    int iterNum { 0 };
    bool notDone { false };
    do {
        flushAssignment();
        runAuctionPhase();
        iterNum++;
        //std::cout << "Iteration " << iterNum << " completed. " << std::endl; 
        // result is d^q
        double currentResult = getDistanceToQthPowerInternal();
        double denominator = currentResult - numBidders * oracle->getEpsilon();
        currentResult = pow(currentResult, 1.0 / wassersteinPower);
        //std::cout << "Current result is " << currentResult << std::endl;
        if ( denominator <= 0 ) {
            //std::cout << "Epsilon is too big." << std::endl;
            notDone = true;
        } else {
            denominator = pow(denominator, 1.0 / wassersteinPower);
            double numerator = currentResult - denominator;
            //std::cout << " numerator: " << numerator << " denominator: " << denominator << std::endl;
            //std::cout << " error bound: " << numerator / denominator << std::endl;
            // if relative error is greater than delta, continue
            notDone = ( numerator / denominator > delta );
        }
        // decrease epsilon for the next iteration
        oracle->setEpsilon( oracle->getEpsilon() / epsilonCommonRatio );
        if (iterNum > maxIterNum) {
            std::cerr << "Maximum iteration number exceeded, exiting. Current result is:"; 
            std::cerr << wassersteinDistance << std::endl;
            std::exit(1);
        }
    } while ( notDone );
    //printMatching();
}

void AuctionRunner::runAuctionPhase(void)
{
    //std::cout << "Entered runAuctionPhase" << std::endl;
    int numUnassignedBidders { 0 };

    // at the beginning of a phase all bidders are unassigned
    unassignedBidders.clear();
    unassignedBiddersIterators.clear();
    for(size_t bidderIdx = 0; bidderIdx < allIndices.size(); ++bidderIdx) {
        unassignedBidders.push_back(bidderIdx);
        unassignedBiddersIterators.push_back( std::prev( unassignedBidders.end() ));
    }

    do {
        // get unassigned bidders
        //unassignedBidders.clear();
        //std::copy_if(allIndices.begin(), allIndices.end(),
                     //std::back_inserter(unassignedBidders), 
                     //[this](IdxType i) { return this->biddersToItems[i] == -1; });

        //std::cout << "Number of unassignedBidders: " << unassignedBidders.size() << std::endl;
        ////for(auto ub : unassignedBidders) { std::cout << bidders[ub] << std::endl; }

        //// bidding phase
        //clearBidTable();
        //std::for_each(unassignedBidders.begin(), unassignedBidders.end(),
                      //[ this]
                      //(IdxType bidderIdx) 
                      //{ this->submitBid(bidderIdx, this->oracle->getOptimalBid(bidderIdx)); }
                //);
        ////std::cout << "Bidding phase done" << std::endl;
        // get items with bids
        // todo: do we need this? maybe, better check in assignToBestBidder and
        // loop over all indices?
        
        //std::vector<IdxType> itemsWithBids;
        //itemsWithBids.clear();
        //std::copy_if(allIndices.begin(), allIndices.end(), 
                     //std::back_inserter(itemsWithBids),
                     //[this](IdxType idx) { return this->bidTable[idx].first != -1; } );
        ////std::cout << "Got items with bids " << itemsWithBids.size() << std::endl;
        //// assignment phase
        //std::for_each(itemsWithBids.begin(), itemsWithBids.end(), 
                      //std::bind(&AuctionRunner::assignToBestBidder, this, std::placeholders::_1) );
        //std::cout << "Assignment phase done" << std::endl;

        // bidding phase
        
        //std::cout << "HI" << std::endl;
        clearBidTable();

        for(const auto bidderIdx : unassignedBidders) {
            submitBid(bidderIdx, oracle->getOptimalBid(bidderIdx));
        }

        //std::cout << "Number of unassignedBidders: " << unassignedBidders.size() << std::endl;

        // todo: maintain list of items that received a bid
        for(auto itemIdx : itemsWithBids ) {
            assignToBestBidder(itemIdx);
        }

        //std::cout << "Assignment phase done" << std::endl;
        
        sanityCheck();
        //printDebug();
    } while (unassignedBidders.size() > 0);
    //std::cout << "runAuctionPhase finished" << std::endl;

#ifdef DEBUG_AUCTION
    for(size_t bidderIdx = 0; bidderIdx < numBidders; ++bidderIdx) {
        if ( biddersToItems[bidderIdx] < 0) {
            std::cerr << "After auction terminated bidder " << bidderIdx;
            std::cerr << " has no items assigned" << std::endl;
            throw "Auction did not give a perfect matching";
        }
    }
#endif

}
 
// assertion: the matching must be perfect
double AuctionRunner::getDistanceToQthPowerInternal(void)
{
    sanityCheck();
    double result = 0.0;
    for(size_t bIdx = 0; bIdx < numBidders; ++bIdx) {
        auto pA = bidders[bIdx];
        assert( 0 <= biddersToItems[bIdx] and biddersToItems[bIdx] < items.size() );
        auto pB = items[biddersToItems[bIdx]];
        result += pow(distLInf(pA, pB), wassersteinPower);
    }
    wassersteinDistance = pow(result, 1.0 / wassersteinPower);
    return result;
}

double AuctionRunner::getWassersteinDistance(void)
{
    runAuction();
    return wassersteinDistance;
}

void AuctionRunner::sanityCheck(void)
{
#ifdef DEBUG_AUCTION
    if (biddersToItems.size() != numBidders) {
        std::cerr << "Wrong size of biddersToItems, must be " << numBidders << ", is " << biddersToItems.size() << std::endl;
        throw "Wrong size of biddersToItems";
    }

    if (itemsToBidders.size() != numBidders) {
        std::cerr << "Wrong size of itemsToBidders, must be " << numBidders << ", is " << itemsToBidders.size() << std::endl;
        throw "Wrong size of itemsToBidders";
    }

    for(size_t bidderIdx = 0; bidderIdx < numBidders; ++bidderIdx) {
        if ( biddersToItems[bidderIdx] >= 0) {

            if ( std::count(biddersToItems.begin(),
                        biddersToItems.end(),
                        biddersToItems[bidderIdx]) > 1 ) {
                std::cerr << "Good " << biddersToItems[bidderIdx];
                std::cerr << " appears in biddersToItems more than once" << std::endl;
                throw "Duplicate in biddersToItems";
            }

            if (itemsToBidders.at(biddersToItems[bidderIdx]) != static_cast<int>(bidderIdx)) {
                std::cerr << "Inconsitency: bidderIdx = " << bidderIdx;
                std::cerr << ", itemIdx in biddersToItems = ";
                std::cerr << biddersToItems[bidderIdx];
                std::cerr << ", bidderIdx in itemsToBidders = ";
                std::cerr << itemsToBidders[biddersToItems[bidderIdx]] << std::endl;
                throw "inconsistent mapping";
            }
        }
    }

    for(IdxType itemIdx = 0; itemIdx < static_cast<IdxType>(numBidders); ++itemIdx) {
        if ( itemsToBidders[itemIdx] >= 0) {

            // check for uniqueness
            if ( std::count(itemsToBidders.begin(),
                        itemsToBidders.end(),
                        itemsToBidders[itemIdx]) > 1 ) {
                std::cerr << "Bidder " << itemsToBidders[itemIdx];
                std::cerr << " appears in itemsToBidders more than once" << std::endl;
                throw "Duplicate in itemsToBidders";
            }
            // check for consistency
            if (biddersToItems.at(itemsToBidders[itemIdx]) != static_cast<int>(itemIdx)) {
                std::cerr << "Inconsitency: itemIdx = " << itemIdx;
                std::cerr << ", bidderIdx in itemsToBidders = ";
                std::cerr << itemsToBidders[itemIdx];
                std::cerr << ", itemIdx in biddersToItems= ";
                std::cerr << biddersToItems[itemsToBidders[itemIdx]] << std::endl;
                throw "inconsistent mapping";
            }
        }
    }
#endif
}

void AuctionRunner::printMatching(void)
{
//#ifdef DEBUG_AUCTION
    sanityCheck();
    for(size_t bIdx = 0; bIdx < biddersToItems.size(); ++bIdx) {
        if (biddersToItems[bIdx] >= 0) {
            auto pA = bidders[bIdx];
            auto pB = items[biddersToItems[bIdx]];
            std::cout <<  pA << " <-> " << pB << "+" << pow(distLInf(pA, pB), wassersteinPower) << std::endl;
        } else {
            assert(false);
        }
    }
//#endif
}


AuctionOracleAbstract::AuctionOracleAbstract(const std::vector<DiagramPoint>& _bidders, const std::vector<DiagramPoint>& _items, const double _wassersteinPower) :
    bidders(_bidders),
    items(_items),
    prices(items.size(), 0.0),
    wassersteinPower(_wassersteinPower)
{
}

double AuctionOracleAbstract::getValueForBidder(size_t bidderIdx, size_t itemIdx)
{
    return pow(distLInf(bidders[bidderIdx], items[itemIdx]), wassersteinPower) + prices[itemIdx];
}

// *****************************
// AuctionOracleLazyHeap
// *****************************

AuctionOracleLazyHeap::AuctionOracleLazyHeap(const std::vector<DiagramPoint>& b, 
                                     const std::vector<DiagramPoint>& g, 
                                     double _wassersteinPower) :
    AuctionOracleAbstract(b, g, _wassersteinPower),
    maxVal(std::numeric_limits<double>::min()),
    biddersUpdateMoments(b.size(), 0),
    updateCounter(0)
{
    assert(b.size() == g.size() );
    assert(b.size() > 1);

    weightMatrix.reserve(b.size());
    const double maxDistUpperBound = 3 * getFurthestDistance3Approx(b, g);
    weightAdjConst = pow(maxDistUpperBound, wassersteinPower);
    //std::cout << "3getFurthestDistance3Approx = " << getFurthestDistance3Approx(b, g) << std::endl;
    //std::cout << "in AuctionOracleLazyHeap weightAdjConst = " << weightAdjConst << std::endl;
    // init weight matrix
    for(const auto& pointA : bidders) {
        std::vector<double> weightVec;
        weightVec.clear();
        weightVec.reserve(b.size());
        for(const auto& pointB : items) {
            double val = weightAdjConst - pow(distLInf(pointA, pointB), wassersteinPower);
            weightVec.push_back( val );
            if ( val > maxVal )
                maxVal = val;
        }
        maxVal = weightAdjConst;
        weightMatrix.push_back(weightVec);
    }
    fillInProfitHeap();
    for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        updateList.push_back(std::make_pair(static_cast<IdxType>(itemIdx), 0));
    }
    for(auto updateListIter = updateList.begin(); updateListIter != updateList.end(); ++updateListIter) {
        itemsIterators.push_back(updateListIter);
    }
}

void AuctionOracleLazyHeap::updateQueueForBidder(IdxType bidderIdx)
{
    assert(0 <= bidderIdx and bidderIdx < static_cast<int>(bidders.size()));
    assert(bidderIdx < static_cast<int>(biddersUpdateMoments.size()));

    int bidderLastUpdateTime = biddersUpdateMoments[bidderIdx];
    auto iter = updateList.begin();
    while (iter != updateList.end() and iter->second >= bidderLastUpdateTime) {
        IdxType itemIdx = iter->first;
        IdxValPair newVal { itemIdx, weightMatrix[bidderIdx][itemIdx] - prices[itemIdx]};
        // to-do: change indexing of profitHeapHandles
        profitHeap[bidderIdx]->decrease(profitHeapHandles[bidderIdx][itemIdx], newVal);
        iter++;
    }
    biddersUpdateMoments[bidderIdx] = updateCounter;
}

void AuctionOracleLazyHeap::fillInProfitHeap(void)
{
    for(size_t bidderIdx = 0; bidderIdx < bidders.size(); ++bidderIdx) {
        profitHeap.push_back( new AuctionHeap() );
        std::vector<AuctionHeap::handle_type> handlesVec;
        profitHeapHandles.push_back(handlesVec);
        for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
            IdxValPair vp { itemIdx, weightMatrix[bidderIdx][itemIdx] - prices[itemIdx] };
            profitHeapHandles[bidderIdx].push_back(  profitHeap[bidderIdx]->push(vp) );
        }
    }
}

AuctionOracleLazyHeap::~AuctionOracleLazyHeap()
{
    for(auto h : profitHeap) {
        delete h;
    }
}

void AuctionOracleLazyHeap::setPrice(IdxType itemIdx, double newPrice)
{
    assert( prices.at(itemIdx) < newPrice );
#ifdef DEBUG_AUCTION
    std::cout << "price incremented by " <<  prices.at(itemIdx) - newPrice << std::endl;
#endif
    prices[itemIdx] = newPrice;
    // lazy: record the moment we updated the price of the items,
    // do not update queues.
    // 1. move the items with updated price to the front of the updateList,
    updateList.splice(updateList.begin(), updateList, itemsIterators[itemIdx]);
    // 2. record the moment we updated the price and increase the counter
    updateList.front().second = updateCounter++;
}

// subtract min. price from all prices
void AuctionOracleLazyHeap::adjustPrices(void)
{
    double minPrice = *(std::min_element(prices.begin(), prices.end()));
    std::transform(prices.begin(), prices.end(), prices.begin(), [minPrice](double a) { return a - minPrice; });
}

DebugOptimalBid AuctionOracleLazyHeap::getOptimalBidDebug(IdxType bidderIdx)
{
    assert(bidderIdx >=0 and bidderIdx < static_cast<IdxType>(bidders.size()) );
    assert(profitHeap.at(bidderIdx) != nullptr);
    assert(profitHeap[bidderIdx]->size() >= 2);

    updateQueueForBidder(bidderIdx);
    DebugOptimalBid result;
    
    auto pHeap = profitHeap[bidderIdx];
    auto topIter = pHeap->ordered_begin(); 
    result.bestItemIdx = topIter->first;
    result.bestItemValue = topIter->second;
    ++topIter; // now points to the second-best items
    result.secondBestItemValue = topIter->second;
    result.secondBestItemIdx = topIter->first;

    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemValue = " << bestItemValue << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << topIter->first << "; secondBestValue = " << secondBestItemValue << "; secondBestPrice = " << prices[topIter->first] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;
    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemsDist= " << (weightAdjConst -  bestItemValue) << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << topIter->first << "; secondBestDist= " << (weightAdjConst - secondBestItemValue) << "; secondBestPrice = " << prices[topIter->first] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;

    return result;
}

IdxValPair AuctionOracleLazyHeap::getOptimalBid(const IdxType bidderIdx) 
{
    assert(bidderIdx >=0 and bidderIdx < static_cast<IdxType>(bidders.size()) );
    assert(profitHeap.at(bidderIdx) != nullptr);
    assert(profitHeap[bidderIdx]->size() >= 2);

    updateQueueForBidder(bidderIdx);
    
    auto pHeap = profitHeap[bidderIdx];
    auto topIter = pHeap->ordered_begin(); 
    IdxType bestItemIdx = topIter->first;
    double bestItemValue = topIter->second;
    ++topIter; // now points to the second-best items
    double secondBestItemValue = topIter->second;

    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemValue = " << bestItemValue << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << topIter->first << "; secondBestValue = " << secondBestItemValue << "; secondBestPrice = " << prices[topIter->first] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;
    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemsDist= " << (weightAdjConst -  bestItemValue) << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << topIter->first << "; secondBestDist= " << (weightAdjConst - secondBestItemValue) << "; secondBestPrice = " << prices[topIter->first] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;

    // bid value: price + value difference + epsilon
    return std::make_pair(bestItemIdx, 
                          prices[bestItemIdx] + 
                          ( bestItemValue - secondBestItemValue ) +
                          epsilon );
}

// *****************************
// AuctionOracleLazyHeapRestricted
// *****************************

AuctionOracleLazyHeapRestricted::AuctionOracleLazyHeapRestricted(const std::vector<DiagramPoint>& b, 
                                     const std::vector<DiagramPoint>& g, 
                                     double _wassersteinPower) :
    AuctionOracleAbstract(b, g, _wassersteinPower),
    maxVal(std::numeric_limits<double>::min()),
    biddersUpdateMoments(b.size(), 0),
    updateCounter(0),
    biddersToProjItems(bidders.size(), std::numeric_limits<size_t>::max()),
    heapHandlesIndices(items.size(), std::numeric_limits<size_t>::max()),
    bestDiagonalItemsComputed(false)
{
    assert(b.size() == g.size() );
    assert(b.size() > 1);

    weightMatrix.reserve(b.size());
    const double maxDistUpperBound = 3 * getFurthestDistance3Approx(b, g);
    weightAdjConst = pow(maxDistUpperBound, wassersteinPower);
    //std::cout << "3getFurthestDistance3Approx = " << getFurthestDistance3Approx(b, g) << std::endl;
    //std::cout << "in AuctionOracleLazyHeapRestricted weightAdjConst = " << weightAdjConst << std::endl;
    // init weight matrix
    for(const auto& pointA : bidders) {
        std::vector<double> weightVec;
        weightVec.clear();
        weightVec.reserve(b.size());
        for(const auto& pointB : items) {
            double val = pow(distLInf(pointA, pointB), wassersteinPower);
            weightVec.push_back( val );
            if ( val > maxVal )
                maxVal = val;
        }
        maxVal = weightAdjConst;
        weightMatrix.push_back(weightVec);
    }
    fillInProfitHeap();
    for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        updateList.push_back(std::make_pair(static_cast<IdxType>(itemIdx), 0));
    }
    for(auto updateListIter = updateList.begin(); updateListIter != updateList.end(); ++updateListIter) {
        itemsIterators.push_back(updateListIter);
    }

    size_t handleIdx {0};
    for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        if (items[itemIdx].isDiagonal() ) {
            heapHandlesIndices[itemIdx] = handleIdx++;
            diagHeapHandles.push_back(diagItemsHeap.push(std::make_pair(itemIdx, 0)));
        }
    }
     // todo: this must be done in readFiles procedure
    //std::cout << "started wasting time..." << std::endl;
    for(size_t bidderIdx = 0; bidderIdx < bidders.size(); ++bidderIdx) {
        for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
            auto bidder = bidders[bidderIdx];
            auto item = items[itemIdx];
            if (bidder.id == item.projId) {
                biddersToProjItems[bidderIdx] = itemIdx;
            }
        }
    }
    //std::cout << "stopped wasting time." << std::endl;
}

void AuctionOracleLazyHeapRestricted::updateQueueForBidder(IdxType bidderIdx)
{
    assert(0 <= bidderIdx and bidderIdx < static_cast<int>(bidders.size()));
    assert(bidderIdx < static_cast<int>(biddersUpdateMoments.size()));
    assert(profitHeap[bidderIdx] != nullptr );

    int bidderLastUpdateTime = biddersUpdateMoments[bidderIdx];
    auto iter = updateList.begin();
    while (iter != updateList.end() and iter->second >= bidderLastUpdateTime) {
        IdxType itemIdx = iter->first;
        size_t handleIdx = itemsIndicesForHeapHandles[bidderIdx][itemIdx];
        if (handleIdx  < items.size() ) {
            IdxValPair newVal { itemIdx, weightMatrix[bidderIdx][itemIdx] + prices[itemIdx]};
            // to-do: change indexing of profitHeapHandles
            profitHeap[bidderIdx]->decrease(profitHeapHandles[bidderIdx][handleIdx], newVal);
        }
        iter++;
    }
    biddersUpdateMoments[bidderIdx] = updateCounter;
}

void AuctionOracleLazyHeapRestricted::fillInProfitHeap(void)
{
    for(size_t bidderIdx = 0; bidderIdx < bidders.size(); ++bidderIdx) {
        DiagramPoint bidder { bidders[bidderIdx] };
        // no heap for diagonal bidders
        if ( bidder.isDiagonal() ) {
            profitHeap.push_back( nullptr );
            profitHeapHandles.push_back(std::vector<AuctionHeap::handle_type>());
            itemsIndicesForHeapHandles.push_back( std::vector<size_t>() );
            continue;
        } else {
            profitHeap.push_back( new AuctionHeap() );
            assert( profitHeap.at(bidderIdx) != nullptr );
            itemsIndicesForHeapHandles.push_back( std::vector<size_t>(items.size(), std::numeric_limits<size_t>::max() ) );

            std::vector<AuctionHeap::handle_type> handlesVec;
            profitHeapHandles.push_back(handlesVec);
            size_t handleIdx { 0 };
            for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
                assert( itemsIndicesForHeapHandles.at(bidderIdx).at(itemIdx) > 0 );
                DiagramPoint item { items[itemIdx] };
                if ( item.isNormal() ) {
                    // item can be assigned to bidder, store in heap 
                    IdxValPair vp { itemIdx, weightMatrix[bidderIdx][itemIdx] + prices[itemIdx] };
                    profitHeapHandles[bidderIdx].push_back(  profitHeap[bidderIdx]->push(vp) );
                    // keep corresponding index in itemsIndicesForHeapHandles
                    itemsIndicesForHeapHandles[bidderIdx][itemIdx] = handleIdx++;
                }
            }
        }
    }
}

AuctionOracleLazyHeapRestricted::~AuctionOracleLazyHeapRestricted()
{
    for(auto h : profitHeap) {
        delete h;
    }
}

void AuctionOracleLazyHeapRestricted::setPrice(IdxType itemIdx, double newPrice)
{
    assert( prices.at(itemIdx) < newPrice );
#ifdef DEBUG_AUCTION
    std::cout << "price incremented by " <<  prices.at(itemIdx) - newPrice << std::endl;
#endif
    prices[itemIdx] = newPrice;
    if (items[itemIdx].isNormal() ) {
        // lazy: record the moment we updated the price of the items,
        // do not update queues.
        // 1. move the items with updated price to the front of the updateList,
        updateList.splice(updateList.begin(), updateList, itemsIterators[itemIdx]);
        // 2. record the moment we updated the price and increase the counter
        updateList.front().second = updateCounter++;
    } else {
        // diagonal items are stored in one heap
        diagItemsHeap.decrease(diagHeapHandles[heapHandlesIndices[itemIdx]], std::make_pair(itemIdx, newPrice));
        bestDiagonalItemsComputed = false;
    }
}

// subtract min. price from all prices
void AuctionOracleLazyHeapRestricted::adjustPrices(void)
{
}

DebugOptimalBid AuctionOracleLazyHeapRestricted::getOptimalBidDebug(IdxType bidderIdx)
{
    DebugOptimalBid result;
    assert(bidderIdx >=0 and bidderIdx < static_cast<IdxType>(bidders.size()) );

    DiagramPoint bidder = bidders[bidderIdx];
    std::vector<IdxValPair> candItems;
    // corresponding point is always considered as a candidate

    size_t projItemIdx = biddersToProjItems[bidderIdx];
    assert( 0 <= projItemIdx and projItemIdx < items.size() );
    DiagramPoint projItem = items[projItemIdx];
    assert(projItem.type != bidder.type);
    assert(projItem.projId == bidder.id);
    assert(projItem.id == bidder.projId);
    // todo: store precomputed distance?
    double projItemValue = pow(distLInf(bidder, projItem), wassersteinPower) + prices[projItemIdx];
    candItems.push_back( std::make_pair(projItemIdx, projItemValue) );
 
    if (bidder.isNormal()) {
        assert(profitHeap.at(bidderIdx) != nullptr);
        assert(profitHeap[bidderIdx]->size() >= 2);
        updateQueueForBidder(bidderIdx);
        auto pHeap = profitHeap[bidderIdx];
        assert( pHeap != nullptr );
        auto topIter = pHeap->ordered_begin(); 
        candItems.push_back( *topIter );
        ++topIter; // now points to the second-best items
        candItems.push_back( *topIter );
        std::sort(candItems.begin(), candItems.end(), CompPairsBySecondStruct());
        assert(candItems[1].second >= candItems[0].second);
    } else {
        // for diagonal bidder the only normal point has already been added
        // the other 2 candidates are diagonal items only, get from the heap
        // with prices
        assert(diagItemsHeap.size() > 1);
        auto topDiagIter = diagItemsHeap.ordered_begin();
        auto topDiag1 = *topDiagIter++;
        auto topDiag2 = *topDiagIter;
        candItems.push_back(topDiag1);
        candItems.push_back(topDiag2);
        std::sort(candItems.begin(), candItems.end(), CompPairsBySecondStruct());
        assert(candItems.size() == 3);
        assert(candItems[2].second >= candItems[1].second);
        assert(candItems[1].second >= candItems[0].second);
    }
    
    result.bestItemIdx = candItems[0].first;
    result.secondBestItemIdx = candItems[1].first;
    result.bestItemValue = candItems[0].second;
    result.secondBestItemValue = candItems[1].second;

    // checking code

    //DebugOptimalBid debugMyResult(result);
    //DebugOptimalBid debugNaiveResult;
    //debugNaiveResult.bestItemValue = 1e20;
    //debugNaiveResult.secondBestItemValue = 1e20;
    //double currItemValue;
    //for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        //if ( bidders[bidderIdx].type != items[itemIdx].type and
                //bidders[bidderIdx].projId != items[itemIdx].id)
            //continue;

        //currItemValue = pow(distLInf(bidders[bidderIdx], items[itemIdx]), wassersteinPower) + prices[itemIdx];
        //if (currItemValue < debugNaiveResult.bestItemValue) {
            //debugNaiveResult.bestItemValue = currItemValue;
            //debugNaiveResult.bestItemIdx  = itemIdx;
        //}
    //}

    //for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        //if (itemIdx == debugNaiveResult.bestItemIdx) {
            //continue;
        //}
        //if ( bidders[bidderIdx].type != items[itemIdx].type and
                //bidders[bidderIdx].projId != items[itemIdx].id)
            //continue;

        //currItemValue = pow(distLInf(bidders[bidderIdx], items[itemIdx]), wassersteinPower) + prices[itemIdx];
        //if (currItemValue < debugNaiveResult.secondBestItemValue) {
            //debugNaiveResult.secondBestItemValue = currItemValue;
            //debugNaiveResult.secondBestItemIdx = itemIdx;
        //}
    //}
    ////std::cout << "got naive result" << std::endl;

    //if ( fabs( debugMyResult.bestItemValue - debugNaiveResult.bestItemValue ) > 1e-6 or
            //fabs( debugNaiveResult.secondBestItemValue - debugMyResult.secondBestItemValue) > 1e-6 ) {
        //std::cerr << "bidderIdx = " << bidderIdx << "; ";
        //std::cerr << bidders[bidderIdx] << std::endl;
        //for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
            //std::cout << itemIdx << ": " << items[itemIdx] << "; price = " << prices[itemIdx] << std::endl;
        //}
        //std::cerr << "debugMyResult: " << debugMyResult << std::endl;
        //std::cerr << "debugNaiveResult: " << debugNaiveResult << std::endl;
        //auto pHeap = profitHeap[bidderIdx];
        //assert( pHeap != nullptr );
        //for(auto topIter = pHeap->ordered_begin(); topIter != pHeap->ordered_end(); ++topIter) {
            //std::cerr << "in heap: " << topIter->first << ": " << topIter->second << "; real value = " << distLInf(bidder, items[topIter->first]) + prices[topIter->first] << std::endl;
        //}
        //for(auto ci : candItems) {
            //std::cout << "ci.idx = " << ci.first << ", value = " << ci.second << std::endl;
        //}

        ////std::cerr << "twoBestItems: " << twoBestItems[0].d << " " << twoBestItems[1].d << std::endl;
        //assert(false);
    //}


    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemValue = " << bestItemValue << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << topIter->first << "; secondBestValue = " << secondBestItemValue << "; secondBestPrice = " << prices[topIter->first] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;
    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemsDist= " << (weightAdjConst -  bestItemValue) << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << topIter->first << "; secondBestDist= " << (weightAdjConst - secondBestItemValue) << "; secondBestPrice = " << prices[topIter->first] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;

    return result;
}

IdxValPair AuctionOracleLazyHeapRestricted::getOptimalBid(const IdxType bidderIdx) 
{
    IdxType bestItemIdx;
    IdxType secondBestItemIdx;
    double bestItemValue;
    double secondBestItemValue;

    auto& bidder = bidders[bidderIdx];
    IdxType projItemIdx = biddersToProjItems[bidderIdx];
    assert( 0 <= projItemIdx and projItemIdx < items.size() );
    DiagramPoint projItem = items[projItemIdx];
    assert(projItem.type != bidder.type);
    assert(projItem.projId == bidder.id);
    assert(projItem.id == bidder.projId);
    // todo: store precomputed distance?
    double projItemValue = pow(distLInf(bidder, projItem), wassersteinPower) + prices[projItemIdx];
   
    if (bidder.isDiagonal()) {
        // for diagonal bidder the only normal point has already been added
        // the other 2 candidates are diagonal items only, get from the heap
        // with prices
        assert(diagItemsHeap.size() > 1);
        if (!bestDiagonalItemsComputed) {
            auto topDiagIter = diagItemsHeap.ordered_begin();
            bestDiagonalItemIdx = topDiagIter->first;
            bestDiagonalItemValue = topDiagIter->second;
            topDiagIter++;
            secondBestDiagonalItemIdx = topDiagIter->first;
            secondBestDiagonalItemValue = topDiagIter->second;
            bestDiagonalItemsComputed = true;
        }

        if ( projItemValue < bestDiagonalItemValue) {
            bestItemIdx = projItemIdx;
            bestItemValue = projItemValue;
            secondBestItemValue = bestDiagonalItemValue;
            secondBestItemIdx = bestDiagonalItemIdx;
        } else if (projItemValue < secondBestDiagonalItemValue) {
            bestItemIdx = bestDiagonalItemIdx;
            bestItemValue = bestDiagonalItemValue;
            secondBestItemValue = projItemValue;
            secondBestItemIdx = projItemIdx;
        } else {
            bestItemIdx = bestDiagonalItemIdx;
            bestItemValue = bestDiagonalItemValue;
            secondBestItemValue = secondBestDiagonalItemValue;
            secondBestItemIdx = secondBestDiagonalItemIdx;
        }
    } else {
        // for normal bidder get 2 best items among non-diagonal (=normal) points 
        // from the corresponding heap 
        assert(diagItemsHeap.size() > 1);
        updateQueueForBidder(bidderIdx);
        auto topNormIter = profitHeap[bidderIdx]->ordered_begin();
        IdxType bestNormalItemIdx { topNormIter->first };
        double bestNormalItemValue { topNormIter->second };
        topNormIter++;
        double secondBestNormalItemValue { topNormIter->second };
        IdxType secondBestNormalItemIdx { topNormIter->first };
 
        if ( projItemValue < bestNormalItemValue) {
            bestItemIdx = projItemIdx;
            bestItemValue = projItemValue;
            secondBestItemValue = bestNormalItemValue;
            secondBestItemIdx = bestNormalItemIdx;
        } else if (projItemValue < secondBestNormalItemValue) {
            bestItemIdx = bestNormalItemIdx;
            bestItemValue = bestNormalItemValue;
            secondBestItemValue = projItemValue;
            secondBestItemIdx = projItemIdx;
        } else {
            bestItemIdx = bestNormalItemIdx;
            bestItemValue = bestNormalItemValue;
            secondBestItemValue = secondBestNormalItemValue;
            secondBestItemIdx = secondBestNormalItemIdx;
        }
    }

    IdxValPair result;

    assert( secondBestItemValue >= bestItemValue );

    result.first = bestItemIdx;
    result.second = ( secondBestItemValue - bestItemValue ) + prices[bestItemIdx] + epsilon;


    // checking code

    //DebugOptimalBid debugMyResult;
    //debugMyResult.bestItemIdx = bestItemIdx;
    //debugMyResult.bestItemValue = bestItemValue;
    //debugMyResult.secondBestItemIdx = secondBestItemIdx;
    //debugMyResult.secondBestItemValue = secondBestItemValue;
    //DebugOptimalBid debugNaiveResult;
    //debugNaiveResult.bestItemValue = 1e20;
    //debugNaiveResult.secondBestItemValue = 1e20;
    //double currItemValue;
    //for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        //if ( bidders[bidderIdx].type != items[itemIdx].type and
                //bidders[bidderIdx].projId != items[itemIdx].id)
            //continue;

        //currItemValue = pow(distLInf(bidders[bidderIdx], items[itemIdx]), wassersteinPower) + prices[itemIdx];
        //if (currItemValue < debugNaiveResult.bestItemValue) {
            //debugNaiveResult.bestItemValue = currItemValue;
            //debugNaiveResult.bestItemIdx  = itemIdx;
        //}
    //}

    //for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        //if (itemIdx == debugNaiveResult.bestItemIdx) {
            //continue;
        //}
        //if ( bidders[bidderIdx].type != items[itemIdx].type and
                //bidders[bidderIdx].projId != items[itemIdx].id)
            //continue;

        //currItemValue = pow(distLInf(bidders[bidderIdx], items[itemIdx]), wassersteinPower) + prices[itemIdx];
        //if (currItemValue < debugNaiveResult.secondBestItemValue) {
            //debugNaiveResult.secondBestItemValue = currItemValue;
            //debugNaiveResult.secondBestItemIdx = itemIdx;
        //}
    //}
    ////std::cout << "got naive result" << std::endl;

    //if ( fabs( debugMyResult.bestItemValue - debugNaiveResult.bestItemValue ) > 1e-6 or
            //fabs( debugNaiveResult.secondBestItemValue - debugMyResult.secondBestItemValue) > 1e-6 ) {
        //std::cerr << "bidderIdx = " << bidderIdx << "; ";
        //std::cerr << bidders[bidderIdx] << std::endl;
        //for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
            //std::cout << itemIdx << ": " << items[itemIdx] << "; price = " << prices[itemIdx] << std::endl;
        //}
        //std::cerr << "debugMyResult: " << debugMyResult << std::endl;
        //std::cerr << "debugNaiveResult: " << debugNaiveResult << std::endl;
        //auto pHeap = profitHeap[bidderIdx];
        //if ( pHeap != nullptr ) {
            //for(auto topIter = pHeap->ordered_begin(); topIter != pHeap->ordered_end(); ++topIter) {
                //std::cerr << "in heap: " << topIter->first << ": " << topIter->second << "; real value = " << distLInf(bidder, items[topIter->first]) + prices[topIter->first] << std::endl;
            //}
        //}
        ////for(auto ci : candItems) {
            ////std::cout << "ci.idx = " << ci.first << ", value = " << ci.second << std::endl;
        ////}

        ////std::cerr << "twoBestItems: " << twoBestItems[0].d << " " << twoBestItems[1].d << std::endl;
        //assert(false);
    // }
    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemValue = " << bestItemValue << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << topIter->first << "; secondBestValue = " << secondBestItemValue << "; secondBestPrice = " << prices[topIter->first] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;
    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemsDist= " << (weightAdjConst -  bestItemValue) << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << topIter->first << "; secondBestDist= " << (weightAdjConst - secondBestItemValue) << "; secondBestPrice = " << prices[topIter->first] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;

    return result;
}


// *****************************
// AuctionOracleKDTree
// *****************************

AuctionOracleKDTree::AuctionOracleKDTree(const std::vector<DiagramPoint>& _bidders, 
        const std::vector<DiagramPoint>& _items, 
        double _wassersteinPower) :
    AuctionOracleAbstract(_bidders, _items, _wassersteinPower),
    heapHandlesIndices(items.size(), std::numeric_limits<size_t>::max()),
    kdtreeItems(items.size(), std::numeric_limits<size_t>::max())
{
    //assert(wassersteinPower == 1.0); // temporarily, to-do: update dnn to search with any q
    size_t dnnItemIdx { 0 };
    size_t trueIdx { 0 };
    dnnPoints.clear();
    // store normal items in kd-tree
    for(const auto& g : items) {
        if (g.isNormal()) {
            kdtreeItems[trueIdx] = dnnItemIdx;
            // index of items is id of dnn-point
            DnnPoint p(trueIdx);
            p[0] = g.x;
            p[1] = g.y;
            dnnPoints.push_back(p);
            assert(dnnItemIdx == dnnPoints.size() - 1);
            dnnItemIdx++;
        }
        trueIdx++;
    }

    assert(dnnPoints.size() < items.size() );
    for(size_t i = 0; i < dnnPoints.size(); ++i) {
        dnnPointHandles.push_back(&dnnPoints[i]);
    }
    DnnTraits traits;
    //std::cout << "kdtree: " << dnnPointHandles.size() << " points" << std::endl;
    kdtree = new dnn::KDTree<DnnTraits>(traits, dnnPointHandles, wassersteinPower);

    size_t dnnItemIdxAll { 0 };
    dnnPointsAll.clear();
    // store all items in kd-tree
    for(const auto& g : items) {
        DnnPoint p(dnnItemIdxAll++);
        p[0] = g.getRealX();
        p[1] = g.getRealY();
        //std::cout << "to all tree: " << p[0] << ", " << p[1] << std::endl;
        dnnPointsAll.push_back(p);
        assert(dnnItemIdxAll == dnnPointsAll.size());
    }

    for(size_t i = 0; i < dnnPointsAll.size(); ++i) {
        dnnPointHandlesAll.push_back(&dnnPointsAll[i]);
    }
    //std::cout << "kdtreeAll: " << dnnPointHandlesAll.size() << " points" << std::endl;
    kdtreeAll = new dnn::KDTree<DnnTraits>(traits, dnnPointHandlesAll, wassersteinPower);
    
    size_t handleIdx {0};
    for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        if (items[itemIdx].isDiagonal() ) {
            heapHandlesIndices[itemIdx] = handleIdx++;
            diagHeapHandles.push_back(diagItemsHeap.push(std::make_pair(itemIdx, 0)));
        }
    }
    //to-do: remove maxVal from 
    //std::cout << "3getFurthestDistance3Approx = " << getFurthestDistance3Approx(_bidders, _items) << std::endl;
    maxVal = 3*getFurthestDistance3Approx(_bidders, _items);
    maxVal = pow(maxVal, wassersteinPower);
    weightAdjConst = maxVal;
    //std::cout << "AuctionOracleKDTree: weightAdjConst = " << weightAdjConst << std::endl;
    //std::cout << "AuctionOracleKDTree constructor done" << std::endl;
    // for debug
}

DebugOptimalBid AuctionOracleKDTree::getOptimalBidDebug(IdxType bidderIdx)
{
    DebugOptimalBid result;
    DiagramPoint bidder = bidders[bidderIdx];
    DnnPoint bidderDnn;
    bidderDnn[0] = bidder.getRealX();
    bidderDnn[1] = bidder.getRealY();

    //std::cout << "bidder.x = " << bidderDnn[0] << std::endl;
    //std::cout << "bidder.y = " << bidderDnn[1] << std::endl;

    std::vector<IdxValPair> candItems;
    

    if ( bidder.isDiagonal() ) {
        // 
        auto twoBestItems = kdtree->findK(bidderDnn, 2);
        //std::cout << "twoBestItems for non-diag: " << twoBestItems[0].d << " " << twoBestItems[1].d << std::endl;
        candItems.push_back( std::make_pair(twoBestItems[0].p->id(), twoBestItems[0].d) );
        candItems.push_back( std::make_pair(twoBestItems[1].p->id(), twoBestItems[1].d) );
        assert(diagItemsHeap.size() > 1);
        auto topDiagIter = diagItemsHeap.ordered_begin();
        auto topDiag1 = *topDiagIter++;
        auto topDiag2 = *topDiagIter;
        candItems.push_back(topDiag1);
        candItems.push_back(topDiag2);
        assert(candItems.size() == 4);
        std::sort(candItems.begin(), candItems.end(), CompPairsBySecondStruct());
        assert(candItems[3].second >= candItems[2].second);
        assert(candItems[2].second >= candItems[1].second);
        assert(candItems[1].second >= candItems[0].second);
    } else {
        auto twoBestItems = kdtreeAll->findK(bidderDnn, 2);
        //std::cout << "twoBestItems for all: " << twoBestItems[0].d << " " << twoBestItems[1].d << std::endl;
        candItems.push_back( std::make_pair(twoBestItems[0].p->id(), twoBestItems[0].d) );
        candItems.push_back( std::make_pair(twoBestItems[1].p->id(), twoBestItems[1].d) );
        //size_t projItemIdx { biddersProjIndices.at(bidderIdx) };
        //assert(items[projItemIdx].projId == bidder.id);
        //double projItemValue { pow(distLInf(bidder, items[projItemIdx]), wassersteinPower) + prices.at(projItemIdx) };
        //candItems.push_back( std::make_pair(projItemIdx, projItemValue) );
        assert(candItems.size() == 2);
        assert(candItems[1].second >= candItems[0].second);
    }

    result.bestItemIdx = candItems[0].first;
    result.secondBestItemIdx = candItems[1].first;
    result.bestItemValue = candItems[0].second;
    result.secondBestItemValue = candItems[1].second;
    //std::cout << "got result: " << result << std::endl;
    //double bestItemsPrice = prices[bestItemIdx];
    //if (items[result.bestItemIdx].type == DiagramPoint::DIAG) {
        //double bestItemValue1 = pow( distLInf(bidder, items[result.bestItemIdx]), q) + prices[result.bestItemIdx];
        //if ( fabs(result.bestItemValue - bestItemValue1) > 1e-6 ) {
            //std::cerr << "XXX: " << result.bestItemValue << " vs " << bestItemValue1 << std::endl;
            //result.bestItemValue = bestItemValue1;
        //}

    //}


    // checking code
    /*
    
    DebugOptimalBid debugMyResult(result);
    DebugOptimalBid debugNaiveResult;
    debugNaiveResult.bestItemValue = 1e20;
    debugNaiveResult.secondBestItemValue = 1e20;
    double currItemValue;
    for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        //if ( bidders[bidderIdx].type == DiagramPoint::NORMAL and
                //items[itemIdx].type == DiagramPoint::DIAG and
                //bidders[bidderIdx].projId != items[itemIdx].id)
            //continue;

        currItemValue = pow(distLInf(bidders[bidderIdx], items[itemIdx]), wassersteinPower) + prices[itemIdx];
        if (currItemValue < debugNaiveResult.bestItemValue) {
            debugNaiveResult.bestItemValue = currItemValue;
            debugNaiveResult.bestItemIdx  = itemIdx;
        }
    }

    for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        if (itemIdx == debugNaiveResult.bestItemIdx) {
            continue;
        }
        currItemValue = pow(distLInf(bidders[bidderIdx], items[itemIdx]), wassersteinPower) + prices[itemIdx];
        if (currItemValue < debugNaiveResult.secondBestItemValue) {
            debugNaiveResult.secondBestItemValue = currItemValue;
            debugNaiveResult.secondBestItemIdx = itemIdx;
        }
    }
    //std::cout << "got naive result" << std::endl;

    if ( fabs( debugMyResult.bestItemValue - debugNaiveResult.bestItemValue ) > 1e-6 or
            fabs( debugNaiveResult.secondBestItemValue - debugMyResult.secondBestItemValue) > 1e-6 ) {
        kdtreeAll->printWeights();
        std::cerr << "bidderIdx = " << bidderIdx << "; ";
        std::cerr << bidders[bidderIdx] << std::endl;
        for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
            std::cout << itemIdx << ": " << items[itemIdx] << "; price = " << prices[itemIdx] << std::endl;
        }
        std::cerr << "debugMyResult: " << debugMyResult << std::endl;
        std::cerr << "debugNaiveResult: " << debugNaiveResult << std::endl;
        //std::cerr << "twoBestItems: " << twoBestItems[0].d << " " << twoBestItems[1].d << std::endl;
        assert(false);
    }
    //std::cout << "returning" << std::endl;

    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemValue = " << bestItemValue << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << secondBestItemIdx << "; secondBestValue = " << secondBestItemValue << "; secondBestPrice = " << prices[secondBestItemIdx] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;
    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemsDist= " << (weightAdjConst - bestItemValue) << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << secondBestItemIdx << "; secondBestDist= " << (weightAdjConst - secondBestItemValue) << "; secondBestPrice = " << prices[secondBestItemIdx] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;
    */

    return result;
}

IdxValPair AuctionOracleKDTree::getOptimalBid(IdxType bidderIdx)
{
    IdxValPair result;
    DebugOptimalBid debugMyResult = getOptimalBidDebug(bidderIdx);
    result.first = debugMyResult.bestItemIdx;
    result.second = ( debugMyResult.secondBestItemValue - debugMyResult.bestItemValue ) + prices[debugMyResult.bestItemIdx] + epsilon;
    return result;
}
/*
a_{ij} = d_{ij} 
value_{ij} = a_{ij} + price_j
*/
void AuctionOracleKDTree::setPrice(IdxType itemIdx, double newPrice)
{
    assert(prices.size() == items.size());
    assert( 0 < diagHeapHandles.size() and diagHeapHandles.size() <= items.size());
    assert(newPrice > prices.at(itemIdx));
    prices[itemIdx] = newPrice;
    if ( items[itemIdx].isNormal() ) {
        //std::cout << "before increasing weight in kdtree " << std::endl;
        //std::cout << kdtreeItems.at(itemIdx) << std::endl;
        assert(0 <= itemIdx and itemIdx < kdtreeItems.size());
        assert(0 <= kdtreeItems[itemIdx] and kdtreeItems[itemIdx] < dnnPointHandles.size());
        kdtree->increase_weight( dnnPointHandles[kdtreeItems[itemIdx]], newPrice);
        kdtreeAll->increase_weight( dnnPointHandlesAll[itemIdx], newPrice);
        //std::cout << "after increasing weight in kdtree" << std::endl;
    } else {
        //std::cout << "before decreasing weight in heap" << std::endl;
        //std::cout << "diagHeapHandles.size = " << diagHeapHandles.size() << std::endl;
        kdtreeAll->increase_weight( dnnPointHandlesAll[itemIdx], newPrice);
        //std::cout << "after increasing weight in kdtree" << std::endl;
        assert(diagHeapHandles.size() > heapHandlesIndices.at(itemIdx));
        diagItemsHeap.decrease(diagHeapHandles[heapHandlesIndices[itemIdx]], std::make_pair(itemIdx, newPrice));
    } 
}

void AuctionOracleKDTree::adjustPrices(void)
{
}

AuctionOracleKDTree::~AuctionOracleKDTree()
{
    delete kdtree;
    delete kdtreeAll;
}

void AuctionOracleKDTree::setEpsilon(double newVal) 
{
    assert(newVal >= 0.0);
    epsilon = newVal;
}

// *****************************
// AuctionOracleRestricted
// *****************************
AuctionOracleRestricted::AuctionOracleRestricted(const std::vector<DiagramPoint>& b, 
                                         const std::vector<DiagramPoint>& g, 
                                         double _wassersteinPower) :
    AuctionOracleAbstract(b, g, _wassersteinPower),
    maxVal(0.0)
{
    assert(b.size() == g.size() );
    assert(b.size() > 1);

    weightMatrix.reserve(b.size());
    for(const auto& pointA : bidders) {
        std::vector<double> weightVec;
        weightVec.clear();
        weightVec.reserve(b.size());
        for(const auto& pointB : items) {
            double val = pow(distLInf(pointA, pointB), wassersteinPower);
            if (val > maxVal) {
                maxVal = val;
            }
            weightVec.push_back( val );
        }
        weightMatrix.push_back(weightVec);
    }
}

IdxValPair AuctionOracleRestricted::getOptimalBid(const IdxType bidderIdx) 
{
    assert(bidderIdx >=0 and bidderIdx < static_cast<IdxType>(bidders.size()) );

    const auto bidder = bidders[bidderIdx];
    
    IdxType bestItemIdx { -1 };
    double bestItemValue { std::numeric_limits<double>::max() };
    IdxType secondBestItemIdx { -1 };
    double secondBestItemValue { std::numeric_limits<double>::max() };

    // find best items idx
    for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        // non-diagonal point should be matched either to another
        // non-diagonal point or to its own projection
        if (isRestricted and bidder.isNormal() ) {
            auto item = items[itemIdx];
            if (item.isDiagonal() and item.projId != bidder.id)
                continue;
        }
        auto currItemValue = weightMatrix[bidderIdx][itemIdx] + prices[itemIdx];
        if ( currItemValue < bestItemValue ) {
            bestItemValue = currItemValue;
            bestItemIdx = itemIdx;
        }
    }

    // find second best items idx and value

    for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        // non-diagonal point should be matched either to another
        // non-diagonal point or to its own projection
        if (isRestricted and bidder.isNormal() ) {
            auto itemsItem = items[itemIdx];
            if (DiagramPoint::DIAG == itemsItem.type and itemsItem.projId != bidder.id)
                continue;
        }

        if (static_cast<IdxType>(itemIdx) == bestItemIdx)
            continue;

        auto currItemValue = weightMatrix[bidderIdx][itemIdx] + prices[itemIdx];
        if ( currItemValue < secondBestItemValue ) {
            secondBestItemValue = currItemValue;
            secondBestItemIdx = itemIdx;
        }
    }

    assert(bestItemValue <= secondBestItemValue);

    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemValue = " << bestItemValue << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << topIter->first << "; secondBestValue = " << secondBestItemValue << "; secondBestPrice = " << prices[topIter->first] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;
    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemValue = " << bestItemValue << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << secondBestItemIdx << "; secondBestValue = " << secondBestItemValue << "; secondBestPrice = " << prices[secondBestItemIdx] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;
    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemsDist= " << (weightAdjConst -  bestItemValue) << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << topIter->first << "; secondBestDist= " << (weightAdjConst - secondBestItemValue) << "; secondBestPrice = " << prices[topIter->first] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;

    // bid value: price + value difference + epsilon
    
    return std::make_pair(bestItemIdx, 
                          prices[bestItemIdx] + 
                          ( -bestItemValue + secondBestItemValue ) +
                          epsilon );
}

void AuctionOracleRestricted::setPrice(const IdxType itemIdx, const double newPrice)
{
    assert(prices.at(itemIdx) < newPrice );
    prices[itemIdx] = newPrice;
}

// *****************************
// AuctionOracleKDTreeRestricted
// *****************************

AuctionOracleKDTreeRestricted::AuctionOracleKDTreeRestricted(const std::vector<DiagramPoint>& _bidders, 
        const std::vector<DiagramPoint>& _items, 
        double _wassersteinPower) :
    AuctionOracleAbstract(_bidders, _items, _wassersteinPower),
    heapHandlesIndices(items.size(), std::numeric_limits<size_t>::max()),
    kdtreeItems(items.size(), std::numeric_limits<size_t>::max()),
    biddersToProjItems(bidders.size(), std::numeric_limits<size_t>::max()),
    bestDiagonalItemsComputed(false)
{
    size_t dnnItemIdx { 0 };
    size_t trueIdx { 0 };
    dnnPoints.clear();
    // store normal items in kd-tree
    for(const auto& g : items) {
        if (g.isNormal() ) {
            kdtreeItems[trueIdx] = dnnItemIdx;
            // index of items is id of dnn-point
            DnnPoint p(trueIdx);
            p[0] = g.x;
            p[1] = g.y;
            dnnPoints.push_back(p);
            assert(dnnItemIdx == dnnPoints.size() - 1);
            dnnItemIdx++;
        }
        trueIdx++;
    }

    assert(dnnPoints.size() < items.size() );
    for(size_t i = 0; i < dnnPoints.size(); ++i) {
        dnnPointHandles.push_back(&dnnPoints[i]);
    }
    DnnTraits traits;
    //std::cout << "kdtree: " << dnnPointHandles.size() << " points" << std::endl;
    kdtree = new dnn::KDTree<DnnTraits>(traits, dnnPointHandles, wassersteinPower);
    
    size_t handleIdx {0};
    for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        if (items[itemIdx].type == DiagramPoint::DIAG) {
            heapHandlesIndices[itemIdx] = handleIdx++;
            diagHeapHandles.push_back(diagItemsHeap.push(std::make_pair(itemIdx, 0)));
        }
    }
    //to-do: remove maxVal from 
    //std::cout << "3getFurthestDistance3Approx = " << getFurthestDistance3Approx(_bidders, _items) << std::endl;
    maxVal = 3*getFurthestDistance3Approx(_bidders, _items);
    maxVal = pow(maxVal, wassersteinPower);
    weightAdjConst = maxVal;
    //std::cout << "AuctionOracleKDTreeRestricted: weightAdjConst = " << weightAdjConst << std::endl;
    //std::cout << "AuctionOracleKDTreeRestricted constructor done" << std::endl;
    // todo: this must be done in readFiles procedure
    //std::cout << "started wasting time..." << std::endl;
    for(size_t bidderIdx = 0; bidderIdx < bidders.size(); ++bidderIdx) {
        for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
            auto bidder = bidders[bidderIdx];
            auto item = items[itemIdx];
            if (bidder.id == item.projId) {
                biddersToProjItems[bidderIdx] = itemIdx;
            }
        }
    }
    //std::cout << "stopped wasting time." << std::endl;
}

DebugOptimalBid AuctionOracleKDTreeRestricted::getOptimalBidDebug(IdxType bidderIdx)
{
    DebugOptimalBid result;
    DiagramPoint bidder = bidders[bidderIdx];

    //std::cout << "bidder.x = " << bidderDnn[0] << std::endl;
    //std::cout << "bidder.y = " << bidderDnn[1] << std::endl;

    // corresponding point is always considered as a candidate
    // if bidder is a diagonal point, projItem is a normal point, 
    // and vice versa.

    size_t projItemIdx = biddersToProjItems[bidderIdx];
    assert( 0 <= projItemIdx and projItemIdx < items.size() );
    DiagramPoint projItem = items[projItemIdx];
    assert(projItem.type != bidder.type);
    assert(projItem.projId == bidder.id);
    assert(projItem.id == bidder.projId);
    // todo: store precomputed distance?
    double projItemValue = pow(distLInf(bidder, projItem), wassersteinPower) + prices[projItemIdx];
   
    if (bidder.isDiagonal()) {
        // for diagonal bidder the only normal point has already been added
        // the other 2 candidates are diagonal items only, get from the heap
        // with prices
        assert(diagItemsHeap.size() > 1);
        if (!bestDiagonalItemsComputed) {
            auto topDiagIter = diagItemsHeap.ordered_begin();
            bestDiagonalItemIdx = topDiagIter->first;
            bestDiagonalItemValue = topDiagIter->second;
            topDiagIter++;
            secondBestDiagonalItemIdx = topDiagIter->first;
            secondBestDiagonalItemValue = topDiagIter->second;
            bestDiagonalItemsComputed = true;
        }

        if ( projItemValue < bestDiagonalItemValue) {
            result.bestItemIdx = projItemIdx;
            result.bestItemValue = projItemValue;
            result.secondBestItemIdx = bestDiagonalItemIdx;
            result.secondBestItemValue = bestDiagonalItemValue;
        } else if (projItemValue < secondBestDiagonalItemValue) {
            result.bestItemIdx = bestDiagonalItemIdx;
            result.bestItemValue = bestDiagonalItemValue;
            result.secondBestItemIdx = projItemIdx;
            result.secondBestItemValue = projItemValue;
        } else {
            result.bestItemIdx = bestDiagonalItemIdx;
            result.bestItemValue = bestDiagonalItemValue;
            result.secondBestItemIdx = secondBestDiagonalItemIdx;
            result.secondBestItemValue = secondBestDiagonalItemValue;
        }
    } else {
        // for normal bidder get 2 best items among non-diagonal points from
        // kdtree
        DnnPoint bidderDnn;
        bidderDnn[0] = bidder.getRealX();
        bidderDnn[1] = bidder.getRealY();
        auto twoBestItems = kdtree->findK(bidderDnn, 2);
        //std::cout << "twoBestItems for all: " << twoBestItems[0].d << " " << twoBestItems[1].d << std::endl;
        size_t bestNormalItemIdx { twoBestItems[0].p->id() };
        double bestNormalItemValue { twoBestItems[0].d };
        size_t secondBestNormalItemIdx { twoBestItems[1].p->id() };
        double secondBestNormalItemValue { twoBestItems[1].d };

        if ( projItemValue < bestNormalItemValue) {
            result.bestItemIdx = projItemIdx;
            result.bestItemValue = projItemValue;
            result.secondBestItemIdx = bestNormalItemIdx;
            result.secondBestItemValue = bestNormalItemValue;
        } else if (projItemValue < secondBestNormalItemValue) {
            result.bestItemIdx = bestNormalItemIdx;
            result.bestItemValue = bestNormalItemValue;
            result.secondBestItemIdx = projItemIdx;
            result.secondBestItemValue = projItemValue;
        } else {
            result.bestItemIdx = bestNormalItemIdx;
            result.bestItemValue = bestNormalItemValue;
            result.secondBestItemIdx = secondBestNormalItemIdx;
            result.secondBestItemValue = secondBestNormalItemValue;
        }
    }

    return result;

    //std::cout << "got result: " << result << std::endl;
    //double bestItemsPrice = prices[bestItemIdx];
    //if (items[result.bestItemIdx].type == DiagramPoint::DIAG) {
        //double bestItemValue1 = pow( distLInf(bidder, items[result.bestItemIdx]), wassersteinPower) + prices[result.bestItemIdx];
        //if ( fabs(result.bestItemValue - bestItemValue1) > 1e-6 ) {
            //std::cerr << "XXX: " << result.bestItemValue << " vs " << bestItemValue1 << std::endl;
            //result.bestItemValue = bestItemValue1;
        //}

    //}


    // checking code
    
    /*
    DebugOptimalBid debugMyResult(result);
    DebugOptimalBid debugNaiveResult;
    debugNaiveResult.bestItemValue = 1e20;
    debugNaiveResult.secondBestItemValue = 1e20;
    double currItemValue;
    for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        //if ( bidders[bidderIdx].type == DiagramPoint::NORMAL and
                //items[itemIdx].type == DiagramPoint::DIAG and
                //bidders[bidderIdx].projId != items[itemIdx].id)
            //continue;

        currItemValue = pow(distLInf(bidders[bidderIdx], items[itemIdx]), wassersteinPower) + prices[itemIdx];
        if (currItemValue < debugNaiveResult.bestItemValue) {
            debugNaiveResult.bestItemValue = currItemValue;
            debugNaiveResult.bestItemIdx  = itemIdx;
        }
    }

    for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
        if (itemIdx == debugNaiveResult.bestItemIdx) {
            continue;
        }
        currItemValue = pow(distLInf(bidders[bidderIdx], items[itemIdx]), wassersteinPower) + prices[itemIdx];
        if (currItemValue < debugNaiveResult.secondBestItemValue) {
            debugNaiveResult.secondBestItemValue = currItemValue;
            debugNaiveResult.secondBestItemIdx = itemIdx;
        }
    }
    //std::cout << "got naive result" << std::endl;

    if ( fabs( debugMyResult.bestItemValue - debugNaiveResult.bestItemValue ) > 1e-6 or
            fabs( debugNaiveResult.secondBestItemValue - debugMyResult.secondBestItemValue) > 1e-6 ) {
        kdtreeAll->printWeights();
        std::cerr << "bidderIdx = " << bidderIdx << "; ";
        std::cerr << bidders[bidderIdx] << std::endl;
        for(size_t itemIdx = 0; itemIdx < items.size(); ++itemIdx) {
            std::cout << itemIdx << ": " << items[itemIdx] << "; price = " << prices[itemIdx] << std::endl;
        }
        std::cerr << "debugMyResult: " << debugMyResult << std::endl;
        std::cerr << "debugNaiveResult: " << debugNaiveResult << std::endl;
        //std::cerr << "twoBestItems: " << twoBestItems[0].d << " " << twoBestItems[1].d << std::endl;
        assert(false);
    }
    //std::cout << "returning" << std::endl;

    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemValue = " << bestItemValue << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << secondBestItemIdx << "; secondBestValue = " << secondBestItemValue << "; secondBestPrice = " << prices[secondBestItemIdx] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;
    //std::cout << "getOptimalBid: bidderIdx = " << bidderIdx << "; bestItemIdx = " << bestItemIdx << "; bestItemsDist= " << (weightAdjConst - bestItemValue) << "; bestItemsPrice = " << prices[bestItemIdx] << "; secondBestItemIdx = " << secondBestItemIdx << "; secondBestDist= " << (weightAdjConst - secondBestItemValue) << "; secondBestPrice = " << prices[secondBestItemIdx] <<  "; bid = " << prices[bestItemIdx] + ( bestItemValue - secondBestItemValue ) + epsilon << "; epsilon = " << epsilon << std::endl;
    */
    return result;
}

IdxValPair AuctionOracleKDTreeRestricted::getOptimalBid(IdxType bidderIdx)
{

    
    DiagramPoint bidder = bidders[bidderIdx];

    //std::cout << "bidder.x = " << bidderDnn[0] << std::endl;
    //std::cout << "bidder.y = " << bidderDnn[1] << std::endl;

    // corresponding point is always considered as a candidate
    // if bidder is a diagonal point, projItem is a normal point, 
    // and vice versa.
    
    size_t bestItemIdx;
    double bestItemValue;
    double secondBestItemValue;


    size_t projItemIdx = biddersToProjItems[bidderIdx];
    assert( 0 <= projItemIdx and projItemIdx < items.size() );
    DiagramPoint projItem = items[projItemIdx];
    assert(projItem.type != bidder.type);
    assert(projItem.projId == bidder.id);
    assert(projItem.id == bidder.projId);
    // todo: store precomputed distance?
    double projItemValue = pow(distLInf(bidder, projItem), wassersteinPower) + prices[projItemIdx];
   
    if (bidder.isDiagonal()) {
        // for diagonal bidder the only normal point has already been added
        // the other 2 candidates are diagonal items only, get from the heap
        // with prices
        assert(diagItemsHeap.size() > 1);
        if (!bestDiagonalItemsComputed) {
            auto topDiagIter = diagItemsHeap.ordered_begin();
            bestDiagonalItemIdx = topDiagIter->first;
            bestDiagonalItemValue = topDiagIter->second;
            topDiagIter++;
            secondBestDiagonalItemIdx = topDiagIter->first;
            secondBestDiagonalItemValue = topDiagIter->second;
            bestDiagonalItemsComputed = true;
        }

        if ( projItemValue < bestDiagonalItemValue) {
            bestItemIdx = projItemIdx;
            bestItemValue = projItemValue;
            secondBestItemValue = bestDiagonalItemValue;
        } else if (projItemValue < secondBestDiagonalItemValue) {
            bestItemIdx = bestDiagonalItemIdx;
            bestItemValue = bestDiagonalItemValue;
            secondBestItemValue = projItemValue;
        } else {
            bestItemIdx = bestDiagonalItemIdx;
            bestItemValue = bestDiagonalItemValue;
            secondBestItemValue = secondBestDiagonalItemValue;
        }
    } else {
        // for normal bidder get 2 best items among non-diagonal points from
        // kdtree
        DnnPoint bidderDnn;
        bidderDnn[0] = bidder.getRealX();
        bidderDnn[1] = bidder.getRealY();
        auto twoBestItems = kdtree->findK(bidderDnn, 2);
        //std::cout << "twoBestItems for all: " << twoBestItems[0].d << " " << twoBestItems[1].d << std::endl;
        size_t bestNormalItemIdx { twoBestItems[0].p->id() };
        double bestNormalItemValue { twoBestItems[0].d };
        double secondBestNormalItemValue { twoBestItems[1].d };

        if ( projItemValue < bestNormalItemValue) {
            bestItemIdx = projItemIdx;
            bestItemValue = projItemValue;
            secondBestItemValue = bestNormalItemValue;
        } else if (projItemValue < secondBestNormalItemValue) {
            bestItemIdx = bestNormalItemIdx;
            bestItemValue = bestNormalItemValue;
            secondBestItemValue = projItemValue;
        } else {
            bestItemIdx = bestNormalItemIdx;
            bestItemValue = bestNormalItemValue;
            secondBestItemValue = secondBestNormalItemValue;
        }
    }

    IdxValPair result;

    assert( secondBestItemValue >= bestItemValue );

    result.first = bestItemIdx;
    result.second = ( secondBestItemValue - bestItemValue ) + prices[bestItemIdx] + epsilon;
    return result;
}
/*
a_{ij} = d_{ij} 
value_{ij} = a_{ij} + price_j
*/
void AuctionOracleKDTreeRestricted::setPrice(IdxType itemIdx, double newPrice)
{
    assert(prices.size() == items.size());
    assert( 0 < diagHeapHandles.size() and diagHeapHandles.size() <= items.size());
    assert(newPrice > prices.at(itemIdx));
    prices[itemIdx] = newPrice;
    if ( items[itemIdx].isNormal() ) {
        //std::cout << "before increasing weight in kdtree " << std::endl;
        //std::cout << kdtreeItems.at(itemIdx) << std::endl;
        assert(0 <= itemIdx and itemIdx < kdtreeItems.size());
        assert(0 <= kdtreeItems[itemIdx] and kdtreeItems[itemIdx] < dnnPointHandles.size());
        kdtree->increase_weight( dnnPointHandles[kdtreeItems[itemIdx]], newPrice);
        //std::cout << "after increasing weight in kdtree" << std::endl;
    } else {
        //std::cout << "before decreasing weight in heap" << std::endl;
        //std::cout << "diagHeapHandles.size = " << diagHeapHandles.size() << std::endl;
        assert(diagHeapHandles.size() > heapHandlesIndices.at(itemIdx));
        diagItemsHeap.decrease(diagHeapHandles[heapHandlesIndices[itemIdx]], std::make_pair(itemIdx, newPrice));
        bestDiagonalItemsComputed = false;
    }
}

void AuctionOracleKDTreeRestricted::adjustPrices(void)
{
}

AuctionOracleKDTreeRestricted::~AuctionOracleKDTreeRestricted()
{
    delete kdtree;
}

void AuctionOracleKDTreeRestricted::setEpsilon(double newVal) 
{
    assert(newVal >= 0.0);
    epsilon = newVal;
}
