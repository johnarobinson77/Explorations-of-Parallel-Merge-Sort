// bfMergeSort.hpp

#pragma once
/**
 * Copyright (c) 2023 John Robinson.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

#include <cstring>
#include <stdint.h>
#include <algorithm>    // std::swap
#include <cmath>        // log2
#include "parallelFor.hpp"

#define sortFor false
#define sortRev true


template <typename T>
class bfMergeSort
{


 
 size_t threads;
 // typedef  T (*comp_t)(T& a, T& b);
 // comp_t comp;

 static inline
   int64_t comp(int64_t& a, int64_t& b) {
   return a - b;
 }

public:

  typedef  void (bfMergeSort<T>::* bfMergeSortMemFn)(T arry[], int64_t len);

  bfMergeSort() {
    setThreads(1);
  };

  bfMergeSort(size_t threads)
  {
    setThreads(threads);
  };

  void setThreads(size_t threads) {
    this->threads = threads;
    //    if (mpi != nullptr) delete mpi;
    //    mpi = new size_t(threads + 1);
  }

  ~bfMergeSort() {
    //    if (mpi != nullptr) delete mpi;
  }

private:
  // the next series of methods are utility fuctions that are used by the various sort routines.

#define min(a,b) ((a)<(b) ? (a) : (b))
#define max(a,b) ((a)>(b) ? (a) : (b))

  //this parity generator routine was was first attributed to https://graphics.stanford.edu/~seander/bithacks.html#ParityWith64Bits
  inline bool getParity32(size_t n)
  {
    n = n ^ (n >> 1);
    n = n ^ (n >> 2);
    n = (n & 0x11111111) * 0x11111111;
    return (n & 0x10000000) == 0x10000000;
  }

  // this routine does an interger devision rounded up.
  inline size_t iDivUp(size_t a, size_t b)
  {
    return ((a % b) == 0) ? (a / b) : (a / b + 1);
  }

  // ParallelMerge and it's supporting functions getMergPaths and mergePath are a CPU implementation of
  // parallel merge function developed for GPUs described in "GPU Merge Path: A GPU Merging Algorithm" by Greenand, McColl, and Bader.
  // Proceedings of the 26th ACM International Conference on Supercomputing

  // For a particular output element of a merge, find all the elements in valA and ValB that will be below it.  
  template <bool DIR>
  size_t mergePath(T valA[], int64_t aCount, T valB[], int64_t bCount, int64_t diag, int64_t N) {

    size_t begin = max(0, diag - bCount);
    size_t end = min(diag, aCount);

    while (begin < end) {
      size_t mid = begin + ((end - begin) >> 1);
      T  aVal = valA[mid];
      T  bVal = valB[diag - 1 - mid];
      bool pred = (DIR && (comp(aVal, bVal) > 0)) ||
        (!DIR && (comp(aVal, bVal) < 0));
      if (pred) begin = mid + 1;
      else end = mid;
    }
    return begin;
  }

   // Divide the output of the merge into #-of-threads equal segments and determine the elements of valA and valB below each segment boundary. 
  template <bool DIR>
  inline void getMergePaths(size_t mpi[], T valA[], int64_t aCount, T valB[], int64_t bCount, int64_t spacing, int64_t N) {

    int64_t partitions = threads;
    mpi[0] = 0;  
    mpi[partitions] = aCount;
    
    parallelFor(1, partitions, [=](int64_t i) {
      int64_t diag = i * spacing;
      mpi[i] = mergePath<DIR>(valA, aCount, valB, bCount, diag, N);
      }, threads);
  }

  // parallelMerge() merges two sorted ranges in the src array into a single sorted range in the dst array
  // aBeg and aEnd inclusive indicate one sorted range in the src array 
  // bBeg and bEnd inclusive indicate the other sorted range in the src array
  // dBeg indicates where in the dst array the merge list should start
  // After determining the ranges of the src segments that will go nto each output segment, use the mergeFF() function to 
  // merge those segments in parallel.
  template <bool DIR>
  inline void parallelMerge(T dst[], T src[], size_t aBeg, size_t aEnd, size_t bBeg, size_t bEnd, size_t dBeg, size_t N) {
    
    size_t mpi[1024];  // an array that hold the intermediate ranges during a parallel merge.

    int64_t aCount = aEnd - aBeg + 1;
    int64_t bCount = bEnd - bBeg + 1;
    int64_t spacing = iDivUp(aCount + bCount, threads);
    getMergePaths<DIR>(mpi, &(src[aBeg]), aCount, &(src[bBeg]), bCount, spacing, N);

    parallelFor(0, threads, [=](int64_t thread) {
      size_t grid = thread * spacing;		// Calculate the relevant index ranges into the source array
      size_t a0 = mpi[thread] + aBeg;				// for this partition
      size_t a1 = mpi[thread + 1] + aBeg;
      size_t b0 = (grid - mpi[thread]) + bBeg;
      size_t b1 = (min(aCount + bCount, grid + spacing) - mpi[thread + 1]) + bBeg;
      size_t wtid = dBeg + thread * spacing;  //Place where this thread will start writing the data
 
      if (a0 == a1) {							// If no a data just copy b
        for (size_t b = b0; b < b1; b++) dst[wtid++] = src[b];
      }
      else if (b0 == b1) {					// If no b data just copy a
        for (size_t a = a0; a < a1; a++) dst[wtid++] = src[a];
      }
      else { // else do a merge using the forward-forward merge
        mergeFF<DIR>(dst, src, a0, a1 - 1, b0, b1 - 1, wtid);
      }
    }, threads);
  }

  // this function is only for debug purposes.
  void print(std::string t, T * in, int n) {
    std::cout << t << " ";
    for (size_t i = 0; i < n; i++) std::cout << in[i] << " ";
    std::cout << std::endl;
  }

  // this function does a test and swap of two adjacent elements.
  template<bool DIR> // DIR=0 is forward order, DIR=1 is reverse order
  inline void testAndSwap(T src[], int64_t idx) {
    if (DIR) {
      if (comp(src[idx], src[idx + 1]) < 0) std::swap(src[idx], src[idx + 1]);  // reverse
    } 
    else {
      if (comp(src[idx], src[idx + 1]) > 0) std::swap(src[idx], src[idx + 1]);  // forward
    }
  }
   // same as the above function but the direction needs to be runtime selectable.
  inline void testAndSwap(bool dir, T src[], int64_t idx) {
    if (dir) {
      testAndSwap<sortRev>(src, idx);  // reverse
    }  else { 
      testAndSwap<sortFor>(src, idx);  // forward
    }
  }
  
  // this function does a swap on two adjacent elements when copying to another array.
  template<bool DIR>
  inline void testAndCopy(T dst[], T src[], int64_t idx) {
    if (DIR) {
      if (comp(src[idx], src[idx + 1]) < 0) { // reverse 
        dst[idx + 1] = src[idx]; dst[idx] = src[idx + 1];
      }
      else { 
        dst[idx] = src[idx]; dst[idx + 1] = src[idx + 1];
      }
    } else {
      if (comp(src[idx], src[idx + 1]) > 0) { // forward 
        dst[idx + 1] = src[idx]; dst[idx] = src[idx + 1];
      }
      else {
        dst[idx] = src[idx]; dst[idx + 1] = src[idx + 1];
      }
    }
  }
  // and this is the runtime selectable version.
  inline void testAndCopy(bool dir, T dst[], T src[], int64_t idx) {
    if (dir)
      testAndCopy<sortRev>(dst, src, idx);  // reverse
    else
      testAndCopy<sortFor>(dst, src, idx);  // forward
  }

  // This is the forward-reverse merge function.  It does a merge of two sorted segments.  
  // The twe segments are assumed to be adjacent in the src array but the segment
  // in the upper portion of the array is in reverse order.
  // The  source indices start at the bottom of the lower segment and top of the
  // upper segments.  The merge progresses by moving the smaller element of the lower or upper
  // segment and advancing it's index towards the middle.  
  template<bool DIR>
  inline void mergeFR(T dst[], T src[], int64_t bot, int64_t top) {
    if (DIR) {
      for (size_t ic = bot, jc = top, kc = bot; kc <= top; kc++) {  // reverse
        dst[kc] = comp(src[ic], src[jc]) > 0 ? src[ic++] : src[jc--];
      }
    } else {
      for (size_t ic = bot, jc = top, kc = bot; kc <= top; kc++) {  //forward
        dst[kc] = comp(src[ic], src[jc]) < 0 ? src[ic++] : src[jc--];
      }
    }
  }
  // Runtime selectable version of the direction.
  inline void mergeFR(T dst[], T src[], int64_t bot, int64_t top, bool dir) {
    if (dir)
      mergeFR<sortRev>(dst, src, bot, top);
    else 
      mergeFR<sortFor>(dst, src, bot, top);
  }

  // mergeFF merges two sorted ranges in the src array into a single sorted range in the dst array
  // aBeg and aEnd inclusive indicate one sorted range 
  // bBeg and bEnd inclusive indicate the other sorted range
  // dBeg indicates where in the dst array the merge list should start
  // First, the elements at aEnd and bEnd are compared.  The end index with the smaller value is used.
  // to cap the merge function.  After than, the reset of the other array is just copied to the dst array
  // The two ranges being merged do not have to be adjacent in memory.
  template<bool DIR>
  inline void mergeFF(T dst[], T src[], size_t aBeg, size_t aEnd, size_t bBeg, size_t bEnd, size_t dBeg) {
    if (DIR) {
      if (comp(src[aEnd], src[bEnd]) > 0) {  // determine which range will be completed first during a compare and copy loop
        while (aBeg <= aEnd) { // the a range will be completely copied first so only compare up the end of a
          dst[dBeg++] = comp(src[aBeg], src[bBeg]) > 0 ? src[aBeg++] : src[bBeg++];
        }
        while (bBeg <= bEnd) {  // then copy the rest of b
          dst[dBeg++] = src[bBeg++];
        }
      }
      else {
        while (bBeg <= bEnd) { // the b range will be completely copied first so only compare up the end of b
          dst[dBeg++] = comp(src[aBeg], src[bBeg]) >= 0 ? src[aBeg++] : src[bBeg++];
        }
        while (aBeg <= aEnd) {  // then copy the rest of a
          dst[dBeg++] = src[aBeg++];
        }
      }
    }
    else {
      if (comp(src[aEnd], src[bEnd]) < 0) { // determine which range will be completed first during a compare and copy loop
        while (aBeg <= aEnd) {  // the a range will be completely copied first so only compare up the end of a
          dst[dBeg++] = comp(src[aBeg], src[bBeg]) <= 0 ? src[aBeg++] : src[bBeg++];
        }
        while (bBeg <= bEnd) { // then copy the rest of b
          dst[dBeg++] = src[bBeg++];
        }
      }
      else {
        while (bBeg <= bEnd) {  // the b range will be completely copied first so only compare up the end of b
          dst[dBeg++] = comp(src[aBeg], src[bBeg]) < 0 ? src[aBeg++] : src[bBeg++];
        }
        while (aBeg <= aEnd) { // then copy the rest of a
          dst[dBeg++] = src[aBeg++];
        }
      }
    }
  }

  // The next two functions are simple insertion sorts used by most of the sort routines to speed up the finest level of sorting.
  // copyInsertionSort() does an insertion sort as part of a copy from the src to the dst array.
  template<bool DIR >
  inline void copyInsertionSort(T dst[], T src[], int64_t lower, int64_t upper) {
    dst[lower] = src[lower];
    for (size_t is = lower + 1; is <= upper; is++) { // for each item to be copied...
      int64_t jd;
      for (jd = lower; jd < is; jd++) {
        if (DIR ^ (comp(src[is], dst[jd]) < 0)) { // find the place where next copy should be placed
          for (int64_t k = is; k > jd; k--) { // move the ones further out of the way.
            dst[k] = dst[k - 1];
          }
          break;
        }
      }
      dst[jd] = src[is]; //and do the copy.
    }
  }

  // runtime selectable direction version of the above.
  inline void copyInsertionSort(T dst[], T src[],  int64_t lower, int64_t upper, bool flip) {
    if (flip)
      copyInsertionSort<sortRev>(dst, src, lower, upper);
    else
      copyInsertionSort<sortFor>(dst, src, lower, upper);
  }

  // inPlaceInsertionSort does an insertion sort of the data leaving at the same spot in the source array.
  template<bool DIR >
  inline void inPlaceInsertionSort(T src[], int64_t lower, int64_t upper) {
    for (size_t is = lower + 1; is <= upper; is++) {
      T tmp = src[is];
      int64_t js;
      if (DIR) for (js = lower; (comp(tmp, src[js]) < 0); js++); // find the place where temp should be placed
      else     for (js = lower; (comp(tmp, src[js]) > 0); js++); // find the place where temp should be placed
      for (int64_t k = is; k > js; k--) src[k] = src[k - 1]; // move the others out of the way
      src[js] = tmp; // place the current 
    }
  }

  // runtime selectable direction version of the above.
  inline void inPlaceInsertionSort(T src[], int64_t lower, int64_t upper, bool flip) {
    if (flip)
      inPlaceInsertionSort<sortRev>(src, lower, upper);
    else
      inPlaceInsertionSort<sortFor>(src, lower, upper);
  }

// pivotPicker() is used as part of the quick sort function qSort.  It makes an attempt to mick a the median 
// the data elements in the range upper to lower inclusive.  It is needed only to avoid the worst case 
// pivots that could lead to poor performance.
  int64_t pivotPicker(T arry[], int64_t lower, int64_t upper) {
    const size_t maxIdxs = 11;
    int64_t idxs = min(maxIdxs, (upper - lower));
    int64_t medianIdxs[maxIdxs];
    
    // pick idxs equally spaced indices in the range lower to upper and sort them into 
    // medianIdxs[] according the the value they point to.
    double delta = double(upper - lower) / double(idxs - 1);
    medianIdxs[0] = lower;
    for (size_t is = 1; is < idxs; is++) {
      int64_t l = lower + llround(is * delta);
      int64_t jd;
      for (jd = 0; jd < is; jd++) {
        if (comp(arry[l], arry[medianIdxs[jd]]) < 0) { // find the place where next idx should be placed
          for (int64_t k = is; k > jd; k--) { // move the ones further out of the way.
            medianIdxs[k] = medianIdxs[k - 1];
          }
          break;
        }
      }
      medianIdxs[jd] = l; //and do the copy.
    }
    // return the middle of the sorted indices
    auto rv = medianIdxs[idxs >> 1];
    return rv;
  }

private:
  // This constant is used to determine the number lowest levels of sorting to be covered by insertion sort.
  // for the mSortFF algorithm.
  const int64_t insertionSortDepth = 4;

  // mSortFF sorts the data in arry using a breadth first algorithm.  It uses the forward-forward merge function.
public:
  template<bool DIR>
  void mSortFF(T arry[], const int64_t len) {

    T* swap;  // pointer to array that that the data will be swapped to during a merge sort
    swap = new T[len];

    const int64_t depth = ceil(log2(len)); // calculate the depth or number of levels of merging that need to be done.
    const bool depthOdd = ((depth - insertionSortDepth) & 0x1) != 1; // and whether the depth - insertionSortDepth is odd

    // numIterations iterations need at the lowest level where the insertion sort is done.
    // The number of iterations needs tp be a power of 2.
    int64_t numIterations = 1 << (depth - insertionSortDepth);
    // The spacing is the fractional number of elements that will be sorted at each level, starting with the insertion sort size.
    double spacing = double(len) / double(numIterations);

    // first level of sort is done with insertions sorts to cut down on the number of lowest level iterations where merges are inefficient.
    // Do an inPlaceInsertionSort() of a copyInsertion sort depending on whether the remaining loops are even or odd.  The reason is to
    // make sure that after all merges between the input arry and swap are done, the final sorted data is in arry.
    if (depthOdd) {
      parallelFor(0L, numIterations, [&](uint64_t start) {
        int64_t lb = llround(start * spacing);
        int64_t le = llround((start + 1) * spacing) - 1;
        inPlaceInsertionSort<sortFor>(arry, lb, le);
        }, threads);
    }
    else {
      parallelFor(0L, numIterations, [&](size_t start) {
        int64_t lb = llround(start * spacing);
        int64_t le = llround((start + 1) * spacing) - 1;
        copyInsertionSort<sortFor>(swap, arry, lb, le);
        }, threads);
    }

    //print("arry", arry, 32);
    //print("swap", swap, 32);

    // merge the next depth-insertionSortDepth levels.  For each successive level, the spacing doubles and the number of iteration halves
    for (size_t d = depth - insertionSortDepth; d > 0; d--, spacing *= 2.0) {
      numIterations /= 2;
      // set up the from/to arrays according to which direction the swap needs to happen
      const bool dOdd = (d & 0x1) == 1;
      T* fromPtr;
      T* toPtr;
      if (dOdd) { fromPtr = swap;  toPtr = arry; }
      else { fromPtr = arry;  toPtr = swap; }

      // determine at which level the merges are distributed among the threads or each merge is divided up by the number of threads
      if (d > 4) {
        parallelFor(0L, numIterations, [&](size_t start) {
          int64_t i = start * 2;
          int64_t lb = llround(i * spacing);
          int64_t lm = llround((i + 1) * spacing);
          int64_t le = llround((i + 2) * spacing) - 1;
          mergeFF<sortFor>(toPtr, fromPtr, lb, lm - 1, lm, le, lb);
          }, threads);
      }
      else {
        for (size_t start = 0; start < numIterations; start++) {
          int64_t i = start * 2;
          int64_t lb = llround(i * spacing);
          int64_t lm = llround((i + 1) * spacing);
          int64_t le = llround((i + 2) * spacing) - 1;
          parallelMerge<sortFor>(toPtr, fromPtr, lb, lm - 1, lm, le, lb, len);
        }
      }
      //print("to", toPtr, 16);
      //print("swap", swap, 16);

    }
    delete[] swap;
  }


  // mSortFR sorts the data in arry using a breadth first algorithm.  It uses the forward-reverse merge function
  // which requires a lot of extra complexity to insure that on every merge, the upper portion is in the 
  // reverse order.  THis is why some sorts and merges are done in reverse order.
public:
  template<bool DIR>
  void mSortFR(T arry[], const int64_t len) {

    T* swap;  // pointer to array that that the data will be swapped to during a merge function
    swap = new T[len];

    const int64_t depth = ceil(log2(len)); // calculate the depth or number of merge levels to be done
    const bool depthOdd = (depth & 0x1) == 1; // and whether it's odd

    // first level of sort.  Swap the adjacent values to put them in the right order. for the next merge level
    int64_t lo = len % 8;  // lo represents the number that elements that are left over after the for loop is complete.
    if (depthOdd) {
      int64_t i;
      int64_t sel;
      parallelFor(0L, ((len - lo) >> 3), [=](uint64_t start) {
        int64_t i = start << 3;
        if (!getParity32(start << 2)) {
          testAndSwap <sortFor>(arry, i + 0);
          testAndSwap <sortRev>(arry, i + 2);
          testAndSwap <sortRev>(arry, i + 4);
          testAndSwap <sortFor>(arry, i + 6);
        } else {
          testAndSwap <sortRev>(arry, i + 0);
          testAndSwap <sortFor>(arry, i + 2);
          testAndSwap <sortFor>(arry, i + 4);
          testAndSwap <sortRev>(arry, i + 6);
        }
        }, threads);
      // handle some overflow if the number of elements is not evenly divisible by 8
      i = len - lo;
      sel = i / 2;
      for (; i <= len - 2; i += 2, sel++) {
        testAndSwap(getParity32(sel), arry, i);
      }
    } else {
      int64_t i;
      int64_t sel;
      parallelFor(0L, ((len - lo) >> 3), [=](uint64_t strt) {
        int64_t i = strt << 3;
        if (!getParity32(strt << 2)) {
          testAndCopy<sortFor>(swap, arry, i + 0);
          testAndCopy<sortRev>(swap, arry, i + 2);
          testAndCopy<sortRev>(swap, arry, i + 4);
          testAndCopy<sortFor>(swap, arry, i + 6);
        }
        else {
          testAndCopy<sortRev>(swap, arry, i + 0);
          testAndCopy<sortFor>(swap, arry, i + 2);
          testAndCopy<sortFor>(swap, arry, i + 4);
          testAndCopy<sortRev>(swap, arry, i + 6);
        }
        }, threads);

      // handle some overflow if the number of elements is not evenly divisible by 8
      i = len - lo;
      sel = i / 2;
      for (; i <= len - 2; i += 2, sel++) {
        testAndCopy(getParity32(sel), swap, arry, i);
      }
      if (i < len) swap[i] = arry[i];
    }

    //print("arry", arry, 32);
    //print("swap", swap, 32);

    // merge the next depth-1 levels
    for (int64_t d = depth - 1, logIncr = 2; d > 0; d--, logIncr += 1) {
      // set up the from/to arrays according to which direction the swap needs to happen
      const bool dOdd = (d & 0x1) == 1;
      T* fromPtr;
      T* toPtr;
      if (dOdd) {
        fromPtr = swap;  toPtr = arry;
      }
      else {
        fromPtr = arry;  toPtr = swap;
      }
      int64_t incr = 1LL << logIncr;
      lo = len % incr; // calculate the number of left overs.

      int64_t sel;
      auto* futures = parallelForNoWait(0L, len >> logIncr, [=](uint64_t strt) {
        int64_t i = strt << logIncr;
        mergeFR(toPtr, fromPtr,i, i + incr - 1, getParity32(strt));
        }, threads);
      if (lo != 0) { // check for left overs
        int64_t i = len - lo;
        int64_t sel = i / incr;
        mergeFR(toPtr, fromPtr,i, len - 1, getParity32(sel));
      }
      parallelForFinish(futures);
      //print("to", toPtr, 36);
      //print("swap", swap, 16);
    }
    delete swap;
    swap = nullptr;
  }

  // mSortRecFR is the closest to the textbook merge sort algorithm.  It uses the forward-reverse mergeFR function
  // for the merges.  mSortRecFR calculates the total recursive depth and multithreading depth.  Then calls the
  // recursive subdivideAndMergeFR function to do the actual sorting.  Note that it only can utilize a power-of-two 
  // number of threads
public:
  template<bool DIR>
  void mSortRecFR(T arry[], const int64_t len) {

    T* swap;  // pointer to array that that the data will be swapped to during a merge function
    swap = new T[len];

    int64_t depth = floor(log2(len)) - 1; // calculate the depth
    const bool depthOdd = (depth & 0x1) == 1; // and whether it's odd

    int64_t threadLevel = floor(log2(threads));

    if (!depthOdd) depth--;
    subdivideAndMergeFR<true>(arry, swap, 0, len - 1, depth, DIR, threadLevel);

    delete swap;
    swap = nullptr;
  }

private:
  // this is the recursive function that does the sort for mSortRecFR.  
  // Note, to work correctly, the initial value of depth must be odd so at the recursion termination point, arry is in the src position.    
  template<bool MTHREAD>
  void subdivideAndMergeFR(T dst[], T src[], int64_t lower, int64_t upper, int64_t depth, bool flip, int64_t threadLevel = -1) {


    // check to see if we are done recusing and do a copy insertion sort
    if (depth == 0) { // at the lowest level, do a 
      copyInsertionSort(dst, src, lower, upper, flip);
      // and return;
      return;
    }

    // calculate the midpoint for the next call
    int64_t mid = (lower + upper) >> 1;
    // subdivide the array at near the midpoint and descend
    // note that at the next level the src and dst flip and the direction of the sort flips only for the upper half.
    // also the src become the dst and the dst becomes the sort.
    if (threadLevel > 0) { // if the thread level has not been reached, start the lower subdivision on a separate thread. 
      std::future<void> f = std::async(std::launch::async, [&] {
        subdivideAndMergeFR<true>(src, dst, lower, mid - 1, depth - 1, flip, threadLevel - 1); // descend down the lower half
        });
      subdivideAndMergeFR<true>(src, dst, mid, upper, depth - 1, !flip, threadLevel - 1); // descend down the upper half half.  Reverse the sort order
      f.wait();
    }
    else {
      subdivideAndMergeFR<false>(src, dst, lower, mid - 1, depth - 1, flip, threadLevel - 1); // descend down the lower half
      subdivideAndMergeFR<false>(src, dst, mid, upper, depth - 1, !flip, threadLevel - 1); // descend down the upper half half  Reverse the sort order
    }
    // merge the results of the lower level recursion.
    mergeFR(dst, src, lower, upper, flip);
    // and return;
    return;

  }

  // mSortRecFR is similar to the textbook merge sort algorithm.  It uses the forward-forward mergeFR function
  // for the merges.  mSortRecFF calculates the total recursive depth and multithreading depth.  Then calls the
  // recursive subdivideAndMergeFF function to do the actual sorting.  Note that it only can utilize a power-of-two 
  // number of threads.
public:
  template<bool DIR>
  void mSortRecFF(T arry[], const int64_t len) {
    T* swap;  // pointer to array that that the data will be swapped to during a merge function
    swap = new T[len];

    int64_t depth = floor(log2(len)) - 1; // calculate the depth
    const bool depthOdd = (depth & 0x1) == 1; // and whether it's odd

    int64_t threadLevel = floor(log2(threads));

    if (!depthOdd) depth--;
    subDivideAndMergeFF<DIR>(arry, swap, 0, len - 1, depth, threadLevel);

    delete swap;
    swap = nullptr;
  }

  // this is the recursive function that does the sort for mSortRecFF.  
  // Note, to work correctly, the initial value of depth must be odd so at the recursion termination point, arry is in the src position.    
private:
  template<bool DIR>
  void subDivideAndMergeFF(T dst[], T src[], int64_t lower, int64_t upper, int64_t depth, int64_t threadLevel = -1) {

    // check to see if we are done recusing and do a copy insertion sort and return;
    if (depth == 0) {
      copyInsertionSort<DIR>(dst, src, lower, upper);
      return;
    }

    // calculate the midpoint for the next call
    int64_t mid = (lower + upper) >> 1;
    // subdivide the array at near the midpoint and descend
    // note that at the next level the src and dst flip but the direction of the sort does not flip for the upper half.
    if (threadLevel > 0) {
      std::future<void> f = std::async(std::launch::async, [&] {
        subDivideAndMergeFF<DIR>(src, dst, lower, mid - 1, depth - 1, threadLevel - 1);
        });
      subDivideAndMergeFF<DIR>(src, dst, mid, upper, depth - 1, threadLevel - 1); // descend down the upper half half
      f.wait();
    }
    else {
      lowLevelMsortFF<DIR, true>(src, dst, lower, mid - 1, depth - 1); // descend down the lower half
      lowLevelMsortFF<DIR, true>(src, dst, mid, upper, depth - 1); // descend down the upper half half
    }

    // merge the results of the lower level recursion and return;
    mergeFF<DIR>(dst, src,lower, mid-1, mid, upper, lower);
    return;

  }

  // this is the recursive function that does the non-threaded sort for mSortRecFF.      
private:
  template<bool DIR, bool doCopyInsertionSort>
  void lowLevelMsortFF(T dst[], T src[], int64_t lower, int64_t upper, int64_t depth) {
    // check to see if we are done recusing and do an insertion sort
    if (depth == 0) {
      if (doCopyInsertionSort)
        copyInsertionSort<DIR>(dst, src, lower, upper);
      else
        inPlaceInsertionSort<DIR>(dst, lower, upper);
      return;
    }

    // calculate the midpoint for the next call
    int64_t mid = (lower + upper) >> 1;
    lowLevelMsortFF<DIR, doCopyInsertionSort>(src, dst, lower, mid - 1, depth - 1); // descend down the lower half
    lowLevelMsortFF<DIR, doCopyInsertionSort>(src, dst, mid, upper, depth - 1); // descend down the upper half half
    // merge the results of the lower level recursion and return.
    mergeFF<DIR>(dst, src, lower, mid - 1, mid, upper, lower);
    return;   // and return;
  }

  // this is a public wrapper around the quick sort function below.
public:
  const int64_t insertionSortSize = 21;  // Size of sort region below which insertion sort is used 
  template<bool DIR>
  void qSort(T arry[], const int64_t len) {

    subdivideAndPartition<sortFor>( arry, 0, len - 1);

  }

  // this is a standard non-treaded implementation of quick sort. It only sorts in the forward direction 
private:
  template<bool DIR>
  void subdivideAndPartition(T arry[], int64_t lower, int64_t upper) {

    // only use the elaborate picker above a certain segment size to cut down on it's
    // overhead to performance.  Otherwise just use the one at the top of the segment.
    if (upper - lower > 41) {
      auto idx = pivotPicker(arry, lower, upper);
      std::swap(arry[upper], arry[idx]);
    }
    T& pivot = arry[upper];
    int64_t il = lower;
    int64_t iu = upper - 1;
    while (true) {  // part ion the data about the pivot
      while (comp(pivot, arry[il]) > 0) il++;  // search bottom up for > pivot
      while (comp(pivot, arry[iu]) < 0 && (iu > il)) iu--;  // seacrh top down for < pivot,
      if (il < iu) {
        std::swap(arry[il++], arry[iu]);  // swap those values.
        //if (il > iu) iu++;
      }
      else break; // but if  il >= iu exit the loop
    }
    if (il < upper) std::swap(arry[iu], pivot); // If il has not gotten all the way to the pivot, swap the pivot into it's proper place.  
    else iu++;  //If it has, that indicates everything is <= pivot.  Increment iu so pivot is included in the next level. 
    if (iu - lower > insertionSortSize) subdivideAndPartition<sortFor>(arry, lower, iu - 1);
    else if (iu - 1 > lower) inPlaceInsertionSort<sortFor>(arry, lower, iu - 1);
    if (upper - iu > insertionSortSize) subdivideAndPartition<sortFor>(arry, iu + 1, upper);
    else if (upper > iu + 1) inPlaceInsertionSort<sortFor>(arry, iu + 1, upper);
    return;
  }

// cSort is just a wrapper around the c++ std::sort, used for comparison purposes.  it only sorts in the forward direction.
public:
  template<bool DIR>
  void cSort(T arry[], int64_t len) {
    std::sort(arry, arry + len);
  }

  // sortPerThread() is a function that can parallelize any of the 6 public sort functions above.  It does so by
  // dividing the input data into segments.  The number of segmens is equal to the number of threads.  Each segment
  // is sorted independently but at the same time by the function provided in singleSort.
  // Then it merges those segments using breadth first merging and the parallelMerge function iteratively until all
  // segments have been merged back into arry.
public:
  void sortPerThread(T arry[], const int64_t len, bfMergeSortMemFn singleSort) {
    size_t lThreads = threads; // save a local copy of the class threads and 
    threads = 1; //temporarily set class threads to 1.

    // calculate the fractional size of each segment to sort.
    // calculating the arry segments using doubles results in segment sizes where the max segment size
    // is only one bigger than the min.
    double delta = double(len) / double(lThreads);

    //sort lThread segments of the input arry using the sort method provided in the function pointer
    parallelFor(0, lThreads, [this, arry, delta, singleSort](int64_t i) {
      int64_t lb = llround(i * delta);
      int64_t le = llround((i + 1) * delta);
      (this->*singleSort)(arry + lb, le - lb);
      }, lThreads);

    //restore the class threads
    threads = lThreads;

    if (threads <= 1) return;

    // create a local swap space
    T* swap;  // pointer to array that that the data will be swapped to during a merge function
    swap = new T[len];

    const int64_t depth = ceil(log2(threads)) - 1; // calculate the number of depth iterations
    bool depthOdd = (depth & 0x1) == 1; // and whether it's odd
    // if it is not odd, then we need to copy arry to swap
    if (!depthOdd) parallelFor(0, threads, [arry, swap, delta](int64_t i) {
      int64_t lb = llround(i * delta);
      int64_t le = llround((i + 1) * delta);
      memcpy(swap + lb, arry + lb, (le - lb) * sizeof(T));
      }, threads);

    double next;
    double start;
    for (int64_t d = depth; d > 0; d--) {
      // set up the from/to arrays according to which direction the swap needs to happen
      const bool Odd = (d & 0x1) == 1;
      T* fromPtr, * toPtr;
      if (Odd) { fromPtr = arry;  toPtr = swap; }
      else { fromPtr = swap;  toPtr = arry; }
      start = 0.0;
      // merged the full sized pairs of of segments for this level.
      for (next = 2.0 * delta; llround(next) < len; next += 2.0 * delta) {
        parallelMerge<sortFor>(toPtr, fromPtr, llround(start), llround(start + delta) - 1, llround(start + delta), llround(next) - 1, llround(start), len);
        start = next;
      }
      // if there is an odd number segments to be merged, take care of the left overs.
      if (llround(next) >= len) {
        if (llround(start + delta) < len) {  // do a merge if there is a segment plus a partial segment
          parallelMerge<sortFor>(toPtr, fromPtr, llround(start), llround(start + delta) - 1, llround(start + delta), len - 1, llround(start), len);
        } else { // use just copy the partial segment.
          double incr = double(len - start) / double(lThreads);
          parallelFor(0, threads, [fromPtr, toPtr, start, incr](int64_t i) {
            int64_t lb = llround(start + i * incr);
            int64_t le = llround(start + (i + 1) * incr);
            memcpy(toPtr+lb, fromPtr+lb, (le - lb) * sizeof(T));
            }, threads);
        }
      }
      delta *= 2.0;
    }
    if (delta < len)  //THe last level my require an additional merge
      parallelMerge<sortFor>(arry, swap, 0, llround(delta) - 1, llround(delta), len - 1, 0, len);

    // clean up
    delete[] swap;

  }

};


