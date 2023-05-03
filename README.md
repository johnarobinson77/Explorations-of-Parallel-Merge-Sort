# Explorations-of-Parallel-Merge-Sort
C++ code and a writeup of various forms of merge sort and how well they parallelize.
# Explorations of Parallel Merge Sort


## Abstract

Merge sort, described by Sedgewick[1][2], is the most versatile and flexible sort algorithm.  Besides being a general sort algorithm, it is deterministic because its performance is O(nlogn) for all data input and is straightforward to parallelize.

I want to explore various ways to parallelize merge sort in this project.  The standard recursive binary subdivision algorithm has some limitations in this regard, and I wanted to see if I could overcome those limitations.  I do not claim that the ideas presented here and in the code are original, only that I have yet to see some implemented.  Those that I have are referenced below.

The first limitation of the recursive or depth-first algorithm is that when implementing parallel versions, the work distributed to different threads only balances well for powers-of-2 numbers of threads due to the binary subdivision of the data at each level of recursion.  The first exploration is to use a breadth-first implementation in which for loops replace recursion.  It is much easier to balance the work among any number of threads when parallelizing a for loop.

The second limitation is that the recursive merge functions only use at most two threads to improve performance at the upper levels where there are fewer regions to merge than threads.  Implemented here is a CPU version of the Greenand, McColl, and Bader [3] algorithm for parallel implementation of the merge function adapted from their GPU algorithm.


## Introduction, 

Four variations of the merge sort algorithm are implemented in the code.  They use combinations of depth-first vs. breadth-first and the two optimized versions of the merge function.



1. mSortRecFR() is the closest to the original algorithm described by Sedgewick.
2. mSortRecFF() is a variation that uses a different optimized merge function.
3. mSortFR() is a breadth-first algorithm in which for loops replace recursion
4. mSortFF() is a variation of the breadth-first algorithm that uses the alternate merge function.

Thrown in for good measure are



5. qSort() is a single-thread-based baseline algorithm with only a few optimizations.
6. cSort() is a wrapper around the std::sort() function.

Finally, there is a sortPerThread() function, which runs any of the above five sort algorithms in parallel on subsections of the input data and then merges the results.

The Performace section will discuss the performance of each of these.


## The Code

The above sorting functions are implemented in a class called bfMergeSort&lt;> in file bfMergeSort.hpp.


### The Merge Functions 

The essential low-level operation of all the merge sort algorithms described here is the merge operation in which two sorted regions of data in a src array are merged into a single region in a dst array.


#### mergeFR()

The mergeFR() function implements the one described by Sedgwick.  It assumes that the two regions in the src array are adjacent and that the upper region is written in reverse order.  The index of the lower region starts at the bottom and increments as the merge operation proceeds.  The index of the upper region starts, and the top of that region in decrements.  The merge loop continues until the output index reaches the size of the upper plus the lower region.  

The input indices will naturally stop incrementing or decrementing at the maximum value, either at the top of the lower region or the bottom of the upper region.  Therefore capping the two input indices to ensure they stay within their regions is unnecessary.  This saves two compare and branches in the inner loop of the algorithm.  Refer to the code in mergeFR() function for the implementation of this merge function.


#### mergeFF()

The downside to Sedgewick’s merge is that it requires some complexity when arranging for the upper region to be placed in memory in the reverse order.  An alternate merge function is implemented here that allows the data in the upper and lower regions to be arranged in the forward direction.  

Before starting the merge function, the last value of the upper and lower regions are compared.  If the lower region has the smaller value, a merge loop is started that ends when the lower region's index reaches that region's end.  After the merge loop is complete, a simple copy of the rest of the upper region is performed.  Similar loops are executed with the opposite indices if the upper region has a smaller value.  This merge function is also optimal because there is only one test of the indices in the inner loop.  It has the advantage that the merged regions do not have to be adjacent in memory.  This merge operation is implemented in the mergeFF() function.


#### parallelMerge()

The parallelMerge() function is a CPU implementation of the algorithm first described by Greenand McColl and Bader [3] for using multiple threads to perform a merge GPUs.  In this algorithm, there is a pre-step that, given equal segments of the output equal to the number of threads, determines the bounding indices of the two input regions required to merge into those segments.  Once these bounding indices are determined, each thread can independently merge that segment.

As mentioned, there is a pre-step required for the merge function.  This pre-step is done in roughly O(logn) time, where n, in this case, is the total number of elements to be merged.  When the number of elements to be merged is large, the pre-step time is inconsequential.  However, it should be avoided for smaller merges.

There is one other two-thread merge in which one thread starts at the bottom and merges halfway, while the other thread merges from the top of the output region and merges halfway down.  That code is not implemented here, but an implementation can be found in the sorting portion of a kdTree builder algorithm by Brown[4] and specifically in this file https://github.com/chezruss/kd-tree/blob/master/kdTreeKmapKnlogn.cpp. 


### Depth-first Algorithms

The algorithm described by Sedgewick can be described as a depth-first algorithm.  The basic algorithm is to subdivide the data into two halves and recurse with both halves.  Recursion continues until there is only one element.  Then on returning from both sides of the recursion, the subdivided regions are merged into a single region using one of the merge functions above.


#### mSortRecFR()

The mSortRecFF sort implementation is closest to Sedgewick's published algorithm.  It uses the optimized forward-reverse merger coded in mergeFR().  1) It merges the data from one array into another but does not copy it back.  Instead, it counts on the merge at the next level up to it back in the original; array.  2) To accomplish the reverse ordering of the data in the upper portion, the direction of the sort is flipped for each successive call to the recursive function on the upper half only.  3) Recursion ends at the level where there are ~8 to 16 elements and performs an insertion sort on that data.  The insertion sort avoids many recursive function calls and short inefficient merges.  Depending on the total depth, the insertion sort is either in-place or part of a copy to the swap array such that the last merge will end up in the original array.


#### mSortRecFF()

mSortRecFF uses the same recursive algorithm as the mSortFR() but uses the optimized forward-forward merge coded in mergeFF() and therefore does not require reversing the order of the upper data.


#### Depth-first Restrictions

Both of the depth-first algorithms have the same multi-threaded restrictions.  At the first several levels of recursion, one of the two calls to the next level is done on a separate thread.  If four threads are available, a separate thread will be started in the first two levels of recursion.  However, if six threads are available, the multithreading will still only occurs at the first two levels, starting with four threads.  Using all six would cause an imbalance in the threads applied to each half of the data, and the performance is limited to the side with fewer threads.


### Breadth-first Algorithms

The following two algorithms are breadth-first because they calculate all of the merges that need to be calculated at a given level, starting with the finest level first and progressing to the course levels.  These algorithms are implemented with two nested for loops instead of recursion.  The inner loop iterates through the merges that must be done at a given level.  The outer loop iterates through the levels.  Ignoring optimizations for the moment, the first iteration of the outer does N/2 "merges" 1 data element each.  The next iteration does N/4 regions of 2 elements each.  The last iteration does one merge of N/2 elements each.

The advantage of the breadth-first algorithm is that it is easy to parallelize the inner loop with any number of threads available and maintain nearly equal loading on each thread.  The number of usable threads is not restricted to a power of two like the depth-first algorithms.  All threads must synchronize before the next level of merging starts.

mSortFF()

This function implements a bread-first merge sort with the following optimizations.  The first depth iteration is a loop on insertion sorts of data regions between 8 and 16 elements.  This saves the outer loop's first 3 or 4 iterations and many short merges.   Depending on the total depth iteration required, the insertion sort is either in place or part of a copy to the swap array such that the last iteration merge will end up in the original array.

Note that there are not enough merges to be given to the available threads at the final levels of iterations.  So at the higher levels where the data to be merged is large, it is better to use the parallelMerge() function which does a better job distributing the work to the available threads.

mSortFR()

mSortFR is also a breadth-first algorithm but uses the forward-reverse merger in mergeFR.  Determining which data needs to be sorted in the forward direction and which in the reverse order is complicated and results in messy code.  After the development of mSortFF(), further work on this algorithm was not justified.


### Miscellaneous Sort Algorithms


#### qSort()

qSort() implements the popular quicksort algorithm.  It is not a variation of the merge sort algorithm.  It is included here only for demonstration purposes.

Two standard optimizations are included in this code.  The first is to use insertion sort at the finest level of sorting.  The second is to use a slightly intelligent pivot picker to avoid the O(n^2) performance issues with certain data sets.  There is no attempt to use multithreading. 


#### cSort()

cSort() is a wrapper function around std::sort().  It is provided for comparison purposes.  It is not, to my knowledge, parallelized.


### The Mixed Breadth-Depth-First Algorithm

While writing the code for these various parallelized merge sorts, it occurred to me that the following algorithm is possible:



1.  Break up the input data into equal size segments where the number of segments equals the number of threads.  
2. Sort those segments on individual threads using any sorting algorithm you want.  
3. Using the breadth-first technique, use parallelMerge() to merge those segments into a single sorted data set.

The advantage of this algorithm is that any number of threads/cores can be applied efficiently, even if the initial sorting is done with a depth-first sort.  As of this writing, the single-threaded depth-first merge sorts are faster than the single-threaded breadth-first algorithms.    

This algorithm is implemented in sortPerThread().  Note that the sortPerThread() parameter list includes a function pointer that can take any of the above sorting functions as an argument.


### Verification

The main() function in bfMergeSortTest.cpp sets up, runs, and verifies any sort algorithms described above.  It loops over the number of threads and the number of tests to run at each thread.

The data can either be a random size or a fixed size.  Also, the data can be random, ordered from smallest to largest or largest to smallest.  If random data is chosen, some amount of that data is copied through the array to ensure equal cases are tested.


## If verification is enabled, the data is checked to make sure the data is in the correct order.  A checksum is also checked to provide the minimum verification that the data was not corrupted.


## Performance

The code provided has been compiled under MSVS, g++-11, and clang on various machines.  However, the best performance was achieved using clang on all machines.


### Merge Sort Performance

The first 2 bar charts below show the performance of the four forms of merge sort by themselves.  The code is:

DF FR	- Depth-First or recursive algorithm using the forward-reverse mergeFR() function

DF FF	- Depth-First or recursive algorithm using the forward-forward mergeFF() function

BF FR	- Breadth-First or iterative algorithm using the forward-reverse mergeFR() function

BF FF	- Breadth-First or iterative algorithm using the forward-forward mergeFF() function

The first bar chart shows the performance on a Mac mini with the 3200 GHz core i7 Intel processor with six cores and 12 hyper-threads.

![/assets/macSpeedChart1.png](https://github.com/johnarobinson77/Explorations-of-Parallel-Merge-Sort/blob/main/macSpeedChart1.png)

The two things to note are that 1) the performance of the recursive algorithms only improves at power-of-two thread counts, which is expected.  2) The performance breadth-first algorithms saturate at four threads.  This is unexpected, and I have not been able to determine why.

Bar chart 2 is the performance data from a 48-core Graviton3 processor on AWS.  The cores run at 2600 GHz.  

![/assets/macSpeedChart1.png](https://github.com/johnarobinson77/Explorations-of-Parallel-Merge-Sort/blob/main/gvtSpeedChart1.png)

Note here that the performance of the breadth-first algorithms 1) improves with every increase in the core count, not just at power-of-two boundaries, and generally outperforming the depth-first algorithms.  2) The breadth-first algorithms' performance does not saturate the way they do in the i7.


### Parallel Top Sort Performance

The following two charts show six sorts run under the sortPerThread() function.  These two show how this method can parallelize any sort algorithm with any thread count.

![/assets/macSpeedChart1.png](https://github.com/johnarobinson77/Explorations-of-Parallel-Merge-Sort/blob/main/macSpeedChart2.png)

The above chart shows incremental performance with any thread count for all sort algorithms, including qSort, (quick sort), and cSort, (std::sort) on the Mac mini i7.

![/assets/macSpeedChart1.png](https://github.com/johnarobinson77/Explorations-of-Parallel-Merge-Sort/blob/main/gvtSpeedChart2.png)

Chart 4 shows the same as Chart 3 but on the Graviton3 processor.

I did try these sort algorithms up to 40 cores on the Graviton processor, but performance is pretty much saturated above 20.


## References

[1] R. Sedgewick. Mergesort. In Algorithms in C++,

     pages 165–166.  Addison-Wesley, Reading, MA, 1992.

[2] R. Sedgewick. Mergesort. In Algorithms in C++,

     pages 173–174.  Addison-Wesley, Reading, MA, 1992

[3] Greenand, McColl, and Bader.  GPU Merge Path: A GPU Merging Algorithm

     Proceedings of the 26th ACM International Conference on Supercomputing

[4] Brown R.  Building a Balanced k-d Tree in O(kn log n) Time

     https://arxiv.org/abs/1410.5420
