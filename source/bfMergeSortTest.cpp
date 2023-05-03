// bfMergeSortTest.cpp 
 
/**
* Copyright(c) 2023 John Robinson.
*
*SPDX - License - Identifier: BSD - 3 - Clause
*/

//

#include <iostream>
#include "bfMergeSort.hpp"
#include <random>
#include <chrono>
#include "parallelFor.hpp"
#include <vector>
#include <cstring>


// a slight rewrite of the Romdomer class from
// https://stackoverflow.com/questions/13445688/how-to-generate-a-random-number-in-c/53887645#53887645
// Usage example
// auto riTestData = RandomInterval<int64_t>(-10000000000LL, 10000000000LL, 1);
// int64_t ri = riTestData();

template <typename RDT>
class RandomInterval {
  // random seed by default
  std::mt19937 gen_;
  std::uniform_int_distribution<RDT> dist_;

public:
  RandomInterval(RDT min, RDT max, unsigned int seed = std::random_device{}())
    : gen_{ seed }, dist_{ min, max } {}

  // if you want predictable numbers
  void SetSeed(unsigned int seed) { gen_.seed(seed); }

  RDT operator()() { return dist_(gen_); }
};

// names for the available test data types
const int64_t dtRandom = 0;
const int64_t dtOrdered = 1;
const int64_t dtReverseOrdered = 2;

int main(int argc, char *argv[]) {

  // default settings
  bool fixedTestSize = true;
  size_t fixedTestSizeNum = 1024 * 1024 * 16;
  size_t minThreads = 1;
  size_t maxThreads = 4;
  bool verifySort = true;
  size_t algorithmSel = 1;
  size_t num_tests = 5;
  size_t dataType = dtRandom;

  bool argError = false;
  for (size_t arg = 1; arg < argc; arg++) {
    bool oneMore = arg < (argc - 1);
    if (strcmp(argv[arg], "-minT") == 0) {
      arg++;
      if (!oneMore || 0 == (minThreads = atoi(argv[arg]))) {
        std::cout << "-minT requires a non-zero integer argument." << std::endl;
        argError = true;
      }
    }
    else if (strcmp(argv[arg], "-maxT") == 0) {
      arg++;
      if (!oneMore || 0 == (maxThreads = atoi(argv[arg]))) {
        std::cout << "-maxT requires a non-zero integer argument." << std::endl;
        argError = true;
      }
    }
    else if (strcmp(argv[arg], "-n") == 0) {
      arg++;
      if (!oneMore || 0 == (fixedTestSizeNum = atoi(argv[arg]))) {
        std::cout << "-t requires a non-zero integer argument." << std::endl;
        argError = true;
      }
    }
    else if (strcmp(argv[arg], "-a") == 0) {
      arg++;
      if (!oneMore || 0 == (algorithmSel = atoi(argv[arg]))) {
        std::cout << "-a requires a non-zero integer argument." << std::endl;
        argError = true;
      }
    }
    else if (strcmp(argv[arg], "-l") == 0) {
      arg++;
      if (!oneMore || 0 == (num_tests = atoi(argv[arg]))) {
        std::cout << "-l requires a non-zero integer argument." << std::endl;
        argError = true;
      }
    }
    else if (strcmp(argv[arg], "-rs") == 0) {
      fixedTestSize = false;
    }
    else if (strcmp(argv[arg], "-dr") == 0) {
      dataType = dtRandom;
    }
    else if (strcmp(argv[arg], "-do") == 0) {
      dataType = dtOrdered;
    }
    else if (strcmp(argv[arg], "-db") == 0) {
      dataType = dtReverseOrdered;
    }
    else if (strcmp(argv[arg], "-nv") == 0) {
      verifySort = false;
    }
    else if (strcmp(argv[arg], "-v") == 0) {
      verifySort = true;
    }
    else {
      std::cout << "Argument " << argv[arg] << " not recognized" << std::endl;
    }
  }
  if (argError) {
    return 1;
  }

  switch (dataType) {
  case dtRandom: {
    std::cout << "Random data: "; break;
  }
  case dtOrdered: {
    std::cout << "Ordered data: "; break;
  }
  case dtReverseOrdered: {
    std::cout << "ReverseOrdered data: "; break;
  }
  default: {
    std::cout << "No such data type: " << dataType << std::endl;
    exit(1);
  }
  }


  switch (algorithmSel) {
  case 1: {
    std::cout << "Sort algorithm " << algorithmSel << ", merge sort BF FF" << std::endl; break;
  }
  case 2: {
    std::cout << "Sort algorithm " << algorithmSel << ", merge sort BF FR" << std::endl; break;
  }
  case 3: {
    std::cout << "Sort algorithm " << algorithmSel << ", merge sort DF FR" << std::endl; break;
  }
  case 4: {
    std::cout << "Sort algorithm " << algorithmSel << ", merge sort DF FF" << std::endl; break;
  }
  case 5: {
    std::cout << "Sort algorithm " << algorithmSel << ", quick sort" << std::endl; break;
  }
  case 6: {
    std::cout << "Sort algorithm " << algorithmSel << ", c std sort" << std::endl; break;
  }
  case 11: {
    std::cout << "Sort algorithm " << algorithmSel << ", merge sort BF FF Parallel Top" << std::endl; break;
  }
  case 12: {
    std::cout << "Sort algorithm " << algorithmSel << ", merge sort BF FR Parallel Top" << std::endl; break;
  }
  case 13: {
    std::cout << "Sort algorithm " << algorithmSel << ", merge sort DF FR Parallel Top" << std::endl; break;
  }
  case 14: {
    std::cout << "Sort algorithm " << algorithmSel << ", merge sort DF FF Parallel Top" << std::endl; break;
  }
  case 15: {
    std::cout << "Sort algorithm " << algorithmSel << ", quick sort Parallel Top" << std::endl; break;
  }
  case 16: {
    std::cout << "Sort algorithm " << algorithmSel << ", c std sort Parallel Top" << std::endl; break;
  }
  default: {
    std::cout << "No such algorithm: " << algorithmSel << std::endl;
    exit(1);
  }
  }

  // set up the merge sort class
  bfMergeSort<int64_t> testSort(1);

  // set up the random number ranges for test data and test size.
  auto riTestData = RandomInterval<int64_t>(-10000000000LL, 10000000000LL, 1);
  auto riSize = RandomInterval<int64_t>(1024LL, 1048576LL, 1);

  // set up a vector of number of treads to test and fill with a range of numbers.
  std::vector<size_t> threadList;
  for (size_t i = minThreads; i <= maxThreads; i++) threadList.push_back(i);
  
  // set up a pointer an array that will be the source for unsorted data. 
  // The same data will be used for multiple iterations of a test when a fixed test size is specified.
  int64_t* source_data = nullptr;
  // this is the checksum of the source data that will be used in the verification process.
  int64_t checksum = 0;

  // set up a counter for number of test tht failed.
  size_t testsFailed = 0;

  // start looping through the vector of the number of threads.
  for (size_t tc = 0; tc < threadList.size(); tc++) {
    size_t threads = threadList[tc];
    // set the number of threads
    testSort.setThreads(threads);

    // loop through the number of tests to run for each thread count.
    for (size_t test_num = 0; test_num < num_tests; test_num++) {
      int64_t test_size;
      if (fixedTestSize)  // set up either a fixed or random test size. 
        test_size = fixedTestSizeNum;
      else
        test_size = riSize();

      // Generate the test data but only once if it's a fixed data size
      if ((test_num == 0 && tc == 0) || !fixedTestSize) {
        if (source_data != nullptr) delete[] source_data;
        source_data = new int64_t[test_size];

        // create the requested data type.
        switch (dataType) {
        case dtRandom: { // generate random data
          for (size_t i = 0; i < test_size; i++) source_data[i] = riTestData();
          // make about 5% of the data equal at prime number spacings
          for (size_t i = 0; i < test_size; i += 19) source_data[i] = source_data[test_size - i];
          break;
        }
        case dtOrdered: {  // generate ordered data
          for (size_t i = 0; i < test_size; i++) source_data[i] = (i); 
          break;
        }
        case dtReverseOrdered: { // generate reverse ordered data
          for (size_t i = 0; i < test_size; i++) source_data[i] = (test_size - i);
          break;
        }
        default: {
          std::cout << "No such data type: " << dataType << std::endl;
          exit(1);
        }
        }

        if (verifySort) {  // if the verify option is set, calculate the checksum;
          checksum = 0;
          for (size_t i = 0; i < test_size; i++) checksum += source_data[i];
        }
      }

      // create the array to be sorted and copy the source data to it.
      int64_t* test_data = new int64_t[test_size];
      memcpy(test_data, source_data, test_size * sizeof(int64_t));

      // time the sort of one of the sort algorithms.
      // Get starting timepoint
      auto start = std::chrono::high_resolution_clock::now();
      // call the selected sort algorithm
      switch (algorithmSel) {
      case 1: {
        testSort.mSortFF<sortFor>(test_data, test_size); break;
      }
      case 2: {
        testSort.mSortFR<sortFor>(test_data, test_size); break;
      }
      case 3: {
        testSort.mSortRecFR<sortFor>(test_data, test_size); break;
      }
      case 4: {
        testSort.mSortRecFF<sortFor>(test_data, test_size); break;
      }
      case 5: {
        testSort.qSort<sortFor>(test_data, test_size); break;
      }
      case 6: {
        testSort.cSort<sortFor>(test_data, test_size); break;
      }
      case 11: {
        bfMergeSort<int64_t>::bfMergeSortMemFn f = &bfMergeSort<int64_t>::mSortFF<sortFor>;
        testSort.sortPerThread(test_data, test_size, f); break;
      }
      case 12: {
        bfMergeSort<int64_t>::bfMergeSortMemFn f = &bfMergeSort<int64_t>::mSortFR<sortFor>;
        testSort.sortPerThread(test_data, test_size, f); break;
      }
      case 13: {
        bfMergeSort<int64_t>::bfMergeSortMemFn f = &bfMergeSort<int64_t>::mSortRecFR<sortFor>;
        testSort.sortPerThread(test_data, test_size, f); break;
      }
      case 14: {
        bfMergeSort<int64_t>::bfMergeSortMemFn f = &bfMergeSort<int64_t>::mSortRecFF<sortFor>;
        testSort.sortPerThread(test_data, test_size, f); break;
      }
      case 15: {
        bfMergeSort<int64_t>::bfMergeSortMemFn f = &bfMergeSort<int64_t>::qSort<sortFor>;
        testSort.sortPerThread(test_data, test_size, f); break;
      }
      case 16: {
        bfMergeSort<int64_t>::bfMergeSortMemFn f = &bfMergeSort<int64_t>::cSort<sortFor>;
        testSort.sortPerThread(test_data, test_size, f); break;
      }
      default: {
        std::cout << "No such algorithm: " << algorithmSel << std::endl;
        exit(1);
      }
      }
      auto stop = std::chrono::high_resolution_clock::now();

      // calculate and print the execution time.
      auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
      std::cout << "Total sort time of " << std::setw(8) << test_size << " elements ";
      std::cout << "using " << threads << " threads = ";
      std::cout << std::fixed << std::setprecision(3) << duration.count()/1000000.0 << " seconds" << std::endl;

      if (verifySort) {
        // run the verify function in parallel because you can.  Needs and atomic error counter and checksum register.
        std::atomic<size_t> errCnt{ 0 };
        std::atomic<int64_t> postChecksum{ test_data[0] };

        // The verify portion makes sure the data is sorted correctly by checking that each array element is >= the one before it.
        // It also computes a checksum which is compared against the source data's checksum to check for data corruption.
        parallelFor(1, test_size, [&errCnt, &postChecksum, test_data](size_t lc) {
          static std::mutex lock;
          std::lock_guard<std::mutex> guard(lock);
          if (!(test_data[lc] >= test_data[lc - 1])) {
            if (errCnt == 0) {
              std::cout << "First error at " << lc << std::endl;
              std::cout << "[" << lc << "] " << test_data[lc] << ": [" << lc - 1 << "] " << test_data[lc - 1] << std::endl;
            }
            errCnt++;
          }
          postChecksum += test_data[lc];
          }, 8);
        bool thisTestFailed = false; 
        if (errCnt > 0) {
          std::cout << "Total of " << errCnt << " errors out of " << test_size << std::endl;
          thisTestFailed = true;
        }
        if (checksum != postChecksum) {
          std::cout << " Sorted data is different than source data." << std::endl;
          std::cout << postChecksum << " " << checksum << std::endl;
          thisTestFailed = true;
        }
        if (thisTestFailed) testsFailed++;
      }
      delete[] test_data;
    }
  }
  if (source_data != nullptr) delete[] source_data;
  std::cout << "Completed " << num_tests * threadList.size() << " tests";
  if (verifySort) std::cout << " with " << testsFailed << " test failures." << std::endl;
  else std::cout << "." << std::endl;
}
