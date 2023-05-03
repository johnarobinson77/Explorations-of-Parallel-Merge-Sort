#ifndef PARALLEL_FOR
#define PARALLEL_FOR

/**
 * Copyright (c) 2023 John Robinson.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */


#include <iostream>
#include <future>
#include <mutex>
#include <vector>
#include <numeric>
#include <atomic>
#include <iomanip>
#include <random>

#define NDEBUG

// parallel_for runs the equivalent to the statement "for (int i = begin; i < end; i++) fn(i);"
// where fn is a lambda.  if num_segs is greater than 1, it divides the the begin-end number range 
// into that many  segments and start starts separate threads to execute the first num_segs-1 segments
// The last segment is run on the current thread.  
// parallel_for returns after all treads segments have completed.
template < typename FN >
void parallelFor(size_t begin, size_t end, FN fn, int64_t num_segs = 1) {

  // compute the number iterations and if 0, then return
  const std::size_t n = end - begin;
  if (n == 0) return;
  // if the number of elements is less than the number of segments or threads, reduce the number of segments
  if (n < num_segs) num_segs = n;
  // create the array of futures.
  std::vector<std::future<void>> futures;
  // Compute the number of iterations to preform on each segment using a double to preserve the fraction.
  // using the method of iterating through bounds applied to each thread using a double and then rounding to
  // the nearest integer ensures that the difference in the number of iterations that each thread will do is 
  // no greater than one.
  double seg_size = (double)n / (double)num_segs;
  size_t sBeg = 0;
  double sEnd = seg_size;
  // start threads for one less than the number of segments
  for (int64_t seg = 0; seg < num_segs-1; seg++) { 
    futures.push_back(std::async(std::launch::async, [&](size_t lb, size_t le) {
#ifndef NDEBUG
      {
        static std::mutex lock;
        std::lock_guard<std::mutex> guard(lock);
        std::cout << "parallel_for: n == " << le - lb << " thread: " << std::this_thread::get_id() << '\n';
      }
#endif // NDEBUG
      for (size_t i = lb; i < le; i++) fn(i);
      }, begin + sBeg, begin + llround(sEnd)));
    sBeg = llround(sEnd);
    sEnd += seg_size;
  }
  // execute the last segment in this thread
  for (size_t i = begin + sBeg; i < begin + llround(sEnd); i++) fn(i);

  // wait for the segments to complete
  for (int seg = 0; seg < futures.size(); seg++)
    futures[seg].wait();
}

// parallelFor runs the equivalent to the statement "for (int i = begin; i < end; i++) fn(i);"
// where fn is a lambda.  if num_segs is greater than 1, it divides the the begin-end number range 
// into that many  segments and start starts separate threads to execute the first num_segs-1 segments
// The last segment is run on the current thread.  
// parallelForBoWait returns before checking that all thread are done.  parallelForFinish() must be
// called make sure all of the threads have finished.  
template < typename FN >
std::vector<std::future<void>>*  parallelForNoWait(size_t begin, size_t end, FN fn, int64_t num_segs = 1) {
  
  auto* futuresNW = new std::vector<std::future<void>>;

  // compute the number iterations and if 0, then return
  const std::size_t n = end - begin;
  if (n == 0) return nullptr;
  // if the number of elements is less than the number of segments or threads, reduce the number of segments
  if (n < num_segs) num_segs = n;
  // Compute the number of iterations to preform on each segment using a double to preserve the fraction.
  // using the method of iterating through bounds applied to each thread using a double and then rounding to
  // the nearest integer ensures that the difference in the number of iterations that each thread will do is 
  // no greater than one.
  double seg_size = (double)n / (double)num_segs;
  size_t sBeg = 0;
  double sEnd = seg_size;
  // start threads for one less than the number of segments
  for (int64_t seg = 0; seg < num_segs - 1; seg++) {
    futuresNW->push_back(std::async(std::launch::async, [&](size_t lb, size_t le) {
      for (size_t i = lb; i < le; i++) fn(i);
      }, begin + sBeg, begin + llround(sEnd)));
    sBeg = llround(sEnd);
    sEnd += seg_size;
    }
  // execute the last segment in this thread
  for (size_t i = begin + sBeg; i < begin + llround(sEnd); i++) fn(i);
  return futuresNW;
  }

void parallelForFinish(std::vector<std::future<void>>* futuresNW) {
  if (futuresNW == nullptr) return;
  // wait for the segments to complete
  for (size_t seg = 0; seg < futuresNW->size(); seg++)
    futuresNW->at(seg).wait();
  delete futuresNW;
}



#endif // PARALLEL_FOR
