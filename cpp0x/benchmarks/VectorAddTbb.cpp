#include "gaalet.h"
#include <sys/time.h>
#include <cmath>
#include <tbb/task.h>
#include <google/profiler.h>

typedef gaalet::mv<1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12>::type Multivector;

template<typename T>
T na(const T& l, const T& r)
{
   T result;
   for(int i=0; i<T::size; ++i) {
      result[i] = l[i] + r[i];
   }
   return result;
}
template<typename T>
T ns(const T& l, const T& r)
{
   T result;
   for(int i=0; i<T::size; ++i) {
      result[i] = l[i] - r[i];
   }
   return result;
}

class SumTask : public tbb::task
{
public:
   SumTask(Multivector& c, const Multivector& a, const Multivector& b, const int& idx)
      :  c_(c),
         a_(a),
         b_(b),
         idx_(idx)
   {  }

   task* execute() {
      c_[idx_] = c_[idx_] + a_[idx_] + b_[idx_] - a_[idx_] - b_[idx_] + a_[idx_] + b_[idx_] - a_[idx_] - b_[idx_] - c_[idx_];
      return NULL;
   }

protected:
   Multivector& c_;
   const Multivector& a_;
   const Multivector& b_;
   const int idx_;
};

int main()
{
   timeval start, end;
   double solveTime;

   Multivector a = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
   Multivector b = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};

   Multivector c;

   gettimeofday(&start, 0);
   for(int i = 0; i<1e7; ++i) {
      c = c + a + b - a - b + a + b - a - b - c;
   }
   gettimeofday(&end, 0);
   solveTime = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)*1e-6;

   std::cout << "a: " << a << std::endl;
   std::cout << "b: " << b << std::endl;
   std::cout << "c: " << c << std::endl;

   std::cout << "operator=(): add solve time: " << solveTime << std::endl;


   gettimeofday(&start, 0);
   for(int i = 0; i<1e7; ++i) {
      c = eval(c + a + b - a - b + a + b - a - b - c);
   }
   gettimeofday(&end, 0);
   solveTime = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)*1e-6;

   std::cout << "a: " << a << std::endl;
   std::cout << "b: " << b << std::endl;
   std::cout << "c: " << c << std::endl;

   std::cout << "eval(): add solve time: " << solveTime << std::endl;


   gettimeofday(&start, 0);
   for(int i = 0; i<1e7; ++i) {
      c.assign(c + a + b - a - b + a + b - a - b - c);
   }
   gettimeofday(&end, 0);
   solveTime = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)*1e-6;

   std::cout << "a: " << a << std::endl;
   std::cout << "b: " << b << std::endl;
   std::cout << "c: " << c << std::endl;

   std::cout << "assign(): add solve time: " << solveTime << std::endl;


   gettimeofday(&start, 0);
   for(int i = 0; i<1e7; ++i) {
      c[0] = c[0] + a[0] + b[0] - a[0] - b[0] + a[0] + b[0] - a[0] - b[0] - c[0];
      c[1] = c[1] + a[1] + b[1] - a[1] - b[1] + a[1] + b[1] - a[1] - b[1] - c[1];
      c[1] = c[1] + a[1] + b[1] - a[1] - b[1] + a[1] + b[1] - a[1] - b[1] - c[1];
      c[2] = c[2] + a[2] + b[2] - a[2] - b[2] + a[2] + b[2] - a[2] - b[2] - c[2];
      c[3] = c[3] + a[3] + b[3] - a[3] - b[3] + a[3] + b[3] - a[3] - b[3] - c[3];
      c[4] = c[4] + a[4] + b[4] - a[4] - b[4] + a[4] + b[4] - a[4] - b[4] - c[4];
      c[5] = c[5] + a[5] + b[5] - a[5] - b[5] + a[5] + b[5] - a[5] - b[5] - c[5];
      c[6] = c[6] + a[6] + b[6] - a[6] - b[6] + a[6] + b[6] - a[6] - b[6] - c[6];
      c[7] = c[7] + a[7] + b[7] - a[7] - b[7] + a[7] + b[7] - a[7] - b[7] - c[7];
      c[8] = c[8] + a[8] + b[8] - a[8] - b[8] + a[8] + b[8] - a[8] - b[8] - c[8];
      c[9] = c[9] + a[9] + b[9] - a[9] - b[9] + a[9] + b[9] - a[9] - b[9] - c[9];
      c[10] = c[10] + a[10] + b[10] - a[10] - b[10] + a[10] + b[10] - a[10] - b[10] - c[10];
      c[11] = c[11] + a[11] + b[11] - a[11] - b[11] + a[11] + b[11] - a[11] - b[11] - c[11];
   }
   gettimeofday(&end, 0);
   solveTime = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)*1e-6;

   std::cout << "a: " << a << std::endl;
   std::cout << "b: " << b << std::endl;
   std::cout << "c: " << c << std::endl;

   std::cout << "handcoded: add solve time: " << solveTime << std::endl;

   //ProfilerStart("VectorAddTbb.prof");
   gettimeofday(&start, 0);
   for(int i = 0; i<1e7; ++i) {
      tbb::task_list tl;
      for(int j=0; j<12; ++j) tl.push_back(*new(tbb::task::allocate_root()) SumTask(c,a,b,j));
      tbb::task::spawn_root_and_wait(tl);
   }
   gettimeofday(&end, 0);
   //ProfilerStop();
   solveTime = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)*1e-6;

   std::cout << "a: " << a << std::endl;
   std::cout << "b: " << b << std::endl;
   std::cout << "c: " << c << std::endl;

   std::cout << "tbb handcoded: add solve time: " << solveTime << std::endl;


   gettimeofday(&start, 0);
   for(int i = 0; i<1e7; ++i) {
      //c = na(c, na(a, ns(b, ns(a, na(b, na(a, ns(b, ns(a, ns(b, c)))))))));
      c = ns(ns(ns(na(na(ns(ns(na(na(c, a), b), a), b), a), b), a), b), c);
   }
   gettimeofday(&end, 0);
   solveTime = (double)(end.tv_sec - start.tv_sec) + (double)(end.tv_usec - start.tv_usec)*1e-6;

   std::cout << "a: " << a << std::endl;
   std::cout << "b: " << b << std::endl;
   std::cout << "c: " << c << std::endl;

   std::cout << "naive_addition::operator+(): add solve time: " << solveTime << std::endl;
}
