/**
 * @file ndarray.h
 * @brief A simple N-dimensional array in Column Major
 * @author Weiqi Shen weiqishen1994@ufl.edu
 * @version 0.1
 * @date 2019-01-26
 * 
 * @copyright Copyright (c) 2019
 * 
 */

#pragma once
#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <iomanip>
#include "err.h"

/**
 * @brief A simple N-dimensional array in Column Major
 * @class ndarray
 * @tparam T data type
 */
template <typename T>
class ndarray
{
public:
  // #### constructors ####

  /// Default constructor to construct an empty ndarray
  ndarray(void);

  /**
   * @brief N-Dimension constructor
   * 
   * @param list A list to specify dimension of the ndarray
   */
  ndarray(std::initializer_list<size_t> list);

  /**
   * @brief 1-D constructor
   * 
   * @param nele number of element
   */
  ndarray(size_t nele);

  /// Copy constructor
  ndarray(const ndarray<T> &in_array);

  /// destructor
  ~ndarray(void);

  // #### methods ####

  /**
   * @brief Overload stream output operator
   * 
   * @param out ostream object
   * @param s array to output
   * @return ostream& reference to the ostream object
   */
  template<class U>
  friend std::ostream& operator<<(std::ostream& out, ndarray<U>& s);

  /**
   * @brief setup a new ndarray object
   * 
   * @param list A list to specify dimension of the ndarray
   */
  void setup(std::initializer_list<size_t> list);

  /**
   * @brief setup a new 1-D ndarray object
   * 
   * @param list A list to specify dimension of the ndarray
   */
  void setup(size_t nele);

  ///Assignment
  ndarray<T> &operator=(const ndarray<T> &in_array);

  ///fill with value
  ndarray<T> &operator=(const T &other);

  /// Access/set ndarray element
  T &operator()(std::initializer_list<size_t> list);

  /// 1-D access/set
  T &operator()(size_t idx);

  /// Return pointer of ndarray element
  T *get_ptr(std::initializer_list<size_t> list);

  /// Return pointer of ndarray element in 1-D
  T *get_ptr(size_t idx = 0);

  /// Get number of elements along one axis
  size_t get_dim(size_t in_dim);

  /// Get number of dimension of the ndarray
  size_t get_n_dim(void);

  ///get the length of the ndarray
  size_t get_len(void);

  /// Method to get maximum value of ndarray
  T get_max(void);

  /// Method to get minimum value of ndarray
  T get_min(void);

  /// Reshape the array
  void reshape(std::initializer_list<size_t> list);

  ///inplace transpose
  void trans(void);

protected:
  size_t *shape;
  T *data;
  size_t n_dim;
  size_t len;

private:
  /// helper method to calculate length of ndarray
  void calc_len(void);
};

// definitions

using namespace std;

// #### constructors ####

// default constructor

template <typename T>
ndarray<T>::ndarray()
{
  len = 0;
  n_dim = 0;
  shape = NULL;
  data = NULL;
}

// constructor 1

template <typename T>
ndarray<T>::ndarray(initializer_list<size_t> list)
{
  //store dimension array
  n_dim = list.size();
  shape = new size_t[n_dim];

  size_t i = 0;
  for (auto l : list)
    shape[i++] = l;

  calc_len();
  data = new T[len];
}

// constructor 2

template <typename T>
ndarray<T>::ndarray(size_t nele)
{
  //store dimension array
  n_dim = 1;
  shape = new size_t[n_dim];
  shape[0] = nele;
  calc_len();
  data = new T[len];
}

// copy constructor

template <typename T>
ndarray<T>::ndarray(const ndarray<T> &in_array)
{
  n_dim = in_array.n_dim;
  shape = new size_t[n_dim];
  copy(in_array.shape, in_array.shape + n_dim, this->shape);

  calc_len();
  data = new T[len];
  copy(in_array.data, in_array.data + len, this->data);
}

// assignment

template <typename T>
ndarray<T> &ndarray<T>::operator=(const ndarray<T> &in_array)
{

  if (this == &in_array)
  {
    return (*this);
  }
  else
  {
    delete[] data;
    delete[] shape;

    n_dim = in_array.n_dim;
    shape = new size_t[n_dim];
    copy(in_array.shape, in_array.shape + n_dim, this->shape);

    calc_len();
    data = new T[len];
    copy(in_array.data, in_array.data + len, this->data);
    return (*this);
  }
}

template <typename T>
ndarray<T> &ndarray<T>::operator=(const T &other)
{
  fill_n(this->data, len, other);
  return *this;
}
// destructor

template <typename T>
ndarray<T>::~ndarray()
{
  delete[] data;
  delete[] shape;
}

// #### methods ####

// setup

template <typename T>
void ndarray<T>::setup(initializer_list<size_t> list)
{
  //delete previous data
  delete[] shape;
  delete[] data;

  //store dimension array
  n_dim = list.size();
  shape = new size_t[n_dim];

  size_t i = 0;
  for (auto l : list)
    shape[i++] = l;

  calc_len();
  data = new T[len];
}

//1-D setup
template <typename T>
void ndarray<T>::setup(size_t nele)
{
  //delete previous data
  delete[] shape;
  delete[] data;

  //store dimension array
  n_dim = 1;
  shape = new size_t[n_dim];
  shape[0]=nele;
  calc_len();
  data = new T[len];
}

template <typename T>
T &ndarray<T>::operator()(size_t idx)
{
  return data[idx];
}

template <typename T>
T &ndarray<T>::operator()(initializer_list<size_t> list)
{
  size_t idx = 0, acc = 1;
  size_t i = 0;

    for (auto l : list)
    {
      idx += acc * l;
      acc *= shape[i++];
    }
#ifdef _DEBUG
      if (idx >= len)
        Fatal_Error("Out of bound");
#endif
  return data[idx];
}

// return pointer

template <typename T>
T *ndarray<T>::get_ptr(size_t idx)
{
  return data + idx;
}

template <typename T>
T *ndarray<T>::get_ptr(initializer_list<size_t> list)
{
  size_t idx = 0, acc = 1;
  size_t i = 0;

  for (auto l : list)
  {
    idx += acc * l;
    acc *= shape[i++];
  }

  return data + idx;
}

// obtain dimension

template <typename T>
size_t ndarray<T>::get_dim(size_t in_dim)
{
  if (in_dim < n_dim)
    return shape[in_dim];
  else
    Fatal_Error("Dimension not supported");
}

template <typename T>
size_t ndarray<T>::get_n_dim()
{
  return n_dim;
}

template <typename T>
void ndarray<T>::calc_len(void)
{
  len = 1;
  for (size_t i = 0; i < n_dim; i++)
    len *= shape[i];
}

template <typename T>
size_t ndarray<T>::get_len(void)
{
  return len;
}
// method to calculate maximum value of ndarray
// Template specialization
template <typename T>
T ndarray<T>::get_max(void)
{
  return *max_element(data, data + len);
}

// method to calculate minimum value of ndarray
// Template specialization
template <typename T>
T ndarray<T>::get_min(void)
{
  return *min_element(data, data + len);
}

template <typename T>
void ndarray<T>::reshape(initializer_list<size_t> list)
{

  size_t i = 0;

  if (n_dim != list.size())
  {
    n_dim = list.size();
    //delete shape
    delete[] shape;
    shape=new size_t[n_dim];
  }
  
#ifdef _DEBUG
  size_t acc = 1;
#endif

  for (auto l : list)
  {
#ifdef _DEBUG
    acc *= l;
#endif
    shape[i++] = l;
  }

#ifdef _DEBUG
  if (acc != len)
    Fatal_Error("Total number of element doesn't agree");
#endif
}

template <typename T>
void ndarray<T>::trans()
{
  if (n_dim == 2)
  {
    size_t c_s = shape[0];
    size_t r_s = shape[1];

    for (size_t k = 0; k < shape[0] * shape[1]; k++) //loop over the target array
    {
      size_t idx = k;
      do
      {
        size_t id_0t = idx % r_s;
        size_t id_1t = idx / r_s;
        idx = id_0t * c_s + id_1t; //idx in original array
      } while (idx < k);
      swap(data[k], data[idx]);
    }
    swap(shape[0], shape[1]);
  }
  else
  {
    Fatal_Error("ERROR: Array transpose function only accepts a 2-dimensional square ndarray");
  }
}

//-------------friend---------------------
//output
template <typename U>
ostream &operator<<(ostream &out, ndarray<U> &s)
{
  if (s.n_dim == 1) //1-D output
  {
    for (size_t i = 0; i < s.shape[0]; i++)
      out << setprecision(4) << setw(10) << s(i);
    out << endl;
  }
  else if (s.n_dim == 2) //2d output
  {
    for (size_t i = 0; i < s.shape[0]; i++)
    {
      for (size_t j = 0; j < s.shape[1]; j++)
      {
        out << setprecision(4) << setw(10) << s({i, j});
      }
      out << endl;
    }
  }
  else if (s.n_dim == 3) //3d output
  {
    for (size_t k = 0; k < s.shape[2]; k++)
    {
      out << "slice: (:,:," << k << ")" << endl;
      for (size_t i = 0; i < s.shape[0]; i++)
      {
        for (size_t j = 0; j < s.shape[1]; j++)
        {
          out << setprecision(4) << setw(10) << s({i, j, k});
        }
        out << endl;
      }
    }
  }
  else
  {
    Fatal_Error("Output not supported");
  }

  return out;
}
