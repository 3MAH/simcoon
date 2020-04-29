/* Takes a derivative of a Tensor2_symmetric, yielding a Dg. */

#pragma once

namespace FTensor
{
  template <class T, int Dim01, int Dim2, char i, char j, char k>
  class dTensor2_symmetric
  {
    const Tensor2_symmetric<T *, Dim01> &a;
    const Tensor1<int, Dim2> &d_ijk;
    const Tensor1<double, Dim2> &d_xyz;

  public:
    typename promote<T, double>::V
    operator()(const int N1, const int N2, const int N3) const
    {
      return (*(a.ptr(N1, N2) + d_ijk(N3)) - *(a.ptr(N1, N2) - d_ijk(N3)))
             * d_xyz(N3) * 0.5;
    }
    dTensor2_symmetric(const Tensor2_symmetric<T *, Dim01> &A,
                       const Tensor1<int, Dim2> &D_ijk,
                       const Tensor1<double, Dim2> &D_xyz)
        : a(A), d_ijk(D_ijk), d_xyz(D_xyz)
    {}
  };

  template <class T, int Dim01, int Dim2, char i, char j, char k>
  const Dg_Expr<const dTensor2_symmetric<T, Dim01, Dim2, i, j, k>,
                typename promote<T, double>::V, Dim01, Dim2, i, j, k>
  d(const Tensor2_symmetric<T *, Dim01> &a, const Index<i, Dim01> index1,
    const Index<j, Dim01> index2, const Index<k, Dim2> index3,
    const Tensor1<int, Dim2> &d_ijk, const Tensor1<double, Dim2> &d_xyz)
  {
    using TensorExpr = dTensor2_symmetric<T, Dim01, Dim2, i, j, k>;
    return Dg_Expr<TensorExpr, typename promote<T, double>::V, Dim01, Dim2, i,
                   j, k>(TensorExpr(a, d_ijk, d_xyz));
  }
}
