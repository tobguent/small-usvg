
// Gets the xy components of the vector
inline Matrix<Scalar, 2, 1> xy() const { return Matrix<Scalar, 2, 1>(this->x(), this->y()); }
// Gets the xyz components of the vector
inline Matrix<Scalar, 3, 1> xyz() const { return Matrix<Scalar, 3, 1>(this->x(), this->y(), this->z()); }

// Component-wise division
template <typename OtherDerived>
EIGEN_STRONG_INLINE const CwiseBinaryOp<internal::scalar_quotient_op<Scalar>, const Derived, const OtherDerived>
operator/(const MatrixBase<OtherDerived>& other) const
{
    return cwiseQuotient(other);
}
