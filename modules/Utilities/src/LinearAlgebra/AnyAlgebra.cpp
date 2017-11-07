#include "MUQ/Utilities/LinearAlgebra/AnyAlgebra.h"

using namespace muq::Utilities;

AnyAlgebra::AnyAlgebra() {}

bool AnyAlgebra::IsScalar(std::type_info const& obj_type) const {
  // is this a scalar type?
  return typeid(double)==obj_type || typeid(float)==obj_type || typeid(int)==obj_type || typeid(unsigned int)==obj_type;
}

bool AnyAlgebra::IsEigenVector(std::type_info const& obj_type) const {
  // is this an Eigen::Vector type?
  return typeid(Eigen::Vector2d)==obj_type
    || typeid(Eigen::Vector2f)==obj_type
    || typeid(Eigen::Vector2i)==obj_type
    || typeid(Eigen::Vector3d)==obj_type
    || typeid(Eigen::Vector3f)==obj_type
    || typeid(Eigen::Vector3i)==obj_type
    || typeid(Eigen::Vector4d)==obj_type
    || typeid(Eigen::Vector4f)==obj_type
    || typeid(Eigen::Vector4i)==obj_type
    || typeid(Eigen::VectorXd)==obj_type
    || typeid(Eigen::VectorXf)==obj_type
    || typeid(Eigen::VectorXi)==obj_type;
}

bool AnyAlgebra::IsEigenMatrix(std::type_info const& obj_type) const {
  // is this an Eigen::Matrix type?
  return typeid(Eigen::Matrix2d)==obj_type
    || typeid(Eigen::Matrix2f)==obj_type
    || typeid(Eigen::Matrix2i)==obj_type
    || typeid(Eigen::Matrix3d)==obj_type
    || typeid(Eigen::Matrix3f)==obj_type
    || typeid(Eigen::Matrix3i)==obj_type
    || typeid(Eigen::Matrix4d)==obj_type
    || typeid(Eigen::Matrix4f)==obj_type
    || typeid(Eigen::Matrix4i)==obj_type
    || typeid(Eigen::MatrixXd)==obj_type
    || typeid(Eigen::MatrixXf)==obj_type
    || typeid(Eigen::MatrixXi)==obj_type;
}

#if MUQ_HAS_SUNDIALS==1
bool AnyAlgebra::IsSundialsVector(std::type_info const& obj_type) const {
  return typeid(N_Vector)==obj_type;
}
#endif

unsigned int AnyAlgebra::EigenVectorSize(boost::any const& vec) const {
  // fixed size vectors
  if( typeid(Eigen::Vector2d)==vec.type() ) { return 2; }
  if( typeid(Eigen::Vector3d)==vec.type() ) { return 3; }
  if( typeid(Eigen::Vector4d)==vec.type() ) { return 4; }

  // generic Eigen::VectorXd
  const Eigen::VectorXd& eig = boost::any_cast<Eigen::VectorXd const&>(vec);
  return eig.size();
}

unsigned int AnyAlgebra::EigenMatrixSize(boost::any const& mat, int const dim) const {
  // fixed size vectors
  if( typeid(Eigen::Matrix2d)==mat.type() ) { return dim==-1? 4 : 2;}
  if( typeid(Eigen::Matrix3d)==mat.type() ) { return dim==-1? 9 : 3; }
  if( typeid(Eigen::Matrix4d)==mat.type() ) { return dim==-1? 16 : 4; }

  // generic Eigen::MatrixXd
  const Eigen::MatrixXd& eig = boost::any_cast<Eigen::MatrixXd const&>(mat);
  switch( dim ) {
  case 0:
    return eig.rows();
  case 1:
    return eig.cols();
  default:
    return eig.size();
  }
}

#if MUQ_HAS_SUNDIALS==1
unsigned int AnyAlgebra::SundialsVectorSize(boost::any const& vec) const {
  const N_Vector& sun = boost::any_cast<N_Vector const&>(vec);

  return NV_LENGTH_S(sun);
}
#endif

unsigned int AnyAlgebra::Size(boost::any const& obj, int const dim) const {
  // scalars are size one
  if( IsScalar(obj.type()) ) { return 1; }

  // get the size of an Eigen::VectorXd
  if( IsEigenVector(obj.type()) ) { return EigenVectorSize(obj); }

  // get the size of an Eigen::MatrixXd
  if( IsEigenMatrix(obj.type()) ) { return EigenMatrixSize(obj, dim); }

  // get the size of Sundials vectors
  if( IsSundialsVector(obj.type()) ) { return SundialsVectorSize(obj); }

  return SizeImpl(obj);
}

unsigned int AnyAlgebra::SizeImpl(boost::any const& obj) const {
  std::cerr << std::endl << "ERROR: cannot compute the size of an object with type " << boost::core::demangle(obj.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::SizeImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::SizeImpl()" << std::endl << std::endl;
  assert(false);

  return 0;
}

boost::any AnyAlgebra::ZeroEigenVector(std::type_info const& type, unsigned int const size) const {  
  if( typeid(Eigen::Vector2d)==type ) { return (Eigen::Vector2d)Eigen::Vector2d::Zero(); }
  if( typeid(Eigen::Vector2f)==type ) { return (Eigen::Vector2f)Eigen::Vector2f::Zero(); }
  if( typeid(Eigen::Vector2i)==type ) { return (Eigen::Vector2i)Eigen::Vector2i::Zero(); }
  
  if( typeid(Eigen::Vector3d)==type ) { return (Eigen::Vector3d)Eigen::Vector3d::Zero(); }
  if( typeid(Eigen::Vector3f)==type ) { return (Eigen::Vector3f)Eigen::Vector3f::Zero(); }
  if( typeid(Eigen::Vector3i)==type ) { return (Eigen::Vector3i)Eigen::Vector3i::Zero(); }
  
  if( typeid(Eigen::Vector4d)==type ) { return (Eigen::Vector4d)Eigen::Vector4d::Zero(); }
  if( typeid(Eigen::Vector4f)==type ) { return (Eigen::Vector4f)Eigen::Vector4f::Zero(); }
  if( typeid(Eigen::Vector4i)==type ) { return (Eigen::Vector4i)Eigen::Vector4i::Zero(); }

  if( typeid(Eigen::VectorXd)==type ) { return (Eigen::VectorXd)Eigen::VectorXd::Zero(size); }
  if( typeid(Eigen::VectorXf)==type ) { return (Eigen::VectorXf)Eigen::VectorXf::Zero(size); }
  if( typeid(Eigen::VectorXi)==type ) { return (Eigen::VectorXi)Eigen::VectorXi::Zero(size); }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any AnyAlgebra::ZeroEigenMatrix(std::type_info const& type, unsigned int const rows, unsigned int const cols) const {
  if( typeid(Eigen::Matrix2d)==type ) { return (Eigen::Matrix2d)Eigen::Matrix2d::Zero(); }
  if( typeid(Eigen::Matrix2f)==type ) { return (Eigen::Matrix2f)Eigen::Matrix2f::Zero(); }
  if( typeid(Eigen::Matrix2i)==type ) { return (Eigen::Matrix2i)Eigen::Matrix2i::Zero(); }
  
  if( typeid(Eigen::Matrix3d)==type ) { return (Eigen::Matrix3d)Eigen::Matrix3d::Zero(); }
  if( typeid(Eigen::Matrix3f)==type ) { return (Eigen::Matrix3f)Eigen::Matrix3f::Zero(); }
  if( typeid(Eigen::Matrix3i)==type ) { return (Eigen::Matrix3i)Eigen::Matrix3i::Zero(); }
  
  if( typeid(Eigen::Matrix4d)==type ) { return (Eigen::Matrix4d)Eigen::Matrix4d::Zero(); }
  if( typeid(Eigen::Matrix4f)==type ) { return (Eigen::Matrix4f)Eigen::Matrix4f::Zero(); }
  if( typeid(Eigen::Matrix4i)==type ) { return (Eigen::Matrix4i)Eigen::Matrix4i::Zero(); }

  if( typeid(Eigen::MatrixXd)==type ) { return (Eigen::MatrixXd)Eigen::MatrixXd::Zero(rows, cols); }
  if( typeid(Eigen::MatrixXf)==type ) { return (Eigen::MatrixXf)Eigen::MatrixXf::Zero(rows, cols); }
  if( typeid(Eigen::MatrixXi)==type ) { return (Eigen::MatrixXi)Eigen::MatrixXi::Zero(rows, cols); }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any AnyAlgebra::ZeroScalar(std::type_info const& type) const {
  if( typeid(double)==type ) { return (double)0.0; }
  if( typeid(float)==type ) { return (float)0.0; }
  if( typeid(int)==type ) { return (int)0; }
  if( typeid(unsigned int)==type ) { return (unsigned int)0; }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any AnyAlgebra::Zero(std::type_info const& type, unsigned int rows, unsigned int const cols) const {
  if( IsScalar(type) ) { return ZeroScalar(type); }
  
  if( IsEigenVector(type) ) { return ZeroEigenVector(type, rows); }

  if( IsEigenMatrix(type) ) { return ZeroEigenMatrix(type, rows, cols); }

  if( type==typeid(double) ) {
    return 0.0;
  } 
  
  return ZeroImpl(type, rows, cols);
}

boost::any AnyAlgebra::ZeroImpl(std::type_info const& type, unsigned int const rows, unsigned int const cols) const {
  std::cerr << std::endl << "ERROR: cannot compute zero of an object with type " << boost::core::demangle(type.name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::ZeroImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::ZeroImpl()" << std::endl << std::endl;
  assert(false);

  return boost::none;
}

double AnyAlgebra::ScalarNorm(boost::any const& obj) const {
  if( typeid(double)==obj.type() ) { return ScalarMagnitude<double>(obj); }
  if( typeid(float)==obj.type() ) { return ScalarMagnitude<float>(obj); }
  if( typeid(int)==obj.type() ) { return ScalarMagnitude<int>(obj); }
  if( typeid(unsigned int)==obj.type() ) { return ScalarMagnitude<unsigned int>(obj); }

  // something went wrong
  assert(false);
  return -1.0;
}

double AnyAlgebra::EigenVectorNorm(boost::any const& vec) const {
  if( typeid(Eigen::Vector2d)==vec.type() ) { return EigenNorm<Eigen::Vector2d>(vec); }
  if( typeid(Eigen::Vector3d)==vec.type() ) { return EigenNorm<Eigen::Vector3d>(vec); }
  if( typeid(Eigen::Vector4d)==vec.type() ) { return EigenNorm<Eigen::Vector4d>(vec); }
  
  return EigenNorm<Eigen::VectorXd>(vec);
}

double AnyAlgebra::EigenMatrixNorm(boost::any const& mat) const {
  if( typeid(Eigen::Matrix2d)==mat.type() ) { return EigenNorm<Eigen::Matrix2d>(mat); }
  if( typeid(Eigen::Matrix3d)==mat.type() ) { return EigenNorm<Eigen::Matrix3d>(mat); }
  if( typeid(Eigen::Matrix4d)==mat.type() ) { return EigenNorm<Eigen::Matrix4d>(mat); }
  
  return EigenNorm<Eigen::MatrixXd>(mat);
}

double AnyAlgebra::Norm(boost::any const& obj) const {
  if( IsScalar(obj.type()) ) { return ScalarNorm(obj); }

  if( IsEigenVector(obj.type()) ) { return EigenVectorNorm(obj); }

  if( IsEigenMatrix(obj.type()) ) { return EigenMatrixNorm(obj); }
    
  return NormImpl(obj);
}

double AnyAlgebra::NormImpl(boost::any const& obj) const {
  std::cerr << std::endl << "ERROR: Cannot compute the norm of an object with type " << boost::core::demangle(obj.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::NormImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::NormImpl()" << std::endl << std::endl;
  assert(false);

  return -1.0;
}

double AnyAlgebra::ScalarInnerProduct(boost::any const& vec1, boost::any const& vec2) const {
  double ip = ScalarInProd<double>(vec1, vec2);
  if( !std::isnan(ip) ) { return ip; }

  ip = ScalarInProd<float>(vec1, vec2);
  if( !std::isnan(ip) ) { return ip; }

  ip = ScalarInProd<int>(vec1, vec2);
  if( !std::isnan(ip) ) { return ip; }

  ip = ScalarInProd<unsigned int>(vec1, vec2);
  if( !std::isnan(ip) ) { return ip; }

  // something went wrong
  return std::numeric_limits<double>::quiet_NaN();
}

double AnyAlgebra::EigenVectorInnerProduct(boost::any const& vec1, boost::any const& vec2) const {
  double ip = EigenVectorInProd<Eigen::Vector2d>(vec1, vec2);
  if( !std::isnan(ip) ) { return ip; }

  ip = EigenVectorInProd<Eigen::Vector3d>(vec1, vec2);
  if( !std::isnan(ip) ) { return ip; }

  ip = EigenVectorInProd<Eigen::Vector4d>(vec1, vec2);
  if( !std::isnan(ip) ) { return ip; }

  return EigenVectorInProd<Eigen::VectorXd, Eigen::VectorXd>(vec1, vec2);
}

double AnyAlgebra::InnerProduct(boost::any const& vec1, boost::any const& vec2) const {
  if( IsScalar(vec1.type()) && IsScalar(vec2.type()) ) { return ScalarInnerProduct(vec1, vec2); }

  if( IsEigenVector(vec1.type()) && IsEigenVector(vec2.type()) ) { return EigenVectorInnerProduct(vec1, vec2); }

  return InnerProductImpl(vec1, vec2);
}

double AnyAlgebra::InnerProductImpl(boost::any const& vec1, boost::any const& vec2) const {
  std::cerr << std::endl << "ERROR: Cannot compute the inner product between vectors with types " << boost::core::demangle(vec1.type().name()) << " and " << boost::core::demangle(vec2.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::InnerProductImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::InnerProductImpl()" << std::endl << std::endl;
  assert(false);

  return 0.0;
}

bool AnyAlgebra::IsScalarZero(boost::any const& obj) const {
  if( obj.type()==typeid(double) ) { return boost::any_cast<double const>(obj)==0.0; }
  if( obj.type()==typeid(float) ) { return boost::any_cast<float const>(obj)==0.0; }
  if( obj.type()==typeid(int) ) { return boost::any_cast<int const>(obj)==0; }
  if( obj.type()==typeid(unsigned int) ) { return boost::any_cast<unsigned int const>(obj)==0; }

  // something when wrong
  assert(false);
  return false;
}

bool AnyAlgebra::IsEigenVectorZero(boost::any const& obj) const {
  if( typeid(Eigen::Vector2d)==obj.type() ) { return IsEigVecZero<Eigen::Vector2d>(obj); }
  if( typeid(Eigen::Vector2f)==obj.type() ) { return IsEigVecZero<Eigen::Vector2f>(obj); }
  if( typeid(Eigen::Vector2i)==obj.type() ) { return IsEigVecZero<Eigen::Vector2i>(obj); }

  if( typeid(Eigen::Vector3d)==obj.type() ) { return IsEigVecZero<Eigen::Vector3d>(obj); }
  if( typeid(Eigen::Vector3f)==obj.type() ) { return IsEigVecZero<Eigen::Vector3f>(obj); }
  if( typeid(Eigen::Vector3i)==obj.type() ) { return IsEigVecZero<Eigen::Vector3i>(obj); }

  if( typeid(Eigen::Vector4d)==obj.type() ) { return IsEigVecZero<Eigen::Vector4d>(obj); }
  if( typeid(Eigen::Vector4f)==obj.type() ) { return IsEigVecZero<Eigen::Vector4f>(obj); }
  if( typeid(Eigen::Vector4i)==obj.type() ) { return IsEigVecZero<Eigen::Vector4i>(obj); }

  if( typeid(Eigen::VectorXd)==obj.type() ) { return IsEigVecZero<Eigen::VectorXd>(obj); }
  if( typeid(Eigen::VectorXf)==obj.type() ) { return IsEigVecZero<Eigen::VectorXf>(obj); }
  if( typeid(Eigen::VectorXi)==obj.type() ) { return IsEigVecZero<Eigen::VectorXi>(obj); }

  // something went wrong
  assert(false);
  return false;
}

bool AnyAlgebra::IsEigenMatrixZero(boost::any const& obj) const {
  if( typeid(Eigen::Matrix2d)==obj.type() ) { return IsEigMatZero<Eigen::Matrix2d>(obj); }
  if( typeid(Eigen::Matrix2f)==obj.type() ) { return IsEigMatZero<Eigen::Matrix2f>(obj); }
  if( typeid(Eigen::Matrix2i)==obj.type() ) { return IsEigMatZero<Eigen::Matrix2i>(obj); }

  if( typeid(Eigen::Matrix3d)==obj.type() ) { return IsEigMatZero<Eigen::Matrix3d>(obj); }
  if( typeid(Eigen::Matrix3f)==obj.type() ) { return IsEigMatZero<Eigen::Matrix3f>(obj); }
  if( typeid(Eigen::Matrix3i)==obj.type() ) { return IsEigMatZero<Eigen::Matrix3i>(obj); }

  if( typeid(Eigen::Matrix4d)==obj.type() ) { return IsEigMatZero<Eigen::Matrix4d>(obj); }
  if( typeid(Eigen::Matrix4f)==obj.type() ) { return IsEigMatZero<Eigen::Matrix4f>(obj); }
  if( typeid(Eigen::Matrix4i)==obj.type() ) { return IsEigMatZero<Eigen::Matrix4i>(obj); }

  if( typeid(Eigen::MatrixXd)==obj.type() ) { return IsEigMatZero<Eigen::MatrixXd>(obj); }
  if( typeid(Eigen::MatrixXf)==obj.type() ) { return IsEigMatZero<Eigen::MatrixXf>(obj); }
  if( typeid(Eigen::MatrixXi)==obj.type() ) { return IsEigMatZero<Eigen::MatrixXi>(obj); }

  // something went wrong
  assert(false);
  return false;
}

bool AnyAlgebra::IsZero(boost::any const& obj) const {
  if( IsScalar(obj.type()) ) { return IsScalarZero(obj); }

  if( IsEigenVector(obj.type()) ) { return IsEigenVectorZero(obj); }
  
  if( IsEigenMatrix(obj.type()) ) { return IsEigenMatrixZero(obj); }
  
  return IsZeroImpl(obj);
}

bool AnyAlgebra::IsZeroImpl(boost::any const& obj) const {
  std::cerr << std::endl << "ERROR: No way to determine if an object with type " << boost::core::demangle(obj.type().name()) << " is the zero vector." << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::IsZero()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::IsZero()" << std::endl << std::endl;
  assert(false);

  return false;
}

boost::any AnyAlgebra::AccessEigenMatrix(boost::any const& mat, unsigned int const i, unsigned int const j) const {
  if( typeid(Eigen::Matrix2d)==mat.type() ) { return AccessEigenMat<Eigen::Matrix2d>(mat, i, j); }
  if( typeid(Eigen::Matrix2f)==mat.type() ) { return AccessEigenMat<Eigen::Matrix2f>(mat, i, j); }
  if( typeid(Eigen::Matrix2i)==mat.type() ) { return AccessEigenMat<Eigen::Matrix2i>(mat, i, j); }

  if( typeid(Eigen::Matrix3d)==mat.type() ) { return AccessEigenMat<Eigen::Matrix3d>(mat, i, j); }
  if( typeid(Eigen::Matrix3f)==mat.type() ) { return AccessEigenMat<Eigen::Matrix3f>(mat, i, j); }
  if( typeid(Eigen::Matrix3i)==mat.type() ) { return AccessEigenMat<Eigen::Matrix3i>(mat, i, j); }

  if( typeid(Eigen::Matrix4d)==mat.type() ) { return AccessEigenMat<Eigen::Matrix4d>(mat, i, j); }
  if( typeid(Eigen::Matrix4f)==mat.type() ) { return AccessEigenMat<Eigen::Matrix4f>(mat, i, j); }
  if( typeid(Eigen::Matrix4i)==mat.type() ) { return AccessEigenMat<Eigen::Matrix4i>(mat, i, j); }

  if( typeid(Eigen::MatrixXd)==mat.type() ) { return AccessEigenMat<Eigen::MatrixXd>(mat, i, j); }
  if( typeid(Eigen::MatrixXf)==mat.type() ) { return AccessEigenMat<Eigen::MatrixXf>(mat, i, j); }
  if( typeid(Eigen::MatrixXi)==mat.type() ) { return AccessEigenMat<Eigen::MatrixXi>(mat, i, j); }

  // something went wront
  assert(false);
  return boost::none;
}

boost::any AnyAlgebra::AccessEigenVector(boost::any const& vec, unsigned int const i) const {
  if( typeid(Eigen::Vector2d)==vec.type() ) { return AccessEigenVec<Eigen::Vector2d>(vec, i); }
  if( typeid(Eigen::Vector2f)==vec.type() ) { return AccessEigenVec<Eigen::Vector2f>(vec, i); }
  if( typeid(Eigen::Vector2i)==vec.type() ) { return AccessEigenVec<Eigen::Vector2i>(vec, i); }

  if( typeid(Eigen::Vector3d)==vec.type() ) { return AccessEigenVec<Eigen::Vector3d>(vec, i); }
  if( typeid(Eigen::Vector3f)==vec.type() ) { return AccessEigenVec<Eigen::Vector3f>(vec, i); }
  if( typeid(Eigen::Vector3i)==vec.type() ) { return AccessEigenVec<Eigen::Vector3i>(vec, i); }

  if( typeid(Eigen::Vector4d)==vec.type() ) { return AccessEigenVec<Eigen::Vector4d>(vec, i); }
  if( typeid(Eigen::Vector4f)==vec.type() ) { return AccessEigenVec<Eigen::Vector4f>(vec, i); }
  if( typeid(Eigen::Vector4i)==vec.type() ) { return AccessEigenVec<Eigen::Vector4i>(vec, i); }

  if( typeid(Eigen::VectorXd)==vec.type() ) { return AccessEigenVec<Eigen::VectorXd>(vec, i); }
  if( typeid(Eigen::VectorXf)==vec.type() ) { return AccessEigenVec<Eigen::VectorXf>(vec, i); }
  if( typeid(Eigen::VectorXi)==vec.type() ) { return AccessEigenVec<Eigen::VectorXi>(vec, i); }

  // something went wront
  assert(false);
  return boost::none;
}

#if MUQ_HAS_SUNDIALS==1
boost::any AnyAlgebra::AccessSundialsVector(N_Vector const& vec, unsigned int const i) const {
  // check the size
    assert(i<NV_LENGTH_S(vec));

    // return the ith element
    return NV_Ith_S(vec, i);
}
#endif

boost::any AnyAlgebra::AccessElement(boost::any const& obj, unsigned int const i, unsigned int const j) const {
  if( IsScalar(obj.type()) ) { return obj; }

  if( IsEigenVector(obj.type()) ) { return AccessEigenVector(obj, i); }

  if( IsEigenMatrix(obj.type()) ) { return AccessEigenMatrix(obj, i, j); }

#if MUQ_HAS_SUNDIALS==1
  if( IsSundialsVector(obj.type()) ) { return AccessSundialsVector(boost::any_cast<const N_Vector&>(obj), i); }
#endif
  
  return AccessElementImpl(obj, i);
}

boost::any AnyAlgebra::AccessElementImpl(boost::any const& vec, unsigned int const i) const {
  std::cerr << std::endl << "ERROR: No way to access element " << i << " of a vector with type " << boost::core::demangle(vec.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::AccessElement()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::AccessElement()" << std::endl << std::endl;
  assert(false);

  return boost::none;
}

boost::any AnyAlgebra::ScalarIdentity(std::type_info const& type) const {
  if( type==typeid(double) ) { return (double)1.0; }
  if( type==typeid(float) ) { return (float)1.0; }
  if( type==typeid(int) ) { return (int)1; }
  if( type==typeid(unsigned int) ) { return (unsigned int)1; }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any AnyAlgebra::EigenMatrixIdentity(std::type_info const& type, unsigned int const rows, unsigned int const cols) const {
  if( type==typeid(Eigen::Matrix2d) ) { return (Eigen::Matrix2d)Eigen::Matrix2d::Identity(); }
  if( type==typeid(Eigen::Matrix2f) ) { return (Eigen::Matrix2f)Eigen::Matrix2f::Identity(); }
  if( type==typeid(Eigen::Matrix2i) ) { return (Eigen::Matrix2i)Eigen::Matrix2i::Identity(); }
  
  if( type==typeid(Eigen::Matrix3d) ) { return (Eigen::Matrix3d)Eigen::Matrix3d::Identity(); }
  if( type==typeid(Eigen::Matrix3f) ) { return (Eigen::Matrix3f)Eigen::Matrix3f::Identity(); }
  if( type==typeid(Eigen::Matrix3i) ) { return (Eigen::Matrix3i)Eigen::Matrix3i::Identity(); }
  
  if( type==typeid(Eigen::Matrix4d) ) { return (Eigen::Matrix4d)Eigen::Matrix4d::Identity(); }
  if( type==typeid(Eigen::Matrix4f) ) { return (Eigen::Matrix4f)Eigen::Matrix4f::Identity(); }
  if( type==typeid(Eigen::Matrix4i) ) { return (Eigen::Matrix4i)Eigen::Matrix4i::Identity(); }

  if( type==typeid(Eigen::MatrixXd) ) { return (Eigen::MatrixXd)Eigen::MatrixXd::Identity(rows, cols); }
  if( type==typeid(Eigen::MatrixXf) ) { return (Eigen::MatrixXf)Eigen::MatrixXf::Identity(rows, cols); }
  if( type==typeid(Eigen::MatrixXi) ) { return (Eigen::MatrixXi)Eigen::MatrixXi::Identity(rows, cols); }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any AnyAlgebra::EigenVectorIdentity(std::type_info const& type, unsigned int const rows, unsigned int const cols) const {
  if( type==typeid(Eigen::Vector2d) ) { return (Eigen::Matrix2d)Eigen::Matrix2d::Identity(); }
  if( type==typeid(Eigen::Vector2f) ) { return (Eigen::Matrix2f)Eigen::Matrix2f::Identity(); }
  if( type==typeid(Eigen::Vector2i) ) { return (Eigen::Matrix2i)Eigen::Matrix2i::Identity(); }
  
  if( type==typeid(Eigen::Vector3d) ) { return (Eigen::Matrix3d)Eigen::Matrix3d::Identity(); }
  if( type==typeid(Eigen::Vector3f) ) { return (Eigen::Matrix3f)Eigen::Matrix3f::Identity(); }
  if( type==typeid(Eigen::Vector3i) ) { return (Eigen::Matrix3i)Eigen::Matrix3i::Identity(); }
  
  if( type==typeid(Eigen::Vector4d) ) { return (Eigen::Matrix4d)Eigen::Matrix4d::Identity(); }
  if( type==typeid(Eigen::Vector4f) ) { return (Eigen::Matrix4f)Eigen::Matrix4f::Identity(); }
  if( type==typeid(Eigen::Vector4i) ) { return (Eigen::Matrix4i)Eigen::Matrix4i::Identity(); }

  if( type==typeid(Eigen::VectorXd) ) { return (Eigen::MatrixXd)Eigen::MatrixXd::Identity(rows, cols); }
  if( type==typeid(Eigen::VectorXf) ) { return (Eigen::MatrixXf)Eigen::MatrixXf::Identity(rows, cols); }
  if( type==typeid(Eigen::VectorXi) ) { return (Eigen::MatrixXi)Eigen::MatrixXi::Identity(rows, cols); }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any AnyAlgebra::Identity(std::type_info const& type, unsigned int const rows, unsigned int const cols) const {
  if( IsScalar(type) ) { return ScalarIdentity(type); }

  if( IsEigenVector(type) ) { return EigenVectorIdentity(type, rows, cols); }

  if( IsEigenMatrix(type) ) { return EigenMatrixIdentity(type, rows, cols); }
  
  return IdentityImpl(type, rows, cols);
}

boost::any AnyAlgebra::IdentityImpl(std::type_info const& type, unsigned int const rows, unsigned int const cols) const {
  std::cerr << std::endl << "ERROR: No way to compute identy object with type " << boost::core::demangle(type.name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::IdentityImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::IdentityImpl()" << std::endl << std::endl;
  assert(false);
  
  return boost::none;
}

boost::any AnyAlgebra::AddScalar(boost::any const& in0, boost::any const& in1) const {
  if( in0.type()==typeid(double) ) { return AddScalar<double>(in0, in1); }
  if( in0.type()==typeid(float) ) { return AddScalar<float>(in0, in1); }
  if( in0.type()==typeid(int) ) { return AddScalar<int>(in0, in1); }
  if( in0.type()==typeid(unsigned int) ) { return AddScalar<unsigned int>(in0, in1); }
  
  // something went wrong
  assert(false);
  return boost::none;
}

boost::any AnyAlgebra::AddEigenVector(boost::any const& in0, boost::any const& in1) const {
  // 2D vectors
  if( in0.type()==typeid(Eigen::Vector2d) ) {
    if( in1.type()==typeid(Eigen::Vector2d) ) { return AddEigenVector<Eigen::Vector2d, Eigen::Vector2d>(in0, in1); }
    return AddEigenVector<Eigen::Vector2d, Eigen::VectorXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector2f) ) {
    if( in1.type()==typeid(Eigen::Vector2f) ) { return AddEigenVector<Eigen::Vector2f, Eigen::Vector2f>(in0, in1); }
    return AddEigenVector<Eigen::Vector2f, Eigen::VectorXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector2i) ) {
    if( in1.type()==typeid(Eigen::Vector2i) ) { return AddEigenVector<Eigen::Vector2i, Eigen::Vector2i>(in0, in1); }
    return AddEigenVector<Eigen::Vector2i, Eigen::VectorXi>(in0, in1); 
  }

  // 3D vectors
  if( in0.type()==typeid(Eigen::Vector3d) ) {
    if( in1.type()==typeid(Eigen::Vector3d) ) { return AddEigenVector<Eigen::Vector3d, Eigen::Vector3d>(in0, in1); }
    return AddEigenVector<Eigen::Vector3d, Eigen::VectorXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector3f) ) {
    if( in1.type()==typeid(Eigen::Vector3f) ) { return AddEigenVector<Eigen::Vector3f, Eigen::Vector3f>(in0, in1); }
    return AddEigenVector<Eigen::Vector3f, Eigen::VectorXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector3i) ) {
    if( in1.type()==typeid(Eigen::Vector3i) ) { return AddEigenVector<Eigen::Vector3i, Eigen::Vector3i>(in0, in1); }
    return AddEigenVector<Eigen::Vector3i, Eigen::VectorXi>(in0, in1); 
  }

  // 4D vectors
  if( in0.type()==typeid(Eigen::Vector4d) ) {
    if( in1.type()==typeid(Eigen::Vector4d) ) { return AddEigenVector<Eigen::Vector4d, Eigen::Vector4d>(in0, in1); }
    return AddEigenVector<Eigen::Vector4d, Eigen::VectorXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector4f) ) {
    if( in1.type()==typeid(Eigen::Vector4f) ) { return AddEigenVector<Eigen::Vector4f, Eigen::Vector4f>(in0, in1); }
    return AddEigenVector<Eigen::Vector4f, Eigen::VectorXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector4i) ) {
    if( in1.type()==typeid(Eigen::Vector4i) ) { return AddEigenVector<Eigen::Vector4i, Eigen::Vector4i>(in0, in1); }
    return AddEigenVector<Eigen::Vector4i, Eigen::VectorXi>(in0, in1); 
  }

  // XD vectors
  if( in0.type()==typeid(Eigen::VectorXd) ) {
    if( in1.type()==typeid(Eigen::Vector2d) ) { return AddEigenVector<Eigen::VectorXd, Eigen::Vector2d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector3d) ) { return AddEigenVector<Eigen::VectorXd, Eigen::Vector3d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector4d) ) { return AddEigenVector<Eigen::VectorXd, Eigen::Vector4d>(in0, in1); }
    return AddEigenVector<Eigen::VectorXd, Eigen::VectorXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::VectorXf) ) {
    if( in1.type()==typeid(Eigen::Vector2f) ) { return AddEigenVector<Eigen::VectorXf, Eigen::Vector2f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector3f) ) { return AddEigenVector<Eigen::VectorXf, Eigen::Vector3f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector4f) ) { return AddEigenVector<Eigen::VectorXf, Eigen::Vector4f>(in0, in1); }
    return AddEigenVector<Eigen::VectorXf, Eigen::VectorXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::VectorXi) ) {
    if( in1.type()==typeid(Eigen::Vector2i) ) { return AddEigenVector<Eigen::VectorXi, Eigen::Vector2i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector3i) ) { return AddEigenVector<Eigen::VectorXi, Eigen::Vector3i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector4i) ) { return AddEigenVector<Eigen::VectorXi, Eigen::Vector4i>(in0, in1); }
    return AddEigenVector<Eigen::VectorXi, Eigen::VectorXi>(in0, in1); 
  }
  
  // something went wrong
  assert(false);
  return boost::none;
}

boost::any AnyAlgebra::AddEigenMatrix(boost::any const& in0, boost::any const& in1) const {
  // 2D matrices
  if( in0.type()==typeid(Eigen::Matrix2d) ) {
    if( in1.type()==typeid(Eigen::Matrix2d) ) { return AddEigenMatrix<Eigen::Matrix2d, Eigen::Matrix2d>(in0, in1); }
    return AddEigenMatrix<Eigen::Matrix2d, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix2f) ) {
    if( in1.type()==typeid(Eigen::Matrix2f) ) { return AddEigenMatrix<Eigen::Matrix2f, Eigen::Matrix2f>(in0, in1); }
    return AddEigenMatrix<Eigen::Matrix2f, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix2i) ) {
    if( in1.type()==typeid(Eigen::Matrix2i) ) { return AddEigenMatrix<Eigen::Matrix2i, Eigen::Matrix2i>(in0, in1); }
    return AddEigenMatrix<Eigen::Matrix2i, Eigen::MatrixXi>(in0, in1); 
  }

  // 3D matrices
  if( in0.type()==typeid(Eigen::Matrix3d) ) {
    if( in1.type()==typeid(Eigen::Matrix3d) ) { return AddEigenMatrix<Eigen::Matrix3d, Eigen::Matrix3d>(in0, in1); }
    return AddEigenMatrix<Eigen::Matrix3d, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix3f) ) {
    if( in1.type()==typeid(Eigen::Matrix3f) ) { return AddEigenMatrix<Eigen::Matrix3f, Eigen::Matrix3f>(in0, in1); }
    return AddEigenMatrix<Eigen::Matrix3f, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix3i) ) {
    if( in1.type()==typeid(Eigen::Matrix3i) ) { return AddEigenMatrix<Eigen::Matrix3i, Eigen::Matrix3i>(in0, in1); }
    return AddEigenMatrix<Eigen::Matrix3i, Eigen::MatrixXi>(in0, in1); 
  }

  // 4D matrices
  if( in0.type()==typeid(Eigen::Matrix4d) ) {
    if( in1.type()==typeid(Eigen::Matrix4d) ) { return AddEigenMatrix<Eigen::Matrix4d, Eigen::Matrix4d>(in0, in1); }
    return AddEigenMatrix<Eigen::Matrix4d, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix4f) ) {
    if( in1.type()==typeid(Eigen::Matrix4f) ) { return AddEigenMatrix<Eigen::Matrix4f, Eigen::Matrix4f>(in0, in1); }
    return AddEigenMatrix<Eigen::Matrix4f, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix4i) ) {
    if( in1.type()==typeid(Eigen::Matrix4i) ) { return AddEigenMatrix<Eigen::Matrix4i, Eigen::Matrix4i>(in0, in1); }
    return AddEigenMatrix<Eigen::Matrix4i, Eigen::MatrixXi>(in0, in1); 
  }

  // XD matrices
  if( in0.type()==typeid(Eigen::MatrixXd) ) {
    if( in1.type()==typeid(Eigen::Matrix2d) ) { return AddEigenMatrix<Eigen::MatrixXd, Eigen::Matrix2d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3d) ) { return AddEigenMatrix<Eigen::MatrixXd, Eigen::Matrix3d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4d) ) { return AddEigenMatrix<Eigen::MatrixXd, Eigen::Matrix4d>(in0, in1); }
    return AddEigenMatrix<Eigen::MatrixXd, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::MatrixXf) ) {
    if( in1.type()==typeid(Eigen::Matrix2f) ) { return AddEigenMatrix<Eigen::MatrixXf, Eigen::Matrix2f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3f) ) { return AddEigenMatrix<Eigen::MatrixXf, Eigen::Matrix3f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4f) ) { return AddEigenMatrix<Eigen::MatrixXf, Eigen::Matrix4f>(in0, in1); }
    return AddEigenMatrix<Eigen::MatrixXf, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::MatrixXi) ) {
    if( in1.type()==typeid(Eigen::Matrix2i) ) { return AddEigenMatrix<Eigen::MatrixXi, Eigen::Matrix2i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3i) ) { return AddEigenMatrix<Eigen::MatrixXi, Eigen::Matrix3i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4i) ) { return AddEigenMatrix<Eigen::MatrixXi, Eigen::Matrix4i>(in0, in1); }
    return AddEigenMatrix<Eigen::MatrixXi, Eigen::MatrixXi>(in0, in1); 
  }
  
  // something went wrong
  assert(false);
  return boost::none;
}

boost::any AnyAlgebra::Add(boost::any const& in0, boost::any const& in1) const {
  if( IsScalar(in0.type()) && IsScalar(in1.type()) ) { return AddScalar(in0, in1); }

  if( IsEigenVector(in0.type()) && IsEigenVector(in1.type()) ) { return AddEigenVector(in0, in1); }

  if( IsEigenMatrix(in0.type()) && IsEigenMatrix(in1.type()) ) { return AddEigenMatrix(in0, in1); }
  
  // the first type is boost::none --- return the second
  if( in0.type()==typeid(boost::none) ) { return in1; }

  // the second type is boost::none --- return the first
  if( in1.type()==typeid(boost::none) ) { return in0; }

  return AddImpl(in0, in1);
}

boost::any AnyAlgebra::AddImpl(boost::any const& in0, boost::any const& in1) const {
  std::cerr << std::endl << "ERROR: No way to add type " << boost::core::demangle(in0.type().name()) << " and type " << boost::core::demangle(in1.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::AddImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::AddImpl()" << std::endl << std::endl;
  assert(false);
  
  return boost::none;
}

boost::any AnyAlgebra::SubtractScalar(boost::any const& in0, boost::any const& in1) const {
  if( in0.type()==typeid(double) ) { return SubtractScalar<double>(in0, in1); }
  if( in0.type()==typeid(float) ) { return SubtractScalar<float>(in0, in1); }
  if( in0.type()==typeid(int) ) { return SubtractScalar<int>(in0, in1); }
  if( in0.type()==typeid(unsigned int) ) { return SubtractScalar<unsigned int>(in0, in1); }
  
  // something went wrong
  assert(false);
  return boost::none;
}

boost::any AnyAlgebra::SubtractEigenVector(boost::any const& in0, boost::any const& in1) const {
  // 2D vectors
  if( in0.type()==typeid(Eigen::Vector2d) ) {
    if( in1.type()==typeid(Eigen::Vector2d) ) { return SubtractEigenVector<Eigen::Vector2d, Eigen::Vector2d>(in0, in1); }
    return SubtractEigenVector<Eigen::Vector2d, Eigen::VectorXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector2f) ) {
    if( in1.type()==typeid(Eigen::Vector2f) ) { return SubtractEigenVector<Eigen::Vector2f, Eigen::Vector2f>(in0, in1); }
    return SubtractEigenVector<Eigen::Vector2f, Eigen::VectorXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector2i) ) {
    if( in1.type()==typeid(Eigen::Vector2i) ) { return SubtractEigenVector<Eigen::Vector2i, Eigen::Vector2i>(in0, in1); }
    return SubtractEigenVector<Eigen::Vector2i, Eigen::VectorXi>(in0, in1); 
  }

  // 3D vectors
  if( in0.type()==typeid(Eigen::Vector3d) ) {
    if( in1.type()==typeid(Eigen::Vector3d) ) { return SubtractEigenVector<Eigen::Vector3d, Eigen::Vector3d>(in0, in1); }
    return SubtractEigenVector<Eigen::Vector3d, Eigen::VectorXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector3f) ) {
    if( in1.type()==typeid(Eigen::Vector3f) ) { return SubtractEigenVector<Eigen::Vector3f, Eigen::Vector3f>(in0, in1); }
    return SubtractEigenVector<Eigen::Vector3f, Eigen::VectorXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector3i) ) {
    if( in1.type()==typeid(Eigen::Vector3i) ) { return SubtractEigenVector<Eigen::Vector3i, Eigen::Vector3i>(in0, in1); }
    return SubtractEigenVector<Eigen::Vector3i, Eigen::VectorXi>(in0, in1); 
  }

  // 4D vectors
  if( in0.type()==typeid(Eigen::Vector4d) ) {
    if( in1.type()==typeid(Eigen::Vector4d) ) { return SubtractEigenVector<Eigen::Vector4d, Eigen::Vector4d>(in0, in1); }
    return SubtractEigenVector<Eigen::Vector4d, Eigen::VectorXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector4f) ) {
    if( in1.type()==typeid(Eigen::Vector4f) ) { return SubtractEigenVector<Eigen::Vector4f, Eigen::Vector4f>(in0, in1); }
    return SubtractEigenVector<Eigen::Vector4f, Eigen::VectorXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Vector4i) ) {
    if( in1.type()==typeid(Eigen::Vector4i) ) { return SubtractEigenVector<Eigen::Vector4i, Eigen::Vector4i>(in0, in1); }
    return SubtractEigenVector<Eigen::Vector4i, Eigen::VectorXi>(in0, in1); 
  }

  // XD vectors
  if( in0.type()==typeid(Eigen::VectorXd) ) {
    if( in1.type()==typeid(Eigen::Vector2d) ) { return SubtractEigenVector<Eigen::VectorXd, Eigen::Vector2d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector3d) ) { return SubtractEigenVector<Eigen::VectorXd, Eigen::Vector3d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector4d) ) { return SubtractEigenVector<Eigen::VectorXd, Eigen::Vector4d>(in0, in1); }
    return SubtractEigenVector<Eigen::VectorXd, Eigen::VectorXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::VectorXf) ) {
    if( in1.type()==typeid(Eigen::Vector2f) ) { return SubtractEigenVector<Eigen::VectorXf, Eigen::Vector2f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector3f) ) { return SubtractEigenVector<Eigen::VectorXf, Eigen::Vector3f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector4f) ) { return SubtractEigenVector<Eigen::VectorXf, Eigen::Vector4f>(in0, in1); }
    return SubtractEigenVector<Eigen::VectorXf, Eigen::VectorXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::VectorXi) ) {
    if( in1.type()==typeid(Eigen::Vector2i) ) { return SubtractEigenVector<Eigen::VectorXi, Eigen::Vector2i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector3i) ) { return SubtractEigenVector<Eigen::VectorXi, Eigen::Vector3i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector4i) ) { return SubtractEigenVector<Eigen::VectorXi, Eigen::Vector4i>(in0, in1); }
    return SubtractEigenVector<Eigen::VectorXi, Eigen::VectorXi>(in0, in1); 
  }
  
  // something went wrong
  assert(false);
  return boost::none;
}

boost::any AnyAlgebra::SubtractEigenMatrix(boost::any const& in0, boost::any const& in1) const {
  // 2D matrices
  if( in0.type()==typeid(Eigen::Matrix2d) ) {
    if( in1.type()==typeid(Eigen::Matrix2d) ) { return SubtractEigenMatrix<Eigen::Matrix2d, Eigen::Matrix2d>(in0, in1); }
    return SubtractEigenMatrix<Eigen::Matrix2d, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix2f) ) {
    if( in1.type()==typeid(Eigen::Matrix2f) ) { return SubtractEigenMatrix<Eigen::Matrix2f, Eigen::Matrix2f>(in0, in1); }
    return SubtractEigenMatrix<Eigen::Matrix2f, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix2i) ) {
    if( in1.type()==typeid(Eigen::Matrix2i) ) { return SubtractEigenMatrix<Eigen::Matrix2i, Eigen::Matrix2i>(in0, in1); }
    return SubtractEigenMatrix<Eigen::Matrix2i, Eigen::MatrixXi>(in0, in1); 
  }

  // 3D matrices
  if( in0.type()==typeid(Eigen::Matrix3d) ) {
    if( in1.type()==typeid(Eigen::Matrix3d) ) { return SubtractEigenMatrix<Eigen::Matrix3d, Eigen::Matrix3d>(in0, in1); }
    return SubtractEigenMatrix<Eigen::Matrix3d, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix3f) ) {
    if( in1.type()==typeid(Eigen::Matrix3f) ) { return SubtractEigenMatrix<Eigen::Matrix3f, Eigen::Matrix3f>(in0, in1); }
    return SubtractEigenMatrix<Eigen::Matrix3f, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix3i) ) {
    if( in1.type()==typeid(Eigen::Matrix3i) ) { return SubtractEigenMatrix<Eigen::Matrix3i, Eigen::Matrix3i>(in0, in1); }
    return SubtractEigenMatrix<Eigen::Matrix3i, Eigen::MatrixXi>(in0, in1); 
  }

  // 4D matrices
  if( in0.type()==typeid(Eigen::Matrix4d) ) {
    if( in1.type()==typeid(Eigen::Matrix4d) ) { return SubtractEigenMatrix<Eigen::Matrix4d, Eigen::Matrix4d>(in0, in1); }
    return SubtractEigenMatrix<Eigen::Matrix4d, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix4f) ) {
    if( in1.type()==typeid(Eigen::Matrix4f) ) { return SubtractEigenMatrix<Eigen::Matrix4f, Eigen::Matrix4f>(in0, in1); }
    return SubtractEigenMatrix<Eigen::Matrix4f, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix4i) ) {
    if( in1.type()==typeid(Eigen::Matrix4i) ) { return SubtractEigenMatrix<Eigen::Matrix4i, Eigen::Matrix4i>(in0, in1); }
    return SubtractEigenMatrix<Eigen::Matrix4i, Eigen::MatrixXi>(in0, in1); 
  }

  // XD matrices
  if( in0.type()==typeid(Eigen::MatrixXd) ) {
    if( in1.type()==typeid(Eigen::Matrix2d) ) { return SubtractEigenMatrix<Eigen::MatrixXd, Eigen::Matrix2d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3d) ) { return SubtractEigenMatrix<Eigen::MatrixXd, Eigen::Matrix3d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4d) ) { return SubtractEigenMatrix<Eigen::MatrixXd, Eigen::Matrix4d>(in0, in1); }
    return SubtractEigenMatrix<Eigen::MatrixXd, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::MatrixXf) ) {
    if( in1.type()==typeid(Eigen::Matrix2f) ) { return SubtractEigenMatrix<Eigen::MatrixXf, Eigen::Matrix2f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3f) ) { return SubtractEigenMatrix<Eigen::MatrixXf, Eigen::Matrix3f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4f) ) { return SubtractEigenMatrix<Eigen::MatrixXf, Eigen::Matrix4f>(in0, in1); }
    return SubtractEigenMatrix<Eigen::MatrixXf, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::MatrixXi) ) {
    if( in1.type()==typeid(Eigen::Matrix2i) ) { return SubtractEigenMatrix<Eigen::MatrixXi, Eigen::Matrix2i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3i) ) { return SubtractEigenMatrix<Eigen::MatrixXi, Eigen::Matrix3i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4i) ) { return SubtractEigenMatrix<Eigen::MatrixXi, Eigen::Matrix4i>(in0, in1); }
    return SubtractEigenMatrix<Eigen::MatrixXi, Eigen::MatrixXi>(in0, in1); 
  }
  
  // something went wrong
  assert(false);
  return boost::none;
}

boost::any AnyAlgebra::Subtract(boost::any const& in0, boost::any const& in1) const {
  if( IsScalar(in0.type()) && IsScalar(in1.type()) ) { return SubtractScalar(in0, in1); }

  if( IsEigenVector(in0.type()) && IsEigenVector(in1.type()) ) { return SubtractEigenVector(in0, in1); }

  if( IsEigenMatrix(in0.type()) && IsEigenMatrix(in1.type()) ) { return SubtractEigenMatrix(in0, in1); }
  
  // the first type is boost::none --- return the second
  if( in0.type()==typeid(boost::none) ) { return in1; }

  // the second type is boost::none --- return the first
  if( in1.type()==typeid(boost::none) ) { return in0; }

  return SubtractImpl(in0, in1);
}

boost::any AnyAlgebra::SubtractImpl(boost::any const& in0, boost::any const& in1) const {
  std::cerr << std::endl << "ERROR: No way to subtract type " << boost::core::demangle(in0.type().name()) << " and type " << boost::core::demangle(in1.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::SubtractImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::SubtractImpl()" << std::endl << std::endl;
  assert(false);
  
  return boost::none;
}

boost::any AnyAlgebra::MultiplyScalar(boost::any const& in0, boost::any const& in1) const {
  if( in0.type()==typeid(double) ) { return MultiplyScalar<double>(in0, in1); }
  if( in0.type()==typeid(float) ) { return MultiplyScalar<float>(in0, in1); }
  if( in0.type()==typeid(int) ) { return MultiplyScalar<int>(in0, in1); }
  if( in0.type()==typeid(unsigned int) ) { return MultiplyScalar<unsigned int>(in0, in1); }
  
  // something went wrong
  assert(false);
  return boost::none;
}

boost::any AnyAlgebra::MultiplyEigenMatrix(boost::any const& in0, boost::any const& in1) const {
  // 2D matrices
  if( in0.type()==typeid(Eigen::Matrix2d) ) {
    if( in1.type()==typeid(Eigen::Matrix2d) ) { return MultiplyEigenMatrix<Eigen::Matrix2d, Eigen::Matrix2d>(in0, in1); }
    return MultiplyEigenMatrix<Eigen::Matrix2d, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix2f) ) {
    if( in1.type()==typeid(Eigen::Matrix2f) ) { return MultiplyEigenMatrix<Eigen::Matrix2f, Eigen::Matrix2f>(in0, in1); }
    return MultiplyEigenMatrix<Eigen::Matrix2f, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix2i) ) {
    if( in1.type()==typeid(Eigen::Matrix2i) ) { return MultiplyEigenMatrix<Eigen::Matrix2i, Eigen::Matrix2i>(in0, in1); }
    return MultiplyEigenMatrix<Eigen::Matrix2i, Eigen::MatrixXi>(in0, in1); 
  }

  // 3D matrices
  if( in0.type()==typeid(Eigen::Matrix3d) ) {
    if( in1.type()==typeid(Eigen::Matrix3d) ) { return MultiplyEigenMatrix<Eigen::Matrix3d, Eigen::Matrix3d>(in0, in1); }
    return MultiplyEigenMatrix<Eigen::Matrix3d, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix3f) ) {
    if( in1.type()==typeid(Eigen::Matrix3f) ) { return MultiplyEigenMatrix<Eigen::Matrix3f, Eigen::Matrix3f>(in0, in1); }
    return MultiplyEigenMatrix<Eigen::Matrix3f, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix3i) ) {
    if( in1.type()==typeid(Eigen::Matrix3i) ) { return MultiplyEigenMatrix<Eigen::Matrix3i, Eigen::Matrix3i>(in0, in1); }
    return MultiplyEigenMatrix<Eigen::Matrix3i, Eigen::MatrixXi>(in0, in1); 
  }

  // 4D matrices
  if( in0.type()==typeid(Eigen::Matrix4d) ) {
    if( in1.type()==typeid(Eigen::Matrix4d) ) { return MultiplyEigenMatrix<Eigen::Matrix4d, Eigen::Matrix4d>(in0, in1); }
    return MultiplyEigenMatrix<Eigen::Matrix4d, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix4f) ) {
    if( in1.type()==typeid(Eigen::Matrix4f) ) { return MultiplyEigenMatrix<Eigen::Matrix4f, Eigen::Matrix4f>(in0, in1); }
    return MultiplyEigenMatrix<Eigen::Matrix4f, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::Matrix4i) ) {
    if( in1.type()==typeid(Eigen::Matrix4i) ) { return MultiplyEigenMatrix<Eigen::Matrix4i, Eigen::Matrix4i>(in0, in1); }
    return MultiplyEigenMatrix<Eigen::Matrix4i, Eigen::MatrixXi>(in0, in1); 
  }

  // XD matrices
  if( in0.type()==typeid(Eigen::MatrixXd) ) {
    if( in1.type()==typeid(Eigen::Matrix2d) ) { return MultiplyEigenMatrix<Eigen::MatrixXd, Eigen::Matrix2d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3d) ) { return MultiplyEigenMatrix<Eigen::MatrixXd, Eigen::Matrix3d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4d) ) { return MultiplyEigenMatrix<Eigen::MatrixXd, Eigen::Matrix4d>(in0, in1); }
    return MultiplyEigenMatrix<Eigen::MatrixXd, Eigen::MatrixXd>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::MatrixXf) ) {
    if( in1.type()==typeid(Eigen::Matrix2f) ) { return MultiplyEigenMatrix<Eigen::MatrixXf, Eigen::Matrix2f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3f) ) { return MultiplyEigenMatrix<Eigen::MatrixXf, Eigen::Matrix3f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4f) ) { return MultiplyEigenMatrix<Eigen::MatrixXf, Eigen::Matrix4f>(in0, in1); }
    return MultiplyEigenMatrix<Eigen::MatrixXf, Eigen::MatrixXf>(in0, in1); 
  }
  if( in0.type()==typeid(Eigen::MatrixXi) ) {
    if( in1.type()==typeid(Eigen::Matrix2i) ) { return MultiplyEigenMatrix<Eigen::MatrixXi, Eigen::Matrix2i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3i) ) { return MultiplyEigenMatrix<Eigen::MatrixXi, Eigen::Matrix3i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4i) ) { return MultiplyEigenMatrix<Eigen::MatrixXi, Eigen::Matrix4i>(in0, in1); }
    return MultiplyEigenMatrix<Eigen::MatrixXi, Eigen::MatrixXi>(in0, in1); 
  }
  
  // something went wrong
  assert(false);
  return boost::none;
}

boost::any AnyAlgebra::MultiplyEigenVectorScalar(boost::any const& in0, boost::any const& in1) const {
  if( in0.type()==typeid(double) ) {
    if( in1.type()==typeid(Eigen::Vector2d) ) { return MultiplyEigenScalar<double, Eigen::Vector2d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector3d) ) { return MultiplyEigenScalar<double, Eigen::Vector3d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector4d) ) { return MultiplyEigenScalar<double, Eigen::Vector4d>(in0, in1); }
    
    return MultiplyEigenScalar<double, Eigen::VectorXd>(in0, in1); 
  }

  if( in0.type()==typeid(float) ) {
    if( in1.type()==typeid(Eigen::Vector2f) ) { return MultiplyEigenScalar<float, Eigen::Vector2f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector3f) ) { return MultiplyEigenScalar<float, Eigen::Vector3f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector4f) ) { return MultiplyEigenScalar<float, Eigen::Vector4f>(in0, in1); }

    return MultiplyEigenScalar<float, Eigen::VectorXf>(in0, in1); 
  }

  if( in0.type()==typeid(int) ) {
    if( in1.type()==typeid(Eigen::Vector2i) ) { return MultiplyEigenScalar<int, Eigen::Vector2i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector3i) ) { return MultiplyEigenScalar<int, Eigen::Vector3i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector4i) ) { return MultiplyEigenScalar<int, Eigen::Vector4i>(in0, in1); }

    return MultiplyEigenScalar<int, Eigen::VectorXi>(in0, in1); 
  }

  if( in0.type()==typeid(unsigned int) ) {
    if( in1.type()==typeid(Eigen::Vector2i) ) { return MultiplyEigenScalar<unsigned int, Eigen::Vector2i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector3i) ) { return MultiplyEigenScalar<unsigned int, Eigen::Vector3i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Vector4i) ) { return MultiplyEigenScalar<unsigned int, Eigen::Vector4i>(in0, in1); }

    return MultiplyEigenScalar<unsigned int, Eigen::VectorXi>(in0, in1); 
  }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any AnyAlgebra::MultiplyEigenMatrixScalar(boost::any const& in0, boost::any const& in1) const {
  if( in0.type()==typeid(double) ) {
    if( in1.type()==typeid(Eigen::Matrix2d) ) { return MultiplyEigenScalar<double, Eigen::Matrix2d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3d) ) { return MultiplyEigenScalar<double, Eigen::Matrix3d>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4d) ) { return MultiplyEigenScalar<double, Eigen::Matrix4d>(in0, in1); }

    return MultiplyEigenScalar<double, Eigen::MatrixXd>(in0, in1); 
  }

  if( in0.type()==typeid(float) ) {
    if( in1.type()==typeid(Eigen::Matrix2f) ) { return MultiplyEigenScalar<float, Eigen::Matrix2f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3f) ) { return MultiplyEigenScalar<float, Eigen::Matrix3f>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4f) ) { return MultiplyEigenScalar<float, Eigen::Matrix4f>(in0, in1); }

    return MultiplyEigenScalar<float, Eigen::MatrixXf>(in0, in1); 
  }

  if( in0.type()==typeid(int) ) {
    if( in1.type()==typeid(Eigen::Matrix2i) ) { return MultiplyEigenScalar<int, Eigen::Matrix2i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3i) ) { return MultiplyEigenScalar<int, Eigen::Matrix3i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4i) ) { return MultiplyEigenScalar<int, Eigen::Matrix4i>(in0, in1); }

    return MultiplyEigenScalar<int, Eigen::MatrixXi>(in0, in1); 
  }

  if( in0.type()==typeid(unsigned int) ) {
    if( in1.type()==typeid(Eigen::Matrix2i) ) { return MultiplyEigenScalar<unsigned int, Eigen::Matrix2i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix3i) ) { return MultiplyEigenScalar<unsigned int, Eigen::Matrix3i>(in0, in1); }
    if( in1.type()==typeid(Eigen::Matrix4i) ) { return MultiplyEigenScalar<unsigned int, Eigen::Matrix4i>(in0, in1); }

    return MultiplyEigenScalar<unsigned int, Eigen::MatrixXi>(in0, in1); 
  }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any AnyAlgebra::Multiply(boost::any const& in0, boost::any const& in1) const {
  if( IsScalar(in0.type()) && IsScalar(in1.type()) ) { return MultiplyScalar(in0, in1); }

  if( IsScalar(in0.type()) && IsEigenVector(in1.type()) ) { return MultiplyEigenVectorScalar(in0, in1); }
  if( IsScalar(in1.type()) && IsEigenVector(in0.type()) ) { return MultiplyEigenVectorScalar(in1, in0); }

  if( IsEigenMatrix(in0.type()) && IsEigenMatrix(in1.type()) ) { return MultiplyEigenMatrix(in0, in1); }

  if( IsScalar(in0.type()) && IsEigenMatrix(in1.type()) ) { return MultiplyEigenMatrixScalar(in0, in1); }
  if( IsScalar(in1.type()) && IsEigenMatrix(in0.type()) ) { return MultiplyEigenMatrixScalar(in1, in0); }

  // the first type is boost::none --- return the second
  if( in0.type()==typeid(boost::none) ) { return in1; }

  // the second type is boost::none --- return the first
  if( in1.type()==typeid(boost::none) ) { return in0; }

  return MultiplyImpl(in0, in1);
}

boost::any AnyAlgebra::MultiplyImpl(boost::any const& in0, boost::any const& in1) const {
  std::cerr << std::endl << "ERROR: No way to multiply type " << boost::core::demangle(in0.type().name()) << " and type " << boost::core::demangle(in1.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::MultiplyImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::MultiplyImpl()" << std::endl << std::endl;
  assert(false);
  
  return boost::none;
}

boost::any AnyAlgebra::ApplyScalarInverse(boost::any const& A, boost::any const& x) const {
  if( A.type()==typeid(double) ) { return Multiply(1.0/boost::any_cast<double const>(A), x); }
  if( A.type()==typeid(float) ) { return Multiply((float)1.0/boost::any_cast<float const>(A), x); }
  if( A.type()==typeid(int) ) { return Multiply(1.0/boost::any_cast<int const>(A), x); }
  if( A.type()==typeid(unsigned int) ) { return Multiply(1.0/boost::any_cast<unsigned int const>(A), x); }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any AnyAlgebra::ApplyEigenVectorInverse(boost::any const& A, boost::any const& x) const {
  if( A.type()==typeid(Eigen::Vector2d) ) {
    if( x.type()==typeid(Eigen::Vector2d) ) { return ApplyEigenVectorInverse<Eigen::Vector2d, Eigen::Vector2d>(A, x); }
    return ApplyEigenVectorInverse<Eigen::Vector2d, Eigen::VectorXd>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector2f) ) {
    if( x.type()==typeid(Eigen::Vector2f) ) { return ApplyEigenVectorInverse<Eigen::Vector2f, Eigen::Vector2f>(A, x); }
    return ApplyEigenVectorInverse<Eigen::Vector2f, Eigen::VectorXf>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector2i) ) {
    if( x.type()==typeid(Eigen::Vector2i) ) { return ApplyEigenVectorInverse<Eigen::Vector2i, Eigen::Vector2i>(A, x); }
    return ApplyEigenVectorInverse<Eigen::Vector2i, Eigen::VectorXi>(A, x);
  }

  if( A.type()==typeid(Eigen::Vector3d) ) {
    if( x.type()==typeid(Eigen::Vector3d) ) { return ApplyEigenVectorInverse<Eigen::Vector3d, Eigen::Vector3d>(A, x); }
    return ApplyEigenVectorInverse<Eigen::Vector3d, Eigen::VectorXd>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector3f) ) {
    if( x.type()==typeid(Eigen::Vector3f) ) { return ApplyEigenVectorInverse<Eigen::Vector3f, Eigen::Vector3f>(A, x); }
    return ApplyEigenVectorInverse<Eigen::Vector3f, Eigen::VectorXf>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector3i) ) {
    if( x.type()==typeid(Eigen::Vector3i) ) { return ApplyEigenVectorInverse<Eigen::Vector3i, Eigen::Vector3i>(A, x); }
    return ApplyEigenVectorInverse<Eigen::Vector3i, Eigen::VectorXi>(A, x);
  }

  if( A.type()==typeid(Eigen::Vector4d) ) {
    if( x.type()==typeid(Eigen::Vector4d) ) { return ApplyEigenVectorInverse<Eigen::Vector4d, Eigen::Vector4d>(A, x); }
    return ApplyEigenVectorInverse<Eigen::Vector4d, Eigen::VectorXd>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector4f) ) {
    if( x.type()==typeid(Eigen::Vector4f) ) { return ApplyEigenVectorInverse<Eigen::Vector4f, Eigen::Vector4f>(A, x); }
    return ApplyEigenVectorInverse<Eigen::Vector4f, Eigen::VectorXf>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector4i) ) {
    if( x.type()==typeid(Eigen::Vector4i) ) { return ApplyEigenVectorInverse<Eigen::Vector4i, Eigen::Vector4i>(A, x); }
    return ApplyEigenVectorInverse<Eigen::Vector4i, Eigen::VectorXi>(A, x);
  }

  if( A.type()==typeid(Eigen::VectorXd) ) {
    if( x.type()==typeid(Eigen::Vector2d) ) { return ApplyEigenVectorInverse<Eigen::VectorXd, Eigen::Vector2d>(A, x); }
    if( x.type()==typeid(Eigen::Vector3d) ) { return ApplyEigenVectorInverse<Eigen::VectorXd, Eigen::Vector3d>(A, x); }
    if( x.type()==typeid(Eigen::Vector4d) ) { return ApplyEigenVectorInverse<Eigen::VectorXd, Eigen::Vector4d>(A, x); }

    return ApplyEigenVectorInverse<Eigen::VectorXd, Eigen::VectorXd>(A, x);
  }
  if( A.type()==typeid(Eigen::VectorXf) ) {
    if( x.type()==typeid(Eigen::Vector2f) ) { return ApplyEigenVectorInverse<Eigen::VectorXf, Eigen::Vector2f>(A, x); }
    if( x.type()==typeid(Eigen::Vector3f) ) { return ApplyEigenVectorInverse<Eigen::VectorXf, Eigen::Vector3f>(A, x); }
    if( x.type()==typeid(Eigen::Vector4f) ) { return ApplyEigenVectorInverse<Eigen::VectorXf, Eigen::Vector4f>(A, x); }

    return ApplyEigenVectorInverse<Eigen::VectorXf, Eigen::VectorXf>(A, x);
  }
  if( A.type()==typeid(Eigen::VectorXi) ) {
    if( x.type()==typeid(Eigen::Vector2i) ) { return ApplyEigenVectorInverse<Eigen::VectorXi, Eigen::Vector2i>(A, x); }
    if( x.type()==typeid(Eigen::Vector3i) ) { return ApplyEigenVectorInverse<Eigen::VectorXi, Eigen::Vector3i>(A, x); }
    if( x.type()==typeid(Eigen::Vector4i) ) { return ApplyEigenVectorInverse<Eigen::VectorXi, Eigen::Vector4i>(A, x); }

    return ApplyEigenVectorInverse<Eigen::VectorXi, Eigen::VectorXi>(A, x);
  }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any AnyAlgebra::ApplyInverse(boost::any const& A, boost::any const& x) const {
  if( IsScalar(A.type()) ) { return ApplyScalarInverse(A, x); }

  if( IsEigenVector(A.type()) ) { return ApplyEigenVectorInverse(A, x); }
  
  return ApplyInverseImpl(A, x);
}

boost::any AnyAlgebra::ApplyInverseImpl(boost::any const& A, boost::any const& x) const {
  std::cerr << std::endl << "ERROR: No way to apply the inverse " << boost::core::demangle(A.type().name()) << " type to type " << boost::core::demangle(x.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::ApplyInverseImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::ApplyInverseImpl()" << std::endl << std::endl;
  assert(false);
  
  return boost::none;
}

boost::any AnyAlgebra::ApplyScalar(boost::any const& A, boost::any const& x) const {
  if( A.type()==typeid(double) ) { return Multiply(boost::any_cast<double const>(A), x); }
  if( A.type()==typeid(float) ) { return Multiply(boost::any_cast<float const>(A), x); }
  if( A.type()==typeid(int) ) { return Multiply(boost::any_cast<int const>(A), x); }
  if( A.type()==typeid(unsigned int) ) { return Multiply(boost::any_cast<unsigned int const>(A), x); }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any AnyAlgebra::ApplyEigenVector(boost::any const& A, boost::any const& x) const {
  if( A.type()==typeid(Eigen::Vector2d) ) {
    if( x.type()==typeid(Eigen::Vector2d) ) { return ApplyEigenVector<Eigen::Vector2d, Eigen::Vector2d>(A, x); }
    return ApplyEigenVector<Eigen::Vector2d, Eigen::VectorXd>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector2f) ) {
    if( x.type()==typeid(Eigen::Vector2f) ) { return ApplyEigenVector<Eigen::Vector2f, Eigen::Vector2f>(A, x); }
    return ApplyEigenVector<Eigen::Vector2f, Eigen::VectorXf>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector2i) ) {
    if( x.type()==typeid(Eigen::Vector2i) ) { return ApplyEigenVector<Eigen::Vector2i, Eigen::Vector2i>(A, x); }
    return ApplyEigenVector<Eigen::Vector2i, Eigen::VectorXi>(A, x);
  }

  if( A.type()==typeid(Eigen::Vector3d) ) {
    if( x.type()==typeid(Eigen::Vector3d) ) { return ApplyEigenVector<Eigen::Vector3d, Eigen::Vector3d>(A, x); }
    return ApplyEigenVector<Eigen::Vector3d, Eigen::VectorXd>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector3f) ) {
    if( x.type()==typeid(Eigen::Vector3f) ) { return ApplyEigenVector<Eigen::Vector3f, Eigen::Vector3f>(A, x); }
    return ApplyEigenVector<Eigen::Vector3f, Eigen::VectorXf>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector3i) ) {
    if( x.type()==typeid(Eigen::Vector3i) ) { return ApplyEigenVector<Eigen::Vector3i, Eigen::Vector3i>(A, x); }
    return ApplyEigenVector<Eigen::Vector3i, Eigen::VectorXi>(A, x);
  }

  if( A.type()==typeid(Eigen::Vector4d) ) {
    if( x.type()==typeid(Eigen::Vector4d) ) { return ApplyEigenVector<Eigen::Vector4d, Eigen::Vector4d>(A, x); }
    return ApplyEigenVector<Eigen::Vector4d, Eigen::VectorXd>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector4f) ) {
    if( x.type()==typeid(Eigen::Vector4f) ) { return ApplyEigenVector<Eigen::Vector4f, Eigen::Vector4f>(A, x); }
    return ApplyEigenVector<Eigen::Vector4f, Eigen::VectorXf>(A, x);
  }
  if( A.type()==typeid(Eigen::Vector4i) ) {
    if( x.type()==typeid(Eigen::Vector4i) ) { return ApplyEigenVector<Eigen::Vector4i, Eigen::Vector4i>(A, x); }
    return ApplyEigenVector<Eigen::Vector4i, Eigen::VectorXi>(A, x);
  }

  if( A.type()==typeid(Eigen::VectorXd) ) {
    if( x.type()==typeid(Eigen::Vector2d) ) { return ApplyEigenVector<Eigen::VectorXd, Eigen::Vector2d>(A, x); }
    if( x.type()==typeid(Eigen::Vector3d) ) { return ApplyEigenVector<Eigen::VectorXd, Eigen::Vector3d>(A, x); }
    if( x.type()==typeid(Eigen::Vector4d) ) { return ApplyEigenVector<Eigen::VectorXd, Eigen::Vector4d>(A, x); }

    return ApplyEigenVector<Eigen::VectorXd, Eigen::VectorXd>(A, x);
  }
  if( A.type()==typeid(Eigen::VectorXf) ) {
    if( x.type()==typeid(Eigen::Vector2f) ) { return ApplyEigenVector<Eigen::VectorXf, Eigen::Vector2f>(A, x); }
    if( x.type()==typeid(Eigen::Vector3f) ) { return ApplyEigenVector<Eigen::VectorXf, Eigen::Vector3f>(A, x); }
    if( x.type()==typeid(Eigen::Vector4f) ) { return ApplyEigenVector<Eigen::VectorXf, Eigen::Vector4f>(A, x); }

    return ApplyEigenVector<Eigen::VectorXf, Eigen::VectorXf>(A, x);
  }
  if( A.type()==typeid(Eigen::VectorXi) ) {
    if( x.type()==typeid(Eigen::Vector2i) ) { return ApplyEigenVector<Eigen::VectorXi, Eigen::Vector2i>(A, x); }
    if( x.type()==typeid(Eigen::Vector3i) ) { return ApplyEigenVector<Eigen::VectorXi, Eigen::Vector3i>(A, x); }
    if( x.type()==typeid(Eigen::Vector4i) ) { return ApplyEigenVector<Eigen::VectorXi, Eigen::Vector4i>(A, x); }

    return ApplyEigenVector<Eigen::VectorXi, Eigen::VectorXi>(A, x);
  }

  // something went wrong
  assert(false);
  return boost::none;
}

boost::any AnyAlgebra::Apply(boost::any const& A, boost::any const& x) const {
  if( IsScalar(A.type()) ) { return ApplyScalar(A, x); }

  if( IsEigenVector(A.type()) ) { return ApplyEigenVector(A, x); }
  
  return ApplyImpl(A, x);
}

boost::any AnyAlgebra::ApplyImpl(boost::any const& A, boost::any const& x) const {
  std::cerr << std::endl << "ERROR: No way to apply " << boost::core::demangle(A.type().name()) << " type to type " << boost::core::demangle(x.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::ApplyImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::ApplyImpl()" << std::endl << std::endl;
  assert(false);
  
  return boost::none;
}
