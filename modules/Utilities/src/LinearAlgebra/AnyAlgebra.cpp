#include "MUQ/Utilities/LinearAlgebra/AnyAlgebra.h"

using namespace muq::Utilities;

AnyAlgebra::AnyAlgebra() {}

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
  if( ScalarAlgebra::IsScalar(obj.type()) ) { return 1; }

  // get the size of an Eigen::VectorXd
  if( EigenVectorAlgebra::IsEigenVector(obj.type()) ) { return EigenVectorAlgebra::Size(obj); }

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
  if( ScalarAlgebra::IsScalar(type) ) { return ZeroScalar(type); }
  
  if( EigenVectorAlgebra::IsEigenVector(type) ) { return EigenVectorAlgebra::Zero(type, rows); }

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

double AnyAlgebra::EigenMatrixNorm(boost::any const& mat) const {
  if( typeid(Eigen::Matrix2d)==mat.type() ) { return EigenNorm<Eigen::Matrix2d>(mat); }
  if( typeid(Eigen::Matrix3d)==mat.type() ) { return EigenNorm<Eigen::Matrix3d>(mat); }
  if( typeid(Eigen::Matrix4d)==mat.type() ) { return EigenNorm<Eigen::Matrix4d>(mat); }
  
  return EigenNorm<Eigen::MatrixXd>(mat);
}

double AnyAlgebra::Norm(boost::any const& obj) const {
  if( ScalarAlgebra::IsScalar(obj.type()) ) { return ScalarAlgebra::Norm(obj); }

  if( EigenVectorAlgebra::IsEigenVector(obj.type()) ) { return EigenVectorAlgebra::Norm(obj); }

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

double AnyAlgebra::InnerProduct(boost::any const& vec1, boost::any const& vec2) const {
  if( ScalarAlgebra::IsScalar(vec1.type()) && ScalarAlgebra::IsScalar(vec2.type()) ) { return ScalarAlgebra::InnerProduct(vec1, vec2); }

  if( EigenVectorAlgebra::IsEigenVector(vec1.type()) && EigenVectorAlgebra::IsEigenVector(vec2.type()) ) { return EigenVectorAlgebra::InnerProduct(vec1, vec2); }

  return InnerProductImpl(vec1, vec2);
}

double AnyAlgebra::InnerProductImpl(boost::any const& vec1, boost::any const& vec2) const {
  std::cerr << std::endl << "ERROR: Cannot compute the inner product between vectors with types " << boost::core::demangle(vec1.type().name()) << " and " << boost::core::demangle(vec2.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::InnerProductImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::InnerProductImpl()" << std::endl << std::endl;
  assert(false);

  return 0.0;
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
  if( ScalarAlgebra::IsScalar(obj.type()) ) { return ScalarAlgebra::IsZero(obj); }

  if( EigenVectorAlgebra::IsEigenVector(obj.type()) ) { return EigenVectorAlgebra::IsZero(obj); }
  
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

#if MUQ_HAS_SUNDIALS==1
boost::any AnyAlgebra::AccessSundialsVector(N_Vector const& vec, unsigned int const i) const {
  // check the size
    assert(i<NV_LENGTH_S(vec));

    // return the ith element
    return NV_Ith_S(vec, i);
}
#endif

boost::any AnyAlgebra::AccessElement(boost::any const& obj, unsigned int const i, unsigned int const j) const {
  if( ScalarAlgebra::IsScalar(obj.type()) ) { return obj; }

  if( EigenVectorAlgebra::IsEigenVector(obj.type()) ) { return EigenVectorAlgebra::AccessElement(obj, i); }

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

boost::any AnyAlgebra::Identity(std::type_info const& type, unsigned int const rows, unsigned int const cols) const {
  if( ScalarAlgebra::IsScalar(type) ) { return ScalarAlgebra::Identity(type); }

  if( EigenVectorAlgebra::IsEigenVector(type) ) { return EigenVectorAlgebra::Identity(type, rows, cols); }

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
  if( ScalarAlgebra::IsScalar(in0.type()) && ScalarAlgebra::IsScalar(in1.type()) ) { return ScalarAlgebra::Add(in0, in1); }

  if( EigenVectorAlgebra::IsEigenVector(in0.type()) && EigenVectorAlgebra::IsEigenVector(in1.type()) ) { return EigenVectorAlgebra::Add(in0, in1); }

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
  if( ScalarAlgebra::IsScalar(in0.type()) && ScalarAlgebra::IsScalar(in1.type()) ) { return ScalarAlgebra::Subtract(in0, in1); }

  if( EigenVectorAlgebra::IsEigenVector(in0.type()) && EigenVectorAlgebra::IsEigenVector(in1.type()) ) { return EigenVectorAlgebra::Subtract(in0, in1); }

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
  if( ScalarAlgebra::IsScalar(in0.type()) && ScalarAlgebra::IsScalar(in1.type()) ) { return ScalarAlgebra::Multiply(in0, in1); }

  if( ScalarAlgebra::IsScalar(in0.type()) && EigenVectorAlgebra::IsEigenVector(in1.type()) ) { return EigenVectorAlgebra::ScalarMultiply(in0, in1); }
  if( ScalarAlgebra::IsScalar(in1.type()) && EigenVectorAlgebra::IsEigenVector(in0.type()) ) { return EigenVectorAlgebra::ScalarMultiply(in1, in0); }

  if( IsEigenMatrix(in0.type()) && IsEigenMatrix(in1.type()) ) { return MultiplyEigenMatrix(in0, in1); }

  if( ScalarAlgebra::IsScalar(in0.type()) && IsEigenMatrix(in1.type()) ) { return MultiplyEigenMatrixScalar(in0, in1); }
  if( ScalarAlgebra::IsScalar(in1.type()) && IsEigenMatrix(in0.type()) ) { return MultiplyEigenMatrixScalar(in1, in0); }

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

boost::any AnyAlgebra::ApplyInverse(boost::any const& A, boost::any const& x) const {
  if( ScalarAlgebra::IsScalar(A.type()) ) { return Multiply(Inverse(A), x); }

  if( EigenVectorAlgebra::IsEigenVector(A.type()) ) { return EigenVectorAlgebra::ApplyInverse(A, x); }
  
  return ApplyInverseImpl(A, x);
}

boost::any AnyAlgebra::ApplyInverseImpl(boost::any const& A, boost::any const& x) const {
  std::cerr << std::endl << "ERROR: No way to apply the inverse " << boost::core::demangle(A.type().name()) << " type to type " << boost::core::demangle(x.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::ApplyInverseImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::ApplyInverseImpl()" << std::endl << std::endl;
  assert(false);
  
  return boost::none;
}

boost::any AnyAlgebra::Apply(boost::any const& A, boost::any const& x) const {
  if( ScalarAlgebra::IsScalar(A.type()) ) { return Multiply(A, x); }

  if( EigenVectorAlgebra::IsEigenVector(A.type()) ) { return EigenVectorAlgebra::Apply(A, x); }
  
  return ApplyImpl(A, x);
}

boost::any AnyAlgebra::ApplyImpl(boost::any const& A, boost::any const& x) const {
  std::cerr << std::endl << "ERROR: No way to apply " << boost::core::demangle(A.type().name()) << " type to type " << boost::core::demangle(x.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::ApplyImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::ApplyImpl()" << std::endl << std::endl;
  assert(false);
  
  return boost::none;
}

boost::any AnyAlgebra::Inverse(boost::any const& obj) const {
  if( ScalarAlgebra::IsScalar(obj.type()) ) { return ScalarAlgebra::Inverse(obj); }
  
  return InverseImpl(obj);
}

boost::any AnyAlgebra::InverseImpl(boost::any const& obj) const {
  std::cerr << std::endl << "ERROR: No way to compute the inverse of type " << boost::core::demangle(obj.type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::inverseImpl()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::InverseImpl()" << std::endl << std::endl;
  assert(false);
  
  return boost::none;
}
