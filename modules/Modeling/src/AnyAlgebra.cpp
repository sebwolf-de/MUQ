#include "MUQ/Modeling/AnyAlgebra.h"

using namespace muq::Modeling;

AnyAlgebra::AnyAlgebra() {}

boost::any AnyAlgebra::IdentityBase(std::reference_wrapper<const boost::any> const& in) const {
  // Eigen::VectorXd type
  if( eigenVecType.compare(in.get().type().name())==0 ) {
    const Eigen::VectorXd& inVec = boost::any_cast<const Eigen::VectorXd&>(in);
    
    return (Eigen::MatrixXd)Eigen::MatrixXd::Identity(inVec.size(), inVec.size());
  }
  
  // double type
  if( doubleType.compare(in.get().type().name())==0 ) {
    return 1.0;
  }
  
  return Identity(in);
}

boost::any AnyAlgebra::Identity(std::reference_wrapper<const boost::any> const& in) const {
  std::cerr << std::endl << "ERROR: No way to compute identy object with type " << boost::core::demangle(in.get().type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::Identity()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::Identity()" << std::endl << std::endl;
  assert(false);
  
  return boost::none;
}

boost::any AnyAlgebra::AddBase(std::reference_wrapper<const boost::any> const& in0, std::reference_wrapper<const boost::any> const& in1) const {
  // the first type is boost::none
  if( noneType.compare(in0.get().type().name())==0 ) {
    // return the second
    return in1.get();
  }

  // both in/out are double
  if( doubleType.compare(in0.get().type().name())==0 && doubleType.compare(in1.get().type().name())==0 ) {
    const double in0d = boost::any_cast<const double>(in0);
    const double in1d = boost::any_cast<const double>(in1);

    return in0d+in1d;
  }

  // both in/out are Eigen::MatrixXd
  if( eigenMatType.compare(in0.get().type().name())==0 && eigenMatType.compare(in1.get().type().name())==0 ) {
    const Eigen::MatrixXd& in0Mat = boost::any_cast<const Eigen::MatrixXd&>(in0);
    const Eigen::MatrixXd& in1Mat = boost::any_cast<const Eigen::MatrixXd&>(in1);

    // make sure the sizes match
    assert(in0Mat.rows()==in1Mat.rows());
    assert(in0Mat.cols()==in1Mat.cols());
	
    return (Eigen::MatrixXd)(in0Mat+in1Mat);
  }

  // both in/out are Eigen::VectorXd
  if( eigenVecType.compare(in0.get().type().name())==0 && eigenVecType.compare(in1.get().type().name())==0 ) {
    const Eigen::VectorXd& in0Vec = boost::any_cast<const Eigen::VectorXd&>(in0);
    const Eigen::VectorXd& in1Vec = boost::any_cast<const Eigen::VectorXd&>(in1);

    // make sure the sizes match
    assert(in0Vec.size()==in1Vec.size());
	
    return (Eigen::VectorXd)(in0Vec+in1Vec);
  }

  return Add(in0, in1);
}

boost::any AnyAlgebra::Add(std::reference_wrapper<const boost::any> const& in0, std::reference_wrapper<const boost::any> const& in1) const {
  std::cerr << std::endl << "ERROR: No way to add type " << boost::core::demangle(in0.get().type().name()) << " and type " << boost::core::demangle(in1.get().type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::Add()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::Add()" << std::endl << std::endl;
  assert(false);
  
  return boost::none;
}

boost::any AnyAlgebra::MultiplyBase(std::reference_wrapper<const boost::any> const& in0, std::reference_wrapper<const boost::any> const& in1) const {
  // both in/out are Eigen::MatrixXd
  if( eigenMatType.compare(in0.get().type().name())==0 && eigenMatType.compare(in1.get().type().name())==0 ) {
    const Eigen::MatrixXd& in0Mat = boost::any_cast<const Eigen::MatrixXd&>(in0);
    const Eigen::MatrixXd& in1Mat = boost::any_cast<const Eigen::MatrixXd&>(in1);

    // make sure the sizes match
    assert(in0Mat.cols()==in1Mat.rows());
	
    return (Eigen::MatrixXd)(in0Mat*in1Mat);
  }

  // double times Eigen::MatrixXd
  if( doubleType.compare(in0.get().type().name())==0 && eigenMatType.compare(in1.get().type().name())==0 ) {
    const double in0Mat = boost::any_cast<double>(in0);
    const Eigen::MatrixXd& in1Mat = boost::any_cast<const Eigen::MatrixXd&>(in1);
    
    // make sure the sizes match
    assert(in1Mat.rows()==1);
	
    return (Eigen::MatrixXd)(in0Mat*in1Mat);
  }

  // Eigen::MatrixXd times double 
  if( eigenMatType.compare(in0.get().type().name())==0 && doubleType.compare(in1.get().type().name())==0 ) {
    const Eigen::MatrixXd& in0Mat = boost::any_cast<const Eigen::MatrixXd&>(in1);
    const double in1Mat = boost::any_cast<double>(in0);
    
    // make sure the sizes match
    assert(in0Mat.cols()==1);
	
    return (Eigen::MatrixXd)(in0Mat*in1Mat);
  }

  return Multiply(in0, in1);
}

boost::any AnyAlgebra::Multiply(std::reference_wrapper<const boost::any> const& in0, std::reference_wrapper<const boost::any> const& in1) const {
  std::cerr << std::endl << "ERROR: No way to multiply type " << boost::core::demangle(in0.get().type().name()) << " and type " << boost::core::demangle(in1.get().type().name()) << std::endl;
  std::cerr << "\tTry overloading boost::any AnyAlgebra::Multiply()" << std::endl << std::endl;
  std::cerr << "\tError in AnyAlgebra::Multiply()" << std::endl << std::endl;
  assert(false);
  
  return boost::none;
}
