#include "MUQ/Modeling/LinearAlgebra/LOBPCG.h"

#include "MUQ/Utilities/RandomGenerator.h"

#include <boost/property_tree/ptree.hpp>

using namespace muq::Modeling;
using namespace muq::Utilities;

LOBPCG::LOBPCG(int    numEigsIn,
               double tolIn,
               int    maxItsIn,
               bool   largestIn,
               int    verbosityIn) : numEigs(numEigsIn),
                                     tol(tolIn),
                                     maxIts(maxItsIn),
                                     largest(largestIn),
                                     verbosity(verbosityIn)
{
  assert(numEigs>0);
}


LOBPCG::LOBPCG(boost::property_tree::ptree const& opts) : LOBPCG(opts.get<int>("NumEigs"),
                                                                 opts.get("Tolerance",-1.0),
                                                                 opts.get("MaxIts", -1),
                                                                 opts.get("Largest",true),
                                                                 opts.get("Verbosity",0))
{}

Eigen::MatrixXd LOBPCG::Orthonormalizer::Compute(Eigen::Ref<const Eigen::MatrixXd> const& V)
{
  Eigen::MatrixXd output(V);
  if(B!=nullptr){
    ComputeInPlace(output,B->Apply(V));
  }else{
    ComputeInPlace(output,V);
  }
  return output;
}

Eigen::MatrixXd LOBPCG::Orthonormalizer::Compute(Eigen::Ref<const Eigen::MatrixXd> const& V, Eigen::Ref<const Eigen::MatrixXd> const& BVin)
{
  Eigen::MatrixXd output(V);
  ComputeInPlace(output, BVin);
  return output;
}

void LOBPCG::Orthonormalizer::ComputeInPlace(Eigen::Ref<Eigen::MatrixXd> V)
{
  if(B!=nullptr){
    ComputeInPlace(V, B->Apply(V));
  }else{
    ComputeInPlace(V,V);
  }
}

void LOBPCG::Orthonormalizer::ComputeInPlace(Eigen::Ref<Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXd> const& BVin)
{
  vDim = V.cols();

  VBV_Chol = (V.transpose()*BVin).eval().llt().matrixL();
  V = VBV_Chol.triangularView<Eigen::Lower>().solve(V.transpose()).transpose();

  if(B!=nullptr){
    BV = VBV_Chol.triangularView<Eigen::Lower>().solve(BV.transpose()).transpose();
  }
}

Eigen::MatrixXd LOBPCG::Orthonormalizer::InverseVBV() const
{
  return VBV_Chol.triangularView<Eigen::Lower>().solve(Eigen::MatrixXd::Identity(vDim,vDim));
}


LOBPCG::Constraints::Constraints(std::shared_ptr<LinearOperator> const& B,
                                 Eigen::MatrixXd                 const& constMat) : Y(constMat)
{
  if(B!=nullptr){
    BY = B->Apply(constMat);
  }else{
    BY = constMat;
  }

  YBY_llt = (constMat.transpose() * BY).eval().llt();
}

void LOBPCG::Constraints::ApplyInPlace(Eigen::MatrixXd& x)
{
    Eigen::MatrixXd YBX = BY.transpose()*x;
    x -= Y*YBY_llt.solve(YBX);
}

void LOBPCG::SortVec(std::vector<std::pair<int,int>> const& swapInds,
                     Eigen::Ref<Eigen::VectorXd>            matrix)
{
  for(auto& swap : swapInds)
    std::swap(matrix(swap.first),matrix(swap.second));
}

void LOBPCG::SortCols(std::vector<std::pair<int,int>> const& swapInds,
                      Eigen::Ref<Eigen::MatrixXd>            matrix)
{
  for(auto& swap : swapInds)
    matrix.col(swap.first).swap(matrix.col(swap.second));
}


std::vector<std::pair<int,int>> LOBPCG::GetSortSwaps(Eigen::Ref<const Eigen::VectorXd> const& resids)
{
  Eigen::VectorXd newResids = resids;
  const unsigned int size = resids.size();

  std::vector<std::pair<int,int>> swaps;
  unsigned int maxInd;

  for(unsigned int i=0; i<size-1; ++i)
  {
    // Find the maximum index
    maxInd = std::distance(newResids.data(), std::max_element(&newResids(i), newResids.data()+size));

    // Swap indices if needed
    if(maxInd!=i){
      swaps.push_back( std::make_pair(i, maxInd) );
      std::swap(newResids(i), newResids(maxInd));
    }
  }

  return swaps;
}

LOBPCG& LOBPCG::compute(std::shared_ptr<LinearOperator> const& A,
                        Eigen::MatrixXd                 const& constMat,
                        std::shared_ptr<LinearOperator>        B,
                        std::shared_ptr<LinearOperator>        M)
{

  const int dim = A->rows();
  if(numEigs>int(std::ceil(0.2*dim))){
    std::cerr << "ERROR: LOBPCG can fail when more than 20\% of the eigenvalues are computed." << std::endl;
    assert(numEigs>int(std::ceil(0.2*dim)));
  }

  // Initial message
  if(verbosity>0){
    std::cout << "Solving ";
    if(B==nullptr){
      std::cout << "standard ";
    }else{
      std::cout << "generalized ";
    }
    std::cout << "eigenvalue problem." << std::endl;

    if(verbosity>1)
      std::cout << "  matrix size = " << A->rows() << " x " << A->cols() << std::endl;
  }

  // If this is the first call to compute, initialize everything
  if(eigVals.size()==0)
    reset(dim); // This will initialize the eigenvectors to random values

  Eigen::MatrixXd X = eigVecs;

  // To start with, all of the indices are active
  unsigned int numActive = numEigs;

  // Get ready to apply constraints (e.g., build Y, BY, )
  Constraints consts(B, constMat);
  if(consts.size()>0)
    consts.ApplyInPlace(X); // this will remove components of the eigenvectors that lie in the constraint subspace

  // Make the current vectors orthogonal wrt B
  Orthonormalizer bOrtho(B);
  bOrtho.ComputeInPlace(X);

  Eigen::MatrixXd BX;
  if(B != nullptr){
    BX = bOrtho.GetBV();
  }else{
    BX = X;
  }

  // Compute the initial Ritz vectors by solving the reduced eigenproblem.
  Eigen::MatrixXd AX = A->Apply(X);

  Eigen::MatrixXd XAX = X.transpose() * AX;

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ritzSolver(XAX);
  Eigen::VectorXd subEigVals = ritzSolver.eigenvalues();
  Eigen::MatrixXd const& subEigVecs = ritzSolver.eigenvectors();

  X = (X*subEigVecs).eval();
  AX = (AX*subEigVecs).eval();
  if(B!=nullptr){
    BX = (BX*subEigVecs).eval();
  }else{
    BX = X;
  }

  //
  Eigen::MatrixXd resids(X.rows(),X.cols());
  Eigen::VectorXd residNorms;
  Eigen::MatrixXd AR, P, BP, AP;
  AR.resize(X.rows(),X.cols());

  for(unsigned it=0; it<maxIts; ++it){

    if(verbosity>1)
      std::cout << "Iteration " << it << std::endl;

    // Compute the residuals
    resids.leftCols(numActive) = AX.leftCols(numActive) - BX.leftCols(numActive)*subEigVals.head(numActive).asDiagonal();

    //
    // Compute the norm of the residuals
    residNorms = resids.colwise().norm();

    // Figure out how many indices are active
    numActive=0;
    for(unsigned int i=0; i<residNorms.size(); ++i){
      if(residNorms(i)>tol){
        numActive++;
      }
    }


    if(verbosity>2)
      std::cout << "  eigenvalues = " << subEigVals.transpose() << std::endl;

    if(verbosity>1)
      std::cout << "  residNorms = " << residNorms.transpose() << std::endl;

      
    // If we've converged for all the vectors, break
    if(numActive==0){
      eigVecs = X;
      eigVals = subEigVals;
      return *this;
    }

    auto swaps = GetSortSwaps(residNorms);

    // Sort everything so the first columns corresponds to the largest residuals
    SortVec(swaps, subEigVals.head(numActive));
    SortCols(swaps, X.leftCols(numActive));
    SortCols(swaps, AX.leftCols(numActive));
    SortCols(swaps, BX.leftCols(numActive));
    SortCols(swaps, resids);

    if(it!=0){
      SortCols(swaps, P.leftCols(numActive));
      SortCols(swaps, AP.leftCols(numActive));
      SortCols(swaps, BP.leftCols(numActive));
    }

    // Make residuals B-orthogonal to current values
    //std::cout << "  X.dot(X) = \n" << X.transpose()*X << std::endl;
    resids -= (X.leftCols(numActive)*(BX.leftCols(numActive).transpose()*resids)).eval();


    // Orthonormalize residuals wrt B-inner product
    bOrtho.ComputeInPlace(resids.leftCols(numActive));
    Eigen::MatrixXd BR;
    if(B != nullptr){
      BR = bOrtho.GetBV();
    } else {
      BR = resids.leftCols(numActive);
    }

    AR.leftCols(numActive) = A->Apply(resids.leftCols(numActive));

    if(it!=0){
      // Orthonormalize P
      bOrtho.ComputeInPlace(P.leftCols(numActive),BP.leftCols(numActive));
      AP.leftCols(numActive) = bOrtho.VBV_Chol.triangularView<Eigen::Lower>().solve(AP.leftCols(numActive).transpose()).transpose();

      if(B!=nullptr){
        BP = B->Apply(P);
      }else{
        BP = P;
      }
    }

    Eigen::MatrixXd XAR = X.transpose()*AR.leftCols(numActive);
    Eigen::MatrixXd RAR = resids.leftCols(numActive).transpose()*AR.leftCols(numActive);
    Eigen::MatrixXd XAX = X.transpose()*AX;
    Eigen::MatrixXd XBX, RBR, XBR;

    XBX = BX.transpose() * X;
    RBR = BR.transpose() * resids.leftCols(numActive);
    XBR = X.transpose() * BR.leftCols(numActive);

    Eigen::MatrixXd gramA(3*numEigs,3*numEigs);
    Eigen::MatrixXd gramB(3*numEigs,3*numEigs);

    if(it!=0){
      gramA.resize(numEigs+2*numActive,numEigs+2*numActive);
      gramB.resize(numEigs+2*numActive,numEigs+2*numActive);

      Eigen::MatrixXd XAP = X.transpose() * AP.leftCols(numActive);
      Eigen::MatrixXd RAP = resids.leftCols(numActive).transpose()*AP.leftCols(numActive);
      Eigen::MatrixXd PAP = P.leftCols(numActive).transpose()*AP.leftCols(numActive);
      Eigen::MatrixXd XBP = X.transpose()*BP.leftCols(numActive);
      Eigen::MatrixXd RBP = resids.leftCols(numActive).transpose()*BP.leftCols(numActive);

      gramA.block(0,0,numEigs,numEigs) = subEigVals.asDiagonal();
      gramA.block(0,numEigs,numEigs,numActive) = XAR;
      gramA.block(0,numEigs+numActive,numEigs,numActive) = XAP;
      gramA.block(numEigs,0,numActive,numEigs) = XAR.transpose();
      gramA.block(numEigs,numEigs,numActive,numActive) = RAR;
      gramA.block(numEigs,numEigs+numActive,numActive,numActive) = RAP;
      gramA.block(numEigs+numActive,0,numActive,numEigs) = XAP.transpose();
      gramA.block(numEigs+numActive,numEigs,numActive,numActive) = RAP.transpose();
      gramA.block(numEigs+numActive,numEigs+numActive,numActive,numActive) = PAP;

      gramB.block(0,0,numEigs,numEigs) = Eigen::MatrixXd::Identity(numEigs,numEigs);
      gramB.block(0,numEigs,numEigs,numActive) = XBR;
      gramB.block(0,numEigs+numActive,numEigs,numActive) = XBP;
      gramB.block(numEigs,0,numActive,numEigs) = XBR.transpose();
      gramB.block(numEigs,numEigs,numActive,numActive) = Eigen::MatrixXd::Identity(numActive,numActive);
      gramB.block(numEigs,numEigs+numActive,numActive,numActive) = RBP;
      gramB.block(numEigs+numActive,0,numActive,numEigs) = XBP.transpose();
      gramB.block(numEigs+numActive,numEigs,numActive,numActive) = RBP.transpose();
      gramB.block(numEigs+numActive,numEigs+numActive,numActive,numActive) = Eigen::MatrixXd::Identity(numActive,numActive);

    }else{

      gramA.resize(2*numEigs,2*numEigs);
      gramB.resize(2*numEigs,2*numEigs);

      gramA.block(0,0,numEigs,numEigs) = subEigVals.asDiagonal();
      gramA.block(0,numEigs,numEigs,numEigs) = XAR;
      gramA.block(numEigs,0,numEigs,numEigs) = XAR.transpose();
      gramA.block(numEigs,numEigs,numEigs,numEigs) = RAR;

      gramB.block(0,0,numEigs,numEigs) = Eigen::MatrixXd::Identity(numEigs,numEigs);
      gramB.block(0,numEigs,numEigs,numEigs) = XBR;
      gramB.block(numEigs,0,numEigs,numEigs) = XBR.transpose();
      gramB.block(numEigs,numEigs,numEigs,numEigs) = Eigen::MatrixXd::Identity(numEigs,numEigs);
    }

    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ritzSolver(gramA, gramB);
    subEigVals = ritzSolver.eigenvalues().tail(numEigs);
    assert(subEigVals(0)<subEigVals(1));
    Eigen::Ref<const Eigen::MatrixXd> subEigVecs = ritzSolver.eigenvectors().rightCols(numEigs);


    Eigen::Ref<const Eigen::MatrixXd> eigX = subEigVecs.topRows(numEigs);
    Eigen::Ref<const Eigen::MatrixXd> eigR = subEigVecs.block(numEigs,0,numEigs,subEigVecs.cols());

    Eigen::MatrixXd pp  = resids.leftCols(numActive) * eigR;
    Eigen::MatrixXd app = AR.leftCols(numActive)*eigR;
    Eigen::MatrixXd bpp = BR.leftCols(numActive)*eigR;

    if(it!=0){
      Eigen::Ref<const Eigen::MatrixXd> eigP = subEigVecs.bottomRows(numEigs);

      pp  += P.leftCols(numActive)*eigP;
      app += AP.leftCols(numActive)*eigP;
      bpp += BP.leftCols(numActive)*eigP;
    }

    X = (X*eigX + pp).eval();
    AX = (AX*eigX + app).eval();
    BX = (BX*eigX + bpp).eval();

    //if(B!=nullptr){
    //  BX = (BX*eigX + bpp).eval();
    //}else{
    //  BX = X;
    //}


    P = pp;
    AP = app;
    BP = bpp;

  }

  return *this;
}


LOBPCG& LOBPCG::compute(std::shared_ptr<LinearOperator> const& A,
                        std::shared_ptr<LinearOperator>        B,
                        std::shared_ptr<LinearOperator>        M)
{
  return compute(A,Eigen::MatrixXd(), B,M);
};

LOBPCG& LOBPCG::reset(int dim)
{
  if(tol<0.0)
    tol = dim * std::sqrt(std::numeric_limits<double>::epsilon());

  if(maxIts<0)
    maxIts = std::min(dim,20);

  eigVals.resize(numEigs);
  eigVecs = RandomGenerator::GetNormal(dim,numEigs);

  return *this;
}
