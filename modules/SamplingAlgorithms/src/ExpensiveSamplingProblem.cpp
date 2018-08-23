#include "MUQ/SamplingAlgorithms/ExpensiveSamplingProblem.h"

using namespace muq::SamplingAlgorithms;

ExpensiveSamplingProblem::ExpensiveSamplingProblem(std::shared_ptr<muq::Modeling::ModPiece> target) : SamplingProblem(target) {}

ExpensiveSamplingProblem::~ExpensiveSamplingProblem() {}
