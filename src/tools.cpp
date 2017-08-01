#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if(estimations.size() == 0 || estimations.size() != ground_truth.size()){
    return rmse;
  }

  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
    VectorXd se = estimations[i] - ground_truth[i];
    se = se.array() * se.array();
    rmse += se;
  }

  //calculate the mean
  rmse = rmse.array() / estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float pxysq = px * px + py * py;

  //check division by zero
  if(fabs(pxysq) < .0001){
    cout << "Division by zero." << endl;
    return Hj;
  }
  //compute the Jacobian matrix
  float absD = sqrt(pxysq);
  float vpDiff = vx * py - vy * px;
  float absCube = pow(absD, 3);
  //cout << vpDiff << endl;
  Hj << px / absD, py / absD, 0, 0,
    -py / pxysq, px / pxysq, 0, 0,
    py * vpDiff / absCube, px * vpDiff / absCube, px / absD, py / absD;

  return Hj;
}
