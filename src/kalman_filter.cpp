#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
  
}

void KalmanFilter::Update(const VectorXd &z) {

  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  UpdateCommon(y);

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  VectorXd z_pred = VectorXd(3); 
  double rho = sqrt(x_(0) * x_(0) + x_(1) * x_(1));
  double phi = atan2(x_(1), x_(0));
  double rho_dot;
  if(fabs(rho) > .00001){
    rho_dot = (x_(0) * x_(2) + x_(1) * x_(3)) / rho;
  }else{
    rho_dot = sqrt(x_(2) * x_(2) + x_(3) * x_(3)); // for now
  }
  z_pred << rho, phi, rho_dot;
  VectorXd y = z - z_pred;  
  UpdateCommon(y);
  
}

void KalmanFilter::UpdateCommon(const Eigen::VectorXd &y){
  MatrixXd Ht = H_.transpose(); 
  MatrixXd S = H_ * P_ * Ht + R_; 
  MatrixXd K = (P_ * Ht) * S.inverse(); 
  
  x_ = x_ + K * y;
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I - K * H_) * P_;
}
