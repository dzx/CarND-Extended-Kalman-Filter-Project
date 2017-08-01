#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  H_laser_ << 1, 0, 0, 0,
	      0, 1, 0, 0;
  
  VectorXd x = VectorXd(4);
  MatrixXd p = MatrixXd(4,4);
  p << 1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1000, 0,
      0, 0, 0, 1000;
  MatrixXd f = MatrixXd(4, 4);
  f << 1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1;

  MatrixXd q = MatrixXd(4, 4);
  
  ekf_.Init(x, p, f, H_laser_, R_laser_, q);
  
  


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 0, 0;
    float x_pos, y_pos;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      x_pos = measurement_pack.raw_measurements_(0) * cos(measurement_pack.raw_measurements_(1));
      y_pos = measurement_pack.raw_measurements_(0) * sin(measurement_pack.raw_measurements_(1));

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_pos = measurement_pack.raw_measurements_(0);
      y_pos = measurement_pack.raw_measurements_(1);
    }
    previous_timestamp_ = measurement_pack.timestamp_;
    ekf_.x_ << x_pos, y_pos, 0, 0;
    

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  // update transition matrix with dt
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  // update process covariance matrix
  float dt4 = pow(dt, 4) / 4;
  float dt3 = pow(dt, 3) / 2;
  float dt2 = pow(dt, 2);
  float dt4x = dt4 * noise_ax;
  float dt3x = dt3 * noise_ax;
  float dt4y = dt3 * noise_ay;
  float dt3y = dt3 * noise_ay;
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt4x, 0, dt3x, 0,
	    0, dt4y, 0, dt3y,
	    dt3x, 0, (dt2 * noise_ax), 0,
	    0, dt3y, 0, (dt2 * noise_ay);
  
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/


  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.R_ = R_radar_;
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    
  } else {
    // Laser updates
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
