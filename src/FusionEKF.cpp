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

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  ekf_.Q_  = MatrixXd(4, 4);
  ekf_.x_ = VectorXd::Ones(4);
  ekf_.F_ = MatrixXd::Identity(4,4);
  ekf_.P_ = MatrixXd::Identity(4,4) * 10e6;
  previous_timestamp_ = 0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  /*
     NOTE: Initialization template provided by Udacity is ignored.

	 Initialization can be done by simply executing only measurement update.
	 As the initial prediction of x_ is unreliable, initial diagnal elements of P_ is set to large value (e.g. 10e6).
    
  */

#if false
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
	if (!is_initialized_) {
		/**
		TODO:
		  * Initialize the state ekf_.x_ with the first measurement.
		  * Create the covariance matrix.
		  * Remember: you'll need to convert radar from polar to cartesian coordinates.
		*/
		// first measurement
		cout << "EKF: " << endl;

		if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
			/**
			Convert radar from polar to cartesian coordinates and initialize state.
			*/
			float rho = measurement_pack.raw_measurements_(0);
			float phi = measurement_pack.raw_measurements_(1);
			float rho_dot = measurement_pack.raw_measurements_(2);
			ekf_.x_ << rho*cos(phi), rho*sin(phi), rho_dot*cos(phi), rho_dot*sin(phi);
		}
		else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
			/**
			Initialize state.
			*/
			ekf_.x_(0) = measurement_pack.raw_measurements_(0);
			ekf_.x_(1) = measurement_pack.raw_measurements_(1);
		}

		previous_timestamp_ = measurement_pack.timestamp_;

		// done initializing, no need to predict or update
		is_initialized_ = true;
		return;
	}
#endif



  /*****************************************************************************
   *  Time Update
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  if (is_initialized_) {
		float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
		previous_timestamp_ = measurement_pack.timestamp_;
		float dt2 = dt*dt;
		ekf_.F_(0, 2) = dt;
		ekf_.F_(1, 3) = dt;

		MatrixXd Q_nu = MatrixXd(2, 2);
		Q_nu << 9, 0,
			0, 9;

		MatrixXd G = MatrixXd(4, 2);
		G << dt2 / 2, 0,
			0, dt2 / 2,
			dt, 0,
			0, dt;

		ekf_.Q_ = G * Q_nu * G.transpose();

		ekf_.Predict();
  }
  else {
		previous_timestamp_ = measurement_pack.timestamp_;
		is_initialized_ = true;
  }
  /*****************************************************************************
   *  Measurement Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.R_ = R_radar_;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
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
