#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.7;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  is_initialized_ = false;

  n_x_ = 5;

  n_aug_ = 7;

  lambda_ = 3 - n_aug_;

  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {
    weights_(i) = 0.5 / (n_aug_ + lambda_);
  }

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (is_initialized_ == false) {
    P_ = MatrixXd::Identity(5, 5);
    time_us_ = meas_package.timestamp_;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) { 
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      x_ << rho * cos(phi), rho * sin(phi), 0, 0, 0;
    } else {
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    };
    is_initialized_ = true;
    return;
  };

  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  Prediction(delta_t);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else {
    UpdateLidar(meas_package);
  };

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  //Predict sigma points

  for (int col = 0; col < 2 * n_aug_ + 1; col++) {
      VectorXd x = Xsig_aug.col(col);
      double px = x(0);
      double py = x(1);
      double v = x(2);
      double psi = x(3);
      double psi_dot = x(4);
      double nu_a = x(5);
      double nu_psi_dotdot = x(6);
      VectorXd dx = VectorXd(n_x_);
      VectorXd noise_x = VectorXd(n_x_);
      
      if (psi_dot < 0.15 && psi_dot > -0.15) {        
          dx << v * cos(psi) * delta_t, v * sin(psi) * delta_t, 0, psi_dot * delta_t, 0;
      } else {
          dx << v / psi_dot * (sin(psi + psi_dot * delta_t) - sin(psi)),
                v / psi_dot * (-cos(psi + psi_dot * delta_t) + cos(psi)),
                0, psi_dot * delta_t, 0;
      }
      
      noise_x << 0.5 * delta_t * delta_t * cos(psi) * nu_a,
                 0.5 * delta_t * delta_t * sin(psi) * nu_a,
                 delta_t * nu_a,
                 0.5 * delta_t * delta_t * nu_psi_dotdot,
                 delta_t * nu_psi_dotdot;
      
      Xsig_pred_.col(col) = x.head(n_x_) + dx + noise_x;
  };

  VectorXd x_pred = VectorXd(n_x_);
  MatrixXd P_pred = MatrixXd(n_x_, n_x_);

  //predicted state mean
  x_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_pred = x_pred + weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_pred;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_pred = P_pred + weights_(i) * x_diff * x_diff.transpose() ;
  }

 x_ = x_pred;
 P_ = P_pred;

  //cout << "Prediction: " << endl;
  //cout << "x_ = " << x_ << endl;
  //cout << "P_ = " << P_ << endl;


}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  VectorXd z = VectorXd(2);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];

  MatrixXd H = MatrixXd(2, n_x_);
  H.fill(0);
  H(0, 0) = 1;
  H(1, 1) = 1;

  MatrixXd R_las = MatrixXd(2, 2);
  R_las.fill(0);
  R_las(0, 0) = std_laspx_ * std_laspx_;
  R_las(1, 1) = std_laspy_ * std_laspy_;

  // KF Measurement update step
  VectorXd y = z - H * x_;
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P_ * Ht + R_las;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  //new state
  x_ = x_ + (K * y);
  P_ = (I - K * H) * P_;

  //cout << "Update Lidar: " << endl;
  //cout << "x_ = " << x_ << endl;
  //cout << "P_ = " << P_ << endl;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  MatrixXd R_rad = MatrixXd(3, 3);
  R_rad.fill(0);
  R_rad(0, 0) = std_radr_ * std_radr_;
  R_rad(1, 1) = std_radphi_ * std_radphi_;
  R_rad(2, 2) = std_radrd_ * std_radrd_;

  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      double px = Xsig_pred_(0, i);
      double py = Xsig_pred_(1, i);
      double v = Xsig_pred_(2, i);
      double psi = Xsig_pred_(3, i);
      
      double rho = sqrt(px*px + py*py);
      if (rho < 1) {
        rho += 1;
      };
      double phi = atan2(py, px);
      double rho_dot = (px * cos(psi) * v + py * sin(psi) * v) / rho;
      
      Zsig(0, i) = rho;
      Zsig(1, i) = phi;
      Zsig(2, i) = rho_dot;
  };
  
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  };

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      MatrixXd diff_z = Zsig.col(i) - z_pred;
      
      //angle normalization
      while (diff_z(1)> M_PI) diff_z(1)-=2.*M_PI;
      while (diff_z(1)<-M_PI) diff_z(1)+=2.*M_PI;

      S = S + weights_(i) * diff_z * diff_z.transpose();
  };
  
  S = S + R_rad;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      Tc = Tc + weights_(i) * (Xsig_pred_.col(i) - x_) * (Zsig.col(i) - z_pred).transpose();
  };
  
  MatrixXd K = Tc * S.inverse();  

  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];

  x_ = x_ + K * (z - z_pred);

  P_ = P_ - K * S * K.transpose();  

  /*
  if (x_(1) < -1000 || P_(0, 0) < -10000 || x_(0) != x_(0)) {
    cout << "Xsig_pred_ = " << Xsig_pred_ << endl;
    cout << "K = " << K << endl;
    cout << "S = " << S << endl;
    cout << "Tc = " << Tc << endl;
    cout << "Zsig = " << Zsig << endl;
    cout << "z_pred = " << z_pred << endl;
  }

  cout << "Update Radar: " << endl;
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
  */

}
