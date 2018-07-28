#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = false;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  //initialize to identity (CAN BE MODIFIED FOR TUNING)
  P_ = MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  is_initialized_ = false;

  n_x_ = 5 ;

  n_aug_ = 7 ;

  lambda_ = 3 - n_aug_;
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

  /**
  Initialize state
  */
  if(!is_initialized_){
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
      //initialize using radar
      float x = meas_package.raw_measurements_[0]*cos(meas_package.raw_measurements_[1]);
      float y = -meas_package.raw_measurements_[0]*sin(meas_package.raw_measurements_[1]);
      x_ << x,y,0,0,0; //note: we cannot initialize v with the radar velocity as it's a different quantity.
      time_us_ = meas_package.timestamp_;
      is_initialized_ = true;


    } else if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
      //initialize using lidar
      x_ << meas_package.raw_measurements_[0],meas_package.raw_measurements_[1],0,0,0; //TODO
      time_us_ = meas_package.timestamp_;
      is_initialized_ = true;
    }

  }


  //predict and update
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
    float dt = meas_package.timestamp_ - time_us_;
    time_us_ = meas_package.timestamp_; 

    Prediction(dt);
    UpdateRadar(meas_package);

  } else if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_){
    float dt = meas_package.timestamp_ - time_us_;
    time_us_ = meas_package.timestamp_; 

    Prediction(dt);
    UpdateLidar(meas_package);
  }

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


  /**
  Generate sigma points
  */

  //preliminary calculations and initializations
  double sq = sqrt(lambda_ + n_x_);  //spread factor (square root of lambda + number of dim of state vector)
  MatrixXd A = P_.llt().matrixL();  //square root of the covariance matrix P
  MatrixXd mean_columns(n_x_,n_x_);    //term corresponding to mean in next 2*n_x_ sigma points
  mean_columns<<x_, x_, x_, x_, x_ ;

  //actual sigma points generation
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);
  Xsig << x_,mean_columns + sq * A,mean_columns - sq * A;

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(n_x_) = x_;
  x_aug(n_x_) = 0;  // mean of the noise is always zero
  x_aug(n_x_ + 1) = 0; // mean of the noise is always zero 

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(n_aug_-2,n_aug_-2) = std_a_*std_a_;
  P_aug(n_aug_-1,n_aug_-1) = std_yawdd_*std_yawdd_;

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  MatrixXd B = P_aug.llt().matrixL(); //square root of P_Aug
  Xsig_aug.col(0)  = x_aug;  //augmented sigma points

  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1)     = x_aug + sqrt(lambda_+n_aug_) * B.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * B.col(i);
  }

  /**
  predict mean and covariance
   */

  //predict sigma points
  //avoid division by zero
  //write predicted sigma points into right column

  double px; 
  double py;
  double v;
  double sigma;
  double sigma_dot;
  double va;
  double vsigma;
  
  for(int i = 0; i<15; i++){ 


    px = Xsig_aug(0,i);
    py = Xsig_aug(1,i);
    v  = Xsig_aug(2,i);
    sigma = Xsig_aug(3,i);
    sigma_dot = Xsig_aug(4,i);
    va = Xsig_aug(5,i);
    vsigma = Xsig_aug(6,i);
    
    if (sigma_dot > 0.001){

      VectorXd update(5);
      update << (v/sigma_dot) * (sin(sigma + sigma_dot*delta_t) - sin(sigma) ) + 0.5*delta_t*delta_t*cos(sigma)*va , (v/sigma_dot) * (-cos(sigma + sigma_dot*delta_t) + cos(sigma) ) + 0.5*delta_t*delta_t*sin(sigma)*va, delta_t * va , sigma_dot*delta_t  + 0.5*delta_t*delta_t*vsigma, delta_t*vsigma;
      Xsig_pred_.col(i) = Xsig_aug.col(i).head(5) + update;

    } else {
      VectorXd update(5);
      update << v*cos(sigma)*delta_t + 0.5*delta_t*delta_t*cos(sigma)*va , v*sin(sigma)*delta_t + 0.5*delta_t*delta_t*sin(sigma)*va, delta_t * va , sigma_dot*delta_t  + 0.5*delta_t*delta_t*vsigma, delta_t*vsigma;
      Xsig_pred_.col(i) = Xsig_aug.col(i).head(5) + update;
    }
  }

  // recover mean and covariance from predicted sigma points
  // set weights
  weights_(0) =  lambda_/(lambda_ + n_aug_);
  for (int i=1; i<2*n_aug_+1; i++) {  
    weights_(i) = 0.5/(n_aug_+lambda_);
  }

  //predicted state : 
  //Predict mean
  x_.fill(0.0);
  for(int i=0; i<2*n_aug_+1; i++){
      x_ = x_+ weights_(i) * Xsig_pred_.col(i);
  }

  //Predict covariance
  //normalize angles
  P_.fill(0.0);
  for(int i=0; i<2*n_aug_+1; i++){

    VectorXd x_diff = Xsig_pred_.col(i)- x_;
      
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
      
    P_ = P_+ weights_(i)*x_diff*x_diff.transpose();
  }

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

  //transform state vector into measurement space
  
  //sigma points transformation 

  int n_z = 3; //number of variables in radar space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1); //used to store sigma points in radar measurement space

  //transform each sigma point into radar measurement space applying geometrical transformations
  for(int i = 0; i < 2 * n_aug_ + 1; i++){
      double px = Xsig_pred_.col(i)(0);
      double py = Xsig_pred_.col(i)(1);
      double v = Xsig_pred_.col(i)(2);
      double sig = Xsig_pred_.col(i)(3);
      
      double rho = sqrt(px*px + py*py);
      double phi = atan2(py,px);
      double drho = (px*cos(sig)*v+py*sin(sig)*v)/rho;
      
      VectorXd elements = VectorXd(3); 
      elements << rho,phi,drho;

      //store sigma points for radar space
      Zsig.col(i) = elements;
  } 
  
  //calculate mean predicted measurement based on radar space sigma points
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for(int i=0; i<2*n_aug_+1; i++){
      z_pred = z_pred+ weights_(i) * Zsig.col(i);
  }  
  
  //calculate measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for(int i=0; i<2*n_aug_+1; i++){

    VectorXd z_diff = Zsig.col(i)- z_pred;
    //normalize angles  
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
      
    S = S+ weights_(i)*z_diff*z_diff.transpose();
  }
  
  //add the noise term
  VectorXd diag = VectorXd(n_z);
  diag << std_radr_*std_radr_, std_radphi_*std_radphi_, std_radrd_*std_radrd_;
  MatrixXd R = MatrixXd(n_z,n_z);
  R.fill(0.0);
  R.diagonal() = diag;
  S = S + R;

  /**
  Once we have the prediction in measurement space we can update the UKF using the measurement
  */

  //calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for(int i =0; i<2*n_aug_+1; i++){
      MatrixXd diffx= Xsig_pred_.col(i)-x_;
      
      while (diffx(3)> M_PI) diffx(3)-=2.*M_PI;
      while (diffx(3)<-M_PI) diffx(3)+=2.*M_PI;
      
      MatrixXd diffz= Zsig.col(i)-z_pred;
      
      while (diffz(1)> M_PI) diffz(1)-=2.*M_PI;
      while (diffz(1)<-M_PI) diffz(1)+=2.*M_PI;
      
      Tc += weights_(i)*diffx*diffz.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd S_inv = S.inverse();
  MatrixXd K = Tc*S_inv;
  //update state mean and covariance matrix
  x_= x_+ K*(meas_package.raw_measurements_ - z_pred);
  P_= P_- K*S*K.transpose();
}
