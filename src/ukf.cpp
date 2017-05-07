
#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

//Define sizes of states
#define N_STATES (5)
#define N_AUG_STATES (7)
#define N_RADAR (3)
#define N_LIDAR (2)
#define N_SIGMA_POINTS (2 * N_AUG_STATES + 1)

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  //states not initialized
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd::Zero(N_STATES);

  // initial covariance matrix
  P_ = MatrixXd::Identity(N_STATES, N_STATES);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.25;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.1*M_PI;
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.05;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // initial covariance matrix
  P_ << std_radr_, 0, 0, 0, 0,
        0, std_radr_, 0, 0, 0,
        0, 0, 2*std_radr_, 0, 0, 
        0, 0, 0, 2*M_PI, 0,
        0, 0, 0, 0, 5*M_PI*M_PI;

  //initialize weights
  weights_ = VectorXd(2*N_STATES + 1); 

  //tools
  tools = Tools();

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  
  */
  //compute the time elapsed between the current and previous measurements
  
  double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;  //dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;

  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    cout << "Initializing...\n";
    VectorXd x_previous = x_;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      VectorXd z = meas_package.raw_measurements_;
      x_(0) = z(0)*cos(z(1)); 
      x_(1) = z(0)*sin(z(1));
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
    }
    
    //avoid initializing to exactly zero, as this does not make sense
    if (x_(0)<0.001){
      x_(0)=0.001;
    }
    if (x_(1)<0.001){
      x_(1)=0.001;
    }

    is_initialized_ = true;

    // done initializing, no need to predict or update
    // cout << x_ << "\n";
    return;
    
  }


  //Predict ahead of processing measurement 
  Prediction(dt);

  //process  measurement 
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    if (use_radar_){
      // cout << "Processing Radar Measurement...\n";
      UpdateRadar(meas_package);
    }
    else{
      // cout <<"Skipping Radar...\n";
    }
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    if(use_laser_){
      // cout << "Processing Laser Measurement...\n";
      UpdateLidar(meas_package);
    }
    else{
      // cout << "Skipping Laser\n...";
    }
    
  }
  // cout <<x_ <<f "\n";

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

  /*****************************************************
  * Compute Augmented Sigma Points
  *****************************************************/
  // useful parameters 
  double lambda = 3 - N_AUG_STATES;//spreading parameter
  //declare augmented variables
  VectorXd x_aug = VectorXd(N_AUG_STATES);
  MatrixXd P_aug = MatrixXd(N_AUG_STATES, N_AUG_STATES);
  MatrixXd Xsig_aug = MatrixXd(N_AUG_STATES, N_SIGMA_POINTS);

  //create augmented mean state
  x_aug.head(x_.size()) = x_;
  x_aug.tail(2) << 0,0;
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd A_aug = P_aug.llt().matrixL();
  //create augmented sigma points

  Xsig_aug.col(0) = x_aug;
  double w = sqrt(lambda+N_AUG_STATES);
  for (int i=0; i<N_AUG_STATES; i++){
     Xsig_aug.col(i+1) = x_aug + w*A_aug.col(i);
     // Xsig_aug.col(i+1)(3) = tools.wrapPi(Xsig_aug.col(i+1)(3));
     Xsig_aug.col(i+N_AUG_STATES+1) = x_aug - w*A_aug.col(i);
     // Xsig_aug.col(i+N_AUG_STATES+1)(3) = tools.wrapPi(Xsig_aug.col(i+N_AUG_STATES+1)(3));

  }
  // cout <<Xsig_aug << "\n";

/*******************************************************************************
 * Predict Sigma Points
******************************************************************************/

  Xsig_pred_ = MatrixXd(N_STATES, N_SIGMA_POINTS);
  for (int i = 0; i< N_SIGMA_POINTS; i++){
      VectorXd pnt = Xsig_aug.col(i);
      
      // Add in process noise first (always added in)
      Xsig_pred_.col(i) << 0.5*pow(delta_t,2)*cos(pnt(3))*pnt(5),
                         0.5*pow(delta_t,2)*sin(pnt(3))*pnt(5),
                         delta_t*pnt(5),
                         0.5*pow(delta_t,2)*pnt(6),
                         delta_t*pnt(6);
      // Add previous states:
      Xsig_pred_.col(i) +=pnt.head(N_STATES);

      //avoid division by zero when integrating state. chose correct integration formula
      if (pnt(4)<0.0001){
          Xsig_pred_.col(i)(0) = Xsig_pred_.col(i)(0) + pnt(2)*cos(pnt(3))*delta_t;
          Xsig_pred_.col(i)(1) = Xsig_pred_.col(i)(1) + pnt(2)*sin(pnt(3))*delta_t;
      }
      else{
          Xsig_pred_.col(i)(0) = Xsig_pred_.col(i)(0) + pnt(2)/pnt(4)*
          (sin(pnt(3) + pnt(4)*delta_t) - sin(pnt(3)));
          Xsig_pred_.col(i)(1) = Xsig_pred_.col(i)(1) + pnt(2)/pnt(4)*
          (-cos(pnt(3) + pnt(4)*delta_t) + cos(pnt(3)));
      }
      //third state has no divide by zero worries
      Xsig_pred_.col(i)(3) = Xsig_pred_.col(i)(3) + pnt(4)*delta_t;
      // Xsig_pred_.col(i)(3) = tools.wrapPi(Xsig_pred_.col(i)(3));
  }


/*******************************************************************************
 * Compute  predicted mean and covariance
 ******************************************************************************/
  // Initialize to zero
  x_.setZero();
  P_.setZero();
  
  
  //Calculate Weights:
  weights_=  Eigen::VectorXd::Constant(2*N_AUG_STATES + 1, 1/(2*(lambda + N_AUG_STATES)));
  weights_(0) = lambda/(lambda + N_AUG_STATES);
  // Calculate Predicted Mean
  for (int i = 0; i < weights_.size(); i++){
      x_+= weights_(i)*Xsig_pred_.col(i);
  }
  x_(3) = tools.wrapPi(x_(3));
  // Calculate predicted covariance
  for (int i = 0; i < Xsig_pred_.cols(); i++){
      MatrixXd diff = Xsig_pred_.col(i) - x_;
      diff(3) = tools.wrapPi(diff(3));
      P_ += weights_(i)*diff*diff.transpose();
  }




}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  /**
  Complete this function! Use lidar data to update the belief about the object's
  position.
  You'll also need to calculate the lidar NIS.
  */
  

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(N_LIDAR, N_SIGMA_POINTS);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(N_LIDAR);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(N_LIDAR,N_LIDAR);

  //transform sigma points into measurement space
  for (int i=0; i<N_SIGMA_POINTS; i++){
      // breakout state for code clarity
      double px = Xsig_pred_.col(i)(0);
      double py =  Xsig_pred_.col(i)(1);
      //set predicted measurements
      Zsig.col(i) << px,
                     py;
  }
  
  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < N_SIGMA_POINTS; i++) {  //iterate over sigma points
    z_pred += weights_(i) * Zsig.col(i);
  }

  //predicted measurement covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < N_SIGMA_POINTS; i++) {  //iterate over sigma points
    // state difference
    VectorXd z_diff = Zsig.col(i) - z_pred;

    S +=weights_(i) * z_diff * z_diff.transpose() ;
  }
  S(0,0)+=pow(std_laspx_,2);
  S(1,1)+=pow(std_laspy_,2);
/*******************************************************************************
 * Compute Kalman Gain and update state and covaraince
 ******************************************************************************/
 //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(N_STATES, N_LIDAR);
  Tc.fill(0.0);
  //calculate cross correlation matrix
  for (int i = 0; i < N_SIGMA_POINTS; i++) { 
    //measurement difference
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = tools.wrapPi(x_diff(3));

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd S_inv = S.inverse();
  MatrixXd K = Tc*S_inv;
  
  //update state mean and covariance matrix
  VectorXd Innovation = meas_package.raw_measurements_ - z_pred;

  x_ += K*Innovation;
  x_(3) = tools.wrapPi(x_(3));

  P_ += -K*S*K.transpose();

  VectorXd NIS_laser_vector = Innovation.transpose()*S_inv*Innovation;

  NIS_laser_  = NIS_laser_vector(0);
  // cout << x_;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use radar data to update the belief about the object's
  position.
  You'll also need to calculate the radar NIS.
  */
  

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(N_RADAR, N_SIGMA_POINTS);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(N_RADAR);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(N_RADAR,N_RADAR);

  //transform sigma points into measurement space
  for (int i=0; i<N_SIGMA_POINTS; i++){
      // breakout state for code clarity
      double px = Xsig_pred_.col(i)(0);
      double py =  Xsig_pred_.col(i)(1);
      double v = Xsig_pred_.col(i)(2);
      double omega = Xsig_pred_.col(i)(3); 
      // precalculate some values
      double w1 = sqrt(pow(px,2) + pow(py,2));
      //set predicted measurements
      Zsig.col(i) << w1,
                     atan2(py,px),
                     (v*(px*cos(omega) + py*sin(omega)))/w1;
  }
  
  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < N_SIGMA_POINTS; i++) {  //iterate over sigma points
    z_pred += weights_(i) * Zsig.col(i);
  }

  //predicted measurement covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < N_SIGMA_POINTS; i++) {  //iterate over sigma points

    // state difference
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1) = tools.wrapPi(z_diff(1));

    S +=weights_(i) * z_diff * z_diff.transpose() ;
  }
  S(0,0)+=pow(std_radr_,2);
  S(1,1)+=pow(std_radphi_,2);
  S(2,2)+=pow(std_radrd_,2); 
/*******************************************************************************
 * Compute Kalman Gain and update state and covaraince
 ******************************************************************************/
 //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(N_STATES, N_RADAR);
  Tc.fill(0.0);
  //calculate cross correlation matrix
  for (int i = 0; i < N_SIGMA_POINTS; i++) { 
    //measurement difference
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    z_diff(1) = tools.wrapPi(z_diff(1));
    
    //state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    x_diff(3) = tools.wrapPi(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd S_inv = S.inverse();
  MatrixXd K = Tc*S_inv;
  
  //update state mean and covariance matrix
  VectorXd Innovation = meas_package.raw_measurements_ - z_pred;

  x_ += K*Innovation;
  x_(3) = tools.wrapPi(x_(3));

  P_ += -K*S*K.transpose();

  VectorXd NIS_radar_vector = Innovation.transpose()*S_inv*Innovation;

  NIS_radar_ = NIS_radar_vector(0);
  // cout << x_;

}