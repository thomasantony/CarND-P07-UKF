#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

using std::cout;
using std::cerr;
using std::endl;

inline VectorXd CTRVModel(double delta_t, const VectorXd &x)
{
  VectorXd x_pred(5);
  //extract values for better readability
  double p_x = x(0);
  double p_y = x(1);
  double v = x(2);
  double yaw = x(3);
  double yawd = x(4);
  double nu_a = x(5);
  double nu_yawdd = x(6);

  //predicted state values
  double px_p, py_p;

  //avoid division by zero
  if (fabs(yawd) > 0.001) {
    px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
    py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
  }
  else {
    px_p = p_x + v*delta_t*cos(yaw);
    py_p = p_y + v*delta_t*sin(yaw);
  }

  double v_p = v;
  double yaw_p = yaw + yawd*delta_t;
  double yawd_p = yawd;

  //add noise
  px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
  py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
  v_p = v_p + nu_a*delta_t;

  yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
  yawd_p = yawd_p + nu_yawdd*delta_t;

  //write predicted sigma point into right column
  x_pred << px_p, py_p, v_p, yaw_p, yawd_p;
  return x_pred;
}

inline double angle_trim(double ang)
{
  return atan2(sin(ang), cos(ang));
}
inline VectorXd CTRV_postprocess(const VectorXd& x)
{
  VectorXd x_out = x;
  // Performs angle normalization for any state vector for this model
  x_out(3) = angle_trim(x(3));
  return x_out;
}
inline MatrixXd LidarMeasurement(const VectorXd& x){
  VectorXd z(2);
  z << x(0), x(1);
  return z;
}
inline VectorXd Lidar_postprocess(const VectorXd& z)
{
  return z;
}

inline VectorXd RadarMeasurement(const VectorXd& x)
{
  VectorXd z(3);
  // extract values for better readability
  double p_x = x(0);
  double p_y = x(1);
  double v   = x(2);
  double yaw = x(3);

  double v1 = cos(yaw)*v;
  double v2 = sin(yaw)*v;

  // measurement model
  z(0) = sqrt(p_x*p_x + p_y*p_y);                        //r
  z(1) = atan2(p_y,p_x);                                 //phi
  z(2) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot

  return z;
}
inline VectorXd Radar_postprocess(const VectorXd& z)
{
  VectorXd z_out = z;
  z_out(1) = angle_trim(z(1));
  return z_out;
}



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
  std_a_ = 9;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.3;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  n_x_ = 5;
  n_aug_ = n_x_ + 2;
  lambda_ = 3 - n_aug_;

  weights_ = VectorXd(2*n_aug_ + 1);
  weights_.fill(0.5/(lambda_ + n_aug_));
  weights_(0) = lambda_/(lambda_ + n_aug_);

  is_initialized_ = false;

  P_ << 500, 0, 0, 0, 0,
        0, 500, 0, 0, 0,
        0, 0, 1000,0, 0,
        0, 0, 0,  10, 0,
        0, 0, 0,  0, 10;
}


UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {


  if(!is_initialized_)
  {
    double p_x, p_y, rho, phi; //, v, yaw, yawd;
    if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
    {
      p_x = meas_package.raw_measurements_[0];
      p_y = meas_package.raw_measurements_[0];
      rho = sqrt(p_x*p_x + p_y*p_y);
      x_ << p_x, p_y, 0, 0, 0;
    } else if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
    {
      rho = meas_package.raw_measurements_[0];
      phi = meas_package.raw_measurements_[1];

      x_ << rho*cos(phi), rho*sin(phi), 0, 0, 0;
    }else{
      cerr << "Invalid measurement type!! " << endl;
    }
    if(std::fabs(rho) < 0.001)
    {
      cout<<"Measurement too small -- ignoring." <<endl;
      return;
    }
    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    return;
  }


  double delta_t = (meas_package.timestamp_ - time_us_)/ static_cast<double>(1e+6);

  // Predict state
  Prediction(delta_t);

  // Update
  if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
  {
    UKF::UpdateLidar(meas_package);
  }else if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
  {
    UKF::UpdateRadar(meas_package);
  }

  // Update stored timestamp
  time_us_ = meas_package.timestamp_;
}

MatrixXd UKF::GenerateAugmentedSigmaPoints(){

  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //calculate square root of P
  MatrixXd A = P_aug.llt().matrixL();

  //set first column of sigma point matrix
  Xsig_aug.col(0)  = x_aug;

  double sqrtpk = sqrt(lambda_+n_aug_);
  //set remaining sigma points
  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1)     = x_aug + sqrtpk * A.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrtpk * A.col(i);
  }
  return Xsig_aug;
}


MatrixXd UKF::PredictSigmaPoints(double delta_t, const MatrixXd& Xsig_aug)
{
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2*n_aug_ + 1);
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    Xsig_pred.col(i) = CTRVModel(delta_t, Xsig_aug.col(i));
  }
  return Xsig_pred;
}
/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  //create sigma point matrix
  MatrixXd Xsig_aug = GenerateAugmentedSigmaPoints();
  Xsig_pred_ = PredictSigmaPoints(delta_t, Xsig_aug);

  x_ = Xsig_pred_ * weights_;

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = CTRV_postprocess(Xsig_pred_.col(i) - x_);

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
//
//  cout << "x_pred : \n" << x_ << endl;
//  cout << "P_pred : \n" << P_ << endl;
}


/* Kalman update functions */
inline double UKF::UpdateUKF(int n_z, const VectorXd& z, const VectorXd& z_pred, const MatrixXd& Zsig, const MatrixXd& S){
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  //calculate cross correlation matrix
  //predicted state covariance matrix
  Tc.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = CTRV_postprocess(Xsig_pred_.col(i) - x_);

    VectorXd z_diff = Radar_postprocess(Zsig.col(i) - z_pred);

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose() ;
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = Radar_postprocess(z - z_pred);

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  double NIS = z_diff.transpose() * S.inverse() * z_diff;
  return NIS;
}


inline MatrixXd UKF::PredictLidarMeasurement(int n_z) {
  MatrixXd Zsig(n_z, 2 * n_aug_ + 1);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
    Zsig.col(i) = LidarMeasurement(Xsig_pred_.col(i));
  }
  return Zsig;
}
/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  int n_z = 2;
  //transform sigma points into measurement space
  MatrixXd Zsig = PredictLidarMeasurement(n_z);
  VectorXd z_pred = Zsig * weights_;

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Lidar_postprocess(Zsig.col(i) - z_pred);

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;

  S = S + R;

  NIS_laser_ = UpdateUKF(n_z, meas_package.raw_measurements_, z_pred, Zsig, S);

//  cout << "After Laser : \n x : \n" << x_ <<endl;
//  cout << "P : \n" << P_ << endl;

  cout << "NIS_laser : "<< NIS_laser_ << endl;
}


// Radar functions follow
inline MatrixXd UKF::PredictRadarMeasurement(int n_z) {
  MatrixXd Zsig(n_z, 2 * n_aug_ + 1);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
    Zsig.col(i) = RadarMeasurement(Xsig_pred_.col(i));
  }
  return Zsig;
}
/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  int n_z = 3;

  //transform sigma points into measurement space
  MatrixXd Zsig = PredictRadarMeasurement(n_z);
  VectorXd z_pred = Zsig * weights_;

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Radar_postprocess(Zsig.col(i) - z_pred);

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
      0, std_radphi_*std_radphi_, 0,
      0, 0,std_radrd_*std_radrd_;

  S = S + R;

  NIS_radar_ = UpdateUKF(n_z, meas_package.raw_measurements_, z_pred, Zsig, S);
//  cout << "After Radar : \n x : \n" << x_ <<endl;
//  cout << "P : \n" << P_ << endl;

  cout << "NIS_radar : "<< NIS_radar_ << endl;
}
