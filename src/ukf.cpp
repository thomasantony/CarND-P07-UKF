#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <unordered_map>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

using std::cout;
using std::cerr;
using std::endl;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF(const DynamicModel& model, const SensorMap &sensors, const Initializer& initializer) {

  initializer_fn_ = initializer;
  sensors_ = sensors;

  // Store the dynamic model
  std::tie(statePredictor_, Q_, statePostProcessor_) = model;


  n_x_ = 5;
  n_aug_ = n_x_ + 2;
  lambda_ = 3 - n_aug_;

  weights_ = VectorXd(2*n_aug_ + 1);
  weights_.fill(0.5/(lambda_ + n_aug_));
  weights_(0) = lambda_/(lambda_ + n_aug_);

  is_initialized_ = false;

  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);
  Xsig_pred_.fill(0.0);
}


UKF::~UKF() {}

void UKF::AddSensor(SensorType type, const SensorModel &sensor) {
  sensors_[type] = sensor;
}


bool UKF::InitializeFilter(MeasurementPackage meas_package)
{
  is_initialized_ = initializer_fn_(meas_package, x_, P_);
  if (!is_initialized_)
  {
    return false;
  }

  n_x_ = x_.size();
  n_aug_ = x_.size() + Q_.rows();

  lambda_ = 3 - n_aug_;

  weights_ = VectorXd(2*n_aug_ + 1);
  weights_.fill(0.5/(lambda_ + n_aug_));
  weights_(0) = lambda_/(lambda_ + n_aug_);

  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  x_aug_ = VectorXd(n_aug_);

  //create augmented state covariance
  P_aug_ = MatrixXd(n_aug_, n_aug_);

  time_us_ = meas_package.timestamp_;
  return is_initialized_;

}
/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  // If measurements are too tiny, make them slightly bigger
  if(meas_package.raw_measurements_.squaredNorm() < 0.00001)
  {
    meas_package.raw_measurements_.fill(0.01);
  }
  if(!is_initialized_)
  {
    last_measurement_ = meas_package;
    if(!InitializeFilter(meas_package))
      return;
  }
  // Predict state
  double delta_t = (meas_package.timestamp_ - time_us_)/ static_cast<double>(1e+6);
  try {
    Prediction(delta_t);
  } catch (std::range_error e) {
    // If convariance matrix is non positive definite (because of numerical instability),
    // restart the filter using previous measurement as initialiser.
    InitializeFilter(last_measurement_);
    // Redo prediction using the current measurement
    // We don't get exception this time, because initial P (identity) is positive definite.
    Prediction(delta_t);
  }

  // Update UKF with measurement
  UpdateUKF(meas_package);

  // Update stored timestamp
  time_us_ = meas_package.timestamp_;
  last_measurement_ = meas_package;
}

MatrixXd UKF::GenerateAugmentedSigmaPoints(){


  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  x_aug_.head(n_x_) = x_;
  x_aug_.tail(n_aug_-n_x_).fill(0.0);

  //create augmented covariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_, n_x_) = P_;
  P_aug_.bottomRightCorner(n_aug_-n_x_, n_aug_-n_x_) = Q_;

  // Numerical stability fix by Olli Vertanen
  // Source: https://discussions.udacity.com/t/numerical-instability-of-the-implementation/230449/2?u=tantony
  //calculate square root of P
  Eigen::LLT<MatrixXd> lltOfPaug(P_aug_);
  if (lltOfPaug.info() == Eigen::NumericalIssue) {
    // if decomposition fails, we have numerical issues
    // std::cout << "LLT failed!" << std::endl;
    //Eigen::EigenSolver<MatrixXd> es(P_aug);
    //cout << "Eigenvalues of P_aug:" << endl << es.eigenvalues() << endl;
    throw std::range_error("LLT failed");
  }
  // 2. get the lower triangle
  MatrixXd A = lltOfPaug.matrixL();

  //set first column of sigma point matrix
  Xsig_aug.col(0)  = x_aug_;

  double sqrtpk = sqrt(lambda_+n_aug_);
  //set remaining sigma points
  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1)     = x_aug_ + sqrtpk * A.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug_ - sqrtpk * A.col(i);
  }
  return Xsig_aug;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  //create sigma point matrix
  MatrixXd Xsig_aug = GenerateAugmentedSigmaPoints();

  // Predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    Xsig_pred_.col(i) = statePostProcessor_(statePredictor_(delta_t, Xsig_aug.col(i)));
  }
  x_ = Xsig_pred_ * weights_;

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 1; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = statePostProcessor_(Xsig_pred_.col(i) - Xsig_pred_.col(0));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}



/* Kalman update function */
inline void UKF::UpdateUKF(MeasurementPackage meas_package){
  SensorFunc sensorFunc;
  SensorNoise sensorNoise;
  PostProcessor sensorPostProcessor;

  // Sensor type not found. Ignore the measurement
  if(sensors_.count(meas_package.sensor_type_) == 0)
  {
    return;
  }
  std::tie(sensorFunc, sensorNoise, sensorPostProcessor) = sensors_[meas_package.sensor_type_];
  auto n_z = sensorNoise.rows();
  const VectorXd& z = meas_package.raw_measurements_;

  //transform sigma points into measurement space
  MatrixXd Zsig(n_z, 2 * n_aug_ + 1);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points
    Zsig.col(i) = sensorFunc(Xsig_pred_.col(i));
  }
  VectorXd z_pred = Zsig * weights_;

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = sensorPostProcessor(Zsig.col(i) - z_pred);

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  // add measurement covariance
  S = S + sensorNoise;

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  //calculate cross correlation matrix
  Tc.fill(0.0);

  for (int i = 1; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    // state difference
    VectorXd x_diff = statePostProcessor_(Xsig_pred_.col(i) - Xsig_pred_.col(0));

    VectorXd z_diff = sensorPostProcessor(Zsig.col(i) - z_pred);

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose() ;
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = sensorPostProcessor(z - z_pred);

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  // Save NIS
  NIS_[meas_package.sensor_type_] = z_diff.transpose() * S.inverse() * z_diff;
}

