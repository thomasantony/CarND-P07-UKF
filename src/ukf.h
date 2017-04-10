#ifndef UKF_H
#define UKF_H

#include <unordered_map>
#include "Eigen/Dense"

enum SensorType {
  LASER,
  RADAR
};

#include "measurement_package.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

using ProcessNoise = MatrixXd;
using DynamicModelFunc = std::function<VectorXd(double, const VectorXd&)>;
using PostProcessor = std::function<VectorXd(const VectorXd&)>;
// Dynamic model includes number of states, dynamic model function,
// process noise matrix and post processor function.
using DynamicModel = std::tuple<DynamicModelFunc, ProcessNoise, PostProcessor>;

using SensorNoise = MatrixXd;
using SensorFunc = std::function<VectorXd (const VectorXd &x)>;
using SensorModel = std::tuple<SensorFunc, SensorNoise, PostProcessor>;

using SensorMap = std::unordered_map<SensorType, SensorModel>;

using Initializer = std::function<bool(MeasurementPackage, VectorXd&, MatrixXd&)>;
class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* dynamic model used for state propagation
  DynamicModelFunc statePredictor_;

  ///* post processor called after any state operation
  PostProcessor statePostProcessor_;

  ///* initializer function for boot-strapping the filter
  Initializer initializer_fn_;

  //* Map of sensor types and sensors to be used by the filter
  SensorMap sensors_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* process noise covariance matrix
  MatrixXd Q_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* augmented state vector
  VectorXd x_aug_;

  ///* augmented covariance matrix
  MatrixXd P_aug_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* the current NIS for all sensors
  std::unordered_map<SensorType, double> NIS_;

  // Last measurement
  MeasurementPackage last_measurement_;
  /**
   * Constructor
   */
  UKF(const DynamicModel& model, const Initializer& initializer);

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * Adds a new sensor model to the UKF
   *
   * @param type Type of sensor
   * @param sensor Sensor model definition
   */
  void AddSensor(SensorType type, const SensorModel& sensor);

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using sensor models
   * @param meas_package The measurement at k+1
   */
  void UpdateUKF(MeasurementPackage meas_package);
private:
  /**
  * Initializes the UKF from a single measurement by calling pre-defined initializer
  *
  * @param meas_package
  */
  bool InitializeFilter(MeasurementPackage meas_package);

  inline MatrixXd GenerateAugmentedSigmaPoints();

};

#endif /* UKF_H */
