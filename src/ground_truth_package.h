#ifndef GROUND_TRUTH_PACKAGE_H_
#define GROUND_TRUTH_PACKAGE_H_

#include "ukf.h"
#include "Eigen/Dense"

class GroundTruthPackage {
public:
  long long timestamp_;

  SensorType sensor_type_;

  Eigen::VectorXd gt_values_;

};

#endif /* GROUND_TRUTH_PACKAGE_H_ */
