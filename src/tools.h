#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"


inline double angle_normalize(double ang)
{
  return atan2(sin(ang), cos(ang));
}
inline double clamp(double val, double min_val, double max_val)
{
  return std::max(min_val, std::min(val, max_val));
}

class Tools {
public:
  /**
  * A helper method to calculate RMSE.
  */
  static Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

};

#endif /* TOOLS_H_ */
