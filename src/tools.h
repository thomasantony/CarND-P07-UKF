#ifndef TOOLS_H_
#define TOOLS_H_
#include <iostream>
#include <vector>
#include "Eigen/Dense"

/**
 * Normalizes an angle to be between -pi and +pi
 * @param ang Angle to be normalized
 * @return Normalized angle
 */
inline double angle_normalize(double ang)
{
  return atan2(sin(ang), cos(ang));
}
/**
 * Clamps a value between a maximum and a minimum bound.
 *
 * @param val
 * @param min_val
 * @param max_val
 * @return
 */
template <typename T>
inline T clamp(T val, T min_val, T max_val)
{
  return std::max(min_val, std::min(val, max_val));
}

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cerr;
using std::endl;

inline VectorXd CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4), dx(4);
  rmse.fill(0.0);

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  if(estimations.size() == 0)
  {
    cerr << "Need at least one estimate, boss!" << endl;
    return rmse;
  }else if(estimations.size() != ground_truth.size())
  {
    cerr << "estimations.size() != ground_truth.size()" <<endl;
    return rmse;
  }
  //  * the estimation vector size should equal ground truth vector size
  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
    dx = (estimations[i] - ground_truth[i]);
    rmse += dx.cwiseAbs2();
  }

  //calculate the mean
  rmse /= estimations.size();

  //calculate the squared root
  rmse = rmse.cwiseSqrt();

  //return the result
  return rmse;
}

#endif /* TOOLS_H_ */
