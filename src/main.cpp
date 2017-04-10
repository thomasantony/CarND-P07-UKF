
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include "Eigen/Dense"
#include "ukf.h"
#include "tools.h"
#include "ground_truth_package.h"
#include "measurement_package.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

namespace{
  /* Function that propagates the state based on the CTRV model */
  VectorXd CTRV_ModelFunc(double delta_t, const VectorXd &x)
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
    if (std::fabs(yawd) > 0.01) {
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
  // Process noise standard deviation longitudinal acceleration in m/s^2
  const auto std_a = 3.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  const auto std_yawdd = 0.55;

  const MatrixXd Process_Noise = (MatrixXd(2,2) << std_a*std_a, 0,
                                           0, std_yawdd*std_yawdd).finished();
  inline VectorXd CTRV_Postprocessor(const VectorXd& x)
  {
    auto max_yawrate = 60*M_PI/180; // rad/s, bound yawrate to a reasonable number

    VectorXd x_out = x;

    x_out(2) = clamp(x(2), -25, 25);
    x_out(3) = angle_normalize(x(3));  // Normalize angle to between -pi and +pi
    x_out(4) = clamp(x(4), -max_yawrate, max_yawrate);
    return x_out;
  }
  const auto CTRV_model = DynamicModel(CTRV_ModelFunc, Process_Noise, CTRV_Postprocessor);

  /***************************************************************************/
  /*                            Sensor Models                                */
  /***************************************************************************/

  /**********************/
  /* Lidar Sensor model */
  /**********************/
  inline MatrixXd Lidar_Measurement(const VectorXd& x){
    VectorXd z(2);
    z << x(0), x(1);
    return z;
  }
  // Laser measurement noise standard deviation position1 in m
  const auto std_laspx = 0.15;

  // Laser measurement noise standard deviation position2 in m
  const auto std_laspy = 0.15;

  const MatrixXd Lidar_Noise = (MatrixXd(2,2) << std_laspx*std_laspx, 0,
                                                  0, std_laspy*std_laspy).finished();
  inline VectorXd Lidar_Postprocessor(const VectorXd &z)
  {
    return z; // No post-processing
  }
  const auto Lidar_Sensor = SensorModel(Lidar_Measurement, Lidar_Noise, Lidar_Postprocessor);

  /**********************/
  /* Radar Sensor model */
  /**********************/
  inline VectorXd Radar_Measurement(const VectorXd& x)
  {
    VectorXd z(3);
    // extract values for better readability
    double p_x = x(0);
    double p_y = x(1);
    double v   = x(2);
    double yaw = x(3);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    double rho = sqrt(p_x*p_x + p_y*p_y);
    if(rho < 0.01)
    {
      rho = 0.01;
    }
    // measurement model
    z(0) = rho;                        //r
    z(1) = atan2(p_y,p_x);             //phi
    z(2) = (p_x*v1 + p_y*v2 ) / rho;   //r_dot

    return z;
  }
  // Radar measurement noise standard deviation radius in m
  const auto std_radr = 0.3;
  // Radar measurement noise standard deviation angle in rad
  const auto std_radphi = 0.03;
  // Radar measurement noise standard deviation radius change in m/s
  const auto std_radrd = 0.3;
  const MatrixXd Radar_Noise = (MatrixXd(3,3) << std_radr*std_radr, 0, 0,
                                                0, std_radphi*std_radphi, 0,
                                                0, 0,std_radrd*std_radrd).finished();

  inline VectorXd Radar_Postprocessor(const VectorXd &z)
  {
    VectorXd z_out = z;
    z_out(1) = angle_normalize(z(1));
    return z_out;
  }
  const auto Radar_Sensor = SensorModel(Radar_Measurement, Radar_Noise, Radar_Postprocessor);

  /* UKF Initializer */
  bool InitUKF(MeasurementPackage first_measurement, VectorXd& x_in, MatrixXd& P_in)
  {
    // initial state
    x_in = VectorXd(5);

    // initial covariance matrix
    P_in = MatrixXd(5,5);
    double p_x, p_y, rho, rhodot, phi; //, v, yaw, yawd;
    if(first_measurement.sensor_type_ == SensorType::LASER)
    {
      p_x = first_measurement.raw_measurements_[0];
      p_y = first_measurement.raw_measurements_[1];
      rho = sqrt(p_x*p_x + p_y*p_y);
      x_in << p_x, p_y, 0, 0, 0;
    } else if(first_measurement.sensor_type_ == SensorType::RADAR)
    {
      rho = first_measurement.raw_measurements_[0];
      phi = first_measurement.raw_measurements_[1];
      rhodot = first_measurement.raw_measurements_[2];
      x_in << rho*cos(phi), rho*sin(phi), rhodot, 0, 0;
    }
    P_in << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

    return true;
  }

}
void check_arguments(int argc, char* argv[]) {
  string usage_instructions = "Usage instructions: ";
  usage_instructions += argv[0];
  usage_instructions += " path/to/input.txt output.txt";

  bool has_valid_args = false;

  // make sure the user has provided input and output files
  if (argc == 1) {
    cerr << usage_instructions << endl;
  } else if (argc == 2) {
    cerr << "Please include an output file.\n" << usage_instructions << endl;
  } else if (argc == 3) {
    has_valid_args = true;
  } else if (argc > 3) {
    cerr << "Too many arguments.\n" << usage_instructions << endl;
  }

  if (!has_valid_args) {
    exit(EXIT_FAILURE);
  }
}

void check_files(ifstream& in_file, string& in_name,
                 ofstream& out_file, string& out_name) {
  if (!in_file.is_open()) {
    cerr << "Cannot open input file: " << in_name << endl;
    exit(EXIT_FAILURE);
  }

  if (!out_file.is_open()) {
    cerr << "Cannot open output file: " << out_name << endl;
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char* argv[]) {

  check_arguments(argc, argv);

  string in_file_name_ = argv[1];
  ifstream in_file_(in_file_name_.c_str(), ifstream::in);

  string out_file_name_ = argv[2];
  ofstream out_file_(out_file_name_.c_str(), ofstream::out);

  check_files(in_file_, in_file_name_, out_file_, out_file_name_);

  /**********************************************
   *  Set Measurements                          *
   **********************************************/

  vector<MeasurementPackage> measurement_pack_list;
  vector<GroundTruthPackage> gt_pack_list;

  string line;

  // prep the measurement packages (each line represents a measurement at a
  // timestamp)
  while (getline(in_file_, line)) {
    string sensor_type;
    MeasurementPackage meas_package;
    GroundTruthPackage gt_package;
    istringstream iss(line);
    long long timestamp;

    // reads first element from the current line
    iss >> sensor_type;

    if (sensor_type.compare("L") == 0) {
      // laser measurement

      // read measurements at this timestamp
      meas_package.sensor_type_ = SensorType::LASER;
      meas_package.raw_measurements_ = VectorXd(2);
      float px;
      float py;
      iss >> px;
      iss >> py;
      meas_package.raw_measurements_ << px, py;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    } else if (sensor_type.compare("R") == 0) {
      // radar measurement

      // read measurements at this timestamp
      meas_package.sensor_type_ = SensorType::RADAR;
      meas_package.raw_measurements_ = VectorXd(3);
      float ro;
      float phi;
      float ro_dot;
      iss >> ro;
      iss >> phi;
      iss >> ro_dot;
      meas_package.raw_measurements_ << ro, phi, ro_dot;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    }

      // read ground truth data to compare later
      float x_gt;
      float y_gt;
      float vx_gt;
      float vy_gt;
      iss >> x_gt;
      iss >> y_gt;
      iss >> vx_gt;
      iss >> vy_gt;
      gt_package.gt_values_ = VectorXd(4);
      gt_package.gt_values_ << x_gt, y_gt, vx_gt, vy_gt;
      gt_pack_list.push_back(gt_package);
  }

  // Create a UKF instance and pass in the dynamic model and initializer
  UKF ukf(CTRV_model, InitUKF);

  // Add sensor models
  ukf.AddSensor(SensorType::LASER, Lidar_Sensor);
  ukf.AddSensor(SensorType::RADAR, Radar_Sensor);

  // used to compute the RMSE later
  vector<VectorXd> estimations;
  vector<VectorXd> ground_truth;

  // start filtering from the second frame (the speed is unknown in the first
  // frame)

  size_t number_of_measurements = measurement_pack_list.size();

  // column names for output file
  out_file_ << "px" << "\t";
  out_file_ << "py" << "\t";
  out_file_ << "v" << "\t";
  out_file_ << "yaw_angle" << "\t";
  out_file_ << "yaw_rate" << "\t";
  out_file_ << "px_measured" << "\t";
  out_file_ << "py_measured" << "\t";
  out_file_ << "px_true" << "\t";
  out_file_ << "py_true" << "\t";
  out_file_ << "vx_true" << "\t";
  out_file_ << "vy_true" << "\t";
  out_file_ << "NIS" << "\n";

  for (size_t k = 0; k < number_of_measurements; ++k) {
    // Call the UKF-based fusion
    ukf.ProcessMeasurement(measurement_pack_list[k]);
    // output the estimation
    out_file_ << ukf.x_(0) << "\t"; // pos1 - est
    out_file_ << ukf.x_(1) << "\t"; // pos2 - est
    out_file_ << ukf.x_(2) << "\t"; // vel_abs -est
    out_file_ << ukf.x_(3) << "\t"; // yaw_angle -est
    out_file_ << ukf.x_(4) << "\t"; // yaw_rate -est

    // output the measurements
    if (measurement_pack_list[k].sensor_type_ == SensorType::LASER) {
      // output the estimation

      // p1 - meas
      out_file_ << measurement_pack_list[k].raw_measurements_(0) << "\t";

      // p2 - meas
      out_file_ << measurement_pack_list[k].raw_measurements_(1) << "\t";
    } else if (measurement_pack_list[k].sensor_type_ == SensorType::RADAR) {
      // output the estimation in the cartesian coordinates
      float ro = measurement_pack_list[k].raw_measurements_(0);
      float phi = measurement_pack_list[k].raw_measurements_(1);
      out_file_ << ro * cos(phi) << "\t"; // p1_meas
      out_file_ << ro * sin(phi) << "\t"; // p2_meas
    }

    // output the ground truth packages
    out_file_ << gt_pack_list[k].gt_values_(0) << "\t";
    out_file_ << gt_pack_list[k].gt_values_(1) << "\t";
    out_file_ << gt_pack_list[k].gt_values_(2) << "\t";
    out_file_ << gt_pack_list[k].gt_values_(3) << "\t";

    // output the NIS values
    out_file_ << ukf.NIS_[measurement_pack_list[k].sensor_type_] << "\n";

    // convert ukf x vector to cartesian to compare to ground truth
    VectorXd ukf_x_cartesian_ = VectorXd(4);

    float x_estimate_ = ukf.x_(0);
    float y_estimate_ = ukf.x_(1);
    float vx_estimate_ = ukf.x_(2) * cos(ukf.x_(3));
    float vy_estimate_ = ukf.x_(2) * sin(ukf.x_(3));
    
    ukf_x_cartesian_ << x_estimate_, y_estimate_, vx_estimate_, vy_estimate_;
    
    estimations.push_back(ukf_x_cartesian_);
    ground_truth.push_back(gt_pack_list[k].gt_values_);

  }

  // compute the accuracy (RMSE)
  cout << "Accuracy - RMSE:" << endl << Tools::CalculateRMSE(estimations, ground_truth) << endl;

  // close files
  if (out_file_.is_open()) {
    out_file_.close();
  }

  if (in_file_.is_open()) {
    in_file_.close();
  }

  cout << "Done!" << endl;
  return 0;
}
