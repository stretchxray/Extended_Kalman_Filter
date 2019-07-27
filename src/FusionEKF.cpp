#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;
    
    previous_timestamp_ = 0;
    
    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);
    
    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
    0, 0.0225;
    
    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
    0, 0.0009, 0,
    0, 0, 0.09;
    
    /**
     TODO:
     * Finish initializing the FusionEKF.
     * Set the process and measurement noises
     */
    
    H_laser_ << 1,0,0,0,
    0,1,0,0;
    
    //initial state transition matrix F
    ekf_.F_ = Eigen::MatrixXd(4,4);
    ekf_.F_ << 1,0,1,0,
    0,1,0,1,
    0,0,1,0,
    0,0,0,1;
    
    //state covariance matrix P
    ekf_.P_ = Eigen::MatrixXd(4,4);
    ekf_.P_ << 1,0,0,0,
    0,1,0,0,
    0,0,1000,0,
    0,0,0,1000;
    
    //process covariance matrix Q
    ekf_.Q_ = Eigen::MatrixXd(4,4);
    
    //set acceleration noise 3^2
    noise_ax_ = 9;
    noise_ay_ = 9;
    
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    
    
    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
        /**
         TODO:
         * Initialize the state ekf_.x_ with the first measurement.
         * Create the covariance matrix.
         * Remember: you'll need to convert radar from polar to cartesian coordinates.
         */
        // first measurement
        cout << "EKF: " << endl;
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 1, 1;
        
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
             Convert radar from polar to cartesian coordinates and initialize state.
             */
            
            double rho = measurement_pack.raw_measurements_[0];
            double phi = measurement_pack.raw_measurements_[1];
            double rho_dot = measurement_pack.raw_measurements_[2];
            
            //convert polar to cartesian coords
            double px = rho * cos(phi);
            double py = rho * sin(phi);
            
            /*
             while we can perfectly calculate px and py from phi, we cannot compute vx and vy from phi. We will need
             yaw (which is introduced in UKF) to compute vx and vy. So even from radar measurement, we can only compute px and py.
             */
            double vx = rho_dot * cos(phi);
            double vy = rho_dot * sin(phi);
            
            ekf_.x_ << px,py,0,0;
            
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            /**
             Initialize state.
             */
            ekf_.x_ << measurement_pack.raw_measurements_[0],measurement_pack.raw_measurements_[1],0,0;
        }
        
        previous_timestamp_ = measurement_pack.timestamp_;
        
        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }
    
    /*****************************************************************************
     *  Prediction
     ****************************************************************************/
    
    /**
     TODO:
     * Update the state transition matrix F according to the new elapsed time.
     - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
     */
    
    // convert to secs
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;
    
    //add t to F
    ekf_.F_(0,2) = dt;
    ekf_.F_(1,3) = dt;
    
    //calculate t for Q
    float dt2 = dt * dt;
    float dt3 = dt2 * dt;
    float dt4 = dt3 * dt;
    float dt4_4 = dt4 / 4;
    float dt3_2 = dt3 / 2;
    
    
    ekf_.Q_ << dt4_4 * noise_ax_, 0, dt3_2 * noise_ax_,0,
    0,dt4_4 * noise_ay_,0,dt3_2 * noise_ay_,
    dt3_2 * noise_ax_,0,dt2 * noise_ax_,0,
    0,dt3_2 * noise_ay_,0,dt2 * noise_ay_;
    
    ekf_.Predict();
    
    /*****************************************************************************
     *  Update
     ****************************************************************************/
    
    /**
     TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
     */
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        Hj_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.H_ = Hj_;
        ekf_.R_ = R_radar_;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } else {
        // Laser updates
        ekf_.H_ = H_laser_;
        ekf_.R_ = R_laser_;
        
        ekf_.Update(measurement_pack.raw_measurements_);
        
    }
    
    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
