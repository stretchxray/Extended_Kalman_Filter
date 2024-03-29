#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
}

void KalmanFilter::Predict() {
    /**
     TODO:
     * predict the state
     */
    
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    /**
     TODO:
     * update the state by using Kalman Filter equations
     */
    
    VectorXd y = z - H_ * x_;
    
    MeasurementUpdate(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    /**
     TODO:
     * update the state by using Extended Kalman Filter equations
     */
    
    VectorXd h_proj = VectorXd(3);
    h_proj << 0,0,0;
    
    if (!(x_[0] == 0 && x_[1] == 0)) {
        double rho = sqrt(x_[0] * x_[0] + x_[1] * x_[1]);
        double phi = atan2(x_[1],x_[0]);
        double rho_dot = (x_[0]*x_[2] + x_[1]*x_[3]) / rho;
        
        h_proj << rho, phi, rho_dot;
    }
    
    VectorXd y = z - h_proj;
    
    if(y[1] < -M_PI){
        y[1] = M_PI + fmod(y[1] + M_PI, 2 * M_PI);
    }else{
        y[1] = fmod(y[1] + M_PI, 2 * M_PI) - M_PI;
    }
    
    MeasurementUpdate(y);
}

void KalmanFilter::MeasurementUpdate(const VectorXd &y) {
    
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K = P_ * Ht * Si;
    
    x_ = x_ + (K * y);
    MatrixXd I = MatrixXd::Identity(x_.size(),x_.size());
    
    P_ = (I - K * H_) * P_;
}
