#ifndef COSTFUNCTION_H_
#define COSTFUNCTION_H_


#include "MapPoint.h"
#include "Camera.h"
#include <iostream>

class CostFunction
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    CostFunction(Camera* camera, MapPoint* mappoint, double fx, double fy, double cx, double cy, const Eigen::Vector2d& ob_z);
    ~CostFunction();

    void ComputeJacobian(Eigen::Matrix<double, 2, 6>& J_pose, Eigen::Matrix<double, 2, 3>& J_position);
    void ComputeResiduals(Eigen::Vector2d& e, Eigen::Matrix2d& weight_info, double& weight_e2);

    MapPoint* p_mappoint_;
	Camera* p_camera_;
private:

    double fx_, fy_, cx_, cy_;
    Eigen::Vector2d ob_z_;
    double huber_delta_;  
    Eigen::Matrix2d info_matrix_;
    
};




#endif // !COSTFUNCTION_H_
