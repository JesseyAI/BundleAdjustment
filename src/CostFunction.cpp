#include "../include/CostFunction.h"


CostFunction::CostFunction(Camera* camera, MapPoint* mappoint, double fx, double fy, double cx, double cy, const Eigen::Vector2d& ob_z)
    :p_camera_(camera), p_mappoint_(mappoint), fx_(fx), fy_(fy), cx_(cx), cy_(cy), ob_z_(ob_z),huber_delta_(1.0), info_matrix_(Eigen::Matrix2d::Identity())
{

}

CostFunction::~CostFunction()
{
}

void CostFunction::ComputeResiduals(Eigen::Vector2d& e, Eigen::Matrix2d& weight_info, double& weight_e2)
{
    Eigen::Vector3d Pc = p_camera_->GetPose()*p_mappoint_->GetPosition();   //这里优化用的位姿是 Tcw

    // std::cout << p_camera_->GetPose().rotationMatrix() << std::endl;
    // std::cout <<"translation \n" <<p_camera_->GetPose().translation() << std::endl;

    // std::cout <<"位置\n" <<p_mappoint_->GetPosition() << std::endl;

	const Eigen::Vector2d Puv(
		fx_ * Pc.x() / Pc.z() + cx_,
		fy_ * Pc.y() / Pc.z() + cy_
	);

    //观测减预测
    e = ob_z_ - Puv;

    double et_info_e = e.transpose() * info_matrix_ * e;

    //compute huber weight
    double weight = 1.0;
    double sqrt_e2 = sqrt(et_info_e);
    if(sqrt_e2 > huber_delta_)
    {
        weight = 2 * huber_delta_ * sqrt_e2 -huber_delta_ * huber_delta_;
        weight = weight / et_info_e;
    }
    
    // compute huber weighted info matrix.
	weight_info = weight * info_matrix_;
	
	// compute weighted et_weighte_info_e
	weight_e2 = weight * et_info_e;

    //std::cout << "weighted_e2:\n" << weight_e2 << std::endl;

	//std::cout << "weighted_info \n" << weight_info << std::endl;

}

void CostFunction::ComputeJacobian(Eigen::Matrix<double, 2, 6>& J_pose, Eigen::Matrix<double, 2, 3>& J_position)
{
    Eigen::Vector3d Pc = p_camera_->GetPose()*p_mappoint_->GetPosition();	

    const double X = Pc.x();
    const double Y = Pc.y();
    const double Z = Pc.z();
    
    const double Z2 = Z * Z;

    J_pose.setZero();
    // Jobian matrix of camera pose
    J_pose(0,0) = fx_ / Z;  
    J_pose(0,2) = - fx_ * X / Z2;
    J_pose(0,3) = - fx_ * X * Y / Z2;
    J_pose(0,4) = fx_ + fx_ * X * X / Z2;
    J_pose(0,5) = - fx_ * Y / Z;

    J_pose(1,1) = fy_ / Z;
    J_pose(1,2) = - fy_ * Y / Z2;
    J_pose(1,3) = - fy_ - fy_ * Y * Y / Z2;
    J_pose(1,4) = fy_ * X * Y / Z2;
    J_pose(1,5) = fy_ * X / Z;

    J_pose = - J_pose;

    // Jobian matrix of 3d points
    Eigen::Matrix<double, 2, 3> temp;
    temp.setZero();
    temp(0,0) = fx_ / Z;
    temp(0,2) = -fx_ * X / Z2;
    temp(1,1) = fy_ / Z;
    temp(1,2) = -fy_ * Y / Z2;

    J_position = -temp * p_camera_->GetPose().rotationMatrix();

}