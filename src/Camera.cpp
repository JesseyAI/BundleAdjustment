#include "../include/Camera.h"

Camera::Camera(const Sophus::SE3d& pose, int id, bool isfixed):
  pose_(pose), id_(id), isfixed_(isfixed), state_index_(-1)
{

}

Camera::~Camera()
{

}

int Camera::GetCamId()
{
    return id_;
}

void Camera::UpdateCamPose(const Sophus::SE3d& delta_pose)
{
    pose_ = delta_pose * pose_;
}

void Camera::UpdateCamPose(const Eigen::Matrix<double, 6, 1>& delta_pose)
{
    Sophus::SE3d se3_delta = Sophus::SE3d::exp(delta_pose);
    pose_ = se3_delta * pose_;
}

void Camera::SetFixed()
{
	isfixed_ = true;
} // setFixed

bool Camera::IsFixed()
{
	return isfixed_;
} // isFixed

void Camera::SetPose (const Sophus::SE3d& pose)
{
	pose_ = pose;
}

Sophus::SE3d Camera::GetPose()
{
  return pose_;
}

