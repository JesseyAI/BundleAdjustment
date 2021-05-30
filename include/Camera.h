#ifndef CAMERA_H_
#define CAMERA_H_

#include <Eigen/Core>
#include <sophus/se3.hpp>

class Camera
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Camera(const Sophus::SE3d& pose, int id, bool isfixed = false);
    ~Camera();

    void SetPose(const Sophus::SE3d& pose);
    void UpdateCamPose(const Sophus::SE3d& delta_pose);
    void UpdateCamPose(const Eigen::Matrix<double, 6, 1>& delta_pose);
    void SetFixed();
    bool IsFixed();
    Sophus::SE3d GetPose();

    int GetCamId();

    int state_index_;
private:
    
    Sophus::SE3d pose_;
    int id_;
    bool isfixed_;
};


#endif // !CAMERA_H_
