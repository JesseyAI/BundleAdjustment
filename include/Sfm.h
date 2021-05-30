//Including Fumdmental, camera pose, triangulate
#ifndef SFM_H_
#define SFM_H_

#include <iostream>
#include <Eigen/Core>
#include <memory>
#include <sophus/se3.hpp>
#include <opencv2/opencv.hpp>
#include <opencv2/calib3d.hpp>
#include <vector>
#include <unordered_map>
#include <opencv2/core/eigen.hpp>
#include <fstream>

struct sFrame
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    sFrame(Eigen::Matrix4d Twc): Twc_(Twc){}

    Eigen::Matrix4d Twc_;
    std::unordered_map<int, Eigen::Vector2d> per_feature_;
};

class SFM
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    typedef std::shared_ptr<SFM> Ptr;
    SFM();
    SFM(double fx, double fy,double cx,double cy);
    ~SFM();

    bool Initializer(std::vector<sFrame>& frames);

    std::vector<Eigen::Vector3d> Triangulate(const std::vector<sFrame>& v_sframes);
private:
    
    cv::Mat_<double> K_;
    double fx_,fy_, cx_, cy_;
    
    cv::Mat ComputeF(const std::vector<cv::Point2d>& vpts1, const std::vector<cv::Point2d>& vpts2);
    void DecomposeF(cv::Mat F, cv::Mat_<double> &R1, cv::Mat_<double> &R2, cv::Mat_<double> &t1, cv::Mat_<double> &t2);
    void RecoverPose(const std::unordered_map<int,Eigen::Vector2d>& umpts1, const std::unordered_map<int, Eigen::Vector2d>& umpts2, Eigen::Matrix3d& R, Eigen::Vector3d& t);
    double TestTriangulate(const std::vector<cv::Point2d>& vpts1, const std::vector<cv::Point2d>& vpts2, cv::Mat_<double> R, cv::Mat_<double> t);
    cv::Mat calEpipolar();

    
};

#endif // SFM_H_
