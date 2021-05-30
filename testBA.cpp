#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <unordered_map>
#include "include/Sfm.h"
#include "include/BundleAdjustment.h"
#include <sophus/se3.hpp>


using namespace std;
using namespace Eigen;

const int num_cameras = 4;
const int num_map_points = 100;
const double height = 1024;
const double width = 768;


const double fx = 750;
const double fy = 500;
const double cx = 250;
const double cy = 250;

void GenerateCameraData(int nCams);
void AddNoise(std::vector<Eigen::Vector3d>& mptpoints, std::vector<sFrame>& t_sframes, double mpt_noise, double cam_noise, double ob_noise);

void createMaPoints(const std::vector<sFrame>& v_frame, std::vector<Eigen::Vector3d>& gt_points)
{
    std::mt19937 gen(std::chrono::system_clock::now().time_since_epoch().count());
    
    std::uniform_real_distribution<double> d_x{0.0, 4.0};
    std::uniform_real_distribution<double> d_y{0.0, 5.0};
    std::uniform_real_distribution<double> d_z{0.1, 5.0};

    int cnt = 0;
    while(true)
    {
        Eigen::Vector3d pt;
        pt[0] = d_x(gen);
        pt[1] = d_y(gen);
        pt[2] = d_z(gen);

        int n_obervation_times = 0;

        for(int j = 0; j < num_cameras; j++)
        {
            Eigen::Matrix3d Rwc = v_frame[j].Twc_.block(0, 0, 3, 3);
            Eigen::Vector3d twc = v_frame[j].Twc_.block(0, 3, 3, 1);
            Eigen::Matrix3d Rcw = Rwc.transpose();
            Eigen::Vector3d tcw = -Rcw * twc;

            Eigen::Vector3d Pc = Rcw * pt + tcw;
            Eigen::Vector2d uv (
                fx*Pc[0]/Pc[2] + cx,
                fy*Pc[1]/Pc[2] + cy
            );

            if( Pc[2] < 0)
                continue;
            n_obervation_times ++;
        }

        if(n_obervation_times == num_cameras)
        {
            gt_points.push_back(pt);
            cnt ++;
        }
            
        if(cnt == num_map_points)
            break;
    }
}


void createCameraPose(std::vector<sFrame>& frame_pose, std::vector<Sophus::SE3d>& gt_pose, int n_cam)
{

    Eigen::Matrix4d T;
    T.setIdentity();

    Eigen::Matrix3d Rx, Ry, Rz;
    Eigen::Matrix3d R; 
    Eigen::Vector3d t;
    cv::RNG rng ( cv::getTickCount() );

    double radius = 8;
    for ( int i = 0; i < n_cam; i ++ ) {
        // Rotation.
        if(i == 0)
        {
            sFrame frame(T);
            frame_pose.push_back(frame);
            R = T.block(0,0,3,3);
            t = T.block(0,3,3,1);
            Sophus::SE3d cam ( R, t );
            gt_pose.push_back ( cam );
            //std::cout << "R \n" << T  << std::endl;
            continue;
        }

        double theta = i * 2 * M_PI / (6 * 4);

        Eigen::Matrix3d R;
        R = Eigen::AngleAxisd(theta, Eigen::Vector3d::UnitZ());
        Eigen::Vector3d t = Eigen::Vector3d(radius * cos(theta) - radius, radius * sin(theta), 1 * sin(2 * theta));

        T.block(0,0,3,3) = R;
        T.block(0,3,3,1) = t;
        // SE3
        sFrame frame(T);
        frame_pose.push_back(frame);

        Sophus::SE3d cam ( R, t );
        gt_pose.push_back ( cam );
		
    }
    
}

void detectFeatures(const Eigen::Matrix3d& K, const std::vector<Eigen::Vector3d> &mappoints, std::vector<sFrame> &v_frames, bool is_add_noise = true)
{
    std::mt19937 gen{12345};
    const double pixel_sigma = 2;
    std::normal_distribution<> d{0.0, pixel_sigma};

    for(size_t m = 0; m < mappoints.size(); m++)
    {
         Eigen::Vector3d Pw = mappoints[m];
        for(size_t nc = 0; nc < v_frames.size(); nc ++)
        {
            Eigen::Matrix3d Rwc = v_frames[nc].Twc_.block(0, 0, 3, 3);
            Eigen::Vector3d twc = v_frames[nc].Twc_.block(0, 3, 3, 1);
            Eigen::Matrix3d Rcw = Rwc.transpose();
            Eigen::Vector3d tcw = -Rcw * twc;

            Eigen::Vector3d Pc = Rcw * Pw + tcw;
        
            double noise_u = is_add_noise ? d(gen) : 0.0;
            double noise_v = is_add_noise ? d(gen) : 0.0;

            Eigen::Vector3d pi = K * Pc;
            double u = pi[0]/pi[2] + noise_u;
            double v = pi[1]/pi[2] + noise_v;

            Eigen::Vector2d ob(u,v);

            v_frames[nc].per_feature_.insert(make_pair(m, ob));
        }    
    }
}

int main(int argc, char** argv)
{
    std::vector<sFrame> v_sframes;
    std::vector<Eigen::Vector3d> gt_mappoints;
    std::vector<Sophus::SE3d> gt_cam_pose;
    createCameraPose(v_sframes, gt_cam_pose, num_cameras);
    createMaPoints(v_sframes, gt_mappoints);

    cv::Mat cv_K = (cv::Mat_<double>(3,3) << fx, 0,  cx, 0, fy, cy,0,  0,  1);
    Eigen::Matrix3d eigen_K;
    cv::cv2eigen(cv_K, eigen_K);
   
    detectFeatures(eigen_K, gt_mappoints, v_sframes, true);

    SFM::Ptr pSfm(new SFM(fx,fy,cx,cy));
    pSfm->Initializer(v_sframes);

    //recover scale
    for(size_t i = 0; i < v_sframes.size(); i++)
    {
        double gt_scale = gt_cam_pose[i].translation().norm();
        v_sframes[i].Twc_.block(0,3,3,1) *= gt_scale;
    }

    std::vector<Eigen::Vector3d> estimate_mappoints;
    estimate_mappoints = pSfm->Triangulate(v_sframes);

    //noise test
    // std::vector<Eigen::Vector3d> estimate_mappoints = gt_mappoints;
    // double mpt_noise = 0.05;
    // double cam_noise = 0.1;
    // double obs_noise = 1;

    // AddNoise(estimate_mappoints, v_sframes, mpt_noise, cam_noise, obs_noise);

    //【4】bundle adjustment
    BundleAdjustment ba;
    ba.SetConvergenceCondition(20, 1e-5,1e-10);

    //add camera
    for(size_t ci = 0; ci < v_sframes.size(); ci++)
    {
        Eigen::Matrix3d Rwc = v_sframes[ci].Twc_.block(0, 0, 3, 3);
        Eigen::Vector3d twc = v_sframes[ci].Twc_.block(0, 3, 3, 1);
        Eigen::Matrix3d Rcw = Rwc.transpose();

        Eigen::Vector3d tcw = -Rcw * twc;

        Sophus::SE3d se3_rt(Rcw, tcw);

        Camera* cam = new Camera(se3_rt, ci, ci==0);
        ba.AddCameraPose(cam);
    }



    //add mappoint
    for(size_t pi = 0; pi < estimate_mappoints.size(); pi++)
    {
        Eigen::Vector3d ptw = estimate_mappoints[pi];
        MapPoint* mpt = new MapPoint(ptw, pi);

        ba.AddMappoints(mpt);
    }

    //add costfunction
    for(size_t np = 0; np < estimate_mappoints.size(); np++)
    {
        MapPoint* mpt = ba.GetMapPoint(np);
        for(size_t nc = 0; nc < v_sframes.size(); nc++)
        {
            Camera* cam = ba.GetCamera(nc);
            CostFunction* cost_func = new CostFunction(cam, mpt, fx, fy, cx, cy, v_sframes[nc].per_feature_.find(np)->second);
            ba.AddCostFunctions(cost_func);
        }   
    }

     double sum_mappoint_errors1 = 0.0;
    for(size_t i = 0; i < gt_mappoints.size(); i++)
    {
        MapPoint* mpt = ba.GetMapPoint(i);
        Eigen::Vector3d opt_mpt = mpt->GetPosition();
        Eigen::Vector3d origin_mpt = gt_mappoints[i];
        sum_mappoint_errors1 += (opt_mpt - origin_mpt).norm();

    }
    std::cout << "Mean point error: " << sum_mappoint_errors1 / (double)(num_map_points) << std::endl;

    double sum_rot_error = 0.0;
	double sum_trans_error = 0.0;
	for(size_t i = 0; i < v_sframes.size(); i ++)
	{
		Camera* cam = ba.GetCamera(i);
		const Sophus::SE3d& opt_pose_cw = cam->GetPose();   //优化后的位姿态
        Eigen::Matrix3d Rwc = opt_pose_cw.rotationMatrix().transpose();
        Eigen::Vector3d twc = -Rwc * opt_pose_cw.translation();
        Sophus::SE3d opt_pose(Rwc, twc); 
		const Sophus::SE3d& org_pose = gt_cam_pose[i];
		Sophus::SE3d pose_err = opt_pose * org_pose.inverse();
		sum_rot_error += pose_err.so3().log().norm();
		sum_trans_error += pose_err.translation().norm();
	}
	std::cout << "Mean rot error: " << sum_rot_error / (double)(v_sframes.size())
	<< "\tMean trans error: " <<  sum_trans_error / (double)(v_sframes.size()) << std::endl;
    
    //optimized
    ba.Optimized();

    // compute pose rmse
    double rmse_R_error = 0.0;
    double rmse_t_error = 0.0;
    for(size_t i = 0; i < v_sframes.size(); i++)
    {
        Camera* cam = ba.GetCamera(i);
		const Sophus::SE3d& opt_pose_cw = cam->GetPose();
        Eigen::Matrix3d Rwc = opt_pose_cw.rotationMatrix().transpose();
        Eigen::Vector3d twc = -Rwc * opt_pose_cw.translation();
        Sophus::SE3d opt_pose(Rwc, twc); 
		const Sophus::SE3d& org_pose = gt_cam_pose[i];
		Sophus::SE3d pose_err = opt_pose * org_pose.inverse();
		rmse_R_error += pose_err.so3().log().norm();
		rmse_t_error += pose_err.translation().norm();
    }
    
    std::cout << "\nRMSE of rotation: " <<  rmse_R_error/(double)num_cameras << std::endl;
    std::cout << "RMSE of translation: " <<  rmse_t_error / (double)num_cameras << std::endl;

    double sum_mappoint_errors = 0.0;
    for(size_t i = 0; i < gt_mappoints.size(); i++)
    {
        MapPoint* mpt = ba.GetMapPoint(i);
        Eigen::Vector3d opt_mpt = mpt->GetPosition();
        Eigen::Vector3d origin_mpt = gt_mappoints[i];
        sum_mappoint_errors += (opt_mpt - origin_mpt).norm();      
    }

    std::cout << "Mean point error: " << sum_mappoint_errors / (double)(num_map_points) << std::endl;
 
    return 0;
}

void AddNoise(std::vector<Eigen::Vector3d>& mptpoints, std::vector<sFrame>& t_sframes, double mpt_noise, double cam_noise, double ob_noise)
{

     cv::RNG rng ( cv::getTickCount() );
    // add noise for mappoint
    for(size_t i = 0; i < mptpoints.size(); i++)
    {
        double noise_x = rng.gaussian ( mpt_noise );
        double noise_y = rng.gaussian ( mpt_noise );
        double noise_z = rng.gaussian ( mpt_noise );

        mptpoints[i] += Eigen::Vector3d(noise_x, noise_y, noise_z);
    }

    // add noise for camera
    std::normal_distribution<> noise_cam{0.0, cam_noise};
    Eigen::Matrix3d Rx, Ry, Rz;
	Eigen::Matrix3d R; 
	Eigen::Vector3d t;
    Eigen::Matrix4d T;
    for(size_t i = 0; i < t_sframes.size();i++)
    {
        // skip the first camera.
		if(i == 0)
			continue;
		
		double tz = rng.gaussian ( cam_noise );
		double ty = rng.gaussian ( cam_noise );
		double tx = rng.gaussian ( cam_noise );
		
		Rz << cos ( tz ), -sin ( tz ), 0.0,
		sin ( tz ), cos ( tz ), 0.0,
		0.0, 0.0, 1.0;
		Ry << cos ( ty ), 0.0, sin ( ty ),
		0.0, 1.0, 0.0,
		-sin ( ty ), 0.0, cos ( ty );
		Rx << 1.0, 0.0, 0.0,
		0.0, cos ( tx ), -sin ( tx ),
		0.0, sin ( tx ), cos ( tx );
		R = Rz * Ry * Rx;
		
		// translation.
		double x = rng.gaussian ( cam_noise );
		double y = rng.gaussian ( cam_noise );
		double z = rng.gaussian ( cam_noise );
		t << x, y, z;
		
		// SE3
		Sophus::SE3d cam_noise ( R, t );

        Eigen::Matrix3d Rwc = t_sframes[i].Twc_.block(0,0,3,3);
        Eigen::Vector3d twc = t_sframes[i].Twc_.block(0,3,3,1);
        Sophus::SE3d cam_temp(Rwc, twc);
        Sophus::SE3d newcamera = cam_temp * cam_noise;
		t_sframes[i].Twc_.block(0,0,3,3) = newcamera.rotationMatrix();
        t_sframes[i].Twc_.block(0,3,3,1) = newcamera.translation();

    }

    // add noise for observation
    std::normal_distribution<> noise_obs{0.0, ob_noise};
    for(size_t i = 0; i < mptpoints.size(); i++)
    {
        for(size_t j = 0; j < t_sframes.size(); j++)
        {
            double noise_u = rng.gaussian ( ob_noise );
            double noise_v = rng.gaussian ( ob_noise );

            t_sframes[j].per_feature_.find(i)->second += Eigen::Vector2d(noise_u, noise_v);
        }
    }
}



