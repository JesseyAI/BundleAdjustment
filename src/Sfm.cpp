#include "../include/Sfm.h"

SFM::SFM()
{

}

SFM::SFM(double fx, double fy, double cx, double cy)
        : fx_(fx), fy_(fy), cx_(cx), cy_(cy)
{
    K_ = cv::Mat(3,3,CV_64F);
    K_ << fx, 0,  cx,
	       0, fy, cy,
		   0, 0,  1;
}

SFM::~SFM()
{

}

bool SFM::Initializer(std::vector<sFrame>& frames)
{
    const int N = frames.size();

	if(N < 2)
	{
		std::cerr << "there is enough cameras !" << std::endl;
		return false;
	}

    for(int i = 0; i < N; i++)
    {
		//fixed the first cameras
        if(i == 0)
		{
			continue;
		}

		Eigen::Matrix3d F;

		Eigen::Matrix3d Rwc;
		Eigen::Vector3d twc;

        RecoverPose(frames[i].per_feature_, frames[0].per_feature_, Rwc ,twc );
        
        //std::cout << "\n 真实的位: \n" << frames[i].Twc_ << std::endl;

        frames[i].Twc_.block(0, 0, 3, 3) = Rwc;
        frames[i].Twc_.block(0, 3, 3, 1) = twc;

        //std::cout << "\n 恢复的姿态: \n" << frames[i].Twc_ << std::endl;

        
    }

    return true;
}

void SFM::RecoverPose(const std::unordered_map<int,Eigen::Vector2d>& umpts1, const std::unordered_map<int, Eigen::Vector2d>& umpts2, Eigen::Matrix3d& R, Eigen::Vector3d& t)
{
	const int N = umpts1.size();

	std::vector<cv::Point2d> vpts1,vpts2;
    for(int i = 0; i < N; i++)
    {
        Eigen::Vector2d temp_pt1, temp_pt2;
        temp_pt1 = umpts1.find(i)->second;
        temp_pt2 = umpts2.find(i)->second;

        vpts1.push_back(cv::Point2d(temp_pt1(0), temp_pt1(1)));
        vpts2.push_back(cv::Point2d(temp_pt2(0), temp_pt2(1)));
    }

    //随机选择八对
    // std::vector<cv::Point2d> vpn1i(8);
    // std::vector<cv::Point2d> vpn2i(8);

    // for(int j = 0; j < 8; j++)
    // {
    //     vpn1i[j] = vpts1[j];
    //     vpn2i[j] = vpts2[j];
    // }

    // cv::Mat F = ComputeF(vpn1i, vpn2i);

    // cv::Mat_<double> R1,R2,t1,t2;
    // DecomposeF(F, R1, R2, t1, t2);

    // if(cv::determinant(R1) + 1.0 < 1e-09)
    // {
    //     F = -F;
    //     DecomposeF(F, R1, R2, t1, t2);
    // }

   
    cv::Mat inline_mask;
    cv::Mat cv_E =  cv::findEssentialMat(vpts1, vpts2, K_, cv::RANSAC, 0.999, 1.0, inline_mask);
    //int num_inlier_E = cv::countNonZero(inline_mask);
    cv::Mat res_R_cv, res_t_cv;
    cv::recoverPose(cv_E, vpts1, vpts2, K_, res_R_cv, res_t_cv, inline_mask);

    cv::cv2eigen(res_R_cv,R);
    cv::cv2eigen(res_t_cv,t);
	
}



cv::Mat SFM::ComputeF(const std::vector<cv::Point2d>& vpts1, const std::vector<cv::Point2d>& vpts2)
{
    const int N = vpts1.size();

    //n*8的系数矩阵
    
	cv::Mat A(N,9,CV_64F);

    for(int i = 0;i < N; i++)
    {
        const double u1 = vpts1[i].x;
        const double v1 = vpts1[i].y;
        const double u2 = vpts2[i].x;
        const double v2 = vpts2[i].y;

        A.at<double>(i,0) = u2 * u1;
        A.at<double>(i,1) = u2 * v1;
        A.at<double>(i,2) = u2;
        A.at<double>(i,3) = v2 * u1;
        A.at<double>(i,4) = v2 * v1;
        A.at<double>(i,5) = v2;
        A.at<double>(i,6) = u1;
        A.at<double>(i,7) = v1;
        A.at<double>(i,8) = 1;
    }
 
	cv::Mat u,w,vt;
 
	cv::SVDecomp(A,w,u,vt,cv::SVD::MODIFY_A | cv::SVD::FULL_UV);
 
	cv::Mat Fpre = vt.row(8).reshape(0, 3);// F 基础矩阵的秩为2 需要在分解 后 取对角矩阵 秩为2 在合成F
 
	cv::SVDecomp(Fpre,w,u,vt,cv::SVD::MODIFY_A | cv::SVD::FULL_UV);
	w.at<double>(2)=0;//  基础矩阵的秩为2，重要的约束条件
 
	return  u * cv::Mat::diag(w)  * vt;// 在合成      
}


void SFM::DecomposeF(cv::Mat F, cv::Mat_<double> &R1, cv::Mat_<double> &R2, cv::Mat_<double> &t1, cv::Mat_<double> &t2)
{

    cv::Mat e = K_.t() * F * K_; 

    cv::SVD svd(e, cv::SVD::MODIFY_A);

    cv::Matx33d W(0, -1, 0,
                  1, 0, 0,
                  0, 0, 1);
    cv::Matx33d Wt(0, 1, 0,
                   -1, 0, 0,
                   0, 0, 1);
    R1 = svd.u * cv::Mat(W) * svd.vt;
    R2 = svd.u * cv::Mat(Wt) * svd.vt;
    t1 = svd.u.col(2);
    t2 = -svd.u.col(2);
}

double SFM::TestTriangulate(const std::vector<cv::Point2d>& vpts1, const std::vector<cv::Point2d>& vpts2, cv::Mat_<double> R, cv::Mat_<double> t)
{
	
	cv::Mat pointcloud;
    // Camera 1 projection matrix K[I | 0]
    cv::Mat P1(3, 4, CV_64F, cv::Scalar(0));
    K_.copyTo(P1.rowRange(0,3).colRange(0,3));

    // Camera 2 projection matrix K[Rcw | tcw]
    cv::Mat P2(3,4, CV_64F);
    cv::Mat Rcw = R.t();
    cv::Mat tcw = -Rcw * t;
    Rcw.copyTo(P2.rowRange(0,3).colRange(0,3));
    tcw.copyTo(P2.rowRange(0,3).col(3));
    P2 = K_ * P2;

    cv::triangulatePoints(P1, P2, vpts1, vpts2, pointcloud);

    int front_count = 0;
    for (int i = 0; i < pointcloud.cols; i++)
    {
        double normal_factor = pointcloud.col(i).at<double>(3);

        cv::Mat_<double> p_3d_l = P1 * (pointcloud.col(i) / normal_factor);
        
        cv::Mat_<double> p_3d_r = P2 * (pointcloud.col(i) / normal_factor);

        if (p_3d_l(2) > 0 && p_3d_r(2) > 0)
            front_count++;
    }
    //ROS_DEBUG("MotionEstimator: %f", 1.0 * front_count / pointcloud.cols);                        
    return 1.0 * front_count / pointcloud.cols;                                         
}

std::vector<Eigen::Vector3d> SFM::Triangulate(const std::vector<sFrame>& v_sframes)
{
	//观测次数
    const size_t n_cam = v_sframes.size();

    //特征点个数
    const size_t n_points = v_sframes[0].per_feature_.size();

    std::vector<Eigen::Vector3d> vpts_world;
    Eigen::Matrix3d eigen_K;
    cv::cv2eigen(K_, eigen_K);

    for(size_t np = 0; np < n_points; np++)
    {
        Eigen::Vector3d Pw;
        Pw.setZero();

        //定义一个动态矩阵，用来存储系数矩阵
        Eigen::MatrixXd D;
        D.resize(2 * n_cam, 4);    //每次观测提供两个方程，因此系数矩阵的维度为2M*4
        D.setZero();

        int p = 0; //系数矩阵D的行
        //系数矩阵填充
        for(size_t nc = 0; nc < n_cam; nc++)
        {
            Eigen::Matrix3d Rwc = v_sframes[nc].Twc_.block(0, 0, 3, 3);
            Eigen::Vector3d twc = v_sframes[nc].Twc_.block(0, 3, 3, 1);
            Eigen::Matrix3d Rcw = Rwc.transpose();
            Eigen::Vector3d tcw = -Rcw * twc;

            Eigen::Matrix<double, 3, 4> Pcw = Eigen::Matrix<double, 3, 4>::Zero();

            Pcw << Rcw, tcw;

            Pcw = eigen_K * Pcw;
            D.row(p) = v_sframes[nc].per_feature_.find(np)->second.x() * Pcw.row(2) - Pcw.row(0);
            D.row(++p) = v_sframes[nc].per_feature_.find(np)->second.y() * Pcw.row(2) - Pcw.row(1);
            p++;
        }

        //求D^t*D
        Eigen::Matrix4d DtD = D.transpose()*D;

        //对D^t*D矩阵作SVD分解
        Eigen::JacobiSVD<Eigen::Matrix4d> svd(DtD, Eigen::ComputeFullV);
        Eigen::Matrix4d V = svd.matrixV();

        //最后的三维点坐标
        Pw = Eigen::Vector3d(V(0,3), V(1,3), V(2,3)) / V(3,3);

        vpts_world.push_back(Pw);
    }

    
    return vpts_world;
}

// std::vector<Eigen::Vector3d> SFM::Triangulate(const std::vector<sFrame>& v_sframes)
// {
// 	//观测次数
//     const size_t n_cam = v_sframes.size();

//     //特征点个数
//     const size_t n_points = v_sframes[0].per_feature_.size();

//     std::vector<Eigen::Vector3d> vpts_world;
//     Eigen::Matrix3d eigen_K;
//     cv::cv2eigen(K_, eigen_K);

//     for(size_t np = 0; np < n_points; np++)
//     {
//         Eigen::Vector3d Pw;
//         Pw.setZero();

//         //定义一个动态矩阵，用来存储系数矩阵
//         Eigen::MatrixXd D;
//         D.resize(2 * n_cam, 4);    //每次观测提供两个方程，因此系数矩阵的维度为2M*4
//         D.setZero();

//         //系数矩阵填充
//         for(size_t nc = 0; nc < n_cam; nc++)
//         {
//             Eigen::Matrix3d Rwc = v_sframes[nc].Twc_.block(0, 0, 3, 3);
//             Eigen::Vector3d twc = v_sframes[nc].Twc_.block(0, 3, 3, 1);
//             Eigen::Matrix3d Rcw = Rwc.transpose();
//             Eigen::Vector3d tcw = -Rcw * twc;

//             Eigen::Matrix<double, 3, 4> Pcw = Eigen::Matrix<double, 3, 4>::Zero();

//             Pcw << Rcw, tcw;

//             Pcw = eigen_K * Pcw;

//             D.block(2 * nc, 0, 1, 4) = v_sframes[nc].per_feature_.find(np)->second.x() * Pcw.row(2) - Pcw.row(0);
//             D.block(2 * nc + 1, 0, 1, 4) = v_sframes[nc].per_feature_.find(np)->second.y() * Pcw.row(2) - Pcw.row(1);
//         }

//         //找出矩阵D中绝对值最大的元素
//         double max_element = D.cwiseAbs().maxCoeff();
//         Eigen::MatrixXd DtD((D/max_element).transpose() * (D/max_element));

//         //对DtD作奇异值分解
//         Eigen::JacobiSVD<Eigen::MatrixXd> svd(DtD, Eigen::ComputeThinU | Eigen::ComputeThinV);

//         // //判断第四个特征值是不是足够小
//         // if( std::abs(svd.singularValues()[3] / svd.singularValues()[2]) > 1e-2 )
//         // {
//         //     std::cout << " the forth sigularValues is not enough small ! " << std::endl;
//         // }

//         Eigen::Vector4d u4 = max_element * svd.matrixU().rightCols(1);
//         Pw = (u4 / u4[3]).head(3);

//         vpts_world.push_back(Pw);
//     }

    
//     return vpts_world;
// }



