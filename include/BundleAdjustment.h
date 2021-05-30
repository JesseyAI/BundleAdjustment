#ifndef BUNDLEADJUSTMENT_H_
#define BUNDLEADJUSTMENT_H_

#include "Camera.h"
#include <sophus/se3.hpp>
#include "MapPoint.h"
#include "CostFunction.h"
#include "TicToc.h"
#include <memory>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <iomanip>


class BundleAdjustment
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    BundleAdjustment();
    ~BundleAdjustment();

    bool Optimized();
    
    
    void SetConvergenceCondition(int max_iters, double min_delta, double min_error);  
    void AddCameraPose(Camera* cam);
    void AddMappoints(MapPoint* mappoint);
    void AddCostFunctions(CostFunction* costfunction);
    MapPoint* GetMapPoint(int mpt_id);
    Camera* GetCamera(int cam_id);

private:
    
    std::map<int, Camera*> mCameras_;
    std::map<int, MapPoint*> mMapPoints_;
    std::set<CostFunction*> mCostFunctions_;

    /* Convergence condition */
	int max_iters_;
	double min_delta_;
	double min_error_;
	double current_error_;
    double current_lamda_;
    double sum_error2_;
	double last_sum_error2_;

    double ni_;            //控制 lambda 缩放大小

    /* */

    Eigen::MatrixXd Hessian_;     // Hessian Matrix
    Eigen::VectorXd b_;     // b_ = J^T r
    Eigen::MatrixXd I_;
    Eigen::VectorXd delta_X_;	  // Delta_X_

    ulong n_cam_state_;       //number of cameras in the state vector
    ulong n_mpt_state_;       //number of mappoints in the state vector
    ulong n_state_size_;   

    void SetOrdering();    // set the dimensional of matrix 
    void MakeHessian();    // H = J^T * J
    void SolveNormalEquation();
    void ComputeLambdaInitLM();
    bool IsGoodStepInLM();
    void UpdateSates();
    void RollbackStates();
    Eigen::VectorXd PCGSolver(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, int maxIter);
	
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatXX;

};

#endif // !BUNDLEADJUSTMENT_H_
