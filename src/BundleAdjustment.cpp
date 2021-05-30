#include "../include/BundleAdjustment.h"

BundleAdjustment::BundleAdjustment()
     :max_iters_ ( 20 ), min_delta_ ( 1e-10 ), min_error_ ( 1e-10 )
{

}

BundleAdjustment::~BundleAdjustment()
{
    //delete pointer
    std::map<int, Camera*>::iterator it_cam;
    for ( it_cam = mCameras_.begin(); it_cam != mCameras_.end(); it_cam ++ ) {
        Camera* cam = it_cam->second;
        delete cam;
    } // for all cameras

    std::map<int, MapPoint*>::iterator it_mpt;
    for ( it_mpt = mMapPoints_.begin(); it_mpt != mMapPoints_.end(); it_mpt++ ) {
        MapPoint* mpt = it_mpt->second;
        delete mpt;
    }
}

void BundleAdjustment::AddCameraPose(Camera* cam)
{
    mCameras_[cam->GetCamId()] = cam;
}

void BundleAdjustment::AddMappoints(MapPoint* mappoint)
{
    mMapPoints_[mappoint->GetMapPointId()] = mappoint;
}

void BundleAdjustment::AddCostFunctions(CostFunction* costfunction)
{
    mCostFunctions_.insert( costfunction);
}

MapPoint* BundleAdjustment::GetMapPoint(int mpt_id)
{
    std::map<int, MapPoint*>::iterator iter = mMapPoints_.find ( mpt_id );
    if ( iter != mMapPoints_.end() ) {
        return iter->second;
    }
    return nullptr;
}

Camera* BundleAdjustment::GetCamera(int cam_id)
{
    std::map<int, Camera*>::iterator iter = mCameras_.find ( cam_id );
    if ( iter != mCameras_.end() ) {
        Camera* cam = iter->second;
        // cam->GetPose().rotationMatrix() = cam->GetPose().rotationMatrix().transpose();
        // cam->GetPose().translation() = cam->GetPose().rotationMatrix() * cam->GetPose().translation();
        return cam;  //return Rwc and twc
    }
    return nullptr;
}


//set convergence condition
void BundleAdjustment::SetConvergenceCondition(int max_iters, double min_delta, double min_error)
{
    max_iters_ = max_iters;
    min_delta_ = min_delta;
    min_error_ = min_error;
}

//LM
bool BundleAdjustment::Optimized()
{
    TicToc t_solve;

    SetOrdering();

    MakeHessian();

    //Initilization for LM
    ComputeLambdaInitLM();

     bool stop = false;
    TicToc perIteration;
    // Optimized by LM
    int nIter = 0;
    while(!stop && (nIter < max_iters_))
    {
        perIteration.Tic();

        std::cout << "iter:" << nIter << ", chi = " << current_error_ << ", Lambda = " << current_lamda_;

        bool isOneStepSucess = false;
        int false_cnt = 0;

        while(!isOneStepSucess)
        {
            SolveNormalEquation();

            //condition 1
            //std::cout <<"\n 更新量： \n" << delta_X_ << std::endl;
            if (delta_X_.squaredNorm() <= 1e-6 || false_cnt > 10)
            {
                stop = true;
                break;
            }

            UpdateSates();

            isOneStepSucess = IsGoodStepInLM();

            if(isOneStepSucess)
            {
                MakeHessian();
                false_cnt = 0;
            }
            else
            {
                false_cnt ++;
                RollbackStates();
            }

        }
        nIter++;
        std::cout <<", perIteration = "<< perIteration.Toc() <<" ms" << std::endl;

        //优化退出条件
        if(sqrt(current_error_) <= min_error_)
            stop = true;

    }
    std::cout << "problem solve cost: " << t_solve.Toc() << " ms" << std::endl;
    return true;
}

/*
 * 为了构造hessian的方便，这里可以考虑将cam和mpt的state_index
 */
void BundleAdjustment::SetOrdering()
{
    int INDEX = 0;
    n_cam_state_ = 0;
    std::map<int, Camera*>::iterator it_cam;
    for(it_cam = mCameras_.begin(); it_cam != mCameras_.end(); it_cam ++)
    {
        Camera* cam = it_cam->second;
        if( cam->IsFixed() == false)
        {
            cam->state_index_ = INDEX;
            INDEX ++;
            n_cam_state_ ++;
        }
        else
        {
            cam->state_index_ = 0;
        }
    }

    n_mpt_state_ = 0;
    INDEX = 0;         // 这里为了合并状态量，应该不用INDEX=0
    std::map<int, MapPoint*>::iterator it_mpt;
    for(it_mpt = mMapPoints_.begin(); it_mpt != mMapPoints_.end(); it_mpt ++)
    {
        MapPoint* mpt = it_mpt->second;
        mpt->state_index_ = INDEX;
            INDEX++;
        n_mpt_state_ ++;
    }

    n_state_size_ = 6 * n_cam_state_ + 3 * n_mpt_state_;

    // set ordering
    Hessian_.resize(n_state_size_, n_state_size_);
    b_.resize(n_state_size_);
    delta_X_.resize(n_state_size_);

    I_.resize ( n_state_size_, n_state_size_ );
    for ( size_t i = 0; i < n_state_size_; i ++ ) {
        I_ ( i,i ) = 1.0;
    }
}

void BundleAdjustment::ComputeLambdaInitLM()
{
    ni_ = 2;
    current_lamda_ = -1;
    current_error_ = 0.0;

    std::set<CostFunction*>::iterator iter;
    for ( iter = mCostFunctions_.begin(); iter != mCostFunctions_.end(); iter ++ )
    {
        CostFunction* cost_func = *iter;

        Eigen::Vector2d e;
        Eigen::Matrix2d weighted_info;
        double weighted_e2;
        cost_func->ComputeResiduals(e, weighted_info, weighted_e2 );
        current_error_ += weighted_e2;
    }

    double max_diagonal = 0;
    long size = Hessian_.cols();
    assert(Hessian_.rows() == Hessian_.cols() && "Hessian is not square");
    for(long i = 0; i < size; i++)
    {
        max_diagonal = std::max(std::fabs(Hessian_(i,i)), max_diagonal);
    }

    double tau = 1e-8;
    current_lamda_ = tau * max_diagonal;
}

void BundleAdjustment::MakeHessian()
{
    //直接构造大的hessian矩阵
    sum_error2_ = 0.0;
    Eigen::MatrixXd H(Eigen::MatrixXd::Zero(n_state_size_, n_state_size_));
    //MatXX H1(MatXX::Zero(n_state_size_, n_state_size_));
    Eigen::VectorXd b(Eigen::VectorXd::Zero(n_state_size_));

    std::set<CostFunction*>::iterator iter;

    for ( iter = mCostFunctions_.begin(); iter != mCostFunctions_.end(); iter ++ )
    {
        CostFunction* cost_func = *iter;

        //Compute Jacobian
        Camera* cam = cost_func->p_camera_;

        MapPoint* mpt = cost_func->p_mappoint_;

        Eigen::Matrix<double, 2, 6> Jtemp1;  //关于pose
        Eigen::Matrix<double, 2, 3> Jtemp2;  //关于三维点

        cost_func->ComputeJacobian(Jtemp1, Jtemp2);

        Eigen::Vector2d e;
        Eigen::Matrix2d weighted_info;
        double weighted_e2;
        cost_func->ComputeResiduals(e, weighted_info, weighted_e2 );

        if(cam->IsFixed() == true)
            Jtemp1.setZero();

        int dim_cam = Jtemp1.cols();
        int dim_mpt = Jtemp2.cols();

        Eigen::MatrixXd hessian11 = Jtemp1.transpose() * weighted_info * Jtemp1;
        Eigen::MatrixXd hessian12 = Jtemp1.transpose() * weighted_info * Jtemp2;
        Eigen::MatrixXd hessian21 = Jtemp2.transpose() * weighted_info * Jtemp1;
        Eigen::MatrixXd hessian22 = Jtemp2.transpose() * weighted_info * Jtemp2;

        H.block(cam->state_index_ * 6, cam->state_index_ * 6, dim_cam, dim_cam).noalias() += hessian11;
        H.block(n_cam_state_ * 6 + mpt->state_index_ * 3, n_cam_state_ * 6 + mpt->state_index_ * 3, dim_mpt,dim_mpt).noalias() += hessian22;
        H.block(cam->state_index_ * 6, n_cam_state_ * 6 + mpt->state_index_ * 3, dim_cam, dim_mpt).noalias() += hessian12;
        H.block(n_cam_state_ * 6 + mpt->state_index_ * 3, cam->state_index_ * 6, dim_mpt, dim_cam).noalias() += hessian21;

        b.segment(cam->state_index_ * 6 , dim_cam).noalias() -= Jtemp1.transpose() * weighted_info *e;
        b.segment(n_cam_state_ * 6 + mpt->state_index_ * 3, dim_mpt).noalias() -= Jtemp2.transpose() * weighted_info * e;

        sum_error2_ += weighted_e2; // add all error2.
    }

    Hessian_ = H;
    b_ = b;

     //std::cout << "Hessian_\n" << Hessian_ << std::endl;
    //std::cout << "b_\n" << b_ << std::endl;

}

/*
 * solve by schur
 */
void BundleAdjustment::SolveNormalEquation()
{
    //Hessian_ += current_lamda_ * I_;

    Eigen::MatrixXd U = Hessian_.block(0, 0, n_cam_state_ * 6, n_cam_state_ * 6);
    Eigen::MatrixXd W = Hessian_.block(0, n_cam_state_ * 6, n_cam_state_ * 6, n_mpt_state_ * 3);
    Eigen::MatrixXd Wt = Hessian_.block(n_cam_state_ * 6, 0, n_mpt_state_ * 3, n_cam_state_ * 6);
    Eigen::MatrixXd V = Hessian_.block(n_cam_state_ * 6, n_cam_state_ * 6, n_mpt_state_ * 3, n_mpt_state_ * 3);

    Eigen::VectorXd bc = b_.segment(0, n_cam_state_ * 6);
    Eigen::VectorXd bm = b_.segment(n_cam_state_ * 6,  n_mpt_state_ * 3);

    Eigen::MatrixXd V_inv = V;

    for(ulong i = 0; i < n_mpt_state_; i++)
    {
        V_inv.block(3 * i, 3 * i, 3, 3) = V_inv.block(3*i, 3*i, 3, 3).inverse();
    }


    Eigen::MatrixXd temp = W*V_inv;
    Eigen::MatrixXd V_schur = U - temp*Wt;
    Eigen::VectorXd bc_schur = bc - temp * bm;

    for(ulong i = 0; i < n_cam_state_ * 6; i++)
    {
        V_schur(i, i) += current_lamda_;
    }

    Eigen::VectorXd delta_X_cam(Eigen::VectorXd::Zero(n_cam_state_ * 6));
    delta_X_cam = V_schur.ldlt().solve(bc_schur);
    delta_X_.head(n_cam_state_ * 6) = delta_X_cam;

    //solve delta_b
    Eigen::VectorXd delta_X_mpt(Eigen::VectorXd::Zero(n_mpt_state_ * 3));
    delta_X_mpt = V_inv * (bm - Wt*delta_X_cam);
    delta_X_.tail(n_mpt_state_ * 3) = delta_X_mpt;

}

/*
 * Nielsen's Strategy(g2o)
 */
bool BundleAdjustment::IsGoodStepInLM()
{
    double scale = 0.0;
    scale = delta_X_.transpose() *(current_lamda_ * delta_X_ + b_);   //这里是否有问题？？？？
    scale += 1e-3;     //make sure it is non-zero

    double temp_error = 0.0;
    std::set<CostFunction*>::iterator iter;
    for ( iter = mCostFunctions_.begin(); iter != mCostFunctions_.end(); iter ++ )
    {
        CostFunction* cost_func = *iter;
        Eigen::Vector2d e;
        Eigen::Matrix2d weighted_info;
        double weighted_e2;
        cost_func->ComputeResiduals(e, weighted_info, weighted_e2 );
        temp_error += weighted_e2;
    }

    double rho = (current_error_ - temp_error) / scale;

    if(rho > 0)
    {
        double alpha = 1 - pow((2 * rho - 1),3);
        alpha = std::min(alpha, 2. / 3.);
        double scaleFactor = std::max(1. / 3., alpha);
        current_lamda_ *= scaleFactor;
        ni_ = 2;
        current_error_ = temp_error;
        return true;
    }
    else
    {
        current_lamda_ *= ni_;
        ni_ *= 2;
        return false;
    }
}

void BundleAdjustment::UpdateSates()
{
    std::map<int, Camera*>::iterator it_cam;
    for(it_cam = mCameras_.begin(); it_cam != mCameras_.end(); it_cam++)
    {
        Camera* cam = it_cam->second;
        if(cam->IsFixed() == false)
        {
            int index = cam->state_index_;
            cam->UpdateCamPose(delta_X_.block(index * 6, 0, 6, 1));
        }
    }

    std::map<int, MapPoint*>::iterator it_mpt;
    for(it_mpt = mMapPoints_.begin(); it_mpt != mMapPoints_.end(); it_mpt++)
    {
        MapPoint* mpt = it_mpt->second;

        int index = mpt->state_index_;

        mpt->UpdatePosition(delta_X_.block(n_cam_state_ * 6 + index * 3, 0, 3, 1));
    }
}

void BundleAdjustment::RollbackStates()
{
    std::map<int, Camera*>::iterator it_cam;
    for(it_cam = mCameras_.begin(); it_cam != mCameras_.end(); it_cam++)
    {
        Camera* cam = it_cam->second;
        if(cam->IsFixed() == false)
        {
            int index = cam->state_index_;
            cam->UpdateCamPose(- delta_X_.block(index * 6, 0, 6, 1));
        }     
    }

    std::map<int, MapPoint*>::iterator it_mpt;
    for(it_mpt = mMapPoints_.begin(); it_mpt != mMapPoints_.end(); it_mpt++)
    {
        MapPoint* mpt = it_mpt->second;

        int index = mpt->state_index_;

        mpt->UpdatePosition(-delta_X_.block(n_cam_state_ * 6 + index * 3, 0, 3, 1));
    }
}

Eigen::VectorXd BundleAdjustment::PCGSolver(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, int maxIter = -1)
{
     assert(A.rows() == A.cols() && "PCG solver ERROR: A is not a square matrix");
    int rows = b.rows();
    int n = maxIter < 0 ? rows : maxIter;
    Eigen::VectorXd x(Eigen::VectorXd::Zero(rows));
    Eigen::MatrixXd M_inv = A.diagonal().asDiagonal().inverse();
    Eigen::VectorXd r0(b);  // initial r = b - A*0 = b
    Eigen::VectorXd z0 = M_inv * r0;
    Eigen::VectorXd p(z0);
    Eigen::VectorXd w = A * p;
    double r0z0 = r0.dot(z0);
    double alpha = r0z0 / p.dot(w);
    Eigen::VectorXd r1 = r0 - alpha * w;
    int i = 0;
    double threshold = 1e-6 * r0.norm();
    while (r1.norm() > threshold && i < n) {
        i++;
        Eigen::VectorXd z1 = M_inv * r1;
        double r1z1 = r1.dot(z1);
        double belta = r1z1 / r0z0;
        z0 = z1;
        r0z0 = r1z1;
        r0 = r1;
        p = belta * p + z1;
        w = A * p;
        alpha = r1z1 / p.dot(w);
        x += alpha * p;
        r1 -= alpha * w;
    }
    return x;
}
