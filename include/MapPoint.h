#ifndef MAPPOINT_H_
#define MAPPOINT_H_

#include <Eigen/Core>


class MapPoint
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    MapPoint(const Eigen::Vector3d& position, int id);
    ~MapPoint();

    int GetMapPointId();
    Eigen::Vector3d GetPosition();
    void UpdatePosition(const Eigen::Vector3d& delta_position);

    int state_index_;
private:
    Eigen::Vector3d position_;
    int id_;
};




#endif // !MAPPOINT_H_
