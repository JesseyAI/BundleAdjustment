#include "../include/MapPoint.h"

MapPoint::MapPoint(const Eigen::Vector3d& position, int id):
  position_(position), id_(id), state_index_(-1)
{

}

MapPoint::~MapPoint()
{
}

int MapPoint::GetMapPointId()
{
    return id_;
}

Eigen::Vector3d MapPoint::GetPosition()
{
  return position_;
}

void MapPoint::UpdatePosition(const Eigen::Vector3d& delta_position)
{
    position_ += delta_position;
} // addDeltaPosition