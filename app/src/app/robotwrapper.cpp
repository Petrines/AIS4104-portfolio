#include "app/robotwrapper.h"

#include <utility/math.h>
#include <utility/vectors.h>

using namespace AIS4104;

RobotWrapper::RobotWrapper(std::shared_ptr<threepp::Robot> robot, std::shared_ptr<Simulation::KinematicsSolver> solver)
    : m_tool_transform(Eigen::Matrix4d::Identity())
    , m_robot(std::move(robot))
    , m_solver(std::move(solver))
{
}

threepp::Robot& RobotWrapper::threepp_robot()
{
    return *m_robot;
}

std::shared_ptr<threepp::Robot> RobotWrapper::threepp_robot_ptr()
{
    return m_robot;
}

const threepp::Robot& RobotWrapper::threepp_robot() const
{
    return *m_robot;
}

uint8_t RobotWrapper::joint_count() const
{
    return m_solver->joint_count();
}

//TASK FIN: Implement the function to calculate the joint positions for the desired tool pose
// a) Use m_tool_transform to calculate the flange pose required by m_solver.ik_solve()
// b) Use the m_solver.ik_solve() overload with the solution selector lambda to choose the most desirable IK solution.
Eigen::VectorXd RobotWrapper::ik_solve_pose(const Eigen::Matrix4d &eef_pose, const Eigen::VectorXd &j0) const
{
    Eigen::Matrix4d flange_pose = eef_pose * m_tool_transform.inverse();
    return ik_solve_flange_pose(flange_pose,j0);
//add det ej fjerna?
}

//FINISHED: Implement the function to calculate the joint positions for the desired flange pose
// a) Use m_tool_transform to calculate the flange pose required by m_solver.ik_solve()
// b) Use the m_solver.ik_solve() overload with the solution selector lambda to choose the most desirable IK solution.
Eigen::VectorXd RobotWrapper::ik_solve_flange_pose(const Eigen::Matrix4d &flange_pose, const Eigen::VectorXd &j0) const
{
//Bør kanskljke adde det ej fjerna??
    return m_solver->ik_solve(flange_pose, j0);
    //return joint_positions();
}

Eigen::Matrix4d RobotWrapper::tool_transform() const
{
    return m_tool_transform;
}

void RobotWrapper::set_tool_transform(Eigen::Matrix4d transform)
{
    m_tool_transform = std::move(transform);
}

//TASK FIN: Calculate the pose of the end effector using forward kinematics;
// Relevant variables are m_solver and m_tool_transform.
Eigen::Matrix4d RobotWrapper::current_pose() const
{
    Eigen::Matrix4d T_flange = current_flange_pose();
    return T_flange * m_tool_transform;
    //return Eigen::Matrix4d::Identity();
}

//TASK FIN: Calculate the position of the end effector using forward kinematics.
// Relevant variables are m_solver and m_tool_transform (or possibly another function of RobotWrapper?).
Eigen::Vector3d RobotWrapper::current_position() const
{
    Eigen::Matrix4d T = current_pose();
    return T.block<3, 1>(0, 3);
    //return Eigen::Vector3d::Zero();
}

//TASK FIN: Calculate the orientation of the end effector using forward kinematics and m_solver (or rely on another function of RobotWrapper?).
Eigen::Vector3d RobotWrapper::current_orientation_zyx() const
{
    Eigen::Matrix4d T = current_pose();
    Eigen::Matrix3d R = T.block<3, 3>(0, 0);
    return utility::euler_zyx_from_rotation_matrix(R);
    //return Eigen::Vector3d::Zero();
}

//FINISHED: Calculate the pose of the end effector using forward kinematics and m_solver.
Eigen::Matrix4d RobotWrapper::current_flange_pose() const
{
    Eigen::VectorXd current_joints = joint_positions(); //litt unødvendig steg jaja
    return m_solver->fk_solve(current_joints);
    //return Eigen::Matrix4d::Identity();
}

//FINISHED: Based on the flange pose, return its linear position.
Eigen::Vector3d RobotWrapper::current_flange_position() const
{
    //Eigen::Matrix4d T = current_flange_pose();
    //return T.block<3, 1>(0, 3);
    //return Eigen::Vector3d::Zero();

    return current_flange_pose().block<3, 1>(0, 3); // kan vere ej endra
}

//FINISHED: Based on the flange pose, return its orientation in the Euler ZYX representation.
Eigen::Vector3d RobotWrapper::current_flange_orientation_zyx() const
{
    Eigen::Matrix4d T = current_flange_pose();
    Eigen::Matrix3d R = T.block<3, 3>(0, 0);
    return utility::euler_zyx_from_rotation_matrix(R);
    //return Eigen::Vector3d::Zero();
}

const Simulation::JointLimits& RobotWrapper::joint_limits() const
{
    return m_solver->joint_limits();
}

Eigen::VectorXd RobotWrapper::joint_positions() const
{
    return utility::to_eigen_vectord(m_robot->jointValues());
}

void RobotWrapper::set_joint_positions(const Eigen::VectorXd &joint_positions)
{
    m_robot->setJointValues(utility::to_std_vectorf(joint_positions));
}
