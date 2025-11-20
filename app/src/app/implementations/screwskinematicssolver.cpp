#include "app/implementations/screwskinematicssolver.h"

#include <utility/math.h>

using namespace AIS4104;

ScrewsKinematicsSolver::ScrewsKinematicsSolver(Eigen::Matrix4d m, std::vector<Eigen::VectorXd> screws, Simulation::JointLimits limits)
    : ScrewsKinematicsSolver(std::move(m), std::move(screws), 4.e-3, 4.e-3, std::move(limits))
{
}

ScrewsKinematicsSolver::ScrewsKinematicsSolver(Eigen::Matrix4d m, std::vector<Eigen::VectorXd> space_screws, double v_e, double w_e, Simulation::JointLimits limits)
    : KinematicsSolver(std::move(limits))
    , m_ve(v_e)
    , m_we(w_e)
    , m_m(std::move(m))
    , m_screws(std::move(space_screws))
{
}

void ScrewsKinematicsSolver::set_epsilons(double v_e, double w_e)
{
    m_ve = v_e;
    m_we = w_e;
}

uint32_t ScrewsKinematicsSolver::joint_count() const
{
    return m_screws.size();
}

//TASK: FIN Implement fk_solve using screws.
Eigen::Matrix4d ScrewsKinematicsSolver::fk_solve(const Eigen::VectorXd &joint_positions)
{
    auto [M, S] = space_chain();

    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();

    for (int i = 0; i < S.size(); ++i)
    {
        const Eigen::VectorXd screw = S[i];
        double theta = joint_positions[i];

        Eigen::Matrix4d T_i = utility::matrix_exponential(screw, theta);

        T = T * T_i;
    }

    T = T * M;

    return T;
    //return Eigen::Matrix4d::Identity();
}

Eigen::VectorXd ScrewsKinematicsSolver::ik_solve(const Eigen::Matrix4d &t_sd, const Eigen::VectorXd &j0)
{
    return ik_solve(t_sd, j0, [&](const std::vector<Eigen::VectorXd> &) { return 0u; });
}

//TASK: FIN Implement ik_solve using screws.
Eigen::VectorXd ScrewsKinematicsSolver::ik_solve(const Eigen::Matrix4d &t_sd, const Eigen::VectorXd &j0, const std::function<uint32_t(const std::vector<Eigen::VectorXd> &)> &solution_selector)
{
    return j0;
}

std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> ScrewsKinematicsSolver::space_chain()
{
    return {m_m, m_screws};
}

//TASK: Implement body_chain(). You can obtain the variables to transform to body frame from space_chain().
std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> ScrewsKinematicsSolver::body_chain()
{

    // auto [m, screws] = space_chain();
    return space_chain();
}

//TASK: Implement space_jacobian() using space_chain()
Eigen::MatrixXd ScrewsKinematicsSolver::space_jacobian(const Eigen::VectorXd &current_joint_positions)
{
    auto [M, S] = space_chain();

    uint32_t n = S.size();
    Eigen::MatrixXd Js(6, n);

    Eigen::Matrix4d T_i_minus_1 = Eigen::Matrix4d::Identity();

    for (uint32_t i = 0; i < n; ++i)
    {
        const Eigen::VectorXd S_i = S[i];
        double theta_i = current_joint_positions[i];

        Eigen::VectorXd Js_i = utility::adjoint_map(S_i, T_i_minus_1);

        Js.col(i) = Js_i;

        Eigen::Matrix4d T_i = utility::matrix_exponential(S_i, theta_i);
        T_i_minus_1 = T_i_minus_1 * T_i;
    }

    return Js;
    //return Eigen::MatrixXd::Identity(current_joint_positions.size(), current_joint_positions.size());
}

//TASK: Implement body_jacobian() using body_chain()
Eigen::MatrixXd ScrewsKinematicsSolver::body_jacobian(const Eigen::VectorXd &current_joint_positions)
{
    return Eigen::MatrixXd::Identity(current_joint_positions.size(), current_joint_positions.size());
}
