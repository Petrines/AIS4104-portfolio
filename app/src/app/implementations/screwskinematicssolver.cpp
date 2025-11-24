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

//FINISHED: Implement fk_solve using screws.
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

//FINISHED: Implement ik_solve using screws.
Eigen::VectorXd ScrewsKinematicsSolver::ik_solve(const Eigen::Matrix4d &t_sd, const Eigen::VectorXd &j0, const std::function<uint32_t(const std::vector<Eigen::VectorXd> &)> &solution_selector)
{
    const int MAX_ITER = 20;

    const double WE_TOL = m_we;
    const double VE_TOL = m_ve;

    Eigen::VectorXd current_joints = j0;

    for (int iter = 0; iter < MAX_ITER; ++iter)
    {
        Eigen::Matrix4d T_sb = fk_solve(current_joints);
        Eigen::Matrix4d T_bs = T_sb.inverse();
        Eigen::Matrix4d T_bd = T_bs * t_sd;
        auto log_T_pair = utility::matrix_logarithm(T_bd);
        Eigen::VectorXd S_e = log_T_pair.first;
        double theta_e = log_T_pair.second;

        Eigen::VectorXd V_e = S_e * theta_e;
        double w_error = V_e.head<3>().norm();
        double v_error = V_e.tail<3>().norm();

        if (w_error < WE_TOL && v_error < VE_TOL)
        {
            return current_joints;
        }
        Eigen::MatrixXd Jb = body_jacobian(current_joints);
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(Jb, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::VectorXd singular_values = svd.singularValues();
        double tolerance = std::numeric_limits<double>::epsilon() * std::max(Jb.rows(), Jb.cols()) * singular_values.array().abs().maxCoeff();
        Eigen::VectorXd singular_values_inv = (singular_values.array() > tolerance).select(singular_values.array().inverse(), 0.0);

        Eigen::MatrixXd Jb_dagger = svd.matrixV() * singular_values_inv.asDiagonal() * svd.matrixU().transpose();
        Eigen::VectorXd Delta_theta = Jb_dagger * V_e;
        current_joints += Delta_theta;
    }

    return current_joints;

    //return j0;
}

std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> ScrewsKinematicsSolver::space_chain()
{
    return {m_m, m_screws};
}

//FINISHED: Implement body_chain(). You can obtain the variables to transform to body frame from space_chain().
std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> ScrewsKinematicsSolver::body_chain()
{
    auto [M_space, S_space] = space_chain();
    Eigen::Matrix4d M_body = M_space.inverse();
    Eigen::MatrixXd Ad_M_inv = utility::adjoint_matrix(M_body);
    std::vector<Eigen::VectorXd> B_body;
    B_body.reserve(S_space.size());

    for (const auto& S_i : S_space)
    {
        Eigen::VectorXd B_i = Ad_M_inv * S_i;
        B_body.push_back(B_i);
    }
    return {M_body, B_body};
}

//FINISHED: Implement space_jacobian() using space_chain()
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
}

//FINISHED: Implement body_jacobian() using body_chain()
Eigen::MatrixXd ScrewsKinematicsSolver::body_jacobian(const Eigen::VectorXd &current_joint_positions)
{
        Eigen::MatrixXd Js = space_jacobian(current_joint_positions);
        Eigen::Matrix4d T_sb = fk_solve(current_joint_positions);
        Eigen::Matrix4d T_bs = T_sb.inverse();
        Eigen::MatrixXd Ad_T_bs = utility::adjoint_matrix(T_bs);
        Eigen::MatrixXd Jb = Ad_T_bs * Js;
        return Jb;
}
