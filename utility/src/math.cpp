#include "utility/math.h"
#include "utility/vectors.h"
#include <numbers>

namespace AIS4104::utility {
    //Equation (x) page x, MR pre-print 2019 digital
    //Page X MR physical book Equation (X)
    //Page 577, MR pre-print 2019 digital
    Eigen::Vector3d euler_zyx_from_rotation_matrix(const Eigen::Matrix3d &r)
    {
        double gamma, beta, alfa;
        if (r(2,0) != 1 || r(2,0) != -1  ) {
            gamma = atan2(r(2,1), r(2,2));                  // yaw (Z)
            beta = atan2(-r(2,0), sqrt(pow(r(0,0),2) + pow(r(1,0),2))); // pitch (Y)
            alfa = atan2(r(1,0), r(0,0));                  // roll (X)
        }
        else if (r(2,0) == -1) {
            beta = std::numbers::pi/2;
            alfa = 0;
            gamma = atan2(r(0,1), r(1,1));

        }
        else if (r(2,0) == 1 ) {
            beta = -std::numbers::pi/2;
            alfa = 0;
            gamma = -atan2(r(0,1), r(1,1));
        }


        Eigen::Vector3d eulerZYX(alfa, beta, gamma);
        return eulerZYX;
    }
//Equation (3.30) page 75, MR pre-print 2019 digital
    Eigen::Matrix3d skew_symmetric(const Eigen::Vector3d &v)
    {
        Eigen::Matrix3d S;
        S <<      0.0,     -v.z(),   v.y(),
               v.z(),      0.0,     -v.x(),
              -v.y(),     v.x(),   0.0;
        return S;
    }
    //new: page 65 MR physical book Equation (3.30)
    Eigen::Vector3d from_skew_symmetric(const Eigen::Matrix3d &m)
    {
        Eigen::Vector3d v;
        v << m(2,1), m(0,2), m(1,0);
        return v;
    }
    //new: Equation (3.20) page 98, MR pre-print 2019 digital
    Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p)
    {
        Eigen::Matrix3d p_hat=skew_symmetric(p);
        Eigen::Matrix3d p_hat_times_R = p_hat * r;
        Eigen::MatrixXd Ad(6, 6);
        Ad.setZero();
        Ad.block<3, 3>(0, 0) = r;
        Ad.block<3, 3>(3, 0) = p_hat_times_R;
        Ad.block<3, 3>(3, 3) = r;

        return Ad;
    }
    //Equation (3.20) page 98, MR pre-print 2019 digital
    Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix4d &tf)
    {
        Eigen::Matrix3d R = tf.block<3,3>(0,0);
        Eigen::Vector3d p = tf.block<3,1>(0,3);

        Eigen::Matrix3d p_hat=skew_symmetric(p);

        Eigen::MatrixXd AdT(6,6);
        AdT.setZero();
        AdT.block<3,3>(0,0) = R;
        AdT.block<3,3>(3,3) = R;
        AdT.block<3,3>(3,0) = p_hat * R;

        return AdT;
    }
    //Equation (3.20) page 98, MR pre-print 2019 digital
    Eigen::VectorXd adjoint_map(const Eigen::VectorXd &twist, const Eigen::Matrix4d &tf)
    {
        Eigen::MatrixXd Ad_g = adjoint_matrix(tf);
        Eigen::VectorXd transformed_twist = Ad_g * twist;

        return transformed_twist;
    }
    //Equation (3.70) page 96, MR pre-print 2019 digital
    Eigen::VectorXd twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v)
    {
        Eigen::VectorXd V(6);
        V << w, v;
        return V;
    }
    //new: Page 87 MR physical book
    Eigen::VectorXd twist(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h, double angular_velocity)
    {
        double dot_theta = angular_velocity;
        Eigen::Vector3d omega = s * dot_theta;
        Eigen::Vector3d term1 = -omega.cross(q);
        Eigen::Vector3d term2 = s * h * dot_theta;
        Eigen::Vector3d v = term1 + term2;
        Eigen::VectorXd V(6);
        V << omega,
                 v;

        return V;
    }
    //new: page 83 MR physical book Equation (3.71)
    Eigen::Matrix4d twist_matrix(const Eigen::Vector3d &w, const Eigen::Vector3d &v)
    {
        Eigen::Matrix3d omega_skew = skew_symmetric(w);
        Eigen::Matrix4d V_matrix = Eigen::Matrix4d::Zero();
        V_matrix.block<3, 3>(0, 0) = omega_skew;
        V_matrix.block<3, 1>(0, 3) = v;
        return V_matrix;
    }
    //new: page 83 MR physical book Equation (3.71)
    Eigen::Matrix4d twist_matrix(const Eigen::VectorXd &twist)
    {
        Eigen::Vector3d w = twist.segment<3>(0);
        Eigen::Vector3d v = twist.segment<3>(3);

        return twist_matrix(w,v);
    }
    //new: page 88 MR physical book --> note to self: samesame som dei andre sånn sett
    Eigen::VectorXd screw_axis(const Eigen::Vector3d &w, const Eigen::Vector3d &v)
    {
        Eigen::VectorXd S(6);
        S << w, v;
        return S;
    }
    //Page 102, MR pre-print 2019 digital
    Eigen::VectorXd screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h)
    {
        Eigen::VectorXd S(6);
        Eigen::Vector3d v = -s.cross(q) + h*s;
        S << s, v;
        return S;
    }
    //Equation (3.51) page 82, MR pre-print 2019 digital
    Eigen::Matrix3d matrix_exponential(const Eigen::Vector3d &w, double theta)
    {
        double rad = theta;
        Eigen::Matrix3d R = Eigen::Matrix3d::Identity() + sin(rad) * skew_symmetric(w) + (1-cos(rad)) * skew_symmetric(w) * skew_symmetric(w);
        return R;
    }
    //Equation (3.25) page 103, MR pre-print 2019 digital
    Eigen::Matrix4d matrix_exponential(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double theta)
    {
        if (w.norm() < 1e-6)
        {
            return transformation_matrix(Eigen::Matrix3d::Identity(), v * theta);
        }
        else
        {
            double rad = theta;
            Eigen::Matrix3d w_skew = skew_symmetric(w);
            Eigen::Matrix3d R = matrix_exponential(w, theta);
            Eigen::Matrix3d G = (Eigen::Matrix3d::Identity() * rad) +
                                ((1 - cos(rad)) * w_skew) +
                                ((rad - sin(rad)) * (w_skew * w_skew));
            Eigen::Vector3d v_theta = G * v;
            return transformation_matrix(R, v_theta);
        }
    }
    //new: Equation (3.25) page 103, MR pre-print 2019 digital
    Eigen::Matrix4d matrix_exponential(const Eigen::VectorXd &screw, double theta)
    {
        Eigen::Vector3d w = screw.segment<3>(0);
        Eigen::Vector3d v = screw.segment<3>(3);
        Eigen::Matrix4d T = matrix_exponential(w, v, theta);
        return T;
    }
    //Algorithm page 85 and 86, MR pre-print 2019 digital
    std::pair<Eigen::Vector3d, double> matrix_logarithm(const Eigen::Matrix3d &r)
    {
        double theta;
        Eigen::Vector3d w;

        if (r == Eigen::Matrix3d::Identity()) {
            theta = 0;
        }
        else {
            double trace_r = r.trace();
            if (trace_r == -1) {
                theta = std::numbers::pi;
                w = (1 / sqrt(2 * (1 + r(2,2)))) * Eigen::Vector3d (r(0,2), r(1,2), 1 + r(2,2));
            }
            else {
                theta = acos(0.5 * (trace_r - 1));
                double w_n = 1 / (2*sin(theta));

                double w_1 = w_n * (r(2, 1) - r(1, 2));
                double w_2 = w_n * (r(0, 2) - r(2, 0));
                double w_3 = w_n * (r(1, 0) - r(0, 1));

                w = Eigen::Vector3d(w_1, w_2, w_3);
            }
        }

        return std::pair<Eigen::Vector3d, double>(w, theta);
    }
    //new: page 90 MR physical book Algorithm on this page
    std::pair<Eigen::VectorXd, double> matrix_logarithm(const Eigen::Matrix3d &r, const Eigen::Vector3d &p)
    {
    const double epsilon = 1e-6;
    Eigen::Vector3d w;
    Eigen::Vector3d v;
    double theta;

    if (r.isApprox(Eigen::Matrix3d::Identity(), epsilon))
    {
        w = Eigen::Vector3d::Zero();
        theta = p.norm();

        if (theta < epsilon) {
            v = Eigen::Vector3d::Zero();
        } else {
            v = p / theta;
        }
    }
    else
    {
        auto log_R_pair = matrix_logarithm(r);
        w = log_R_pair.first;
        theta = log_R_pair.second;
        Eigen::Matrix3d w_skew = skew_symmetric(w);
        Eigen::Matrix3d G_inv = Eigen::Matrix3d::Identity();
        G_inv -= 0.5 * w_skew;
        if (std::abs(theta) > epsilon) {
            double term_scalar = (1.0 / theta) - (0.5 * (1.0 / std::tan(theta / 2.0)));
            G_inv += term_scalar * w_skew * w_skew;
        }
        v = G_inv * p;
    }
    Eigen::VectorXd S(6);
    S << w, v;

    return std::make_pair(S, theta);
    }
    //--> note to self: assignment 2
    double cot(double x) {
        return 1 / (std::sin(x) / std::cos(x)); // fordi tan = sin/cos
    }
    //--> note to self:Is also from assignment 2
    //new: page 90 MR physical book Algorithm on this page
    std::pair<Eigen::VectorXd, double> matrix_logarithm(const Eigen::Matrix4d &tf)
    {
        const Eigen::Matrix3d R = tf.topLeftCorner<3,3>();
        const Eigen::Vector3d p = tf.topRightCorner<3,1>();

        Eigen::Vector3d w;
        Eigen::Vector3d v;
        double theta;

        if (R == Eigen::Matrix3d::Identity()) {
            w = Eigen::Vector3d::Zero();
            v = p / p.norm();
            theta = p.norm();
        }
        else {
            std::pair<Eigen::Vector3d, double> m_log = matrix_logarithm(R);
            w = m_log.first;
            theta = m_log.second;
            const Eigen::Matrix3d skew_w = skew_symmetric(w);
            v = (((1/theta) * Eigen::Matrix3d::Identity()) - 0.5 * skew_w + ((1/theta) - 0.5* cot(theta/2)) * skew_w*skew_w) * p;
        }

        return std::pair<Eigen::VectorXd, double>(twist(w,v), theta);
        //return std::make_pair(Eigen::VectorXd::Zero(6),0);

    }

    //Equation found in page 72, MR pre-print 2019 digital
    Eigen::Matrix3d rotate_x(double radians)
    {
        Eigen::Matrix3d matrix;
        double theta = radians;

        matrix << 1,          0,           0,
                  0,  std::cos(theta), -std::sin(theta),
                  0,  std::sin(theta),  std::cos(theta);

        return matrix;
    }
    //Equation found in page 72, MR pre-print 2019 digital
    Eigen::Matrix3d rotate_y(double radians)
    {
        Eigen::Matrix3d matrix;
        double theta = radians;

        matrix <<  std::cos(theta),  0,  std::sin(theta),
                   0,                1,  0,
                  -std::sin(theta),  0,  std::cos(theta);
        return matrix;
    }
    //Equation found in page 72, MR pre-print 2019 digital
    Eigen::Matrix3d rotate_z(double radians)
    {
        Eigen::Matrix3d matrix;
        double theta = radians;


        matrix <<  std::cos(theta), -std::sin(theta), 0,
                   std::sin(theta),  std::cos(theta), 0,
                   0,                0,               1;

        return matrix;
    }
    //Found in page 67, MR pre-print 2019 digital
    Eigen::Matrix3d rotation_matrix_from_frame_axes(const Eigen::Vector3d &x, const Eigen::Vector3d &y, const Eigen::Vector3d &z)
    {

        Eigen::Matrix3d Matrix;

        Matrix.col(0) = x.normalized();
        Matrix.col(1) = y.normalized();
        Matrix.col(2) = z.normalized();

        return Matrix;
    }
    //Found in page 577, MR pre-print 2019 digital
    Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e)
    {
        double alpha = e(0); // rotasjon om z
        double beta  = e(1); // rotasjon om y
        double gamma = e(2); // rotasjon om x


        Eigen::Matrix3d Rz = rotate_z(alpha);
        Eigen::Matrix3d Ry = rotate_y(beta);
        Eigen::Matrix3d Rx = rotate_x(gamma);

        Eigen::Matrix3d R = Eigen::Matrix3d::Identity() * Rz * Ry * Rx;  // rekkefølgen Z * Y * X
        return R;
    }
    //Equation (3.51) page 82, MR pre-print 2019 digital
    Eigen::Matrix3d rotation_matrix_from_axis_angle(const Eigen::Vector3d &axis, double radians)
    {
        Eigen::Vector3d u = axis.normalized();
        double theta = radians;

        Eigen::Matrix3d u_skew = skew_symmetric(u);

        // Rodrigues rotasjonsformel
        Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
        Eigen::Matrix3d matrix =
            I * std::cos(theta) +
            (1 - std::cos(theta)) * (u * u.transpose()) +
            u_skew * std::sin(theta);

        return matrix;
    }
    //new: Page 87, MR pre-print 2019 digital
    Eigen::Matrix3d rotation_matrix(const Eigen::Matrix4d &tf)
    {
        return tf.block<3,3>(0,0);
    }
    //new: Page 87, MR pre-print 2019 digital
    Eigen::Matrix4d transformation_matrix(const Eigen::Vector3d &p)
    {
        Eigen::Matrix4d matrix = Eigen::Matrix4d::Identity();
        matrix.block<3,1>(0,3) = p;
        return matrix;
    }
    //new: Page 87, MR pre-print 2019 digital
    Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r)
    {
        Eigen::Matrix4d matrix = Eigen::Matrix4d::Identity();
        matrix.block<3,3>(0,0) = r;
        return matrix;
    }
    //new: Page 87, MR pre-print 2019 digital
    Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p)
    {
        Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
        T.block<3,3>(0,0) = r;   // øvre venstre 3x3 = R
        T.block<3,1>(0,3) = p;   // øvre høyre 3x1 = p
        return T;
    }
}