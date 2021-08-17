/**
* Kalman filter implementation using Eigen. Based on the following
* introductory paper:
*
*     http://www.cs.unc.edu/~welch/media/pdf/kalman_intro.pdf
*
* @author: Hayk Martirosyan
* @date: 2014.11.15
*/

#include <Eigen/Dense>
#include "physics.hh"
#include "globals.hh"
#include <iostream>
#include <fstream>

#pragma once

class KalmanFilter
{

public:
  /**
  * Create a Kalman filter with the specified matrices.
  *   A - System dynamics matrix
  *   C - Output matrix
  *   Q - Process noise covariance
  *   R - Measurement noise covariance
  *   P - Estimate error covariance
  */
  KalmanFilter(
      double dy,
      const Eigen::MatrixXd &A,
      const Eigen::MatrixXd &C,
      const Eigen::MatrixXd &Q,
      const Eigen::MatrixXd &R,
      const Eigen::MatrixXd &P);

  KalmanFilter();

  void init_gain(const Eigen::VectorXd &x0, std::vector<physics::digi_hit *> first_layer);
  void init_means(const Eigen::VectorXd x0, const Eigen::VectorXd q, const Eigen::MatrixXd B, const Eigen::MatrixXd D);
  void init_smooth_gain();
  void init_smooth_means();

  // find index of nearest hit to position
  std::vector<int> find_nearest(std::vector<physics::digi_hit *> hits, Eigen::VectorXd position)
  {
    std::ofstream file;
    file.open("print.txt", std::ios_base::app);

    int min_index = -1;
    int second_min_index = -1;
//    std::vector<int> chi_inds = {};

    int j = 0;

    double min_val;

    if (initialized)
      min_val = cuts::kalman_chi_add;
    else
      min_val = 1e6; // if the filter hasn't been initialised take hit
                     // with lowest chi regardless of the value

    // R_k^{k-1} in Fruhwirth
//    Eigen::MatrixXd err_metric = R + C * (A * P * A.transpose() + Q) * C.transpose();
    Eigen::MatrixXd err_metric = R + C * P * C.transpose();
//    file << " Err metric is " << err_metric << std::endl;

//    file << "chi's and hits are " << std::endl;

    for (auto hit : hits)
    {
      Eigen::VectorXd hit_eig(3);
      hit_eig << hit->x, hit->t, hit->z;

      double del_chi = (hit_eig - C * position).transpose() * err_metric.inverse() * (hit_eig - C * position);

      file << del_chi << ", and hit is " << hit_eig.transpose() << std::endl;

      Eigen::VectorXd x_temp(6);
      x_temp = position + K * (hit_eig - C * position);

      double v = 0;
      for (int i=3; i < x_temp.size(); i++) v += std::pow(x_temp[i], 2);
      v = std::sqrt(v);

      file << "speed is " << v << std::endl;
      file << "beta is " << v / constants::c << std::endl;

      if (del_chi < min_val && cuts::kalman_v_add[0] < v / constants::c && v / constants::c < cuts::kalman_v_add[1])
      {
        file << "using this hit " << std::endl;

        second_min_index = min_index; // keep track of second lowest for king moves
        min_index = j;
//        chi_inds.insert(chi_inds.begin(),j); // store all inds that pass the cuts (could have made the track)
        min_val = del_chi;
      }
      j++;
    }
    //file << std::endl;

    file << "layer index is " << min_index << std::endl;
    file << "King moves layer index is " << second_min_index << std::endl;

    if (min_index != -1)
      file << "King moves index is " << hits[min_index]->index << std::endl;
//    if (chi_inds.size() > 0)
//      file << "hit index is " << hits[chi_inds[0]]->index << std::endl;

    if (second_min_index != -1)
      file << "King moves index is " << hits[second_min_index]->index << std::endl;

    file.close();

//    return chi_inds;
    return {min_index, second_min_index};
  }


  // find index of nearest track to current vertex best estimate
  int find_nearest_vertex(std::vector<physics::track *> tracks)
  {
//    std::cout << "kalaman hh " << std::endl;

    std::ofstream file;
    file.open("print.txt", std::ios_base::app);

    int min_index = -1;
    int j = 0;

    double min_val;
//    double min_val = 1e8; // no minumum value for acceptance, we only use beta to cut

    if (initialized)
      min_val = cuts::kalman_vertex_chi_add;
    else
      min_val = 1e8; // if the filter hasn't been initialised take hit
                     // with lowest chi regardless of the value

   // std::cout << "B is " << B_f.back() << std::endl;
   // std::cout << "C is " << C << std::endl;


    //std::cout << "G is " << G << std::endl;
    //std::cout << "W is " << W << std::endl;


    //std::cout << "G_B is " << G_B << std::endl;
    //std::cout << "P_new is " << P_new << std::endl;

    //std::cout << "x_hat is " << x_hat.transpose() << std::endl;

    for (auto track : tracks)
    {
      //this->update_matrices_means(track);

      file << "NEW TRACK" << std::endl;

  //    std::cout << "kalman -" << std::endl;

      Eigen::MatrixXd R = track->P_s;

      double dt = x_f.back()[3] - track->t0;

//      std::cout << "kalman -" << std::endl;

      Eigen::MatrixXd B(6, 3);
      B << dt, 0, 0,
           0, dt / track->vy, 0,
           0, 0, dt,
           1, 0, 0,
           0, 1, 0,
           0, 0, 1;

//      file << " B is \n" << B << std::endl;

      Eigen::MatrixXd G = R.inverse(); // R is V in Fruhwirth
//      Eigen::MatrixXd W = (B_f.back().transpose() * G * B_f.back()).inverse();
      Eigen::MatrixXd W = (B.transpose() * G * B).inverse();

//      file << "G is \n" << G << std::endl;
//      file << "W is \n" << W << std::endl;

//      Eigen::MatrixXd G_B = G - G * B_f.back() * W * B_f.back().transpose() * G;
//      Eigen::MatrixXd P_new = (P_f.back().inverse() + C.transpose() * G_B * C).inverse();
      Eigen::MatrixXd G_B = G - G * B * W * B.transpose() * G;
      Eigen::MatrixXd P_new = (P.inverse() + C.transpose() * G_B * C).inverse();

//      file << "G_B is \n" << G_B << std::endl;
      file << "P_new is \n" << P << std::endl;

      Eigen::VectorXd y(6);
      y << track->x0, track->t0, track->z0, track->vx, track->vy, track->vz;

      //std::cout << "kalman -" << std::endl;

      file << " y is " << y.transpose() << std::endl;

//      Eigen::VectorXd x_hat_new = P_new * (P_f.back().inverse() * x_hat + C.transpose() * G_B * y);
//      Eigen::VectorXd q = W * B_f.back().transpose() * G * (y - C * x_hat_new);
      Eigen::VectorXd x_hat_new = P_new * (P.inverse() * x_hat + C.transpose() * G_B * y);
      Eigen::VectorXd q = W * B.transpose() * G * (y - C * x_hat_new);

      file << " x_hat_new is " << x_hat_new.transpose() << std::endl;
      file << " q is " << q.transpose() << std::endl;

    //  std::cout << "B is \n" << B << std::endl;
//      std::cout << "kalaman hh 0 \n" << "W\n" << W << std::endl << "G\n" << G << std::endl << "C\n" <<
//		C << std::endl << "P_new\n"<< P_new << std::endl;

//      std::cout << "part 1 " << std::endl << W * B.transpose() * G * C * P_new << std::endl;
//      std::cout << "part 2 " << std::endl << C.transpose() * G * B * W << std::endl;
//      std::cout << "part 3 " << std::endl << W * B.transpose() * G * C * P_new * C.transpose() * G * B * W << std::endl;

      Eigen::MatrixXd D = W + W * B.transpose() * G * C * P_new * C.transpose() * G * B * W;

  //    std::cout << "kalaman hh 1 " << D << std::endl;
/*
      double dt_v = x_hat_new[3] - track->t0; // difference between initial and vertex time
      file << "vertex time " << x_hat_new[3] - track->t0 << std::endl;;
      file << "predicted position at vertex time " << track->x0 + track->vx * dt_v << ", " <<
						      track->y0 + track->vy * dt_v << ", " <<
						      track->z0 + track->vz * dt_v << std::endl;
*/
//      Eigen::VectorXd res_t = (y - C * x_hat_new - B_f.back() * q); // track residual
      Eigen::VectorXd res_t = (y - C * x_hat_new - B * q); // track residual
      Eigen::VectorXd res_v = (x_hat_new - x_hat); // vertex residual



      file << "res_t is " << res_t.transpose() << std::endl;
      file << "vertex residual is " << res_v.transpose() << "\n for vertex " << x_hat_new.transpose() << std::endl;

      double del_chi = res_t.transpose() * G * res_t;
//      del_chi = del_chi + res_v.transpose() * P_f.back().inverse() * res_v;
      del_chi = del_chi + res_v.transpose() * P.inverse() * res_v;

      file << "chi increment is " << del_chi << std::endl;
      file << "Beta is " << q.norm() / constants::c << ", for velocity " << q.transpose() << std::endl;

//      double v_2 =  q.squaredNorm();

//      std::cout << "kalman 2" << std::endl;

      double sigma_v = (q.transpose() * D * q)(0) / q.squaredNorm(); // error in v
//      sigma_v = sigma_v / v_2;
      //std::cout << "kalman 3 " << sigma_v << std::endl;

      pulls_v_f.push_back((q.norm() - constants::c) / sigma_v);

    //  std::cout << "kalman 4 " << pull_v << std::endl;

      file << "pull is " << pulls_v_f.back() << std::endl;
//      file << "for D \n" << D << std::endl;

//      if (del_chi < min_val && std::abs(pulls_v_f.back()) < kalman::pull_cut_add)
      if (del_chi < min_val
          && kalman::v_cut_add[0] < q.norm() / constants::c && q.norm() / constants::c < kalman::v_cut_add[1])
      {
        min_index = j;
        min_val = del_chi;
      }

  //    std::cout << "kalman 5 " << std::endl;

      j++;
    }
    //file << std::endl;

    file << "minimum index is " << min_index << std::endl;
    if (min_index != -1)
      file << "track index is " << tracks[min_index]->index << std::endl;

    file.close();
//    std::cout << "kalaman hh " << std::endl;

    return min_index;
  }

  double update_gain(const std::vector<physics::digi_hit *> y);
  double update_gain(const std::vector<physics::digi_hit *> y, double dy);
  double update_means(const std::vector<physics::track*> tracks);
  double smooth_gain(const physics::digi_hit *y, int k); //, bool dropping);
  double smooth_means(int k);

  Eigen::VectorXd state() { return x_hat; };
  //double time() { return t; };
  Eigen::MatrixXd cov_matrix() { return P; };

  void push_q(Eigen::VectorXd q) { q_f.push_back(q); };

  std::vector<Eigen::VectorXd> x_f_list() { return x_f; };
  std::vector<Eigen::VectorXd> x_p_list() { return x_p; };

  std::vector<Eigen::VectorXd> x_s; // smoothed State
  std::vector<Eigen::MatrixXd> P_s; // smoothed Cov
  std::vector<Eigen::VectorXd> q_s; // smooth vertex velocity

  Eigen::VectorXd x_n;
  Eigen::MatrixXd P_n;

  // record index of hits added by find_nearest()
  std::vector<int> added_inds;
  std::vector<int> king_move_inds;
  std::vector<double> pulls_v_f;
  std::vector<double> pulls_v_s;

  std::vector<physics::digi_hit *> added_hits;
  std::vector<physics::digi_hit *> unadded_hits;
  std::vector<physics::track *> added_tracks;
  std::vector<physics::track *> unadded_tracks;

//  std::ofstream file;

  // chi square value
  double chi;

  // Discrete height step
  double dy;

  // run options
  bool dropping;
  //  bool adding;
  bool alrdy_drop;


  // scattering per y m
  double x_scat;
  double z_scat;

//  Eigen::MatrixXd B, R;

private:
  // Matrices for computation
  Eigen::MatrixXd A, C, Q, P, K, P0, err_metric, B, D, R;

  // System dimensions
  int m, n;

  // Initial and current time
  double t0, t;

  // Is the filter initialized?
  bool initialized;

  // n-size identity
  Eigen::MatrixXd I;

  // Estimated states
  Eigen::VectorXd x_hat, x_hat_new;

  // Storage for update() data
  std::vector<Eigen::VectorXd> x_p; // Prediction State
  std::vector<Eigen::VectorXd> x_f; // Filter State
  std::vector<Eigen::VectorXd> y_v; // Vertexer Track Parameters
  std::vector<Eigen::MatrixXd> P_p; // Prediction Cov
  std::vector<Eigen::MatrixXd> P_f; // Filter Cov
  std::vector<Eigen::MatrixXd> A_k; // Action Matrix
  std::vector<Eigen::VectorXd> q_f; // Filter State

  std::vector<Eigen::MatrixXd> W_f;
  std::vector<Eigen::MatrixXd> G_f;
  std::vector<Eigen::MatrixXd> G_B_f;
  std::vector<Eigen::MatrixXd> B_f;
  std::vector<Eigen::MatrixXd> C_f;
  std::vector<Eigen::MatrixXd> D_f;
  std::vector<Eigen::MatrixXd> D_s;

  //Storage for smooth() data
  std::vector<double> chi_2_s; // smoothed chi^2

  void update_matrices(physics::digi_hit *a_hit);
  void update_matrices_means(const physics::track* new_track);

  void Q_update(double dy, double a, double b, double c);
  void Q_propagate(const physics::track* track);
};
