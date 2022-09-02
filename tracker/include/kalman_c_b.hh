/**
* Kalman filter implementation using Eigen. Based on the following
* introductory paper:
*
*     http://www.cs.unc.edu/~welch/media/pdf/kalman_intro.pdf
*
* @author: Hayk Martirosyan
* @date: 2014.11.15
*/
#ifndef KM_C_B_DEFINE
#define KM_C_B_DEFINE

#include <Eigen/Dense>
#include "physics.hh"
#include "globals.hh"
#include <iostream>
#include <fstream>
#include "par_handler.hh"
#include "statistics.hh"

#pragma once

class KalmanFilter_c_b
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
  KalmanFilter_c_b(
      double dy,
      const Eigen::MatrixXd &A,
      const Eigen::MatrixXd &C,
      const Eigen::MatrixXd &Q,
      const Eigen::MatrixXd &R,
      const Eigen::MatrixXd &P);

  KalmanFilter_c_b();

  void init_gain(const Eigen::VectorXd &x0, std::vector<physics::digi_hit *> first_layer);
  void init_means(const Eigen::VectorXd x0, const Eigen::VectorXd q, const Eigen::MatrixXd B, const Eigen::MatrixXd D);
  void init_smooth_gain();
  void init_smooth_means();

  std::vector<int> find_nearest(std::vector<physics::digi_hit *> hits, Eigen::VectorXd position)
  { // find index of nearest hit to position

    int min_index = -1;
    int second_min_index = -1;

    int j = 0;

    double min_val;
    bool first_hit;
    if (initialized) {
      min_val = par_handler->par_map["kalman_chi_add"];
//      min_val = par_handler->par_map["kalman_pval_add"];
//      min_val = cuts::kalman_chi_add;
      first_hit = false;
    }

    else {
      min_val = 1e6; // if the filter hasn't been initialised take hit
                     // with lowest chi regardless of the value
      first_hit = true;
    }

    Eigen::MatrixXd P_temp = P;

    Stat_Funcs sts;

//    Eigen::MatrixXd err_metric = R + C * P * C.transpose();

    for (auto hit : hits)
    {
      Eigen::VectorXd hit_eig(3);
      hit_eig << hit->x, hit->t, hit->z;
	  std::cout << "find_nearest: update_matrices" << std::endl;
      update_matrices(hit);

//      Eigen::VectorXd x_hat_new_temp = A * x_hat;
      Eigen::VectorXd x_hat_new_temp = A * position;

      if (initialized) {
        P_temp = A * P * A.transpose() + Q;
		std::cout << "find_nearest: P_temp:" << std::endl;
		std::cout << P_temp << std::endl;
//        Eigen::MatrixXd K_temp = P_temp * C.transpose() * (C * P_temp * C.transpose() + R).inverse();
      }
	  //TODO: Make sure the math here is still valid after the coordinate change
      Eigen::MatrixXd err_metric = R + C * P_temp * C.transpose();
	  std::cout << "find_nearest: err_metric" << std::endl;
      // compute the chi increment
//      double del_chi = (hit_eig - C * position).transpose() * err_metric.inverse() * (hit_eig - C * position);
      double del_chi = (hit_eig - C * x_hat_new_temp).transpose() * err_metric.inverse() * (hit_eig - C * x_hat_new_temp);

      /*
      // calculate p value from the chi increment
      double ndof = added_hits.size() + 1;
      ndof = 4.0 * ndof - 6.0 ? ndof > 1.0 : 1.0;
      //del_chi = sts.chi_prob(del_chi, ndof);
      del_chi = sts.chi_prob_eld(del_chi, ndof);
      */

      if (first_hit) {
        min_val = del_chi; // if not initialized, overwrite min_val with value of first hit (guarantees one of the hits gets used)
        first_hit = false;
        min_index = j;
      }

      if (del_chi < min_val) // && cuts::kalman_v_add[0] < v / constants::c && v / constants::c < cuts::kalman_v_add[1])
      { // if hit has lowest chi so far and meets beta cuts, keep track of its index and the previous best

        // keep track of second lowest for king moves
        second_min_index = min_index;
        min_index = j;

        min_val = del_chi;
      }
      j++;
    }

    return {min_index, second_min_index};
  }

  int find_nearest_vertex(std::vector<physics::track *> tracks)
  { // find index of nearest track to current vertex best estimate
    // CURRENTLY NOT IN USE

    int min_index = -1;
    int j = 0;

    double min_val;

    if (initialized)
//      min_val = cuts::kalman_vertex_chi_add;
      min_val = par_handler->par_map["kalman_vertex_chi_add"];
    else
      min_val = 1e8; // if the filter hasn't been initialised take hit
                     // with lowest chi regardless of the value

    for (auto track : tracks)
    {
      Eigen::MatrixXd R = track->P_s;

      double dt = x_f.back()[3] - track->t0;

      Eigen::MatrixXd B(6, 3);
      B << dt, 0, 0,
          0, dt / track->vy, 0,
          0, 0, dt,
          1, 0, 0,
          0, 1, 0,
          0, 0, 1;

      Eigen::MatrixXd G = R.inverse(); // R is V in Fruhwirth
      Eigen::MatrixXd W = (B.transpose() * G * B).inverse();
      Eigen::MatrixXd G_B = G - G * B * W * B.transpose() * G;

      // state covariance matrix
      Eigen::MatrixXd P_new = (P.inverse() + C.transpose() * G_B * C).inverse();

      Eigen::VectorXd y(6);
      y << track->x0, track->t0, track->z0, track->vx, track->vy, track->vz;

      // filtered state and vertex velocity
      Eigen::VectorXd x_hat_new = P_new * (P.inverse() * x_hat + C.transpose() * G_B * y);
      Eigen::VectorXd q = W * B.transpose() * G * (y - C * x_hat_new);

      // vertex velocity covariance matrix
      Eigen::MatrixXd D = W + W * B.transpose() * G * C * P_new * C.transpose() * G * B * W;

      // track and vertex residuals
      Eigen::VectorXd res_t = (y - C * x_hat_new - B * q);
      Eigen::VectorXd res_v = (x_hat_new - x_hat);

      double del_chi = res_t.transpose() * G * res_t;
      del_chi = del_chi + res_v.transpose() * P.inverse() * res_v;

      double sigma_v = (q.transpose() * D * q)(0) / q.squaredNorm(); // error in v
      pulls_v_f.push_back((q.norm() - constants::c) / sigma_v);      // pull of v from speed of light

      //      if (del_chi < min_val && std::abs(pulls_v_f.back()) < kalman::pull_cut_add)
//      if (del_chi < min_val && kalman::v_cut_add[0] < q.norm() / constants::c && q.norm() / constants::c < kalman::v_cut_add[1])
      if (del_chi < min_val && par_handler->par_map["v_cut_add[0]"] < q.norm() / constants::c
	  && q.norm() / constants::c < par_handler->par_map["v_cut_add[1]"])
      {
        min_index = j;
        min_val = del_chi;
      }
      j++;
    }
    return min_index;
  }

  Eigen::VectorXd to_cartesian_v(double tanx, double tanz) {

    Eigen::VectorXd vel(3);
	vel << tanx, 1, tanz;
	double N = std::sqrt(1 + tanx*tanx + tanz*tanz);
    vel = vel * beta * constants::c / N;

    return vel;
  }
/*  Eigen::VectorXd to_cartesian_v(double theta, double phi) {

    Eigen::VectorXd vel(3);
	vel << 
    vel << std::sin(theta) * std::cos(phi), std::cos(theta), std::sin(theta) * std::sin(phi);
    vel = vel * beta * constants::c;

    return vel;
 }
*/
  double update_gain(const std::vector<physics::digi_hit *> y);
  double update_gain(const std::vector<physics::digi_hit *> y, double dy);
  double update_means(const std::vector<physics::track *> tracks);

  void king_moves_algorithm(const std::vector<physics::digi_hit *> y_list, std::vector<int> hit_inds);

  double smooth_gain(const physics::digi_hit *y, int k);
  double smooth_means(int k);

  Eigen::VectorXd state() { return x_hat; };
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

  // chi square value
  double chi;

  // Discrete height step
  double dy;
  double y_val;

  // run options
  bool dropping;
  bool alrdy_drop;

  // scattering per y m
  double x_scat;
  double z_scat;

  // constrained beta
  double beta;

  ParHandler* par_handler;

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

  std::vector<double> chi_2_s; // smoothed chi^2

  void update_matrices(physics::digi_hit *a_hit);
  void update_matrices_means(const physics::track *new_track);

  void Q_update(double dy, double a, double b, double c);
  void Q_propagate(const physics::track *track);
};


#endif
