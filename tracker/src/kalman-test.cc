/**
* Implementation of KalmanFilter. Based on code by 
* 
* Hayk Martirosyan (2014.11.15)
* 
* and the following paper
* 
* Frühwirth, R. (1987). Application of Kalman filtering to track and vertex fitting.
* Nuclear Instruments and Methods in Physics Research Section A:
* Accelerators, Spectrometers, Detectors and Associated Equipment,
* 262(2-3), 444-450. doi:10.1016/0168-9002(87)90887-4
*
* @authors: Keegan Humphrey and Hao Liao
* @date: May 2021
*/
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <map>
#include "kalman.hh"
#include <unsupported/Eigen/MatrixFunctions>
#include "physics.hh"
#include "kalman-test.hh"
#include "globals.hh"
#include <fstream>
#include <algorithm>

void kalman_track::kalman_all(std::vector<physics::digi_hit *> trackhits, seed *current_seed)
{ // tracking algorithm using the kalman filter

  // sort hits by layer
  layer_sort(trackhits);

  // status == 1 => currently working
  status = 1;

  if (layers.size() < cuts::track_nlayers)
  {
    // status == -1 => not enough layers in the event
    status = -1;
  }

  while (status == 1)
  {

    init_seed_info(current_seed);

    init_matrices(current_seed);

    init_first_state();

    // initialise the KalmanFilter object for fitting
    // and assign initialised object to class variable
    KalmanFilter kf_init(0, A, C, Q, R, P);
    kf = kf_init;

    filter();

    for (int i = 1; i < kf.added_hits.size(); i++)
    {
      found_hits.push_back(kf.added_hits[i]);
    }

    if (status == 0)
      break;

    if (kf.x_f_list().size() < cuts::track_nlayers && !finding)
    {
      status = -2;
      break;
    }

    if (!finding)
    {
      smooth();
    }

    prepare_output();

    // status == 2 => seed succeeded
    status = 2;

  } // while loop
}

void kalman_track::init_matrices(seed *current_seed)
{ // initialise matrices for the filter

  A = Eigen::MatrixXd::Zero(n, n);
  C = Eigen::MatrixXd::Zero(m, n);
  Q = Eigen::MatrixXd::Zero(n, n);
  R = Eigen::MatrixXd::Zero(m, m);
  P = Eigen::MatrixXd::Zero(n, n);

  // Projection Matrix
  C << 1, 0, 0, 0, 0, 0,
      0, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0;

  // seed position vectors
  auto P1 = current_seed->hits.first->PosVector();
  auto P2 = current_seed->hits.second->PosVector();
  auto dr = P2 - P1;

  // differences between the vector components
  double dt = current_seed->hits.second->t - current_seed->hits.first->t;
  double dx = dr.x;
  double dy = dr.y;
  double dz = dr.z;

  // propagate velocity errors
  double v_x_err = 2.0 * dx / dt * std::sqrt(std::pow(first_hit->ex, 2) / (dx * dx) + 1.0 / (dt * dt));
  double v_y_err = 2.0 * dy / dt * std::sqrt(9.0 / (dy * dy * 12.0) + 1.0 / (dt * dt));
  double v_z_err = 2.0 * dz / dt * std::sqrt(std::pow(first_hit->ez, 2) / (dz * dz) + 1.0 / (dt * dt));

  Eigen::VectorXd errs(n);
  errs << first_hit->ex, first_hit->et, first_hit->ez, v_x_err, v_y_err, v_z_err;

  // track covariance matrix
  P << errs[0], 0, 0, 0, 0, 0,
      0, errs[1], 0, 0, 0, 0,
      0, 0, errs[2], 0, 0, 0,
      0, 0, 0, errs[3], 0, 0,
      0, 0, 0, 0, errs[4], 0,
      0, 0, 0, 0, 0, errs[5];
  P = P * P;
}

void kalman_track::init_first_state()
{ // init function to find the first state vector for the filter

  // find first hit by running the filter backwards from seed
  if (finding)
  {
    int next_layer;

    if (seed_was_used)
    { // if first seed hit was used linearly backpropagate to first layer with hits

      // get highest layer with hits below seed layer
      // or lowest layer with hits above seed layer
      next_layer = find_next_layer();

      x0 = find_guess(seedguess, seedguess[1] - layer_hits[next_layer][0]->y);

      first_hit_list = layer_hits[next_layer];
    }
    find_first();
  }
  else
  { // we are only fitting, not looking for hits

    // unambiguous since only 1 hit per layer on second filter
    lowest_hit = layer_hits[layers[0]][0];

    // pass previous fit velocity instead of seedguess in tr>
    velocity = {seedguess[3], seedguess[4], seedguess[5]};

    added_layers = layers;
    filter_start_layer = layers[0];
  }
}

void kalman_track::find_first()
{ // find the first hit to start the filter

  KalmanFilter kf_find_init(0, A, C, Q, R, P);
  kf_find = kf_find_init;

  kf_find.init_gain(x0, first_hit_list);

  // find layer index to start filter
  int start_ind;

  // layers index of first_hit_list
  if (seed_was_used)
    start_ind = det_ind_to_layers[(first_hit_list[0]->det_id).layerIndex];
  else
    start_ind = det_ind_to_layers[seed_layer];

  added_layers.push_back(layers[start_ind]);

  for (int i = start_ind; i > 0; i--)
  {

    double y_step = layer_hits[layers[i - 1]][0]->y - layer_hits[layers[i]][0]->y;

    if (skipped)
    { // if no hits are found add previous y_step to current one (skip a layer)
      y_step += kf_find.dy;
      skipped = false;
    }

    double chi = kf_find.update_gain(layer_hits[layers[i - 1]], y_step);

    if (chi == -1.0)
    {
      skipped = true;
      continue;
    }
    else
    {
      added_layers.push_back(layers[i - 1]);
    }
  }

  //reverse so bottom to top
  std::reverse(kf_find.added_hits.begin(), kf_find.added_hits.end());

  // take last hit in backward seed filter and use it to start forward filter
  found_hits = kf_find.added_hits;

  lowest_hit = kf_find.added_hits.back();
  velocity = {kf_find.x_f_list().back()[3], kf_find.x_f_list().back()[4], kf_find.x_f_list().back()[5]};

  filter_start_layer = layers[start_ind];
}

void kalman_track::filter()
{ // filter algorithm

  Eigen::VectorXd x_filter(6);
  x_filter << lowest_hit->x, lowest_hit->t, lowest_hit->z, velocity[0], velocity[1], velocity[2];

  std::vector<physics::digi_hit *> lowest_hit_list = {lowest_hit};

  kf.init_gain(x_filter, lowest_hit_list);
  kf.dropping = dropping;

  skipped = false;

  for (int i = det_ind_to_layers[filter_start_layer]; i < layers.size() - 1; i++)
  {

    double y_step = layer_hits[layers[i + 1]][0]->y - layer_hits[layers[i]][0]->y;

    if (skipped)
    {
      y_step += kf.dy;
      skipped = false;
    }

    double chi = kf.update_gain(layer_hits[layers[i + 1]], y_step);

    // no hit was found skip this layer
    if (chi == -1.0)
    {
      skipped = true;
      continue;
    }
    else
    {
      chi_f.push_back(chi);
      x_scat.push_back(kf.x_scat);
      z_scat.push_back(kf.z_scat);
    }
  }
}

void kalman_track::smooth()
{ // smoothing algorithm

  kf.init_smooth_gain();

  // last filter best estimate is unaffected by smoother (see Fruhwirth paper)
  chi_s = {chi_f.back()};

  for (int i = kf.added_hits.size() - 1; i > 0; i--)
  {
    double chi = kf.smooth_gain(kf.added_hits[i - 1], i - 1);

    chi_s.insert(chi_s.begin(), chi);
  }
}

void kalman_vertex::vertexer(std::vector<physics::track *> tracks_list, vertex_seed *seed)
{ // vertexer algorithm using kalman filter weighted means formalism

  int i = 0;

  status = 1;

  while (status == 1)
  {
    init_seed_info(seed);

    init_matrices(seed);

    KalmanFilter kfv_init(0, A, C, Q, R, P);
    kfv = kfv_init;
    kfv.dropping = dropping;

    filter(tracks_list);

    if (kfv.added_tracks.size() < 2)
    {
      status = 0;
      break;
    }

    smooth();

    added_tracks = kfv.added_tracks;
    unadded_tracks = kfv.unadded_tracks;

    prepare_output();

    status = 2;
  } // while loop
}

void kalman_vertex::init_matrices(vertex_seed *seed)
{

  A = Eigen::MatrixXd::Zero(m, m);
  C = Eigen::MatrixXd::Zero(m, s);
  Q = Eigen::MatrixXd::Zero(m, m);
  P = Eigen::MatrixXd::Zero(s, s);

  B = Eigen::MatrixXd::Zero(m, v);
  R = Eigen::MatrixXd::Zero(m, m);
  D = Eigen::MatrixXd::Zero(v, v);

  R = seed->tracks.first->P_s;

  // projects x = {x_v,y_v,z_v,t} onto p = {x_v,t_v,z_v,vx,vy,vz}
  // A in Fruhwirth paper (eq 24)
  C << 1, 0, 0, 0,
      0, 0, 0, 1,
      0, 0, 1, 0,
      0, 0, 0, 0,
      0, 0, 0, 0,
      0, 0, 0, 0;

  physics::track *tr1 = seed->tracks.first;
  physics::track *tr2 = seed->tracks.second;

  double dt = x0[3] - tr1->t0; // negative

  std::vector<double> dr = {tr2->x0 - tr1->x0, tr2->y0 - tr1->y0, tr2->z0 - tr1->z0};
  std::vector<double> dv = {tr1->vx - tr2->vx, tr1->vy - tr2->vy, tr1->vz - tr2->vz};

  // temporary variable
  std::vector<double> tmp;
  for (int i = 0; i < dr.size(); i++)
    tmp.push_back(dr[i] - 2 * dv[i] * x0[3]);

  // find norm of dv
  double v_2 = 0;
  for (auto di : dv)
    v_2 += di * di;

  Eigen::VectorXd q(3);
  //  q = (q0 + q1) / 2;
  q = q0;

  // jacobian for track to vertex parameter transformation
  // see overleaf doc for calculation
  Eigen::MatrixXd jac(4, 6);
  jac << v_2 + dv[0] * q[0], -q[0], dv[2] * q[0], v_2 * dt + q[0] * tmp[0], q[0] * tmp[1], q[0] * tmp[2],
      dv[0] * q[1], -q[1], dv[2] * dv[1], q[1] * tmp[0], v_2 * dt + q[1] * tmp[1], q[1] * tmp[2],
      dv[0] * q[2], -q[2], v_2 + dv[2] * dv[2], q[2] * tmp[0], q[2] * tmp[1], v_2 * dt + q[2] * tmp[2],
      dv[0], 0, dv[2], tmp[0], tmp[1], tmp[2];

  jac = jac / v_2;

  //  P = jac * (tr1->P_s + tr2->P_s) / 2 * jac.transpose();
  P = jac * tr1->P_s * jac.transpose();

  // positive
  dt = -dt;

  // projects q = {vx,vy,vz} onto p = {dx,dt,dz,vx,vy,vz}
  B << dt, 0, 0,
      0, dt / q0[1], 0,
      0, 0, dt,
      1, 0, 0,
      0, 1, 0,
      0, 0, 1;

  // vertex velocity covariance matrix
  D << R(3, 3), R(3, 4), R(3, 5),
      R(4, 3), R(4, 4), R(4, 5),
      R(5, 3), R(5, 4), R(5, 5);
}

void kalman_vertex::filter(std::vector<physics::track *> tracks_list)
{ // filter algorithm for the vertexing algorithm

  kfv.init_means(x0, q0, B, D);

  bool done = false;
  int j = 0;

  chi_v = 0;

  while (!done && tracks_list.size() != 0)
  {
    double chi = kfv.update_means(tracks_list);

    // no new tracks added
    if (chi == -1.0)
    {
      done = true;
      break;
    }
    else
    {
      tracks_list = kfv.unadded_tracks;

      // at most 1 added track
      added_tracks.push_back(kfv.added_tracks.back());

      chi_v += chi;
    }

    j++;
  }
}

void kalman_vertex::smooth()
{
  kfv.init_smooth_means();

  for (int j = kfv.added_tracks.size() - 1; j >= 0; j--)
  {
    kfv.smooth_means(j);
  }
}
