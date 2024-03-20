/*****************************************************************************
 * Multiscale Universal Interface Code Coupling Library                       *
 *                                                                            *
 * Copyright (C) 2019 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis, *
 *                    A. Skillen, W. Liu, S. Longshaw, O. Mahfoze             *
 *                                                                            *
 * This software is jointly licensed under the Apache License, Version 2.0    *
 * and the GNU General Public License version 3, you may use it according     *
 * to either.                                                                 *
 *                                                                            *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License");            *
 * you may not use this file except in compliance with the License.           *
 * You may obtain a copy of the License at                                    *
 *                                                                            *
 * http://www.apache.org/licenses/LICENSE-2.0                                 *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 *                                                                            *
 * ** GNU General Public License, version 3 **                                *
 *                                                                            *
 * This program is free software: you can redistribute it and/or modify       *
 * it under the terms of the GNU General Public License as published by       *
 * the Free Software Foundation, either version 3 of the License, or          *
 * (at your option) any later version.                                        *
 *                                                                            *
 * This program is distributed in the hope that it will be useful,            *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU General Public License for more details.                               *
 *                                                                            *
 * You should have received a copy of the GNU General Public License          *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
 *****************************************************************************/

/**
 * @file sampler_rbf.h
 * @author A. Skillen, W. Liu, S. Longshaw, O. Mahfoze
 * @date November 2018
 * @brief Spatial sampler using Radial Basis Function interpolation.
 */

#ifndef MUI_SAMPLER_RBF_H_
#define MUI_SAMPLER_RBF_H_

#include "../../general/util.h"
#include "../../uniface.h"
#include "../../linear_algebra/solver.h"
#include <iterator>
#include <ctime>
#include <sys/types.h>
#include <sys/stat.h>
#include <mpi.h>

namespace mui {

template<typename CONFIG = default_config,
typename O_TP = typename CONFIG::REAL, typename I_TP = O_TP>
class sampler_rbf {

public:
    using OTYPE = O_TP;
    using ITYPE = I_TP;
    using REAL = typename CONFIG::REAL;
    using INT = typename CONFIG::INT;
    using point_type = typename CONFIG::point_type;
    using EXCEPTION = typename CONFIG::EXCEPTION;

    static const bool QUIET = CONFIG::QUIET;
    static const bool DEBUG = CONFIG::DEBUG;

    /**
     * Input parameters:
     * 1. REAL r:
     *       The search radius used to construct each RBF
     * 2. std::vector<point_type>& pts:
     *       Vector of points that pre-set for RBF interpolation
     * 3. INT basisFunc:
     *       Parameter for basis function selection. Implemented functions are as follows:
     *       Gaussian(default): basisFunc_=0
     *       Wendland's C0:     basisFunc_=1
     *       Wendland's C2:     basisFunc_=2
     *       Wendland's C4:     basisFunc_=3
     *       Wendland's C6:     basisFunc_=4
     * 4. bool conservative:
     *       Switch for the mode of RBF interpolation:
     *       consistent mode(default): conservative=false
     *       conservative mode:        conservative=true
     * 5. bool smoothFunc:
     *       Switch for the smoothing function of the transformation matrix:
     *       without smoothing function(default): smoothFunc=false
     *       with smoothing function:             smoothFunc=true
     * 6. bool generateMatrix:
     *       Switch for whether to generate the transformation matrix:
     *       Don't generate the matrix (readRBFMatrix() should be used): generateMatrix=false
     *       Generate the matrix: generateMatrix=true
     * 7. const std::string& writeFileAddress:
     *       The address that the transformation matrix I/O uses.
     *       The default value of is an empty string - empty means no files will be written.
     *       The directory will be created on object construction if it doesn't already exist.
     * 8. REAL cutOff:
     *       Parameter to set the cut-off of the Gaussian basis function (only valid for basisFunc_=0).
     *       The default value of cutoff is 1e-9
     * 9. REAL cgSolveTol:
     *       The tolerance used to determine convergence for the ConjugateGradient solver
     *       The default value of cgSolveTol is 1e-6
     * 10. INT cgMaxIter:
     *       The maximum number of iterations each Eigen ConjugateGradient solve can take
     *       The default value of cgMaxIter is 0, which means the solver decides
     * 11. INT pouSize:
     *       The size of each partition used within the RBF-POU approach
     *       The default value of pouSize is 50, setting to 0 disables the partitioned approach
     * 12. INT precond:
     *       The Preconditioner of the Conjugate Gradient solver. Implemented as follows:
     *       No Preconditioner:                         precond_=0
     *       Diagonal Preconditioner (default):         precond_=1
     * 13. MPI_Comm local_comm:
     *       The MPI communicator from the local application. Used for ghost cell construction.
     *       The default value of local_comm is MPI_COMM_NULL, i.e. no MPI communicator.
     */

    sampler_rbf(REAL r, const std::vector<point_type> &pts, INT basisFunc = 0,
            bool conservative = false, bool smoothFunc = false, bool generateMatrix = true,
            const std::string &writeFileAddress = std::string(), REAL cutOff = 1e-9,
            REAL cgSolveTol = 1e-6, INT cgMaxIter = 0, INT pouSize = 50,
            INT precond = 1, MPI_Comm local_comm = MPI_COMM_NULL) :
            r_(r),
            pts_(pts),
            basisFunc_(basisFunc),
            conservative_(conservative),
            consistent_(!conservative),
            smoothFunc_(smoothFunc),
            generateMatrix_(generateMatrix),
            writeFileAddress_(writeFileAddress),
            precond_(precond),
            initialised_(false),
            local_mpi_comm_world_(local_comm),
            CABrow_(0),
            CABcol_(0),
            Hrow_(0),
            Hcol_(0),
            pouEnabled_(pouSize == 0 ? false : true),
            cgSolveTol_(cgSolveTol),
            cgMaxIter_(cgMaxIter),
            N_sp_(pouSize),
            M_ap_(pouSize),
            local_rank_(0),
            local_size_(0) {
        //set s to give rbf(r)=cutOff (default 1e-9)
        s_ = std::pow(-std::log(cutOff), 0.5) / r_;
        twor_ = r_ * r_;
        // Ensure CG solver parameters are sensible
        if (cgMaxIter_ < 0)
            cgMaxIter_ = 0;
        if (cgSolveTol_ < 0)
            cgSolveTol_ = 0;

        // Initialise the extened local points (ptsExtend_) by assign local points (pts_)
        ptsExtend_.assign(pts_.begin(), pts_.end());

        if (local_mpi_comm_world_ != MPI_COMM_NULL) {
            MPI_Comm_size(local_mpi_comm_world_, &local_size_);
            MPI_Comm_rank(local_mpi_comm_world_, &local_rank_);
        }

        // Create the matrix output folder if a directory name is given
        if ( !writeFileAddress_.empty() ) {
            int createError = 0;
            #if defined(_WIN32) //Handle creation in Windows
              createError = _mkdir(writeFileAddress_.c_str());
            #else //Handle creation otherwise
              createError = mkdir(writeFileAddress_.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            #endif
            if (createError != 0 && createError != EEXIST)
              EXCEPTION(std::runtime_error("MUI Error [sampler_rbf.h]: Problem creating RBF matrix folder"));
        }
    }

    template<template<typename, typename > class CONTAINER>
    inline OTYPE filter(point_type focus, const CONTAINER<ITYPE, CONFIG> &data_points) const {
      OTYPE sum = 0;

      // RBF matrix not yet created
      if (!initialised_) {
          if (generateMatrix_) { // Generating the matrix
              const clock_t begin_time = std::clock();
              facilitateGhostPoints();
              REAL error = computeRBFtransformationMatrix(data_points, writeFileAddress_);

              if (!QUIET) {
                   std::cout << "MUI [sampler_rbf.h]: Matrices generated in: "
                             << static_cast<double>(std::clock() - begin_time) / CLOCKS_PER_SEC << "s ";
                   if (generateMatrix_) {
                       std::cout << std::endl
                                 << "                     Average CG error: " << error << std::endl;
                   }
               }
          }
          else // Matrix not found and not generating
              EXCEPTION(std::runtime_error("MUI Error [sampler_rbf.h]: RBF matrix not found, call readRBFMatrix() first"));
      }

      //Output for debugging
      if ((!QUIET) && (DEBUG)) {
          std::cout << "MUI [sampler_rbf.h]: Number of remote points: " << data_points.size()
                  << " at rank " << local_rank_ << " out of total ranks "
                  << local_size_ << std::endl;
          std::cout << "MUI [sampler_rbf.h]: Number of local points: " << pts_.size()
                  << " at rank " << local_rank_ << " out of total ranks "
                  << local_size_ << std::endl;
          std::cout << "MUI [sampler_rbf.h]: Number of ghost local points: " << ptsGhost_.size()
                  << " at rank " << local_rank_ << " out of total ranks "
                  << local_size_ << std::endl;
          std::cout << "MUI [sampler_rbf.h]: Number of extended local points: "
                  << ptsExtend_.size() << " at rank " << local_rank_
                  << " out of total ranks " << local_size_ << std::endl;

          for (auto xPtsExtend : ptsExtend_) {
              std::cout << "          ";
              int xPtsExtendSize;
              try {
                  if (sizeof(xPtsExtend) == 0) {
                      throw "MUI Error [sampler_rbf.h]: Error zero xPtsExtend element exception";
                  }
                  else if (sizeof(xPtsExtend) < 0) {
                      throw "MUI Error [sampler_rbf.h]: Error invalid xPtsExtend element exception";
                  }
                  else if (sizeof(xPtsExtend[0]) == 0) {
                      throw "MUI Error [sampler_rbf.h]: Division by zero value of xPtsExtend[0] exception";
                  }
                  else if (sizeof(xPtsExtend[0]) < 0) {
                      throw "MUI Error [sampler_rbf.h]: Division by invalid value of xPtsExtend[0] exception";
                  }
                  xPtsExtendSize = sizeof(xPtsExtend) / sizeof(xPtsExtend[0]);
              }
              catch (const char *msg) {
                  std::cerr << msg << std::endl;
              }
              for (int i = 0; i < xPtsExtendSize; ++i) {
                  std::cout << xPtsExtend[i] << " ";
              }
              std::cout << std::endl;
          }
      }

      auto p = std::find_if(pts_.begin(), pts_.end(), [focus](point_type b) {
                  return normsq(focus - b) < std::numeric_limits<REAL>::epsilon();
              });

      if (p == std::end(pts_))
          EXCEPTION(std::runtime_error("MUI Error [sampler_rbf.h]: Point not found. Must pre-set points for RBF interpolation"));

      auto i = std::distance(pts_.begin(), p);
      REAL tolerance = 1e-5;

      for (size_t j = 0; j < data_points.size(); j++) {
          /// Check whether the order of the remote points coincides with the order when generating the coupling matrix [H]
          ///   and fix it if not. This `check and fix` is essential in the parallel condition as the order of the remote point
          ///   will randomly mixed at the partition boundary and the result will show an randomly oscillating behaviour. The
          ///   below code will ensure a correct match between the remote point and corresponding coupling matrix element.
          OTYPE HElement = 0;

          if (normsq(remote_pts_[j] - data_points[j].first) < (std::numeric_limits<REAL>::epsilon() + tolerance)){
              HElement = H_.get_value(i, j);
          }
          else {

              point_type data_p = data_points[j].first;

              auto remote_p = std::find_if(remote_pts_.begin(), remote_pts_.end(), [data_p, tolerance](point_type b) {
                          return normsq(data_p - b) < (std::numeric_limits<REAL>::epsilon() + tolerance);
                      });

              if (remote_p == std::end(remote_pts_)) {
                  std::cout<< "Missing Remote points: " << data_points[j].first[0] << " " << data_points[j].first[1] << " Size of remote_pts_: " << remote_pts_.size()<< " at rank: " << local_rank_ << std::endl;
                  EXCEPTION(std::runtime_error("MUI Error [sampler_rbf.h]: Remote point not found. Must use the same set of remote points as used to construct the RBF coupling matrix"));
              }

              auto new_j = std::distance(remote_pts_.begin(), remote_p);

              HElement = H_.get_value(i, new_j);
          }

          sum += HElement * data_points[j].second;
      }

      return sum;
    }

    inline geometry::any_shape<CONFIG> support(point_type focus, REAL domain_mag) const {
        return geometry::any_shape<CONFIG>();
    }

    inline void preSetFetchPoints(std::vector<point_type> &pts) {
        pts_ = pts;
        initialised_ = false;
    }

    inline void preSetFetchPointsExtend(std::vector<point_type> &pts) {
        ptsExtend_ = pts;
        initialised_ = false;
    }

    inline void addFetchPoint(point_type pt) {
        pts_.emplace_back(pt);
        initialised_ = false;
    }

    inline void addFetchPointExtend(point_type pt) {
        ptsExtend_.emplace_back(pt);
        initialised_ = false;
    }

    inline void addFetchPointGhost(point_type pt) {
        ptsGhost_.emplace_back(pt);
        initialised_ = false;
    }

    inline void readRBFMatrix(const std::string& readFileAddress) const {
        std::ifstream inputFileMatrixSize(readFileAddress + "/matrixSize.dat");

        if (!inputFileMatrixSize) {
            std::cerr << "MUI Error [sampler_rbf.h]: Could not locate the file address of matrixSize.dat"
                    << std::endl;
        }
        else {
            std::string tempS;
            std::vector<INT> tempV;
            while (std::getline(inputFileMatrixSize, tempS)) {
                // Skips the line if the first two characters are '//'
                if (tempS[0] == '/' && tempS[1] == '/')
                    continue;
                std::stringstream lineStream(tempS);
                std::string tempSS;
                while (std::getline(lineStream, tempSS, ',')) {
                    tempV.emplace_back(std::stoi(tempSS));
                }
            }
            CABrow_ = tempV[0];
            CABcol_ = tempV[1];
            CAArow_ = tempV[2];
            CAAcol_ = tempV[3];
            Hrow_ = tempV[4];
            Hcol_ = tempV[5];
            remote_pts_num_ = tempV[6];
            remote_pts_dim_ = tempV[7];
        }

        std::ifstream inputFileCAB(readFileAddress + "/connectivityAB.dat");

        if (!inputFileCAB) {
            std::cerr
                    << "MUI Error [sampler_rbf.h]: Could not locate the file address on the connectivityAB.dat"
                    << std::endl;
        }
        else {
            connectivityAB_.resize(CABrow_);
            for (INT i = 0; i < CABrow_; i++) {
                connectivityAB_[i].resize(CABcol_, -1);
                std::string tempS;
                while (std::getline(inputFileCAB, tempS)) {
                    // Skips the line if the first two characters are '//'
                    if (tempS[0] == '/' && tempS[1] == '/')
                        continue;
                    std::stringstream lineStream(tempS);
                    std::string tempSS;
                    std::vector<INT> tempV;
                    while (std::getline(lineStream, tempSS, ',')) {
                        tempV.emplace_back(std::stoi(tempSS));
                    }
                    connectivityAB_.emplace_back(tempV);
                }
            }
        }

        if (smoothFunc_) {
            std::ifstream inputFileCAA(readFileAddress + "/connectivityAA.dat");

            if (!inputFileCAA) {
                std::cerr
                        << "MUI Error [sampler_rbf.h]: Could not locate the file address on the connectivityAA.dat"
                        << std::endl;
            }
            else {
                if ((CAArow_ == 0) || (CAAcol_ == 0)) {
                    std::cerr
                            << "MUI Error [sampler_rbf.h]: Error on the size of connectivityAA matrix in matrixSize.dat. Number of rows: "
                            << CAArow_ << " number of columns: " << CAAcol_
                            << ". Make sure matrices were generated with the smoothing function switched on."
                            << std::endl;
                }
                else {
                    connectivityAA_.resize(CAArow_);

                    for (INT i = 0; i < CAArow_; i++) {
                        connectivityAA_[i].resize(CAAcol_, -1);
                        std::string tempS;
                        while (std::getline(inputFileCAA, tempS)) {
                            // Skips the line if the first two characters are '//'
                            if (tempS[0] == '/' && tempS[1] == '/')
                                continue;
                            std::stringstream lineStream(tempS);
                            std::string tempSS;
                            std::vector<INT> tempV;
                            while (std::getline(lineStream, tempSS, ',')) {
                                tempV.emplace_back(std::stoi(tempSS));
                            }
                            connectivityAA_.emplace_back(tempV);
                        }
                    }
                }
            }
        }

        H_.resize(Hrow_, Hcol_);
        H_.set_zero();

        std::ifstream inputFileHMatrix(readFileAddress + "/Hmatrix.dat");

        if (!inputFileHMatrix) {
            std::cerr << "MUI Error [sampler_rbf.h]: Could not locate the file address on the Hmatrix.dat"
                    << std::endl;
        }
        else {
            inputFileHMatrix >> H_;

            if ((H_.get_rows() != Hrow_) || (H_.get_cols() != Hcol_)) {
                std::cerr << "row of H_ (" << H_.get_rows()
                        << ") is not NOT equal to Hrow_ (" << Hrow_ << "), or"
                        << std::endl << "column of H_ (" << H_.get_cols()
                        << ") is not NOT equal to Hcol_ (" << Hcol_ << ")"
                        << std::endl;
            }
        }

        std::ifstream inputFileRemotePoints(readFileAddress + "/remotePoints.dat");

        if (!inputFileRemotePoints) {
            std::cerr << "MUI Error [sampler_rbf.h]: Could not locate the file address on the Hmatrix.dat"
                    << std::endl;
        }
        else {

            if (CONFIG::D != remote_pts_dim_){
                EXCEPTION(std::runtime_error("MUI Error [sampler_rbf.h]: CONFIG::D must equal to remote point dimension in remotePoints.dat"));
            }

            std::string tempS;
            while (std::getline(inputFileRemotePoints, tempS)) {
                // Skips the line if the first two characters are '//'
                if (tempS[0] == '/' && tempS[1] == '/')
                    continue;
                std::stringstream lineStream(tempS);
                std::string tempSS;
                point_type tempP;
                size_t temp_i = 0;
                while (std::getline(lineStream, tempSS, ',')) {
                    assert(temp_i<static_cast<size_t>(remote_pts_dim_));
                    tempP[temp_i] = static_cast<REAL>(std::stold(tempSS));
                    temp_i++;
                }
                remote_pts_.emplace_back(tempP);
            }

            assert(remote_pts_.size() == static_cast<size_t>(remote_pts_num_));
        }

        initialised_ = true;
    }

    // Functions to facilitate ghost points

    // Determine bounding box of local points
    std::pair<point_type, point_type> localBoundingBox(const std::vector<point_type> pt) const {
        point_type lbbMin, lbbMax;
        REAL maxVal = std::numeric_limits<REAL>::max();
        try {
            if (CONFIG::D == 1) {
                lbbMax[0] = -maxVal;
                lbbMin[0] = maxVal;
            }
            else if (CONFIG::D == 2) {
                lbbMax[0] = -maxVal;
                lbbMax[1] = -maxVal;
                lbbMin[0] = maxVal;
                lbbMin[1] = maxVal;
            }
            else if (CONFIG::D == 3) {
                lbbMax[0] = -maxVal;
                lbbMax[1] = -maxVal;
                lbbMax[2] = -maxVal;
                lbbMin[0] = maxVal;
                lbbMin[1] = maxVal;
                lbbMin[2] = maxVal;
            }
            else {
                throw "MUI Error [sampler_rbf.h]: Invalid value of CONFIG::D exception";
            }
        }
        catch (const char *msg) {
            std::cerr << msg << std::endl;
        }

        for (auto xPts : pt) {
            try {
                switch (CONFIG::D) {
                case 1:
                    if (xPts[0] > lbbMax[0])
                        lbbMax[0] = xPts[0];
                    if (xPts[0] < lbbMin[0])
                        lbbMin[0] = xPts[0];
                    break;
                case 2:
                    if (xPts[0] > lbbMax[0])
                        lbbMax[0] = xPts[0];
                    if (xPts[0] < lbbMin[0])
                        lbbMin[0] = xPts[0];
                    if (xPts[1] > lbbMax[1])
                        lbbMax[1] = xPts[1];
                    if (xPts[1] < lbbMin[1])
                        lbbMin[1] = xPts[1];
                    break;
                case 3:
                    if (xPts[0] > lbbMax[0])
                        lbbMax[0] = xPts[0];
                    if (xPts[0] < lbbMin[0])
                        lbbMin[0] = xPts[0];
                    if (xPts[1] > lbbMax[1])
                        lbbMax[1] = xPts[1];
                    if (xPts[1] < lbbMin[1])
                        lbbMin[1] = xPts[1];
                    if (xPts[2] > lbbMax[2])
                        lbbMax[2] = xPts[2];
                    if (xPts[2] < lbbMin[2])
                        lbbMin[2] = xPts[2];
                    break;
                default:
                    throw "MUI Error [sampler_rbf.h]: Invalid value of CONFIG::D exception";
                }
            }
            catch (const char *msg) {
                std::cerr << msg << std::endl;
            }
        }
        return std::make_pair(lbbMin, lbbMax);
    }

    // Determine extended bounding box of local points include ghost area
    std::pair<point_type, point_type> localExtendBoundingBox(std::pair<point_type, point_type> lbb, REAL r) const {
        point_type lbbExtendMin, lbbExtendMax;
        lbbExtendMin = lbb.first;
        lbbExtendMax = lbb.second;
        for (INT i = 0; i < CONFIG::D; ++i) {
            lbbExtendMax[i] += r;
            lbbExtendMin[i] -= r;
        }
        return std::make_pair(lbbExtendMin, lbbExtendMax);
    }

    // Collect points that will be sent to other processors as ghost points
    std::pair<std::vector<std::pair<INT, INT>>, std::vector<std::pair<INT, std::vector<point_type>>>> getGhostPointsToSend(
            std::pair<point_type, point_type> lbbExtend, MPI_Comm local_world,
            int local_rank, int local_size) const {
        std::pair<INT, std::pair<point_type, point_type>> localBB, globalBB[local_size];

        // Assign rank ID and extended bounding box of local points to localBB
        localBB.first = local_rank;
        localBB.second.first = lbbExtend.second;
        localBB.second.second = lbbExtend.first;

        // Gather the bounding boxes from all processes
        MPI_Allgather(&localBB,
                sizeof(std::pair<INT, std::pair<point_type, point_type>>),
                MPI_BYTE, globalBB,
                sizeof(std::pair<INT, std::pair<point_type, point_type>>),
                MPI_BYTE, local_world);

        // output for debugging
        if ((!QUIET) && (DEBUG)) {
            for (auto xGlobalBB : globalBB) {
                std::cout << "MUI [sampler_rbf.h]: globalBBs: " << xGlobalBB.first << " "
                        << xGlobalBB.second.first[0] << " "
                        << xGlobalBB.second.first[1] << " "
                        << xGlobalBB.second.second[0] << " "
                        << xGlobalBB.second.second[1] << " at rank "
                        << local_rank_ << std::endl;
            }
        }

        // Declare vector to gather number of points to send to each processors
        std::vector<std::pair<INT, INT>> ghostPointsCountToSend;
        for (auto xGlobalBB : globalBB) {
            if (xGlobalBB.first == local_rank)
                continue;
            ghostPointsCountToSend.push_back(std::make_pair(xGlobalBB.first, 0));
        }

        // Loop over local points to collect ghost points for other processors
        std::vector<std::pair<INT, std::vector<point_type>>> ghostPointsToSend;
        for (size_t i = 0; i < pts_.size(); ++i) {
            for (auto xGlobalBB : globalBB) {
                if (xGlobalBB.first == local_rank)
                    continue;
                try {
                    switch (CONFIG::D) {
                    case 1: {
                        if ((pts_[i][0] <= xGlobalBB.second.first[0])
                                && (pts_[i][0] >= xGlobalBB.second.second[0])) {

                            auto ghostPointsToSendIter =
                                    std::find_if(ghostPointsToSend.begin(),
                                            ghostPointsToSend.end(),
                                            [xGlobalBB](std::pair<INT, std::vector<point_type>> b) {
                                                return b.first == xGlobalBB.first;
                                            });

                            auto ghostPointsCountToSendIter = std::find_if(
                                    ghostPointsCountToSend.begin(),
                                    ghostPointsCountToSend.end(),
                                    [xGlobalBB](std::pair<INT, INT> b) {
                                        return b.first == xGlobalBB.first;
                                    });

                            if (ghostPointsToSendIter
                                    == std::end(ghostPointsToSend)) {
                                std::vector<point_type> vecPointTemp { pts_[i] };
                                ghostPointsToSend.push_back(std::make_pair(xGlobalBB.first,    vecPointTemp));
                            }
                            else
                                ghostPointsToSendIter->second.push_back(pts_[i]);

                            assert(ghostPointsCountToSendIter != std::end(ghostPointsCountToSend));

                            ++ghostPointsCountToSendIter->second;
                        }
                        break;
                    }
                    case 2: {
                        if ((pts_[i][0] <= xGlobalBB.second.first[0])
                                && (pts_[i][0] >= xGlobalBB.second.second[0])
                                && (pts_[i][1] <= xGlobalBB.second.first[1])
                                && (pts_[i][1] >= xGlobalBB.second.second[1])) {

                            auto ghostPointsToSendIter =
                                    std::find_if(ghostPointsToSend.begin(),
                                            ghostPointsToSend.end(),
                                            [xGlobalBB](std::pair<INT, std::vector<point_type>> b) {
                                                return b.first == xGlobalBB.first;
                                            });

                            auto ghostPointsCountToSendIter = std::find_if(
                                    ghostPointsCountToSend.begin(),
                                    ghostPointsCountToSend.end(),
                                    [xGlobalBB](std::pair<INT, INT> b) {
                                        return b.first == xGlobalBB.first;
                                    });

                            if (ghostPointsToSendIter
                                    == std::end(ghostPointsToSend)) {
                                std::vector<point_type> vecPointTemp { pts_[i] };
                                ghostPointsToSend.push_back(std::make_pair(xGlobalBB.first,    vecPointTemp));
                            }
                            else
                                ghostPointsToSendIter->second.push_back(pts_[i]);

                            assert(ghostPointsCountToSendIter != std::end(ghostPointsCountToSend));

                            ++ghostPointsCountToSendIter->second;
                        }
                        break;
                    }
                    case 3: {
                        if ((pts_[i][0] <= xGlobalBB.second.first[0])
                                && (pts_[i][0] >= xGlobalBB.second.second[0])
                                && (pts_[i][1] <= xGlobalBB.second.first[1])
                                && (pts_[i][1] >= xGlobalBB.second.second[1])
                                && (pts_[i][2] <= xGlobalBB.second.first[2])
                                && (pts_[i][2] >= xGlobalBB.second.second[2])) {

                            auto ghostPointsToSendIter =
                                    std::find_if(ghostPointsToSend.begin(),
                                            ghostPointsToSend.end(),
                                            [xGlobalBB](std::pair<INT, std::vector<point_type>> b) {
                                                return b.first == xGlobalBB.first;
                                            });

                            auto ghostPointsCountToSendIter = std::find_if(
                                    ghostPointsCountToSend.begin(),
                                    ghostPointsCountToSend.end(),
                                    [xGlobalBB](std::pair<INT, INT> b) {
                                        return b.first == xGlobalBB.first;
                                    });

                            if (ghostPointsToSendIter
                                    == std::end(ghostPointsToSend)) {
                                std::vector<point_type> vecPointTemp { pts_[i] };
                                ghostPointsToSend.push_back(std::make_pair(xGlobalBB.first, vecPointTemp));
                            }
                            else
                                ghostPointsToSendIter->second.push_back(pts_[i]);

                            assert(ghostPointsCountToSendIter != std::end(ghostPointsCountToSend));

                            ++ghostPointsCountToSendIter->second;
                        }
                        break;
                    }
                    default:
                        throw "MUI Error [sampler_rbf.h]: Invalid value of CONFIG::D exception";
                    }
                }
                catch (const char *msg) {
                    std::cerr << msg << std::endl;
                }
            }
        }

        // output for debugging
        if ((!QUIET) && (DEBUG)) {
            std::cout << "MUI [sampler_rbf.h]: Total size of GhostPointsToSend "
                    << ghostPointsToSend.size() << " at rank " << local_rank
                    << std::endl;
            for (auto xGhostPointsToSend : ghostPointsToSend) {
                std::cout << "MUI [sampler_rbf.h]: xGhostPointsToSend to rank "
                        << xGhostPointsToSend.first << " at rank " << local_rank
                        << std::endl;
                for (auto xVectPts : xGhostPointsToSend.second) {
                    std::cout << "MUI [sampler_rbf.h]: " << xVectPts[0] << " " << xVectPts[1]
                            << " at rank " << local_rank << std::endl;
                }
            }
            for (auto xGhostPointsCountToSend : ghostPointsCountToSend) {
                std::cout << "MUI [sampler_rbf.h]: xGhostPointsCountToSend to rank "
                        << xGhostPointsCountToSend.first << " has "
                        << xGhostPointsCountToSend.second
                        << " points to send at rank " << local_rank_
                        << std::endl;
            }
        }
        return std::make_pair(ghostPointsCountToSend, ghostPointsToSend);
    }

    // Distribution of ghost points among processors by all to all
    std::vector<point_type> distributeGhostPoints(
            std::vector<std::pair<INT, INT>> ghostPointsCountToSend,
            std::vector<std::pair<INT, std::vector<point_type>>> ghostPointsToSend,
            MPI_Comm local_world, int local_rank) const {
        std::vector<point_type> ptsGhost;
        std::vector<std::pair<INT, INT>> ghostPointsCountToRecv;
        for (auto xGhostPointsCountToSend : ghostPointsCountToSend) {
            assert(xGhostPointsCountToSend.first != local_rank);
            assert(xGhostPointsCountToSend.second >= 0);

            // Determined of number of points to transfer by pairwise communication
            MPI_Send(&xGhostPointsCountToSend.second, 1, MPI_INT,
                    xGhostPointsCountToSend.first, 0, local_world);
            int pointsCountTemp = -1;
            MPI_Recv(&pointsCountTemp, 1, MPI_INT,
                    xGhostPointsCountToSend.first, 0, local_world,
                    MPI_STATUS_IGNORE);
            assert(pointsCountTemp >= 0);
            ghostPointsCountToRecv.push_back(
                    std::make_pair(xGhostPointsCountToSend.first,
                            pointsCountTemp));

            // Send ghost points by pairwise communication
            if (xGhostPointsCountToSend.second != 0) {
                auto ghostPointsToSendIter = std::find_if(
                        ghostPointsToSend.begin(), ghostPointsToSend.end(),
                        [xGhostPointsCountToSend](
                                std::pair<INT, std::vector<point_type>> b) {
                            return b.first == xGhostPointsCountToSend.first;
                        });

                assert(ghostPointsToSendIter->second.size() == static_cast<size_t>(xGhostPointsCountToSend.second));

                int buffer_size = ghostPointsToSendIter->second.size() * sizeof(point_type);
                char *buffer = new char[buffer_size];
                int position = 0;
                for (auto &pointElement : ghostPointsToSendIter->second) {
                    MPI_Pack(&pointElement, sizeof(point_type), MPI_BYTE,
                            buffer, buffer_size, &position, local_world);
                }
                MPI_Send(buffer, position, MPI_PACKED, xGhostPointsCountToSend.first, 1, local_world);
                delete[] buffer;

                // output for debugging
                if ((!QUIET) && (DEBUG)) {
                    for (auto xsend : ghostPointsToSendIter->second) {
                        std::cout << "MUI [sampler_rbf.h]: Send ghost point to rank "
                                << xGhostPointsCountToSend.first
                                << " with value of " << xsend[0] << " "
                                << xsend[1] << " at rank " << local_rank
                                << std::endl;
                    }
                }
            }

            // Receive ghost points by pairwise communication
            auto ghostPointsCountToRecvIter = std::find_if(
                    ghostPointsCountToRecv.begin(),
                    ghostPointsCountToRecv.end(),
                    [xGhostPointsCountToSend](std::pair<INT, INT> b) {
                        return b.first == xGhostPointsCountToSend.first;
                    });

            if (ghostPointsCountToRecvIter->second != 0) {
                std::vector<point_type> vecPointTemp;
                MPI_Status status;
                MPI_Probe(ghostPointsCountToRecvIter->first, 1, local_world,
                        &status);
                int buffer_size;
                MPI_Get_count(&status, MPI_PACKED, &buffer_size);
                char *buffer = new char[buffer_size];
                MPI_Recv(buffer, buffer_size, MPI_PACKED,
                        ghostPointsCountToRecvIter->first, 1, local_world,
                        MPI_STATUS_IGNORE);
                int position = 0;
                int count;
                MPI_Get_elements(&status, MPI_BYTE, &count);
                vecPointTemp.resize(count / sizeof(point_type));
                for (auto &pointElement : vecPointTemp) {
                    MPI_Unpack(buffer, buffer_size, &position, &pointElement,
                            sizeof(point_type), MPI_BYTE, local_world);
                }
                delete[] buffer;
                for (auto xvecPointTemp : vecPointTemp) {
                    ptsGhost.emplace_back(xvecPointTemp);
                    // output for debugging
                    if ((!QUIET) && (DEBUG)) {
                        std::cout << "MUI [sampler_rbf.h]: Receive ghost point from rank "
                                << ghostPointsCountToRecvIter->first
                                << " with value of " << xvecPointTemp[0] << " "
                                << xvecPointTemp[1] << " at rank " << local_rank
                                << std::endl;
                    }
                }
            }
        }
        return ptsGhost;
    }

    // Facilitate Ghost points
    void facilitateGhostPoints() const {
        std::pair<point_type, point_type> lbb = localBoundingBox(pts_);
        std::pair<point_type, point_type> lbbExtend = localExtendBoundingBox(lbb, r_);
        if (local_mpi_comm_world_ != MPI_COMM_NULL) {
            std::pair<std::vector<std::pair<INT, INT>>, std::vector<std::pair<INT, std::vector<point_type>>>>
            ghostPointsToSendPair = getGhostPointsToSend(lbbExtend, local_mpi_comm_world_, local_rank_, local_size_);
            ptsGhost_ = distributeGhostPoints(ghostPointsToSendPair.first,
                                              ghostPointsToSendPair.second, local_mpi_comm_world_, local_rank_);
            // output for debugging
            if ((!QUIET) && (DEBUG)) {
                std::cout << "MUI [sampler_rbf.h]: Local bounding box: " << lbb.second[0] << " "
                        << lbb.second[1] << " " << lbb.first[0] << " "
                        << lbb.first[1] << " at rank " << local_rank_
                        << " out of total ranks " << local_size_ << std::endl;
                std::cout << "MUI [sampler_rbf.h]: Extended local bounding box: "
                        << lbbExtend.second[0] << " " << lbbExtend.second[1]
                        << " " << lbbExtend.first[0] << " "
                        << lbbExtend.first[1] << " at rank " << local_rank_
                        << " out of total ranks " << local_size_ << std::endl;
            }
        }
        // Construct the extended local points by combine local points with ghost points
        ptsExtend_.insert(ptsExtend_.end(), ptsGhost_.begin(), ptsGhost_.end());
    }

private:
    template<template<typename, typename > class CONTAINER>
    REAL computeRBFtransformationMatrix(const CONTAINER<ITYPE, CONFIG> &data_points, const std::string &fileAddress) const {
      bool writeMatrix = !fileAddress.empty();

      // Refine partition size depending on if PoU enabled, whether conservative or consistent
        // and if problem size smaller than defined patch size
        if (conservative_) { // Conservative RBF, using local point set for size
            if (pouEnabled_) { // PoU enabled so check patch size not larger than point set
                if (ptsExtend_.size() < N_sp_)
                    N_sp_ = ptsExtend_.size();
            }
            else
                N_sp_ = ptsExtend_.size();
        }

        if (consistent_) { // Consistent RBF, using remote point set for size
            if (pouEnabled_) { // PoU enabled so check patch size not larger than point set
                if (data_points.size() < N_sp_)
                    N_sp_ = data_points.size();
            }
            else
                N_sp_ = data_points.size();
        }

        if (smoothFunc_) {
            if (ptsExtend_.size() < M_ap_)
                M_ap_ = ptsExtend_.size() - 1;
        }

        REAL errorReturn = 0;

        if (conservative_)
            buildConnectivityConservative(data_points, N_sp_, writeMatrix, fileAddress);
        else
            buildConnectivityConsistent(data_points, N_sp_, writeMatrix, fileAddress);

        H_.resize(ptsExtend_.size(), data_points.size());
        H_.set_zero();

        if (smoothFunc_) {
            buildConnectivitySmooth(M_ap_, writeMatrix, fileAddress);
            H_toSmooth_.resize(ptsExtend_.size(), data_points.size());
            H_toSmooth_.set_zero();
        }

        if (writeMatrix) {
            std::ofstream outputFileMatrixSize(fileAddress + "/matrixSize.dat");

            if (!outputFileMatrixSize) {
                std::cerr
                        << "MUI Error [sampler_rbf.h]: Could not locate the file address of matrixSize.dat!"
                        << std::endl;
            }
            else {
                outputFileMatrixSize
                        << "// *********************************************************************************************************************************************";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize
                        << "// **** This is the 'matrixSize.dat' file of the RBF spatial sampler of the MUI library";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize
                        << "// **** This file contains the size (number of rows and number of columns) of the Point Connectivity Matrix (N) and the Coupling Matrix (H).";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize
                        << "// **** The file uses the Comma-Separated Values (CSV) format and the ASCII format with the meanings as follows: ";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize
                        << "// ****            The number of rows of the Point Connectivity Matrix (N), ";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize
                        << "// ****            The number of columns of the Point Connectivity Matrix (N),";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize
                        << "// ****            The number of rows of the Point Connectivity Matrix (M) (for smoothing), ";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize
                        << "// ****            The number of columns of the Point Connectivity Matrix (M) (for smoothing),";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize
                        << "// ****            The number of rows of the Coupling Matrix (H),";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize
                        << "// ****            The number of columns of the Coupling Matrix (H)";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize
                        << "// ****            The size of remote points";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize
                        << "// ****            The dimension of remote points";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize
                        << "// *********************************************************************************************************************************************";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize << "//  ";
                outputFileMatrixSize << "\n";
                outputFileMatrixSize << connectivityAB_.size();
                outputFileMatrixSize << ",";
                outputFileMatrixSize << connectivityAB_[0].size();
                outputFileMatrixSize << ",";
                if (smoothFunc_) {
                    outputFileMatrixSize << connectivityAA_.size();
                    outputFileMatrixSize << ",";
                    outputFileMatrixSize << connectivityAA_[0].size();
                    outputFileMatrixSize << ",";
                }
                else {
                    outputFileMatrixSize << "0";
                    outputFileMatrixSize << ",";
                    outputFileMatrixSize << "0";
                    outputFileMatrixSize << ",";
                }
                outputFileMatrixSize << H_.get_rows();
                outputFileMatrixSize << ",";
                outputFileMatrixSize << H_.get_cols();
                outputFileMatrixSize << ",";
                outputFileMatrixSize << data_points.size();
                outputFileMatrixSize << ",";
                outputFileMatrixSize << CONFIG::D;
                outputFileMatrixSize << "\n";
            }
        }

        if (conservative_) { // Build matrix for conservative RBF
            errorReturn = buildMatrixConservative(data_points, N_sp_, M_ap_, smoothFunc_, pouEnabled_);
        }
        else {
            // Build matrix for consistent RBF
            errorReturn = buildMatrixConsistent(data_points, N_sp_, M_ap_, smoothFunc_, pouEnabled_);
        }

        if (writeMatrix) {
            std::ofstream outputFileHMatrix(fileAddress + "/Hmatrix.dat");

            if (!outputFileHMatrix) {
                std::cerr
                        << "MUI Error [sampler_rbf.h]: Could not locate the file address of Hmatrix.dat!"
                        << std::endl;
            }
            else {
                outputFileHMatrix
                        << "// ************************************************************************************************";
                outputFileHMatrix << "\n";
                outputFileHMatrix
                        << "// **** This is the 'Hmatrix.dat' file of the RBF spatial sampler of the MUI library";
                outputFileHMatrix << "\n";
                outputFileHMatrix
                        << "// **** This file contains the entire matrix of the Coupling Matrix (H).";
                outputFileHMatrix << "\n";
                outputFileHMatrix
                        << "// **** The file uses the Comma-Separated Values (CSV) format with ASCII for the entire H matrix";
                outputFileHMatrix << "\n";
                outputFileHMatrix
                        << "// ************************************************************************************************";
                outputFileHMatrix << "\n";
                outputFileHMatrix << "// ";
                outputFileHMatrix << "\n";
                outputFileHMatrix << H_;
            }
        }

        initialised_ = true;

        return errorReturn;
    }

    template<template<typename, typename > class CONTAINER>
    inline REAL buildMatrixConsistent(const CONTAINER<ITYPE, CONFIG> &data_points, const size_t NP,
                const size_t MP, bool smoothing, bool pou) const {
        REAL errorReturn = 0;
        std::pair<INT, REAL> iterErrorReturn(0, 0);
        if( pou ) { // Using PoU approach
            for (size_t row = 0; row < ptsExtend_.size(); row++) {
                linalg::sparse_matrix<INT, REAL> Css; //< Matrix of radial basis function evaluations between prescribed points
                linalg::sparse_matrix<INT, REAL> Aas; //< Matrix of RBF evaluations between prescribed and interpolation points

                Css.resize((1 + NP + CONFIG::D), (1 + NP + CONFIG::D));
                Aas.resize((1 + NP + CONFIG::D), 1);

                // Set matrix Css
                // Define intermediate matrix for performance purpose
                linalg::sparse_matrix<INT, REAL> Css_coo((1 + NP + CONFIG::D),(1 + NP + CONFIG::D),"COO");

                for (size_t i = 0; i < NP; i++) {
                    for (size_t j = i; j < NP; j++) {
                        int glob_i = connectivityAB_[row][i];
                        int glob_j = connectivityAB_[row][j];

                        auto d = norm(
                                data_points[glob_i].first
                                        - data_points[glob_j].first);

                        if (d < r_) {
                            REAL w = rbf(d);
                            Css_coo.set_value(i, j, w, false);
//                            Css.set_value(i, j, w);
                            if (i != j) {
                                Css_coo.set_value(j, i, w, false);
//                                Css.set_value(j, i, w);
                            }
                        }
                    }
                }

                for (size_t i = 0; i < NP; i++) {
                    Css_coo.set_value(i, NP, 1, false);
                    Css_coo.set_value(NP, i, 1, false);
//                    Css.set_value(i, NP, 1);
//                    Css.set_value(NP, i, 1);

                    int glob_i = connectivityAB_[row][i];

                    for (INT dim = 0; dim < CONFIG::D; dim++) {
                        Css_coo.set_value(i, (NP + dim + 1),
                                data_points[glob_i].first[dim], false);
//                        Css.set_value(i, (NP + dim + 1),
//                                data_points[glob_i].first[dim]);
                        Css_coo.set_value((NP + dim + 1), i,
                                data_points[glob_i].first[dim], false);
//                        Css.set_value((NP + dim + 1), i,
//                                data_points[glob_i].first[dim]);
                    }
                }
                Css = Css_coo;

                // Set Aas
                // Define intermediate matrix for performance purpose
                linalg::sparse_matrix<INT, REAL> Aas_coo((1 + NP + CONFIG::D),1,"COO");

                for (size_t j = 0; j < NP; j++) {
                    int glob_j = connectivityAB_[row][j];

                    auto d = norm(ptsExtend_[row] - data_points[glob_j].first);

                    if (d < r_) {
                        Aas_coo.set_value(j, 0, rbf(d), false);
//                        Aas.set_value(j, 0, rbf(d));
                    }
                }
                Aas_coo.set_value(NP, 0, 1, false);
//                Aas.set_value(NP, 0, 1);

                for (int dim = 0; dim < CONFIG::D; dim++) {
                    Aas_coo.set_value((NP + dim + 1), 0, ptsExtend_[row][dim], false);
//                    Aas.set_value((NP + dim + 1), 0, ptsExtend_[row][dim]);
                }
                Aas = Aas_coo;

                linalg::sparse_matrix<INT, REAL> H_i;

                if (precond_ == 0) {
                    linalg::conjugate_gradient<INT, REAL> cg(Css, Aas, cgSolveTol_, cgMaxIter_);
                    iterErrorReturn = cg.solve();
                    H_i = cg.getSolution();
                }
                else if (precond_ == 1) {
                    linalg::diagonal_preconditioner<INT, REAL> M(Css);
                    linalg::conjugate_gradient<INT, REAL> cg(Css, Aas, cgSolveTol_, cgMaxIter_, &M);
                    iterErrorReturn = cg.solve();
                    H_i = cg.getSolution();
                }
                else {
                    std::cerr
                            << "MUI Error [sampler_rbf.h]: Invalid CG Preconditioner function number ("
                            << precond_ << ")" << std::endl
                            << "Please set the CG Preconditioner function number (precond_) as: "
                            << std::endl << "precond_=0 (No Preconditioner); "
                            << std::endl << "precond_=1 (Diagonal Preconditioner); " << std::endl;
                    std::abort();
                }

                if (DEBUG) {
                    std::cout << "MUI [sampler_rbf.h]: #iterations of H_i:     " << iterErrorReturn.first
                            << ". Error of H_i: " << iterErrorReturn.second
                            << std::endl;
                }

                errorReturn += iterErrorReturn.second;

                if (smoothing) {
                    for (size_t j = 0; j < NP; j++) {
                        INT glob_j = connectivityAB_[row][j];
                        H_toSmooth_.set_value(row, glob_j, H_i.get_value(j, 0));
                    }
                }
                else {
                    for (size_t j = 0; j < NP; j++) {
                        INT glob_j = connectivityAB_[row][j];
                        H_.set_value(row, glob_j, H_i.get_value(j, 0));
                    }
                }
            }

            errorReturn /= static_cast<REAL>(pts_.size());

            if (smoothing) {
                for (size_t row = 0; row < ptsExtend_.size(); row++) {
                    for (size_t j = 0; j < NP; j++) {
                        INT glob_j = connectivityAB_[row][j];
                        REAL h_j_sum = 0.;
                        REAL f_sum = 0.;

                        for (size_t k = 0; k < MP; k++) {
                            INT row_k = connectivityAA_[row][k];
                            if (row_k == static_cast<INT>(row)) {
                                std::cerr << "MUI Error [sampler_rbf.h]: Invalid row_k value: " << row_k
                                        << std::endl;
                            }
                            else
                                h_j_sum += std::pow(dist_h_i(row, row_k), -2.);
                        }

                        for (size_t k = 0; k < MP; k++) {
                            INT row_k = connectivityAA_[row][k];
                            if (row_k == static_cast<INT>(row)) {
                                std::cerr << "MUI Error [sampler_rbf.h]: Invalid row_k value: " << row_k
                                        << std::endl;
                            }
                            else {
                                REAL w_i = ((std::pow(dist_h_i(row, row_k), -2.)) / (h_j_sum));
                                f_sum += w_i * H_toSmooth_.get_value(row_k, glob_j);
                            }
                        }
                        H_.set_value(row, glob_j, (0.5 * (f_sum + H_toSmooth_.get_value(row, glob_j))));
                    }
                }
            }
        }
        else { // Not using PoU
            linalg::sparse_matrix<INT,REAL> Css; //< Matrix of radial basis function evaluations between prescribed points
            linalg::sparse_matrix<INT,REAL> Aas; //< Matrix of RBF evaluations between prescribed and interpolation points

            Css.resize((1 + data_points.size() + CONFIG::D), (1 + data_points.size() + CONFIG::D));
            Aas.resize(ptsExtend_.size(), (1 + data_points.size() + CONFIG::D));

            //set Css
            // Define intermediate matrix for performance purpose
            linalg::sparse_matrix<INT, REAL> Css_coo((1 + data_points.size() + CONFIG::D),(1 + data_points.size() + CONFIG::D),"COO");

            for ( size_t i = 0; i < data_points.size(); i++ ) {
                for ( size_t j = i; j < data_points.size(); j++ ) {
                    auto d = norm(data_points[i].first - data_points[j].first);

                    if ( d < r_ ) {
                        REAL w = rbf(d);
                        Css_coo.set_value((i + CONFIG::D + 1), (j + CONFIG::D + 1), w, false);
//                        Css.set_value((i + CONFIG::D + 1), (j + CONFIG::D + 1), w);

                        if ( i != j )
                            Css_coo.set_value((j + CONFIG::D + 1), (i + CONFIG::D + 1), w, false);
//                            Css.set_value((j + CONFIG::D + 1), (i + CONFIG::D + 1), w);
                    }
                }
            }
            Css = Css_coo;

            //set Aas
            // Define intermediate matrix for performance purpose
            linalg::sparse_matrix<INT, REAL> Aas_coo(ptsExtend_.size(),(1 + data_points.size() + CONFIG::D),"COO");

            for ( size_t i = 0; i < ptsExtend_.size(); i++ ) {
                for ( size_t j = 0; j < data_points.size(); j++ ) {
                    auto d = norm(ptsExtend_[i] - data_points[j].first);

                    if ( d < r_ ) {

                        Aas_coo.set_value(i, (j + CONFIG::D + 1), rbf(d), false);
//                        Aas.set_value(i, (j + CONFIG::D + 1), rbf(d));
                    }
                }
            }
            Aas = Aas_coo;

            linalg::sparse_matrix<INT,REAL> Aas_trans = Aas.transpose();

            linalg::sparse_matrix<INT, REAL> H_more;

            if (precond_ == 0) {
                linalg::conjugate_gradient<INT, REAL> cg(Css, Aas_trans, cgSolveTol_, cgMaxIter_);
                iterErrorReturn = cg.solve();
                H_more = cg.getSolution();
            }
            else if (precond_ == 1) {
                linalg::diagonal_preconditioner<INT, REAL> M(Css);
                linalg::conjugate_gradient<INT, REAL> cg(Css, Aas_trans, cgSolveTol_, cgMaxIter_, &M);
                iterErrorReturn = cg.solve();
                H_more = cg.getSolution();
            }
            else {
                std::cerr
                        << "MUI Error [sampler_rbf.h]: Invalid CG Preconditioner function number ("
                        << precond_ << ")" << std::endl
                        << "Please set the CG Preconditioner function number (precond_) as: "
                        << std::endl << "precond_=0 (No Preconditioner); "
                        << std::endl << "precond_=1 (Diagonal Preconditioner); " << std::endl;
                std::abort();
            }

            if (DEBUG) {
                std::cout << "MUI [sampler_rbf.h]: #iterations of H_more:     " << iterErrorReturn.first
                        << ". Error of H_more: " << iterErrorReturn.second
                        << std::endl;
            }

            errorReturn = iterErrorReturn.second;

            if ( smoothing ) {
                for ( size_t i = 0; i < data_points.size(); i++ ) {
                    for (size_t j = 0; j < ptsExtend_.size(); j++ ) {
                        H_toSmooth_.set_value(j, i, H_more.get_value((i + CONFIG::D + 1), j));
                    }
                }
            }
            else {
                for ( size_t i = 0; i < data_points.size(); i++ ) {
                    for ( size_t j = 0; j < ptsExtend_.size(); j++ ) {
                        H_.set_value(j, i, H_more.get_value((i + CONFIG::D + 1), j));
                    }
                }
            }

            if ( smoothing ) {
                for ( size_t row = 0; row < ptsExtend_.size(); row++ ) {
                    for ( size_t j = 0; j < data_points.size(); j++ ) {
                        REAL h_j_sum = 0.;
                        REAL f_sum = 0.;
                        for ( size_t k = 0; k < MP; k++ ) {
                            INT row_k = connectivityAA_[row][k];
                            if ( row_k == static_cast<INT>(row) )
                                std::cerr << "Invalid row_k value: " << row_k << std::endl;
                            else
                                h_j_sum += std::pow(dist_h_i(row, row_k), -2.);
                        }

                        for ( size_t k = 0; k < MP; k++ ) {
                            INT row_k = connectivityAA_[row][k];
                            if ( row_k == static_cast<INT>(row) )
                                std::cerr << "Invalid row_k value: " << row_k << std::endl;
                            else {
                                REAL w_i = ((std::pow(dist_h_i(row, row_k), -2.)) / (h_j_sum));
                                f_sum += w_i * H_toSmooth_.get_value(row_k, j);
                            }
                        }
                        H_.set_value(row, j, 0.5 * (f_sum + H_toSmooth_.get_value(row, j)));
                    }
                }
            }
        }

        return errorReturn;
    }

    template<template<typename, typename > class CONTAINER>
    inline REAL buildMatrixConservative(
            const CONTAINER<ITYPE, CONFIG> &data_points, const size_t NP,
            const size_t MP, bool smoothing, bool pou) const {
        REAL errorReturn = 0;
        std::pair<INT, REAL> iterErrorReturn(0, 0);
        if( pou ) { // Using partitioned approach
            for (size_t row = 0; row < data_points.size(); row++) {
                linalg::sparse_matrix<INT, REAL> Css; //< Matrix of radial basis function evaluations between prescribed points
                linalg::sparse_matrix<INT, REAL> Aas; //< Matrix of RBF evaluations between prescribed and interpolation points

                Css.resize((1 + NP + CONFIG::D), (1 + NP + CONFIG::D));
                Aas.resize((1 + NP + CONFIG::D), 1);

                //set Css
                // Define intermediate matrix for performance purpose
                linalg::sparse_matrix<INT, REAL> Css_coo((1 + NP + CONFIG::D),(1 + NP + CONFIG::D),"COO");

                for (size_t i = 0; i < NP; i++) {
                    for (size_t j = i; j < NP; j++) {
                        INT glob_i = connectivityAB_[row][i];
                        INT glob_j = connectivityAB_[row][j];

                        auto d = norm(ptsExtend_[glob_i] - ptsExtend_[glob_j]);

                        if (d < r_) {
                            REAL w = rbf(d);
                            Css_coo.set_value(i, j, w, false);
//                            Css.set_value(i, j, w);
                            if (i != j) {
                                Css_coo.set_value(j, i, w, false);
//                                Css.set_value(j, i, w);
                            }
                        }
                    }
                }

                for (size_t i = 0; i < NP; i++) {
                    Css_coo.set_value(i, NP, 1, false);
//                    Css.set_value(i, NP, 1);
                    Css_coo.set_value(NP, i, 1, false);
//                    Css.set_value(NP, i, 1);

                    INT glob_i = connectivityAB_[row][i];

                    for (INT dim = 0; dim < CONFIG::D; dim++) {
                        Css_coo.set_value(i, (NP + dim + 1), ptsExtend_[glob_i][dim], false);
//                        Css.set_value(i, (NP + dim + 1), ptsExtend_[glob_i][dim]);
                        Css_coo.set_value((NP + dim + 1), i, ptsExtend_[glob_i][dim], false);
//                        Css.set_value((NP + dim + 1), i, ptsExtend_[glob_i][dim]);
                    }
                }

                Css = Css_coo;

                //set Aas
                // Define intermediate matrix for performance purpose
                linalg::sparse_matrix<INT, REAL> Aas_coo((1 + NP + CONFIG::D),1,"COO");

                for (size_t j = 0; j < NP; j++) {
                    INT glob_j = connectivityAB_[row][j];

                    auto d = norm(data_points[row].first - ptsExtend_[glob_j]);

                    if (d < r_) {
                        Aas_coo.set_value(j, 0, rbf(d), false);
//                        Aas.set_value(j, 0, rbf(d));
                    }
                }

                Aas_coo.set_value(NP, 0, 1, false);
//                Aas.set_value(NP, 0, 1);

                for (int dim = 0; dim < CONFIG::D; dim++) {
                    Aas_coo.set_value((NP + dim + 1), 0, data_points[row].first[dim], false);
//                    Aas.set_value((NP + dim + 1), 0, data_points[row].first[dim]);
                }
                Aas = Aas_coo;

                linalg::sparse_matrix<INT, REAL> H_i;

                if (precond_ == 0) {
                    linalg::conjugate_gradient<INT, REAL> cg(Css, Aas, cgSolveTol_, cgMaxIter_);
                    iterErrorReturn = cg.solve();
                    H_i = cg.getSolution();
                }
                else if (precond_ == 1) {
                    linalg::diagonal_preconditioner<INT, REAL> M(Css);
                    linalg::conjugate_gradient<INT, REAL> cg(Css, Aas, cgSolveTol_, cgMaxIter_, &M);
                    iterErrorReturn = cg.solve();
                    H_i = cg.getSolution();
                }
                else {
                    std::cerr
                            << "MUI Error [sampler_rbf.h]: Invalid CG Preconditioner function number ("
                            << precond_ << ")" << std::endl
                            << "Please set the CG Preconditioner function number (precond_) as: "
                            << std::endl << "precond_=0 (No Preconditioner); "
                            << std::endl << "precond_=1 (Diagonal Preconditioner); " << std::endl;
                    std::abort();
                }

                if (DEBUG) {
                    std::cout << "MUI [sampler_rbf.h]: #iterations of H_i:     " << iterErrorReturn.first
                            << ". Error of H_i: " << iterErrorReturn.second
                            << std::endl;
                }

                errorReturn += iterErrorReturn.second;

                if (smoothing) {
                    for (size_t j = 0; j < NP; j++) {
                        INT glob_j = connectivityAB_[row][j];
                        H_toSmooth_.set_value(glob_j, row, H_i.get_value(j, 0));
                    }
                }
                else {
                    for (size_t j = 0; j < NP; j++) {
                        INT glob_j = connectivityAB_[row][j];
                        H_.set_value(glob_j, row, H_i.get_value(j, 0));
                    }
                }
            }

            errorReturn /= static_cast<REAL>(data_points.size());

            if (smoothing) {
                for (size_t row = 0; row < data_points.size(); row++) {
                    for (size_t j = 0; j < NP; j++) {
                        INT row_i = connectivityAB_[row][j];
                        REAL h_j_sum = 0.;
                        REAL f_sum = 0.;

                        for (size_t k = 0; k < MP; k++) {
                            INT row_k = connectivityAA_[row_i][k];
                            if (row_k == static_cast<INT>(row_i)) {
                                std::cerr << "MUI Error [sampler_rbf.h]: Invalid row_k value: " << row_k
                                        << std::endl;
                            }
                            else
                                h_j_sum += std::pow(dist_h_i(row_i, row_k), -2.);
                        }

                        for (size_t k = 0; k < MP; k++) {
                            INT row_k = connectivityAA_[row_i][k];
                            if (row_k == static_cast<INT>(row_i)) {
                                std::cerr << "MUI Error [sampler_rbf.h]: Invalid row_k value: " << row_k
                                        << std::endl;
                            }
                            else {
                                REAL w_i = ((std::pow(dist_h_i(row_i, row_k), -2.))
                                        / (h_j_sum));
                                f_sum += w_i * H_toSmooth_.get_value(row_k, row);
                            }
                        }
                        H_.set_value(row_i, row,
                                (0.5 * (f_sum + H_toSmooth_.get_value(row_i, row))));
                    }
                }
            }
        }
        else { // Not using partitioned approach
            linalg::sparse_matrix<INT,REAL> Css; //< Matrix of radial basis function evaluations between prescribed points
            linalg::sparse_matrix<INT,REAL> Aas; //< Matrix of RBF evaluations between prescribed and interpolation points

            Css.resize((1 + ptsExtend_.size() + CONFIG::D), (1 + ptsExtend_.size() + CONFIG::D));
            Aas.resize(data_points.size(), (1 + ptsExtend_.size() + CONFIG::D));

            //set Css
            // Define intermediate vectors for performance purpose
            linalg::sparse_matrix<INT, REAL> Css_coo((1 + ptsExtend_.size() + CONFIG::D),(1 + ptsExtend_.size() + CONFIG::D),"COO");

            for ( size_t i = 0; i < ptsExtend_.size(); i++ ) {
                for ( size_t j = i; j < ptsExtend_.size(); j++ ) {
                    auto d = norm(ptsExtend_[i] - ptsExtend_[j]);

                    if ( d < r_ ) {
                        REAL w = rbf(d);
                        Css_coo.set_value((i + CONFIG::D + 1), (j + CONFIG::D + 1), w, false);
//                        Css.set_value((i + CONFIG::D + 1), (j + CONFIG::D + 1), w);

                        if ( i != j ) {
                            Css_coo.set_value((j + CONFIG::D + 1), (i + CONFIG::D + 1), w, false);
//                            Css.set_value((j + CONFIG::D + 1), (i + CONFIG::D + 1), w);
                        }
                    }
                }
            }
            Css = Css_coo;

            //set Aas
            // Define intermediate vectors for performance purpose
            linalg::sparse_matrix<INT, REAL> Aas_coo(data_points.size(),(1 + ptsExtend_.size() + CONFIG::D),"COO");

            for ( size_t i = 0; i < data_points.size(); i++ ) {
                for ( size_t j = 0; j < ptsExtend_.size(); j++ ) {
                    auto d = norm(data_points[i].first - ptsExtend_[j]);

                    if ( d < r_ ) {
                        Aas_coo.set_value(i, (j + CONFIG::D + 1), rbf(d), false);
//                        Aas.set_value(i, (j + CONFIG::D + 1), rbf(d));
                    }
                }
            }
            Aas = Aas_coo;

            linalg::sparse_matrix<INT,REAL> Aas_trans = Aas.transpose();

            linalg::sparse_matrix<INT, REAL> H_more;

            if (precond_ == 0) {
                linalg::conjugate_gradient<INT, REAL> cg(Css, Aas_trans, cgSolveTol_, cgMaxIter_);
                iterErrorReturn = cg.solve();
                H_more = cg.getSolution();
            }
            else if (precond_ == 1) {
                linalg::diagonal_preconditioner<INT, REAL> M(Css);
                linalg::conjugate_gradient<INT, REAL> cg(Css, Aas_trans, cgSolveTol_, cgMaxIter_, &M);
                iterErrorReturn = cg.solve();
                H_more = cg.getSolution();
            }
            else {
                std::cerr
                        << "MUI Error [sampler_rbf.h]: Invalid CG Preconditioner function number ("
                        << precond_ << ")" << std::endl
                        << "Please set the CG Preconditioner function number (precond_) as: "
                        << std::endl << "precond_=0 (No Preconditioner); "
                        << std::endl << "precond_=1 (Diagonal Preconditioner); " << std::endl;
                std::abort();
            }

            if (DEBUG) {
                std::cout << "MUI [sampler_rbf.h]: #iterations of H_more:     " << iterErrorReturn.first
                        << ". Error of H_more: " << iterErrorReturn.second
                        << std::endl;
            }

            errorReturn = iterErrorReturn.second;

            if ( smoothing ) {
                for ( size_t i = 0; i < ptsExtend_.size(); i++ ) {
                    for (size_t j = 0; j < data_points.size(); j++ ) {
                        H_toSmooth_.set_value(i, j, H_more.get_value((i + CONFIG::D + 1), j));
                    }
                }
            }
            else {
                for ( size_t i = 0; i < ptsExtend_.size(); i++ ) {
                    for ( size_t j = 0; j < data_points.size(); j++ ) {
                        H_.set_value(i, j, H_more.get_value((i + CONFIG::D + 1), j));
                    }
                }
            }

            if ( smoothing ) {
                for ( size_t row = 0; row < ptsExtend_.size(); row++ ) {
                    for ( size_t j = 0; j < data_points.size(); j++ ) {
                        REAL h_j_sum = 0.;
                        REAL f_sum = 0.;
                        for ( size_t k = 0; k < MP; k++ ) {
                            INT row_k = connectivityAA_[row][k];
                            if ( row_k == static_cast<INT>(row) )
                                std::cerr << "Invalid row_k value: " << row_k << std::endl;
                            else
                                h_j_sum += std::pow(dist_h_i(row, row_k), -2.);
                        }

                        for ( size_t k = 0; k < MP; k++ ) {
                            INT row_k = connectivityAA_[row][k];
                            if ( row_k == static_cast<INT>(row) )
                                std::cerr << "Invalid row_k value: " << row_k << std::endl;
                            else {
                                REAL w_i = ((std::pow(dist_h_i(row, row_k), -2.)) / (h_j_sum));
                                f_sum += w_i * H_toSmooth_.get_value(row_k, j);
                            }
                        }
                        H_.set_value(row, j, 0.5 * (f_sum + H_toSmooth_.get_value(row, j)));
                    }
                }
            }
        }

        return errorReturn;
    }

    template<template<typename, typename > class CONTAINER>
    inline void buildConnectivityConsistent(const CONTAINER<ITYPE, CONFIG> &data_points, const size_t NP, bool writeMatrix,
            const std::string& fileAddress) const {
        std::ofstream outputFileCAB;
        std::ofstream outputFileRemotePoints;
        if (writeMatrix) {
            outputFileCAB.open(fileAddress + "/connectivityAB.dat");

            if (!outputFileCAB) {
                std::cerr
                        << "MUI Error [sampler_rbf.h]: Could not locate the file address on the connectivityAB.dat"
                        << std::endl;
            }
            else {
                outputFileCAB
                        << "// ************************************************************************************************";
                outputFileCAB << "\n";
                outputFileCAB
                        << "// **** This is the 'connectivityAB.dat' file of the RBF spatial sampler of the MUI library";
                outputFileCAB << "\n";
                outputFileCAB
                        << "// **** This file contains the entire matrix of the Point Connectivity Matrix (N).";
                outputFileCAB << "\n";
                outputFileCAB
                        << "// **** The file uses the Comma-Separated Values (CSV) format with ASCII for the entire N matrix";
                outputFileCAB << "\n";
                outputFileCAB
                        << "// ************************************************************************************************";
                outputFileCAB << "\n";
                outputFileCAB << "// ";
                outputFileCAB << "\n";
            }

            outputFileRemotePoints.open(fileAddress + "/remotePoints.dat");

            if (!outputFileRemotePoints) {
                std::cerr
                        << "MUI Error [sampler_rbf.h]: Could not locate the file address on the remotePoints.dat"
                        << std::endl;
            }
            else {
                outputFileRemotePoints
                        << "// ************************************************************************************************";
                outputFileRemotePoints << "\n";
                outputFileRemotePoints
                        << "// **** This is the 'remotePoints.dat' file of the RBF spatial sampler of the MUI library";
                outputFileRemotePoints << "\n";
                outputFileRemotePoints
                        << "// **** This file contains the remote points in the correct order that build the coupling matrix H.";
                outputFileRemotePoints << "\n";
                outputFileRemotePoints
                        << "// **** The file uses the Comma-Separated Values (CSV) format with ASCII for the entire N matrix";
                outputFileRemotePoints << "\n";
                outputFileRemotePoints
                        << "// ************************************************************************************************";
                outputFileRemotePoints << "\n";
                outputFileRemotePoints << "// ";
                outputFileRemotePoints << "\n";
            }
        }

        connectivityAB_.resize(ptsExtend_.size());

        for (size_t i = 0; i < ptsExtend_.size(); i++) {
            INT pointsCount = 0;
            for (size_t n = 0; n < NP; n++) {
                REAL cur = std::numeric_limits<REAL>::max();
                INT bestj = -1;
                for (size_t j = 0; j < data_points.size(); j++) {
                    auto added = std::find_if(connectivityAB_[i].begin(),
                            connectivityAB_[i].end(), [j](INT k) {
                                return static_cast<size_t>(k) == j;
                            });

                    if (added != connectivityAB_[i].end())
                        continue;

                    auto d = normsq(ptsExtend_[i] - data_points[j].first);
                    if (d < cur) {
                        cur = d;
                        bestj = j;
                    }

                    if (n == 0 && d < twor_)
                        pointsCount++;
                }

                connectivityAB_[i].emplace_back(bestj);

                if (writeMatrix && (n < NP - 1))
                    outputFileCAB << bestj << ",";
                else if (writeMatrix)
                    outputFileCAB << bestj;
            }

            if (writeMatrix && i < ptsExtend_.size() - 1)
                outputFileCAB << '\n';
        }

        for (size_t i = 0; i < data_points.size(); i++) {

            point_type pointTemp;
            for (INT dim = 0; dim < CONFIG::D; dim++) {
                pointTemp[dim] = data_points[i].first[dim];
            }

            remote_pts_.push_back(pointTemp);

            if (writeMatrix) {
                for (INT dim = 0; dim < CONFIG::D; dim++) {
                    outputFileRemotePoints << remote_pts_[i][dim];
                    if (dim == (CONFIG::D - 1)) {
                        outputFileRemotePoints << '\n';
                    } else {
                        outputFileRemotePoints << ",";
                    }
                }
            }
        }

        if (writeMatrix) {
            outputFileCAB.close();
            outputFileRemotePoints.close();
        }
    }

    template<template<typename, typename > class CONTAINER>
    inline void buildConnectivityConservative(const CONTAINER<ITYPE, CONFIG> &data_points, const size_t NP, bool writeMatrix,
            const std::string& fileAddress) const {
        std::ofstream outputFileCAB;
        std::ofstream outputFileRemotePoints;
        if (writeMatrix) {
            outputFileCAB.open(fileAddress + "/connectivityAB.dat");

            if (!outputFileCAB) {
                std::cerr
                        << "MUI Error [sampler_rbf.h]: Could not locate the file address on the connectivityAB.dat"
                        << std::endl;
            }
            else {
                outputFileCAB
                        << "// ************************************************************************************************";
                outputFileCAB << "\n";
                outputFileCAB
                        << "// **** This is the 'connectivityAB.dat' file of the RBF spatial sampler of the MUI library";
                outputFileCAB << "\n";
                outputFileCAB
                        << "// **** This file contains the entire matrix of the Point Connectivity Matrix (N).";
                outputFileCAB << "\n";
                outputFileCAB
                        << "// **** The file uses the Comma-Separated Values (CSV) format with ASCII for the entire N matrix";
                outputFileCAB << "\n";
                outputFileCAB
                        << "// ************************************************************************************************";
                outputFileCAB << "\n";
                outputFileCAB << "// ";
                outputFileCAB << "\n";
            }

            outputFileRemotePoints.open(fileAddress + "/remotePoints.dat");

            if (!outputFileRemotePoints) {
                std::cerr
                        << "MUI Error [sampler_rbf.h]: Could not locate the file address on the remotePoints.dat"
                        << std::endl;
            }
            else {
                outputFileRemotePoints
                        << "// ************************************************************************************************";
                outputFileRemotePoints << "\n";
                outputFileRemotePoints
                        << "// **** This is the 'remotePoints.dat' file of the RBF spatial sampler of the MUI library";
                outputFileRemotePoints << "\n";
                outputFileRemotePoints
                        << "// **** This file contains the remote points in the correct order that build the coupling matrix H.";
                outputFileRemotePoints << "\n";
                outputFileRemotePoints
                        << "// **** The file uses the Comma-Separated Values (CSV) format with ASCII for the entire N matrix";
                outputFileRemotePoints << "\n";
                outputFileRemotePoints
                        << "// ************************************************************************************************";
                outputFileRemotePoints << "\n";
                outputFileRemotePoints << "// ";
                outputFileRemotePoints << "\n";
            }
        }

        connectivityAB_.resize(data_points.size());

        for (size_t i = 0; i < data_points.size(); i++) {
            INT pointsCount = 0;
            for (size_t n = 0; n < NP; n++) {
                REAL cur = std::numeric_limits<REAL>::max();
                INT bestj = -1;
                for (size_t j = 0; j < ptsExtend_.size(); j++) {
                    auto added = std::find_if(connectivityAB_[i].begin(),
                            connectivityAB_[i].end(), [j](INT k) {
                                return static_cast<size_t>(k) == j;
                            });

                    if (added != connectivityAB_[i].end())
                        continue;

                    auto d = normsq(data_points[i].first - ptsExtend_[j]);
                    if (d < cur) {
                        cur = d;
                        bestj = j;
                    }

                    if (n == 0 && d < twor_)
                        pointsCount++;
                }

                connectivityAB_[i].emplace_back(bestj);

                if (writeMatrix && n < NP - 1)
                    outputFileCAB << bestj << ",";
                else if (writeMatrix)
                    outputFileCAB << bestj;
            }

            if (writeMatrix && i < ptsExtend_.size() - 1)
                outputFileCAB << '\n';

            point_type pointTemp;
            for (INT dim = 0; dim < CONFIG::D; dim++) {
                pointTemp[dim] = data_points[i].first[dim];
            }
            remote_pts_.push_back(pointTemp);

            if (writeMatrix) {
                for (INT dim = 0; dim < CONFIG::D; dim++) {
                    outputFileRemotePoints << remote_pts_[i][dim];
                    if (dim == (CONFIG::D - 1)) {
                        outputFileRemotePoints << '\n';
                    } else {
                        outputFileRemotePoints << ",";
                    }
                }
            }
        }

        if (writeMatrix) {
            outputFileCAB.close();
            outputFileRemotePoints.close();
        }
    }

    void buildConnectivitySmooth(const size_t MP, bool writeMatrix, const std::string& fileAddress) const {
        std::ofstream outputFileCAA;

        if (writeMatrix) {
            outputFileCAA.open(fileAddress + "/connectivityAA.dat");

            if (!outputFileCAA) {
                std::cerr
                        << "MUI Error [sampler_rbf.h]: Could not locate the file address on the connectivityAA.dat!"
                        << std::endl;
            }
            else {
                outputFileCAA
                        << "// ************************************************************************************************";
                outputFileCAA << "\n";
                outputFileCAA
                        << "// **** This is the 'connectivityAA.dat' file of the RBF spatial sampler of the MUI library";
                outputFileCAA << "\n";
                outputFileCAA
                        << "// **** This file contains the entire matrix of the Point Connectivity Matrix (M) (for smoothing).";
                outputFileCAA << "\n";
                outputFileCAA
                        << "// **** The file uses the Comma-Separated Values (CSV) format with ASCII for the entire N matrix";
                outputFileCAA << "\n";
                outputFileCAA
                        << "// ************************************************************************************************";
                outputFileCAA << "\n";
                outputFileCAA << "// ";
                outputFileCAA << "\n";
            }
        }

        connectivityAA_.resize(ptsExtend_.size());

        for (size_t i = 0; i < ptsExtend_.size(); i++) {
            for (size_t n = 0; n < MP; n++) {
                REAL cur = std::numeric_limits<REAL>::max();
                INT bestj = -1;
                for (size_t j = 0; j < ptsExtend_.size(); j++) {
                    if (i == j)
                        continue;

                    auto added = std::find_if(connectivityAA_[i].begin(),
                            connectivityAA_[i].end(), [j](INT i) {
                                return static_cast<size_t>(i) == j;
                            });

                    if (added != connectivityAA_[i].end())
                        continue;

                    auto d = normsq(ptsExtend_[i] - ptsExtend_[j]);
                    if (d < cur) {
                        cur = d;
                        bestj = j;
                    }
                }

                connectivityAA_[i].emplace_back(bestj);

                if (writeMatrix && n < MP - 1)
                    outputFileCAA << bestj << ",";
                else if (writeMatrix)
                    outputFileCAA << bestj;
            }

            if (writeMatrix && i < ptsExtend_.size() - 1)
                outputFileCAA << '\n';
        }

        if (writeMatrix)
            outputFileCAA.close();
    }

    //Radial basis function for two points
    inline REAL rbf(point_type x1, point_type x2) const {
        auto d = norm(x1 - x2);
        return rbf(d);
    }

    //Radial basis function for calculated distance
    inline REAL rbf(REAL d) const {
        switch (basisFunc_) {
        case 0:
            //Gaussian
            return (d < r_) ? std::exp(-(s_ * s_ * d * d)) : 0.;
        case 1:
            //Wendland's C0
            return std::pow((1. - d), 2.);
        case 2:
            //Wendland's C2
            return (std::pow((1. - d), 4.)) * ((4. * d) + 1.);
        case 3:
            //Wendland's C4
            return (std::pow((1. - d), 6)) * ((35. * d * d) + (18. * d) + 3.);
        case 4:
            //Wendland's C6
            return (std::pow((1. - d), 8.))
                    * ((32. * d * d * d) + (25. * d * d) + (8. * d) + 1.);
        default:
            std::cerr
                    << "MUI Error [sampler_rbf.h]: invalid RBF basis function number ("
                    << basisFunc_ << ")" << std::endl
                    << "Please set the RBF basis function number (basisFunc_) as: "
                    << std::endl << "basisFunc_=0 (Gaussian); " << std::endl
                    << "basisFunc_=1 (Wendland's C0); " << std::endl
                    << "basisFunc_=2 (Wendland's C2); " << std::endl
                    << "basisFunc_=3 (Wendland's C4); " << std::endl
                    << "basisFunc_=4 (Wendland's C6); " << std::endl;
            return 0;
        }
    }

    ///Distances function
    inline REAL dist_h_i(INT ptsExtend_i, INT ptsExtend_j) const {
        switch (CONFIG::D) {
        case 1:
            return std::sqrt((std::pow((ptsExtend_[ptsExtend_i][0]
                             - ptsExtend_[ptsExtend_j][0]), 2.)));
        case 2:
            return std::sqrt((std::pow((ptsExtend_[ptsExtend_i][0]
                             - ptsExtend_[ptsExtend_j][0]), 2.))
                             + (std::pow((ptsExtend_[ptsExtend_i][1]
                                          - ptsExtend_[ptsExtend_j][1]), 2.)));
        case 3:
            return std::sqrt((std::pow((ptsExtend_[ptsExtend_i][0]
                             - ptsExtend_[ptsExtend_j][0]), 2.))
                             + (std::pow((ptsExtend_[ptsExtend_i][1]
                             - ptsExtend_[ptsExtend_j][1]), 2.))
                             + (std::pow((ptsExtend_[ptsExtend_i][2]
                              - ptsExtend_[ptsExtend_j][2]), 2.)));
        default:
            std::cerr << "MUI Error [sampler_rbf.h]: CONFIG::D must equal 1-3" << std::endl;
            return 0.;
        }
    }

protected:
    REAL r_;
    const std::vector<point_type> pts_; //< Local points
    const INT basisFunc_;
    const bool conservative_;
    const bool consistent_;
    const bool smoothFunc_;
    const bool generateMatrix_;
    const std::string writeFileAddress_;
    INT precond_;
    mutable bool initialised_;
    MPI_Comm local_mpi_comm_world_;
    mutable INT CABrow_;
    mutable INT CABcol_;
    mutable INT Hrow_;
    mutable INT Hcol_;
    bool pouEnabled_;
    REAL cgSolveTol_;
    INT cgMaxIter_;

    mutable size_t N_sp_;
    mutable size_t M_ap_;
    int local_rank_;
    int local_size_;

    REAL twor_;
    REAL s_;
    mutable INT CAArow_;
    mutable INT CAAcol_;
    mutable INT remote_pts_num_;
    mutable INT remote_pts_dim_;
    mutable std::vector<point_type> ptsGhost_; //< Local ghost points
    mutable std::vector<point_type> ptsExtend_; //< Extended local points, i.e. local points and ghost local points
    mutable std::vector<std::vector<INT> > connectivityAB_;
    mutable std::vector<std::vector<INT> > connectivityAA_;
    mutable linalg::sparse_matrix<INT, REAL> H_; //< Transformation Matrix
    mutable linalg::sparse_matrix<INT, REAL> H_toSmooth_;
    mutable std::vector<point_type> remote_pts_; //< Remote points
};
} // mui

#endif /* MUI_SAMPLER_RBF_H_ */
