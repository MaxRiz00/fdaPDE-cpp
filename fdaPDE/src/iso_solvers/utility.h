// This file is part of fdaPDE, a C++ library for physics-informed
// spatial and functional data analysis.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __FDAPDE_DRIVERS_UTILITY_H__
#define __FDAPDE_DRIVERS_UTILITY_H__

namespace fdapde {
namespace internals {

// checks if supplied penalty is valid
template <typename Penalty> class is_valid_penalty_pair {
    using Penalty_ = std::decay_t<Penalty>;
   public:
    static constexpr bool value =
      is_pair_v<Penalty> &&
      requires(std::tuple_element_t<0, Penalty_> t) {   // first pair element: bilinear form
          { t.assemble() } -> std::same_as<Eigen::SparseMatrix<double>>;
      } &&
      requires(std::tuple_element_t<1, Penalty_> t) {   // second pair element: linear form
          { t.assemble() } -> std::same_as<Eigen::Matrix<double, Dynamic, 1>>;
      } &&   // linear and bilinear form have same discretization category
      std::is_same_v<typename std::tuple_element_t<0, Penalty_>::discretization_category FDAPDE_COMMA
                     typename std::tuple_element_t<1, Penalty_>::discretization_category>;
};
template <typename Penalty> constexpr bool is_valid_penalty_pair_v = is_valid_penalty_pair<Penalty>::value;
  
// efficient left multiplication Q*x, with Q = W * (I - X * (X^\top * W * X)^{-1} * X^\top * W)
template <typename WeightMatrix, typename DesignMatrix, typename InvDesignMatrix>
Eigen::Matrix<double, Dynamic, Dynamic> lmbQ(
  const WeightMatrix& W, const DesignMatrix& X, const InvDesignMatrix& invXtWX,
  const Eigen::Matrix<double, Dynamic, Dynamic>& x) {
    using MatrixType = Eigen::Matrix<double, Dynamic, Dynamic>;
    if (X.cols() == 0) return W * x;
    MatrixType v = X.transpose() * W * x;   // X^\top * W * x
    MatrixType z = invXtWX.solve(v);        // (X^\top * W * X)^{-1} * X^\top * W * x
    // compute W*x - W*X*z = W*x - (W*X*(X^\top*W*X)^{-1}*X^\top*W)*x = W(I - H)*x = Q*x
    return W * (x - X * z);
}
template <typename DesignMatrix, typename InvDesignMatrix>
Eigen::Matrix<double, Dynamic, Dynamic>
lmbQ(const DesignMatrix& X, const InvDesignMatrix& invXtX, const Eigen::Matrix<double, Dynamic, Dynamic>& x) {
    using MatrixType = Eigen::Matrix<double, Dynamic, Dynamic>;
    if (X.cols() == 0) return x;
    MatrixType v = X.transpose() * x;   // X^\top * x
    MatrixType z = invXtX.solve(v);     // (X^\top * X)^{-1} * X^\top * x
    // compute x - X*z = x - (X*(X^\top*X)^{-1}*X^\top)*x = (I - H)*x = Q*x
    return x - X * z;
}

// pointwise basis evaluation for finite element basis system
template <typename IsoMesh_, typename CoordsMatrix_>
    requires(internals::is_eigen_dense_xpr_v<CoordsMatrix_>)
Eigen::SparseMatrix<double> point_basis_eval(IsoSpace<IsoMesh_> iso_space, CoordsMatrix_&& coords) {
    //std::cout << "POIT Evaluating basis at " << coords.rows() << " locations." << std::endl;
    static constexpr int local_dim = IsoMesh_::local_dim;
    static constexpr int embed_dim = IsoMesh_::embed_dim;
    fdapde_assert(coords.rows() > 0 && coords.cols() == embed_dim);

    auto& dof_handler = iso_space.dof_handler();

    int n_shape_functions = iso_space.n_shape_functions();
    int n_mapped_dofs = dof_handler.n_mapped_dofs();
    int n_extended_dofs = dof_handler.n_dofs();
    std::vector<int> dof_map = dof_handler.dof_map(); // for periodicity

    int n_locs = coords.rows();
    Eigen::SparseMatrix<double> psi_(n_locs, n_extended_dofs);
    // evaluate basis system at locations
    std::vector<fdapde::Triplet<double>> triplet_list;
    triplet_list.reserve(n_locs * n_extended_dofs);

    //std::cout << "Number of shape functions: " << n_shape_functions << std::endl;
    //std::cout << "Number of mapped dofs: " << n_mapped_dofs << std::endl;
    //std::cout << "Points: "<<coords << std::endl;



    Eigen::Matrix<int, Dynamic, 1> cell_id(coords.rows());

    std::vector<Eigen::Matrix<double, local_dim, 1>> param_coords(coords.rows());
    for(int i = 0; i<coords.rows();i++) {
        //std::cout << "Locating point " << i << ": " << coords.row(i) << std::endl;
        param_coords[i] = iso_space.mesh().invert_point(coords.row(i),5);
        cell_id[i] = iso_space.mesh().locate_param(param_coords[i]); // to improve
        //std::cout << "Point " << coords.row(i) << " located in cell ID: " << cell_id[i] << std::endl;
    }
    

    //std::cout<<"shape of cell_id: " << cell_id.rows() << " x " << cell_id.cols() << std::endl;
    //std::cout<<"shape of coords: " << coords.rows() << " x " << coords.cols() << std::endl;
    // build basis evaluation matrix
    for (int i = 0; i < n_locs; ++i) {
        if (cell_id[i] != -1) {   // point falls inside domain
            Eigen::Matrix<double, embed_dim, 1> p_i(coords.row(i));
            auto param_p_i = param_coords[i];
            //std::cout<< "Evaluating basis at point " << p_i.transpose() << " in cell ID: " << cell_id[i] << std::endl;
            //std::cout << "Parametric coordinates: " << param_p_i.transpose() << std::endl;
            auto cell = dof_handler.cell(cell_id[i]);
            for (int h = 0; h < cell.dofs().size(); ++h) {
                int active_dof = cell.dofs()[h];
                triplet_list.emplace_back(
                  i, active_dof, iso_space.eval_shape_value(active_dof, param_p_i));   // \psi_j(p_i) da capire che se fare la trasparmaziona linearecell.invJ() * (p_i - cell.node(0))
            }
        }
    }
    //std::cout << "Triplet list size: " << triplet_list.size() << std::endl;
    //std::cout <<"n_locs x n_dofs: " << n_locs << " x " << n_dofs << std::endl;
    // finalize construction
    psi_.setFromTriplets(triplet_list.begin(), triplet_list.end()); 

    Eigen::SparseMatrix<double> psi_reduced(n_locs, n_mapped_dofs);

    std::vector<Eigen::Triplet<double>> reduced_triplet_list;
    reduced_triplet_list.reserve(triplet_list.size());  // same size, roughly

    for (int i = 0; i < n_locs; ++i) {
        // Copy the sparse row into a dense vector
        Eigen::VectorXd row_dense = psi_.row(i);

        // Apply periodic constraints
        Eigen::VectorXd constrained_row = row_dense;
        dof_handler.enforce_periodic_constraints(constrained_row);

        // Insert nonzero entries into the reduced triplet list
        for (int j = 0; j < constrained_row.size(); ++j) {
            if (constrained_row(j) != 0.0) {
                reduced_triplet_list.emplace_back(i, j, constrained_row(j));
            }
        }
    }

    // Now build psi_reduced
    psi_reduced.setFromTriplets(reduced_triplet_list.begin(), reduced_triplet_list.end());


    psi_.makeCompressed();
    psi_reduced.makeCompressed();

    //std::cout << "Basis evaluation completed." << std::endl;
    //std::cout << "psi_reduced size: " << psi_reduced.rows() << " x " << psi_reduced.cols() << std::endl;
    std::cout << "psi_reduced size: " << psi_.rows() << " x " << psi_.cols() << std::endl;
    return psi_reduced;
}
template <typename IsoMesh_, typename GeoIndex_>
    requires(!internals::is_eigen_dense_xpr_v<GeoIndex_>)
Eigen::SparseMatrix<double>
point_basis_eval(const IsoSpace<IsoMesh_>& iso_space, const GeoIndex_& geo_index) {
    //std::cout << "Evaluating basis at " << geo_index.rows() << " locations." << std::endl;
    if (geo_index.points_at_dofs()) {
        int n_dofs = iso_space.n_dofs();
        int n_locs = geo_index.rows();
        //std::cout << "Evaluating basis at " << n_locs << " locations, with " << n_dofs << " dofs." << std::endl;
        Eigen::SparseMatrix<double> psi_(n_locs, n_dofs);
        psi_.setIdentity();   // \psi_i(p_j) = 1 \iff i == j, otherwise \psi_i(p_j) = 0
        return psi_;
    }
    //std::cout << "Evaluating basis at " << geo_index.rows() << " locations." << std::endl;
    return point_basis_eval(iso_space, geo_index.coordinates());
}


// areal basis evaluation for finite element basis system
/*
template <typename IsoMesh_>
std::pair<Eigen::SparseMatrix<double>, Eigen::Matrix<double, Dynamic, 1>> areal_basis_eval(
  const IsoSpace<IsoMesh_>& iso_space, const BinaryMatrix<Dynamic, Dynamic>& incidence_mat) {
    using IsoSpace_ = IsoSpace<IsoMesh_>;
    fdapde_assert(incidence_mat.rows() > 0 && incidence_mat.cols() == iso_space.mesh().n_cells());
    static constexpr int local_dim = IsoMesh_::local_dim;
    static constexpr int embed_dim = IsoMesh_::embed_dim;
    using IsoType = typename IsoSpace_::IsoType;
    using cell_dof_descriptor = typename IsoSpace_::cell_dof_descriptor;
    using BasisType = typename cell_dof_descriptor::BasisType;
    using Quadrature = typename IsoType::template cell_quadrature_t<local_dim>;
    static constexpr int n_quadrature_nodes = Quadrature::order;
    static constexpr int n_shape_functions = iso_space.n_shape_functions();
    // compile time evaluation of \int_{\hat K} \psi_i on reference element \hat K
    static constexpr Matrix<double, n_shape_functions, 1> int_table_ {[]() {
        std::array<double, n_shape_functions> int_table_ {};
        BasisType basis {cell_dof_descriptor().dofs_phys_coords()};
        for (int i = 0; i < n_shape_functions; ++i) {
            for (int k = 0; k < n_quadrature_nodes; ++k) {
                int_table_[i] += Quadrature::weights[k] * basis[i](Quadrature::nodes.row(k).transpose());
            }
        }
        return int_table_;
    }};

    int n_dofs = iso_space.n_dofs();
    int n_regions = incidence_mat.rows();
    Eigen::SparseMatrix<double> psi_(n_regions, n_dofs);
    Eigen::Matrix<double, Dynamic, 1> D(n_regions);
    std::vector<fdapde::Triplet<double>> triplet_list;
    triplet_list.reserve(n_regions * n_shape_functions);

    const auto& dof_handler = iso_space.dof_handler();
    int tail = 0;
    for (int k = 0; k < n_regions; ++k) {
        int head = 0;
        double Di = 0;   // measure of region D_i
        for (int l = 0, n_cells = incidence_mat.cols(); l < n_cells; ++l) {
            if (incidence_mat(k, l)) {   // element with ID l belongs to k-th region
                auto cell = dof_handler.cell(l);
                for (int h = 0; h < n_shape_functions; ++h) {
                    // compute \int_e \psi_h on physical element e
                    triplet_list.emplace_back(k, cell.dofs()[h], int_table_[h] * cell.parametric_measure());
                    head++;
                }
                Di += cell.parametric_measure();
            }
        }
        // divide each \int_{D_i} \psi_j by the measure of region D_i
        for (int j = 0; j < head; ++j) { triplet_list[tail + j].value() /= Di; }
        D[k] = Di;
        tail += head;
    }
    // finalize construction
    psi_.setFromTriplets(triplet_list.begin(), triplet_list.end());
    psi_.makeCompressed();
    return std::make_pair(std::move(psi_), std::move(D));
}
template <typename IsoMesh_, typename GeoIndex_>
std::pair<Eigen::SparseMatrix<double>, Eigen::Matrix<double, Dynamic, 1>>
areal_basis_eval(const IsoSpace<IsoMesh_>& iso_space, const GeoIndex_& geo_index) {
    return areal_basis_eval(iso_space, geo_index.incidence_matrix());
}
    */


}   // namespace internals
}   // namespace fdapde

#endif // __FE_ELLIPTIC_DRIVER_H__
