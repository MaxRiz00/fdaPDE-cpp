#include "iso_solvers.h"
#include "regression.h"
#include "../helpers.h"
#include <random>
#include <unordered_map>
#include <fstream>
#include <Eigen/Dense>
#include <vector>
#include <sstream>
#include <limits>

typedef Eigen::Matrix<double, 3, 1> Vec;

struct WindPoint {
    Vec pos;
    Vec wind;
};

// Read 3D points from CSV
auto read_csv_points(const std::string& filename) {
    std::vector<Vec> points;
    std::ifstream file(filename);
    std::string line;
    std::getline(file, line); // skip header

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string val;
        std::vector<double> coords;
        while (std::getline(ss, val, ',')) {
            coords.push_back(std::stod(val));
        }
        if (coords.size() == 3)
            points.emplace_back(coords[0], coords[1], coords[2]);
    }
    return points;
}

// Read scalar data (temperature)
auto read_csv_scalar(const std::string& filename, const std::string& colname) {
    std::ifstream file(filename);
    std::string line, header;
    std::getline(file, header); // skip header
    std::vector<double> data;
    while (std::getline(file, line)) {
        data.push_back(std::stod(line));
    }
    return data;
}

// Read wind data (x,y,z,wx,wy,wz)
auto read_wind_data(const std::string& filename) {
    std::vector<WindPoint> wind;
    std::ifstream file(filename);
    std::string line;
    std::getline(file, line); // skip header

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string val;
        std::vector<double> vals;
        while (std::getline(ss, val, ',')) {
            vals.push_back(std::stod(val));
        }
        if (vals.size() == 6)
            wind.push_back({Vec(vals[0], vals[1], vals[2]), Vec(vals[3], vals[4], vals[5])});
    }
    return wind;
}

// Find nearest wind point
auto nearest_wind(const Vec& point, const std::vector<decltype(read_wind_data(""))::value_type>& wind_data) {
    double min_dist = std::numeric_limits<double>::max();
    size_t min_idx = 0;
    for (size_t i = 0; i < wind_data.size(); ++i) {
        double dist = (point - wind_data[i].pos).norm();
        if (dist < min_dist) {
            min_dist = dist;
            min_idx = i;
        }
    }
    /*
    std::cout << "Nearest wind point found at point: " 
              << wind_data[min_idx].pos.transpose() 
              << " with wind vector: " 
              << wind_data[min_idx].wind.transpose() 
              << " at distance: " 
              << min_dist << std::endl;
              */
    return wind_data[min_idx];
}


Vec tangent_idw_interpolator(const Vec& p, const std::vector<WindPoint>& wind_data, double eps = 1e-7) {
    double wsum = 0.0;
    Vec interp = Vec::Zero();
    
    // Ensure p is normalized (assuming unit sphere)
    Vec p_norm = p.normalized();

    for (const auto& wp : wind_data) {
        // Ensure wind point position is normalized
        Vec wp_pos_norm = wp.pos.normalized();
        
        // Use geodesic distance on unit sphere for better accuracy
        double cos_angle = std::clamp(p_norm.dot(wp_pos_norm), -1.0, 1.0);
        double geodesic_dist = std::acos(cos_angle);
        
        if (geodesic_dist < eps) {
            // Exact match - project wind to tangent space at p
            Vec exact_wind = wp.wind - wp.wind.dot(p_norm) * p_norm;
            return exact_wind;
        }

        // IDW weight using geodesic distance
        double w = 1.0 / (geodesic_dist * geodesic_dist);
        
        // Ensure the wind vector is tangent at its source point
        Vec wind_tangent_at_source = wp.wind - wp.wind.dot(wp_pos_norm) * wp_pos_norm;
        
        // Simple parallel transport: project to tangent space at target point p
        Vec wind_at_p = wind_tangent_at_source - wind_tangent_at_source.dot(p_norm) * p_norm;
        
        interp += w * wind_at_p;
        wsum += w;
    }

    if (wsum > eps) {
        interp /= wsum;
    }

    // Final projection to ensure result is tangent at p
    interp -= (interp.dot(p_norm)) * p_norm;

    return interp;
}

using namespace fdapde;

int main() {
    auto mesh = IsoMesh<2, 3>::sphere_patch(1.0, 20.0, 60.0, -40.0, 50.0);  //60 EUROPE
    mesh.refine_knots({25,25});
    std::filesystem::create_directories("../results/");

    std::vector<double> c_values = {0,1,10,100,200}; //,100, 200, 300, 400, 500, 600, 700, 800
    int num_seeds = 10;

    std::string root_folder = "../results/";
    std::string wind_file = "../data/wind.csv";

    auto wind_data = read_wind_data(wind_file);

    std::ofstream lambda_out(root_folder + "lambdas.csv");
    lambda_out << "c,lambda_mean,lambda_std\n";

    for (double c : c_values) {
        std::vector<double> rmse_list;
        std::vector<double> lambda_list;

        std::cout<<"c = "<<c<<std::endl;

        for (int seed = 0; seed < num_seeds; ++seed) {
            std::stringstream train_file, resp_file, test_file, test_resp_file;
            train_file << "../data/train_locs_" << seed << ".csv";
            //train_file << "../dataRAW/train_locs.csv";
            resp_file << "../data/train_response_" << seed << ".csv";
            test_file << "../data/test_locs.csv";
            test_resp_file << "../data/test_response.csv";

            auto points = read_csv_points(train_file.str());
            GeoFrame data(mesh);
            auto& locs = data.insert_scalar_layer<POINT>("locs", train_file.str());
            locs.load_csv<double>(resp_file.str());
            
            
            MatrixField<3, 3, 1> b;

            // Define wind field using nearest neighbor interpolation
            
            b[0] = [&](const Vec& p) { return nearest_wind(p, wind_data).wind(0); };
            b[1] = [&](const Vec& p) { return nearest_wind(p, wind_data).wind(1); };
            b[2] = [&](const Vec& p) { return nearest_wind(p, wind_data).wind(2); };


            //for (int i = 0; i < 3; ++i)
            //  b[i] = [&, i](const Vec& p) { return tangent_idw_interpolator(p, wind_data)(i); };
               



            IsoSpace Vh(mesh);
            TrialFunction f(Vh);
            TestFunction v(Vh);

            auto a = integral(mesh, QGL2DP9)((laplacian(f) + c * dot(b, grad(f))) * (laplacian(v) +  c * dot(b, grad(v))));
            ZeroField<3> u;
            auto F = integral(mesh, QGL2DP9)(v * u);

            SRPDE m("target ~ f", data, iso_ls_elliptic(a, F));

            std::vector<double> lambda_grid;
            int n_lambdas = 100;
            double log_min = std::log10(1e-15);
            double log_max = std::log10(1e-5);
            for (int i = 0; i < n_lambdas; ++i) {
                double log_lambda = log_min + i * (log_max - log_min) / (n_lambdas - 1);
                lambda_grid.push_back(std::pow(10.0, log_lambda));
            }

            GridOptimizer<1> opt;
            opt.optimize(m.gcv(), lambda_grid);
            std::cout << "Optimal lambda: " << opt.optimum()[0] << "\n";
            double lambda = opt.optimum()[0];
            lambda_list.push_back(lambda);
            m.fit(opt.optimum()[0]);

            IsoFunction sol(Vh);
            sol = Vh.dof_handler().expand_solution(m.f());

            std::vector<Vec> test_points = read_csv_points(test_file.str());
            std::vector<double> true_values_vec = read_csv_scalar(test_resp_file.str(), "target");
            Eigen::VectorXd true_values = Eigen::Map<Eigen::VectorXd>(true_values_vec.data(), true_values_vec.size());

            Eigen::VectorXd predictions(test_points.size());
            for (size_t i = 0; i < test_points.size(); ++i)
                predictions(i) = sol(mesh.invert_point(test_points[i], 5));

            double mse = (predictions - true_values).squaredNorm() / predictions.size();
            double rmse = std::sqrt(mse);
            rmse_list.push_back(rmse);

            std::cout<<"RMSE: "<<rmse<<std::endl;

            std::ostringstream c_folder;
            c_folder << "c_" << std::scientific << std::setprecision(0) << c;
            std::string base_folder = root_folder + c_folder.str() + "/";
            std::filesystem::create_directories(base_folder);

            if (seed ==5 ) helpers::export_results(mesh, sol, base_folder, 10);

            /*
            ScalarField<3> phys_constr{[&](const Vec& p) {
                Vec grad = sol.phys_grad(mesh.invert_point(p, 5));
                auto b_val = b(p);
                Vec wind;
                for (int i = 0; i<3 ; ++i) wind(i) = b_val(i,0);
                // return inner product of gradient and wind vector
                return grad.dot(wind) * grad.dot(wind);
            }};

            double constraint = integral(mesh, QGL2DP9)(phys_constr );
            std::cout << "Constraint value for c=" << c << ", seed=" << seed << ": " << constraint << "\n";
            */
        }

        double lambda_sum = 0.0;
        for (double l : lambda_list) lambda_sum += l;
        double lambda_mean = lambda_sum / lambda_list.size();

        double lambda_sq_sum = 0.0;
        for (double l : lambda_list) lambda_sq_sum += (l - lambda_mean) * (l - lambda_mean);
        double lambda_std = std::sqrt(lambda_sq_sum / lambda_list.size());

        // Save to CSV
        lambda_out << c << "," << lambda_mean << "," << lambda_std << "\n";

        std::ostringstream name;
        name << std::scientific << std::setprecision(1) << c;
        std::ofstream fout(root_folder + "RMSE_" + name.str() + ".csv");
        fout << "seed,rmse\n";
        for (int i = 0; i < rmse_list.size(); ++i) {
            fout << i << "," << rmse_list[i] << "\n";
        }
        fout.close();
    }

    lambda_out.close();

    std::cout << "Grid search completed and RMSEs saved.\n";
}