#include "../sources/flow/momentum_iterator.hpp"
#include "../sources/L.hpp"
#include "../sources/mean_field/Model.hpp"
#include "../sources/flow/FlowContainer.hpp"

#include <Eigen/Dense>
#include <mrock/utility/InputFileReader.hpp>
#include <mrock/utility/OutputConvenience.hpp>
#include <nlohmann/json.hpp>

#include <array>
#include <iostream>

#ifndef OUTPUT_DATA_DIR
#define OUTPUT_DATA_DIR "build/" //"../../data/nickel_cut/"
#endif

using namespace NickelCUT;

typedef flow::momentum_iterator<L> mt;

Eigen::MatrixXd compute_M(const flow::FlowContainer& flow_result, 
    const std::array<double, N>& occupation_numbers, 
    const mt& x) 
{
    Eigen::MatrixXd dynamical_matrix = Eigen::MatrixXd::Zero(N, N);

    for (mt k = mt::begin(); k != mt::end(); ++k) {
        for (mt K = mt::begin(); K != mt::end(); ++K) {
            dynamical_matrix(k.get_position(), k.get_position()) += (flow_result.interactions_differing_spin(k+x,K,flow::Gamma<L>)
                + flow_result.interactions_differing_spin(k,K,flow::Gamma<L>)
                - flow_result.interactions_same_spin(k+x,K,K-k-x)
                + flow_result.interactions_same_spin(k+x,K,flow::Gamma<L>)
                - flow_result.interactions_same_spin(k,K,K-k-x)
                + flow_result.interactions_same_spin(k,K,flow::Gamma<L>)) * occupation_numbers[K];
        }
        dynamical_matrix(k.get_position(), k.get_position()) += 0.5 * (flow_result.epsilon_tilde[k] + flow_result.epsilon_tilde[k+x]);
        dynamical_matrix(k.get_position(), k.get_position()) *= 2. * (1. - occupation_numbers[k] - occupation_numbers[k+x]);

        for (mt l = mt::begin(); l != mt::end(); ++l) {
            dynamical_matrix(k.get_position(), l.get_position()) += 2 * flow_result.interactions_differing_spin(-k-x,k,k-l) * (
                (1. - occupation_numbers[k] - occupation_numbers[k+x]) * (1. - occupation_numbers[l] - occupation_numbers[l+x])
            );
        }
    }
    return dynamical_matrix;
} 

Eigen::DiagonalMatrix<double, Eigen::Dynamic> compute_N(const std::array<double, N>& occupation_numbers, 
    const mt& x) 
{
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> norm_matrix;
    norm_matrix.setZero(N);

    for (mt k = mt::begin(); k != mt::end(); ++k) {
        norm_matrix.diagonal()(k.get_position()) = occupation_numbers[k] + occupation_numbers[k+x] - 1.;
    }
    return norm_matrix;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Invalid number of arguments: Use <path_to_executable> <configfile>" << std::endl;
        return -1;
    }
    mrock::utility::InputFileReader input(argv[1]);

    const std::string load_binary_dir = std::string(OUTPUT_DATA_DIR) 
        + (std::string(OUTPUT_DATA_DIR).back() == '/' ? "" : "/") // ensures that OUTPUT_DATA_DIR ends in "/"
        + input.getString("output_dir") + "/"
        + "binaries/";
    mean_field::Model model(load_binary_dir, input);

    const std::string output_dir = std::string(OUTPUT_DATA_DIR) 
        + (std::string(OUTPUT_DATA_DIR).back() == '/' ? "" : "/") // ensures that OUTPUT_DATA_DIR ends in "/"
        + input.getString("output_dir") + "/"
        + model.data_dir_name();
    
    const nlohmann::json j_metadata = model.generate_meta_data_json();

    
}