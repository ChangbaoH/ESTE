#ifndef GMM_H
#define GMM_H

#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>

const double EPS = 0.00001;  // convergence threshold


// GMM for one dim data with two label
class GaussianMixtureModel {
public:
    int k; // cluster num
    std::vector<double> pi; // mix parameters
    std::vector<double> mu; // means
    std::vector<double> sigma_2; // variance


public:
    // Constructor
    GaussianMixtureModel(std::vector<double>& my_mu, std::vector<double>& my_sigma_2) {
        k = my_mu.size();
        for (int i = 0; i < k; ++i) {
            pi.push_back(1.00/k);
            mu.push_back(my_mu[i]);
            sigma_2.push_back(my_sigma_2[i]);
        }
    }

    double gaussianPDF(double x, double mean, double variance) {
        double coeff = 1.0 / std::sqrt(2 * M_PI * variance);
        double exponent = std::exp(- (std::pow(x - mean, 2)) / (2 * variance));
        return coeff * exponent;
    }


    // E-Step, updata responsibilities
    std::vector<std::vector<double>> expectationStep(std::vector<double>& data) {
        std::vector<std::vector<double>> responsibilities(data.size(), std::vector<double>(k));

        for (int i = 0; i < data.size(); ++i) {
            std::vector<double> pi_temp;
            for (int j = 0; j < k; ++j) {
                pi_temp.push_back(pi[j] * gaussianPDF(data[i], mu[j], sigma_2[j]));
            }

            double sum = std::accumulate(pi_temp.begin(),pi_temp.end(),0.0);
            for (int j = 0; j < k; ++j) {
                responsibilities[i][j] = pi_temp[j] / sum;
            }
        }

        return responsibilities;
    }

    // M-Step, update parameters
    void maximizationStep(std::vector<double>& data, std::vector<std::vector<double>>& responsibilities) {
        std::vector<double> sumGamma(k, 0.0);
        std::vector<double> weightedSum(k, 0.0);
        std::vector<double> varianceSum(k, 0.0);

        // updata means
        for (size_t i = 0; i < data.size(); ++i) {
            for (int j = 0; j < k; ++j) {
                sumGamma[j] += responsibilities[i][j];
                weightedSum[j] += responsibilities[i][j] * data[i];
            }
        }

        for (int i = 0; i < k; ++i) {
            mu[i] = weightedSum[i]/sumGamma[i];
        }

        // update variance, pi
        for (size_t i = 0; i < data.size(); ++i) {
            for (int j = 0; j < k; ++j) {
                varianceSum[j] += responsibilities[i][j] * std::pow(data[i] - mu[j] , 2);
            }
        }

        for (int i = 0; i < k; ++i) {
            sigma_2[i] = varianceSum[i]/sumGamma[i];
            pi[i] = sumGamma[i]/data.size();
        }

    }

    // run EM
    std::vector<std::vector<double>> runEM(std::vector<double>& data, int maxIterations = 10000) {
//        double threshold = 0;
        std::vector<std::vector<double>> cluster_res(k);

        for (int iter = 0; iter < maxIterations; ++iter) {
            // E-Step, updata responsibilities
            std::vector<std::vector<double>> responsibilities = expectationStep(data);

            // save old_value
            std::vector<double> old_mu(mu); // means
            std::vector<double> old_sigma_2(sigma_2); // variance

            // M-Step, update parameters
            maximizationStep(data, responsibilities);

            //Check for convergence
            bool is_convergence = true;
            for (int i = 0; i < k; ++i) {
                if(std::fabs(mu[i] - old_mu[i])/old_mu[i] > EPS | std::fabs(sigma_2[i] - old_sigma_2[i])/old_sigma_2[i] > EPS){
                    is_convergence = false;
                    break;
                }
            }
            if (is_convergence){
                // std::cout << "Converged after " << iter + 1 << " iterations." << std::endl;

                std:vector<std::vector<double>> data_cluster(k);

                for (int i = 0; i < responsibilities.size(); ++i) {
                    int max_Index = -1;
                    double max_Value = -1;
                    for (int j = 0; j < k; ++j) {
                        if(responsibilities[i][j] > max_Value){
                            max_Index = j;
                            max_Value = responsibilities[i][j];
                        }
                    }

                    data_cluster[max_Index].push_back(data[i]);
                }

                for (int i = 0; i < k; ++i) {
                    double min_temp = *std::min_element(data_cluster[i].begin(), data_cluster[i].end());
                    double max_temp = *std::max_element(data_cluster[i].begin(), data_cluster[i].end());
                    cluster_res[i].push_back(min_temp);
                    cluster_res[i].push_back(max_temp);
                }

                break;
            }

        }

        return cluster_res;
    }

    // print current model
    void printParameters() const {
        cout << "pi: ";
        for (int i = 0; i < k; ++i) {
            cout << pi[i] << " ";
        }
        cout << endl << "mu: ";
        for (int i = 0; i < k; ++i) {
            cout << mu[i] << " ";
        }
        cout << endl << "sigma_2: ";
        for (int i = 0; i < k; ++i) {
            cout << sigma_2[i] << " ";
        }
        cout << endl;
    }
};


#endif //GMM_H
