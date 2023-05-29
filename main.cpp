#include <iostream>
#include<vector>
#include <cmath>
#include <string>
#include <fstream>
#include <cassert>

template <typename T>
class StormletVerlet{

    StormletVerlet(int p_dim, std::vector<T> p_q0, std::vector<T> p_p0, std::vector<T> p_m, int p_iter_number = 1000, T p_step=1e-3) :
    dim(p_dim), q0(p_q0), p0(p_p0), m(p_m), iter_number(p_iter_number), step(p_step)
    {
         assert(p_q0.size() == p_p0.size());
         assert(p_p0.size() % p_m.size() == 0);
    }

    private:
        int dim;  //dimension of configuration space
        std::vector<T> q0;  //initial coordinates
        std::vector<T> p0;  //initial momentum vector
        std::vector<T> m; // vector of masses of bodies
        int iter_number = 1000;  // amount of iterations for integrator
        double step = 1e-3;

        static T get_norm(std::vector<T>& x, size_t  norm_type=2){
            T norm = 0.0;
            for (const auto &x_i : x){
                norm += pow(x_i, norm_type);
            }
            return norm;
        }

        T get_U(std::vector<T>& q){
            /* Potential Energy
             * q - vector in configuration space
             * m - mass vector
             * */
            T u = 0.0;
            std::vector<T> q_tmp;
            for (size_t i = 0; i != q.size() - 1; ++i){
                for (size_t j = i + 1; j != q.size(); ++j){
                    q_tmp = q[i] - q[j];
                    u += m[i] * m[j] / get_norm(q_tmp);
                }
            }
            return u;
        }

        T get_T(std::vector<T>& p){
            /* Kinetic Energy
            * p - momenta vector
            * m - mass vector
            * */
            T t_kin = 0.0;
            for(size_t i = 0; i != p.size(); ++i){
                t_kin += get_norm(p[i]) / (2.0 * m[i]);
            }
            return t_kin;
        }

        std::vector<T> get_dT(std::vector<T>& p){
            /*Derivative of kinetic energy by momentum
             * p - momemtum vector
             * m - mass vector
             * return: Vector (dT / dp_i), where i = 1, .., dim
             * */
            std::vector<T> dT(dim);
            for (size_t i; i != p.size(); ++p){
                dT[i] = p[i] / m[i];
            }
            return dT;
        }

        std::vector<T> get_dU(std::vector<T>& q){
            /*Derivative of potential energy by coords
            * q - coords vector
            * m - mass vector
            * return: Vector (dU / dq_i), where i = 1, .., dim
            * */
            std::vector<T> dU(dim);
            for(size_t i = 0; i != q.size() - 1; ++i){
                for (size_t j = 0; j != q.size(); ++j){

                }
            }
            return dU;
        }

        T ** integrate(bool to_file = false, std::string const & filename="result.dat"){
            std::vector<T> p = p0;
            std::vector<T> q = q0;
            std::vector<T> p_half2;
            T ** result = new T *[iter_number];
            for (size_t i = 0; i != iter_number; ++i) {
                result[i] = new T[dim + 1];
            }
            for (int i = 0; i < iter_number; ++i){
                p_half2 = p - step / 2.0 * get_dU(q, m);
                q += step * get_dT(p_half2);
                p =  p_half2 - step / 2.0 * get_dU(q);
                for (int j = 0; j < q.size(); ++j){result[i][j] = q[j];}
                for (int j = q.size(); j < dim; ++j){result[i][j] = p[j - q.size()];}
                result[i][dim] = i * step;
            }
            return result;
        }



};



int main() {

    return 0;
}
