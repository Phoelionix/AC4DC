#include <iostream>
#include <fstream>
#include "src/FreeDistribution.h"
#include "src/Constant.h"


using namespace std;

class BasisTester : public SplineIntegral{
    public:
    BasisTester(size_t F_size, double min_e, double max_e, GridSpacing grid_type) {
        
        Distribution::set_elec_points(F_size, min_e, max_e, grid_type);
        // HACK: There are two distinct BasisSet-inheriting things, the Distribution static BasisIntegral
        // and this object. 
        grid_type.zero_degree=1;
        this->set_parameters(F_size, min_e, max_e, grid_type);
    }

    void check_ee(string fstem)
    {
        
        // Index runs along final states
        double dQ = 0;        
        string Q_loc = fstem + "_Q.csv";

        ofstream qout(Q_loc);

        // Make Q_mat the correct size to store all the summed Q^J_{KL} coefficients 
        std::vector<std::vector<std::vector<double>>> Q_mat;
        
        Q_mat.resize(Distribution::size);
        for (auto&& v : Q_mat) {
            v.resize(Distribution::size);
            for (auto& vi : v) {
                vi.resize(Distribution::size, 0);
            }
        }
        
        // Read all QEE params into file
        SplineIntegral::pair_list ls;
        for (size_t J = 0; J < Distribution::size; J++) {
            for (size_t K=0; K < Distribution::size; K++) {
                ls = calc_Q_ee(J, K);
                for (auto&& entry : ls) {
                    Q_mat[K][entry.idx][J] = entry.val;
                }
            }
        } 

        // Print some helpful labels
        
        qout<<"#L= ";
        for (size_t L=0; L<Distribution::size; L++) {
            qout <<(supp_max(L)+supp_min(L))*0.5<<" ";
        }
        qout<<endl;
        // Print the whole QEE tensor to stdout
        for (size_t J = 0; J < Distribution::size; J++)
        {
            cout<<"J = "<<J<<endl;
            for (size_t K = 0; K < Distribution::size; K++)
            {
                for (size_t L = 0; L < Distribution::size; L++)
                {
                    cout<<Q_mat[K][L][J]<<" ";
                }
                cout<<endl;
            }
        }

        for (size_t K = 0; K < Distribution::size; K++)
        {
            qout <<(supp_max(K)+supp_min(K))*0.5;
            for (size_t L=0; L<Distribution::size; L++) {

                // Aggregate charges in each (K,L) combination, summed over J.
                // This should sum to zero.
                double dQ_tmp=0;
                for (size_t J = 0; J < Distribution::size; J++)
                {
                    dQ_tmp += Q_mat[K][L][J];
                }
                qout<<" "<<dQ_tmp;
                dQ += dQ_tmp;
            }
            qout<<endl;
        }
        qout.close();
        
        cerr<<" total dQ:"<<endl;
        cerr<< dQ <<endl;
    }
};


int main(int argc, char const *argv[]) {
    if (argc < 6) cerr << "usage: "<<argv[0]<<" [filestem] [min_e (Ha)] [max_e (Ha)] [num_basis] [ grid type = [l]inear, [q]uadratic, [e]xponential ]";
    string fstem(argv[1]);
    double min_e = atof(argv[2]);
    double max_e = atof(argv[3]);
    int N = atoi(argv[4]);
    GridSpacing gt;
    istringstream is(argv[5]);
    is >> gt;
    BasisTester bt(N, min_e, max_e, gt);
    bt.check_ee(fstem);
    
    return 0;
}
