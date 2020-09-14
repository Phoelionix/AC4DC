#include <iostream>
#include "src/FreeDistribution.h"
#include "src/Constant.h"
#include <fstream>


using namespace std;

RateData::EIIdata get_fake_eii() {
    RateData::EIIdata tmp;
    tmp.init = 0;
    tmp.push_back(3, 2, 20, 1.9);
    tmp.push_back(5, 2, 44.3, 8.9);
    tmp.push_back(1, 1, 1.3, 0.9);
    tmp.push_back(7, 1, 5, 0.9);
    return tmp;
}

const auto& one = [](double e) -> double {return 1.;};

class BasisTester : public SplineIntegral{
    public:
    BasisTester(size_t F_size, double min_e, double max_e, GridSpacing grid_type) {
        
        Distribution::set_elec_points(F_size, min_e, max_e, grid_type);
        // HACK: There are two distinct BasisSet-inheriting things, the Distribution static BasisIntegral
        // and this object. 
        this->set_parameters(F_size, min_e, max_e, 0, grid_type);
    }

    void check_eii(string filestem)
    {
        
        RateData::EIIdata eii_process = get_fake_eii();
        std::vector<double> tot_gamma;
        int num_fin_states = eii_process.size();
        tot_gamma.resize(num_fin_states, 0);
        std::vector<SparsePair> eg;
        double dQ = 0;
        ofstream out(filestem + "_qeii.csv");

        out <<"#K energy sum_eta gamma_K   sum_J Q_K "<<endl;
        for (size_t K = 0; K < Distribution::size; K++)
        {
            Gamma_eii(eg, eii_process, K);
            double tot = 0;
            for (size_t j = 0; j < num_fin_states; j++)
            {
                tot_gamma[j] += eg[j].val;
                tot += eg[j].val;
            }

            double dQ_tmp = 0;
            for (size_t J = 0; J < Distribution::size; J++)
            {
                dQ_tmp += calc_Q_eii(eii_process, J, K);
            }

            out << K <<" "<<(supp_max(K)+supp_min(K))*0.5<<" "<<tot <<" "<< dQ_tmp <<endl;
            dQ += dQ_tmp;

        }
        out.close();
        cout<<endl<<"init -> fin   gamma_total   "<<endl;
        double tot=0;
        for (size_t j = 0; j < num_fin_states; j++)
        {
            cout<<"[ "<<eii_process.init<<"->"<<eii_process.fin[j]<<" ] "<<tot_gamma[j]<<endl;
            tot += tot_gamma[j];
        }
        cout<<tot<<endl;

        cout<<" total dQ/dP_init:"<<endl;
        cout<< dQ <<endl;
        cout<<" ratio: "<<endl;

        cout<<"Fraction: " << tot/ dQ <<endl;

        
    }
};


int main(int argc, char const *argv[]) {
    if (argc < 6) cerr << "usage: "<<argv[0]<<" [filestem] [min_e (Ha)] [max_e (Ha)] [num_basis] [grid type = [l]inear, [q]uadratic, [e]xponential ]";
    double density = 1;
    string filestem(argv[1]);
    double min_e = atof(argv[2]);
    double max_e = atof(argv[3]);
    int N = atoi(argv[4]);
    GridSpacing gt;
    istringstream inp(argv[5]);
    inp >> gt;
    
    BasisTester bt(N, min_e, max_e, gt);
    bt.check_eii(filestem);
    
    return 0;
}
