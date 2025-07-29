
using namespace disk;

template<typename T, typename Mesh, typename TestCase>
class interface_method {

    using Mat  = Eigen::Matrix<T, Dynamic, Dynamic>;
    using Vect = Eigen::Matrix<T, Dynamic, 1>;

protected:
    interface_method() = default;

    virtual void
    make_contrib_cut(const Mesh& msh, const typename Mesh::cell_type& cl, const TestCase& test_case, const hho_degree_info hdi) {};

public:

    std::pair<Mat, Vect>
    make_contrib_uncut(const Mesh& msh, const typename Mesh::cell_type& cl, const hho_degree_info hdi, const TestCase& test_case) {

        T kappa;
        if (location(msh, cl) == location::IN_NEGATIVE_SIDE)
            kappa = test_case.parms.kappa_1;
        else
            kappa = test_case.parms.kappa_2;

        auto gr = make_hho_gradrec_vector(msh, cl, hdi);
        Mat stab = make_hho_naive_stabilization(msh, cl, hdi);
        Mat lc = kappa * (gr.second + stab);
        Mat f = make_rhs(msh, cl, hdi.cell_degree(), test_case.rhs_fun);
        return std::make_pair(lc, f);
    }
};


template<typename T, typename Mesh, typename TestCase>
class gradrec_interface_method : public interface_method<T, Mesh, TestCase> {

    using Mat  = Matrix<T, Dynamic, Dynamic>;
    using Vect = Matrix<T, Dynamic, 1>;
    using Tuple = std::tuple<double, location, std::vector<double>>;

public:
    TestCase test_case_;  

    gradrec_interface_method(const TestCase& test_case) : test_case_(test_case) {}


};

template<typename T, typename Mesh, typename TestCase>
auto make_gradrec_interface_method(const TestCase& test_case) {
    return gradrec_interface_method<T, Mesh, TestCase>(test_case);
}
