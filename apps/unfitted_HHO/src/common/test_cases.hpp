
#include "preprocessor.hpp"

template<typename T, typename Function, typename Mesh>
class test_case {
public:
    Function level_set_;
    std::function<T(const typename Mesh::point_type&)> sol_fun;
    std::function<T(const typename Mesh::point_type&)> rhs_fun;
    std::function<T(const typename Mesh::point_type&)> bcs_fun;

    test_case() = default;

    test_case(Function level_set__,
              std::function<T(const typename Mesh::point_type&)> sol_fun_,
              std::function<T(const typename Mesh::point_type&)> rhs_fun_,
              std::function<T(const typename Mesh::point_type&)> bcs_fun_)
        : level_set_(level_set__), sol_fun(sol_fun_), rhs_fun(rhs_fun_), bcs_fun(bcs_fun_)
    {}
};

template<typename Mesh>
auto select_test_case(std::size_t id, const level_set<typename Mesh::coordinate_type>* level_set_function) {
    
    using T = typename Mesh::coordinate_type;
    switch (id) {
        case 0:
            std::cout << "        Test case                  -t : 0 - sin(pi x).sin(pi y)" << std::endl;
            return make_test_case_laplacian_sin_sin<Mesh>(level_set_function);
        case 1:
            std::cout << "        Test case                  -t : 1 - sin(pi x).sin(pi y)" << std::endl;
            return make_test_case_laplacian_sin_sin<Mesh>(level_set_function);
        default:
            throw std::runtime_error("Unknown test case ID: " + std::to_string(id));
    }
}

template<typename Mesh, typename Function>
auto make_test_case_laplacian_sin_sin(Function level_set_function) {

    using T = typename Mesh::coordinate_type;
    using Point = typename Mesh::point_type;
    
    typename preprocessor<Mesh>::template params<T> default_params;
    return test_case<T, Function, Mesh>(
        level_set_function,
        // solution
        [](const Point& pt) -> T {
            return std::sin(M_PI * pt.x()) * std::sin(M_PI * pt.y());
        },
        // rhs
        [](const Point& pt) -> T {
            return 2.0 * M_PI * M_PI * std::sin(M_PI * pt.x()) * std::sin(M_PI * pt.y());
        },
        // boundary conditions
        [](const Point& pt) -> T {
            return std::sin(M_PI * pt.x()) * std::sin(M_PI * pt.y());
        }
    );
}

