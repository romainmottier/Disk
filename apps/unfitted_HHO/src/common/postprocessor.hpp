#pragma once
#ifndef postprocessor_hpp
#define postprocessor_hpp

#include <iomanip>

template<typename Mesh>
class postprocessor {
    
public:

    using Tuple = std::tuple<double, disk::location, std::vector<double>>;
    using VecTuple = std::vector<std::tuple<double, disk::location, std::vector<double>>>;

//     // PICK CELLS & FIND CELLS   
//     /////////////////////////////////////////////////////////////////////////////
//     /////////////////////////////////////////////////////////////////////////////
//     /////////////////////////////////////////////////////////////////////////////

//     /// Find the cells associated to the requested point
//     static std::set<size_t> 
//     find_cells(typename Mesh::point_type & pt, Mesh & msh, bool verbose_Q = false){
        
//         using RealType = double;
//         auto norm =  [](const typename Mesh::point_type& a, const typename Mesh::point_type& b ) -> RealType {
//             RealType dx = (b.x() - a.x());
//             RealType dy = (b.y() - a.y());
//             RealType norm = std::sqrt(dx*dx + dy*dy);
//             return norm;
//         };
        
//         // find minimum distance to the requested point
//         size_t np = msh.points.size();
//         std::vector<RealType> distances(np);
        
//         size_t ip = 0;
//         for (auto& point : msh.points)
//         {
//             RealType dist = norm(pt,point);
//             distances[ip] = dist;
//             ip++;
//         }
        
//         size_t index = std::min_element(distances.begin(),distances.end()) - distances.begin();
//         if(verbose_Q){
//             RealType min_dist = *std::min_element(distances.begin(), distances.end());
//             typename Mesh::point_type nearest_point = msh.points.at(index);
//             std::cout << "Nearest point detected : " << std::endl;
//             std::cout << "  x =  " << nearest_point.x() << std::endl;
//             std::cout << "  y =  " << nearest_point.y() << std::endl;
//             std::cout << "Distance = " << min_dist << std::endl;
//             std::cout << "Global index = " << index << std::endl;
//         }
        
//         std::set<size_t> cell_indexes;
//         size_t cell_i = 0;
//         for (auto& cell : msh.cells)
//         {
//             auto cell_nodes = nodes(msh,cell);
//             size_t n_p = cell_nodes.size();
//             for (size_t l = 0; l < n_p; l++)
//             {
//                 auto node = cell_nodes[l];
//                 if(index == node.ptid){
//                     cell_indexes.insert(cell_i);
//                 }
//             }
//             cell_i++;
//         }
        
//         if(verbose_Q){
//             std::cout << "Detected cells indexes : " << std::endl;
//             for(auto index : cell_indexes){
//                 std::cout << index << std::endl;
//             }
//         }

//         return cell_indexes;
//     }
    
//     /// Pick the cell that contains the requested point
//     static size_t 
//     pick_cell(typename Mesh::point_type & pt, Mesh & msh, std::set<size_t> & cell_indexes, bool verbose_Q = false){
        
//         using RealType = double;
        
//         auto triangle_member_Q = [] (typename Mesh::point_type & p, typename Mesh::point_type & p0, typename Mesh::point_type & p1, typename Mesh::point_type & p2)
//         {
//             RealType dx = p.x()-p2.x();
//             RealType dy = p.y()-p2.y();
//             RealType dx21 = p2.x()-p1.x();
//             RealType dy12 = p1.y()-p2.y();
//             RealType d = dy12*(p0.x()-p2.x()) + dx21*(p0.y()-p2.y());
//             RealType s = dy12*dx + dx21*dy;
//             RealType t = (p2.y()-p0.y())*dx + (p0.x()-p2.x())*dy;
//             if (d < 0.0) {
//                 return s<=0.0 && t<=0.0 && s+t>=d;
//             }
//             return s>=0 && t>=0 && s+t<=d;
//         };
        
//         size_t n_cells = cell_indexes.size();
//         if (n_cells == 1) {
//             size_t first_index = *cell_indexes.begin();
//             return first_index;
//         }
//         bool is_member_Q = false;
//         for(auto index : cell_indexes){
//             auto& cell = msh.cells.at(index);
//             auto bar = barycenter(msh, cell);
//             auto cell_nodes = nodes(msh,cell);
//             size_t n_p = cell_nodes.size();
            
//             // building teselation
//             std::vector<std::vector<typename Mesh::point_type>> triangles(n_p);
//             for (size_t l = 0; l < n_p; l++)
//             {

//                 std::vector<typename Mesh::point_type> chunk(3);
//                 if( l == n_p - 1){
//                     chunk[0] = msh.points.at(cell_nodes[l].ptid);
//                     chunk[1] = msh.points.at(cell_nodes[0].ptid);
//                     chunk[2] = bar;
//                 }else{
//                     chunk[0] = msh.points.at(cell_nodes[l].ptid);
//                     chunk[1] = msh.points.at(cell_nodes[l+1].ptid);
//                     chunk[2] = bar;
//                 }
//                 triangles[l] = chunk;
//             }
            
//             // check whether the point is memeber of any triangle
//             for (auto triangle : triangles) {
//                 is_member_Q = triangle_member_Q(pt,triangle[0],triangle[1],triangle[2]);
//                 if (is_member_Q) {
//                     std::cout << "Detected cell index = " << index << std::endl;
//                     return index;
//                 }
//             }

//         }
        
//         if(!is_member_Q){
//             if(verbose_Q){
//                 std::cout << "Point is not member of cells set. Returning cell_indexes[0] " << std::endl;
//             }
//             size_t first_index = *cell_indexes.begin();
//             return first_index;
//         }
        
//         return -1;
//     }

//     // SILO CONDITIONING 
//     /////////////////////////////////////////////////////////////////////////////
//     /////////////////////////////////////////////////////////////////////////////
//     /////////////////////////////////////////////////////////////////////////////
    
//     static void 
//     write_silo_conditioning(std::string silo_file_name, Mesh & msh, hho_degree_info & hho_di, Matrix<double, Dynamic, 1> & conditioning, interface_assembler<Mesh, std::function<double(const typename Mesh::point_type& )>> & assembler) {

//         using RealType = double;
//         timecounter tc;
//         tc.tic();

//         auto dim = 2;
//         std::vector<RealType> cond;

//         auto cell_table = assembler.get_cell_table();
//         for (auto& cl : msh.cells) {
//             auto offset_cl = cell_table.at(offset(msh, cl));
//             if (!is_cut(msh, cl)) 
//                 cond.push_back(conditioning(offset_cl));
//             else 
//                 cond.push_back(std::max(conditioning(offset_cl), conditioning(offset_cl+1)));
//         }
        
//         silo_database silo;
//         silo_file_name += ".silo";
//         silo.create(silo_file_name.c_str());
//         silo.add_mesh(msh, "mesh");
//         silo.add_variable("mesh", "cond",  cond.data(),  cond.size(), zonal_variable_t);
//         silo.close();
//         tc.toc();
//         std::cout << bold << yellow << "         Silo file - CONDITIONING - rendered in : " << tc << " seconds" << reset << std::endl;
//     }

//     // ERRORS
//     /////////////////////////////////////////////////////////////////////////////
//     /////////////////////////////////////////////////////////////////////////////
//     /////////////////////////////////////////////////////////////////////////////
    
//     static std::vector<double> 
//     compute_error_elliptic_second_order_agglo(Mesh & msh, hho_degree_info & hho_di, interface_assembler<Mesh, std::function<double(const typename Mesh::point_type& )>> & assembler, Matrix<double, Dynamic, 1> & x_dof,std::function<double(const typename Mesh::point_type& )> sol_fun, std::function<Matrix<double, 1, 2>(const typename Mesh::point_type& )> sol_grad, double previous_h, double previous_L2, double previous_H1, std::ostream & error_file = std::cout){

//        timecounter tc;
//        tc.tic();

//        using RealType = double;

//        RealType H1_error = 0.0;
//        RealType L2_error = 0.0;
//        size_t   cell_i   = 0;
//        RealType h = 10;
//        for (auto& cl : msh.cells) {
       
//            // Diameter
//            RealType h_l = diameter(msh, cl);
//            if (h_l < h) 
//                h = h_l;
               
//             // UNCUT CELLS 
//             if (!is_cut(msh, cl)) {
//                 // Bases & dofs infos
//                 cell_basis<cuthho_poly_mesh<RealType>, RealType> cb(msh, cl, hho_di.cell_degree());
//                 auto cbs = cb.size();
//                 auto fcs = faces(msh, cl);
//                 auto num_faces = fcs.size();
//                 auto fbs = face_basis<cuthho_poly_mesh<RealType>,RealType>::size(hho_di.face_degree());
//                 Matrix<RealType, Dynamic, 1> locdata_n, locdata_p, locdata;
//                 Matrix<RealType, Dynamic, 1> cell_dofs_n, cell_dofs_p, cell_dofs;
//                 // Compute errors
//                 locdata = assembler.take_local_data(msh, cl, x_dof, location::IN_NEGATIVE_SIDE);
//                 cell_dofs = locdata.head(cbs);
//                 auto qps = integrate(msh, cl, 2*hho_di.cell_degree());
//                 for (auto& qp : qps) {
//                     /* Compute H1-error */
//                     auto t_dphi = cb.eval_gradients( qp.first );
//                     Matrix<RealType, 1, 2> grad = Matrix<RealType, 1, 2>::Zero();
//                     for (size_t i = 1; i < cbs; i++ )
//                         grad += cell_dofs(i) * t_dphi.block(i, 0, 1, 2);
//                     H1_error += qp.second * (sol_grad(qp.first) - grad).dot(sol_grad(qp.first) - grad);
//                     auto t_phi = cb.eval_basis( qp.first );
//                     auto v = cell_dofs.dot(t_phi);
//                     /* Compute L2-error */
//                     L2_error += qp.second * (sol_fun(qp.first) - v) * (sol_fun(qp.first) - v);
//                 }
//             }
//             // CUT CELLS
//             else {
//                 // DISCRETIZATION INFOS NEGATIVE SIDE
//                 #ifndef subcell_centering
//                 cell_basis<cuthho_poly_mesh<RealType>, RealType> neg_cell_basis(msh, cl, hho_di.cell_degree());
//                 #else
//                 cut_cell_basis<cuthho_poly_mesh<RealType>, RealType> neg_cell_basis(msh, cl, hho_di.cell_degree(), location::IN_NEGATIVE_SIDE);
//                 #endif
//                 auto cbs = neg_cell_basis.size();
//                 auto locdata_n = assembler.take_local_data(msh, cl, x_dof, location::IN_NEGATIVE_SIDE);
//                 auto cell_dofs_n = locdata_n.head(cbs);
//                 auto qps_n = integrate(msh, cl, 2*hho_di.cell_degree(), location::IN_NEGATIVE_SIDE);
//                 for (auto& qp : qps_n) {
//                     /* Compute H1-error */
//                     auto t_dphi = neg_cell_basis.eval_gradients( qp.first );
//                     Matrix<RealType, 1, 2> grad = Matrix<RealType, 1, 2>::Zero();
//                     for (size_t i = 1; i < cbs; i++ )
//                         grad += cell_dofs_n(i) * t_dphi.block(i, 0, 1, 2);
//                     H1_error += qp.second * (sol_grad(qp.first) - grad).dot(sol_grad(qp.first) - grad);
//                     auto t_phi = neg_cell_basis.eval_basis( qp.first );
//                     auto v = cell_dofs_n.dot(t_phi);
//                     /* Compute L2-error */
//                     L2_error += qp.second * (sol_fun(qp.first) - v) * (sol_fun(qp.first) - v);
//                 }  
//                 // DISCRETIZATION INFOS POSITIVE SIDE
//                 #ifndef subcell_centering
//                 cell_basis<cuthho_poly_mesh<RealType>, RealType> pos_cell_basis(msh, cl, hho_di.cell_degree());
//                 #else
//                 cut_cell_basis<cuthho_poly_mesh<RealType>, RealType> pos_cell_basis(msh, cl, hho_di.cell_degree(), location::IN_POSITIVE_SIDE);           
//                 #endif
//                 cbs = pos_cell_basis.size();
//                 auto locdata_p = assembler.take_local_data(msh, cl, x_dof, location::IN_POSITIVE_SIDE);
//                 auto cell_dofs_p = locdata_p.head(cbs);
//                 auto qps_p = integrate(msh, cl, 2*hho_di.cell_degree(), location::IN_POSITIVE_SIDE);
//                 for (auto& qp : qps_p) {
//                     /* Compute H1-error */
//                     auto t_dphi = pos_cell_basis.eval_gradients( qp.first );
//                     Matrix<RealType, 1, 2> grad = Matrix<RealType, 1, 2>::Zero();
//                     for (size_t i = 1; i < cbs; i++ )
//                         grad += cell_dofs_p(i) * t_dphi.block(i, 0, 1, 2);
//                     H1_error += qp.second * (sol_grad(qp.first) - grad).dot(sol_grad(qp.first) - grad);
//                     auto t_phi = pos_cell_basis.eval_basis( qp.first );
//                     auto v = cell_dofs_p.dot(t_phi);
//                     /* Compute L2-error */
//                     L2_error += qp.second * (sol_fun(qp.first) - v) * (sol_fun(qp.first) - v);
//                 }
//             }
//        }
//        H1_error = std::sqrt(H1_error);
//        L2_error = std::sqrt(L2_error);
//        RealType orderH = log(previous_H1 / H1_error) / log(previous_h / h);
//        RealType orderL = log(previous_L2 / L2_error) / log(previous_h / h);
//        error_file << "Characteristic h size = " << h << std::endl;
//        error_file << "L2-norm error = " << L2_error << std::endl;
//        error_file << "H1-norm error = " << H1_error << std::endl;
//        error_file << "order L2 = " << orderL << std::endl;
//        error_file << "order H1 = " << orderH << std::endl << std::endl;
//        std::vector<RealType> vec = {h, H1_error, L2_error};
//        tc.toc();

//        std::cout << bold << yellow << "         H1-Error: " << H1_error << reset << std::endl;
//        std::cout << bold << yellow << "         L2-Error: " << L2_error << reset << std::endl;
//        std::cout << bold << yellow << "         order H1 = " << orderH << std::endl;
//        std::cout << bold << yellow << "         order L2 = " << orderL << std::endl;
//        std::cout << bold << yellow << "         Error completed: " << tc << " seconds" << reset << std::endl;

//        return vec;

//     }

//     template<typename testType>
//     static std::vector<double> 
//     compute_error_elliptic_second_order_poly_ext(Mesh & msh, VecTuple POK, hho_degree_info & hho_di, interface_assembler<Mesh, std::function<double(const typename Mesh::point_type& )>> & assembler, Matrix<double, Dynamic, 1> & x_dof, const testType &test_case, double previous_h, double previous_L2, double previous_H1, std::ostream & error_file = std::cout) {

//         timecounter tc;
//         tc.tic();
        
//         using RealType = double;
        
//         // TEST CASE
//         auto sol_fun = test_case.sol_fun;
//         auto sol_grad = test_case.sol_grad;
        
//         RealType H1_error = 0.0;
//         RealType L2_error = 0.0;
//         RealType h = 10;
        
//         for (auto& p_ok : POK) {
            
//             // CELL INFOS 
//             auto cell_index = std::get<0>(p_ok);
//             auto cl = msh.cells[cell_index];
//             auto loc = std::get<1>(p_ok);
//             double kappa;
//             if (loc == location::IN_NEGATIVE_SIDE)
//                 kappa = test_case.parms.kappa_1;
//             else 
//                 kappa = test_case.parms.kappa_2;

//             // DIAMETER
//             RealType h_l = diameter(msh, cl);
//             if (h_l < h) 
//                 h = h_l;

//             // BASES & DOFS INFOS  
//             #ifndef subcell_centering
//             cell_basis<cuthho_poly_mesh<RealType>, RealType> cb(msh, cl, hho_di.cell_degree());
//             #else
//             cut_cell_basis<cuthho_poly_mesh<RealType>, RealType> cb(msh, cl, hho_di.cell_degree(), loc);
//             #endif
//             auto cbs = cb.size();
//             Matrix<RealType, Dynamic, 1> locdata_n, locdata_p, locdata;
//             Matrix<RealType, Dynamic, 1> cell_dofs_n, cell_dofs_p, cell_dofs;
            
//             // COMPUTE ERROR OF ONE SIDE
//             locdata = assembler.take_local_data(msh, cl, x_dof, loc);
//             cell_dofs = locdata.head(cbs);
            
//             // UNCUT CELLS 
//             if (!is_cut(msh, cl)) {
//                 auto qps = integrate(msh, cl, 2*hho_di.cell_degree());
//                 for (auto& qp : qps) {
//                     /* Compute H1-error */
//                     auto t_dphi = cb.eval_gradients( qp.first );
//                     Matrix<RealType, 1, 2> grad = Matrix<RealType, 1, 2>::Zero();
//                     for (size_t i = 1; i < cbs; i++ ) 
//                         grad += cell_dofs(i) * t_dphi.block(i, 0, 1, 2);
//                     H1_error += kappa * qp.second * (sol_grad(qp.first) - grad).dot(sol_grad(qp.first) - grad);
//                     /* Compute L2-error */
//                     auto t_phi = cb.eval_basis( qp.first );
//                     auto v = cell_dofs.dot(t_phi);
//                     L2_error += kappa * qp.second * (sol_fun(qp.first) - v) * (sol_fun(qp.first) - v);
//                 }
//             }
//             // CUT CELLS
//             else {
//                 auto qps = integrate(msh, cl, 2*hho_di.cell_degree(), loc);
//                 for (auto& qp : qps) {
//                     /* Compute H1-error */
//                     auto t_dphi = cb.eval_gradients( qp.first );
//                     Matrix<RealType, 1, 2> grad = Matrix<RealType, 1, 2>::Zero();
//                     for (size_t i = 1; i < cbs; i++ ) 
//                         grad += cell_dofs(i) * t_dphi.block(i, 0, 1, 2);
//                     H1_error += kappa * qp.second * (sol_grad(qp.first) - grad).dot(sol_grad(qp.first) - grad);
//                     /* Compute L2-error */
//                     auto t_phi = cb.eval_basis( qp.first );
//                     auto v = cell_dofs.dot(t_phi);
//                     L2_error += kappa * qp.second * (sol_fun(qp.first) - v) * (sol_fun(qp.first) - v);
//                 }
//             }
//             // DEPENDENT CELLS 
//             for (auto& dp_cl : std::get<2>(p_ok)) {
//                 // CELL INFOS 
//                 auto dp_cell = msh.cells[dp_cl];
//                 auto locdata = assembler.take_local_data(msh, cl, x_dof, loc);
//                 auto cell_dofs = locdata.head(cbs);
//                 // COMPUTE ERRORS
//                 auto qps = integrate(msh, dp_cell, 2*hho_di.cell_degree(), loc);
//                 for (auto& qp : qps) {
//                     /* Compute H1-error */
//                     auto t_dphi = cb.eval_gradients( qp.first );
//                     Matrix<RealType, 1, 2> grad = Matrix<RealType, 1, 2>::Zero();
//                     for (size_t i = 1; i < cbs; i++ )
//                         grad += cell_dofs(i) * t_dphi.block(i, 0, 1, 2);
//                     H1_error += kappa * qp.second * (sol_grad(qp.first) - grad).dot(sol_grad(qp.first) - grad);
//                     /* Compute L2-error */
//                     auto t_phi = cb.eval_basis( qp.first );
//                     auto v = cell_dofs.dot(t_phi);
//                     L2_error += kappa * qp.second * (sol_fun(qp.first) - v) * (sol_fun(qp.first) - v);
//                 }
//             }
//         }
//         H1_error = std::sqrt(H1_error);
//         L2_error = std::sqrt(L2_error);
//         RealType orderH = log(previous_H1 / H1_error) / log(previous_h / h);
//         RealType orderL = log(previous_L2 / L2_error) / log(previous_h / h);
//         error_file << "Characteristic h size = " << h << std::endl;
//         error_file << "L2-norm error = " << L2_error << std::endl;
//         error_file << "H1-norm error = " << H1_error << std::endl;
//         error_file << "order L2 = " << orderL << std::endl;
//         error_file << "order H1 = " << orderH << std::endl << std::endl;
//         std::vector<RealType> vec = {h, H1_error, L2_error};
//         tc.toc();

//         std::cout << bold << yellow << "         H1-Error: " << H1_error << reset << std::endl;
//         std::cout << bold << yellow << "         L2-Error: " << L2_error << reset << std::endl;
//         std::cout << bold << yellow << "         order H1: " << orderH << reset << std::endl;
//         std::cout << bold << yellow << "         order L2: " << orderL << reset << std::endl;
//         std::cout << bold << yellow << "         Error completed: " << tc << " seconds" << reset << std::endl;

//        return vec;

//     }

//     #ifndef subcell_centering
//     static void 
//     compute_errors_grad_one_field(Mesh & msh, hho_degree_info & hho_di, interface_assembler<Mesh, std::function<double(const typename Mesh::point_type& )>> & assembler, Matrix<double, Dynamic, 1> & grad_dof, std::function<Matrix<double, 1, 2>(const typename Mesh::point_type& )> flux_fun, std::ostream & error_file = std::cout){

//         timecounter tc;
//         tc.tic();

//         using RealType = double;
        
//         RealType grad_l2_error = 0.0;
//         size_t cell_i = 0;
//         RealType h = 10.0;
//         std::vector<RealType> l2_error_grad(msh.cells.size());
        
//         auto cell_table = assembler.get_cell_table();
//         for (auto& cl : msh.cells ) {
            
//             // CELL INFOS
//             l2_error_grad[cell_i] = 0.0;
//             RealType h_l = diameter(msh, cl);
//             if (h_l < h)
//                 h = h_l;

//             vector_cell_basis<cuthho_poly_mesh<RealType>, RealType> vector_cell_basis(msh, cl, hho_di.grad_degree());
//             auto grad_dofs_size = vector_cell_basis.size();
//             auto offset_cl = cell_table.at(offset(msh, cl));
           
//             if (!is_cut(msh, cl)) {
//                 auto cell_GRAD_offset = offset_cl * grad_dofs_size;
//                 Eigen::VectorXd grad_dofs = grad_dof.block(cell_GRAD_offset, 0, grad_dofs_size, 1);
//                 auto qps = integrate(msh, cl, 2*hho_di.grad_degree());
//                 for (auto& qp : qps) {
//                     auto vec_t_phi = vector_cell_basis.eval_basis(qp.first);
//                     auto grad = grad_dofs.transpose() * vec_t_phi;
//                     l2_error_grad[cell_i] += qp.second * (flux_fun(qp.first) - grad).dot(flux_fun(qp.first) - grad);
//                 }
//             }
//             else {
//                 {   // NEGATIVE SIDE
//                     auto cell_GRAD_offset = offset_cl * grad_dofs_size;
//                     Eigen::VectorXd grad_dofs = grad_dof.block(cell_GRAD_offset, 0, grad_dofs_size, 1);
//                     auto qps = integrate(msh, cl, 2*hho_di.grad_degree(), location::IN_NEGATIVE_SIDE);
//                     for (auto& qp : qps) {
//                         auto vec_t_phi = vector_cell_basis.eval_basis(qp.first);
//                         auto grad = grad_dofs.transpose() * vec_t_phi;
//                         l2_error_grad[cell_i] += qp.second * (flux_fun(qp.first) - grad).dot(flux_fun(qp.first) - grad);
//                     }
//                 }
//                 {   // POSITIVE SIDE
//                     auto cell_GRAD_offset = offset_cl*grad_dofs_size + grad_dofs_size;
//                     Eigen::VectorXd grad_dofs = grad_dof.block(cell_GRAD_offset, 0, grad_dofs_size, 1);
//                     auto qps = integrate(msh, cl, 2*hho_di.grad_degree(), location::IN_POSITIVE_SIDE);
//                     for (auto& qp : qps) {
//                         auto vec_t_phi = vector_cell_basis.eval_basis(qp.first);
//                         auto grad = grad_dofs.transpose() * vec_t_phi;
//                         l2_error_grad[cell_i] += qp.second * (flux_fun(qp.first) - grad).dot(flux_fun(qp.first) - grad);
//                     }
//                 }
//             }
//             cell_i++;
//         }
        
//         grad_l2_error = std::accumulate(l2_error_grad.begin(), l2_error_grad.end(), 0.0);
//         tc.toc();
       
//         std::cout << bold << yellow << "         Test Gradient completed" << reset << std::endl;
//         error_file << "Characteristic h size = " << std::setprecision(16) << h << std::endl;
//         error_file << "L2-norm grad error = " << std::setprecision(16) << std::sqrt(grad_l2_error) << std::endl;
//         error_file << std::endl;
//         error_file.flush();
       
//     }
//     #else
//     static void 
//     compute_errors_grad_one_field(Mesh & msh, hho_degree_info & hho_di, interface_assembler<Mesh, std::function<double(const typename Mesh::point_type& )>> & assembler, Matrix<double, Dynamic, 1> & grad_dof, std::function<Matrix<double, 1, 2>(const typename Mesh::point_type& )> flux_fun, std::ostream & error_file = std::cout){

//         timecounter tc;
//         tc.tic();

//         using RealType = double;
        
//         RealType grad_l2_error = 0.0;
//         size_t cell_i = 0;
//         RealType h = 10.0;
//         std::vector<RealType> l2_error_grad(msh.cells.size());
                    
//         auto cell_table = assembler.get_cell_table();
//         for (auto& cl : msh.cells ) {

//             // CELL INFOS
//             RealType h_l = diameter(msh, cl);
//             if (h_l < h) 
//                 h = h_l;

//             l2_error_grad[cell_i] = 0.0;
//             auto offset_cl = cell_table.at(offset(msh, cl));
//             if (!is_cut(msh, cl)) {
//                 cut_vector_cell_basis<cuthho_poly_mesh<RealType>, RealType> cut_vector_cell_basis(msh, cl, hho_di.grad_degree(), locate(msh, cl));
//                 auto grad_dofs_size = cut_vector_cell_basis.size();
//                 auto cell_GRAD_offset = offset_cl * grad_dofs_size;
//                 Eigen::VectorXd grad_dofs = grad_dof.block(cell_GRAD_offset, 0, grad_dofs_size, 1);
//                 auto qps = integrate(msh, cl, 2*hho_di.grad_degree());
//                 for (auto& qp : qps) {
//                     // auto vec_t_phi = vector_cell_basis.eval_basis(qp.first);
//                     auto vec_t_phi = cut_vector_cell_basis.eval_basis(qp.first);
//                     auto grad = grad_dofs.transpose() * vec_t_phi;
//                     l2_error_grad[cell_i] += qp.second * (flux_fun(qp.first) - grad).dot(flux_fun(qp.first) - grad);
//                 }
//             }
//             else {
//                 {   // NEGATIVE SIDE                 
//                     cut_vector_cell_basis<cuthho_poly_mesh<RealType>, RealType> cut_vector_cell_basis(msh, cl, hho_di.grad_degree(), location::IN_NEGATIVE_SIDE);
//                     auto grad_dofs_size = cut_vector_cell_basis.size();
//                     auto cell_GRAD_offset = offset_cl * grad_dofs_size;
//                     Eigen::VectorXd grad_dofs = grad_dof.block(cell_GRAD_offset, 0, grad_dofs_size, 1);
//                     auto qps = integrate(msh, cl, 2*hho_di.grad_degree(), location::IN_NEGATIVE_SIDE);
//                     for (auto& qp : qps) {
//                         // auto vec_t_phi = vector_cell_basis.eval_basis(qp.first);
//                         auto vec_t_phi = cut_vector_cell_basis.eval_basis(qp.first);
//                         auto grad = grad_dofs.transpose() * vec_t_phi;
//                         l2_error_grad[cell_i] += qp.second * (flux_fun(qp.first) - grad).dot(flux_fun(qp.first) - grad);
//                     }
//                 }
//                 {   // POSITIVE SIDE
//                     cut_vector_cell_basis<cuthho_poly_mesh<RealType>, RealType> cut_vector_cell_basis(msh, cl, hho_di.grad_degree(), location::IN_POSITIVE_SIDE);
//                     auto grad_dofs_size = cut_vector_cell_basis.size();
//                     auto cell_GRAD_offset = offset_cl*grad_dofs_size + grad_dofs_size;
//                     Eigen::VectorXd grad_dofs = grad_dof.block(cell_GRAD_offset, 0, grad_dofs_size, 1);
//                     auto qps = integrate(msh, cl, 2*hho_di.grad_degree(), location::IN_POSITIVE_SIDE);
//                     for (auto& qp : qps) {
//                         // auto vec_t_phi = vector_cell_basis.eval_basis(qp.first);
//                         auto vec_t_phi = cut_vector_cell_basis.eval_basis(qp.first);
//                         auto grad = grad_dofs.transpose() * vec_t_phi;
//                         l2_error_grad[cell_i] += qp.second * (flux_fun(qp.first) - grad).dot(flux_fun(qp.first) - grad);
//                     }
//                 }
//             }
//             cell_i++;
//         }
        
//         grad_l2_error = std::accumulate(l2_error_grad.begin(), l2_error_grad.end(), 0.0);
//         tc.toc();
       
//         std::cout << bold << yellow << "         Gradient error completed: " << tc << " seconds" << reset << std::endl;
//         error_file << "Characteristic h size = " << std::setprecision(16) << h << std::endl;
//         error_file << "L2-norm grad error = " << std::setprecision(16) << std::sqrt(grad_l2_error) << std::endl;
//         error_file << std::endl;
//         error_file.flush();
       
//     }
//     #endif 

//     /// Compute L2 and H1 errors for one field approximation
//     static void compute_errors_grad_grad(Mesh & msh, hho_degree_info & hho_di, double grad_grad_dofs, std::ostream & error_file = std::cout){

//         using RealType = double;

//         RealType h = 10.0;
//         for (auto& cell : msh.cells) {
//             RealType h_l = diameter(msh, cell);
//             if (h_l < h) 
//                 h = h_l;
//         }

//         auto error = std::abs(grad_grad_dofs - M_PI*M_PI/2.0);
//         error_file << "Characteristic h size = " << std::setprecision(16) << h << std::endl;
//         error_file << "L2-norm error = " << std::setprecision(16) << error << std::endl;
//         error_file << std::endl;
       
//     }
    
//     /// Compute L2 and H1 errors for one field approximation
//     static void 
//     compute_errors_one_field_bis(Mesh & msh, hho_degree_info & hho_di, interface_assembler<Mesh, std::function<double(const typename Mesh::point_type& )>> & assembler, Matrix<double, Dynamic, 1> & x_dof,std::function<double(const typename Mesh::point_type& )> scal_fun, std::function<Matrix<double, 1, 2>(const typename Mesh::point_type& )> flux_fun, std::ostream & error_file = std::cout){

//        timecounter tc;
//        tc.tic();

//        using RealType = double;

//        RealType scalar_l2_error = 0.0;
//        RealType flux_l2_error = 0.0;
//        size_t cell_i = 0;
//        RealType h = 10.0;
//        std::vector<RealType> l2_error_vec(msh.cells.size());
//        std::vector<RealType> flux_l2_error_vec(msh.cells.size());
//        for (auto& cell : msh.cells) {

//            l2_error_vec[cell_i] = 0.0;
//             RealType h_l = diameter(msh, cell);
//            if (h_l < h) 
//                h = h_l;
        
           
//            cell_basis<cuthho_poly_mesh<RealType>, RealType> cell_basis(msh, cell, hho_di.cell_degree());
//            auto cbs = cell_basis.size();
//            if ( locate(msh, cell) == location::ON_INTERFACE )
//            {
               
//                auto dofs_n = assembler.take_local_data(msh, cell, x_dof, location::IN_NEGATIVE_SIDE);
//                auto dofs_p = assembler.take_local_data(msh, cell, x_dof, location::IN_POSITIVE_SIDE);

//                auto cell_dofs_n = dofs_n.head(cbs);
//                auto cell_dofs_p = dofs_p.head(cbs);
               
//                // negative side
//                auto qps_n = integrate(msh, cell, 2*hho_di.cell_degree(), location::IN_NEGATIVE_SIDE);
//                for (auto& qp : qps_n)
//                {
//                    /* Compute H1-error */
//                    auto t_dphi = cell_basis.eval_gradients( qp.first );
//                    Matrix<RealType, 1, 2> grad = Matrix<RealType, 1, 2>::Zero();

//                    for (size_t i = 1; i < cbs; i++ )
//                        grad += cell_dofs_n(i) * t_dphi.block(i, 0, 1, 2);

//                    flux_l2_error_vec[cell_i] += qp.second * (flux_fun(qp.first) - grad).dot(flux_fun(qp.first) - grad);
                   

//                    auto t_phi = cell_basis.eval_basis( qp.first );
//                    auto v = cell_dofs_n.dot(t_phi);
                   
//                    /* Compute L2-error */
//                    l2_error_vec[cell_i] += qp.second * (scal_fun(qp.first) - v) * (scal_fun(qp.first) - v);
//                }
               
//                // positive side
//                auto qps_p = integrate(msh, cell, 2*hho_di.cell_degree(), location::IN_POSITIVE_SIDE);
//                for (auto& qp : qps_p)
//                {
//                    /* Compute H1-error */
//                    auto t_dphi = cell_basis.eval_gradients( qp.first );
//                    Matrix<RealType, 1, 2> grad = Matrix<RealType, 1, 2>::Zero();

//                    for (size_t i = 1; i < cbs; i++ )
//                        grad += cell_dofs_n(i) * t_dphi.block(i, 0, 1, 2);

//                    flux_l2_error_vec[cell_i] += qp.second * (flux_fun(qp.first) - grad).dot(flux_fun(qp.first) - grad);

//                    auto t_phi = cell_basis.eval_basis( qp.first );
//                    auto v = cell_dofs_n.dot(t_phi);
                   
//                    /* Compute L2-error */
//                    l2_error_vec[cell_i] += qp.second * (scal_fun(qp.first) - v) * (scal_fun(qp.first) - v);
//                }

//            }
//            else {
               
//                auto dofs = assembler.take_local_data(msh, cell, x_dof, locate(msh, cell));
//                 auto cell_dofs = dofs.head(cbs);

//                // uncut case
//                auto qps = integrate(msh, cell, 2*hho_di.cell_degree());
//                for (auto& qp : qps)
//                {
//                    /* Compute H1-error */
//                    auto t_dphi = cell_basis.eval_gradients( qp.first );
//                    Matrix<RealType, 1, 2> grad = Matrix<RealType, 1, 2>::Zero();

//                    for (size_t i = 1; i < cbs; i++ )
//                        grad += cell_dofs(i) * t_dphi.block(i, 0, 1, 2);

//                    flux_l2_error_vec[cell_i] += qp.second * (flux_fun(qp.first) - grad).dot(flux_fun(qp.first) - grad);

//                    auto t_phi = cell_basis.eval_basis( qp.first );
//                    auto v = cell_dofs.dot(t_phi);
                   
//                    /* Compute L2-error */
//                    l2_error_vec[cell_i] += qp.second * (scal_fun(qp.first) - v) * (scal_fun(qp.first) - v);
                   
//                }
//            }
//            cell_i++;
//        }
       
//        scalar_l2_error = std::accumulate(l2_error_vec.begin(), l2_error_vec.end(),0.0);
//        flux_l2_error = std::accumulate(flux_l2_error_vec.begin(), flux_l2_error_vec.end(),0.0);
//        tc.toc();
       
//     //    std::cout << bold << cyan << "Error completed: " << tc << " seconds" << reset << std::endl;
//        error_file << "Characteristic h size = " << std::setprecision(16) << h << std::endl;
//        error_file << "L2-norm error = " << std::setprecision(16) << std::sqrt(scalar_l2_error) << std::endl;
//        error_file << "H1-norm error = " << std::setprecision(16) << std::sqrt(flux_l2_error) << std::endl;
//        std::cout << "Characteristic h size = " << std::setprecision(16) << h << std::endl;
//        std::cout << "L2-norm error = " << std::setprecision(16) << std::sqrt(scalar_l2_error) << std::endl;
//        std::cout << "H1-norm error = " << std::setprecision(16) << std::sqrt(flux_l2_error) << std::endl;
//     //    error_file << std::endl;
//        error_file.flush();
       
//     }
    
//     /// Compute L2 and H1 errors for one field approximation
//     static void 
//     compute_errors_one_field(Mesh & msh, hho_degree_info & hho_di, one_field_interface_assembler<Mesh, std::function<double(const typename Mesh::point_type& )>> & assembler, Matrix<double, Dynamic, 1> & x_dof,std::function<double(const typename Mesh::point_type& )> scal_fun, std::function<Matrix<double, 1, 2>(const typename Mesh::point_type& )> flux_fun, std::ostream & error_file = std::cout){

//        timecounter tc;
//        tc.tic();

//        using RealType = double;

//        RealType scalar_l2_error = 0.0;
//        RealType flux_l2_error = 0.0;
//        size_t cell_i = 0;
//        RealType h = 10.0;
//        std::vector<RealType> l2_error_vec(msh.cells.size());
//        std::vector<RealType> flux_l2_error_vec(msh.cells.size());
//        for (auto& cell : msh.cells)
//        {
//            l2_error_vec[cell_i] = 0.0;
//             RealType h_l = diameter(msh, cell);
//            if (h_l < h) {
//                h = h_l;
//            }
           
//            cell_basis<cuthho_poly_mesh<RealType>, RealType> cell_basis(msh, cell, hho_di.cell_degree());
//            auto cbs = cell_basis.size();
//            if ( locate(msh, cell) == location::ON_INTERFACE )
//            {
               
//                auto dofs_n = assembler.take_local_data(msh, cell, x_dof, location::IN_NEGATIVE_SIDE);
//                auto dofs_p = assembler.take_local_data(msh, cell, x_dof, location::IN_POSITIVE_SIDE);

//                auto cell_dofs_n = dofs_n.head(cbs);
//                auto cell_dofs_p = dofs_p.head(cbs);
               
//                // negative side
//                auto qps_n = integrate(msh, cell, 2*hho_di.cell_degree(), location::IN_NEGATIVE_SIDE);
//                for (auto& qp : qps_n)
//                {
//                    /* Compute H1-error */
//                    auto t_dphi = cell_basis.eval_gradients( qp.first );
//                    Matrix<RealType, 1, 2> grad = Matrix<RealType, 1, 2>::Zero();

//                    for (size_t i = 1; i < cbs; i++ )
//                        grad += cell_dofs_n(i) * t_dphi.block(i, 0, 1, 2);

//                    flux_l2_error_vec[cell_i] += qp.second * (flux_fun(qp.first) - grad).dot(flux_fun(qp.first) - grad);
                   

//                    auto t_phi = cell_basis.eval_basis( qp.first );
//                    auto v = cell_dofs_n.dot(t_phi);
                   
//                    /* Compute L2-error */
//                    l2_error_vec[cell_i] += qp.second * (scal_fun(qp.first) - v) * (scal_fun(qp.first) - v);
//                }
               
//                // positive side
//                auto qps_p = integrate(msh, cell, 2*hho_di.cell_degree(), location::IN_POSITIVE_SIDE);
//                for (auto& qp : qps_p)
//                {
//                    /* Compute H1-error */
//                    auto t_dphi = cell_basis.eval_gradients( qp.first );
//                    Matrix<RealType, 1, 2> grad = Matrix<RealType, 1, 2>::Zero();

//                    for (size_t i = 1; i < cbs; i++ )
//                        grad += cell_dofs_n(i) * t_dphi.block(i, 0, 1, 2);

//                    flux_l2_error_vec[cell_i] += qp.second * (flux_fun(qp.first) - grad).dot(flux_fun(qp.first) - grad);

//                    auto t_phi = cell_basis.eval_basis( qp.first );
//                    auto v = cell_dofs_n.dot(t_phi);
                   
//                    /* Compute L2-error */
//                    l2_error_vec[cell_i] += qp.second * (scal_fun(qp.first) - v) * (scal_fun(qp.first) - v);
//                }

//            }else{
               
//                auto dofs = assembler.take_local_data(msh, cell, x_dof);
//                 auto cell_dofs = dofs.head(cbs);

//                // uncut case
//                auto qps = integrate(msh, cell, 2*hho_di.cell_degree());
//                for (auto& qp : qps)
//                {
//                    /* Compute H1-error */
//                    auto t_dphi = cell_basis.eval_gradients( qp.first );
//                    Matrix<RealType, 1, 2> grad = Matrix<RealType, 1, 2>::Zero();

//                    for (size_t i = 1; i < cbs; i++ )
//                        grad += cell_dofs(i) * t_dphi.block(i, 0, 1, 2);

//                    flux_l2_error_vec[cell_i] += qp.second * (flux_fun(qp.first) - grad).dot(flux_fun(qp.first) - grad);

//                    auto t_phi = cell_basis.eval_basis( qp.first );
//                    auto v = cell_dofs.dot(t_phi);
                   
//                    /* Compute L2-error */
//                    l2_error_vec[cell_i] += qp.second * (scal_fun(qp.first) - v) * (scal_fun(qp.first) - v);
                   
//                }
//            }
//            cell_i++;
//        }
       
//        scalar_l2_error = std::accumulate(l2_error_vec.begin(), l2_error_vec.end(),0.0);
//        flux_l2_error = std::accumulate(flux_l2_error_vec.begin(), flux_l2_error_vec.end(),0.0);
//        tc.toc();
       
//        std::cout << bold << cyan << "Error completed: " << tc << " seconds" << reset << std::endl;
//        error_file << "Characteristic h size = " << std::setprecision(16) << h << std::endl;
//        error_file << "L2-norm error = " << std::setprecision(16) << std::sqrt(scalar_l2_error) << std::endl;
//        error_file << "H1-norm error = " << std::setprecision(16) << std::sqrt(flux_l2_error) << std::endl;
//        error_file << std::endl;
//        error_file.flush();
       
//     }
    
//     /// Compute L2 and H1 errors for two fields approximation
//     static void 
//     compute_errors_two_fields(Mesh & msh, hho_degree_info & hho_di, two_fields_interface_assembler<Mesh, std::function<double(const typename Mesh::point_type& )>> & assembler, Matrix<double, Dynamic, 1> & x_dof,std::function<double(const typename Mesh::point_type& )> scal_fun, std::function<Matrix<double, 1, 2>(const typename Mesh::point_type& )> flux_fun, std::ostream & error_file = std::cout){

//         timecounter tc;
//         tc.tic();

//         using RealType = double;
 
//         RealType scalar_l2_error = 0.0;
//         RealType flux_l2_error = 0.0;
//         size_t cell_i = 0;
//         RealType h = 10.0;
//         std::vector<RealType> l2_error_vec(msh.cells.size());
//         std::vector<RealType> flux_l2_error_vec(msh.cells.size());
//         for (auto& cell : msh.cells)
//         {
//             l2_error_vec[cell_i] = 0.0;
//              RealType h_l = diameter(msh, cell);
//             if (h_l < h) {
//                 h = h_l;
//             }

            
//             {
//                 cell_basis<cuthho_poly_mesh<RealType>, RealType> cell_basis(msh, cell, hho_di.cell_degree());
//                 vector_cell_basis<cuthho_poly_mesh<RealType>, RealType> vec_cell_basis(msh, cell, hho_di.grad_degree());
                
//                 auto cbs = cell_basis.size();
//                 auto gbs = vec_cell_basis.size();
//                 if ( locate(msh, cell) == location::ON_INTERFACE )
//                 {
                    
//                     Matrix<RealType, Dynamic, 1> cell_dof_n = assembler.gather_cell_dof(msh,cell,x_dof,location::IN_NEGATIVE_SIDE);
//                     Matrix<RealType, Dynamic, 1> cell_dof_p = assembler.gather_cell_dof(msh,cell,x_dof,location::IN_POSITIVE_SIDE);
                    
//                     // negative side
//                     auto qps_n = integrate(msh, cell, 2*hho_di.cell_degree(), location::IN_NEGATIVE_SIDE);
//                     for (auto& qp : qps_n)
//                     {
//                         // scalar evaluation
//                         Matrix<RealType, Dynamic, 1> scal_cell_dof = cell_dof_n.tail(cbs);
//                         auto t_phi = cell_basis.eval_basis( qp.first );
//                         RealType uh = scal_cell_dof.dot( t_phi );
//                         l2_error_vec[cell_i] += qp.second * (scal_fun(qp.first) - uh) * (scal_fun(qp.first) - uh);
                    
//                         // flux evaluation
//                         Matrix<RealType, Dynamic, 1> vec_cell_dof = cell_dof_n.head(gbs);
//                         auto t_phi_v = vec_cell_basis.eval_basis( qp.first );
//                         Matrix<RealType, 1, 2> grad_uh = Matrix<RealType, 1, 2>::Zero();
//                         for (size_t i = 0; i < t_phi_v.rows(); i++){
//                           grad_uh = grad_uh + vec_cell_dof(i)*t_phi_v.block(i, 0, 1, 2);
//                         }
//                         auto grad_u_exact = flux_fun(qp.first);
//                         flux_l2_error_vec[cell_i] += qp.second * (grad_u_exact - grad_uh).dot(grad_u_exact - grad_uh);
                        

//                     }
                    
//                     // positive side
//                     auto qps_p = integrate(msh, cell, 2*hho_di.cell_degree(), location::IN_POSITIVE_SIDE);
//                     for (auto& qp : qps_p)
//                     {
//                         // scalar evaluation
//                         Matrix<RealType, Dynamic, 1> scal_cell_dof = cell_dof_p.tail(cbs);
//                         auto t_phi = cell_basis.eval_basis( qp.first );
//                         RealType uh = scal_cell_dof.dot( t_phi );
//                         l2_error_vec[cell_i] += qp.second * (scal_fun(qp.first) - uh) * (scal_fun(qp.first) - uh);
                        
//                         // flux evaluation
//                         Matrix<RealType, Dynamic, 1> vec_cell_dof = cell_dof_p.head(gbs);
//                         auto t_phi_v = vec_cell_basis.eval_basis( qp.first );
//                         Matrix<RealType, 1, 2> grad_uh = Matrix<RealType, 1, 2>::Zero();
//                         for (size_t i = 0; i < t_phi_v.rows(); i++){
//                           grad_uh = grad_uh + vec_cell_dof(i)*t_phi_v.block(i, 0, 1, 2);
//                         }
//                         auto grad_u_exact = flux_fun(qp.first);
//                         flux_l2_error_vec[cell_i] += qp.second * (grad_u_exact - grad_uh).dot(grad_u_exact - grad_uh);
                        
//                     }

//                 }else{
                    
//                     // uncut case
//                     // scalar evaluation
//                     Matrix<RealType, Dynamic, 1> cell_dof = assembler.gather_cell_dof(msh,cell,x_dof,locate(msh, cell));
//                     auto qps = integrate(msh, cell, 2*hho_di.cell_degree());
//                     for (auto& qp : qps)
//                     {
//                         Matrix<RealType, Dynamic, 1> scal_cell_dof = cell_dof.tail(cbs);
//                         auto t_phi = cell_basis.eval_basis( qp.first );
//                         RealType uh = scal_cell_dof.dot( t_phi );
//                         l2_error_vec[cell_i] += qp.second * (scal_fun(qp.first) - uh) * (scal_fun(qp.first) - uh);
                        
//                         // flux evaluation
//                         Matrix<RealType, Dynamic, 1> vec_cell_dof = cell_dof.head(gbs);
//                         auto t_phi_v = vec_cell_basis.eval_basis( qp.first );
//                         Matrix<RealType, 1, 2> grad_uh = Matrix<RealType, 1, 2>::Zero();
//                         for (size_t i = 0; i < t_phi_v.rows(); i++){
//                           grad_uh = grad_uh + vec_cell_dof(i)*t_phi_v.block(i, 0, 1, 2);
//                         }
//                         auto grad_u_exact = flux_fun(qp.first);
//                         flux_l2_error_vec[cell_i] += qp.second * (grad_u_exact - grad_uh).dot(grad_u_exact - grad_uh);
                        
//                     }
//                 }

//             }
//             cell_i++;
//         }
//         scalar_l2_error = std::accumulate(l2_error_vec.begin(), l2_error_vec.end(),0.0);
//         flux_l2_error = std::accumulate(flux_l2_error_vec.begin(), flux_l2_error_vec.end(),0.0);
//         tc.toc();
        
//         std::cout << bold << yellow << "         Error completed: " << tc << " seconds" << reset << std::endl;
//         error_file << "Characteristic h size = " << std::setprecision(16) << h << std::endl;
//         error_file << "L2-norm error = " << std::setprecision(16) << std::sqrt(scalar_l2_error) << std::endl;
//         error_file << "H1-norm error = " << std::setprecision(16) << std::sqrt(flux_l2_error) << std::endl;
//         error_file.flush();
        
//     }

    // WRITE SCRIPT PYTHON CV TESTS
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////

    static void 
    write_conv_sol(const std::string& txtFilename) {

        if (txtFilename.size() < 4 || txtFilename.substr(txtFilename.size() - 4) != ".txt") {
            std::cerr << "Erreur : Le fichier d'entrée n'a pas une extension .txt valide." << std::endl;
            return;
        }
        
        std::string pyFilename = txtFilename.substr(0, txtFilename.size() - 4) + ".py";
        
        std::ofstream pyFile(pyFilename);
        if (!pyFile.is_open()) {
            std::cerr << "Erreur : Impossible de créer le fichier Python." << std::endl;
            return;
        }
        
      pyFile << "import matplotlib.pyplot as plt\n";
      pyFile << "import numpy as np\n";
      pyFile << "import sys\n\n";
      
      pyFile << "error_type = \"L2\"\n";
      pyFile << "if len(sys.argv) > 1:\n";
      pyFile << "    if sys.argv[1] in [\"L2\", \"H1\"]:\n";
      pyFile << "        error_type = sys.argv[1]\n";
      pyFile << "    else:\n";
      pyFile << "        print(\"Argument invalide. Utilisez 'L2' ou 'H1'.\")\n";
      pyFile << "        sys.exit(1)\n\n";
      
      pyFile << "filename = \"" << txtFilename << "\"\n";
      pyFile << "with open(filename, \"r\") as file:\n";
      pyFile << "    lines = file.readlines()\n\n";
      
      pyFile << "results_L2 = {}\n";
      pyFile << "results_H1 = {}\n";
      pyFile << "h_values = {}\n\n";
      
      pyFile << "current_degree = None\n";
      pyFile << "for line in lines:\n";
      pyFile << "    line = line.strip()\n";
      pyFile << "    if line.startswith(\"Polynomial degree k :\"):\n";
      pyFile << "        current_degree = int(line.split(\":\")[1].strip())\n";
      pyFile << "        if current_degree not in results_L2:\n";
      pyFile << "            results_L2[current_degree] = []\n";
      pyFile << "            results_H1[current_degree] = []\n";
      pyFile << "            h_values[current_degree] = []\n";
      pyFile << "    elif line.startswith(\"Characteristic h size =\"):\n";
      pyFile << "        h = float(line.split(\"=\")[1].strip())\n";
      pyFile << "        h_values[current_degree].append(h)\n";
      pyFile << "    elif line.startswith(\"L2-norm error =\"):\n";
      pyFile << "        L2_error = float(line.split(\"=\")[1].strip())\n";
      pyFile << "        results_L2[current_degree].append(L2_error)\n";
      pyFile << "    elif line.startswith(\"H1-norm error =\"):\n";
      pyFile << "        H1_error = float(line.split(\"=\")[1].strip())\n";
      pyFile << "        results_H1[current_degree].append(H1_error)\n\n";
      
      pyFile << "plt.figure(figsize=(10, 6))\n\n";
      
      pyFile << "if error_type == \"L2\":\n";
      pyFile << "    results = results_L2\n";
      pyFile << "    ylabel = \"L2-norm Error\"\n";
      pyFile << "    title = \"L2-norm Error vs h\"\n";
      pyFile << "elif error_type == \"H1\":\n";
      pyFile << "    results = results_H1\n";
      pyFile << "    ylabel = \"H1-norm Error\"\n";
      pyFile << "    title = \"H1-norm Error vs h\"\n\n";
      
      pyFile << "for degree in sorted(results.keys()):\n";
      pyFile << "    h_values_np = np.array(h_values[degree])\n";
      pyFile << "    results_np = np.array(results[degree])\n";
      pyFile << "    log_h = np.log(h_values_np)\n";
      pyFile << "    log_error = np.log(results_np)\n";
      pyFile << "    slope, _ = np.polyfit(log_h, log_error, 1)\n";
      pyFile << "    plt.loglog(h_values_np, results_np, marker='o', label=f\"k={degree} {ylabel} (rate={slope:.2f})\")\n\n";
      
      pyFile << "plt.xlabel(\"h\", fontsize=12)\n";
      pyFile << "plt.ylabel(ylabel, fontsize=12)\n";
      pyFile << "plt.title(title, fontsize=14)\n";
      pyFile << "plt.legend()\n";
      pyFile << "plt.grid(which=\"both\", linestyle=\"--\", linewidth=0.5)\n";
      pyFile << "plt.tight_layout()\n\n";
      
      pyFile << "plt.show()\n";
      pyFile.close();
      
    }

//     static void 
//     write_conv_grad(const std::string& txtFilename) {

//         if (txtFilename.size() < 4 || txtFilename.substr(txtFilename.size() - 4) != ".txt") {
//             std::cerr << "Erreur : Le fichier d'entrée n'a pas une extension .txt valide." << std::endl;
//             return;
//         }

//         std::string pyFilename = txtFilename.substr(0, txtFilename.size() - 4) + ".py";

//         std::ofstream pyFile(pyFilename);
//         if (!pyFile.is_open()) {
//             std::cerr << "Erreur : Impossible de créer le fichier Python." << std::endl;
//             return;
//         }

//         pyFile << "import matplotlib.pyplot as plt\n";
//         pyFile << "import numpy as np\n\n";

//         pyFile << "filename = \"" << txtFilename << "\"\n";
//         pyFile << "with open(filename, \"r\") as file:\n";
//         pyFile << "    lines = file.readlines()\n\n";

//         pyFile << "results = {}\n";
//         pyFile << "degree = None\n\n";

//         pyFile << "# Lire les données du fichier\n";
//         pyFile << "for line in lines:\n";
//         pyFile << "    line = line.strip()\n";
//         pyFile << "    if line.startswith(\"Polynomial degree k :\"):\n";
//         pyFile << "        degree = int(line.split(\":\")[1].strip())\n";
//         pyFile << "        if degree not in results:\n";
//         pyFile << "            results[degree] = {\"h\": [], \"error\": []}\n";
//         pyFile << "    elif line.startswith(\"Characteristic h size =\"):\n";
//         pyFile << "        h = float(line.split(\"=\")[1].strip())\n";
//         pyFile << "        results[degree][\"h\"].append(h)\n";
//         pyFile << "    elif line.startswith(\"L2-norm grad error =\"):\n";
//         pyFile << "        error = float(line.split(\"=\")[1].strip())\n";
//         pyFile << "        results[degree][\"error\"].append(error)\n\n";

//         pyFile << "plt.figure(figsize=(10, 6))\n\n";

//         pyFile << "for degree, values in sorted(results.items()):\n";
//         pyFile << "    h_values = np.array(values[\"h\"])\n";
//         pyFile << "    error_values = np.array(values[\"error\"])\n\n";

//         pyFile << "    log_h = np.log(h_values)\n";
//         pyFile << "    log_error = np.log(error_values)\n";
//         pyFile << "    slope, _ = np.polyfit(log_h, log_error, 1)  # Pente (taux de convergence)\n\n";

//         pyFile << "    plt.loglog(h_values, error_values, marker='o', label=f\"k={degree} (rate={slope:.2f})\")\n\n";

//         pyFile << "plt.xlabel(\"h\", fontsize=12)\n";
//         pyFile << "plt.ylabel(\"L2-norm grad error\", fontsize=12)\n";
//         pyFile << "plt.legend()\n";
//         pyFile << "plt.grid(which=\"both\", linestyle=\"--\", linewidth=0.5)\n";
//         pyFile << "plt.tight_layout()\n\n";

//         pyFile << "# Afficher le graphique\n";
//         pyFile << "plt.show()\n";

//         pyFile.close();
//     }

//     static void 
//     write_conv_grad_grad(const std::string& txtFilename) {
        
//     if (txtFilename.size() < 4 || txtFilename.substr(txtFilename.size() - 4) != ".txt") {
//         std::cerr << "Erreur : Le fichier d'entrée n'a pas une extension .txt valide." << std::endl;
//         return;
//     }

//     std::string pyFilename = txtFilename.substr(0, txtFilename.size() - 4) + ".py";

//     std::ofstream pyFile(pyFilename);
//     if (!pyFile.is_open()) {
//         std::cerr << "Erreur : Impossible de créer le fichier Python." << std::endl;
//         return;
//     }

//     pyFile << "import matplotlib.pyplot as plt\n";
//     pyFile << "import numpy as np\n\n";

//     pyFile << "filename = \"" << txtFilename << "\"\n";
//     pyFile << "with open(filename, \"r\") as file:\n";
//     pyFile << "    lines = file.readlines()\n\n";

//     pyFile << "results = {}\n";
//     pyFile << "degree = None\n\n";

//     pyFile << "# Lire les données du fichier\n";
//     pyFile << "for line in lines:\n";
//     pyFile << "    line = line.strip()\n";
//     pyFile << "    if line.startswith(\"Polynomial degree k :\"):\n";
//     pyFile << "        degree = int(line.split(\":\")[1].strip())\n";
//     pyFile << "        if degree not in results:\n";
//     pyFile << "            results[degree] = {\"h\": [], \"error\": []}\n";
//     pyFile << "    elif line.startswith(\"Characteristic h size =\"):\n";
//     pyFile << "        h = float(line.split(\"=\")[1].strip())\n";
//     pyFile << "        results[degree][\"h\"].append(h)\n";
//     pyFile << "    elif line.startswith(\"L2-norm error =\"):\n";
//     pyFile << "        error = float(line.split(\"=\")[1].strip())\n";
//     pyFile << "        results[degree][\"error\"].append(error)\n\n";

//     pyFile << "plt.figure(figsize=(10, 6))\n\n";

//     pyFile << "for degree, values in sorted(results.items()):\n";
//     pyFile << "    h_values = np.array(values[\"h\"])\n";
//     pyFile << "    error_values = np.array(values[\"error\"])\n\n";

//     pyFile << "    # Vérifiez si les longueurs correspondent\n";
//     pyFile << "    if len(h_values) != len(error_values):\n";
//     pyFile << "        print(f\"Attention : Les longueurs de 'h' et 'error' ne correspondent pas pour le degré {degree}.\")\n";
//     pyFile << "        continue\n\n";

//     pyFile << "    log_h = np.log(h_values)\n";
//     pyFile << "    log_error = np.log(error_values)\n";
//     pyFile << "    slope, _ = np.polyfit(log_h, log_error, 1)  # Pente (taux de convergence)\n\n";

//     pyFile << "    plt.loglog(h_values, error_values, marker='o', label=f\"k={degree} (rate={slope:.2f})\")\n\n";

//     pyFile << "plt.xlabel(\"h\", fontsize=12)\n";
//     pyFile << "plt.ylabel(\"L2-norm error\", fontsize=12)\n";
//     pyFile << "plt.legend()\n";
//     pyFile << "plt.grid(which=\"both\", linestyle=\"--\", linewidth=0.5)\n";
//     pyFile << "plt.tight_layout()\n\n";

//     pyFile << "# Afficher le graphique\n";
//     pyFile << "plt.show()\n";
        
//         pyFile.close();
//     }

//     // SILO
//     /////////////////////////////////////////////////////////////////////////////
//     /////////////////////////////////////////////////////////////////////////////
//     /////////////////////////////////////////////////////////////////////////////
    
//     static void write_silo_one_field(std::string silo_file_name, size_t it, Mesh & msh, hho_degree_info & hho_di, interface_assembler<Mesh, std::function<double(const typename Mesh::point_type& )>> & assembler, Matrix<double, Dynamic, 1> & x_dof,
//     std::function<double(const typename Mesh::point_type& )> scal_fun, bool cell_centered_Q = false){

//         timecounter tc;
//         tc.tic();
        
//         auto dim = 2;
//         auto num_cells = msh.cells.size();
//         auto num_points = msh.points.size();
//         using RealType = double;
//         std::vector<RealType> exact_u, approx_u;
        
//         if (cell_centered_Q) {
//             exact_u.reserve( num_cells );
//             approx_u.reserve( num_cells );

//             size_t cell_i = 0;
//             for (auto& cell : msh.cells)
//             {
//                 auto bar = barycenter(msh, cell);
//                 exact_u.push_back( scal_fun(bar) );
                
//                 // scalar evaluation
//                 {
//                     cell_basis<cuthho_poly_mesh<RealType>, RealType> cell_basis(msh, cell, hho_di.cell_degree());
//                     if ( locate(msh, cell) == location::ON_INTERFACE )
//                     {
//                         // negative side
//                         {
//                             Matrix<RealType, Dynamic, 1> scalar_cell_dof = assembler.gather_cell_dof(msh,cell,x_dof,location::IN_NEGATIVE_SIDE);
//                             auto t_phi = cell_basis.eval_basis( bar );
//                             RealType uh = scalar_cell_dof.dot( t_phi );
//                             approx_u.push_back(uh);
//                         }
                        
// //                        // positive side
// //                        {
// //                            Matrix<RealType, Dynamic, 1> scalar_cell_dof = assembler.gather_cell_dof(msh,cell,x_dof,location::IN_POSITIVE_SIDE);
// //                            auto t_phi = cell_basis.eval_basis( bar );
// //                            RealType uh = scalar_cell_dof.dot( t_phi );
// //                            approx_u.push_back(uh);
// //                        }
                        
//                     }else{
//                         Matrix<RealType, Dynamic, 1> scalar_cell_dof = assembler.gather_cell_dof(msh,cell,x_dof,locate(msh, cell));
//                         auto t_phi = cell_basis.eval_basis( bar );
//                         RealType uh = scalar_cell_dof.dot( t_phi );
//                         approx_u.push_back(uh);
//                     }
                    

//                 }
//                 cell_i++;
//             }

//         }else{

//             exact_u.reserve( num_points );
//             approx_u.reserve( num_points );

//             // scan for selected cells, common cells are discardable
//             std::map<size_t, size_t> node_to_cell;
//             size_t cell_i = 0;
//             for (auto& cell : msh.cells)
//             {
//                 auto cell_nodes = nodes(msh,cell);
//                 size_t n_p = cell_nodes.size();
//                 for (size_t l = 0; l < n_p; l++)
//                 {
//                     auto node = cell_nodes[l];
//                     node_to_cell[node.ptid] = cell_i;
//                 }
//                 cell_i++;
//             }
            

//             for (auto& node_id : node_to_cell)
//             {
//                 auto bar = msh.points.at(node_id.first);
//                 exact_u.push_back( scal_fun(bar) );

//                 cell_i = node_id.second;
//                 auto cell = msh.cells.at(cell_i);

//                 // scalar evaluation
//                 {
//                     cell_basis<cuthho_poly_mesh<RealType>, RealType> cell_basis(msh, cell, hho_di.cell_degree());
//                     if ( locate(msh, cell) == location::ON_INTERFACE )
//                     {
//                         auto node = msh.nodes.at(node_id.first);
                        
//                         if (locate(msh, node) == location::IN_NEGATIVE_SIDE)
//                         // negative side
//                         {
//                             Matrix<RealType, Dynamic, 1> scalar_cell_dof = assembler.gather_cell_dof(msh,cell,x_dof,location::IN_NEGATIVE_SIDE);
//                             auto t_phi = cell_basis.eval_basis( bar );
//                             RealType uh = scalar_cell_dof.dot( t_phi );
//                             approx_u.push_back(uh);
//                         }else
//                         // positive side
//                         {
//                             Matrix<RealType, Dynamic, 1> scalar_cell_dof = assembler.gather_cell_dof(msh,cell,x_dof,location::IN_POSITIVE_SIDE);
//                             auto t_phi = cell_basis.eval_basis( bar );
//                             RealType uh = scalar_cell_dof.dot( t_phi );
//                             approx_u.push_back(uh);
//                         }

//                     }else{
//                         Matrix<RealType, Dynamic, 1> scalar_cell_dof = assembler.gather_cell_dof(msh,cell,x_dof,locate(msh, cell));
//                         auto t_phi = cell_basis.eval_basis( bar );
//                         RealType uh = scalar_cell_dof.dot( t_phi );
//                         approx_u.push_back(uh);
//                     }


//                 }

//             }

//         }

//         silo_database silo;
//         silo_file_name += std::to_string(it) + ".silo";
//         silo.create(silo_file_name.c_str());
//         silo.add_mesh(msh, "mesh");
//         if (cell_centered_Q) {
//             silo.add_variable("mesh", "v", exact_u.data(), exact_u.size(), zonal_variable_t);
//             silo.add_variable("mesh", "vh", approx_u.data(), approx_u.size(), zonal_variable_t);
//         }else{
//             silo.add_variable("mesh", "v", exact_u.data(), exact_u.size(), nodal_variable_t);
//             silo.add_variable("mesh", "vh", approx_u.data(), approx_u.size(), nodal_variable_t);
//         }

//         silo.close();
//         tc.toc();
//         // std::cout << std::endl;
//         // std::cout << bold << cyan << "Silo file rendered in : " << tc << " seconds" << reset << std::endl;
//     }
    
//     static void 
//     write_silo_two_fields(std::string silo_file_name, size_t it, Mesh & msh, hho_degree_info & hho_di, two_fields_interface_assembler<Mesh, std::function<double(const typename Mesh::point_type& )>> & assembler, Matrix<double, Dynamic, 1> & x_dof, std::function<double(const typename Mesh::point_type& )> scal_fun, bool cell_centered_Q = false) {

//         timecounter tc;
//         tc.tic();
        
//         auto dim = 2;
//         auto num_cells = msh.cells.size();
//         auto num_points = msh.points.size();
//         using RealType = double;
//         std::vector<RealType> exact_u, approx_u;
        
//         exact_u.reserve( num_points );
//         approx_u.reserve( num_points );
//         // scan for selected cells, common cells are discardable
//         std::map<size_t, size_t> node_to_cell;
//         size_t cell_i = 0;
//         for (auto& cell : msh.cells) {
//             auto cell_nodes = nodes(msh,cell);
//             size_t n_p = cell_nodes.size();
//             for (size_t l = 0; l < n_p; l++) {
//                 auto node = cell_nodes[l];
//                 node_to_cell[node.ptid] = cell_i;
//             }
//             cell_i++;
//         }    
        
//         for (auto& node_id : node_to_cell) {
//             auto bar = msh.points.at(node_id.first);
//             exact_u.push_back( scal_fun(bar) );
//             cell_i = node_id.second;
//             auto cell = msh.cells.at(cell_i);
//             // scalar evaluation
//             cell_basis<cuthho_poly_mesh<RealType>, RealType> cell_basis(msh, cell, hho_di.cell_degree());
//             auto cbs = cell_basis.size();
//             if ( locate(msh, cell) == location::ON_INTERFACE ) {
//                 auto node = msh.nodes.at(node_id.first);
//                 if (locate(msh, node) == location::IN_NEGATIVE_SIDE) { // negative side
//                     Matrix<RealType, Dynamic, 1> cell_dof = assembler.gather_cell_dof(msh,cell,x_dof,location::IN_NEGATIVE_SIDE);
//                     Matrix<RealType, Dynamic, 1> scal_cell_dof = cell_dof.tail(cbs);
//                     auto t_phi = cell_basis.eval_basis( bar );
//                     RealType uh = scal_cell_dof.dot( t_phi );
//                     approx_u.push_back(uh);
//                 }
//                 else { // positive side
//                     Matrix<RealType, Dynamic, 1> cell_dof = assembler.gather_cell_dof(msh,cell,x_dof,location::IN_POSITIVE_SIDE);
//                     Matrix<RealType, Dynamic, 1> scal_cell_dof = cell_dof.tail(cbs);
//                     auto t_phi = cell_basis.eval_basis( bar );
//                     RealType uh = scal_cell_dof.dot( t_phi );
//                     approx_u.push_back(uh);
//                 }
//             }
//             else {
//                 Matrix<RealType, Dynamic, 1> cell_dof = assembler.gather_cell_dof(msh,cell,x_dof,locate(msh, cell));
//                 Matrix<RealType, Dynamic, 1> scal_cell_dof = cell_dof.tail(cbs);
//                 auto t_phi = cell_basis.eval_basis( bar );
//                 RealType uh = scal_cell_dof.dot( t_phi );
//                 approx_u.push_back(uh);
//             }
//         }
        
//         silo_database silo;
//         silo_file_name += std::to_string(it) + ".silo";
//         silo.create(silo_file_name.c_str());
//         silo.add_mesh(msh, "mesh");
//         silo.add_variable("mesh", "v", exact_u.data(), exact_u.size(), nodal_variable_t);
//         silo.add_variable("mesh", "vh", approx_u.data(), approx_u.size(), nodal_variable_t);
//         silo.close();
//         tc.toc();
//         std::cout << bold << yellow << "         Silo file rendered in : " << tc << " seconds" << reset << std::endl;
//     }
        
//     #ifndef subcell_centering
//     template<typename testType>
//     static void write_silo_poly_ext(std::string silo_file_name, size_t it, Mesh & msh, 
//     hho_degree_info & hho_di, Matrix<double, Dynamic, 1> &x_dof, testType &test_case,
//     interface_assembler<Mesh, std::function<double(const typename Mesh::point_type& )>> & assembler) {

//         using RealType = double;
//         timecounter tc;
//         tc.tic();

//         auto dim = 2;
//         auto num_points = msh.nodes.size();
//         auto msh_nodes = msh.nodes;
        
//         std::vector<RealType> exact_u, approx_u;
//         std::vector<RealType> exact_Gx, exact_Gy;
//         std::vector<RealType> approx_Gx, approx_Gy;
        
//         exact_u.reserve(num_points);
//         approx_u.reserve(num_points);
//         exact_Gx.reserve(num_points);
//         exact_Gy.reserve(num_points);
//         approx_Gx.reserve(num_points);
//         approx_Gy.reserve(num_points);

//         auto scal_fun = test_case.sol_fun;
//         auto flux_fun = test_case.sol_grad;

//         auto cell_table = assembler.get_cell_table();
//         for (auto& node : msh_nodes) {

//             // Coords of the node
//             auto pt = msh.points.at(node.ptid); 
            
//             // Exact functions
//             exact_u.push_back( scal_fun(pt) );
//             exact_Gx.push_back( flux_fun(pt)(0) );
//             exact_Gy.push_back( flux_fun(pt)(1) );
            
//             RealType uh = 0.0;
//             RealType Gx = 0.0;
//             RealType Gy = 0.0;
//             auto cells = find_cells(pt, msh);
//             auto size_cells = cells.size();
//             for (auto& cell : cells ) {
//                 auto cl = msh.cells[cell];
//                 cell_basis<cuthho_poly_mesh<RealType>, RealType> cell_basis(msh, cl, hho_di.cell_degree());
//                 vector_cell_basis<cuthho_poly_mesh<RealType>, RealType> vector_cell_basis(msh, cl, hho_di.grad_degree());
//                 auto local_dofs = cell_basis.size();
//                 auto grad_dofs_size = vector_cell_basis.size();
//                 auto offset_cl = cell_table.at(offset(msh, cl));
//                 if (!is_cut(msh, cl)) {
//                     // Solution evaluation
//                     Matrix<RealType, Dynamic, 1> scalar_cell_dof = assembler.gather_cell_dof(msh, cl, x_dof, locate(msh, cl));
//                     auto t_phi = cell_basis.eval_basis(pt);
//                     auto dt_phi = cell_basis.eval_gradients(pt);
//                     uh += scalar_cell_dof.dot(t_phi);
//                     Gx += scalar_cell_dof.dot(dt_phi.col(0));
//                     Gy += scalar_cell_dof.dot(dt_phi.col(1));
//                 }
//                 else {
//                     if (locate(msh, node) == location::IN_NEGATIVE_SIDE) {
//                         // Solution evaluation
//                         Matrix<RealType, Dynamic, 1> scalar_cell_dof = assembler.gather_cell_dof(msh, cl, x_dof, location::IN_NEGATIVE_SIDE);
//                         auto t_phi = cell_basis.eval_basis(pt);
//                         auto dt_phi = cell_basis.eval_gradients(pt);
//                         uh += scalar_cell_dof.dot(t_phi);
//                         Gx += scalar_cell_dof.dot(dt_phi.col(0));
//                         Gy += scalar_cell_dof.dot(dt_phi.col(1));
//                     }
//                     else if (locate(msh, node) == location::IN_POSITIVE_SIDE) {
//                         // Solution evaluation
//                         Matrix<RealType, Dynamic, 1> scalar_cell_dof = assembler.gather_cell_dof(msh, cl, x_dof, location::IN_POSITIVE_SIDE);
//                         auto t_phi = cell_basis.eval_basis(pt);
//                         auto dt_phi = cell_basis.eval_gradients(pt);
//                         uh += scalar_cell_dof.dot(t_phi);
//                         Gx += scalar_cell_dof.dot(dt_phi.col(0));
//                         Gy += scalar_cell_dof.dot(dt_phi.col(1));
//                     }
//                 }
//             }
//             uh /= size_cells;
//             Gx /= size_cells;
//             Gy /= size_cells;   
//             approx_u.push_back(uh);
//             approx_Gx.push_back(Gx);
//             approx_Gy.push_back(Gy);
//         }
        
//         silo_database silo;
//         silo_file_name += ".silo";
//         silo.create(silo_file_name.c_str());
//         silo.add_mesh(msh, "mesh");
//         silo.add_variable("mesh", "1_v",   exact_u.data(),   exact_u.size(),  nodal_variable_t);
//         silo.add_variable("mesh", "2_vh",  approx_u.data(),  approx_u.size(), nodal_variable_t);
//         silo.add_variable("mesh", "3_Gx",  exact_Gx.data(),  exact_Gx.size(), nodal_variable_t);
//         silo.add_variable("mesh", "4_Gy",  exact_Gy.data(),  exact_Gy.size(), nodal_variable_t);
//         silo.add_variable("mesh", "5_Gxh", approx_Gx.data(), exact_Gx.size(), nodal_variable_t);
//         silo.add_variable("mesh", "6_Gyh", approx_Gy.data(), exact_Gy.size(), nodal_variable_t);
//         silo.close();
//         tc.toc();
//         std::cout << bold << yellow << "         Silo file - SOLUTION - rendered in : " << tc << " seconds" << reset << std::endl;
//     }
//     #else
//     template<typename testType>
//     static void write_silo_poly_ext(std::string silo_file_name, size_t it, Mesh & msh, 
//     hho_degree_info & hho_di, Matrix<double, Dynamic, 1> &x_dof, testType &test_case,
//     interface_assembler<Mesh, std::function<double(const typename Mesh::point_type& )>> & assembler) {

//         using RealType = double;
//         timecounter tc;
//         tc.tic();

//         auto dim = 2;
//         auto num_points = msh.nodes.size();
//         auto msh_nodes = msh.nodes;
        
//         std::vector<RealType> exact_u, approx_u;
//         std::vector<RealType> exact_Gx, exact_Gy;
//         std::vector<RealType> approx_Gx, approx_Gy;
//         std::vector<RealType> Rec_Gx, Rec_Gy;
        
//         exact_u.reserve(num_points);
//         approx_u.reserve(num_points);
//         exact_Gx.reserve(num_points);
//         exact_Gy.reserve(num_points);
//         approx_Gx.reserve(num_points);
//         approx_Gy.reserve(num_points);

//         auto scal_fun = test_case.sol_fun;
//         auto flux_fun = test_case.sol_grad;
//         auto cell_table = assembler.get_cell_table();
//         for (auto& node : msh_nodes) {

//             // Coords of the node
//             auto pt = msh.points.at(node.ptid); 
            
//             // Exact functions
//             exact_u.push_back( scal_fun(pt) );
//             exact_Gx.push_back( flux_fun(pt)(0) );
//             exact_Gy.push_back( flux_fun(pt)(1) );
            
//             RealType uh = 0.0;
//             RealType Gx = 0.0;
//             RealType Gy = 0.0;
//             auto cells = find_cells(pt, msh);
//             auto size_cells = cells.size();
//             for (auto& cell : cells ) {
//                 auto cl = msh.cells[cell];
//                 auto offset_cl = cell_table.at(offset(msh, cl));
//                 if (!is_cut(msh, cl)) {
//                     cell_basis<cuthho_poly_mesh<RealType>, RealType> cell_basis(msh, cl, hho_di.cell_degree());
//                     vector_cell_basis<cuthho_poly_mesh<RealType>, RealType> vector_cell_basis(msh, cl, hho_di.grad_degree());
//                     auto local_dofs = cell_basis.size();
//                     auto grad_dofs_size = vector_cell_basis.size();
//                     // Solution evaluation
//                     Matrix<RealType, Dynamic, 1> scalar_cell_dof = assembler.gather_cell_dof(msh, cl, x_dof, locate(msh, cl));
//                     auto t_phi = cell_basis.eval_basis(pt);
//                     auto dt_phi = cell_basis.eval_gradients(pt);
//                     uh += scalar_cell_dof.dot(t_phi);
//                     Gx += scalar_cell_dof.dot(dt_phi.col(0));
//                     Gy += scalar_cell_dof.dot(dt_phi.col(1));
//                 }
//                 else {
//                     if (locate(msh, node) == location::IN_NEGATIVE_SIDE) {
//                         cut_cell_basis<cuthho_poly_mesh<RealType>, RealType> cell_basis(msh, cl, hho_di.cell_degree(), location::IN_NEGATIVE_SIDE);
//                         cut_vector_cell_basis<cuthho_poly_mesh<RealType>, RealType> vector_cell_basis(msh, cl, hho_di.grad_degree(), location::IN_NEGATIVE_SIDE);
//                         auto local_dofs = cell_basis.size();
//                         auto grad_dofs_size = vector_cell_basis.size();
//                         // Solution evaluation
//                         Matrix<RealType, Dynamic, 1> scalar_cell_dof = assembler.gather_cell_dof(msh, cl, x_dof, location::IN_NEGATIVE_SIDE);
//                         auto t_phi = cell_basis.eval_basis(pt);
//                         auto dt_phi = cell_basis.eval_gradients(pt);
//                         uh += scalar_cell_dof.dot(t_phi);
//                         Gx += scalar_cell_dof.dot(dt_phi.col(0));
//                         Gy += scalar_cell_dof.dot(dt_phi.col(1));
//                     }
//                     else if (locate(msh, node) == location::IN_POSITIVE_SIDE) {
//                         cut_cell_basis<cuthho_poly_mesh<RealType>, RealType> cell_basis(msh, cl, hho_di.cell_degree(), location::IN_POSITIVE_SIDE);
//                         cut_vector_cell_basis<cuthho_poly_mesh<RealType>, RealType> vector_cell_basis(msh, cl, hho_di.grad_degree(), location::IN_POSITIVE_SIDE);
//                         auto local_dofs = cell_basis.size();
//                         auto grad_dofs_size = vector_cell_basis.size();
//                         // Solution evaluation
//                         Matrix<RealType, Dynamic, 1> scalar_cell_dof = assembler.gather_cell_dof(msh, cl, x_dof, location::IN_POSITIVE_SIDE);
//                         auto t_phi = cell_basis.eval_basis(pt);
//                         auto dt_phi = cell_basis.eval_gradients(pt);
//                         uh += scalar_cell_dof.dot(t_phi);
//                         Gx += scalar_cell_dof.dot(dt_phi.col(0));
//                         Gy += scalar_cell_dof.dot(dt_phi.col(1));
//                     }
//                 }
//             }
//             uh /= size_cells;
//             Gx /= size_cells;
//             Gy /= size_cells;   
//             approx_u.push_back(uh);
//             approx_Gx.push_back(Gx);
//             approx_Gy.push_back(Gy);
//         }
        
//         silo_database silo;
//         silo_file_name += ".silo";
//         silo.create(silo_file_name.c_str());
//         silo.add_mesh(msh, "mesh");
//         silo.add_variable("mesh", "1_v",   exact_u.data(),   exact_u.size(),  nodal_variable_t);
//         silo.add_variable("mesh", "2_vh",  approx_u.data(),  approx_u.size(), nodal_variable_t);
//         silo.add_variable("mesh", "3_Gx",  exact_Gx.data(),  exact_Gx.size(), nodal_variable_t);
//         silo.add_variable("mesh", "4_Gy",  exact_Gy.data(),  exact_Gy.size(), nodal_variable_t);
//         silo.add_variable("mesh", "5_Gxh", approx_Gx.data(), exact_Gx.size(), nodal_variable_t);
//         silo.add_variable("mesh", "6_Gyh", approx_Gy.data(), exact_Gy.size(), nodal_variable_t);
//         silo.close();
//         tc.toc();
//         std::cout << bold << yellow << "         Silo file - SOLUTION CENTERED - rendered in : " << tc << " seconds" << reset << std::endl;
//     }
//     #endif 

//     // RECORD SENSORS ACOUSTIC WAVES 
//     /////////////////////////////////////////////////////////////////////////////
//     /////////////////////////////////////////////////////////////////////////////
//     /////////////////////////////////////////////////////////////////////////////

//     /// Record data at provided point for one field approximation
//     static void 
//     record_data_acoustic_one_field(size_t it, std::pair<typename Mesh::point_type,size_t> & pt_cell_index, Mesh & msh, hho_degree_info & hho_di, interface_assembler<Mesh, std::function<double(const typename Mesh::point_type& )>> & assembler, Matrix<double, Dynamic, 1> & x_dof, std::ostream & seismogram_file = std::cout){

//         timecounter tc;
//         tc.tic();

//         using RealType = double;
//         auto dim = 2;

//         Matrix<double, Dynamic, 1> vh = Matrix<double, Dynamic, 1>::Zero(2, 1);

//         typename Mesh::point_type pt = pt_cell_index.first;
        
//         if(pt_cell_index.second == -1){
//             std::set<size_t> cell_indexes = find_cells(pt, msh, true);
//             size_t cell_index = pick_cell(pt, msh, cell_indexes, true);
//             assert(cell_index != -1);
//             pt_cell_index.second = cell_index;
//             seismogram_file << "\"Time\"" << "," << "\"vhx\"" << "," << "\"vhy\"" << std::endl;
//         }

//         {
//             size_t cell_ind = pt_cell_index.second;
//             // scalar evaluation
//             auto cell = msh.cells.at(cell_ind);

//             // scalar evaluation
//             {
//                 cell_basis<cuthho_poly_mesh<RealType>, RealType> cell_basis(msh, cell, hho_di.cell_degree());
//                 if ( locate(msh, cell) == location::ON_INTERFACE )
//                 {
//                     auto node = msh.nodes.at(0);
//                     throw std::invalid_argument("Recoding at cut cell. Not implemented.");
//                     if (locate(msh, node) == location::IN_NEGATIVE_SIDE)
//                     // negative side
//                     {
//                         Matrix<RealType, Dynamic, 1> scalar_cell_dof = assembler.gather_cell_dof(msh,cell,x_dof,location::IN_NEGATIVE_SIDE);
//                         auto t_dphi = cell_basis.eval_gradients( pt );
//                         size_t cbs = cell_basis.size();
//                            Matrix<RealType, 1, 2> grad = Matrix<RealType, 1, 2>::Zero();
//                            for (size_t i = 1; i < cbs; i++ )
//                                grad += scalar_cell_dof(i) * t_dphi.block(i, 0, 1, 2);
//                         vh = grad;
//                     }else
//                     // positive side
//                     {
//                         Matrix<RealType, Dynamic, 1> scalar_cell_dof = assembler.gather_cell_dof(msh,cell,x_dof,location::IN_POSITIVE_SIDE);
//                         auto t_dphi = cell_basis.eval_gradients( pt );
//                         size_t cbs = cell_basis.size();
//                            Matrix<RealType, 1, 2> grad = Matrix<RealType, 1, 2>::Zero();
//                            for (size_t i = 1; i < cbs; i++ )
//                                grad += scalar_cell_dof(i) * t_dphi.block(i, 0, 1, 2);
//                         vh = grad;
//                     }

//                 }else{
//                     Matrix<RealType, Dynamic, 1> scalar_cell_dof = assembler.gather_cell_dof(msh,cell,x_dof,locate(msh, cell));
//                     auto t_dphi = cell_basis.eval_gradients( pt );
//                     size_t cbs = cell_basis.size();
//                        Matrix<RealType, 1, 2> grad = Matrix<RealType, 1, 2>::Zero();
//                        for (size_t i = 1; i < cbs; i++ )
//                            grad += scalar_cell_dof(i) * t_dphi.block(i, 0, 1, 2);
                    
//                     vh = grad;
//                 }


//             }
//         }
//         tc.toc();
//         std::cout << bold << cyan << "Value recorded: " << tc << " seconds" << reset << std::endl;
//         seismogram_file << it << "," << std::setprecision(16) <<  vh(0,0) << "," << std::setprecision(16) <<  vh(1,0) << std::endl;
//         seismogram_file.flush();

//     }
    
//     /// Record data at provided point for two fields approximation
//     static void 
//     record_data_acoustic_two_fields(size_t it, std::pair<typename Mesh::point_type,size_t> & pt_cell_index, Mesh & msh, hho_degree_info & hho_di, two_fields_interface_assembler<Mesh, std::function<double(const typename Mesh::point_type& )>> & assembler, Matrix<double, Dynamic, 1> & x_dof, std::ostream & seismogram_file = std::cout){

//         timecounter tc;
//         tc.tic();

//         using RealType = double;
//         auto dim = 2;

//         Matrix<double, Dynamic, 1> vh = Matrix<double, Dynamic, 1>::Zero(2, 1);

//         typename Mesh::point_type pt = pt_cell_index.first;
        
//         if(pt_cell_index.second == -1){
//             std::set<size_t> cell_indexes = find_cells(pt, msh, true);
//             size_t cell_index = pick_cell(pt, msh, cell_indexes, true);
//             assert(cell_index != -1);
//             pt_cell_index.second = cell_index;
//             seismogram_file << "\"Time\"" << "," << "\"vhx\"" << "," << "\"vhy\"" << std::endl;
//         }

//         {
//             size_t cell_ind = pt_cell_index.second;
//             auto cell = msh.cells.at(cell_ind);

//             // flux evaluation
//             {
//                 vector_cell_basis<cuthho_poly_mesh<RealType>, RealType> vec_cell_basis(msh, cell, hho_di.grad_degree());
//                 auto gbs = vec_cell_basis.size();
                
//                 if ( locate(msh, cell) == location::ON_INTERFACE )
//                 {
//                     auto node = msh.nodes.at(0);
//                     throw std::invalid_argument("Recoding at cut cell. Not implemented.");
//                     Matrix<RealType, Dynamic, 1> cell_dof_n = assembler.gather_cell_dof(msh,cell,x_dof,location::IN_NEGATIVE_SIDE);
//                     Matrix<RealType, Dynamic, 1> cell_dof_p = assembler.gather_cell_dof(msh,cell,x_dof,location::IN_POSITIVE_SIDE);
                    
//                     if (locate(msh, node) == location::IN_NEGATIVE_SIDE)
//                     // negative side
//                     {
//                         Matrix<RealType, Dynamic, 1> vec_cell_dof = cell_dof_n.head(gbs);
//                         auto t_phi_v = vec_cell_basis.eval_basis( pt );
//                         Matrix<RealType, 1, 2> grad_uh = Matrix<RealType, 1, 2>::Zero();
//                         for (size_t i = 0; i < t_phi_v.rows(); i++){
//                           grad_uh = grad_uh + vec_cell_dof(i)*t_phi_v.block(i, 0, 1, 2);
//                         }
//                         vh = grad_uh;
//                     }else
//                     // positive side
//                     {
//                         Matrix<RealType, Dynamic, 1> vec_cell_dof = cell_dof_p.head(gbs);
//                         auto t_phi_v = vec_cell_basis.eval_basis( pt );
//                         Matrix<RealType, 1, 2> grad_uh = Matrix<RealType, 1, 2>::Zero();
//                         for (size_t i = 0; i < t_phi_v.rows(); i++){
//                           grad_uh = grad_uh + vec_cell_dof(i)*t_phi_v.block(i, 0, 1, 2);
//                         }
//                         vh = grad_uh;
//                     }

//                 }else{
//                     Matrix<RealType, Dynamic, 1> cell_dof = assembler.gather_cell_dof(msh,cell,x_dof,locate(msh, cell));
                    
//                     Matrix<RealType, Dynamic, 1> vec_cell_dof = cell_dof.head(gbs);
//                     auto t_phi_v = vec_cell_basis.eval_basis( pt );
//                     Matrix<RealType, 1, 2> grad_uh = Matrix<RealType, 1, 2>::Zero();
//                     for (size_t i = 0; i < t_phi_v.rows(); i++){
//                       grad_uh = grad_uh + vec_cell_dof(i)*t_phi_v.block(i, 0, 1, 2);
//                     }
//                     vh = grad_uh;
//                 }


//             }
//         }
//         tc.toc();
//         std::cout << bold << cyan << "Value recorded: " << tc << " seconds" << reset << std::endl;
//         seismogram_file << it << "," << std::setprecision(16) <<  vh(0,0) << "," << std::setprecision(16) <<  vh(1,0) << std::endl;
//         seismogram_file.flush();

//     }

};


#endif /* postprocessor_hpp */
