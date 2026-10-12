
//  Created by Omar Durán
//  Contributor: Romain Mottier

#pragma once
#ifndef scal_vec_analytic_functions_hpp
#define scal_vec_analytic_functions_hpp

class scal_vec_analytic_functions {
    
    public:
    
    /// Enumerate defining the function type
    enum EFunctionType { 
      EFunctionNonPolynomial            = 0,  
      EFunctionQuadraticInTime          = 1,
      EFunctionQuadraticInSpace         = 2, 
      reproduction_acoustic             = 3,
      reproduction_elastic              = 4,
      EFunctionNonPolynomial_paper      = 5, 
      EFunctionCubicInTimeAcoustic      = 6,
      EFunctionQuarticInTimeAcoustic    = 7,
      EFunctionQuadraticInSpaceAcoustic = 8,
      EFunctionPlaneWaveAcoustic = 9,
      EFunctionInterfaceCoupling = 10,
      // single-medium solutions on (0,1)^2 vanishing on the boundary, rho = 1, c = 1 (acoustic), mu = lambda = 1
      // (elastic); reproduction_acoustic and reproduction_elastic are eqs. (60) and (78) of these papers:
      EFunctionAcousticQuadraticInSpace = 11,   // Burman-Duran-Ern, CAMC 2022, eq. (61)
      EFunctionAcousticNonPolynomial = 12,      // Burman-Duran-Ern, CAMC 2022, eq. (62)
      EFunctionElasticQuadraticInSpace = 13,    // Burman-Duran-Ern-Steins, JSC 2021, eq. (79)
      EFunctionElasticNonPolynomial = 14,       // Burman-Duran-Ern-Steins, JSC 2021, eq. (80)
      // (61) and (62) with cos(sqrt2 pi t) instead of sin(sqrt2 pi t) (p(0) = p, p_t(0) = 0): the HHO-Newmark
      // results of the published CAMC 2022 figures (Fig. 1 right, Fig. 2 left) are those of these solutions
      EFunctionAcousticQuadraticInSpaceCos = 15,
      EFunctionAcousticNonPolynomialCos = 16};
      
      scal_vec_analytic_functions() {
          m_function_type = EFunctionNonPolynomial;
      }

      // EFunctionNonPolynomial_paper (Mottier-Ern-Guillot, sinusoidal test case): spatial frequency omega and temporal
      // frequency theta of u = x^2 sin(omega pi x) sin(omega pi y) sin(theta pi t) (acoustic) and
      // u_x = u_y = x^2 cos(omega pi x / 2) sin(omega pi y) cos(theta pi t) (elastic); default omega = 1, theta = 10
      void set_paper_frequencies(double omega, double theta) {
          m_paper_w = omega;
          m_paper_theta = theta * M_PI;
      }
      
      ~scal_vec_analytic_functions() {  
      }
          
      void set_function_type(EFunctionType function_type) {
          m_function_type = function_type;
      }
      
      std::function<disk::static_vector<double, 2>
      (const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_u(double & t) {
          
          switch (m_function_type) {
              
              case EFunctionNonPolynomial: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,ux,uy;
                      x = pt.x();
                      y = pt.y();
                      ux = x*x*std::sin(M_PI*y)*std::cos((M_PI*x)/2.)*std::cos(std::sqrt(2.0)*M_PI*t);
                      uy = x*x*std::sin(M_PI*y)*std::cos((M_PI*x)/2.)*std::cos(std::sqrt(2.0)*M_PI*t);
                      disk::static_vector<double, 2> u{ux,uy};
                      return u;
                  };
              }
              break;
              
              case EFunctionQuadraticInTime: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,ux,uy;
                      x = pt.x();
                      y = pt.y();
                      ux = x*t*t*std::sin(M_PI*x)*std::sin(M_PI*y);
                      uy = x*t*t*std::sin(M_PI*x)*std::sin(M_PI*y);
                      disk::static_vector<double, 2> u{ux,uy};
                      return u;
                  };
              }
              break;
              
              case EFunctionQuadraticInSpace: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,ux,uy;
                      x = pt.x();
                      y = pt.y();
                      ux = (1 + x)*x*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t);
                      uy = (1 + x)*x*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t);
                      disk::static_vector<double, 2> u{ux,uy};
                      return u;
                  };
              }
              break;
              
              case reproduction_elastic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x,y,ux,uy;
                      x = pt.x();
                      y = pt.y();
                      ux = -t*t*std::cos(M_PI*y)*std::sin(M_PI*x);
                      uy =  t*t*std::cos(M_PI*x)*std::sin(M_PI*y);
                      disk::static_vector<double, 2> u{ux,uy};
                      return u;
                  };
              }
              break;
              
              case EFunctionNonPolynomial_paper: {
                  return [&t, w0 = m_paper_w, th0 = m_paper_theta](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x,y,ux,uy, w, theta;
                      x = pt.x();
                      y = pt.y();
                      w = w0; 
                      theta = th0; 
                      // w = 5.0;
                      // theta = std::sqrt(2)*M_PI;
                      ux = x*x*std::sin(M_PI*w*y)*std::cos(theta*t)*std::cos((M_PI*w*x)/2.);
                      uy = x*x*std::sin(M_PI*w*y)*std::cos(theta*t)*std::cos((M_PI*w*x)/2.);
                      disk::static_vector<double, 2> u{ux,uy};
                      return u;
                  };
              }
              break;
              
              // Coupled solution with non-vanishing interface traces (interface x = 0, elastic x < 0: rho = 1,
              // mu = lambda = 1; acoustic x > 0: rho = 1, c = 1), T(t) = sin(w t), w = sqrt(2) pi:
              //   u   = ((sin(pi x)/(3 pi) + cos(pi x)) sin(pi y) T', -pi x cos(pi y) T'),
              //   phi = (cos(pi x) - 2 pi sin(pi x)) sin(pi y) T.
              // It satisfies the transmission conditions of the coupled formulations, d_x phi / rho = d_t u_x,
              // sigma_xx = d_t phi, sigma_xy = 0 on x = 0, with nonzero traces; the acoustic source vanishes.
              case EFunctionInterfaceCoupling: {
                  // displacement
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x = pt.x(), y = pt.y();
                      double T = std::sqrt(2.0)*M_PI*std::cos(std::sqrt(2.0)*M_PI*t);
                      double vx = (std::sin(M_PI*x)/(3.0*M_PI) + std::cos(M_PI*x))*std::sin(M_PI*y)*T;
                      double vy = -M_PI*x*std::cos(M_PI*y)*T;
                      disk::static_vector<double, 2> v{vx,vy};
                      return v;
                  };
              }
              break;
              
              case reproduction_acoustic: case EFunctionAcousticQuadraticInSpace: case EFunctionAcousticNonPolynomial:
              case EFunctionAcousticQuadraticInSpaceCos: case EFunctionAcousticNonPolynomialCos: {
                  // no elastic part
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      disk::static_vector<double, 2> z{0.0, 0.0};
                      return z;
                  };
              }
              break;
              
              case EFunctionElasticQuadraticInSpace: {
                  // displacement (Burman-Duran-Ern-Steins, JSC 2021, eq. (79))
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x = pt.x(), y = pt.y();
                      double rx = std::sin(std::sqrt(2.0)*M_PI*t)*x*(1.0-x)*y*(1.0-y), ry = rx;
                      disk::static_vector<double, 2> r{rx, ry};
                      return r;
                  };
              }
              break;
              
              case EFunctionElasticNonPolynomial: {
                  // displacement (Burman-Duran-Ern-Steins, JSC 2021, eq. (80))
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x = pt.x(), y = pt.y();
                      double rx = std::sin(std::sqrt(2.0)*M_PI*t)*(-std::sin(M_PI*x)*std::cos(M_PI*y)), ry = std::sin(std::sqrt(2.0)*M_PI*t)*(std::cos(M_PI*x)*std::sin(M_PI*y));
                      disk::static_vector<double, 2> r{rx, ry};
                      return r;
                  };
              }
              break;
              
              default: {
                  std::cout << " Function not implemented " << std::endl;
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      disk::static_vector<double, 2> u;
                      return u;
                  };
              }
              break;
          }   
      }
      
      std::function<disk::static_vector<double, 2>
      (const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_v(double & t){
          
          switch (m_function_type) {
              
              case EFunctionNonPolynomial: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,vx,vy;
                      x = pt.x();
                      y = pt.y();
                      vx = -(std::sqrt(2.0)*M_PI*x*x*std::cos((M_PI*x)/2.)*std::sin(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*y));
                      vy = -(std::sqrt(2.0)*M_PI*x*x*std::cos((M_PI*x)/2.)*std::sin(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*y));
                      disk::static_vector<double, 2> v{vx,vy};
                      return v;
                  };
              }
              break;
              
              case EFunctionQuadraticInTime: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,vx,vy;
                      x = pt.x();
                      y = pt.y();
                      vx = 2.0*t*x*std::sin(M_PI*x)*std::sin(M_PI*y);
                      vy = 2.0*t*x*std::sin(M_PI*x)*std::sin(M_PI*y);
                      disk::static_vector<double, 2> v{vx,vy};
                      return v;
                  };
              }
              break;
              
              case EFunctionQuadraticInSpace: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,vx,vy;
                      x = pt.x();
                      y = pt.y();
                      vx = std::sqrt(2.0)*M_PI*(1 + x)*x*x*(1 - y)*y*std::cos(std::sqrt(2.0)*M_PI*t);
                      vy = std::sqrt(2.0)*M_PI*(1 + x)*x*x*(1 - y)*y*std::cos(std::sqrt(2.0)*M_PI*t);
                      disk::static_vector<double, 2> v{vx,vy};
                      return v;
                  };
              }
              break;
              
              case reproduction_elastic:{
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x,y,vx,vy;
                      x = pt.x();
                      y = pt.y();
                      vx = -2*t*std::cos(M_PI*y)*std::sin(M_PI*x);
                      vy =  2*t*std::cos(M_PI*x)*std::sin(M_PI*y);
                      disk::static_vector<double, 2> v{vx,vy};
                      return v;
                  };
              }
              break;
              
              case EFunctionNonPolynomial_paper:{
                  return [&t, w0 = m_paper_w, th0 = m_paper_theta](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x,y,vx,vy,w,theta;
                      x = pt.x();
                      y = pt.y();
                      w = w0; 
                      theta = th0; 
                      // w = 5.0;
                      // theta = std::sqrt(2)*M_PI;
                      vx = -(theta*x*x*std::sin(theta*t)*std::sin(M_PI*w*y)*std::cos((M_PI*w*x)/2.));
                      vy = -(theta*x*x*std::sin(theta*t)*std::sin(M_PI*w*y)*std::cos((M_PI*w*x)/2.));
                      disk::static_vector<double, 2> v{vx,vy};
                      return v;
                  };
              }
              break;
              
              case EFunctionInterfaceCoupling: {
                  // velocity
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x = pt.x(), y = pt.y();
                      double T = -2.0*M_PI*M_PI*std::sin(std::sqrt(2.0)*M_PI*t);
                      double vx = (std::sin(M_PI*x)/(3.0*M_PI) + std::cos(M_PI*x))*std::sin(M_PI*y)*T;
                      double vy = -M_PI*x*std::cos(M_PI*y)*T;
                      disk::static_vector<double, 2> v{vx,vy};
                      return v;
                  };
              }
              break;
              
              case reproduction_acoustic: case EFunctionAcousticQuadraticInSpace: case EFunctionAcousticNonPolynomial:
              case EFunctionAcousticQuadraticInSpaceCos: case EFunctionAcousticNonPolynomialCos: {
                  // no elastic part
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      disk::static_vector<double, 2> z{0.0, 0.0};
                      return z;
                  };
              }
              break;
              
              case EFunctionElasticQuadraticInSpace: {
                  // velocity
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x = pt.x(), y = pt.y();
                      double rx = std::sqrt(2.0)*M_PI*std::cos(std::sqrt(2.0)*M_PI*t)*x*(1.0-x)*y*(1.0-y), ry = rx;
                      disk::static_vector<double, 2> r{rx, ry};
                      return r;
                  };
              }
              break;
              
              case EFunctionElasticNonPolynomial: {
                  // velocity
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x = pt.x(), y = pt.y();
                      double rx = std::sqrt(2.0)*M_PI*std::cos(std::sqrt(2.0)*M_PI*t)*(-std::sin(M_PI*x)*std::cos(M_PI*y)), ry = std::sqrt(2.0)*M_PI*std::cos(std::sqrt(2.0)*M_PI*t)*(std::cos(M_PI*x)*std::sin(M_PI*y));
                      disk::static_vector<double, 2> r{rx, ry};
                      return r;
                  };
              }
              break;
              
              default: {
                  std::cout << " Function not implemented " << std::endl;
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      disk::static_vector<double, 2> v;
                      return v;
                  };
              }
              break;    
          }
      }
      
      std::function<disk::static_vector<double, 2>
      (const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_a(double & t){
          
          switch (m_function_type) {
              
              case EFunctionNonPolynomial: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,ax,ay;
                      x = pt.x();
                      y = pt.y();
                      ax = -2*M_PI*M_PI*x*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::cos((M_PI*x)/2.)*std::sin(M_PI*y);
                      ay = -2*M_PI*M_PI*x*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::cos((M_PI*x)/2.)*std::sin(M_PI*y);
                      disk::static_vector<double, 2> a{ax,ay};
                      return a;
                  };
              }
              break;
              
              case EFunctionQuadraticInTime: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,ax,ay;
                      x = pt.x();
                      y = pt.y();
                      ax = 2.0*x*std::sin(M_PI*x)*std::sin(M_PI*y);
                      ay = 2.0*x*std::sin(M_PI*x)*std::sin(M_PI*y);
                      disk::static_vector<double, 2> a{ax,ay};
                      return a;
                  };
              }
              break;
              
              case EFunctionQuadraticInSpace: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,ax,ay;
                      x = pt.x();
                      y = pt.y();
                      ax = -2*M_PI*M_PI*(1 + x)*x*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t);
                      ay = -2*M_PI*M_PI*(1 + x)*x*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t);
                      disk::static_vector<double, 2> a{ax,ay};
                      return a;
                  };
              }
              break;
              
              case reproduction_elastic:{
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x,y,ax,ay;
                      x = pt.x();
                      y = pt.y();
                      ax = -2*std::cos(M_PI*y)*std::sin(M_PI*x);
                      ay = 2*std::cos(M_PI*x)*std::sin(M_PI*y);
                      disk::static_vector<double, 2> a{ax,ay};
                      return a;
                  };
              }
              break;
              
              case EFunctionNonPolynomial_paper: {
                  return [&t, w0 = m_paper_w, th0 = m_paper_theta](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,ax,ay,w,theta;
                      x = pt.x();
                      y = pt.y();
                      w = w0; 
                      theta = th0; 
                      // w = 5.0;
                      // theta = std::sqrt(2)*M_PI;
                      ax = -theta*theta*x*x*std::sin(M_PI*w*y)*std::cos(theta*t)*std::cos((M_PI*w*x)/2.);
                      ay = -theta*theta*x*x*std::sin(M_PI*w*y)*std::cos(theta*t)*std::cos((M_PI*w*x)/2.);
                      disk::static_vector<double, 2> a{ax,ay};
                      return a;
                  };
              }
              break;
              
              case EFunctionInterfaceCoupling: {
                  // acceleration
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x = pt.x(), y = pt.y();
                      double T = -2.0*std::sqrt(2.0)*M_PI*M_PI*M_PI*std::cos(std::sqrt(2.0)*M_PI*t);
                      double vx = (std::sin(M_PI*x)/(3.0*M_PI) + std::cos(M_PI*x))*std::sin(M_PI*y)*T;
                      double vy = -M_PI*x*std::cos(M_PI*y)*T;
                      disk::static_vector<double, 2> v{vx,vy};
                      return v;
                  };
              }
              break;
              
              case reproduction_acoustic: case EFunctionAcousticQuadraticInSpace: case EFunctionAcousticNonPolynomial:
              case EFunctionAcousticQuadraticInSpaceCos: case EFunctionAcousticNonPolynomialCos: {
                  // no elastic part
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      disk::static_vector<double, 2> z{0.0, 0.0};
                      return z;
                  };
              }
              break;
              
              case EFunctionElasticQuadraticInSpace: {
                  // acceleration
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x = pt.x(), y = pt.y();
                      double rx = -2.0*M_PI*M_PI*std::sin(std::sqrt(2.0)*M_PI*t)*x*(1.0-x)*y*(1.0-y), ry = rx;
                      disk::static_vector<double, 2> r{rx, ry};
                      return r;
                  };
              }
              break;
              
              case EFunctionElasticNonPolynomial: {
                  // acceleration
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x = pt.x(), y = pt.y();
                      double rx = -2.0*M_PI*M_PI*std::sin(std::sqrt(2.0)*M_PI*t)*(-std::sin(M_PI*x)*std::cos(M_PI*y)), ry = -2.0*M_PI*M_PI*std::sin(std::sqrt(2.0)*M_PI*t)*(std::cos(M_PI*x)*std::sin(M_PI*y));
                      disk::static_vector<double, 2> r{rx, ry};
                      return r;
                  };
              }
              break;
              
              default: {
                  std::cout << " Function not implemented " << std::endl;
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      disk::static_vector<double, 2> f;
                      return f;
                  };
              }
              break;   
          }   
      }
      
      std::function<disk::static_vector<double, 2>
      (const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_f(double & t){
          
          switch (m_function_type) {
              
              case EFunctionNonPolynomial: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,fx,fy;
                      x = pt.x();
                      y = pt.y();
                      fx = -(std::cos(std::sqrt(2.0)*M_PI*t)*(-4*M_PI*x*std::sin((M_PI*x)/2.)*(M_PI*x*std::cos(M_PI*y) + 6*std::sin(M_PI*y)) + std::cos((M_PI*x)/2.)*(16*M_PI*x*std::cos(M_PI*y) + (24 + M_PI*M_PI*x*x)*std::sin(M_PI*y))))/4.0;
                      fy = (std::cos(std::sqrt(2.0)*M_PI*t)*(4*M_PI*x*std::sin((M_PI*x)/2.)*(M_PI*x*std::cos(M_PI*y) + 2*std::sin(M_PI*y)) + std::cos((M_PI*x)/2.)*(-16*M_PI*x*std::cos(M_PI*y) + (-8 + 5*M_PI*M_PI*x*x)*std::sin(M_PI*y))))/4.;
                      disk::static_vector<double, 2> f{fx,fy};
                      return f;
                  };
              }
              break;
              
              case EFunctionQuadraticInTime: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,fx,fy;
                      x = pt.x();
                      y = pt.y();
                      fx = -2 * ( M_PI*t*t*std::cos(M_PI*x)*( M_PI*x*std::cos(M_PI*y) + 3*std::sin(M_PI*y) ) +
                      std::sin(M_PI*x)*(M_PI*t*t*std::cos(M_PI*y) - (1+2*M_PI*M_PI*t*t)*x*std::sin(M_PI*y)));
                      fy = (1 + M_PI*M_PI*t*t)*x*std::cos(M_PI*(x - y)) - (1+3*M_PI*M_PI*t*t)*x*std::cos(M_PI*(x+y))
                      - 2*M_PI*t*t*std::sin(M_PI*(x+y));
                      disk::static_vector<double, 2> f{fx,fy};
                      return f;
                  };
              }
              break;
              
              case EFunctionQuadraticInSpace: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,fx,fy;
                      x = pt.x();
                      y = pt.y();
                      fx = 2 * ( x * (-2+(-2+x)*x) -3*y -x*(5+x*(-6+M_PI*M_PI*(1+x)))*y
                           + (3+x*(9+M_PI*M_PI*x*(1+x)))*y*y)*std::sin(std::sqrt(2.0)*M_PI*t);
                      fy = 2 * ( x*x*(6 + M_PI*M_PI*(-1 + y))*y + (-1 + y)*y + x*x*x*(3 + M_PI*M_PI*(-1 + y)*y)
                           + x*(-2 + y + 3*y*y))*std::sin(std::sqrt(2.0)*M_PI*t);
                      disk::static_vector<double, 2> f{fx,fy};
                      return f;
                  };
              }
              break;
              
              case reproduction_elastic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x,y,fx,fy;
                      x = pt.x();
                      y = pt.y();
                      fx = -2*(1 + M_PI*M_PI*t*t)*std::cos(M_PI*y)*std::sin(M_PI*x);
                      fy = 2*(1 + M_PI*M_PI*t*t)*std::cos(M_PI*x)*std::sin(M_PI*y);
                      disk::static_vector<double, 2> f{fx,fy};
                      return f;
                  };
              }
              break;
              
              case EFunctionNonPolynomial_paper: {
                  return [&t, w0 = m_paper_w, th0 = m_paper_theta](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,fx,fy,w,theta;
                      x = pt.x();
                      y = pt.y();
                      w = w0; 
                      theta = th0; 
                      // w = 5.0;
                      // theta = std::sqrt(2)*M_PI;
                      
                      fx = - theta*theta*x*x*std::sin(M_PI*w*y)*std::cos(theta*t)*std::cos(M_PI*w*x/2.0) 
                           + M_PI*M_PI*w*w*x*x*std::sin(M_PI*w*x/2)*std::cos(theta*t)*std::cos(M_PI*w*y) 
                           + 7*M_PI*M_PI*w*w*x*x*std::sin(M_PI*w*y)*std::cos(theta*t)*cos(M_PI*w*x/2)/4.0  
                           + 6*M_PI*w*x*std::sin(M_PI*w*x/2)*std::sin(M_PI*w*y)*std::cos(theta*t) 
                           - 4*M_PI*w*x*std::cos(theta*t)*std::cos(M_PI*w*x/2)*std::cos(M_PI*w*y) 
                           - 6*std::sin(M_PI*w*y)*std::cos(theta*t)*std::cos(M_PI*w*x/2);
                           
                      fy = - theta*theta*x*x*std::sin(M_PI*w*y)*std::cos(theta*t)*std::cos(M_PI*w*x/2.0) 
                           + M_PI*M_PI*w*w*x*x*std::sin(M_PI*w*x/2)*std::cos(theta*t)*std::cos(M_PI*w*y) 
                           + 13*M_PI*M_PI*w*w*x*x*std::sin(M_PI*w*y)*std::cos(theta*t)*cos(M_PI*w*x/2)/4.0 
                           + 2*M_PI*w*x*std::sin(M_PI*w*x/2)*std::sin(M_PI*w*y)*std::cos(theta*t) 
                           - 4*M_PI*w*x*std::cos(theta*t)*std::cos(M_PI*w*x/2)*std::cos(M_PI*w*y) 
                           - 2*std::sin(M_PI*w*y)*std::cos(theta*t)*std::cos(M_PI*w*x/2);
                      disk::static_vector<double, 2> f{fx,fy};
                      return f;
                  };
              }
              break;
              
              case EFunctionInterfaceCoupling: {
                  // rho d_tt u - div sigma(u)
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x = pt.x(), y = pt.y();
                      double T = std::sqrt(2.0)*M_PI*std::cos(std::sqrt(2.0)*M_PI*t);
                      double fx = ((2.0*M_PI/3.0)*std::sin(M_PI*x) + 2.0*M_PI*M_PI*std::cos(M_PI*x) - 2.0*M_PI*M_PI)*std::sin(M_PI*y)*T;
                      double fy = (-M_PI*M_PI*M_PI*x - (2.0*M_PI/3.0)*std::cos(M_PI*x) + 2.0*M_PI*M_PI*std::sin(M_PI*x))*std::cos(M_PI*y)*T;
                      disk::static_vector<double, 2> f{fx,fy};
                      return f;
                  };
              }
              break;
              
              case reproduction_acoustic: case EFunctionAcousticQuadraticInSpace: case EFunctionAcousticNonPolynomial:
              case EFunctionAcousticQuadraticInSpaceCos: case EFunctionAcousticNonPolynomialCos: {
                  // no elastic part
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      disk::static_vector<double, 2> z{0.0, 0.0};
                      return z;
                  };
              }
              break;
              
              case EFunctionElasticQuadraticInSpace: {
                  // rho u_tt - div sigma(u), rho = mu = lambda = 1
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x = pt.x(), y = pt.y();
                      double X = x*(1.0-x)*y*(1.0-y);
                      double rx = std::sin(std::sqrt(2.0)*M_PI*t)*(-2.0*M_PI*M_PI*X + 6.0*y*(1.0-y) - 2.0*(1.0-2.0*x)*(1.0-2.0*y) + 2.0*x*(1.0-x));
                      double ry = std::sin(std::sqrt(2.0)*M_PI*t)*(-2.0*M_PI*M_PI*X + 2.0*y*(1.0-y) - 2.0*(1.0-2.0*x)*(1.0-2.0*y) + 6.0*x*(1.0-x));
                      disk::static_vector<double, 2> r{rx, ry};
                      return r;
                  };
              }
              break;
              
              case EFunctionElasticNonPolynomial: {
                  // rho u_tt - div sigma(u) = 0
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x = pt.x(), y = pt.y();
                      double rx = 0.0*x*y*t, ry = 0.0;
                      disk::static_vector<double, 2> r{rx, ry};
                      return r;
                  };
              }
              break;
              
              default:
              {
                  std::cout << " Function not implemented " << std::endl;
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      disk::static_vector<double, 2> f;
                      return f;
                  };
              }
              break;
          }    
      }
      
      std::function<disk::static_matrix<double,2,2>
      (const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_sigma(double & t){
          
          switch (m_function_type) {
              
              case EFunctionNonPolynomial: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_matrix<double,2,2> {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      disk::static_matrix<double,2,2> sigma = disk::static_matrix<double,2,2>::Zero(2,2);
                      sigma(0,0) = M_PI*x*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::cos((M_PI*x)/2.)*std::cos(M_PI*y) + 4*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::cos((M_PI*x)/2.)*std::sin(M_PI*y) -
                      M_PI*x*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::sin((M_PI*x)/2.)*std::sin(M_PI*y) +
                      (4*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::cos((M_PI*x)/2.)*std::sin(M_PI*y) - M_PI*x*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::sin((M_PI*x)/2.)*std::sin(M_PI*y))/2.0;
                           
                      sigma(0,1) = M_PI*x*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::cos((M_PI*x)/2.)*std::cos(M_PI*y) + 2*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::cos((M_PI*x)/2.)*std::sin(M_PI*y) -
                      (M_PI*x*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::sin((M_PI*x)/2.)*std::sin(M_PI*y))/2.;
                      sigma(1,0) = sigma(0,1);
                      
                      sigma(1,1) = 3*M_PI*x*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::cos((M_PI*x)/2.)*std::cos(M_PI*y) +
                      (4*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::cos((M_PI*x)/2.)*std::sin(M_PI*y) - M_PI*x*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::sin((M_PI*x)/2.)*std::sin(M_PI*y))/2.;
                      return sigma;
                  };
              }
              break;
              
              case EFunctionQuadraticInTime: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_matrix<double,2,2> {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      disk::static_matrix<double,2,2> sigma = disk::static_matrix<double,2,2>::Zero(2,2);
                      sigma(0,0) = t*t*(M_PI*x*std::cos(M_PI*y)*std::sin(M_PI*x)
                      + 3*(M_PI*x*std::cos(M_PI*x) + std::sin(M_PI*x))*std::sin(M_PI*y));
                      
                      sigma(0,1) = t*t*(M_PI*x*std::cos(M_PI*y)*std::sin(M_PI*x)
                      + (M_PI*x*std::cos(M_PI*x) + std::sin(M_PI*x))*std::sin(M_PI*y));
                      
                      sigma(1,0) = sigma(0,1);
                      
                      sigma(1,1) = t*t*(3*M_PI*x*std::cos(M_PI*y)*std::sin(M_PI*x)
                      + (M_PI*x*std::cos(M_PI*x) + std::sin(M_PI*x))*std::sin(M_PI*y));
                      return sigma;
                  };
              }
              break;
              
              case EFunctionQuadraticInSpace: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_matrix<double,2,2> {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      disk::static_matrix<double,2,2> sigma = disk::static_matrix<double,2,2>::Zero(2,2);
                      sigma(0,0) =   2*x*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t)
                      + 4*x*(1+x)*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t)
                      + (2*x*x*(1 + x)*(1 - y)*std::sin(std::sqrt(2.0)*M_PI*t)
                      - 2*x*x*(1 + x)*y*std::sin(std::sqrt(2.0)*M_PI*t))/2.0
                      + (2*x*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t)
                      + 4*x*(1 + x)*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t))/2.0;
                      
                      sigma(0,1) =   x*x*(1 + x)*(1 - y)*std::sin(std::sqrt(2.0)*M_PI*t)
                      - x*x*(1 + x)*y*std::sin(std::sqrt(2.0)*M_PI*t)
                      + x*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t)
                      + 2*x*(1 + x)*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t);
                      
                      sigma(1,0) = sigma(0,1);
                      
                      sigma(1,1) =   2*x*x*(1 + x)*(1 - y)*std::sin(std::sqrt(2.0)*M_PI*t)
                      - 2*x*x*(1 + x)*y*std::sin(std::sqrt(2.0)*M_PI*t)
                      + (2*x*x*(1 + x)*(1 - y)*std::sin(std::sqrt(2.0)*M_PI*t)
                      - 2*x*x*(1 + x)*y*std::sin(std::sqrt(2.0)*M_PI*t))/2.0
                      + (2*x*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t)
                      + 4*x*(1 + x)*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t))/2.0;
                      return sigma;
                  };
              }
              break;
              
              case reproduction_elastic:
              {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_matrix<double,2,2> {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      disk::static_matrix<double,2,2> sigma = disk::static_matrix<double,2,2>::Zero(2,2);
                      sigma(0,0) = -2*M_PI*t*t*std::cos(M_PI*x)*std::cos(M_PI*y);
                      sigma(1,1) = 2*M_PI*t*t*std::cos(M_PI*x)*std::cos(M_PI*y);
                      return sigma;
                  };
              }
              break;
              
              case EFunctionNonPolynomial_paper: {
                  return [&t, w0 = m_paper_w, th0 = m_paper_theta](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_matrix<double,2,2> {
                      double x,y,w,theta;
                      x = pt.x();
                      y = pt.y();
                      
                      w = w0; 
                      theta = th0; 
                      // w = 5.0;
                      // theta = std::sqrt(2)*M_PI;
                      
                      disk::static_matrix<double,2,2> sigma = disk::static_matrix<double,2,2>::Zero(2,2);
                      sigma(0,0) = - 3.0*M_PI*w*x*x*std::sin(M_PI*w*x/2.0)*std::sin(M_PI*w*y)*std::cos(theta*t)/2.0               
                      + M_PI*w*x*x*std::cos(theta*t)*std::cos(M_PI*w*x/2.0)*std::cos(M_PI*w*y)   
                      + 6.0*x*std::sin(M_PI*w*y)*std::cos(theta*t)*std::cos(M_PI*w*x/2.0);
                      
                      sigma(0,1) = - M_PI*w*x*x*std::sin(M_PI*w*x/2.0)*std::sin(M_PI*w*y)*std::cos(theta*t)/2.0 
                      + M_PI*w*x*x*std::cos(theta*t)*std::cos(M_PI*w*x/2.0)*std::cos(M_PI*w*y)                     
                      + 2.0*x*std::sin(M_PI*w*y)*std::cos(theta*t)*std::cos(M_PI*w*x/2.0);
                      
                      sigma(1,0) = sigma(0,1);
                      
                      sigma(1,1) = - M_PI*w*x*x*std::sin(M_PI*w*x/2.0)*std::sin(M_PI*w*y)*std::cos(theta*t)/2.0 
                      + 3*M_PI*w*x*x*std::cos(theta*t)*std::cos(M_PI*w*x/2.0)*std::cos(M_PI*w*y) 
                      + 2.0*x*std::sin(M_PI*w*y)*std::cos(theta*t)*std::cos(M_PI*w*x/2.0);
                      return sigma;
                  };
              }
              break;
              
              case EFunctionInterfaceCoupling: {
                  // sigma(u) = 2 mu eps(u) + lambda div(u) I, mu = lambda = 1
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_matrix<double,2,2> {
                      double x = pt.x(), y = pt.y();
                      double T = std::sqrt(2.0)*M_PI*std::cos(std::sqrt(2.0)*M_PI*t);
                      double dxux = (std::cos(M_PI*x)/3.0 - M_PI*std::sin(M_PI*x))*std::sin(M_PI*y)*T;
                      double dyuy = M_PI*M_PI*x*std::sin(M_PI*y)*T;
                      disk::static_matrix<double,2,2> sigma = disk::static_matrix<double,2,2>::Zero(2,2);
                      sigma(0,0) = 3.0*dxux + dyuy;
                      sigma(1,1) = dxux + 3.0*dyuy;
                      sigma(0,1) = (std::sin(M_PI*x)/3.0 + M_PI*std::cos(M_PI*x) - M_PI)*std::cos(M_PI*y)*T;
                      sigma(1,0) = sigma(0,1);
                      return sigma;
                  };
              }
              break;
              
              case reproduction_acoustic: case EFunctionAcousticQuadraticInSpace: case EFunctionAcousticNonPolynomial:
              case EFunctionAcousticQuadraticInSpaceCos: case EFunctionAcousticNonPolynomialCos: {
                  // no elastic part
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_matrix<double,2,2> {
                      return disk::static_matrix<double,2,2>::Zero();
                  };
              }
              break;
              
              case EFunctionElasticQuadraticInSpace: {
                  // sigma(u) = 2 mu eps(u) + lambda div(u) I, mu = lambda = 1
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_matrix<double,2,2> {
                      double x = pt.x(), y = pt.y();
                      double gx = (1.0-2.0*x)*y*(1.0-y), gy = x*(1.0-x)*(1.0-2.0*y);
                      double s11 = std::sin(std::sqrt(2.0)*M_PI*t)*(3.0*gx + gy), s22 = std::sin(std::sqrt(2.0)*M_PI*t)*(gx + 3.0*gy), s12 = std::sin(std::sqrt(2.0)*M_PI*t)*(gx + gy);
                      disk::static_matrix<double,2,2> sigma;
                      sigma << s11, s12, s12, s22;
                      return sigma;
                  };
              }
              break;
              
              case EFunctionElasticNonPolynomial: {
                  // sigma(u) = 2 mu eps(u) (div u = 0), mu = 1
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_matrix<double,2,2> {
                      double x = pt.x(), y = pt.y();
                      double s11 = -2.0*M_PI*std::sin(std::sqrt(2.0)*M_PI*t)*std::cos(M_PI*x)*std::cos(M_PI*y), s22 = -s11, s12 = 0.0;
                      disk::static_matrix<double,2,2> sigma;
                      sigma << s11, s12, s12, s22;
                      return sigma;
                  };
              }
              break;
              
              default: {
                  std::cout << " Function not implemented " << std::endl;
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_matrix<double,2,2> {
                      disk::static_matrix<double,2,2> sigma(2,2);
                      return sigma;
                  };
              }
              break;
          }
          
      }
      
      // ACOUSTIC
      std::function<double
      (const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_s_u(double & t) {
          
          switch (m_function_type) {
              
              case EFunctionNonPolynomial: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return x*x*std::sin(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;
              
              case EFunctionQuadraticInTime: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return t*t*x*std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;
              
              case EFunctionQuadraticInSpace: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return (1 - x)*x*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t);
                  };
              }
              break;

              case EFunctionQuadraticInSpaceAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return y*(1-x*x)*(1-y)*std::sin(M_PI*t);
                  };
              }
              break;
              
              case reproduction_acoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return t*t*std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;
              
              case EFunctionNonPolynomial_paper: {
                  return [&t, w0 = m_paper_w, th0 = m_paper_theta](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y,w,theta;
                      x = pt.x();
                      y = pt.y();
                      
                      w = w0; 
                      theta = th0; 
                      // w = 5.0;
                      // theta = std::sqrt(2)*M_PI;
                      
                      return x*x*std::sin(theta*t)*std::sin(w*M_PI*x)*std::sin(w*M_PI*y);
                  };
              }
              
              case EFunctionQuarticInTimeAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return t*t*t*t*std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;
              
              case EFunctionCubicInTimeAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return t*t*t*std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;

              case EFunctionPlaneWaveAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y, theta;
                      x = pt.x();
                      y = pt.y();
                      theta = x + y - std::sqrt(2)*t;
                      return -(1.0/std::sqrt(2)) * std::sin(theta);
                  };
              }
              break;
              
              case EFunctionInterfaceCoupling: {
                  // potential
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x = pt.x(), y = pt.y();
                      return (std::cos(M_PI*x) - 2.0*M_PI*std::sin(M_PI*x))*std::sin(M_PI*y)*std::sin(std::sqrt(2.0)*M_PI*t);
                  };
              }
              break;
              
              case reproduction_elastic: case EFunctionElasticQuadraticInSpace: case EFunctionElasticNonPolynomial: {
                  // no acoustic part
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      return 0.0;
                  };
              }
              break;
              
              case EFunctionAcousticQuadraticInSpace: {
                  // p (Burman-Duran-Ern, CAMC 2022, eq. (61))
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x = pt.x(), y = pt.y();
                      return std::sin(std::sqrt(2.0)*M_PI*t)*x*(1.0-x)*y*(1.0-y);
                  };
              }
              break;
              
              case EFunctionAcousticNonPolynomial: {
                  // p (Burman-Duran-Ern, CAMC 2022, eq. (62))
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x = pt.x(), y = pt.y();
                      return std::sin(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*x)*std::sin(M_PI*y)/(std::sqrt(2.0)*M_PI);
                  };
              }
              break;
              
              case EFunctionAcousticQuadraticInSpaceCos: {
                  // p (eq. (61) with cos)
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x = pt.x(), y = pt.y();
                      return std::cos(std::sqrt(2.0)*M_PI*t)*x*(1.0-x)*y*(1.0-y);
                  };
              }
              break;
              
              case EFunctionAcousticNonPolynomialCos: {
                  // p (eq. (62) with cos)
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x = pt.x(), y = pt.y();
                      return std::cos(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*x)*std::sin(M_PI*y)/(std::sqrt(2.0)*M_PI);
                  };
              }
              break;
              
              default: {
                  std::cout << " Function not implemented " << std::endl;
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      return 0;
                  };
              }
              break;
          }
      }
      
      std::function<double
      (const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_s_v(double & t){
                
          switch (m_function_type) {
              
              case EFunctionNonPolynomial: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return std::sqrt(2.0)*M_PI*x*x*std::cos(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;
              
              case EFunctionQuadraticInTime: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return 2*t*x*std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;
                    
              case EFunctionQuadraticInSpace: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return std::sqrt(2.0)*M_PI*(1 - x)*x*x*(1 - y)*y*std::cos(std::sqrt(2.0)*M_PI*t);
                  };
              }
              break;
              
              case EFunctionQuadraticInSpaceAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return M_PI*y*(1-x*x)*(1-y)*std::cos(M_PI*t);
                  };
              }
              break;

              case reproduction_acoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      return 2.0*t*std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
                  };
              }
              break;
              
              case EFunctionNonPolynomial_paper: {
                  return [&t, w0 = m_paper_w, th0 = m_paper_theta](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y,w,theta;
                      x = pt.x();
                      y = pt.y();
                      w = w0; 
                      theta = th0; 
                      // w = 5.0;
                      // theta = std::sqrt(2)*M_PI;
                      return theta*x*x*std::sin(M_PI*w*x)*std::sin(M_PI*w*y)*std::cos(theta*t);
                  };
              }
              break;
              
              case EFunctionQuarticInTimeAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return 4.0*t*t*t * std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;
              
              case EFunctionCubicInTimeAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return 3.0*t*t * std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;

              case EFunctionPlaneWaveAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y, theta;
                      x = pt.x();
                      y = pt.y();
                      theta = x + y - std::sqrt(2)*t;
                      return std::cos(theta);
                  };
              }
              break;
              
              case EFunctionInterfaceCoupling: {
                  // time derivative of the potential
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x = pt.x(), y = pt.y();
                      return (std::cos(M_PI*x) - 2.0*M_PI*std::sin(M_PI*x))*std::sin(M_PI*y)*std::sqrt(2.0)*M_PI*std::cos(std::sqrt(2.0)*M_PI*t);
                  };
              }
              break;
              
              case reproduction_elastic: case EFunctionElasticQuadraticInSpace: case EFunctionElasticNonPolynomial: {
                  // no acoustic part
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      return 0.0;
                  };
              }
              break;
              
              case EFunctionAcousticQuadraticInSpace: {
                  // p_t
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x = pt.x(), y = pt.y();
                      return std::sqrt(2.0)*M_PI*std::cos(std::sqrt(2.0)*M_PI*t)*x*(1.0-x)*y*(1.0-y);
                  };
              }
              break;
              
              case EFunctionAcousticNonPolynomial: {
                  // p_t
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x = pt.x(), y = pt.y();
                      return std::cos(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;
              
              case EFunctionAcousticQuadraticInSpaceCos: {
                  // p_t (eq. (61) with cos)
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x = pt.x(), y = pt.y();
                      return -std::sqrt(2.0)*M_PI*std::sin(std::sqrt(2.0)*M_PI*t)*x*(1.0-x)*y*(1.0-y);
                  };
              }
              break;
              
              case EFunctionAcousticNonPolynomialCos: {
                  // p_t (eq. (62) with cos)
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x = pt.x(), y = pt.y();
                      return -std::sqrt(2.0)*M_PI*std::sin(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*x)*std::sin(M_PI*y)/(std::sqrt(2.0)*M_PI);
                  };
              }
              break;
              
              default: {
                  std::cout << " Function not implemented " << std::endl;
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      return 0;
                  };
              }
              break;
          }
          
      }
      
      std::function<double
      (const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_s_a(double & t){
          
          switch (m_function_type) {
              case EFunctionNonPolynomial: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return -2*M_PI*M_PI*x*x*std::sin(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;
              
              case EFunctionQuadraticInTime: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return 2*x*std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;
                    
              case EFunctionQuadraticInSpace: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return -2*M_PI*M_PI*(1 - x)*x*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t);
                  };
              }
              break;
              
              case EFunctionQuadraticInSpaceAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return -M_PI*M_PI*y*(1-x*x)*(1-y)*std::sin(M_PI*t);
                  };
              }
              break;

              case reproduction_acoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return 2.0*std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
                  };
              }
              break;
                    
              case EFunctionNonPolynomial_paper: {
                  return [&t, w0 = m_paper_w, th0 = m_paper_theta](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y,w,theta;
                      x = pt.x();
                      y = pt.y();
                      w = w0; 
                      theta = th0; 
                      // w = 5.0;
                      // theta = std::sqrt(2)*M_PI;
                      return -theta*theta*x*x*std::sin(theta*t)*std::sin(M_PI*w*x)*std::sin(M_PI*w*y);
                  };
              }
              break;
              
              case EFunctionQuarticInTimeAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return 12*t*t * std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;
              
              case EFunctionCubicInTimeAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y;
                      x = pt.x();
                      y = pt.y();
                      return 6*t * std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;

              case EFunctionPlaneWaveAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y, theta;
                      x = pt.x();
                      y = pt.y();
                      theta = x + y - std::sqrt(2)*t;
                      return std::sqrt(2) * std::sin(theta);
                  };
              }
              break;
              
              case EFunctionInterfaceCoupling: {
                  // second time derivative of the potential
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x = pt.x(), y = pt.y();
                      return (std::cos(M_PI*x) - 2.0*M_PI*std::sin(M_PI*x))*std::sin(M_PI*y)*-2.0*M_PI*M_PI*std::sin(std::sqrt(2.0)*M_PI*t);
                  };
              }
              break;
              
              case reproduction_elastic: case EFunctionElasticQuadraticInSpace: case EFunctionElasticNonPolynomial: {
                  // no acoustic part
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      return 0.0;
                  };
              }
              break;
              
              case EFunctionAcousticQuadraticInSpace: {
                  // p_tt
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x = pt.x(), y = pt.y();
                      return -2.0*M_PI*M_PI*std::sin(std::sqrt(2.0)*M_PI*t)*x*(1.0-x)*y*(1.0-y);
                  };
              }
              break;
              
              case EFunctionAcousticNonPolynomial: {
                  // p_tt
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x = pt.x(), y = pt.y();
                      return -std::sqrt(2.0)*M_PI*std::sin(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*x)*std::sin(M_PI*y);
                  };
              }
              break;
              
              case EFunctionAcousticQuadraticInSpaceCos: {
                  // p_tt (eq. (61) with cos)
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x = pt.x(), y = pt.y();
                      return -2.0*M_PI*M_PI*std::cos(std::sqrt(2.0)*M_PI*t)*x*(1.0-x)*y*(1.0-y);
                  };
              }
              break;
              
              case EFunctionAcousticNonPolynomialCos: {
                  // p_tt (eq. (62) with cos)
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x = pt.x(), y = pt.y();
                      return -2.0*M_PI*M_PI*std::cos(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*x)*std::sin(M_PI*y)/(std::sqrt(2.0)*M_PI);
                  };
              }
              break;
              
              default: {
                  std::cout << " Function not implemented " << std::endl;
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      return 0;
                  };
              }
              break;
          }    
      }
      
      std::function<double
      (const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_s_f(double & t){
          
          switch (m_function_type) {
              case EFunctionNonPolynomial: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y,f;
                      x = pt.x();
                      y = pt.y();
                      f = -2*std::sin(std::sqrt(2.0)*M_PI*t)*(2*M_PI*x*std::cos(M_PI*x) + std::sin(M_PI*x))*std::sin(M_PI*y);
                      return f;
                  };
              }
              break;
              
              case EFunctionQuadraticInTime: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y,f;
                      x = pt.x();
                      y = pt.y();
                      f = 2*(-(M_PI*t*t*std::cos(M_PI*x)) + (1 + M_PI*M_PI*t*t)*x*std::sin(M_PI*x))*std::sin(M_PI*y);
                      return f;
                  };
              }
              break;
              
              case EFunctionQuadraticInSpace: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y,f;
                      x = pt.x();
                      y = pt.y();
                      f = 2*((-1 + y)*y - 3*x*(-1 + y)*y + x*x*x*(-1 - M_PI*M_PI*(-1 + y)*y)
                      + x*x*(1 + M_PI*M_PI*(-1 + y)*y))*std::sin(std::sqrt(2.0)*M_PI*t);
                      return f;
                  };
              }
              break;

              case EFunctionQuadraticInSpaceAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y,f;
                      x = pt.x();
                      y = pt.y();
                      f = - M_PI*M_PI*y*(1-x*x)*(1-y)*std::sin(M_PI*t)
                          + 2*y*(1-y)*std::sin(M_PI*t)
                          + 2*(1-x*x)*std::sin(M_PI*t);
                      return f;
                  };
              }
              break;
              
              case reproduction_acoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y,f;
                      x = pt.x();
                      y = pt.y();
                      return 2.0*(1.0 + M_PI*M_PI*t*t)*std::sin(M_PI*pt.x())*std::sin(M_PI*pt.y());
                  };
              }
              break;
              
              case EFunctionNonPolynomial_paper: {
                  return [&t, w0 = m_paper_w, th0 = m_paper_theta](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x,y,f,w,theta;
                      x = pt.x();
                      y = pt.y();
                      w = w0; 
                      theta = th0; 
                      // w = 5.0;
                      // theta = std::sqrt(2)*M_PI;
                      f = - theta*theta*x*x*std::sin(theta*t)*std::sin(M_PI*w*x)*std::sin(M_PI*w*y) 
                      + 2*M_PI*M_PI*theta*w*w*x*x*std::sin(theta*t)*std::sin(M_PI*w*x)*std::sin(M_PI*w*y)/theta 
                      - 4*M_PI*theta*w*x*std::sin(theta*t)*std::sin(M_PI*w*y)*std::cos(M_PI*w*x)/theta 
                      - 2*theta*std::sin(theta*t)*std::sin(M_PI*w*x)*std::sin(M_PI*w*y)/theta;
                      return f;
                  };
              }
              break;
              
              case EFunctionQuarticInTimeAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y,f;
                      x = pt.x();
                      y = pt.y();
                      f = 2*M_PI*M_PI*t*t*t*t*std::sin(M_PI*x)*std::sin(M_PI*y) + 12*t*t*std::sin(M_PI*x)*std::sin(M_PI*y);
                      return f;
                  };
              }
              break;
              
              case EFunctionCubicInTimeAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y,f;
                      x = pt.x();
                      y = pt.y();
                      f = 2*M_PI*M_PI*t*t*t*std::sin(M_PI*x)*std::sin(M_PI*y) + 6*t*std::sin(M_PI*x)*std::sin(M_PI*y);
                      return f;
                  };
              }
              break;
              
              case EFunctionPlaneWaveAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double { 
                      double x,y,f;
                      x = pt.x();
                      y = pt.y();
                      f = 0.0;
                      return f;
                  };
              }
              break;

              case EFunctionInterfaceCoupling: {
                  // d_tt phi - Delta phi = 0
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      return 0.0;
                  };
              }
              break;
              
              case reproduction_elastic: case EFunctionElasticQuadraticInSpace: case EFunctionElasticNonPolynomial: {
                  // no acoustic part
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      return 0.0;
                  };
              }
              break;
              
              case EFunctionAcousticQuadraticInSpace: {
                  // p_tt - Delta p (rho = 1, c = 1)
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x = pt.x(), y = pt.y();
                      return std::sin(std::sqrt(2.0)*M_PI*t)*(-2.0*M_PI*M_PI*x*(1.0-x)*y*(1.0-y) + 2.0*(x*(1.0-x) + y*(1.0-y)));
                  };
              }
              break;
              
              case EFunctionAcousticNonPolynomial: {
                  // p_tt - Delta p = 0
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x = pt.x(), y = pt.y();
                      return 0.0*x*y*t;
                  };
              }
              break;
              
              case EFunctionAcousticQuadraticInSpaceCos: {
                  // p_tt - Delta p (rho = 1, c = 1) (eq. (61) with cos)
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x = pt.x(), y = pt.y();
                      return std::cos(std::sqrt(2.0)*M_PI*t)*(-2.0*M_PI*M_PI*x*(1.0-x)*y*(1.0-y) + 2.0*(x*(1.0-x) + y*(1.0-y)));
                  };
              }
              break;
              
              case EFunctionAcousticNonPolynomialCos: {
                  // p_tt - Delta p (rho = 1, c = 1) (eq. (62) with cos)
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      double x = pt.x(), y = pt.y();
                      return 0.0*x*y*t;
                  };
              }
              break;
              
              default: {
                  std::cout << " Function not implemented " << std::endl;
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> double {
                      return 0;
                  };
              }
              break;
          }
      }
      
      std::function<disk::static_vector<double, 2>
      (const typename disk::generic_mesh<double, 2>::point_type& )> Evaluate_s_q(double & t){
          
          switch (m_function_type) {
              case EFunctionNonPolynomial: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,qx,qy;
                      x = pt.x();
                      y = pt.y();
                      qx = M_PI*x*x*std::cos(M_PI*x)*std::sin(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*y) + 2*x*std::sin(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*x)*std::sin(M_PI*y);
                      qy = M_PI*x*x*std::cos(M_PI*y)*std::sin(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*x);
                      disk::static_vector<double, 2> q{qx,qy};
                      return q;
                  };
              }
              break;
              
              case EFunctionQuadraticInTime: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,qx,qy;
                      x = pt.x();
                      y = pt.y();
                      std::vector<double> flux(2);
                      qx = M_PI*t*t*x*std::cos(M_PI*x)*std::sin(M_PI*y) + t*t*std::sin(M_PI*x)*std::sin(M_PI*y);
                      qy = M_PI*t*t*x*std::cos(M_PI*y)*std::sin(M_PI*x);
                      disk::static_vector<double, 2> q{qx,qy};
                      return q;
                  };
              }
              break;
              
              case EFunctionQuadraticInSpace: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,qx,qy;
                      x = pt.x();
                      y = pt.y();
                      std::vector<double> flux(2);
                      qx = 2*(1 - x)*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t)
                           - x*x*(1 - y)*y*std::sin(std::sqrt(2.0)*M_PI*t);
                      qy = (1 - x)*x*x*(1 - y)*std::sin(std::sqrt(2.0)*M_PI*t)
                           - (1 - x)*x*x*y*std::sin(std::sqrt(2.0)*M_PI*t);
                      disk::static_vector<double, 2> q{qx,qy};
                      return q;
                  };
              }
              break;

              case EFunctionQuadraticInSpaceAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,qx,qy;
                      x = pt.x();
                      y = pt.y();
                      std::vector<double> flux(2);
                      qx = -2*x*y*(1-y)*std::sin(M_PI*t);
                      qy = - y*(1-x*x)*std::sin(M_PI*t) 
                           + (1-x*x)*(1-y)*std::sin(M_PI*t);
                      disk::static_vector<double, 2> q{qx,qy};
                      return q;
                  };
              }
              break;
              
              case reproduction_acoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x,y,qx,qy;
                      x = pt.x();
                      y = pt.y();
                      qx = M_PI*t*t*std::cos(M_PI*x)*std::sin(M_PI*y);
                      qy = M_PI*t*t*std::sin(M_PI*x)*std::cos(M_PI*y);
                      disk::static_vector<double, 2> q{qx,qy};
                      return q;
                  };
              }
              break;
              
              case EFunctionNonPolynomial_paper: {
                  return [&t, w0 = m_paper_w, th0 = m_paper_theta](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2> {
                      double x,y,qx,qy,w,theta;
                      x = pt.x();
                      y = pt.y();
                      w = w0; 
                      theta = th0; 
                      // w = 5.0;
                      // theta = std::sqrt(2)*M_PI;
                      qx = M_PI*theta*w*x*x*std::sin(theta*t)*std::sin(M_PI*w*y)*std::cos(M_PI*w*x)/theta + 2*theta*x*std::sin(theta*t)*std::sin(M_PI*w*x)*std::sin(M_PI*w*y)/theta;
                      qy = M_PI*theta*w*x*x*std::sin(theta*t)*std::sin(M_PI*w*x)*std::cos(M_PI*w*y)/theta;
                      disk::static_vector<double, 2> q{qx,qy};
                      return q;
                  };
              }
              break;
              
              case EFunctionQuarticInTimeAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x,y,qx,qy;
                      x = pt.x();
                      y = pt.y();
                      qx = M_PI*t*t*t*t*std::cos(M_PI*x)*std::sin(M_PI*y);
                      qy = M_PI*t*t*t*t*std::sin(M_PI*x)*std::cos(M_PI*y);
                      disk::static_vector<double, 2> q{qx,qy};
                      return q;
                  };
              }
              break;
              
              case EFunctionCubicInTimeAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x,y,qx,qy;
                      x = pt.x();
                      y = pt.y();
                      qx = M_PI*t*t*t*std::cos(M_PI*x)*std::sin(M_PI*y);
                      qy = M_PI*t*t*t*std::sin(M_PI*x)*std::cos(M_PI*y);
                      disk::static_vector<double, 2> q{qx,qy};
                      return q;
                  };
              }
              break;
              
              case EFunctionPlaneWaveAcoustic: {
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x, y, qx, qy, theta;
                      x = pt.x();
                      y = pt.y();
                      theta = x + y - std::sqrt(2)*t;
                      qx = (1.0/std::sqrt(2.0)) * std::sin(theta);
                      qy = (1.0/std::sqrt(2.0)) * std::sin(theta);
                      disk::static_vector<double, 2> q{qx,qy};
                      return q;
                  };
              }
              break;

              case EFunctionInterfaceCoupling: {
                  // grad phi
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x = pt.x(), y = pt.y();
                      double T = std::sin(std::sqrt(2.0)*M_PI*t);
                      double qx = (-M_PI*std::sin(M_PI*x) - 2.0*M_PI*M_PI*std::cos(M_PI*x))*std::sin(M_PI*y)*T;
                      double qy = M_PI*(std::cos(M_PI*x) - 2.0*M_PI*std::sin(M_PI*x))*std::cos(M_PI*y)*T;
                      disk::static_vector<double, 2> q{qx,qy};
                      return q;
                  };
              }
              break;
              
              case reproduction_elastic: case EFunctionElasticQuadraticInSpace: case EFunctionElasticNonPolynomial: {
                  // no acoustic part
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      disk::static_vector<double, 2> z{0.0, 0.0};
                      return z;
                  };
              }
              break;
              
              case EFunctionAcousticQuadraticInSpace: {
                  // grad p
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x = pt.x(), y = pt.y();
                      double rx = std::sin(std::sqrt(2.0)*M_PI*t)*(1.0-2.0*x)*y*(1.0-y), ry = std::sin(std::sqrt(2.0)*M_PI*t)*x*(1.0-x)*(1.0-2.0*y);
                      disk::static_vector<double, 2> r{rx, ry};
                      return r;
                  };
              }
              break;
              
              case EFunctionAcousticNonPolynomial: {
                  // grad p
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x = pt.x(), y = pt.y();
                      double rx = std::sin(std::sqrt(2.0)*M_PI*t)*std::cos(M_PI*x)*std::sin(M_PI*y)/std::sqrt(2.0), ry = std::sin(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*x)*std::cos(M_PI*y)/std::sqrt(2.0);
                      disk::static_vector<double, 2> r{rx, ry};
                      return r;
                  };
              }
              break;
              
              case EFunctionAcousticQuadraticInSpaceCos: {
                  // grad p
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x = pt.x(), y = pt.y();
                      double rx = std::cos(std::sqrt(2.0)*M_PI*t)*(1.0-2.0*x)*y*(1.0-y), ry = std::cos(std::sqrt(2.0)*M_PI*t)*x*(1.0-x)*(1.0-2.0*y);
                      disk::static_vector<double, 2> r{rx, ry};
                      return r;
                  };
              }
              break;
              
              case EFunctionAcousticNonPolynomialCos: {
                  // grad p
                  return [&t](const typename disk::generic_mesh<double, 2>::point_type& pt) -> disk::static_vector<double, 2> {
                      double x = pt.x(), y = pt.y();
                      double rx = std::cos(std::sqrt(2.0)*M_PI*t)*std::cos(M_PI*x)*std::sin(M_PI*y)/std::sqrt(2.0), ry = std::cos(std::sqrt(2.0)*M_PI*t)*std::sin(M_PI*x)*std::cos(M_PI*y)/std::sqrt(2.0);
                      disk::static_vector<double, 2> r{rx, ry};
                      return r;
                  };
              }
              break;
              
              default: {
                  std::cout << " Function not implemented " << std::endl;
                  return [](const typename disk::generic_mesh<double, 2>::point_type& pt)
                  -> disk::static_vector<double, 2>{
                      disk::static_vector<double, 2> f{0,0};
                      return f;
                  };
              }
              break;
          }    
      }
      
      private:
      
      EFunctionType m_function_type;
      double m_paper_w = 1.0;               // EFunctionNonPolynomial_paper: omega and theta pi (set_paper_frequencies)
      double m_paper_theta = 10.0 * M_PI;
      
  };
  
  #endif /* scal_vec_analytic_functions_hpp */
  