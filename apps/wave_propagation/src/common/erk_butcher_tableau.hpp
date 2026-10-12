//
//  erk_butcher_tableau.hpp
//  acoustics
//
//  Created by Omar Durán on 5/8/20.
//

#pragma once
#ifndef erk_butcher_tableau_hpp
#define erk_butcher_tableau_hpp

#include <Eigen/Dense>

class erk_butcher_tableau
{
    public:
    
    static void erk_tables(int s, Matrix<double, Dynamic, Dynamic> &a, Matrix<double, Dynamic, 1> &b, Matrix<double, Dynamic, 1> &c){
        
        // Mathematical Aspects of Discontinuous Galerkin Methods
        // Authors: Di Pietro, Daniele Antonio, Ern, Alexandre
        a = Matrix<double, Dynamic, Dynamic>::Zero(s, s);
        b = Matrix<double, Dynamic, 1>::Zero(s, 1);
        c = Matrix<double, Dynamic, 1>::Zero(s, 1);
        
        switch (s) {
            case 1:
                {
                    a(0,0) = 0.0;
                    b(0,0) = 1.0;
                    c(0,0) = 0.0;
                }
                break;
            case 2:
                {
//                    a(0,0) = 0.0;
//                    a(1,0) = 1.0;
//                    a(1,1) = 0.0;
//
//                    b(0,0) = 0.5;
//                    b(1,0) = 0.5;
//
//                    c(0,0) = 0.0;
//                    c(1,0) = 1.0;
                    
                    a(0,0) = 0.0;
                    a(1,0) = 0.5;
                    a(1,1) = 0.0;

                    b(0,0) = 0.0;
                    b(1,0) = 1.0;

                    c(0,0) = 0.0;
                    c(1,0) = 0.5;
                    
                }
                break;
            case 3:
                {

//                    a(0,0) = 0.0;
//                    a(1,0) = 1.0/3.0;
//                    a(1,1) = 0.0;
//                    a(2,0) = 0.0;
//                    a(2,1) = 2.0/3.0;
//                    a(2,2) = 0.0;
//
//                    b(0,0) = 1.0/4.0;
//                    b(1,0) = 0.0;
//                    b(2,0) = 3.0/4.0;
//
//                    c(0,0) = 0.0;
//                    c(1,0) = 1.0/3.0;
//                    c(2,0) = 2.0/3.0;
                    
                    a(0,0) = 0.0;
                    a(1,0) = 1.0/2.0;
                    a(1,1) = 0.0;
                    a(2,0) = -1.0;
                    a(2,1) = 2.0;
                    a(2,2) = 0.0;
                    
                    b(0,0) = 1.0/6.0;
                    b(1,0) = 2.0/3.0;
                    b(2,0) = 1.0/6.0;
                    
                    c(0,0) = 0.0;
                    c(1,0) = 1.0/2.0;
                    c(2,0) = 1.0;
                    
                }
                break;
                case 4:
                {

                    a(0,0) = 0.0;
                    a(1,0) = 1.0/2.0;
                    a(2,0) = 0.0;
                    a(3,0) = 0.0;
                    a(1,1) = 0.0;
                    a(2,1) = 1.0/2.0;
                    a(3,1) = 0.0;
                    a(2,2) = 0.0;
                    a(3,2) = 1.0;
                    a(3,3) = 0.0;
                    
                    b(0,0) = 1.0/6.0;
                    b(1,0) = 1.0/3.0;
                    b(2,0) = 1.0/3.0;
                    b(3,0) = 1.0/6.0;
                    
                    c(0,0) = 0.0;
                    c(1,0) = 1.0/2.0;
                    c(2,0) = 1.0/2.0;
                    c(3,0) = 1.0;
                    
                }
                break;
            default:
            {
                std::cout << "Error:: Method not implemented." << std::endl;
            }
                break;
        }
        
    }
    
};

#endif /* erk_butcher_tableau_hpp */
