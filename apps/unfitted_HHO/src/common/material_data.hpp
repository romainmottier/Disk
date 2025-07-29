
#ifndef material_data_hpp
#define material_data_hpp

#include <stdio.h>

template<typename T = double>
class material_data {
        
    T m_kappa_1;
    T m_kappa_2;
    
public:
    
    /// Default constructor
    material_data(T kappa_1, T kappa_2){
        m_kappa_1 = kappa_1;
        m_kappa_2 = kappa_2;
    }
    
    /// Copy constructor
    material_data(const material_data & other){
        m_kappa_1 = other.kappa_1;
        m_kappa_2 = other.kappa_2;
    }
    
    /// Assignement constructor
    const material_data & operator=(const material_data & other){
        
        // check for self-assignment
        if(&other == this){
            return *this;
        }
        
        m_kappa_1 = other.m_kappa_1;
        m_kappa_2 = other.m_kappa_2;

        return *this;
        
    }
    
    /// Desconstructor
    virtual ~material_data(){
        
    }
    
};

#endif /* material_data_hpp */
