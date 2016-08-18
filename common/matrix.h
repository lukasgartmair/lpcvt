/*
 *  Copyright (c) 2010, INRIA, Project ALICE
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://alice.loria.fr
 *
 *     ALICE Project
 *     INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#ifndef __GEEX_MATHEMATICS_MATRIX__
#define __GEEX_MATHEMATICS_MATRIX__

#include <iostream>

namespace Geex {

//______________________________________________________________________


    template <class FT, int DIM> class Matrix ;
    template <class FT, int DIM> Matrix<FT, DIM> operator*(
        FT op1, const Matrix<FT, DIM>& op2
    ) ;


  /**
    * A class for representing dense matrices.
    */
    template <class FT, int DIM> class Matrix {
    public:
        Matrix() ;
        void load_zero() ;
        void load_identity() ;

        FT& operator()(int i, int j) ;
        const FT& operator()(int i, int j) const ;

        Matrix<FT, DIM>& operator+=(const Matrix<FT, DIM>& rhs) ;
        Matrix<FT, DIM>& operator-=(const Matrix<FT, DIM>& rhs) ;
        Matrix<FT, DIM>& operator*=(FT rhs) ;
        Matrix<FT, DIM>& operator/=(FT rhs) ;

        Matrix<FT, DIM> operator+(const Matrix<FT, DIM>& op2) const ;
        Matrix<FT, DIM> operator-(const Matrix<FT, DIM>& op2) const ;
        Matrix<FT, DIM> operator*(const Matrix<FT, DIM>& op2) const ;

        Matrix<FT, DIM> inverse() const ;
        Matrix<FT, DIM> transpose() const ;

        void mult(const FT* x, FT* y) const ;

        // Routines for interfacing with Fortran, OpenGL etc...

        const FT* data() const { return &(coeff_[0][0]) ; }
        FT* data() { return &(coeff_[0][0]) ; }

        void get_lower_triangle(FT* store) {
            for(unsigned int i=0; i<DIM; i++) {
                for(unsigned int j=0; j<=i; j++) {
                    *store++ = coeff_[i][j] ;
                }
            }
        }

    private:
        FT coeff_[DIM][DIM] ;
    } ;

//_______________________________________________________________________

    template <class FT, int DIM> inline 
    Matrix<FT, DIM>::Matrix() {
        load_identity() ;
    }


    template <class FT, int DIM> inline 
    FT& Matrix<FT, DIM>::operator()(int i, int j) {
        return coeff_[i][j] ;
    }

    template <class FT, int DIM> inline 
    const FT& Matrix<FT, DIM>::operator()(int i, int j) const {
        return coeff_[i][j] ;
    }

    template <class FT, int DIM> inline 
    Matrix<FT, DIM>& Matrix<FT, DIM>::operator+=(const Matrix<FT, DIM>& rhs) {
        for(int i=0; i<DIM; i++) {
            for(int j=0; j<DIM; j++) {
                coeff_[i][j] += rhs.coeff_[i][j] ;
            }
        }
        return *this ;
    }

    template <class FT, int DIM> inline 
    Matrix<FT, DIM>& Matrix<FT, DIM>::operator-=(const Matrix<FT, DIM>& rhs) {
        for(int i=0; i<DIM; i++) {
            for(int j=0; j<DIM; j++) {
                coeff_[i][j] -= rhs.coeff_[i][j] ;
            }
        }
        return *this ;
    }

    template <class FT, int DIM> inline 
    Matrix<FT, DIM>& Matrix<FT, DIM>::operator*=(FT rhs) {
        for(int i=0; i<DIM; i++) {
            for(int j=0; j<DIM; j++) {
                coeff_[i][j] *= rhs ;
            }
        }
        return *this ;
    }

    template <class FT, int DIM> inline 
    Matrix<FT, DIM>& Matrix<FT, DIM>::operator/=(FT rhs) {
        for(int i=0; i<DIM; i++) {
            for(int j=0; j<DIM; j++) {
                coeff_[i][j] /= rhs ;
            }
        }
        return *this ;
    }

    template <class FT, int DIM> inline 
    Matrix<FT, DIM> Matrix<FT, DIM>::operator+(const Matrix<FT, DIM>& op2) const {
        Matrix<FT, DIM> result = *this ;
        result += op2 ;
        return result ;
    }

    template <class FT, int DIM> inline 
    Matrix<FT, DIM> Matrix<FT, DIM>::operator-(const Matrix<FT, DIM>& op2) const {
        Matrix<FT, DIM> result = *this ;
        result -= op2 ;
        return result ;
    
    }

    template <class FT, int DIM> inline 
    Matrix<FT, DIM> Matrix<FT, DIM>::operator*(const Matrix<FT, DIM>& op2) const {
        Matrix<FT, DIM> result ;
        result.load_zero() ;
        for(int i=0; i<DIM; i++) {
            for(int j=0; j<DIM; j++) {
                for(int k=0; k<DIM; k++) {
                    result.coeff_[i][j] += coeff_[i][k] * op2.coeff_[k][j] ;
                }
            }
        }
        return result ;
    }

    template <class FT, int DIM> inline 
    std::ostream& operator << (std::ostream& output, const Matrix<FT, DIM>& m) {
        for(int i=0; i<DIM; i++) {
            for(int j=0; j<DIM; j++) {
                output << m(i,j) << " " ;
            }
        }    
        return output ;
    }

    template <class FT, int DIM> inline 
    std::istream& operator >> (std::istream& input, Matrix<FT, DIM>& m) {
        for(int i=0; i<DIM; i++) {
            for(int j=0; j<DIM; j++) {
                input >> m(i,j) ;
            }
        }    
        return input ;
    }



//_______________________________________________________________________

    template <class FT, int DIM> void
    Matrix<FT, DIM>::load_zero() {
        for(int i=0; i<DIM; i++) {
            for(int j=0; j<DIM; j++) {
                coeff_[i][j] = FT(0) ;
            }
        }
    }

    template <class FT, int DIM> void
    Matrix<FT, DIM>::load_identity() {
        for(int i=0; i<DIM; i++) {
            for(int j=0; j<DIM; j++) {
                coeff_[i][j] = (i==j) ? FT(1) : FT(0) ;
            }
        }
    }

    template <class FT, int DIM> Matrix<FT, DIM> 
    Matrix<FT, DIM>::inverse() const {
        FT val, val2;
        int i, j, k, ind;
        Matrix<FT, DIM> tmp = (*this);
        Matrix<FT, DIM> result;
    
        result.load_identity();
    
        for (i = 0; i != DIM; i++) {
            val = tmp(i,i);			/* find pivot */
            ind = i;
            for (j = i + 1; j != DIM; j++) {
                if (fabs(tmp(j,i)) > fabs(val)) {
                    ind = j;
                    val = tmp(j,i);
                }
            }
            
            if (ind != i) {			
                for (j = 0; j != DIM; j++) {
                    val2 = result(i,j);
                    result(i,j) = result(ind,j);
                    result(ind,j) = val2;           /* swap columns */
                    val2 = tmp(i,j);
                    tmp(i,j) = tmp(ind,j);
                    tmp(ind,j) = val2;
                }
            }
            
            for (j = 0; j != DIM; j++) {
                tmp(i,j)    /= val;
                result(i,j) /= val;
            }
        
            for (j = 0; j != DIM; j++) {		
                if (j == i)
                    continue;                       /* eliminate column */
                val = tmp(j,i);
                for (k = 0; k != DIM; k++) {
                    tmp(j,k)     -= tmp(i,k)     * val;
                    result(j,k)  -= result(i,k)  * val;
                }
            }
        }
    
        return result;
    }

    template <class FT, int DIM> Matrix<FT, DIM> 
    Matrix<FT, DIM>::transpose() const {
        Matrix<FT, DIM> result ;
        for(int i=0; i<DIM; i++) {
            for(int j=0; j<DIM; j++) {
                result(i,j) = (*this)(j,i) ;
            }
        }
        return result ;
    }

    template <class FT, int N> inline void
    Matrix<FT,N>::mult(const FT* x, FT* y) const {
        for(int i=0; i<N; i++) {
            y[i] = 0 ;
            for(int j=0; j<N; j++) {
                y[i] += (*this)(i,j)*x[j] ;
            }
        }
    }

//_______________________________________________________________________

}

#endif
