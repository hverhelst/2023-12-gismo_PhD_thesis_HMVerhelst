/** @file gsGluingData.h

    @brief Compute the gluing data for one interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Farahat
*/

#pragma once

#include <gsUnstructuredSplines/src/gsC1SurfGD.h>
#include <gsUnstructuredSplines/src/gsC1SurfGluingDataAssembler.h>
#include <gsUnstructuredSplines/src/gsC1SurfGluingDataVisitor.h>
#include <gsUnstructuredSplines/src/gsG1AuxiliaryPatch.h>


namespace gismo
{

template<class T, class Visitor = gsC1SurfGluingDataVisitor<T>>
class gsC1SurfGluingData : public gsC1SurfGD<T>
{

public:
    gsC1SurfGluingData()
    {
        setGDEdge();
    }

    gsC1SurfGluingData(gsMultiPatch<T> const & mp,
    gsMultiBasis<T> & mb)
    :  gsC1SurfGD<T>(mp, mb)
    {
        // Solve the system for alpha_L and alpha_R and beta
        refresh();
        assemble();
        solve();

        // Solve the system for beta_L and beta_R
        refreshBeta();
        assembleBeta();
        solveBeta();
    }

    gsMatrix<T> evalAlpha_R(gsMatrix<T> points)
    {
        gsMatrix<T> ones(1, points.cols());
        ones.setOnes();
        return sol.row(0) * ( ones - points ) + sol.row(1) * points;
    }

    gsMatrix<T> evalAlpha_L(gsMatrix<T> points)
    {
        gsMatrix<T> ones(1, points.cols());
        ones.setOnes();
        return sol.row(2) * ( ones - points ) + sol.row(3) * points;
    }

    gsMatrix<T> evalBeta_R(gsMatrix<T> points)
    {
        gsMatrix<T> ones(1, points.cols());
        ones.setOnes();
        return solBeta.row(0) * ( ones - points ) + solBeta.row(1) * points;
    }

    gsMatrix<T> evalBeta_L(gsMatrix<T> points)
    {
        gsMatrix<T> ones(1, points.cols());
        ones.setOnes();
        return solBeta.row(2) * ( ones - points ) + solBeta.row(3) * points;
    }

    gsMatrix<T> evalBeta(gsMatrix<T> points)
    {
        gsMatrix<T> ones(1, points.cols());
        ones.setOnes();
        return sol.row(4) * ( ones - points ).cwiseProduct( ones - points)
            + sol.row(5) * ( ones - points ).cwiseProduct(points) + sol.row(6) * points;
    }

    gsMatrix<T> getSol(){ return sol; }

    gsMatrix<T> getSolBeta(){ return solBeta; }



protected:

    gsSparseSystem<T> mSys;
    gsSparseSystem<T> mSysBeta;
    gsMatrix<T> dirichletDofs;
    gsMatrix<T> dirichletDofsBeta;

    gsMatrix<T> sol; // In order, it contains: alpha_0L, alpha_1L, alpha_0R, alpha_1R, beta_0, beta_1, beta_2
                    // to construct the linear combination of the GD:
                    // alpha_L = ( 1 - t ) * alpha_0L + alpha_1L * t
                    // alpha_R = ( 1 - t ) * alpha_0R + alpha_1R * t
                    //beta = ( 1 - t )^2 * beta_0 + 2 * t * ( 1 - t ) * beta_1 + t^2 * beta_2

    gsMatrix<T> solBeta;

    void refresh()
    {
        gsVector<index_t> size(1);
        size << 7;

        gsDofMapper map(size);
        map.finalize();

        gsSparseSystem<T> sys(map);
        mSys = sys;
    }

    void refreshBeta()
    {
        gsVector<index_t> size(1);
        size << 4;

        gsDofMapper mapBeta(size);
        mapBeta.finalize();

        gsSparseSystem<T> sysBeta(mapBeta);
        mSysBeta = sysBeta;
    }

    void assemble()
    {
        mSys.reserve(49, 1); // Reserve for the matrix 7x7 values

        dirichletDofs.setZero(mSys.colMapper(0).boundarySize(),1);

        // Assemble volume integrals
        Visitor visitor;
        apply(visitor);

        mSys.matrix().makeCompressed();
    }

    void assembleBeta()
    {
        mSysBeta.reserve(16, 1); // Reserve for the matrix 4x4 values

        dirichletDofsBeta.setZero(mSysBeta.colMapper(0).boundarySize(), 1);

        // Assemble volume integrals
        Visitor visitorBeta;
        applyBeta(visitorBeta);

        mSysBeta.matrix().makeCompressed();
    }

    void apply(Visitor visitor)
    {
    #pragma omp parallel
        {
            Visitor
    #ifdef _OPENMP
            // Create thread-private visitor
            visitor_(visitor);
            const int tid = omp_get_thread_num();
            const int nt  = omp_get_num_threads();
    #else
                    &visitor_ = visitor;
    #endif
            gsQuadRule<T> quRule ; // Quadrature rule
            gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
            gsVector<T> quWeights; // Temp variable for mapped weights

            const gsBasis<T> & basis = this->m_mb[0].basis(0).component(1); // = 0

            // Initialize reference quadrature rule and visitor data
            visitor_.initialize(basis,quRule);

            // Initialize domain element iterator -- using unknown 0
            typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator(boundary::none);

    #ifdef _OPENMP
            for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
    #else
            for (; domIt->good(); domIt->next() )
    #endif
            {
                // Map the Quadrature rule to the element
                quRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

                // Perform required evaluations on the quadrature nodes
                visitor_.evaluate(quNodes, this->m_mp);

                // Assemble on element
                visitor_.assemble(*domIt, quWeights);

                // Push to global matrix and right-hand side vector
    #pragma omp critical(localToGlobal)
                visitor_.localToGlobal( dirichletDofs, mSys); // omp_locks inside
            }
        }//omp parallel
    }

    void applyBeta(Visitor visitorBeta)
    {
    #pragma omp parallel
        {
            Visitor
    #ifdef _OPENMP
            // Create thread-private visitor
            visitor_Beta(visitorBeta);
            const int tid = omp_get_thread_num();
            const int nt  = omp_get_num_threads();
    #else
                    &visitor_Beta = visitorBeta;
    #endif
            gsQuadRule<T> quRule ; // Quadrature rule
            gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
            gsVector<T> quWeights; // Temp variable for mapped weights

            const gsBasis<T> & basis = this->m_mb[0].basis(0).component(1); // = 0

            // Initialize reference quadrature rule and visitor data
            visitor_Beta.initialize(basis,quRule);

            // Initialize domain element iterator -- using unknown 0
            typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator(boundary::none);

    #ifdef _OPENMP
            for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
    #else
            for (; domIt->good(); domIt->next() )
    #endif
            {
                // Map the Quadrature rule to the element
                quRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

                // Perform required evaluations on the quadrature nodes
                visitor_Beta.evaluateBeta(quNodes, this->m_mp, sol);

                // Assemble on element
                visitor_Beta.assembleBeta(*domIt, quWeights);

                // Push to global matrix and right-hand side vector
    #pragma omp critical(localToGlobal)
                visitor_Beta.localToGlobalBeta( dirichletDofsBeta, mSysBeta); // omp_locks inside
            }
        }//omp parallel
    }



    void solve()
    {
        typename gsSparseSolver<T>::SimplicialLDLT solver;

        solver.compute(mSys.matrix());
        sol = solver.solve(mSys.rhs()); // My solution
    }

    void solveBeta()
    {
        typename gsSparseSolver<T>::SimplicialLDLT solver;

        solver.compute(mSysBeta.matrix());
        solBeta = solver.solve(mSysBeta.rhs()); // My solution
    }

    void setGDEdge()
    {
        gsMatrix<T> solTMP(7, 1);
        gsMatrix<T> solBetaTMP(4, 1);

        solTMP(0, 0) = 1;
        solTMP(1, 0) = 1;
        solTMP(2, 0) = 1;
        solTMP(3, 0) = 1;
        solTMP(4, 0) = 0;
        solTMP(5, 0) = 0;
        solTMP(6, 0) = 0;

        solBetaTMP(0, 0) = 0;
        solBetaTMP(1, 0) = 0;
        solBetaTMP(2, 0) = 0;
        solBetaTMP(3, 0) = 0;

        sol = solTMP;
        solBeta = solBetaTMP;
    }
}; // class gsGluingData



} // namespace gismo

