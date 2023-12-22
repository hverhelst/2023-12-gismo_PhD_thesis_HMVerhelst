/** @file gsC1SurfEdge.h

    @brief Creates the (approx) C1 Edge space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Farahat & P. Weinmueller
*/

#pragma once

#include <gsUnstructuredSplines/src/gsG1AuxiliaryPatch.h>
#include <gsUnstructuredSplines/src/gsC1SurfBasisEdge.h>

#include <gsUnstructuredSplines/src/gsC1SurfGluingData.h>

namespace gismo
{
template<short_t d, class T>
class gsC1SurfEdge
{

private:

    /// Shared pointer for gsC1SurfEdge
    typedef memory::shared_ptr<gsC1SurfEdge> Ptr;

    /// Unique pointer for gsC1SurfEdge
    typedef memory::unique_ptr<gsC1SurfEdge> uPtr;


public:
    /// Empty constructor
    ~gsC1SurfEdge() { }

    gsC1SurfEdge(const gsMultiPatch<T> & mp, const boundaryInterface & item){
        auxGeom.push_back(gsG1AuxiliaryPatch<d,T>(mp.patch(item.first().patch), item.first().patch));
        auxGeom.push_back(gsG1AuxiliaryPatch<d,T>(mp.patch(item.second().patch), item.second().patch));
        m_item = item;
    }

    gsC1SurfEdge(const gsMultiPatch<T> & sp, const patchSide & item){
        auxGeom.push_back(gsG1AuxiliaryPatch<d,T>(sp.patch(item.patch), item.patch));
    }

    void computeG1InterfaceBasis()
    {
        basisEdgeResult.clear();

        gsMultiPatch<T> mp_init;
        mp_init.addPatch(auxGeom[0].getPatch());// Right -> 0 ====> v along the interface
        mp_init.addPatch(auxGeom[1].getPatch()); // Left -> 1 ====> u along the interface

        reparametrizeInterface(m_item);
        gsMultiPatch<T> test_mp; // auxGeom contains now the reparametrized geometry
        test_mp.addPatch(auxGeom[0].getPatch());
        test_mp.addPatch(auxGeom[1].getPatch());
        gsMultiBasis<T> test_mb(test_mp);
        gsMultiPatch<T> g1Basis_0, g1Basis_1;

        gsC1SurfGluingData<T> g1BasisEdge(test_mp, test_mb);
        gsC1SurfBasisEdge<T> g1BasisEdge_0(test_mp.patch(0), test_mb.basis(0), 1, false, g1BasisEdge);
        gsC1SurfBasisEdge<T> g1BasisEdge_1(test_mp.patch(1), test_mb.basis(1), 0, false, g1BasisEdge);
        g1BasisEdge_0.setG1BasisEdge(g1Basis_0);
        g1BasisEdge_1.setG1BasisEdge(g1Basis_1);

        gsMatrix<T> points(1,5);
        points << 0.0, 0.25, 0.5, 0.75, 1.0;
        gsMatrix<T> alphaL = g1BasisEdge.evalAlpha_L(points);
        gsMatrix<T> alphaR = g1BasisEdge.evalAlpha_R(points);
        gsMatrix<T> betaL = g1BasisEdge.evalBeta_L(points);
        gsMatrix<T> betaR = g1BasisEdge.evalBeta_R(points);
        gsMatrix<T> beta = g1BasisEdge.evalBeta(points);

//      Patch 0 -> Right
        auxGeom[0].parametrizeBasisBack(g1Basis_0);
//      Patch 1 -> Left
        auxGeom[1].parametrizeBasisBack(g1Basis_1);

        basisEdgeResult.push_back(auxGeom[0].getG1Basis());
        basisEdgeResult.push_back(auxGeom[1].getG1Basis());
    }

    void computeG1BoundaryBasis(const index_t boundaryInd)
    {
        basisEdgeResult.clear();

        reparametrizeBoundary(boundaryInd);
        gsMultiPatch<T> test_mp(auxGeom[0].getPatch());
        gsMultiBasis<T> test_mb(test_mp);
        gsMultiPatch<T> g1Basis_edge;

        gsC1SurfGluingData<T> bdyGD; // Empty constructor creates the sol and solBeta in a suitable way to manage the GD on the boundary
        gsC1SurfBasisEdge<T> g1BasisEdge(test_mp, test_mb, 1, true, bdyGD);
        g1BasisEdge.setG1BasisEdge(g1Basis_edge);

        auxGeom[0].parametrizeBasisBack(g1Basis_edge);

        basisEdgeResult.push_back(auxGeom[0].getG1Basis());
    }

    gsG1AuxiliaryPatch<d,T> & getSinglePatch(const index_t i)
    {
        return auxGeom[i];
    }

    std::vector<gsMultiPatch<T>> getBasis(){return basisEdgeResult;}

protected:

    // Store temp solution
    std::vector<gsMultiPatch<T>> basisEdgeResult;

    std::vector<gsG1AuxiliaryPatch<d,T>> auxGeom;

    // Need for avoiding the computeTopology
    boundaryInterface m_item;

private:

    // Compute topology
    // After computeTopology() the patches will have the same patch-index as the position-index in auxGeom
    // EXAMPLE: global patch-index-order inside auxGeom: [2, 3, 4, 1, 0]
    //          in auxTop: 2->0, 3->1, 4->2, 1->3, 0->4
    void computeAuxTopology();

    void reparametrizeInterface(const boundaryInterface & item);

    void reparametrizeBoundary(index_t side);


}; // Class gsC1SurfEdge

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsC1SurfEdge.hpp)
#endif