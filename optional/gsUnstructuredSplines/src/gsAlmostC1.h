/** @file gsAlmostC1.h

    @brief Creates the D-Patch smoothing matrix.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsCore/gsBoxTopology.h>
#include <gsCore/gsMultiPatch.h>
#include <gsIO/gsOptionList.h>
#include <gsUnstructuredSplines/src/gsDPatchBase.h>
// #include <gsUnstructuredSplines/src/gsDPatchBase.hpp>
namespace gismo
{


/**
 * @brief      Constructs the D-Patch, from which the transformation matrix can be called
 *
 * @tparam     d     parametric dimension
 *
 * @ingroup    UnstructuredSplines
 */
template<short_t d,class T>
class gsAlmostC1 : public gsDPatchBase<d,T>
{
protected:

typedef typename std::vector<std::tuple<index_t,index_t,T>> sparseEntry_t;

public:

    using Base = gsDPatchBase<d,T>;


    /// Shared pointer for gsAlmostC1
    typedef memory::shared_ptr< gsAlmostC1 > Ptr;

    /// Unique pointer for gsAlmostC1
    typedef memory::unique_ptr< gsAlmostC1 > uPtr;

    /// Empty constructor
    gsAlmostC1()
    :
    Base()
    {
    }

    // using Base::compute;
    // void compute();

    /**
     * @brief      Default constructor
     *
     * @param      mp    Multipatch of the geometry
     */
    gsAlmostC1(gsMultiPatch<T> const & mp) ;

    GISMO_CLONE_FUNCTION(gsAlmostC1)

    ~gsAlmostC1();

    using Base::defaultOptions;

    using Base::exportToPatches;

protected:
    /**
     * @brief       Computes the C1 coefficients for pre-multiplication to make the multipatch
     *
     * Takes the coefficients which are tagged as "free" in the modified DoFMapper (m_mapModified) and when a boundary vertex with valence=3 is present, this one is shifted.
     *
     */
    gsMatrix<T> _preCoefficients();
    gsMatrix<T> _preCoefficients(const gsMultiPatch<T> & patches)
    {
        GISMO_UNUSED(patches);
        return _preCoefficients();
    }

    /**
     * @brief       Computes the C1 coefficients for pre-multiplication to make the multipatch
     *
     * Takes the coefficients which are tagged as "free" in the modified DoFMapper (m_mapModified)
     *
     */
    gsMatrix<T> freeCoefficients();

    /**
     * @brief      Set the coefficients of mp to \a coefs
     *
     * @param[in]  coefs  The coefs
     * @param      mp     The multipatch to update
     *
     */
    void setCoefficients(const gsMatrix<T> & coefs, gsMultiPatch<T> & mp) const;


    /**
     * @brief       Exports a single modified patch with index \a patch
     *
     * The patch is obtained by transforming the coefficients of the D-Patch to the original basis, such that the original basis functions can be used to plot the geometry (and the patch structure will remain intact).
     * To construct the geometry, the coefficients for the C1 basis are multiplied with the transpose of the transformation matrix. The C1 coefficients are obtained with \ref preCoefficients().
     *
     */
    // gsGeometry<T>* exportPatch(index_t patch, bool computeCoefs=true);
    using Base::exportPatch;

    /**
     * @brief      Exports the modified geometry to a @a gsMultiPatch object
     *
     * @return     A multipatch with the geometry
     */
    // gsMultiPatch<T> exportToPatches();
    // using Base::exportPatch;

    /**
     * @brief      Returns the smoothing matrix into \a matrix
     *
     * @param      matrix  The matrix
     */
    // const void matrix_into(gsSparseMatrix<T> & matrix) const
    // { matrix = m_matrix; }

    /**
     * @brief      Returns the smoothing matrix
     *
     * The number of columns of the matrix corresponds to the number of basis functions in the local basis; this is the sum of all the basis functions over all the patches.
     * The number of rows of the matrix corresponds to the number of global basis functions, i.e. the number of basis functions corresponding to the D-Patch.
     * Multiplying the basis with the local basis function values gives the values of the basis functions in the global basis.
     *
     * @return     A matrix \a result to transfer local to global coefficients
     */
    // const gsSparseMatrix<T> matrix() const
    // {
    //     gsSparseMatrix<T> result; matrix_into(result);
    //     return result;
    // }

protected:

    using Base::getSharpCorners;

    gsMatrix<T> _getNormals(const std::vector<patchCorner> & corners) const;

    std::tuple<gsMatrix<T>,gsMatrix<T>,gsMatrix<index_t>> _makeTriangle(const patchCorner & corner) const;

    gsMatrix<T,3,3> _getRotationMatrix(const gsVector<T,3> & a, const gsVector<T,3> & b) const;

    using Base::_indexFromSides;
    using Base::_indexFromVert;
    using Base::_vertexData;
    using Base::_sideIndex;
    using Base::_vertIndex;
    using Base::_getLowestCorners;
    using Base::_removeLowestCorners;
    using Base::_getLowestIndices;
    using Base::_removeLowestIndices;
    using Base::_getInterfaceIndices;
    using Base::_getAllInterfaceIndices;

protected:

    /**
     * @brief      Initializes the matrix, the basis and the mappers
     */

    void _countDoFs();

    /**
     * @brief      Computes the modified mapper
     *
     * The modified mapper is computed based on the elimination of different functions with different conditions.
     * 1) For interfaces, it eliminates all the nodes except the first two and the last two
     * 2) For boundaries, there is no elimination
     * 3) For vertices, there are few options
     *  a) Boundary vertices
     *      i)  Valence 1: No eliminations
     *      ii) Valence 2: the two outer-most basis functions on the interface (the one at the vertex and the one next to it on the interface) are both eliminated
     *      iii)Valence 3: On all the patches, the basis functions corresponding to the vertex are coupled to eachother. The basis functions next to this one (on an interface OR on a boundary) are eliminated
     *  b) Interior vertices: all basis functions along the interface are eliminated if not done so
     */
    // void _computeMapper(); // also initialize the mappers!
    using Base::_computeMapper;

    /**
     * @brief      Calculates the smoothing matrix.
     */
    // void _computeSmoothMatrix();
    using Base::_computeSmoothMatrix;

    void _initBasis();

    void _initTHB();

    void _refBoxes(std::vector<std::vector<index_t>> & patchBoxes);

    /**
     * @brief      Prepares the THB basis if needed.
     *
     * This function constructs THB refinements on the places where they are needed, i.e. around EVs. It also constructs the transfer matrix (m_tMatrix) forms the transformation between the original B-spline basis and the THB-Spline basis.
     */
    void _makeTHB();

    /**
     * @brief      Computes D-Patch smoothing
     *
     * Given a basis with THB refinement around the EVs, this function computes the D-Patch smoothing
     */
    void _computeEVs();

protected:

    std::vector<std::vector<patchCorner> > _getSpecialCornerLists(const gsMultiPatch<T> & patches);

    using Base::_push;
    using Base::_pushAndCheck;


protected:
    // void _computeInterfaceMapper(boundaryInterface iface);
    using Base::_computeInterfaceMapper;

    // void _computeBoundaryMapper(patchSide boundary);
    using Base::_computeBoundaryMapper;

    void _computeVertexMapper(patchCorner pcorner);


private:
    // // Boundary vertex of valence 1
    // template<bool _boundary, index_t _v> // valence=2
    // typename std::enable_if<  _boundary && _v==1, void>::type
    // _computeVertexMapperBoundary_v1(patchCorner pcorner, index_t valence);
    using Base::_computeMapperRegularCorner_v1;

    // // Boundary vertex of valence 2 with C1 smoothness
    // template<bool _boundary, index_t _v, bool _smooth> // valence=2
    // typename std::enable_if<  _boundary && _v==2 && _smooth, void>::type
    // _computeVertexMapperBoundarySmooth_v2(patchCorner pcorner, index_t valence);
    using Base::_computeMapperRegularBoundaryVertexSmooth_v2;

    // Boundary vertex of valence 2 with C0 smoothness
    // template<bool _boundary, index_t _v, bool _smooth> // valence=2
    // typename std::enable_if<  _boundary && _v==2 && (!_smooth), void>::type
    void _computeMapperRegularBoundaryVertexNonSmooth_v2(patchCorner pcorner, index_t valence);

    // Boundary vertex of valence !(1,2,3) with C1 smoothness
    // template<bool _boundary, index_t _v, bool _smooth>
    // typename std::enable_if<  _boundary && _v==-1 && _smooth, void>::type
    void _computeMapperIrregularBoundaryVertexSmooth_v(patchCorner pcorner, index_t valence);

    // Boundary vertex of valence !(1,2,3) with C0 smoothness
    // template<bool _boundary, index_t _v, bool _smooth>
    // typename std::enable_if<  _boundary && _v==-1 && (!_smooth), void>::type
    // void _computeVertexMapperBoundaryNonSmooth_v(patchCorner pcorner, index_t valence);
    using Base::_computeMapperIrregularBoundaryVertexNonSmooth_v;

    // Ordinary interior vertex
    // template<bool _boundary, index_t _v> // valence=2
    // typename std::enable_if<  (!_boundary) && _v==4, void>::type
    void _computeMapperInteriorVertex_v4(patchCorner pcorner, index_t valence);

    // Extraordinary interior vertex
    // template<bool _boundary, index_t _v>
    // typename std::enable_if<  (!_boundary) && _v==-1, void>::type
    void _computeMapperInteriorVertex_v(patchCorner pcorner, index_t valence);


protected:
    /**
     * @brief      Handles a vertex in the global matrix
     *
     * We use the following notation convention (per patch!):
     * b00 is the basis function at the vertex
     * b10 is the basis function next to the vertex along the first interface that connects to the vertex
     * b20 is the basis function next to b10 along the first interface that connects to the vertex
     * etc.
     *
     * b01 is the basis function next to the vertex along the second interface that connects to the vertex
     * b02 is the basis function next to b01 along the second interface that connects to the vertex
     * etc.
     *
     * b11 is the basis function with offset 1 from both interfaces and from the vertex itself
     * b22 is the basis function with offset 2 from both interfaces and from the vertex itself
     * etc.
     *
     * There are different options.
     * a) Boundary vertices
     *      i)  Valence 1: b00, b10, b01 and b00 all get weight 1.0 w.r.t the same basis function in the local basis
     *      ii) Valence 2: This case contains an interface between two patches. We use index k to denote the row basis functions along the interface. So k=0 corresponds to the basis functions on the boundary and k=1 corresponds to the basis functions with offset 1 from the boundaries. Using this convention, the functions bk1 in the local basis, are coupled to bk1 in the global basis with weight 1. The functions bk0 in the local basis (on the considered patch) are coupled to bk1 in the global basis with weight 0.5. The functions bk0 in the local basis (on the other patch) are coupled to bk1 (on the considered patch) in the global basis with weight 0.5.
     *      iii)Valence 3: In this case, the matched vertices on all the adjacent patches are treated in one go! Note that all the basis functions corresponding to the vertex (b00) on all patches are matched! We couple the b00 functions of all patches (in the local basis) with weight 1/4 to the b00 of the adjacent patch with the lowest number in the global basis. Then, the b11 on the considered patch is coupled with weight 1 to itself and with weight 0.25 to the b00s of the other patches. Then, we will handle the vertices where an interface and a boundary meet (there are two of these). For the patch corners that are on an interface, we find the b11 and b10 vertices (orthogonal to the interface) and we give all b10s weight 0.5 w.r.t. the b11s in the global basis (on both patches). Lastly, we add weight 0.5 for the b10s along the boundaries (so only for two patches) to the (matched) b00 basis function (all b00s refer to the same dof in the global basis).
     * b) Interior vertices (all valences):
     *      i)  b11 gets weight 1.0 w.r.t the same basis function in the local basis
     *      ii) all associated b00s in the local basis get weight 1/valence w.r.t. b11 in the global basis
     *      iii)the b10 and b01 in the local basis get weight 1/2 w.r.t. b11 in the global basis
     *
     * @param[in]  pcorner  The patchcorner
     */
    // void _handleVertex(patchCorner pcorner);

    // interior vertices
    void _handleInteriorVertex(patchCorner pcorner, index_t valence);

    /**
     * @brief      Handles an interface in the global matrix
     *
     * Gives all the DoFs that have offset 1 (orthogonal) from the interface weight 1.0 w.r.t itself. All the DoFs ON the interface (on both patches) will have weight 0.5 to the DoF with offset 1.
     * This interface handling excludes the indices that are in the 0 and 1 ring around vertices.
     *
     * @param[in]  iface  The interface
     */
    // void _handleInterface(boundaryInterface iface);
    /**
     * @brief      Handles a boundary in the global matrix
     *
     * Handles all DoFs on the boundary with unit-weight, except the ones in the 0 and 1 rings around the vertices.
     *
     * @param[in]  side  The boundary side
     */
    // void _handleBoundary(patchSide side);
    /**
     * @brief      Handles the interior in the global matrix
     *
     * Gives all left-over DoFs, which are in the interior, weight 1 w.r.t. itself
     */
    // void _handleInterior();

private:
    /**
     * @brief      Handles a regular corner
     *
     * @param[in]  pcorner  The pcorner
     */
    // void _handleRegularCorner(patchCorner pcorner);

    // template<bool _regular, bool _smooth> // valence=2
    // typename std::enable_if<  _regular  &&   _smooth   , void>::type
    // _handleBoundaryVertex(patchCorner pcorner, index_t valence);

    // template<bool _regular, bool _smooth> // valence=2
    // typename std::enable_if<  _regular  && (!_smooth)   , void>::type
    // _handleBoundaryVertex(patchCorner pcorner, index_t valence);

    // template<bool _regular, bool _smooth> // valence > 2
    // typename std::enable_if<(!_regular) &&   _smooth    , void>::type
    // _handleBoundaryVertex(patchCorner pcorner, index_t valence);

    // template<bool _regular, bool _smooth> // valence > 1
    // typename std::enable_if<(!_regular) && (!_smooth)   , void>::type
    // _handleBoundaryVertex(patchCorner pcorner, index_t valence);

    // void _handleRegularBoundaryVertexSmooth(patchCorner pcorner, index_t valence);

    void _handleRegularBoundaryVertexNonSmooth(patchCorner pcorner, index_t valence);

    void _handleIrregularBoundaryVertexSmooth(patchCorner pcorner, index_t valence);

    void _handleIrregularBoundaryVertexNonSmooth(patchCorner pcorner, index_t valence);

protected:
    using Base::_whichHandled;

protected:
    using Base::m_patches;
    gsMultiPatch<T> m_RefPatches;
    using Base::m_bases;
    using Base::m_topology;

    using Base::m_Bbases;
    using Base::m_tMatrix;
    using Base::m_sideCheck;
    using Base::m_vertCheck;
    using Base::m_basisCheck;
    using Base::m_C0s;

    using Base::m_mapModified;
    using Base::m_mapOriginal;

    using Base::m_matrix;

    using Base::m_options;

    using Base::m_size;

    using Base::m_coefs;

    using Base::m_nSides;
    using Base::m_nVerts;
};



}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsAlmostC1.hpp)
#endif
