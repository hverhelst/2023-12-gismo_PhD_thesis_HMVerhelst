/** @file gsPatchReparameterized.h
 *
    @brief Reparametrize one Patch

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Farahat & P. Weinmueller
*/
#pragma once

#include <gismo.h>
#include <gsUnstructuredSplines/src/gsContainerBasis.h>

namespace gismo
{

template<short_t d, class T>
class gsPatchReparameterized
{
public:

    gsPatchReparameterized()
    {}

    gsPatchReparameterized(const gsGeometry<T> & patch, gsBasis<T> & singlePatch)
    : m_patchRotated(patch), m_basisRotated(singlePatch)
    {
        rotationNum = 0;
        axisOrientation = 0;

        mapIndex.setZero(4);
        mapIndex[0] = 1;
        mapIndex[1] = 2;
        mapIndex[2] = 3;
        mapIndex[3] = 4;
    };

    void swapAxis()
    {
/*
        gsTensorBSplineBasis<d, T> & temp_L = dynamic_cast<gsTensorBSplineBasis<d, T> &>(m_patchRotated.basis(0));
        gsBSplineBasis<T> temp_basisLU = dynamic_cast<gsBSplineBasis<T> &>(temp_L.component(0));
        gsBSplineBasis<T> temp_basisLV = dynamic_cast<gsBSplineBasis<T> &>(temp_L.component(1));
        index_t dimU = temp_L.size(0);
        index_t dimV = temp_L.size(1);

        // The number of cols has to match the dimension of the space
        gsMatrix<T> mpar(dimU * dimV, m_patchRotated.targetDim());
        for (index_t j = 0; j < dimU; j++)
        {   // Loop over the rows
            for (index_t i = 0; i < dimV; i++)
            {
                mpar.row(i + j * dimV) = m_patchRotated.patch(0).coefs().row(j + dimU * i);
            }
        }
        // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
        gsTensorBSpline<d, T> newgeom1(temp_basisLV.knots(), temp_basisLU.knots(), mpar);

        // Create a new single-patch object
        gsMultiPatch<T> newpatch;
        newpatch.addPatch(newgeom1);
        m_patchRotated.swap(newpatch);
*/
        gsBasis<T> & temp_L_basis = m_patchRotated.basis(0);
        index_t dimU = temp_L_basis.component(0).size();
        index_t dimV = temp_L_basis.component(1).size();

        // The number of cols has to match the dimension of the space
        gsMatrix<T> mpar(dimU * dimV, m_patchRotated.targetDim());
        if (gsTensorBSplineBasis<d, T> * mb =
                dynamic_cast<gsTensorBSplineBasis<d, T> *>(&m_patchRotated.basis(0)) )
        {
            for (index_t j = 0; j < dimU; j++)
               for (index_t i = 0; i < dimV; i++)
                    mpar.row(i + j * dimV) = m_patchRotated.patch(0).coefs().row(j + dimU * i);


            gsTensorBSplineBasis<d, T> & temp_L = dynamic_cast<gsTensorBSplineBasis<d, T> &>(m_patchRotated.basis(0));

            // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
            gsTensorBSpline<d, T> newgeom1(temp_L.knots(1), temp_L.knots(0), mpar);

            // Create a new single-patch object
            gsMultiPatch<T> newpatch(newgeom1);
            m_patchRotated.swap(newpatch);
        }
        else if (gsTensorNurbsBasis<d, T> * mb =
                dynamic_cast<gsTensorNurbsBasis<d, T> *>(&m_patchRotated.basis(0)) )
        {
            gsMatrix<T> weight_rep(dimU * dimV, 1);
            for (index_t j = 0; j < dimU; j++)
            {   // Loop over the rows
                for (index_t i = 0; i < dimV; i++)
                {
                    mpar.row(i + j * dimV) = m_patchRotated.patch(0).coefs().row(j + dimU * i);
                    weight_rep.row(i + j * dimV) = mb->weights().row(j + dimU * i);
                }
            }

            gsTensorNurbsBasis<d, T> & temp_L = dynamic_cast<gsTensorNurbsBasis<d, T> &>(m_patchRotated.basis(0));

            // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
            gsTensorNurbs<d, T> newgeom1(temp_L.knots(1), temp_L.knots(0), mpar, weight_rep);

            // Create a new single-patch object
            gsMultiPatch<T> newpatch(newgeom1);
            m_patchRotated.swap(newpatch);
        }
        else
            gsInfo << "Only gsTensorBSplineBasis<d, T> and gsTensorNurbsBasis<d, T> are implemented \n";


        // BASES
        //m_basisRotated.swapAxis();
        m_basisRotated.basis(0).reverse();

        // Map Index
        mapIndex[0] = 3;
        mapIndex[1] = 4;
        mapIndex[2] = 1;
        mapIndex[3] = 2;

        axisOrientation = 1;
    }


    void rotateParamAntiClock(){
/*
        gsTensorBSplineBasis<d, T> & temp_L = dynamic_cast<gsTensorBSplineBasis<d, T> &>(m_patchRotated.basis(0));
        gsBSplineBasis<T> temp_basisLU = dynamic_cast<gsBSplineBasis<T> &>(temp_L.component(0));
        gsBSplineBasis<T> temp_basisLV = dynamic_cast<gsBSplineBasis<T> &>(temp_L.component(1));
        index_t dimU = temp_L.size(0);
        index_t dimV = temp_L.size(1);

        // The number of cols has to match the dimension of the space
        gsMatrix<T> mpar(dimU * dimV, m_patchRotated.targetDim());

        // Loop over the cols
        for (index_t j = 0; j < dimU; j++)
        {   // Loop over the rows
            for (index_t i = 0; i < dimV; i++)
            {
                mpar.row(i + j * dimV) = m_patchRotated.patch(0).coefs().row((dimU - 1 - j) + dimU * i);
            }
        }
        temp_basisLU.knots().reverse(); // For different regularity
        // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
        gsTensorBSpline<d, T> newgeom1(temp_basisLV.knots(), temp_basisLU.knots(), mpar);

        // Create a new single-patch object
        gsMultiPatch<T> newpatch(newgeom1);
        m_patchRotated.swap(newpatch);
*/
        gsBasis<T> & temp_L_basis = m_patchRotated.basis(0);
        index_t dimU = temp_L_basis.component(0).size();
        index_t dimV = temp_L_basis.component(1).size();

        // The number of cols has to match the dimension of the space
        gsMatrix<T> mpar(dimU * dimV, m_patchRotated.targetDim());
        if (gsTensorBSplineBasis<d, T> * mb =
                dynamic_cast<gsTensorBSplineBasis<d, T> *>(&m_patchRotated.basis(0)) )
        {
            for (index_t j = 0; j < dimU; j++)
                for (index_t i = 0; i < dimV; i++)
                    mpar.row(i + j * dimV) = m_patchRotated.patch(0).coefs().row((dimU - 1 - j) + dimU * i);

            gsTensorBSplineBasis<d, T> & temp_L = dynamic_cast<gsTensorBSplineBasis<d, T> &>(m_patchRotated.basis(0));
            temp_L.knots(0).reverse(); // For different regularity

            // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
            gsTensorBSpline<d, T> newgeom1(temp_L.knots(1), temp_L.knots(0), mpar);

            // Create a new single-patch object
            gsMultiPatch<T> newpatch(newgeom1);
            m_patchRotated.swap(newpatch);
        }
        else if (gsTensorNurbsBasis<d, T> * mb =
                dynamic_cast<gsTensorNurbsBasis<d, T> *>(&m_patchRotated.basis(0)) )
        {
            gsMatrix<T> weight_rep(dimU * dimV, 1);
            for (index_t j = (dimU - 1 ) ; j >= 0; j--)
            {   // Loop over the rows
                for (index_t i = 0; i < dimV; i++)
                {
                    mpar.row(i + j * dimV) = m_patchRotated.patch(0).coefs().row((dimU - 1 - j) + dimU * i);
                    weight_rep.row(i + j * dimV) = mb->weights().row((dimU - 1 - j) + dimU * i);
                }
            }

            gsTensorNurbsBasis<d, T> & temp_L = dynamic_cast<gsTensorNurbsBasis<d, T> &>(m_patchRotated.basis(0));
            temp_L.knots(0).reverse(); // For different regularity

            // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
            gsTensorNurbs<d, T> newgeom1(temp_L.knots(1), temp_L.knots(0), mpar, weight_rep);

            // Create a new single-patch object
            gsMultiPatch<T> newpatch(newgeom1);
            m_patchRotated.swap(newpatch);
        }
        else
            gsInfo << "Only gsTensorBSplineBasis<d, T> and gsTensorNurbsBasis<d, T> are implemented \n";



        // BASES
        //m_basisRotated.swapAxis();
        m_basisRotated.basis(0).reverse();

        // Map Index
        if (getOrient() == 0)
        {
            mapIndex[0] = 4;
            mapIndex[1] = 3;
            mapIndex[2] = 1;
            mapIndex[3] = 2;
        }
        else
        {
            mapIndex[0] = 1;
            mapIndex[1] = 2;
            mapIndex[2] = 4;
            mapIndex[3] = 3;
        }


        // Update the number of rotation of the axis
        rotationNum++;
    }

    void rotateParamClock()
    {
/*
        gsTensorBSplineBasis<d, T> & temp_L = dynamic_cast<gsTensorBSplineBasis<d, T> &>(m_patchRotated.basis(0));
        gsBSplineBasis<T> temp_basisLU = dynamic_cast<gsBSplineBasis<T> &>(temp_L.component(0));
        gsBSplineBasis<T> temp_basisLV = dynamic_cast<gsBSplineBasis<T> &>(temp_L.component(1));
        index_t dimU = temp_L.size(0);
        index_t dimV = temp_L.size(1);
*/

        gsBasis<T> & temp_L_basis = m_patchRotated.basis(0);
        index_t dimU = temp_L_basis.component(0).size();
        index_t dimV = temp_L_basis.component(1).size();

        // The number of cols has to match the dimension of the space
        gsMatrix<T> mpar(dimU * dimV, m_patchRotated.targetDim());
        if (gsTensorBSplineBasis<d, T> * mb =
                dynamic_cast<gsTensorBSplineBasis<d, T> *>(&m_patchRotated.basis(0)) )
        {
            for (index_t j = (dimU - 1 ) ; j >= 0; j--)
                for (index_t i = 0; i < dimV; i++)
                    mpar.row(i + (dimU - j - 1) * dimV) = m_patchRotated.patch(0).coefs().row((dimV * dimU  -1 - j) - dimU * i);


            gsTensorBSplineBasis<d, T> & temp_L = dynamic_cast<gsTensorBSplineBasis<d, T> &>(m_patchRotated.basis(0));
            temp_L.knots(1).reverse(); // For different regularity

            // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
            gsTensorBSpline<d, T> newgeom1(temp_L.knots(1), temp_L.knots(0), mpar);

            // Create a new single-patch object
            gsMultiPatch<T> newpatch(newgeom1);
            m_patchRotated.swap(newpatch);
        }
        else if (gsTensorNurbsBasis<d, T> * mb =
                dynamic_cast<gsTensorNurbsBasis<d, T> *>(&m_patchRotated.basis(0)) )
        {
            gsMatrix<T> weight_rep(dimU * dimV, 1);
            for (index_t j = (dimU - 1 ) ; j >= 0; j--)
            {   // Loop over the rows
                for (index_t i = 0; i < dimV; i++)
                {
                    mpar.row(i + (dimU - j - 1) * dimV) = m_patchRotated.patch(0).coefs().row((dimV * dimU  -1 - j) - dimU * i);
                    weight_rep.row(i + (dimU - j - 1) * dimV) = mb->weights().row((dimV * dimU  -1 - j) - dimU * i);
                }
            }

            gsTensorNurbsBasis<d, T> & temp_L = dynamic_cast<gsTensorNurbsBasis<d, T> &>(m_patchRotated.basis(0));
            temp_L.knots(1).reverse(); // For different regularity

            // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
            gsTensorNurbs<d, T> newgeom1(temp_L.knots(1), temp_L.knots(0), mpar, weight_rep);

            // Create a new single-patch object
            gsMultiPatch<T> newpatch(newgeom1);
            m_patchRotated.swap(newpatch);
        }
        else
            gsInfo << "Only gsTensorBSplineBasis<d, T> and gsTensorNurbsBasis<d, T> are implemented \n";


        // BASES
        //m_basisRotated.swapAxis();
        m_basisRotated.basis(0).reverse();

        // Map Index
        if (getOrient() == 0)
        {
            mapIndex[0] = 3;
            mapIndex[1] = 4;
            mapIndex[2] = 2;
            mapIndex[3] = 1;
        }
        else
        {
            mapIndex[0] = 2;
            mapIndex[1] = 1;
            mapIndex[2] = 3;
            mapIndex[3] = 4;
        }

        // Update the number of rotation of the axis
        rotationNum--;
        checkRotation();
    }

    void rotateParamAntiClockTwice()
    {
/*
        gsTensorBSplineBasis<d, T> & temp_L = dynamic_cast<gsTensorBSplineBasis<d, T> &>(m_patchRotated.basis(0));
        gsBSplineBasis<T> temp_basisLU = dynamic_cast<gsBSplineBasis<T> &>(temp_L.component(0));
        gsBSplineBasis<T> temp_basisLV = dynamic_cast<gsBSplineBasis<T> &>(temp_L.component(1));
        index_t dimU = temp_L.size(0);
        index_t dimV = temp_L.size(1);

        // The number of cols has to match the dimension of the space
        gsMatrix<T> mpar(dimU * dimV, m_patchRotated.targetDim());

        for (index_t i = 0; i < ( dimU * dimV ); i++)
        {
            mpar.row(i) = m_patchRotated.patch(0).coefs().row((dimU * dimV - 1) - i);
        }
        temp_basisLU.knots().reverse(); // For different regularity
        temp_basisLV.knots().reverse(); // For different regularity

        // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
        gsTensorBSpline<d, T> newgeom1(temp_basisLU.knots(), temp_basisLV.knots(), mpar);

        // Create a new single-patch object
        gsMultiPatch<T> newpatch(newgeom1);
        m_patchRotated.swap(newpatch);
*/
        gsBasis<T> & temp_L_basis = m_patchRotated.basis(0);
        index_t dimU = temp_L_basis.component(0).size();
        index_t dimV = temp_L_basis.component(1).size();

        // The number of cols has to match the dimension of the space
        gsMatrix<T> mpar(dimU * dimV, m_patchRotated.targetDim());
        if (gsTensorBSplineBasis<d, T> * mb =
                dynamic_cast<gsTensorBSplineBasis<d, T> *>(&m_patchRotated.basis(0)) )
        {
            for (index_t i = 0; i < ( dimU * dimV ); i++)
                mpar.row(i) = m_patchRotated.patch(0).coefs().row((dimU * dimV - 1) - i);

            gsTensorBSplineBasis<d, T> & temp_L = dynamic_cast<gsTensorBSplineBasis<d, T> &>(m_patchRotated.basis(0));
            temp_L.knots(0).reverse(); // For different regularity
            temp_L.knots(1).reverse(); // For different regularity

            // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
            gsTensorBSpline<d, T> newgeom1(temp_L.knots(0), temp_L.knots(1), mpar);

            // Create a new single-patch object
            gsMultiPatch<T> newpatch(newgeom1);
            m_patchRotated.swap(newpatch);
        }
        else if (gsTensorNurbsBasis<d, T> * mb =
                dynamic_cast<gsTensorNurbsBasis<d, T> *>(&m_patchRotated.basis(0)) )
        {
            gsMatrix<T> weight_rep(dimU * dimV, 1);
            for (index_t i = 0; i < ( dimU * dimV ); i++)
            {
                mpar.row(i) = m_patchRotated.patch(0).coefs().row((dimU * dimV - 1) - i);
                weight_rep.row(i) = mb->weights().row((dimU * dimV - 1) - i);
            }

            gsTensorNurbsBasis<d, T> & temp_L = dynamic_cast<gsTensorNurbsBasis<d, T> &>(m_patchRotated.basis(0));
            temp_L.knots(0).reverse(); // For different regularity
            temp_L.knots(1).reverse(); // For different regularity

            // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
            gsTensorNurbs<d, T> newgeom1(temp_L.knots(0), temp_L.knots(1), mpar, weight_rep);

            // Create a new single-patch object
            gsMultiPatch<T> newpatch(newgeom1);
            m_patchRotated.swap(newpatch);
        }
        else
            gsInfo << "Only gsTensorBSplineBasis<d, T> and gsTensorNurbsBasis<d, T> are implemented \n";


        // Map Index
        if (getOrient() == 0)
        {
            mapIndex[0] = 2;
            mapIndex[1] = 1;
            mapIndex[2] = 4;
            mapIndex[3] = 3;
        }
        else
        {
            mapIndex[0] = 4;
            mapIndex[1] = 3;
            mapIndex[2] = 2;
            mapIndex[3] = 1;
        }

        // Update the number of rotation of the axis (anti-clockwise)
        rotationNum+=2;
        checkRotation();
    }

    void parametrizeBasisBack(gsMultiPatch<T> & g1Basis)
    {
        switch (rotationNum)
        {
            case 2:
                //gsInfo << "Patch Basis rotated twice anticlockwise\n";
                rotateBasisAntiClockTwice(g1Basis);
                break;
            case -1:
                //gsInfo << "Patch Basis rotated anticlockwise\n";
                rotateBasisAntiClock(g1Basis);
                break;
            case 1:
                //gsInfo << "Patch Basis rotated clockwise\n";
                rotateBasisClock(g1Basis);
                break;
            case 0:
                //gsInfo << "Patch Basis not rotated\n";
                break;
            default:
                break;
        }

        if(axisOrientation)
            swapBasisAxis(g1Basis);

    }

    void rotateBasisAntiClockTwice(gsMultiPatch<T> & g1Basis)
    {
        gsMultiPatch<T> newpatch;
        for(size_t np = 0; np < g1Basis.nPatches(); np++)
        {
            gsMultiBasis<T> auxBase(g1Basis.patch(np));
            gsTensorBSplineBasis<d, T>
                & temp_L = dynamic_cast<gsTensorBSplineBasis<d, T> &>(auxBase.basis(0));
            gsBSplineBasis<T> temp_basisLU = dynamic_cast<gsBSplineBasis<T> &>(temp_L.component(0));
            gsBSplineBasis<T> temp_basisLV = dynamic_cast<gsBSplineBasis<T> &>(temp_L.component(1));
            index_t dimU = temp_L.size(0);
            index_t dimV = temp_L.size(1);

            // The number of cols has to match the dimension of the space
            gsMatrix<T> mpar(dimU * dimV, g1Basis.patch(np).targetDim());

            for (index_t i = 0; i < (dimU * dimV); i++)
            {
                mpar.row(i) = g1Basis.patch(np).coefs().row((dimU * dimV - 1) - i);
            }
            // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
            gsTensorBSpline<d, T> newgeom1(temp_basisLU.knots(), temp_basisLV.knots(), mpar);

            newpatch.addPatch(newgeom1);
        }
        g1Basis.swap(newpatch);
    }

    void rotateBasisAntiClock(gsMultiPatch<T> & g1Basis)
    {
        gsMultiPatch<T> newpatch;
        for(size_t np = 0; np < g1Basis.nPatches(); np++)
        {
            gsMultiBasis<T> auxBase(g1Basis.patch(np));
            gsTensorBSplineBasis<d, T>
                & temp_L = dynamic_cast<gsTensorBSplineBasis<d, T> &>(auxBase.basis(0));
            gsBSplineBasis<T> temp_basisLU = dynamic_cast<gsBSplineBasis<T> &>(temp_L.component(0));
            gsBSplineBasis<T> temp_basisLV = dynamic_cast<gsBSplineBasis<T> &>(temp_L.component(1));
            index_t dimU = temp_L.size(0);
            index_t dimV = temp_L.size(1);

            // The number of cols has to match the dimension of the space
            gsMatrix<T> mpar(dimU * dimV, g1Basis.patch(np).targetDim());

            // Loop over the cols
            for (index_t j = 0; j < dimU; j++)
            {   // Loop over the rows
                for (index_t i = 0; i < dimV; i++)
                {
                    mpar.row(i + j * dimV) = g1Basis.patch(np).coefs().row((dimU - 1 - j) + dimU * i);
                }
            }

            // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
            gsTensorBSpline<d, T> newgeom1(temp_basisLV.knots(), temp_basisLU.knots(), mpar);

            newpatch.addPatch(newgeom1);
        }
        g1Basis.swap(newpatch);
    }

    void rotateBasisClock(gsMultiPatch<T> & g1Basis)
    {
        gsMultiPatch<T> newpatch;
        for(size_t np = 0; np < g1Basis.nPatches(); np++)
        {
            gsMultiBasis<T> auxBase(g1Basis.patch(np));
            gsTensorBSplineBasis<d, T>
                & temp_L = dynamic_cast<gsTensorBSplineBasis<d, T> &>(auxBase.basis(0));
            gsBSplineBasis<T> temp_basisLU = dynamic_cast<gsBSplineBasis<T> &>(temp_L.component(0));
            gsBSplineBasis<T> temp_basisLV = dynamic_cast<gsBSplineBasis<T> &>(temp_L.component(1));
            index_t dimU = temp_L.size(0);
            index_t dimV = temp_L.size(1);

            // The number of cols has to match the dimension of the space
            gsMatrix<T> mpar(dimU * dimV, g1Basis.patch(np).targetDim());

            for (index_t j = (dimU - 1); j >= 0; j--)
            {   // Loop over the rows
                for (index_t i = 0; i < dimV; i++)
                {
                    mpar.row(i + (dimU - j - 1) * dimV) =
                        g1Basis.patch(np).coefs().row((dimV * dimU - 1 - j) - dimU * i);
                }
            }

            // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
            gsTensorBSpline<d, T> newgeom1(temp_basisLV.knots(), temp_basisLU.knots(), mpar);

            newpatch.addPatch(newgeom1);
        }
        g1Basis.swap(newpatch);
    }

    void swapBasisAxis(gsMultiPatch<T> & g1Basis)
    {
        gsMultiPatch<T> newpatch;
        for(size_t np = 0; np < g1Basis.nPatches(); np++)
        {
            gsMultiBasis<T> auxBase(g1Basis.patch(np));
            gsTensorBSplineBasis<d, T> & temp_L = dynamic_cast<gsTensorBSplineBasis<d, T> &>(auxBase.basis(0));
            gsBSplineBasis<T> temp_basisLU = dynamic_cast<gsBSplineBasis<T> &>(temp_L.component(0));
            gsBSplineBasis<T> temp_basisLV = dynamic_cast<gsBSplineBasis<T> &>(temp_L.component(1));
            index_t dimU = temp_L.size(0);
            index_t dimV = temp_L.size(1);

            // The number of cols has to match the dimension of the space
            gsMatrix<T> mpar(dimU * dimV, g1Basis.patch(np).targetDim());

            for (index_t j = 0; j < dimU; j++)
            {   // Loop over the rows
                for (index_t i = 0; i < dimV; i++)
                {
                    mpar.row(i + j * dimV) = g1Basis.patch(np).coefs().row(j + dimU * i);
                }
            }

            // Create a new geometry starting from kntot vectors and the matrix of the coefficients reparametrized
            gsTensorBSpline<d, T> newgeom1(temp_basisLV.knots(), temp_basisLU.knots(), mpar);

            newpatch.addPatch(newgeom1);
        }
        g1Basis.swap(newpatch);
    }


    const index_t getOrient(){
        return axisOrientation;
    }

    void checkRotation()
    {
        if(rotationNum == 4)
            rotationNum = 0;
    }

    void setNumberOfRotation(index_t numRot)
    {
        rotationNum = numRot;
    }

    const index_t getNumberOfRotation(){
        return rotationNum;
    }

    gsGeometry<T> & getPatchRotated() { return m_patchRotated.patch(0); }
    gsBasis<T> & getBasisRotated() { return m_basisRotated.basis(0); }

    index_t getMapIndex(index_t glSide) const { return mapIndex[glSide-1]; }

protected:

    gsMultiPatch<T> m_patchRotated;

    gsMultiBasis<T> m_basisRotated;

    // Referenz to initial geometry
    gsVector<index_t> mapIndex;

    // Stores the changing of the axis
    // 0 -> axis not changed
    // 1 -> axis swapped (x, y --> y, x)
    bool axisOrientation;

    // How many rotation of the axis has been executed
    // Positive -> anticlockwise    Negative -> clockwise
    index_t rotationNum;

};

}
