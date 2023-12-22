/** @file gsC1SurfVertex.h

    @brief Creates the (approx.) C1 Vertex space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Farahat
*/

#pragma once


#include <gsUnstructuredSplines/src/gsC1SurfGluingData.h>




namespace gismo
{
    template <class T>
    class gsG1ASVisitorBasisVertex
    {
    public:

        gsG1ASVisitorBasisVertex()
        {
        }

        void initialize(const gsBasis<T>       & basis, //
                        gsQuadRule<T>    & rule)
        {
            gsVector<index_t> numQuadNodes( basis.dim() );
            for (index_t i = 0; i < basis.dim(); ++i) // to do: improve
                numQuadNodes[i] = basis.degree(i) + 1;

            // Setup Quadrature
            rule = gsGaussRule<T>(numQuadNodes);// NB!

            // Set Geometry evaluation flags
            md.flags = NEED_MEASURE ;

            localMat.resize(6);
            localRhs.resize(6);

            rhsVals.resize(6);
        }

        // Evaluate on element.
        inline void evaluate(gsBasis<T>       & basis, //
                             gsBasis<T>       & basis_geo,
                             std::vector<gsBSplineBasis<T>>       & basis_plus,
                             std::vector<gsBSplineBasis<T>>      & basis_minus,
                             const gsGeometry<T>    & geo, // patch
                             gsMatrix<T>            & quNodes,
                             gsMatrix<T>  & gluingData,
                             std::vector<bool> & isBoundary,
                             gsMatrix<T>  & Phi)
        {
            md.points = quNodes;

            // Compute the active basis functions
            // Assumes actives are the same for all quadrature points on the elements
            basis.active_into(md.points.col(0), actives);

            // Evaluate basis functions on element
            basis.eval_into(md.points, basisData);

            // Compute geometry related values
            geo.computeMap(md);


            numActive = actives.rows();

            // Computing c, c+ and c-
            std::vector<gsMatrix<T>> c_0, c_1;
            std::vector<gsMatrix < >> c_0_plus, c_1_plus, c_2_plus;
            std::vector<gsMatrix < >> c_0_plus_deriv, c_1_plus_deriv, c_2_plus_deriv;
            std::vector<gsMatrix < >> c_0_minus, c_1_minus;
            for (index_t i = 0; i < 2; i++) // i == 0 == u , i == 1 == v
            {
                gsMatrix<T> b_0, b_1;
                gsMatrix<T> b_0_plus, b_1_plus, b_2_plus;
                gsMatrix<T> b_0_plus_deriv, b_1_plus_deriv, b_2_plus_deriv;
                gsMatrix<T> b_0_minus, b_1_minus;

//                gsBSplineBasis<T> bsp_temp = dynamic_cast<gsBSplineBasis<T> & >(basis_geo.component(i));
//                T p = bsp_temp.maxDegree();
//                T h_geo = bsp_temp.knots().at(p + 2);

                basis_geo.component(i).evalSingle_into(0, md.points.row(i),b_0); // first
                basis_geo.component(i).evalSingle_into(1, md.points.row(i),b_1); // second

                basis_plus[i].evalSingle_into(0, md.points.row(i),b_0_plus);
                basis_plus[i].evalSingle_into(1, md.points.row(i),b_1_plus);
                basis_plus[i].evalSingle_into(2, md.points.row(i),b_2_plus);

                basis_plus[i].derivSingle_into(0, md.points.row(i),b_0_plus_deriv);
                basis_plus[i].derivSingle_into(1, md.points.row(i),b_1_plus_deriv);
                basis_plus[i].derivSingle_into(2, md.points.row(i),b_2_plus_deriv);

                basis_minus[i].evalSingle_into(0, md.points.row(i),b_0_minus);
                basis_minus[i].evalSingle_into(1, md.points.row(i),b_1_minus);

                // Point zero
                gsMatrix<T> zero;
                zero.setZero(2,1);

                gsMatrix<T> b_1_0, b_1_minus_0;
                basis_geo.component(i).derivSingle_into(1, zero.row(i),b_1_0);
                basis_minus[i].derivSingle_into(1, zero.row(i),b_1_minus_0);
//                c_0.push_back(b_0 + b_1);
//                c_1.push_back((h_geo / p) * b_1);
//
//                c_0_minus.push_back(b_0_minus + b_1_minus);
//                c_1_minus.push_back(h_geo/ (p-1) * b_1_minus);
                T factor_b_1 = 1.0/b_1_0(0,0);
                c_0.push_back(b_0 + b_1);
                c_1.push_back(factor_b_1 * b_1);

                T factor_b_1_minus = 1.0/b_1_minus_0(0,0);
                c_0_minus.push_back(b_0_minus + b_1_minus);
                c_1_minus.push_back(factor_b_1_minus * b_1_minus);

                gsMatrix<T> der_b_1_plus_0, der2_b_1_plus_0, der2_b_2_plus_0;
                basis_plus[i].derivSingle_into(1, zero.row(i), der_b_1_plus_0);
                basis_plus[i].deriv2Single_into(1, zero.row(i), der2_b_1_plus_0);
                basis_plus[i].deriv2Single_into(2, zero.row(i), der2_b_2_plus_0);

                T factor_c_1_plus = 1/der_b_1_plus_0(0,0);
                T factor2_c_1_plus = -der2_b_1_plus_0(0,0)/(der_b_1_plus_0(0,0)*der2_b_2_plus_0(0,0));
                T factor_c_2_plus = 1/der2_b_2_plus_0(0,0);

                c_0_plus.push_back(b_0_plus + b_1_plus + b_2_plus);
                c_1_plus.push_back(factor_c_1_plus * b_1_plus + factor2_c_1_plus * b_2_plus);
                c_2_plus.push_back(factor_c_2_plus * b_2_plus );

                c_0_plus_deriv.push_back(b_0_plus_deriv + b_1_plus_deriv + b_2_plus_deriv);
                c_1_plus_deriv.push_back(factor_c_1_plus * b_1_plus_deriv + factor2_c_1_plus * b_2_plus_deriv);
                c_2_plus_deriv.push_back(factor_c_2_plus * b_2_plus_deriv);
            }

//        if (g1OptionList.getInt("gluingData") == gluingData::global)
//        {
//            alpha.push_back(gluingData[0].get_alpha_tilde().eval(md.points.row(0))); // u
//            alpha.push_back(gluingData[1].get_alpha_tilde().eval(md.points.row(1))); // v
//            alpha_0.push_back(gluingData[0].get_alpha_tilde().eval(zero.row(0))); // u
//            alpha_0.push_back(gluingData[1].get_alpha_tilde().eval(zero.row(0))); // v
//            alpha_deriv.push_back(gluingData[0].get_alpha_tilde().deriv(zero.row(0))); // u
//            alpha_deriv.push_back(gluingData[1].get_alpha_tilde().deriv(zero.row(0))); // v
//
//            beta.push_back(gluingData[0].get_beta_tilde().eval(md.points.row(0))); // u
//            beta.push_back(gluingData[1].get_beta_tilde().eval(md.points.row(1))); // v
//            beta_0.push_back(gluingData[0].get_beta_tilde().eval(zero.row(0))); // u
//            beta_0.push_back(gluingData[1].get_beta_tilde().eval(zero.row(0))); // v
//            beta_deriv.push_back(gluingData[0].get_beta_tilde().deriv(zero.row(0))); // u
//            beta_deriv.push_back(gluingData[1].get_beta_tilde().deriv(zero.row(0))); // v
//
//        }
//        else if (g1OptionList.getInt("gluingData") == gluingData::local)
//        {
//            alpha.push_back(gluingData[0].get_local_alpha_tilde(0).eval(md.points.row(0))); // u
//            alpha.push_back(gluingData[1].get_local_alpha_tilde(0).eval(md.points.row(1))); // v
//            alpha_0.push_back(gluingData[0].get_local_alpha_tilde(0).eval(zero.row(0))); // u
//            alpha_0.push_back(gluingData[1].get_local_alpha_tilde(0).eval(zero.row(0))); // v
//            alpha_deriv.push_back(gluingData[0].get_local_alpha_tilde(0).deriv(zero.row(0))); // u
//            alpha_deriv.push_back(gluingData[1].get_local_alpha_tilde(0).deriv(zero.row(0))); // v
//
//            beta.push_back(gluingData[0].get_local_beta_tilde(0).eval(md.points.row(0))); // u
//            beta.push_back(gluingData[1].get_local_beta_tilde(0).eval(md.points.row(1))); // v
//            beta_0.push_back(gluingData[0].get_local_beta_tilde(0).eval(zero.row(0))); // u
//            beta_0.push_back(gluingData[1].get_local_beta_tilde(0).eval(zero.row(0))); // v
//            beta_deriv.push_back(gluingData[0].get_local_beta_tilde(0).deriv(zero.row(0))); // u
//            beta_deriv.push_back(gluingData[1].get_local_beta_tilde(0).deriv(zero.row(0))); // v
//        }

            std::vector<gsMatrix<T>> alpha, beta, alpha_0, beta_0, alpha_deriv, beta_deriv;

            // Point zero
            gsMatrix<T> zero;
            zero.setZero(2,1);

            // Geo data:
            gsMatrix<T> geo_jac, geo_der2;
            geo_jac = geo.jacobian(zero);
            geo_der2 = geo.deriv2(zero);

            gsMatrix<T> zeros(1, md.points.cols());
            zeros.setZero();

            // Point One
            gsMatrix<T> one;
            one.setOnes(2,1);

            gsMatrix<T> ones(1, md.points.cols());
            ones.setOnes();

            alpha.push_back( gluingData(0, 0) * ( ones - md.points.row(0) ) + gluingData(1, 0) * md.points.row(0) ); // u
            alpha.push_back( gluingData(0, 1) * ( ones - md.points.row(1) ) + gluingData(1, 1) * md.points.row(1) ); // v


            alpha_0.push_back( gluingData(0, 0) * ( one.row(0) - zero.row(0) ) + gluingData(1, 0) * zero.row(0) ); // u
            alpha_0.push_back( gluingData(0, 1) * ( one.row(0) - zero.row(0) ) + gluingData(1, 1) * zero.row(0) ); // v
            alpha_deriv.push_back( ( gluingData(1, 0) - gluingData(0, 0) ) * ones.col(0) ); // u
            alpha_deriv.push_back( ( gluingData(1, 1) - gluingData(0, 1) ) * ones.col(0) ); // v

            beta.push_back( gluingData(2, 0) * ( ones - md.points.row(0) ) + gluingData(3, 0) * md.points.row(0) ); // u
            beta.push_back( gluingData(2, 1) * ( ones - md.points.row(1) ) + gluingData(3, 1) * md.points.row(1) ); // v
            beta_0.push_back( gluingData(2, 0) * ( one.row(0) - zero.row(0) ) + gluingData(3, 0) * zero.row(0) ); // u
            beta_0.push_back( gluingData(2, 1) * ( one.row(0) - zero.row(0) ) + gluingData(3, 1) * zero.row(0) ); // v
            beta_deriv.push_back( ( gluingData(3, 0) - gluingData(2, 0) ) * ones.col(0) ); // u
            beta_deriv.push_back( ( gluingData(3, 1) - gluingData(2, 1) ) * ones.col(0) ); // v

            // Compute dd^^(i_k) and dd^^(i_k-1)
            gsMatrix<T> dd_ik_plus, dd_ik_minus;
            gsMatrix<T> dd_ik_minus_deriv, dd_ik_plus_deriv;

            dd_ik_minus = ( -1 / alpha_0[0](0,0) ) * ( geo_jac.col(1) +
                                                       beta_0[0](0,0) * geo_jac.col(0) );

            dd_ik_plus = ( 1 / alpha_0[1](0,0) ) * ( geo_jac.col(0) +
                                                     beta_0[1](0,0) * geo_jac.col(1) );

            gsMatrix<T> geo_deriv2_12(geo.targetDim(), 1), geo_deriv2_11(geo.targetDim(), 1), geo_deriv2_22(geo.targetDim(), 1);

            geo_deriv2_12.row(0) = geo_der2.row(2);
            geo_deriv2_12.row(1) = geo_der2.row(5);

            geo_deriv2_11.row(0) = geo_der2.row(0);
            geo_deriv2_11.row(1) = geo_der2.row(3);

            geo_deriv2_22.row(0) = geo_der2.row(1);
            geo_deriv2_22.row(1) = geo_der2.row(4);

            if(geo.parDim() + 1 == geo.targetDim()) // In the surface case the dimension of the second derivative vector is 9x1
            {
                geo_deriv2_12.row(2) = geo_der2.row(8);

                geo_deriv2_11.row(2) = geo_der2.row(6);

                geo_deriv2_22.row(2) = geo_der2.row(7);
            }

            gsMatrix<T> alpha_squared_u = alpha_0[0]*alpha_0[0];
            gsMatrix<T> alpha_squared_v = alpha_0[1]*alpha_0[1];

            dd_ik_minus_deriv = -1/(alpha_squared_u(0,0)) * // N^2
                                ( ( geo_deriv2_12 + (beta_deriv[0](0,0) * geo_jac.col(0) + beta_0[0](0,0) * geo_deriv2_11) ) * alpha_0[0](0,0) -
                                  ( geo_jac.col(1) + beta_0[0](0,0) * geo_jac.col(0) ) * alpha_deriv[0](0,0) );


            dd_ik_plus_deriv = 1/(alpha_squared_v(0,0)) *
                               ( ( geo_deriv2_12 + (beta_deriv[1](0,0) * geo_jac.col(1) + beta_0[1](0,0) * geo_deriv2_22) ) * alpha_0[1](0,0) -
                                 ( geo_jac.col(0) + beta_0[1](0,0) * geo_jac.col(1) ) * alpha_deriv[1](0,0) );

            //if (isBoundary[0] == false)
            //    gsInfo << dd_ik_minus_deriv << "\n";
            //if (isBoundary[1] == false)
            //    gsInfo << dd_ik_plus_deriv << "\n";

            // Comupute d_(0,0)^(i_k), d_(1,0)^(i_k), d_(0,1)^(i_k), d_(1,1)^(i_k) ; i_k == 2
            std::vector<gsMatrix<T>> d_ik;
            d_ik.push_back(Phi.row(0).transpose());

            d_ik.push_back(Phi.block(1, 0, geo.targetDim(), 6).transpose() * geo_jac.col(0) ); // deriv into u


//        gsInfo << "d_ik: " << d_ik.back()  << "\n";
//        gsInfo << "======================================= \n";
//        gsInfo << "geo_jac.col(0): " << geo_jac.col(0)  << "\n";
//        gsInfo << "======================================= \n";


            d_ik.push_back(Phi.block(1, 0, geo.targetDim(), 6).transpose() * geo_jac.col(1) ); // deriv into v

            // Hessian
            if(geo.parDim() + 1 == geo.targetDim()) // In the surface case the dimension of the second derivative vector is 9x1
            {
                d_ik.push_back( (geo_jac(0,0) * Phi.row(4).transpose() + geo_jac(1,0) * Phi.row(7).transpose() + geo_jac(2,0) * Phi.row(10).transpose()) * geo_jac(0,1) +
                                (geo_jac(0,0) * Phi.row(5).transpose() + geo_jac(1,0) * Phi.row(8).transpose() + geo_jac(2,0) * Phi.row(11).transpose()) * geo_jac(1,1) +
                                (geo_jac(0,0) * Phi.row(6).transpose() + geo_jac(1,0) * Phi.row(9).transpose() + geo_jac(2,0) * Phi.row(12).transpose()) * geo_jac(2,1) +
                                Phi.block(1, 0, 1, 6).transpose() * geo_der2.row(2) +
                                Phi.block(2, 0, 1, 6).transpose() * geo_der2.row(5) +
                                Phi.block(3, 0, 1, 6).transpose() * geo_der2.row(8) );

            }
            else
            {
                d_ik.push_back( (geo_jac(0,0) * Phi.col(3) + geo_jac(1,0) * Phi.col(4)) * geo_jac(0,1) +
                                (geo_jac(0,0) * Phi.col(4) + geo_jac(1,0) * Phi.col(5)) * geo_jac(1,1) +
                                Phi.block(0,1, 6,1) * geo_der2.row(2) +
                                Phi.block(0,2, 6,1) * geo_der2.row(5)); // Hessian
            }

//        gsInfo << "d_ik: " << d_ik.back() << " : " << Phi << "\n";
            // Compute d_(*,*)^(il,ik)
            std::vector<gsMatrix<T>> d_ilik_minus, d_ilik_plus;

//      d_(*,*)^(ik-1,ik)
            d_ilik_minus.push_back(Phi.row(0).transpose());

            d_ilik_minus.push_back(Phi.block(1, 0, geo.targetDim(), 6).transpose() * geo_jac.col(0));

            if(geo.parDim() + 1 == geo.targetDim()) // In the surface case the dimension of the second derivative vector is 9x1
            {
                d_ilik_minus.push_back( (geo_jac(0,0) * Phi.row(4).transpose() + geo_jac(1,0) * Phi.row(7).transpose() + geo_jac(2,0) * Phi.row(10).transpose()) * geo_jac(0,0) +
                                        (geo_jac(0,0) * Phi.row(5).transpose() + geo_jac(1,0) * Phi.row(8).transpose() + geo_jac(2,0) * Phi.row(11).transpose()) * geo_jac(1,0) +
                                        (geo_jac(0,0) * Phi.row(6).transpose() + geo_jac(1,0) * Phi.row(9).transpose() + geo_jac(2,0) * Phi.row(12).transpose()) * geo_jac(2,0) +
                                        Phi.block(1, 0, 1, 6).transpose() * geo_der2.row(0) +
                                        Phi.block(2, 0, 1, 6).transpose() * geo_der2.row(3) +
                                        Phi.block(3, 0, 1, 6).transpose() * geo_der2.row(6) );
            }
            else
            {
                d_ilik_minus.push_back( (geo_jac(0,0) * Phi.col(3) + geo_jac(1,0) * Phi.col(4))*geo_jac(0,0) +
                                        (geo_jac(0,0) * Phi.col(4) + geo_jac(1,0) * Phi.col(5))*geo_jac(1,0) +
                                        Phi.block(0,1,6,1) * geo_der2.row(0) +
                                        Phi.block(0,2,6,1) * geo_der2.row(3));
            }

            d_ilik_minus.push_back(Phi.block(1, 0, geo.targetDim(), 6).transpose() * dd_ik_minus);

            if(geo.parDim() + 1 == geo.targetDim()) // In the surface case the dimension of the second derivative vector is 9x1
            {
                d_ilik_minus.push_back( (geo_jac(0,0) * Phi.row(4).transpose() + geo_jac(1,0) * Phi.row(7).transpose() + geo_jac(2,0) * Phi.row(10).transpose()) * dd_ik_minus(0,0) +
                                        (geo_jac(0,0) * Phi.row(5).transpose() + geo_jac(1,0) * Phi.row(8).transpose() + geo_jac(2,0) * Phi.row(11).transpose()) * dd_ik_minus(1,0) +
                                        (geo_jac(0,0) * Phi.row(6).transpose() + geo_jac(1,0) * Phi.row(9).transpose() + geo_jac(2,0) * Phi.row(12).transpose()) * dd_ik_minus(2,0) +
                                        Phi.block(1, 0, 1, 6).transpose() * dd_ik_minus_deriv.row(0) +
                                        Phi.block(2, 0, 1, 6).transpose() * dd_ik_minus_deriv.row(1) +
                                        Phi.block(3, 0, 1, 6).transpose() * dd_ik_minus_deriv.row(2) );
            }
            else
            {
                d_ilik_minus.push_back( (geo_jac(0,0) * Phi.col(3) + geo_jac(1,0) * Phi.col(4)) * dd_ik_minus(0,0) +
                                        (geo_jac(0,0) * Phi.col(4) + geo_jac(1,0) * Phi.col(5)) * dd_ik_minus(1,0) +
                                        Phi.block(0,1,6,1) * dd_ik_minus_deriv.row(0) +
                                        Phi.block(0,2,6,1) * dd_ik_minus_deriv.row(1));
            }


//      d_(*,*)^(ik+1,ik)
            d_ilik_plus.push_back(Phi.row(0).transpose());

            d_ilik_plus.push_back(Phi.block(1, 0, geo.targetDim(), 6).transpose() * geo_jac.col(1));

            if(geo.parDim() + 1 == geo.targetDim()) // In the surface case the dimension of the second derivative vector is 9x1
            {
                d_ilik_plus.push_back( (geo_jac(0,1) * Phi.row(4).transpose() + geo_jac(1,1) * Phi.row(7).transpose() + geo_jac(2,1) * Phi.row(10).transpose()) * geo_jac(0,1) +
                                       (geo_jac(0,1) * Phi.row(5).transpose() + geo_jac(1,1) * Phi.row(8).transpose() + geo_jac(2,1) * Phi.row(11).transpose()) * geo_jac(1,1) +
                                       (geo_jac(0,1) * Phi.row(6).transpose() + geo_jac(1,1) * Phi.row(9).transpose() + geo_jac(2,1) * Phi.row(12).transpose()) * geo_jac(2,1) +
                                       Phi.block(1, 0, 1, 6).transpose() * geo_der2.row(1) +
                                       Phi.block(2, 0, 1, 6).transpose() * geo_der2.row(4) +
                                       Phi.block(3, 0, 1, 6).transpose() * geo_der2.row(7) );
            }
            else
            {
                d_ilik_plus.push_back( (geo_jac(0,1) * Phi.col(3) + geo_jac(1,1) * Phi.col(4)) * geo_jac(0,1) +
                                       (geo_jac(0,1) * Phi.col(4) + geo_jac(1,1) * Phi.col(5)) * geo_jac(1,1) +
                                       Phi.block(0,1,6,1) * geo_der2.row(1) +
                                       Phi.block(0,2,6,1) * geo_der2.row(4) );
            }

            d_ilik_plus.push_back(Phi.block(1, 0, geo.targetDim(), 6).transpose() * dd_ik_plus);

            if(geo.parDim() + 1 == geo.targetDim()) // In the surface case the dimension of the second derivative vector is 9x1
            {
                d_ilik_plus.push_back( (geo_jac(0,1) * Phi.row(4).transpose() + geo_jac(1,1) * Phi.row(7).transpose() + geo_jac(2,1) * Phi.row(10).transpose()) * dd_ik_plus(0,0) +
                                       (geo_jac(0,1) * Phi.row(5).transpose() + geo_jac(1,1) * Phi.row(8).transpose() + geo_jac(2,1) * Phi.row(11).transpose()) * dd_ik_plus(1,0) +
                                       (geo_jac(0,1) * Phi.row(6).transpose() + geo_jac(1,1) * Phi.row(9).transpose() + geo_jac(2,1) * Phi.row(12).transpose()) * dd_ik_plus(2,0) +
                                       Phi.block(1, 0, 1, 6).transpose() * dd_ik_plus_deriv.row(0) +
                                       Phi.block(2, 0, 1, 6).transpose() * dd_ik_plus_deriv.row(1) +
                                       Phi.block(3, 0, 1, 6).transpose() * dd_ik_plus_deriv.row(2) );
            }
            else
            {
                d_ilik_plus.push_back( (geo_jac(0,1) * Phi.col(3) + geo_jac(1,1) * Phi.col(4)) * dd_ik_plus(0,0) +
                                       (geo_jac(0,1) * Phi.col(4) + geo_jac(1,1) * Phi.col(5)) * dd_ik_plus(1,0) +
                                       Phi.block(0,1,6,1) * dd_ik_plus_deriv.row(0) +
                                       Phi.block(0,2,6,1) * dd_ik_plus_deriv.row(1) );
            }


            for (index_t i = 0; i < 6; i++)
            {
                rhsVals.at(i) = d_ilik_minus.at(0)(i,0) * (c_0_plus.at(0).cwiseProduct(c_0.at(1)) -
                                                           beta[0].cwiseProduct(c_0_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) +
                                d_ilik_minus.at(1)(i,0) * (c_1_plus.at(0).cwiseProduct(c_0.at(1)) -
                                                           beta[0].cwiseProduct(c_1_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) +
                                d_ilik_minus.at(2)(i,0) * (c_2_plus.at(0).cwiseProduct(c_0.at(1)) -
                                                           beta[0].cwiseProduct(c_2_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) -
                                d_ilik_minus.at(3)(i,0) * alpha[0].cwiseProduct(c_0_minus.at(0).cwiseProduct(c_1.at(1))) -
                                d_ilik_minus.at(4)(i,0) * alpha[0].cwiseProduct(c_1_minus.at(0).cwiseProduct(c_1.at(1))); // f*_(ik-1,ik)

                rhsVals.at(i) += d_ilik_plus.at(0)(i,0) * (c_0_plus.at(1).cwiseProduct(c_0.at(0)) -
                                                           beta[1].cwiseProduct(c_0_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
                                 d_ilik_plus.at(1)(i,0) * (c_1_plus.at(1).cwiseProduct(c_0.at(0)) -
                                                           beta[1].cwiseProduct(c_1_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
                                 d_ilik_plus.at(2)(i,0) * (c_2_plus.at(1).cwiseProduct(c_0.at(0)) -
                                                           beta[1].cwiseProduct(c_2_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
                                 d_ilik_plus.at(3)(i,0) * alpha[1].cwiseProduct(c_0_minus.at(1).cwiseProduct(c_1.at(0))) +
                                 d_ilik_plus.at(4)(i,0) * alpha[1].cwiseProduct(c_1_minus.at(1).cwiseProduct(c_1.at(0))); // f*_(ik+1,ik)

                rhsVals.at(i) -= d_ik.at(0)(i,0) * c_0.at(0).cwiseProduct(c_0.at(1)) + d_ik.at(2)(i,0) * c_0.at(0).cwiseProduct(c_1.at(1)) +
                                 d_ik.at(1)(i,0) * c_1.at(0).cwiseProduct(c_0.at(1)) + d_ik.at(3)(i,0) * c_1.at(0).cwiseProduct(c_1.at(1)); // f*_(ik)

                localMat.at(i).setZero(numActive, numActive);
                localRhs.at(i).setZero(numActive, rhsVals.at(i).rows());//multiple right-hand sides

                localMat.at(i).setZero(numActive, numActive);
                localRhs.at(i).setZero(numActive, rhsVals.at(i).rows());//multiple right-hand sides
            }



        } // evaluate

        inline void assemble(gsDomainIterator<T>    & element,
                             const gsVector<T>      & quWeights)
        {
            gsMatrix<T> & basisVals  = basisData;
            for (index_t i = 0; i < 6; i++)
            {
                // ( u, v)
                localMat.at(i).noalias() = basisData * quWeights.asDiagonal() * md.measures.asDiagonal() * basisData.transpose();

                for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
                {
                    T weight = quWeights[k];
                    gsMatrix<T> Jk = md.jacobian(k);

                    if( Jk.dim().second + 1 == Jk.dim().first)
                    {
                        gsMatrix<T> G = Jk.transpose() * Jk;
                        T detG = G.determinant();
                        weight *= sqrt(detG);
                    }
                    else
                    {
                        weight *=  md.measure(k);
                    }
                    // Multiply weight by the geometry measure

                    localRhs.at(i).noalias() += weight * (basisVals.col(k) * rhsVals.at(i).col(k).transpose());
                }
            }
        }

        inline void localToGlobal(const index_t patchIndex,
                                  const std::vector<gsMatrix<T> >    & eliminatedDofs,
                                  std::vector< gsSparseSystem<T> >     & system)
        {
            gsMatrix<index_t> actives_temp;
            for (size_t i = 0; i < system.size(); i++) // 6
            {
                // Map patch-local DoFs to global DoFs
                system.at(i).mapColIndices(actives, patchIndex, actives_temp);
                // Add contributions to the system matrix and right-hand side
                system.at(i).push(localMat.at(i), localRhs.at(i), actives_temp, eliminatedDofs[0], 0, 0);
            }
        }

    protected:
        gsMatrix<index_t> actives;
        gsMatrix<T> basisData;
        index_t numActive;

    protected:
        // Local values of the right hand side
        std::vector< gsMatrix<T> >  rhsVals;

    protected:
        // Local matrices
        std::vector< gsMatrix<T> > localMat;
        std::vector< gsMatrix<T> > localRhs;

        gsMapData<T> md;

    }; // class gsVisitorG1BasisVertex
} // namespace gismo
