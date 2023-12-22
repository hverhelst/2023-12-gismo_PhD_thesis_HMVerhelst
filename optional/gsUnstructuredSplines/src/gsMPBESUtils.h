#pragma once

#include <gsHSplines/gsHTensorBasis.h>
#include <gsHSplines/gsTHBSplineBasis.h>

#include <gsUnstructuredSplines/src/gsMPBESBasis.h>
#include <gsUnstructuredSplines/src/gsMPBESHSplineBasis.h>
#include <gsUnstructuredSplines/src/gsMPBESBSplineBasis.h>

#include <gsIO/gsFileData.h>
#include <gsIO/gsReadFile.h>

namespace gismo
{

template<short_t d, class T>
gsMPBESBasis<d,T> * getCompBasisFromMultiPatch(const gsMultiPatch<T> & mp,index_t incrSmoothness = -1,index_t minEVDistance = -1 )
{
    gsMPBESBasis<d,T> * compBasis=NULL;
    bool tensorBSpline= true;
    bool hTensor = true;
    std::vector<gsTensorBSplineBasis<d,T>* >tensorBases;
    for(size_t i = 0;i<mp.nPatches();++i)
    {
        tensorBases.push_back(dynamic_cast<gsTensorBSplineBasis<d,T> * >(& mp.basis(i)));
        tensorBSpline = tensorBSpline && tensorBases[i]!=NULL;
    }
    if(tensorBSpline)
        compBasis = (new gsMPBESBSplineBasis<d,T>(tensorBases,mp,incrSmoothness,minEVDistance));
    else
    {
        std::vector<gsHTensorBasis<d,T>* >hBases;
        for(size_t i = 0;i<mp.nPatches();++i)
        {
            hBases.push_back(dynamic_cast<gsHTensorBasis<d,T> * >(& mp.basis(i)));
            hTensor = hTensor && hBases[i]!=NULL;
        }
        if(hTensor)
            compBasis = (new gsMPBESHSplineBasis<d,T>(hBases,mp,incrSmoothness,minEVDistance));
    }
    GISMO_ASSERT(tensorBSpline||hTensor,"No suitable basis for gsMappedBasis found.");
    return compBasis;
}

template<short_t d, class T>
gsMPBESBasis<d,T> * getCompBasisFromMultiPatch_withCoefs(const gsMultiPatch<T> & mp, std::vector<gsMatrix<T>* >&coefs,index_t incrSmoothness = -1,index_t minEVDistance = -1 )
{
    gsMPBESBasis<d,T> * compBasis=NULL;
    bool tensorBSpline= true;
    bool hTensor = true;
    std::vector<gsTensorBSplineBasis<d,T>* >tensorBases;
    for(size_t i = 0;i<mp.nPatches();++i)
    {
        tensorBases.push_back(dynamic_cast<gsTensorBSplineBasis<d,T> * >(& mp.basis(i)));
        tensorBSpline = tensorBSpline && tensorBases[i]!=NULL;
    }
    if(tensorBSpline)
        compBasis = (new gsMPBESBSplineBasis<d,T>(tensorBases,mp,coefs,incrSmoothness,minEVDistance));
    else
    {
        std::vector<gsHTensorBasis<d,T>* >hBases;
        for(size_t i = 0;i<mp.nPatches();++i)
        {
            hBases.push_back(dynamic_cast<gsHTensorBasis<d,T> * >(& mp.basis(i)));
            hTensor = hTensor && hBases[i]!=NULL;
        }
        if(hTensor)
            compBasis = (new gsMPBESHSplineBasis<d,T>(hBases,mp,coefs,incrSmoothness,minEVDistance));
    }
    GISMO_ASSERT(tensorBSpline||hTensor,"No suitable basis for gsMappedBasis found.");
    return compBasis;
}

} // end namespace gismo
