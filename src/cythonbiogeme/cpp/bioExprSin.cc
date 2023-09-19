//-*-c++-*------------------------------------------------------------
//
// File name : bioExprSin.cc
// @date   Tue Apr 17 12:13:00 2018
// @author Michel Bierlaire
// @version Revision 1.0
//
//--------------------------------------------------------------------

#include "bioExprSin.h"
#include "bioExceptions.h"
#include "bioDebug.h"
#include <cmath>
#include <sstream> 

bioExprSin::bioExprSin(bioExpression* c) :
  child(c) {

  listOfChildren.push_back(c) ;
}

bioExprSin::~bioExprSin() {

}

const bioDerivatives* bioExprSin::getValueAndDerivatives(std::vector<bioUInt> literalIds,
							 bioBoolean gradient,
							 bioBoolean hessian) {

  theDerivatives.with_g = gradient ;
  theDerivatives.with_h = hessian ;

  bioUInt n = literalIds.size() ;
  theDerivatives.resize(n) ;

  const bioDerivatives* childResult = child->getValueAndDerivatives(literalIds,gradient,hessian) ;
  bioReal cf = childResult->f ;
  theDerivatives.f = sin(cf) ;
  bioReal d1 = cos(cf) ;
  bioReal d2 = -sin(cf) ;

  if (gradient) {
    for (bioUInt i = 0 ; i < n ; ++i) {
      theDerivatives.g[i] = d1 * childResult->g[i] ;
      if (hessian) {
	for (bioUInt j = 0 ; j < n ; ++j) {
	  theDerivatives.h[i][j] =
	    (
	     d1 * childResult->h[i][j] +
	     d2 * childResult->g[i] * childResult->g[j]
	     );
	}
      }
    }
  }
  return &theDerivatives ;
}

bioString bioExprSin::print(bioBoolean hp) const {
  std::stringstream str ; 
  str << "sin(" << child->print(hp) << ")";
  return str.str() ;

}
