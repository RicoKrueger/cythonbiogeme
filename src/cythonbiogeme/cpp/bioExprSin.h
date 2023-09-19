//-*-c++-*----------------------%--------------------------------------
//
// File name : bioExprSin.h
// @date   Tue Apr 17 12:12:43 2018
// @author Michel Bierlaire
// @version Revision 1.0
//
//--------------------------------------------------------------------

#ifndef bioExprSin_h
#define bioExprSin_h

#include "bioExpression.h"
#include "bioString.h"

class bioExprSin: public bioExpression {
 public:
  bioExprSin(bioExpression* c) ;
  ~bioExprSin() ;
  virtual const bioDerivatives* getValueAndDerivatives(std::vector<bioUInt> literalIds,
						 bioBoolean gradient,
						bioBoolean hessian) ;

  virtual bioString print(bioBoolean hp = false) const ;

 protected:
  bioExpression* child ;
};
#endif
