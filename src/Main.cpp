
#include "DDC.h"

/*************************************/
/*       Main DDCcore function       */
/*************************************/
// [[Rcpp::export]]
Rcpp::List DDC_cpp(arma::mat & X, const double & tolProbCell,
                   const double & tolProbRow, const double & tolProbReg,
                   const double & tolProbCorr , const double & corrlim,
                   const int & combinRule,const int & rowdetect,
                   const int & includeSelf, const int & fastDDC,
                   const int & absCorr, const int & qdim, const int & transFun,
                   const int & treetype, const int & searchtype, 
                   const double & radius, const double & eps,
                   const int & bruteForce, unsigned int & k,
                   const unsigned int & numiter, const double & precScale)
{
  
  try
  {
    
    const double qCell     = std::sqrt(R::qchisq(tolProbCell, 1,true,false));
    const double qRow      = std::sqrt(R::qchisq(tolProbRow, 1,true,false));
    const double qRegr     = std::sqrt(R::qchisq(tolProbReg, 1,true,false));
    const double qCorr     = R::qchisq(tolProbCorr, 2,true,false);
    
    LocScaleEstimators::Xlocscale locscaleX;
    arma::mat Z = X;
    arma::mat Zest;
    arma::mat Zres;
    arma::uvec indcells;

    
    /////////////////////////////////////
    //    STEP 1: STANDARDIZE DATA     //
    /////////////////////////////////////
    
    // Robust standardization
    
    if (fastDDC == 0) {
      locscaleX = LocScaleEstimators::estLocScale(X, 0, precScale);
    } else {
      locscaleX = LocScaleEstimators::estLocScale(X, transFun, precScale);
    }
    Z = X.each_row() - locscaleX.loc.t();
    Z = Z.each_row() / locscaleX.scale.t();
    
 
    /////////////////////////////////////
    //    STEP 2: UNIVARIATE ANALYSIS  //
    /////////////////////////////////////
    
    arma::uvec indNAs = arma::find_nonfinite(X);
    arma::mat U = Z;
    // Univariate analysis : replace outliers by NAs
    
    U.for_each([qCell](arma::mat::elem_type &value) {
      value = std::abs(value) > qCell ? arma::datum::nan : value;
    });
    
    // for (unsigned int i = 0; i < U.n_cols; i++)
    // {
    //   U.col(i) = limitFilt(Z.col(i), qCell);
    // }
    
    arma::uvec UniIndex = DDC::vdiff(arma::find_nonfinite(U), indNAs); //does not include original missings
    

    
    /////////////////////////////////////////////////////
    //    STEP 3: CALCULATE CORRELATIONS AND SLOPES    //
    /////////////////////////////////////////////////////
    
    k = k > (X.n_cols-1) ? (X.n_cols-1) : k; //   - 1 since we do not have to compute the correlation of a column with itself
    
    // For each column j of U, find the k columns h != j of U that
    // it has the highest absolute correlation robCorr with :
    arma::umat ngbrs(U.n_cols, k, arma::fill::zeros);
    arma::mat robcors(U.n_cols, k, arma::fill::zeros);
    
    if (fastDDC == 0) {
      for (unsigned int i = 0; i < U.n_cols; i++)
      {
        DDC::kbestcorr tempresult = DDC::kBestCorr(U.col(i), U, i, k, qCorr, precScale);
        ngbrs.row(i) = tempresult.selected.t();
        robcors.row(i) = tempresult.corrs.t();
      }
    } else {
      
      DDC::fastRobCorout nn2result = DDC::FastRobCorActual(X, locscaleX.loc, 
                                                           locscaleX.scale, k,
                                                           qdim, absCorr,
                                                           transFun, precScale, 
                                                           bruteForce, treetype,
                                                           searchtype, radius,
                                                           eps, 0);
      
      ngbrs = nn2result.ngbrs;
      robcors = nn2result.robcorrs;
    }
    
   
   arma::mat corrweight = arma::abs(robcors); // should have no NAs
    
    if (corrlim > 0) {
      corrweight(arma::find(corrweight < corrlim)).zeros();
    }
    
    arma::umat ngb0 = ngbrs;
    
    ngb0.elem(find(corrweight == 0)).fill(U.n_cols); // out of bounds index for the unused ngbrs 
    arma::mat robslopes(U.n_cols, k, arma::fill::zeros);
    
    for (unsigned int i = 0; i < U.n_cols; i++) 
    {
      robslopes.row(i) = DDC::compSlopes(U.col(i), ngb0.row(i).t(), U, qRegr, precScale).t();
    }
    
    arma::uvec colStandalone = arma::find(sum(corrweight, 1) == 0);
    arma::uvec colConnected = DDC::vdiff(arma::regspace<arma::uvec>(0, X.n_cols - 1), colStandalone);
    arma::uvec indexStandalone = DDC::col2cell(colStandalone, X.n_rows);
    indexStandalone = DDC::vinter(indexStandalone, UniIndex);
    //= list of flagged cells in standalone variables.
    if (includeSelf == 1) {
      // if you want to include column j in its own prediction:
      ngbrs = arma::join_rows(arma::regspace<arma::uvec>(0, (X.n_cols - 1)), ngbrs);
      robcors = arma::join_rows(arma::ones<arma::vec>(X.n_cols), robcors);
      corrweight = arma::join_rows(arma::ones<arma::vec>(X.n_cols), corrweight);
      robslopes = arma::join_rows(arma::ones<arma::vec>(X.n_cols), robslopes);
    }
    
    
  
    for (unsigned int iter = 0; iter < numiter; iter++) {
      
      ////////////////////////////////////
      //    STEP 4 : ESTIMATE CELLS     //
      ////////////////////////////////////
      
      Zest = U; // These values will remain for standalone columns.
      
      // Estimation for connected variables :
      
      for (unsigned int i = 0; i < colConnected.size(); i++) {
        Zest.col(colConnected(i)) = DDC::predictCol(U.col(colConnected(i)), U, colConnected(i),
                 ngbrs, corrweight, robslopes, combinRule);
      }
    
      ////////////////////////////////////
      //    STEP 5 : DESHRINKAGE        //
      ////////////////////////////////////
      
      // Deshrinkage : rescale Zest[, j] using robSlope of Z[, j] on Zest[, j]
      
      for (unsigned int i = 0; i < colConnected.size(); i++) {
        Zest.col(colConnected(i)) = DDC::deShrink(Zest.col(colConnected(i)),
                 Z, colConnected(i), qRegr, precScale);
      }
      
      // Finally, all NAs are replaced by zeroes :
      Zest(find_nonfinite(Zest)).zeros();
      
      ////////////////////////////////////
      //    STEP 6 : FLAGGING CELLS     //
      ////////////////////////////////////
      
      // Compute cell residuals :
      Zres = Z - Zest; // original minus estimated
      Zres.cols(colStandalone) = Z.cols(colStandalone);
      
      arma::vec scalest(colConnected.size(), arma::fill::zeros);
      
      for (unsigned int i = 0; i < colConnected.size(); i++)
      {
        scalest(i) = LocScaleEstimators::scale1StepM(Zres.col(colConnected(i)), 
                LocScaleEstimators::rhoHuber15,
                arma::datum::nan, precScale);
        Zres.col(colConnected(i)) = Zres.col(colConnected(i)) / scalest(i);
      }
      
      // We don't have to scale the standalone columns, as they
      // were already standardized in the beginning.
      // Next, flag outlying cells by their large residuals :
      indcells = find(arma::abs(Zres) > qCell); // does not flag the NAs as cells
      U(indcells).fill(arma::datum::nan);
      
    } // ends the iteration
    
    indcells = DDC::vdiff(indcells, DDC::col2cell(colStandalone, X.n_cols));
    
    // are the indices of outlying cells in connected variables only
    indcells = arma::unique(arma::sort(arma::join_cols(indcells, indexStandalone)));
    // are the indices of both types of outlying cells
    
    ////////////////////////////////////
    //    STEP 7 : FLAGGING ROWS      //
    ////////////////////////////////////
    
    arma::vec Ti(X.n_rows, arma::fill::zeros);
    arma::uvec indrows;
    arma::uvec indall = indcells;
    
    if (rowdetect == 1) {
      
      for (unsigned int i = 0; i < X.n_rows; i++) {
        arma::vec tempres = arma::erf(arma::sqrt(arma::pow(Zres.row(i), 2) / 2)).t();
        //pchisq(x,1) = erf(sqrt(x/2))
        Ti(i) = arma::mean(tempres(arma::find_finite(tempres))) - 0.5;
      }
      
      // calculate the test value(outlyingness) of each row :
      arma::uvec finiteTi = arma::find_finite(Ti);
      double medTi = arma::median(Ti(finiteTi));
      Ti = (Ti - medTi) / (1.4826 * median(arma::abs(Ti(finiteTi) - medTi)));
      indrows = DDC::vinter(find_nonfinite(DDC::limitFilt(Ti, qRow)), arma::find(Ti > 0));
      indall = arma::unique(arma::join_cols(indcells, DDC::row2cell(indrows, X.n_rows, X.n_cols)));
    }
    
    /////////////////////////////////////////////
    //    STEP 8: UNSTANDARDIZE AND IMPUTE     //
    /////////////////////////////////////////////
    
    // compute Xest(storing in the existing matrix Zest to save space)
    Zest = Zest.each_row() % locscaleX.scale.t();
    Zest = Zest.each_row() + locscaleX.loc.t();
    
    // compute Ximp(storing in the existing matrix X to save space)
    X(indcells) = Zest(indcells);  // imputes the outlying cells of X
    X(indNAs) = Zest(indNAs); // imputes the missing values of X
    
    // for conversion to R-indexing
    indcells = indcells +1;
    indNAs   = indNAs +1;
    indrows  = indrows +1;
    ngbrs    = ngbrs + 1;
    indall   = indall +1;
    
    
    return Rcpp::List::create( Rcpp::Named("k") = k,
                               Rcpp::Named("ngbrs") = ngbrs,
                               Rcpp::Named("robcors") = robcors,
                               Rcpp::Named("robslopes") = robslopes,
                               Rcpp::Named("Xest") = Zest,
                               Rcpp::Named("stdResid") = Zres,
                               Rcpp::Named("indcells") = indcells,
                               Rcpp::Named("Ti") = Ti,
                               Rcpp::Named("indrows") = indrows,
                               Rcpp::Named("indNAs") = indNAs,
                               Rcpp::Named("indall") = indall,
                               Rcpp::Named("Ximp") = X,
                               Rcpp::Named("Z") = Z);
    
  } catch( std::exception& __ex__ )
  {
    forward_exception_to_r( __ex__ );
  } catch(...)
  {
    ::Rf_error( "c++ exception " "(unknown reason)" );
  }
  return Rcpp::wrap(NA_REAL);
}




/**********************************/
/*       Main Wrap function       */
/**********************************/


// [[Rcpp::export]]
Rcpp::List Wrap_cpp(arma::mat & X, arma::vec & loc, arma::vec & scale, double precScale) {
  try
  {
    
    arma::mat Xw = X;
    
    for (unsigned int i = 0; i < X.n_cols; i++) {
      arma::uvec finiteinds = arma::find_finite(X.col(i));
      arma::vec u = X.col(i) - loc(i);
      u = u / scale(i);
      arma::vec ufin = u(finiteinds);
      LocScaleEstimators::psiTanh(ufin);
      u(finiteinds) = ufin * scale(i) + loc(i) ;
      if (finiteinds.size() < X.n_rows) {
        arma::uvec infiniteinds = DDC::vdiff(arma::regspace<arma::uvec>(0,(X.n_rows - 1)), finiteinds);
        u(infiniteinds).fill(loc(i));
      }
      Xw.col(i) = u;
    } 
    
    return(Rcpp::List::create( Rcpp::Named("Xw") = Xw));
  } catch( std::exception& __ex__ )
  {
    forward_exception_to_r( __ex__ );
  } catch(...)
  {
    ::Rf_error( "c++ exception " "(unknown reason)" );
  }
  
  return Rcpp::wrap(NA_REAL);
}


/*****************************************/
/*       Main estLocScale function       */
/*****************************************/


// [[Rcpp::export]]
Rcpp::List estLocScale_cpp(arma::mat & X, int type,  double precScale) {
  // type 0 = biweight location + huber 1.5 scale
  // type 1 = huber location + huber scale
  // type 2 =  reweighted unimcd + wrapping location
  // type 3 =  reweighted unimcd
  // type 4 = reweighted raw mcd
  // type 5 = median/mad + wrapping 1 step M
  try
  {
    LocScaleEstimators::Xlocscale locscaleX;
    locscaleX = LocScaleEstimators::estLocScale(X, type, precScale); 
    
    return(Rcpp::List::create( Rcpp::Named("loc") = locscaleX.loc,
                               Rcpp::Named("scale") = locscaleX.scale));
    
  } catch( std::exception& __ex__ )
  {
    forward_exception_to_r( __ex__ );
  } catch(...)
  {
    ::Rf_error( "c++ exception " "(unknown reason)" );
  }
  
  return Rcpp::wrap(NA_REAL);
}


