

#include "DDC.h"

/*************************************/
/*       Main DDCcore function       */
/*************************************/
// [[Rcpp::export]]
Rcpp::List DDC_cpp(arma::mat & X, const double & tolProbCell,
                   const double & tolProbRow, const double & tolProbReg,
                   const double & tolProbCorr , const double & corrlim,
                   const int & combinRule,
                   const int & includeSelf, const int & fastDDC,
                   const int & qdim, const int & transFun,
                   unsigned int & k,
                   const unsigned int & numiter, const double & precScale,
                   const int & standType, const int & corrType,
                   const unsigned int & nCorr, const unsigned int & nLocScale,
                   arma::uvec & goodCols)
{
  
  try
  {
    // hard coded options for NN search in fastDDC
    const int absCorr = 1;
    const int treetype = 0;
    const int searchtype = 1;
    const double radius = 0;
    const double eps = 0;
    
    const double qCell     = std::sqrt(R::qchisq(tolProbCell, 1,true,false));
    const double qRow      = std::sqrt(R::qchisq(tolProbRow, 1,true,false));
    const double qRegr     = std::sqrt(R::qchisq(tolProbReg, 1,true,false));
    const double qCorr     = R::qchisq(tolProbCorr, 2,true,false);
    
    LocScaleEstimators::Xlocscale locscaleX;
    arma::mat Z = X;
    arma::mat Zest;
    arma::mat Zres;
    arma::uvec indcells;
    arma::vec scalest;
    arma::vec deshrinkage;
    
    
    /////////////////////////////////////
    //    STEP 1: STANDARDIZE DATA     //
    /////////////////////////////////////
    
    // Robust standardization
    if (fastDDC == 0) {
      locscaleX = LocScaleEstimators::estLocScale(X, nLocScale,
                                                  standType, precScale, 1);
    } else {
      locscaleX = LocScaleEstimators::estLocScale(X, nLocScale,
                                                  2, precScale, 1);
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
    
    
    arma::uvec UniIndex = DDC::vdiff(arma::find_nonfinite(U), indNAs); //does not include original missings
    
    
    /////////////////////////////////////////////////////
    //    STEP 3: CALCULATE CORRELATIONS AND SLOPES    //
    /////////////////////////////////////////////////////
    
    
    // k = k > (X.n_cols - 1) ? (X.n_cols - 1) : k; //   - 1 since we do not have to compute the correlation of a column with itself
    k = k > (goodCols.size() - 1) ? (goodCols.size() - 1) : k; //   - 1 since we do not have to compute the correlation of a column with itself
    
    
    // For each column j of U, find the k columns h != j of U that
    // it has the highest absolute correlation robCorr with :
    arma::umat ngbrs(U.n_cols, k);
    ngbrs.fill(U.n_cols); // out of bound index used for inactive ngbrs
    arma::mat robcors(U.n_cols, k, arma::fill::zeros);
    
    
    if (fastDDC == 0) {
      switch(corrType) {
      case 1: {// wrap
      arma::mat Xw(U.n_rows, goodCols.size(), arma::fill::zeros);
      for (unsigned int i = 0; i < goodCols.size(); i++) {
        arma::vec u = Z.col(goodCols(i));
        arma::uvec finiteinds = arma::find_finite(u);
        arma::vec ufin = u(finiteinds);
        LocScaleEstimators::psiTanh(ufin);
        u.zeros();
        u(finiteinds) = ufin * locscaleX.scale(goodCols(i)); // no "+center": u filled with zeroes
        Xw.col(i) = u;
      }
      
      arma::mat corW = arma::cor(Xw);
      
      
      
      for (unsigned int i = 0; i < goodCols.size(); i++) {
        arma::vec tempcors = corW.col(i);
        arma::uvec selected = arma::sort_index(arma::abs(tempcors),
                                               "descend");
        selected = selected(arma::find(selected != i)); // removes correlation of j with j
        // selected has length d - 1
        selected = selected.head(k); 
        ngbrs.row(goodCols(i)) = goodCols(selected).t();
        robcors.row(goodCols(i)) = tempcors(selected).t();
      }
      break;
    }
      case 2: {// rank
        arma::mat Zr(U.n_rows, goodCols.size(), arma::fill::zeros);
        for (unsigned int i = 0; i < goodCols.size(); i++) {
          arma::vec u = Z.col(goodCols(i));
          arma::uvec finiteinds = arma::find_finite(u);
          arma::vec ufin = u(finiteinds);
          u.zeros();
          u(finiteinds) = LocScaleEstimators::rank(ufin);
          u(finiteinds) = u(finiteinds) - arma::mean(u(finiteinds));
          Zr.col(i) = u;
        }
        arma::mat corR = arma::cor(Zr);
        
        
        
        for (unsigned int i = 0; i < goodCols.size(); i++) {
          arma::vec tempcors = corR.col(i);
          arma::uvec selected = arma::sort_index(arma::abs(tempcors), "descend");
          // selected = selected.tail(U.n_cols - 1); // removes correlation of j with j
          selected = selected(arma::find(selected != i)); // removes correlation of j with j
          // selected has length d - 1
          selected = selected.head(k); 
          ngbrs.row(goodCols(i)) = goodCols(selected).t();
          robcors.row(goodCols(i)) = tempcors(selected).t();
        }
        break;
      }
      case 3: {// gkwls
        arma::mat Ugood = U.cols(goodCols);
        for (unsigned int i = 0; i < goodCols.size(); i++) {
        DDC::kbestcorr tempresult = DDC::kBestCorr(Ugood.col(i), Ugood, i,
                                                   k, qCorr, precScale);
        ngbrs.row(goodCols(i)) = goodCols(tempresult.selected).t();
        robcors.row(goodCols(i)) = tempresult.corrs.t();
      }
        
      }
      }
    } else {
      arma::mat Xgood = X.cols(goodCols);
      arma::vec locGood = locscaleX.loc(goodCols);
      arma::vec scaleGood = locscaleX.scale(goodCols);
      
      DDC::fastRobCorout nn2result = DDC::FastRobCorActual(Xgood, locGood, 
                                                           scaleGood, k,
                                                           qdim, nCorr, absCorr,
                                                           transFun, precScale,
                                                           treetype,
                                                           searchtype, radius,
                                                           eps, 0);
      
      if (goodCols.size() < X.n_cols) { // correct ngbrs if there are columns with > 50% NAs
        for (unsigned int i = 0; i < goodCols.size(); i++) {
          ngbrs.row(goodCols(i)) = goodCols(nn2result.ngbrs.row(i).t()).t();
          robcors.row(goodCols(i)) = nn2result.robcorrs.row(i);
        }
      } else {
        ngbrs.rows(goodCols) = nn2result.ngbrs;
        robcors.rows(goodCols) = nn2result.robcorrs;
      }
    }
    
    
    arma::mat corrweight = arma::abs(robcors); // should have no NAs
    
    // corrweight.cols(overHalfNA).zeros(); // make columns with over half NA standalone
    
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
    
    arma::uvec colStandalone = arma::find(arma::sum(corrweight, 1) == 0);
    
    arma::uvec colConnected = DDC::vdiff(arma::regspace<arma::uvec>(0, X.n_cols - 1),
                                         colStandalone);
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
        Zest.col(colConnected(i)) = DDC::predictCol(U.col(colConnected(i)),
                 U, colConnected(i),
                 ngbrs, corrweight, robslopes, combinRule);
      }
      
      ////////////////////////////////////
      //    STEP 5 : DESHRINKAGE        //
      ////////////////////////////////////
      
      // Deshrinkage : rescale Zest[, j] using robSlope of Z[, j] on Zest[, j]
      
      deshrinkage = arma::zeros(colConnected.size());
      
      for (unsigned int i = 0; i < colConnected.size(); i++) {
        deshrinkage(i) = DDC::deShrink(Zest.col(colConnected(i)),
                    Z, colConnected(i), qRegr, precScale);
        Zest.col(colConnected(i)) = Zest.col(colConnected(i)) * deshrinkage(i);
      }
      
      // Finally, all NAs are replaced by zeroes :
      Zest(find_nonfinite(Zest)).zeros();
      
      ////////////////////////////////////
      //    STEP 6 : FLAGGING CELLS     //
      ////////////////////////////////////
      
      // Compute cell residuals :
      Zres = Z - Zest; // original minus estimated
      Zres.cols(colStandalone) = Z.cols(colStandalone);
      
      scalest = arma::zeros(colConnected.size());
      
      for (unsigned int i = 0; i < colConnected.size(); i++)
      {
        scalest(i) = LocScaleEstimators::scale1StepM(Zres.col(colConnected(i)), 
                LocScaleEstimators::rhoHuber25,
                arma::datum::nan, precScale);
        Zres.col(colConnected(i)) = Zres.col(colConnected(i)) / scalest(i);
      }
      
      // We don't have to scale the standalone columns, as they
      // were already standardized in the beginning.
      // Next, flag outlying cells by their large residuals :
      indcells = find(arma::abs(Zres) > qCell); // does not flag the NAs as cells
      U(indcells).fill(arma::datum::nan);
      
    } // ends the iteration
    
    
    indcells = DDC::vdiff(indcells, DDC::col2cell(colStandalone, X.n_rows));
    // are the indices of outlying cells in connected variables only
    
    indcells = arma::sort(arma::unique(arma::join_cols(indcells, indexStandalone)));
    // are the indices of both types of outlying cells
    
    ////////////////////////////////////
    //    STEP 7 : FLAGGING ROWS      //
    ////////////////////////////////////
    
    arma::vec Ti(X.n_rows, arma::fill::zeros);
    arma::uvec indrows;
    arma::uvec indall = indcells;
    
    double medTi = 0;
    double madTi = 1;
    
      
      for (unsigned int i = 0; i < X.n_rows; i++) {
        arma::vec tempres = arma::erf(arma::sqrt(arma::pow(Zres.row(i), 2) / 2)).t();
        Ti(i) = arma::mean(tempres(arma::find_finite(tempres))) - 0.5;
      }
      
      // calculate the test value(outlyingness) of each row :
      arma::uvec finiteTi = arma::find_finite(Ti);
      medTi = arma::median(Ti(finiteTi));
      madTi = (1.482602218505602 * median(arma::abs(Ti(finiteTi) - medTi)));
      Ti = (Ti - medTi) / madTi;
      indrows = DDC::vinter(find_nonfinite(DDC::limitFilt(Ti, qRow)), arma::find(Ti > 0));
      indall = arma::unique(arma::join_cols(indcells, DDC::row2cell(indrows, X.n_rows, X.n_cols)));
    
    
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
    indcells = indcells + 1;
    indNAs   = indNAs + 1;
    indrows  = indrows + 1;
    ngbrs    = ngbrs + 1;
    indall   = indall + 1;
    colConnected = colConnected + 1;
    
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
                               Rcpp::Named("Z") = Z,
                               Rcpp::Named("locX") = locscaleX.loc,
                               Rcpp::Named("scaleX") = locscaleX.scale,
                               Rcpp::Named("deshrinkage") = deshrinkage,
                               Rcpp::Named("scalestres") = scalest,
                               Rcpp::Named("medTi") = medTi,
                               Rcpp::Named("madTi") = madTi,
                               Rcpp::Named("colConnected") = colConnected);
    
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
Rcpp::List estLocScale_cpp(arma::mat & X, unsigned int nLocScale, int type,  double precScale,
                           const int center, const double alpha) {
  // type 0 = biweight location + huber 1.5 scale
  // type 1 = huber location + huber scale
  // type 2 =  reweighted unimcd + wrapping location
  // type 3 =  reweighted unimcd
  // type 4 = reweighted raw mcd
  // type 5 = median/mad + wrapping 1 step M
  try
  {
    LocScaleEstimators::Xlocscale locscaleX;
    locscaleX = LocScaleEstimators::estLocScale(X, nLocScale, type,
                                                precScale, center, alpha); 
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



/*************************************/
/*       Unimcd Hidden export       */
/************************************/


// [[Rcpp::export]]
Rcpp::List unimcd_cpp(arma::vec & y, const double alpha) {
  try
  {
    LocScaleEstimators::locscale out;
    out = LocScaleEstimators::uniMcd(y, alpha);
    
    return(Rcpp::List::create( Rcpp::Named("loc") = out.loc,
                               Rcpp::Named("scale") = out.scale,
                               Rcpp::Named("weights") = out.weights));
    
  } catch( std::exception& __ex__ )
  {
    forward_exception_to_r( __ex__ );
  } catch(...)
  {
    ::Rf_error( "c++ exception " "(unknown reason)" );
  }
  
  return Rcpp::wrap(NA_REAL);
}


