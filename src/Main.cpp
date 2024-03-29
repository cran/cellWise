

#include "DDC.h"
#include "cellMCD.h"

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
                   arma::uvec & goodCols,
                   const int & fixedCenter,
                   const arma::vec & center)
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
    if (fixedCenter == 0) {
      if (fastDDC == 0) {
        locscaleX = LocScaleEstimators::estLocScale(X, nLocScale,
                                                    standType, precScale, 1);
      } else {
        locscaleX = LocScaleEstimators::estLocScale(X, nLocScale,
                                                    2, precScale, 1);
      }
      Z = X.each_row() - locscaleX.loc.t();
    } else {
      Z = X.each_row() - center.t();
      if (fastDDC == 0) {
        locscaleX = LocScaleEstimators::estLocScale(Z, nLocScale,
                                                    standType, precScale, 0);
      } else {
        locscaleX = LocScaleEstimators::estLocScale(Z, nLocScale,
                                                    2, precScale, 0);
      }
      locscaleX.loc = center;
    }
    
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
      double newloc = arma::mean(ufin);
      double newsca = arma::stddev(ufin);
      u(finiteinds) = ufin * (scale(i) / newsca) + (loc(i) - newloc * (scale(i) / newsca) ) ;
      
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
Rcpp::List estLocScale_cpp(arma::mat & X, unsigned int nLocScale,
                           int type,  double precScale,
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
Rcpp::List unimcd_cpp(arma::vec & y, const double alpha, 
                      const int center) {
  try
  {
    LocScaleEstimators::locscale out;
    out = LocScaleEstimators::uniMcd(y, alpha, center);
    
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


// [[Rcpp::export]]
Rcpp::List findCellPath_cpp(arma::mat & predictors,
                            arma::vec & response,
                            arma::vec & weights,
                            arma::mat & Sigmai,
                            const arma::uvec & naMask) {
  
  try
  {
    
    // prepare input to lar ols by scaling the predictors/response
    // with the weights in w
    arma::mat x = predictors;
    x = x.each_row() / weights.t();
    arma::vec y = response;
    arma::mat Gram = Sigmai;
    Gram = Gram.each_row() / weights.t();
    Gram = Gram.each_col() / weights;
    
    //// start of lar_ols
    
    const double precScale = 1e-10;
    arma::uvec inactive =  arma::regspace<arma::uvec>(0, y.size() - 1);
    arma::uvec ignores  =  arma::uvec();
    
    arma::vec Cvec = (y.t() * x).t();
    arma::vec residuals = y;
    const int maxSteps = std::min(x.n_rows, x.n_cols);
    
    // containers for used output
    arma::mat beta = arma::mat(maxSteps + 1, x.n_cols, arma::fill::zeros);
    arma::vec RSS(maxSteps + 1 , arma::fill::zeros);
    RSS(0) = arma::conv_to<double>::from(arma::sum(arma::pow(residuals, 2))); // initial MD^2
    arma::cube bmat  = arma::cube(maxSteps + 1, x.n_cols, x.n_cols, arma::fill::zeros);
    arma::uvec path;
    arma::rowvec betaOLS(x.n_cols);
    
    //others
    arma::uvec active;
    arma::ivec actions = arma::ivec(maxSteps);
    int totactions = 0;
    int totPosactions = 0;
    // to indicate multiple variables entering the model at once:
    arma::uvec nbactions = arma::uvec(maxSteps, arma::fill::ones);
    arma::vec Sign;
    
    // matrix holding R of Cholesky decomposition of Sigmai[active, active]
    arma::mat R(x.n_cols,x.n_cols, arma::fill::zeros); 
    int rankR = 0;
    int dimR = 0;
    int k = 0;
    double Cmax = 1;
    double normXnew;
    arma::uvec newVars;
    arma::uvec oldVars;
    
    
    while ((k < maxSteps) && (active.size() < x.n_cols - ignores.size())) {
      arma::ivec action;
      arma::uvec Posaction;
      arma::vec C = Cvec(inactive);//correlations of the inactive variables
      Cmax = arma::max(arma::abs(C));
      
      
      if (Cmax < 1e-12) {
        // correlation of zero, exit the while loop
        break; 
      }
      
      // select variables going active and staying inactive
      if ((k == 0) && (arma::sum(naMask)) > 0) { 
        newVars = arma::find(naMask == 1);
      } else {
        newVars = arma::find(arma::abs(C) >= Cmax - precScale);
      }
      
      oldVars = DDC::vdiff(arma::regspace<arma::uvec>(0, C.size() - 1), newVars);
      
      
      C = C(oldVars);
      newVars = inactive(newVars);
      
      for(unsigned int i = 0; i < newVars.size(); i++) {
        // update R, the Cholesky decomposition of Sigmai[active, active]
        
        double xtx = Gram(newVars(i), newVars(i));
        normXnew = std::sqrt(xtx);
        if (dimR == 0) {
          R(0,0) = normXnew;
          rankR = 1;
          dimR = 1;
        } else {
          arma::vec Xtx = Gram.row(newVars(i)).t();
          Xtx = Xtx.elem(active);
          arma::mat r = arma::solve(arma::trimatl((R.submat(0, 0,dimR - 1,
                                                            dimR - 1)).t()), Xtx);
          
          double rpp;
          rpp = std::pow(normXnew,2) - arma::as_scalar(arma::sum(arma::pow(r,2)));
          
          if (rpp <= precScale){
            rpp = precScale;
          } else {
            rpp = std::sqrt(rpp);
            rankR = rankR + 1;
          }
          dimR = dimR + 1;
          arma::vec newcol(x.n_cols, arma::fill::zeros);
          newcol.head(dimR - 1) = r;
          newcol(dimR - 1) = rpp;
          R.col(dimR - 1) = newcol;
          
        }
        if ((unsigned)rankR == active.size()) { // an ignore variable has entered
          dimR = dimR - 1;
          ignores.insert_rows(ignores.n_rows, newVars.row(i));
          action.insert_rows(action.n_rows, arma::conv_to<arma::ivec>::from(newVars.row(i) * (-1)));
        } else {
          active.insert_rows(active.n_rows, newVars.row(i));
          Posaction.insert_rows(Posaction.n_rows, newVars.row(i));
          path.insert_rows(path.n_rows, newVars.row(i));
          Sign.insert_rows(Sign.n_rows, arma::sign(Cvec.row(newVars(i))));
          action.insert_rows(action.n_rows, arma::conv_to<arma::ivec>::from(newVars.row(i)));
        }
      }
      
      arma::mat Gi1 = arma::solve(arma::trimatu(R.submat(0, 0,dimR - 1, dimR - 1)),
                                  (arma::solve(arma::trimatl((R.submat(0, 0,dimR - 1,
                                                                       dimR - 1)).t()),
                                                                       Sign)));
      
      double A = 1 / std::sqrt(arma::sum(Gi1 % Sign));
      arma::vec  w = A * Gi1; // need to cast matrix to vector
      double gamhat;
      
      // construct the inactive set by taking out actives and ignores
      
      arma::uvec actAndIgn = active;
      actAndIgn.insert_rows(actAndIgn.n_rows, ignores);
      inactive = DDC::vdiff(arma::regspace<arma::uvec>(0, x.n_cols - 1),
                            arma::sort(actAndIgn));
      
      
      if ((active.size() >= x.n_cols - ignores.size())) {
        gamhat = Cmax / A;
      } else {
        arma::vec a;
        a = (w.t() * Gram(active,inactive)).t();
        arma::vec gam(2* C.size());
        
        gam.head(C.size()) = (Cmax - C) / (A - a);
        gam.tail(C.size()) = (Cmax + C) / (A + a);
        
        gam = gam.elem(find(gam > precScale));
        if (gam.size() > 0) {
          gamhat = std::min(Cmax/A, arma::min(gam)); 
        } else {
          gamhat = Cmax/A;
        }
      }
      
      Cvec = Cvec - gamhat * Gram.cols(active) * w;
      
      betaOLS.zeros();
      betaOLS(active) =  (arma::solve(arma::trimatu(R.submat(0, 0,dimR - 1, dimR - 1)),
                          (arma::solve(arma::trimatl((R.submat(0, 0,dimR - 1,
                                                               dimR - 1)).t()),
                                                               (x.cols(active)).t()))) * y).t();
      residuals = y -  (betaOLS.cols(active) *  (x.cols(active)).t()).t();
      
      arma::mat Ip(dimR, dimR);
      Ip.eye();
      arma::mat RI = arma::solve(arma::trimatu(R.submat(0, 0, dimR - 1, dimR - 1)), Ip); //inverse
      arma::mat bmat_replacement(x.n_cols, x.n_cols, arma::fill:: zeros);
      bmat_replacement(active, active) =  RI * RI.t();
      
      bmat_replacement = bmat_replacement.each_row() / weights.t();
      bmat_replacement = bmat_replacement.each_col() / weights;
      for(unsigned int i = 0; i < Posaction.size(); i++) {
        bmat.row(totPosactions + i + 1) = bmat_replacement; // the first one is a zero matrix
        beta.row(totPosactions + i + 1) = betaOLS;
        RSS(totPosactions + i + 1) = arma::conv_to<double>::from(arma::sum(arma::pow(residuals, 2)));
      }
      
      actions.subvec(totactions, totactions+action.size() - 1) = action;
      nbactions(k) = action.size();
      totactions = totactions + action.size();
      totPosactions = totPosactions + Posaction.size();
      
      // Now return the biasmat. 
      // we return the R of the cholesky decomposition instead
      // can be used to construct bias may through Rinverse * Rinverse.t();
      // arma::mat Ip(dimR - 1, dimR - 1);
      // Ip.eye();
      // arma::mat RI = arma::solve(arma::trimatu(R.submat(0, 0, dimR - 1, dimR - 1)), Ip); //inverse
      // arma::mat bmat_replacement(x.n_cols, x.n_cols, arma::fill:: zeros);
      // bmat_replacement(active, active) =  RI * RI.t();
      
      // arma::mat bmat_replacement(x.n_cols, x.n_cols, arma::fill:: zeros);
      // bmat_replacement.submat(0, 0, dimR - 1, dimR - 1) = R.submat(0, 0, dimR - 1, dimR - 1);
      // biasMat.row(k+1) = bmat_replacement; // the first one is a zero matrix
      // biasMat.row(k + 1) = R; // the first one is a zero matrix
      k = k + 1;
      
    }
    
    // now need to check whether beta matrices etc. are complete
    
    arma::uvec missingInds;
    missingInds = DDC::vdiff(arma::regspace<arma::uvec>(0, x.n_cols - 1),
                             arma::sort(path));
    if (missingInds.size() > 0) {
      path.insert_rows(path.size(), missingInds);
      betaOLS.zeros();
      betaOLS   = arma::solve(x, y, arma::solve_opts::likely_sympd).t();
      residuals = y -  (betaOLS *  x.t()).t();
      
      for (unsigned int i = 0; i < missingInds.size(); i++) {
        // bmat.row(bmat.n_rows - i - 1)  = solve(Sigmai); // maybe we don't need these?
        beta.row(beta.n_rows - i - 1) = betaOLS;
        RSS(RSS.size() - i - 1) = arma::conv_to<double>::from(arma::sum(arma::pow(residuals, 2)));
      }
    }
    //// End of lar_ols
    // Now start constructing the output
    
    // arma::mat beta = arma::mat(maxSteps + 1, x.n_cols);
    // arma::vec RSS(maxSteps, arma::fill::zeros);
    // RSS_sw(0) = arma::conv_to<double>::from(arma::sum(arma::pow(residuals, 2))); // initial MD^2
    // arma::cube bmat  = arma::cube(maxSteps + 1, x.n_cols, x.n_cols, arma::fill::zeros);
    // arma::vec ordering;
    // 
    // int nbignores = ignores.size();
    // 
    // arma::uvec path;
    // pathInd = 0
    // for (unsigned int i = 0; i < actions.size(); i++) {
    //   if(actions(i) > 0){
    //     path(pathInd) = (unsigned int) actions(i);
    //     
    //     pathInd = pathInd + 1;
    //   } else if (actions(i) == 0) {
    //     if (arma::find(ignores == 0).size() == 0) {
    //       path(pathInd) = 0;
    //       pathInd = pathInd + 1;
    //     } 
    //   } 
    // }
    
    // unscale with the weights, bmat already has been rescaled
    beta = beta.each_row() / weights.t();
    
    // converg path to 1-based indexing:
    path = path + 1;
    
    return(Rcpp::List::create( Rcpp::Named("beta") = beta,
                               Rcpp::Named("biasMat") = bmat,
                               Rcpp::Named("ordering") = path,
                               Rcpp::Named("RSS") = RSS));
    
  } catch( std::exception& __ex__ )
  {
    forward_exception_to_r( __ex__ );
  } catch(...)
  {
    ::Rf_error( "c++ exception " "(unknown reason)" );
  }
  
  return Rcpp::wrap(NA_REAL);
}


// cellMCD //

// [[Rcpp::export]]
Rcpp::List allpreds_cpp(arma::mat& X,
                        arma::mat& S,
                        arma::vec& mu,
                        arma::umat& W) {
  try
  {
    
    arma::uword n = W.n_rows;
    arma::uword d = W.n_cols;
    //arma::mat Sinv = arma::inv_sympd(arma::symmatu(S));
    arma::mat Sinv = arma::inv(S);
    arma::mat preds(n, d);
    arma::mat cvars(n, d);
    
    
    std::unordered_map<std::string, arma::uvec> Wmap;
    uniqueRows(W, Wmap);
    
    
    for (std::pair<std::string, arma::uvec> s : Wmap) {
      
      
      arma::uvec myrows = s.second; // all indices with rowpattern(string) s
      arma::urowvec w = W.row(myrows(0)); // the actual rowpattern
      
      arma::uvec obs = arma::find(w);
      arma::uvec mis = arma::find(1 - w);
      
      if (mis.n_elem < d) {
        if (mis.n_elem > 0) {
          
          arma::mat Sjinv = subinverse_cpp(S, Sinv, obs);
          
          arma::vec cvar = arma::diagvec(S(mis, mis) -
            S(mis, obs) * Sjinv * S(obs, mis));
          for (arma::uword i = 0; i < myrows.n_elem; i++) {
            arma::uvec rowi(1);
            rowi(0) = myrows(i);
            cvars(rowi, mis) = cvar.t();
            preds(rowi, mis) = mu(mis).t() + (X(rowi, obs) - mu(obs).t()) * Sjinv * S(obs, mis);
          }
          
          // now compute preds and cvars for all available entries
          //
          
          if (obs.n_elem == 1) {
            for (arma::uword i = 0; i < myrows.n_elem; i++) {
              arma::uvec rowi(1);
              rowi(0) = myrows(i);
              preds(rowi, obs) = mu(obs);
              cvars(rowi, obs) = S(obs, obs);
            }
          } else {
            for (arma::uword j = 0; j < obs.n_elem; j++) {
              arma::uvec obsj(1);
              obsj(0) = obs(j);
              arma::uvec helps = obs(arma::find(obs != obs(j)));
              arma::mat M2inv = subinverse_cpp(S, Sinv, helps);
              
              arma::mat cvar = S(obsj, obsj) - S(obsj, helps) *
                M2inv * S(helps, obsj);
              
              for (arma::uword i = 0; i < myrows.n_elem; i++) {
                arma::uvec rowi(1);
                rowi(0) = myrows(i);
                cvars(rowi, obsj) = cvar;
                preds(rowi, obsj) = mu(obsj).t() +
                  (X(rowi, helps) - mu(helps).t()) * M2inv * S(helps, obsj);
              }
            }
          }
          
        } else { //no missings in this pattern
          
          arma:: vec cvar = (1 / arma::diagvec(Sinv));
          
          for (arma::uword i = 0; i < myrows.n_elem; i++) {
            cvars.row(myrows(i)) =  cvar.t();
          }
          
          for (arma::uword j = 0; j < d; j++) {
            arma::uvec jvec(1);
            jvec(0) = j;
            arma::uvec j_neg = arma::regspace<arma::uvec>(0, d - 1);
            j_neg.shed_rows(jvec);
            
            arma::mat Sjinv = subinverse_cpp(S, Sinv, j_neg);
            arma::mat Xtemp = X(myrows, j_neg);
            Xtemp.each_row() -= mu(j_neg).t();
            
            preds(myrows, jvec) = (arma::as_scalar(mu(j)) +
              Xtemp * Sjinv * S(j_neg, jvec));
          }
        }
      } else { // all elements of this pattern are missing
        for (arma::uword i = 0; i < myrows.n_elem; i++) {
          preds.row(myrows(i)) = mu.t();
          cvars.row(myrows(i)) =  arma::diagvec(Sinv).t();
        }
      }
      
    }
    
    
    return Rcpp::List::create(Rcpp::Named("preds") = preds,
                              Rcpp::Named("cvars") = cvars);
  } catch( std::exception& __ex__ )
  {
    forward_exception_to_r( __ex__ );
  } catch(...)
  {
    ::Rf_error( "c++ exception " "(unknown reason)" );
  }
  
  return Rcpp::wrap(NA_REAL);
}



// [[Rcpp::export]]
double Objective_cpp(arma::mat &X,
                     arma::umat &W,
                     arma::vec &mu,
                     arma::mat &Sigma,
                     arma::mat &Sigmai) {
  try
  {
    
    double objective = 0;
    double objective_add = 0;
    
    
    std::unordered_map<std::string, arma::uvec> Wmap;
    uniqueRows(W, Wmap);
    
    
    for (std::pair<std::string, arma::uvec> s : Wmap) {
      
      arma::uvec myrows = s.second; // all indices with rowpattern(string) s
      arma::urowvec w = W.row(myrows(0)); // the actual rowpattern
      
      arma::uvec obs = arma::find(w);
      arma::uvec mis = arma::find(1 - w);
      
      arma::vec subMu = mu(obs);
      arma::mat subSigma = Sigma(obs, obs);
      arma::mat invsubSigma = subinverse_cpp(Sigma, Sigmai, obs);
      
      
      arma::mat Xtemp = X(myrows, obs);
      Xtemp.each_row() -= subMu.t();
      arma::mat MD = arma::sum((Xtemp * invsubSigma) % Xtemp, 1);
      objective_add = arma::as_scalar(arma::sum(MD + 
        arma::log_det_sympd(arma::symmatu(subSigma)) + 
        std::log(2 * arma::datum::pi) * (obs.n_elem)));
      objective    = objective + objective_add;
    } 
    
    return(objective);
  } catch( std::exception& __ex__ )
  {
    forward_exception_to_r( __ex__ );
  } catch(...)
  {
    ::Rf_error( "c++ exception " "(unknown reason)" );
  }
  
  return NA_REAL;
}


// [[Rcpp::export]]
arma::umat updateW_cpp(const arma::mat &X,
                       arma::umat W,
                       const arma::vec &mu,
                       const arma::mat &Sigma,
                       const arma::mat &Sigmai,
                       const arma::vec &lambda,
                       const arma::uword &h) {
  try
  {
    
    arma::uword n = W.n_rows;
    arma::uword d = W.n_cols;
    arma::uvec ord = arma::stable_sort_index(arma::sum(W, 0));
    arma::vec Delta(n);
    arma::uvec goodCells;
    arma::uvec wnew(n);
    arma::uvec Delta_rank(n);
    for (arma::uword j = 0; j < d; j++) {
      arma::uword jtemp = ord(j);
      Delta = Deltacalc_cpp(X, W, Sigma, Sigmai, mu, jtemp + 1);
      goodCells = arma::find(Delta <= lambda(jtemp));
      if (goodCells.n_elem < h) {
        Delta_rank = arma::stable_sort_index(Delta - lambda(jtemp));
        goodCells = Delta_rank.head(h);
      }
      wnew.zeros();
      wnew(goodCells).fill(1);
      W.col(jtemp) = wnew;
    }
    return(W);
  } catch( std::exception& __ex__ )
  {
    forward_exception_to_r( __ex__ );
  } catch(...)
  {
    ::Rf_error( "c++ exception " "(unknown reason)" );
  }
  arma::umat error_out(1,1);
  error_out.fill(arma::datum::nan);
  return error_out;
}
