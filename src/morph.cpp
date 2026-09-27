#include <Rcpp.h>
using namespace Rcpp;

// Morphological operator for gridded objects (c++)
//
// Originally written for loom, inspired by R::mmand::morph by Jon
// Clayden. Ported to geomio to operate on flat vectors (row-major) instead
// of R matrices, avoiding unnecessary memory copies.

// --- Blend functions ---
// Pairwise combination of neighborhood values with kernel values

NumericVector blendIdentity(NumericVector &temp, NumericVector &kernel){
  if(any(is_na(kernel))){
    for(int i = 0; i < temp.size(); i++){
      if(is_na(kernel)[i]){
        temp[i] = NA_REAL;
      }
    }
  }
  return(temp);
}

NumericVector blendEqual(NumericVector &temp, NumericVector &kernel){
  NumericVector out(1);
  out[0] = 1.0;
  for(int i = 0; i < temp.size(); i++){
    if(!is_na(kernel)[i] & (temp[i] != kernel[i])){
      out[0] = 0.0;
      break;
    }
  }
  return(out);
}

// Blends that drop kernel-NA cells collect the kept elements into a new vector.
// Erasing in place while the loop index advances would skip the element that
// slides into the vacated position, one skip per NA.

NumericVector blendLower(NumericVector &temp, NumericVector &kernel){
  NumericVector out(temp.size());
  int n = 0;
  for(int i = 0; i < temp.size(); i++){
    if(NumericVector::is_na(kernel[i])) continue;
    out[n++] = (temp[i] > kernel[i]) ? 0.0 : temp[i];
  }
  return(head(out, n));
}

NumericVector blendGreater(NumericVector &temp, NumericVector &kernel){
  NumericVector out(temp.size());
  int n = 0;
  for(int i = 0; i < temp.size(); i++){
    if(NumericVector::is_na(kernel[i])) continue;
    out[n++] = (temp[i] < kernel[i]) ? 0.0 : temp[i];
  }
  return(head(out, n));
}

NumericVector blendPlus(NumericVector &temp, NumericVector &kernel){
  NumericVector out(temp.size());
  int n = 0;
  for(int i = 0; i < temp.size(); i++){
    if(NumericVector::is_na(kernel[i])) continue;
    out[n++] = temp[i] + kernel[i];
  }
  return(head(out, n));
}

NumericVector blendMinus(NumericVector &temp, NumericVector &kernel){
  NumericVector out(temp.size());
  int n = 0;
  for(int i = 0; i < temp.size(); i++){
    if(NumericVector::is_na(kernel[i])) continue;
    out[n++] = temp[i] - kernel[i];
  }
  return(head(out, n));
}

NumericVector blendProduct(NumericVector &temp, NumericVector &kernel){
  NumericVector out(temp.size());
  int n = 0;
  for(int i = 0; i < temp.size(); i++){
    if(NumericVector::is_na(kernel[i])) continue;
    out[n++] = temp[i] * kernel[i];
  }
  return(head(out, n));
}

// --- Merge functions ---
// Reduce blended values to a single output value

double mergeMin(const NumericVector &values){
  return(min(na_omit(values)));
}

double mergeMax(const NumericVector &values){
  return(max(na_omit(values)));
}

double mergeAll(const NumericVector &values){
  return all(na_omit(values) != 0).is_true() ? 1.0 : 0.0;
}

double mergeAny(const NumericVector &values){
  return any(na_omit(values) != 0).is_true() ? 1.0 : 0.0;
}

double mergeNotAll(const NumericVector &values){
  return all(na_omit(values) != 0).is_true() ? 0.0 : 1.0;
}

double mergeNotAny(const NumericVector &values){
  return any(na_omit(values) != 0).is_true() ? 0.0 : 1.0;
}

double mergeSum(const NumericVector &values){
  return(sum(na_omit(values)));
}

double mergeNA(const NumericVector &values){
  if(sum(na_omit(values)) != 0){
    return(sum(values));
  } else{
    return(NA_REAL);
  }
}

double mergeMean(const NumericVector &values){
  return(mean(na_omit(values)));
}

double mergeMedian(const NumericVector &values){
  return(median(na_omit(values)));
}

double mergeSD(const NumericVector &values){
  return(sd(na_omit(values)));
}

double mergeCV(const NumericVector &values){
  return(sd(na_omit(values)) / mean(na_omit(values)));
}

// --- Function pointer tables ---

// [[Rcpp::plugins(cpp11)]]

using BlendFunctionPtr = NumericVector (*)(NumericVector &, NumericVector &);
using MergeFunctionPtr = double (*)(const NumericVector &);

static BlendFunctionPtr blendFuns[7] = {
  &blendIdentity, &blendEqual, &blendLower, &blendGreater,
  &blendPlus, &blendMinus, &blendProduct
};
static MergeFunctionPtr mergeFuns[12] = {
  &mergeMin, &mergeMax, &mergeAll, &mergeAny, &mergeNotAll, &mergeNotAny,
  &mergeSum, &mergeMean, &mergeMedian, &mergeSD, &mergeCV, &mergeNA
};

//' Morphological operation on a gridded object (c++)
//'
//' C++ function that applies a morphological operation to a gridded object.
//' The function slides a kernel over the grid, blends the neighborhood values
//' with the kernel, and merges the result into a single output value per cell.
//' @param vals [numeric(.)][numeric]\cr flat vector of cell values
//'   (row-major order).
//' @param valRows [integer(1)][integer]\cr number of rows in the grid.
//' @param valCols [integer(1)][integer]\cr number of columns in the grid.
//' @param kernel [matrix(numeric)][matrix]\cr the structuring element. NA
//'   cells in the kernel are ignored.
//' @param value [numeric(.)][numeric]\cr which cell values to process. Cells
//'   whose value is not in this vector are left unchanged.
//' @param blend [integer(1)][integer]\cr blend mode (1-indexed): 1=identity,
//'   2=equal, 3=lower, 4=greater, 5=plus, 6=minus, 7=product.
//' @param merge [integer(1)][integer]\cr merge mode (1-indexed): 1=min,
//'   2=max, 3=all, 4=any, 5=!all, 6=!any, 7=sum, 8=mean, 9=median, 10=sd,
//'   11=cv, 12=sumNa.
//' @param rotateKernel [logical(1)][logical]\cr whether to try all 4 rotations
//'   of the kernel (useful for hit-or-miss transforms).
//' @param strictKernel [logical(1)][logical]\cr whether the kernel must fit
//'   entirely within the grid (TRUE) or is clipped at borders (FALSE).
//' @return A numeric vector of the same length as \code{vals} with the
//'   morphological operation applied.
//' @export
// [[Rcpp::export]]
NumericVector morphCpp(NumericVector &vals, int valRows, int valCols,
                       NumericMatrix &kernel, NumericVector &value,
                       int blend, int merge,
                       bool rotateKernel, bool strictKernel) {

  int kRows = kernel.nrow(), kCols = kernel.ncol();
  int yMar = kRows / 2, xMar = kCols / 2;

  NumericVector out = clone(vals);
  NumericVector tempMat, tempKernel;
  NumericMatrix rotatedKernel, kernelToRotate;
  double merged;

  BlendFunctionPtr blendFun = blendFuns[blend - 1];
  MergeFunctionPtr mergeFun = mergeFuns[merge - 1];

  for(int y = 0; y < valRows; y++){
    for(int x = 0; x < valCols; x++){

      int idx = y * valCols + x;

      // skip if target value is not in the allowed values
      if(NumericVector::is_na(vals[idx])) continue;
      if(!any(vals[idx] == value).is_true()) continue;

      int mL = x - xMar;
      int mR = x + xMar;
      int kL = 0;
      int mT = y - yMar;
      int mB = y + yMar;
      int kT = 0;

      // handle borders
      if(strictKernel){
        // skip cells where the kernel doesn't fit entirely
        if(mL < 0 || mR > valCols - 1 || mT < 0 || mB > valRows - 1) continue;
      } else {
        if(mL < 0){
          kL = xMar - x;
          mL = 0;
        }
        if(mR > valCols - 1){
          mR = valCols - 1;
        }
        if(mT < 0){
          kT = yMar - y;
          mT = 0;
        }
        if(mB > valRows - 1){
          mB = valRows - 1;
        }
      }

      // extract neighborhood from flat vector into a temporary vector
      int nRows = mB - mT + 1;
      int nCols = mR - mL + 1;
      tempMat = NumericVector(nRows * nCols);
      tempKernel = NumericVector(nRows * nCols);
      int ti = 0;
      for(int cx = mL; cx <= mR; cx++){
        for(int ry = mT; ry <= mB; ry++){
          tempMat[ti] = vals[ry * valCols + cx];
          tempKernel[ti] = kernel(ry - mT + kT, cx - mL + kL);
          ti++;
        }
      }

      merged = mergeFun(blendFun(tempMat, tempKernel));

      // rotate kernel if needed
      if(rotateKernel && merged != 1.0){
        rotatedKernel = clone(kernel);
        for(int r = 0; r < 4; r++){
          kernelToRotate = clone(rotatedKernel);
          for(int m = 0; m < kCols; m++){
            rotatedKernel(m, _) = rev(kernelToRotate(_, m));
          }
          // re-extract kernel with rotation
          ti = 0;
          for(int cx = mL; cx <= mR; cx++){
            for(int ry = mT; ry <= mB; ry++){
              tempKernel[ti] = rotatedKernel(ry - mT + kT, cx - mL + kL);
              ti++;
            }
          }
          // re-extract mat (blendFun may have modified it)
          ti = 0;
          for(int cx = mL; cx <= mR; cx++){
            for(int ry = mT; ry <= mB; ry++){
              tempMat[ti] = vals[ry * valCols + cx];
              ti++;
            }
          }
          merged = mergeFun(blendFun(tempMat, tempKernel));
          if(merged == 1.0) break;
        }
      }

      out[idx] = merged;
    }
  }

  return(out);
}
