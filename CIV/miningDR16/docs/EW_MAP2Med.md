# TODO: Switch from MAP to Median/Mean for REW comparison with SALSA

## 1. Extract posterior summaries
- [ ] Compute `b_median`, `logN_median` from posterior samples
- [ ] Compute `b_mean`, `logN_mean` from posterior samples
- [ ] Keep MAP values for comparison (do not delete yet)

## 2. Regenerate absorption profiles
For each absorber:
- [ ] Generate transmission profile A(λ) using MAP parameters
- [ ] Generate transmission profile A(λ) using median parameters
- [ ] Generate transmission profile A(λ) using mean parameters
- [ ] Apply same LSF convolution and same pixel grid in all cases

## 3. Recompute REW consistently
- [ ] Use the same wavelength/velocity integration window for all cases
- [ ] Compute:
  - REW_MAP
  - REW_median
  - REW_mean

## 4. Compare internally
- [ ] Plot REW_median / REW_MAP
- [ ] Plot REW_mean / REW_MAP
- [ ] Check if systematic shift exists (especially for strong/saturated lines)

## 5. Compare to SALSA
- [ ] Compare REW_median vs SALSA REW
- [ ] Compare REW_mean vs SALSA REW
- [ ] Check which summary reduces bias relative to SALSA

## 6. Decide final estimator
- [ ] If posterior is skewed → prefer median
- [ ] If posterior is symmetric → mean is acceptable
- [ ] Avoid MAP unless posterior is clearly Gaussian
