# Equivalent Width (EW) Calculation – TODO / Fix List

## 🚨 Critical Bugs to Fix

- [ ] Fix 1548 integration window condition
  - Ensure window is symmetric in velocity around `z_abs`
  - Remove asymmetric `< z_abs + dz_mid` logic
  - Verify no unintended clipping of line wings

- [ ] Replace `(1+z_qso)` with `(1+z_abs)` where appropriate
  - Especially in Doppler width → redshift conversions
  - All absorber-frame calculations must use `z_abs`

- [ ] Fix 1550 Voigt EW indexing bug
  - Ensure 1550 integration uses `absorptionL1_1550`
  - Avoid reusing 1548 arrays accidentally

- [ ] Fix error array indexing
  - `ErrREW_1548_flux` and `ErrREW_1550_flux` must include `jj` index
  - Currently overwriting 2D slice of 3D array
  - Ensure dimensions match `sigma_width` loop

- [ ] Fix sigma_width dimension mismatch
  - Loop uses `[1,2,3,4,5]`
  - Arrays preallocated with 3rd dimension = 4
  - Make consistent

---

## ⚠️ Methodology Improvements

- [ ] Define integration window in velocity space
  - Use symmetric ±Δv around absorber redshift
  - Apply identical window to both doublet lines
  - Avoid redshift-based clipping logic

- [ ] Ensure EW computed in absorber rest frame
  - Confirm integration uses absorber rest wavelength
  - If integrating in observed frame, divide by `(1+z_abs)`

- [ ] Standardize window selection method
  - Choose between:
    - Fixed ± velocity window (e.g. ±300 km/s)
    - Contiguous detection threshold (e.g. 2σ below continuum)
  - Document choice clearly

---

## 📏 Numerical / Stability Improvements

- [ ] Ensure EW integration is area-preserving
  - Confirm rebinning uses bin-averaging, not interpolation
  - Ensure Δλ used in error matches trapz spacing

- [ ] Fix EW error formula
  - Use:
    ```
    sigma_W^2 = sum( (sigma_F ./ continuum).^2 .* (delta_lambda).^2 )
    ```
  - Avoid using `10^(pixel_spacing)` unless strictly correct

- [ ] Verify EW invariance under LSF convolution
  - Test:
    - EW on fine grid
    - EW after convolution + rebin
  - Should match within ~1%

---

## 🔬 Validation Tests to Add

- [ ] Compare:
  - Intrinsic model EW
  - Post-LSF EW
  - Post-rebin EW
  - Measured EW from noisy spectrum

- [ ] Inject synthetic Voigt lines with known N
  - Verify recovered EW matches analytic value

- [ ] Test sensitivity to integration window width

- [ ] Compare results to SALSA synthetic spectra using identical measurement pipeline

---

## 📊 Science-Level Consistency Checks

- [ ] Confirm EW(1548) / EW(1550) behaves physically
- [ ] Check doublet ratio consistency for unsaturated systems
- [ ] Quantify EW bias vs S/N
- [ ] Quantify EW bias vs resolution

---

## 📄 Documentation

- [ ] Explicitly document:
  - Integration range definition
  - Frame of EW (absorber rest vs QSO rest)
  - LSF handling
  - Error propagation method
  - Saturation flagging logic
